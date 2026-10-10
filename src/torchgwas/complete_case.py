"""Complete-case OLS for phenotypes with missing values, from the scan's full-sample products.

Each trait is tested on its own observed samples, as plink2's --glm does per
phenotype: intercept and covariates refitted on those samples, df their count
less the rank and the genotype. The scan's full-sample products give this
almost for free:

- The phenotype is residualized on [1, covariates] over its own observed rows
  and zero-filled elsewhere (preprocess.residualize_and_standardize,
  missing='complete_case'). Then Z_S^T y_S = 0, so the numerator g^T P_S y is
  the scan's own product g^T y, whatever g's centring.
- The call's residual sum over the subset S = all rows less the missing rows M
  is a downdate of full-sample quantities, with Z = [1/sqrt(n), Q]:

      g^T P_S g = (sum g^2 - sum_M g^2) - u^T (Z^T Z - Z_M^T Z_M)^+ u,
      u = Z^T g - Z_M^T g_M.

  Z^T g is the scan's covariate product, sum g^2 its residual sum's first term;
  only the M rows are gathered, so the work is (missing cells) x (rank + 1)
  per variant. Traits sharing a missingness pattern share it.
- df = (samples in S with the call observed) - rank(Z_S) - 1.

Calls missing inside S keep the scan's convention: centred at the variant's
observed mean (zero), df reduced by one each. That mean is over every scanned
sample, not over S alone, so with missing calls the statistic differs a
little from a run on S's samples only (0.3% of a JAGWAS chi2 at 3% of calls
missing and 2 of 150 subjects dropped; tests/test_jagwas_group_drop.py).
Making it S's mean needs the call-side sums over each variant's missing
calls, which the plan does not gather.

Measured (benchmarks/missing_phenotype_conventions_20260927.py, relative
error of -log10 P against FP64 lstsq on the observed rows): this formula 0
(machine precision); the release convention, mean-imputed t times
sqrt(trait_df / df), -1.3% median on the full-scale panel (up to 2% missing)
and -31% where 30% of samples are missing by a covariate.
"""
from __future__ import annotations

import numpy as np

# Gathered cells per variant in one device block (chunk x BLOCK_CELLS floats).
BLOCK_CELLS = 1 << 14
# Patterns per device group (its (chunk, patterns, rank) FP64 state).
PATTERN_BLOCK = 256
# Missing cells per batch of patterns whose Z_M^T Z_M are summed at once
# (cells x rank^2 FP64: 440 MB at rank 29).
DOWNDATE_CELLS = 1 << 16
# Traits per column block when a panel is read for its missing cells.
PANEL_BLOCK = 4096


def complete_case_basis(n_samples, q_matrix):
    """Z = [1/sqrt(n), Q]: the scan's covariate design, orthonormal (Q is centred)."""
    one = np.full((int(n_samples), 1), 1.0 / np.sqrt(n_samples))
    return one if q_matrix is None else np.column_stack([one, np.asarray(q_matrix, dtype=np.float64)])


def _missing_cells(missing_of_block, n_traits, block):
    """Every missing cell, by trait then row: missing_of_block(first, last) -> (samples, last - first) bool."""
    traits, rows = [np.empty(0, dtype=np.int64)], [np.empty(0, dtype=np.int64)]
    for first in range(0, n_traits, max(1, int(block))):
        last = min(n_traits, first + max(1, int(block)))
        trait, row = np.nonzero(np.ascontiguousarray(missing_of_block(first, last).T))
        traits.append(trait.astype(np.int64) + first)
        rows.append(row.astype(np.int64))
    return np.concatenate(traits), np.concatenate(rows)


class CompleteCasePlan:
    """The missing rows of each trait with missing phenotype values, grouped by pattern."""

    def __init__(self, observed, q_matrix):
        observed = np.asarray(observed, dtype=bool)
        if observed.ndim != 2:
            raise ValueError('observed must be (samples, traits)')
        n, k = observed.shape
        self._build(n, k, *_missing_cells(lambda a, b: ~observed[:, a:b], k, PANEL_BLOCK), q_matrix)

    def _build(self, n, k, cell_traits, cell_rows, q_matrix):
        """The plan from its missing cells, sorted by trait then row."""
        self.n_samples, self.n_traits = n, k
        self.basis = complete_case_basis(n, q_matrix)
        rank = self.basis.shape[1]
        self.traits, starts, counts = np.unique(cell_traits, return_index=True, return_counts=True)
        self.traits = self.traits.astype(np.int64)
        index, rows = {}, []
        self.pattern_of_trait = np.empty(self.traits.size, dtype=np.int64)
        for position, (start, count) in enumerate(zip(starts.tolist(), counts.tolist())):
            missing = cell_rows[start:start + count]
            key = missing.tobytes()
            pattern = index.get(key)
            if pattern is None:
                pattern = index[key] = len(rows)
                rows.append(missing)
            self.pattern_of_trait[position] = pattern
        self.rows = rows
        sizes = np.array([len(missing) for missing in rows], dtype=np.int64)
        offsets = np.concatenate([[0], np.cumsum(sizes)])
        flat = np.concatenate(rows) if rows else np.empty(0, dtype=np.int64)
        self.inverse = np.empty((len(rows), rank, rank), dtype=np.float64)
        self.subset_rank = np.empty(len(rows), dtype=np.int64)
        # Z^T Z itself, not I: a float32 basis is orthonormal only to ~1e-7,
        # which would lift a subset's null direction above the rank cut.
        gram = self.basis.T @ self.basis
        first = 0
        while first < len(rows):
            # Patterns [first, last): at most DOWNDATE_CELLS cells, or one larger pattern.
            last = int(np.searchsorted(offsets, offsets[first] + DOWNDATE_CELLS, side='right')) - 1
            last = min(len(rows), max(first + 1, last))
            zm = self.basis[flat[offsets[first]:offsets[last]]]
            downdate = np.add.reduceat(zm[:, :, None] * zm[:, None, :], offsets[first:last] - offsets[first], axis=0)
            values, vectors = np.linalg.eigh(gram[None] - downdate)
            # A covariate constant on the subset leaves Z_S rank-deficient;
            # its df then counts the rank Z_S actually has.
            keep = values > 1e-9 * np.maximum(values.max(axis=-1, keepdims=True), 1e-300)
            scale = np.where(keep, 1.0 / np.where(keep, values, 1.0), 0.0)
            self.inverse[first:last] = (vectors * scale[:, None, :]) @ vectors.transpose(0, 2, 1)
            self.subset_rank[first:last] = keep.sum(axis=-1)
            first = last
        self.observed_counts = n - sizes
        # Where each sample's column sits in the scan's calls (None: in order).
        self._positions = None
        self._device = {}

    @classmethod
    def from_panel(cls, phenotype, q_matrix, block=PANEL_BLOCK):
        """The plan for a raw panel (NaN missing), read a column block at a time, or None when it is complete.

        phenotype: an array or anything sliced as phenotype[:, first:last]
        (a memory-mapped panel or its views), so no samples x traits mask is
        held at once.
        """
        n, k = phenotype.shape
        traits, rows = _missing_cells(
            lambda a, b: np.isnan(np.asarray(phenotype[:, a:b], dtype=np.float64)), k, block)
        if not traits.size:
            return None
        plan = cls.__new__(cls)
        plan._build(n, k, traits, rows, q_matrix)
        return plan

    def at_positions(self, positions):
        """This plan for calls laid out as whole rows: sample i's column is positions[i]."""
        import copy
        placed = copy.copy(self)
        placed._positions = np.asarray(positions, dtype=np.int64)
        placed._device = {}
        return placed

    @classmethod
    def from_phenotype(cls, phenotype, q_matrix):
        """The plan for a raw panel (NaN missing), or None when it is complete."""
        return cls.from_panel(np.asarray(phenotype, dtype=np.float64), q_matrix)

    def residualize(self, values, trait):
        """Trait `trait`'s complete-case residuals: its own observed rows, zero elsewhere, unit variance there."""
        position = int(np.searchsorted(self.traits, trait))
        pattern = self.pattern_of_trait[position]
        keep = np.ones(self.n_samples, dtype=bool)
        keep[self.rows[pattern]] = False
        z = self.basis[keep]
        y = np.asarray(values, dtype=np.float64)[keep]
        # The same pseudo-inverse as the call's denominator (correct): with a
        # float32 basis a covariate constant on the subset leaves a noise
        # direction that lstsq would keep and the rank cut drops, and the two
        # projections must agree (63% error in t when they did not).
        coef = self.inverse[pattern] @ (z.T @ y)
        out = np.zeros(self.n_samples, dtype=np.float64)
        resid = y - z @ coef
        scale = resid.std()
        out[keep] = resid / (scale if scale > 0 else 1.0)
        return out

    def residualize_block(self, values, traits):
        """residualize for several of the plan's traits at once: values (samples, len(traits)), raw.

        Two products of the block's size and one small solve per trait, where
        residualize column by column copied the basis rows of each trait's
        subset (125 ms a trait at 33,417 samples and rank 29).
        """
        values = np.asarray(values, dtype=np.float64)
        traits = np.asarray(traits, dtype=np.int64)
        patterns = self.pattern_of_trait[np.searchsorted(self.traits, traits)]
        observed = np.ones(values.shape, dtype=bool)
        sizes = self.observed_counts[patterns]
        if traits.size:
            observed[np.concatenate([self.rows[p] for p in patterns]),
                     np.repeat(np.arange(traits.size), self.n_samples - sizes)] = False
        y = np.where(observed, values, 0.0)
        coef = np.einsum('jrs,sj->rj', self.inverse[patterns], self.basis.T @ y)
        resid = (y - self.basis @ coef) * observed
        deviation = (resid - resid.sum(axis=0) / sizes) * observed
        scale = np.sqrt((deviation * deviation).sum(axis=0) / sizes)
        return resid / np.where(scale > 0, scale, 1.0)

    # -- the device side ----------------------------------------------------

    def _state(self, device):
        import torch
        state = self._device.get(device)
        if state is not None:
            return state
        rank = self.basis.shape[1]
        # Whole patterns in order of their missing count, grouped while a
        # group gathers at most BLOCK_CELLS cells per variant and holds at
        # most PATTERN_BLOCK patterns; a larger pattern is its own group, its
        # rows split into segments that accumulate.
        groups, current, width = [], [], 0
        for pattern in np.argsort([len(r) for r in self.rows], kind='stable'):
            size = max(1, len(self.rows[pattern]))
            if size > BLOCK_CELLS:
                if current:
                    groups.append(current)
                    current, width = [], 0
                groups.append([int(pattern)])
                continue
            wider = max(width, size)
            if current and ((len(current) + 1) * wider > BLOCK_CELLS or len(current) >= PATTERN_BLOCK):
                groups.append(current)
                current, wider = [], size
            current.append(int(pattern))
            width = wider
        if current:
            groups.append(current)
        device_groups = []
        for group in groups:
            segments = []
            longest = max(len(self.rows[p]) for p in group)
            for first in range(0, max(1, longest), BLOCK_CELLS):
                parts = [self.rows[p][first:first + BLOCK_CELLS] for p in group]
                width = max(1, max(len(part) for part in parts))
                rows = np.zeros((len(group), width), dtype=np.int64)
                valid = np.zeros((len(group), width), dtype=np.float32)
                # FP64, cast to the scan's dtype at use (float32 would cap an
                # FP64 scan at ~1e-9 relative).
                zm = np.zeros((len(group), width, rank), dtype=np.float64)
                for slot, part in enumerate(parts):
                    rows[slot, :len(part)] = part if self._positions is None else self._positions[part]
                    valid[slot, :len(part)] = 1.0
                    zm[slot, :len(part)] = self.basis[part]
                segments.append(dict(rows=torch.as_tensor(rows, device=device),
                                     valid=torch.as_tensor(valid, device=device),
                                     padded=bool((valid == 0).any()),
                                     zm=torch.as_tensor(zm, device=device)))
            device_groups.append(dict(
                patterns=torch.as_tensor(group, device=device), segments=segments,
                inverse=torch.as_tensor(np.stack([self.inverse[p] for p in group]), dtype=torch.float64,
                                        device=device),
                subset_rank=torch.as_tensor(self.subset_rank[group], dtype=torch.float64, device=device),
                observed=torch.as_tensor(self.observed_counts[group], dtype=torch.float64, device=device)))
        state = dict(groups=device_groups, traits=torch.as_tensor(self.traits, device=device),
                     pattern_of_trait=torch.as_tensor(self.pattern_of_trait, device=device))
        self._device[device] = state
        return state

    def correct(self, centered, sum_squares, projections, products, phenotype_ss, variant_df,
                beta, t, *, calls_observed=None, missing_products=None):
        """Replace the missing traits' beta and t with complete-case OLS; return the (chunk, K) pair df.

        centered (chunk, n): the scan's centred calls (missing calls zero);
        sum_squares (chunk,): their sum of squares; projections (chunk, r):
        Z^T g for Z = complete_case_basis (the scan's covariate products);
        products (chunk, K): g^T y; phenotype_ss (K,); variant_df (chunk,),
        the variant's observed calls less rank and genotype.
        calls_observed(rows) -> (chunk, len(rows)) bool, the calls observed at
        those samples, or None when every call is observed.
        beta and t (chunk, K) are updated in place. missing_products
        (chunk, J), the products of the plan's traits alone, may stand in for
        products.
        """
        import torch
        device = centered.device
        state = self._state(device)
        chunk = centered.shape[0]
        rank = self.basis.shape[1]
        patterns = len(self.rows)
        rss = torch.empty((chunk, patterns), dtype=torch.float64, device=device)
        df = torch.empty((chunk, patterns), dtype=torch.float64, device=device)
        sum_squares = sum_squares.double()
        projections = projections.double()
        # Calls missing anywhere in the variant: n less its observed calls.
        missing_total = (self.n_samples - variant_df.double() - rank - 1)[:, None]
        for group in state['groups']:
            count = group['patterns'].numel()
            s2 = torch.zeros((chunk, count), dtype=torch.float64, device=device)
            u = torch.zeros((chunk, count, rank), dtype=torch.float64, device=device)
            absent = torch.zeros((chunk, count), dtype=torch.float64, device=device)
            for segment in group['segments']:
                width = segment['rows'].shape[1]
                flat = segment['rows'].reshape(-1)
                gm = centered.index_select(1, flat).view(chunk, count, width)
                if segment['padded']:
                    # Padding gathers row 0; its basis rows are zero already.
                    gm = gm * segment['valid']
                # One pass each (a profile of the chunk-wide version: the
                # elementwise passes, not the gather, took 60% of 3.3 ms).
                s2 += torch.linalg.vecdot(gm, gm).double()
                # A strided batch view: bmm reads (p, c, w) without a copy.
                u += torch.bmm(gm.transpose(0, 1), _cast(segment, 'zm', gm.dtype)).transpose(0, 1).double()
                if calls_observed is not None:
                    lost = ~calls_observed(flat).view(chunk, count, width)
                    if segment['padded']:
                        lost &= segment['valid'].bool()
                    absent += lost.sum(-1).double()
            u = projections[:, None, :] - u
            quadratic = torch.einsum('cpr,prs,cps->cp', u, group['inverse'], u)
            samples = group['observed'][None, :].expand(chunk, count)
            if calls_observed is not None:
                # Observed calls in S: the subset less the calls missing there.
                samples = samples - (missing_total - absent)
            rss.index_copy_(1, group['patterns'], sum_squares[:, None] - s2 - quadratic)
            df.index_copy_(1, group['patterns'], samples - group['subset_rank'][None, :] - 1.0)
        which = state['pattern_of_trait']
        rss, df = rss[:, which], df[:, which]
        traits = state['traits']
        gy = (products.index_select(1, traits) if missing_products is None else missing_products).double()
        ss = phenotype_ss.double().index_select(0, traits)[None, :]
        valid = (rss > 1e-10 * sum_squares[:, None].clamp_min(1e-300)) & (df > 0)
        safe_rss = torch.where(valid, rss, torch.ones_like(rss))
        safe_df = torch.where(valid, df, torch.ones_like(df))
        kept_beta = gy / safe_rss
        residual = (ss - gy * gy / safe_rss).clamp_min(1e-12)
        kept_t = kept_beta / torch.sqrt(residual / safe_df / safe_rss)
        nan = torch.tensor(float('nan'), dtype=torch.float64, device=device)
        beta.index_copy_(1, traits, torch.where(valid, kept_beta, nan).to(beta.dtype))
        t.index_copy_(1, traits, torch.where(valid, kept_t, nan).to(t.dtype))
        pair_df = variant_df.to(torch.float32)[:, None].expand(chunk, self.n_traits).clone()
        pair_df.index_copy_(1, traits, df.to(torch.float32))
        return pair_df


def device_significant_pairs_by_pair_df(beta, t, status, pair_df, critical, *, start=0, threshold=None):
    """reduce.device_significant_pairs with a df per (variant, trait) pair, as missing phenotypes give.

    pair_df (rows, traits) float32; critical is reduce.device_significance_critical's
    FP32 lookup by whole-number df. A pair passes when t is finite, its p is
    at most `threshold`, df > 0 and its variant is valid (status 0).

    - Whole-number df (complete-case OLS, CompleteCasePlan.correct): |t| >=
      critical[df], exactly.
    - Fractional df ('impute', ImputedPlan: variant_df x trait_df / df) sits
      between two whole ones, and the critical |t| falls as df grows, so
      |t| >= critical[floor] passes and |t| < critical[ceil] fails for certain.
      Only the pairs in between take the exact tail on the device
      (tails.neg_log10_p_device, FP64) against -log10 threshold.

    Yields the same seven-field chunks, each pair with its own df, so only
    passing pairs leave the device -- not the dense beta, t and df the host
    selector needs (11 GB per 2,048 x 446,000 chunk). Kept out of reduce.py,
    whose source the device-selection census hashes.
    """
    import math

    import torch

    from .host_significance import predicate_block_shape
    from .reduce import _owned_pairs
    from .selection_geometry import DEVICE_SELECTION_MAX_CELLS, device_selection_shape
    from .tails import neg_log10_p_device
    if (t.dtype != torch.float32 or beta.dtype != torch.float32 or pair_df.dtype != torch.float32
            or critical.dtype != torch.float32):
        raise ValueError('Device significance requires native FP32 statistics and df')
    rows, traits = t.shape
    if traits < 1 or beta.shape != t.shape or pair_df.shape != t.shape or status.shape != (rows,):
        raise ValueError('Invalid significance input shapes')
    valid_row = status == 0
    # One sync per chunk: complete-case df are whole numbers, and then no block looks further.
    any_fractional = bool(((pair_df != pair_df.floor()) & torch.isfinite(pair_df)).any())
    width, height, _ = device_selection_shape(rows, traits, DEVICE_SELECTION_MAX_CELLS)
    for first in range(0, rows, height):
        last = min(rows, first + height)
        for left in range(0, traits, width):
            right = min(traits, left + width)
            keep = t.new_empty((last - first, right - left), dtype=torch.bool)
            predicate_height, predicate_width, _ = predicate_block_shape(last - first, right - left)
            for top in range(0, last - first, predicate_height):
                bottom = min(last - first, top + predicate_height)
                for column in range(0, right - left, predicate_width):
                    stop = min(right - left, column + predicate_width)
                    df = pair_df[first + top:first + bottom, left + column:left + stop]
                    usable = valid_row[first + top:first + bottom, None] & (df > 0) & torch.isfinite(df)
                    whole = df.floor()
                    sure = torch.where(usable, critical[whole.to(torch.int64).clamp(0, len(critical) - 1)], torch.inf)
                    magnitude = t[first + top:first + bottom, left + column:left + stop].abs()
                    mask = keep[top:bottom, column:stop]
                    torch.ge(magnitude, sure, out=mask)
                    mask &= magnitude < torch.inf
                    fractional = usable & (df != whole) if any_fractional else None
                    if fractional is not None and bool(fractional.any()):
                        floor = critical[df.ceil().to(torch.int64).clamp(0, len(critical) - 1)]
                        between = fractional & (magnitude >= floor) & (magnitude < sure)
                        if bool(between.any()):
                            if threshold is None:
                                raise ValueError('fractional pair df need the p threshold')
                            where = between.nonzero(as_tuple=True)
                            exact = torch.empty((where[0].numel(), 1), dtype=torch.float64, device=t.device)
                            neg_log10_p_device(magnitude[where].reshape(-1, 1), df[where].reshape(-1, 1), out=exact)
                            mask[where] = exact.reshape(-1) >= -math.log10(threshold)
            ri, ti = keep.nonzero().unbind(1)
            if not ri.numel():
                yield (start + first, start + last, np.empty(0, np.int64), np.empty(0, np.int64),
                       np.empty(0, np.float32), np.empty(0, np.float32), np.empty(0, np.float32))
                continue
            packed = torch.stack((ri.to(torch.int32), ti.to(torch.int32),
                                  beta[first:last, left:right][ri, ti].view(torch.int32),
                                  t[first:last, left:right][ri, ti].view(torch.int32),
                                  pair_df[first:last, left:right][ri, ti].view(torch.int32))).cpu().numpy()
            yield (start + first, start + last, *_owned_pairs(packed, start + first, left))


class ImputedPlan:
    """The release convention for missing phenotypes, with CompleteCasePlan's interface.

    The panel is mean-imputed (preprocess, missing='impute'); each pair's t is
    the imputed panel's t times sqrt(trait_df / df) and its df is variant_df
    x trait_df / df. Nothing is gathered, so correct() needs no calls.
    """

    needs_calls = False

    def __init__(self, observed_counts, covariate_rank, n_samples):
        counts = np.asarray(observed_counts, dtype=np.float64)
        df = float(n_samples - covariate_rank - 2)
        trait_df = counts - covariate_rank - 2
        self.n_samples, self.n_traits = int(n_samples), int(counts.size)
        self.traits = np.flatnonzero(counts < n_samples).astype(np.int64)
        self.scale = np.sqrt(trait_df / df)
        self.factor = trait_df / df
        self._device = {}

    def correct(self, centered, sum_squares, projections, products, phenotype_ss, variant_df,
                beta, t, **_):
        """Rescale t in place and return the (chunk, K) pair df."""
        import torch
        device = t.device
        state = self._device.get(device)
        if state is None:
            # FP64, cast at use: a float32 scale caps an FP64 scan at ~3e-8.
            state = self._device[device] = (torch.as_tensor(self.scale, dtype=torch.float64, device=device),
                                            torch.as_tensor(self.factor, dtype=torch.float64, device=device))
        scale, factor = state
        t.mul_(scale[None, :].to(t.dtype))
        return (variant_df.double().reshape(-1, 1) * factor[None, :]).float()


def _cast(segment, key, dtype):
    """segment[key] at `dtype`, cast once and kept (scans run one dtype)."""
    cached = segment.get((key, dtype))
    if cached is None:
        cached = segment[key].to(dtype)
        segment[(key, dtype)] = cached
    return cached


def call_mask(raw, *, encoding=None, missing_value=None):
    """calls_observed(rows) for CompleteCasePlan.correct, read from the scan's raw calls.

    raw (chunk, n) float with NaN for missing (the torch backend's dosage), or
    the native backend's input: int8 (missing -9 unless given), uint8/float32
    with an optional sentinel, or packed 2-bit rows (pgen_2bit: code 3 missing;
    plink_2bit: code 1 missing).
    """
    import torch
    if encoding in ('pgen_2bit', 'plink_2bit'):
        missing_code = 3 if encoding == 'pgen_2bit' else 1

        def observed(rows):
            byte = raw.index_select(1, rows // 4).to(torch.int32)
            return ((byte >> (2 * (rows % 4)).to(torch.int32)) & 3) != missing_code
        return observed
    if raw.dtype == torch.int8 and missing_value is None:
        missing_value = -9

    def observed(rows):
        values = raw.index_select(1, rows)
        present = ~torch.isnan(values) if values.is_floating_point() else torch.ones_like(values, dtype=torch.bool)
        if missing_value is not None:
            present &= values != missing_value
        return present
    return observed
