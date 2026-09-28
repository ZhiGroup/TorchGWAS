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
observed mean (zero), df reduced by one each.

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


def complete_case_basis(n_samples, q_matrix):
    """Z = [1/sqrt(n), Q]: the scan's covariate design, orthonormal (Q is centred)."""
    one = np.full((int(n_samples), 1), 1.0 / np.sqrt(n_samples))
    return one if q_matrix is None else np.column_stack([one, np.asarray(q_matrix, dtype=np.float64)])


class CompleteCasePlan:
    """The missing rows of each trait with missing phenotype values, grouped by pattern."""

    def __init__(self, observed, q_matrix):
        observed = np.asarray(observed, dtype=bool)
        if observed.ndim != 2:
            raise ValueError('observed must be (samples, traits)')
        n, k = observed.shape
        self.n_samples, self.n_traits = n, k
        self.basis = complete_case_basis(n, q_matrix)
        rank = self.basis.shape[1]
        self.traits = np.flatnonzero(~observed.all(axis=0)).astype(np.int64)
        index, rows = {}, []
        self.pattern_of_trait = np.empty(self.traits.size, dtype=np.int64)
        for position, trait in enumerate(self.traits):
            missing = np.flatnonzero(~observed[:, trait])
            key = missing.tobytes()
            if key not in index:
                index[key] = len(rows)
                rows.append(missing)
            self.pattern_of_trait[position] = index[key]
        self.rows = rows
        self.inverse, self.subset_rank = [], []
        # Z^T Z itself, not I: a float32 basis is orthonormal only to ~1e-7,
        # which would lift a subset's null direction above the rank cut.
        gram = self.basis.T @ self.basis
        for missing in rows:
            zm = self.basis[missing]
            values, vectors = np.linalg.eigh(gram - zm.T @ zm)
            # A covariate constant on the subset leaves Z_S rank-deficient;
            # its df then counts the rank Z_S actually has.
            keep = values > 1e-9 * max(values.max(), 1e-300)
            self.inverse.append((vectors[:, keep] / values[keep]) @ vectors[:, keep].T)
            self.subset_rank.append(int(keep.sum()))
        self.subset_rank = np.asarray(self.subset_rank, dtype=np.int64)
        self.observed_counts = n - np.array([len(missing) for missing in rows], dtype=np.int64)
        # Where each sample's column sits in the scan's calls (None: in order).
        self._positions = None
        self._device = {}

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
        observed = ~np.isnan(np.asarray(phenotype, dtype=np.float64))
        return None if observed.all() else cls(observed, q_matrix)

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
