"""JAGWAS: score statistic, FP64 correlation, eigen truncation (or a rank cutoff), block-triangular projection.

One joint statistic per variant, T = z' R^+ z over the kept eigen-directions
(or traits), chi-square on r (the kept count) degrees of freedom. reduce.JagwasReduction is this class.

Statistic. z_j = sqrt(df) r_j, computed from the scan's t as z = t / sqrt(1 +
t^2 / df) (r = t / sqrt(df + t^2)). z is linear in the phenotype, so a trait
that is a linear combination of others (r_c = a r_1 + b r_2) carries no extra
T. t is not (t = r sqrt(df / (1 - r^2))): for a strong effect the collinear
direction gets v't ~ z^3 / (2N), which its small pivot amplifies. Injected into
near-singular real panels, a t-based T overstated hits by up to 69 at z=10 and
2,300 at z=20 whatever the cutoff; the score form stays within 0.017 of its
exact value at every effect size, and matches t in the null to O(r^2)
(benchmarks/jagwas_nonlinearity_20260926.py).

R. The statistics come from the scanned (residualised, FP32) panel, so their
null correlation is exactly that panel's Gram matrix, formed in FP64 in sample
blocks from the FP32 values (2 N K^2 FP64 FLOPs once per device). An FP32
Gram rounds each entry by ~2.4e-7, which on near-singular real panels (cond
5e7..1e9) moved T by up to 37% (benchmarks/jagwas_rank_error_20260925.py).

Cutoff (default): eigen truncation. With R = U Lambda U', T keeps the
eigen-directions with eigenvalue above rcond x the largest (numpy pinv's rule;
TORCHGWAS_JAGWAS_RCOND, default 1e-3, the reference JAGWAS pipeline's value),
and never more than the rounding target below allows (the sum of 1/lambda over
the kept directions plays tr R_S^-1). df is the kept count. Why not drop
traits: on collinear imaging panels a one-ulp FP32 perturbation of the
phenotype (what another device or kernel rounds differently) changed which
traits greedy pivoting kept, df by 1-2 and T by up to 7% at the rounding
cutoff and 12% at VIF <= 100, because those panels hold near-exact ties; the
eigen-directions above 1e-3 kept the same df and T within 1e-7. The directions
below it were also where a few outlier samples made the panels heavy-tailed.

Rounding cutoff (rcond=0 or TORCHGWAS_JAGWAS_RCOND=0). With R exact, what a
small pivot amplifies is only the rounding of z. To first order, under the
null the rms rounding error of T over a trait set S is 2 eps_z sqrt(tr R_S^-1)
(tr R_S^-1 = sum of the VIFs), with eps_z = u sqrt(N) the random-rounding bound
on each z (u the scan dtype's unit roundoff; measured on H100, A100 and 2080 Ti
at 0.1-0.2 of it). The kept set is the longest prefix of the greedy
pivoted-Cholesky order whose error stays within a target,
TORCHGWAS_JAGWAS_T_ROUNDING (default 0.01, about 0.002 in -log10 p), and never
beyond R's FP64 numerical rank (a pivot at or below K eps64 max diag R is zero:
with FP64 statistics an exactly collinear trait would otherwise pass the
rounding target and add a degree of freedom that carries no chi-square).
tr R_S^-1 does not depend on the order, so a panel whose full set meets the
target keeps its unpivoted factor: the check is ||L^-1||_F, one reduction and
one read of two device scalars (Cholesky info, the norm). Otherwise host LAPACK
dpstrf (its FP64 numerical-rank default) orders the traits; the leading block
of L^-1 is the inverse of the leading block of L, so every prefix's trace is a
cumulative row sum. min_residual (TORCHGWAS_JAGWAS_MIN_RESIDUAL) also stops the
prefix at the first trait with less than that share of its variance its own.
Variant shards share one JagwasRankSelection, so every device keeps the same
count (and, here, the same traits).

Projection. Both factors are cut into row blocks, one GEMM each: the eigen
R (from Lambda_k^-1/2 U_k' = Q R) is upper trapezoidal, block [s, e) using
columns [s, K); L^-1 is lower triangular, block [s, e) using [0, e). Either is
about (B+1)/(2B) of the dense work.

Groups. JagwasGroups runs one such test per trait group from a single scan of
the groups' concatenated panel, so G panels on the same samples and
covariates cost one genotype pass instead of G.
"""
from __future__ import annotations

import math
import threading
import warnings

import numpy as np
import torch

from .jagwas_blocks import (checked_min_residual, checked_rcond, cutoff_settings, gram_rows,  # noqa: F401
                            projection_flops_per_variant, projection_gemm_dimensions,
                            rounding_target_setting, triangular_blocks)

_CUDA_LINALG_LOCK = threading.Lock()
_CUDA_LINALG_READY = False


def _ensure_cuda_linalg_initialized(device):
    """Load CUDA linalg once before independent factors enter from threads.

    PyTorch's lazy CUDA linalg wrapper can race on concurrent first use (2.5.1
    LinearAlgebraStubs.cpp). A tiny Cholesky loads its dispatch library under a
    process-wide lock; the actual phenotype factors remain concurrent. CPU and
    meta execution do not initialize CUDA, and failures never mark it ready.
    """
    global _CUDA_LINALG_READY
    if device.type != 'cuda' or _CUDA_LINALG_READY:
        return
    with _CUDA_LINALG_LOCK:
        if not _CUDA_LINALG_READY:
            with torch.cuda.device(device):
                identity = torch.eye(2, dtype=torch.float64, device=device)
                torch.linalg.cholesky(identity)
            _CUDA_LINALG_READY = True


class JagwasRankSelection:
    """The kept trait set: decided by the first prepared factor, adopted by every other device."""

    def __init__(self, target=None):
        self.target = rounding_target_setting() if target is None else float(target)
        if not self.target > 0:
            raise ValueError('the jagwas rounding target must be positive')
        self._lock = threading.Lock()
        self.record = None

    def trace_limit(self, samples, traits, dtype):
        """(precision of z, largest tr R_S^-1 the kept set may have).

        The rounding target bounds the null rms rounding error of T. R's own
        FP64 numerical rank bounds it too: a pivot at or below K eps64 max
        diag R (the LAPACK dpstrf default) is zero, and an exactly collinear
        trait there adds a degree of freedom but no chi-square. Every greedy
        pivot is at least 1/tr R_S^-1, so tr <= 1/(K eps64) keeps them all
        above it. With FP32 statistics the rounding target is the far smaller
        limit (~2e5 against ~1e13 at K=100, N=35k); with FP64 it is this one.
        """
        precision = 0.5 * torch.finfo(dtype).eps * math.sqrt(samples)
        numerical_rank = 1.0 / (traits * torch.finfo(torch.float64).eps)
        return precision, min((self.target / (2.0 * precision)) ** 2, numerical_rank)

    def agree(self, record):
        with self._lock:
            if self.record is None:
                self.record = record
            return self.record


class JagwasReduction:
    """The joint test over the eigen-directions (default) or traits the cutoff keeps."""

    mode = "jagwas"
    width = 1

    def __init__(self, selection=None, name=None, rcond=None, min_residual=None):
        self._inverse_cholesky = None
        self._n_traits = None
        self._selection = JagwasRankSelection() if selection is None else selection
        self._kept = None
        self._blocks = None
        self.name = name
        # rcond (default 1e-3) drops eigen-directions (_eigen_factor); rcond=0
        # selects the rounding cutoff over traits, which min_residual extends
        # (drop a trait while less than that fraction of its variance is its own).
        if rcond not in (None, 0) and min_residual is not None:
            raise ValueError("pass rcond (drop eigen-directions) or min_residual (drop traits), not both")
        if rcond == 0:
            self.rcond = None
            self.min_residual = None if min_residual is None else checked_min_residual(min_residual)
        elif rcond is not None:
            self.rcond, self.min_residual = checked_rcond(rcond), None
        elif min_residual is not None:
            self.rcond, self.min_residual = None, checked_min_residual(min_residual)
        else:
            self.rcond, self.min_residual = cutoff_settings()

    def spawn(self):
        """A reduction for another device that shares this run's kept-trait decision."""
        return type(self)(self._selection, self.name, 0 if self.rcond is None else self.rcond, self.min_residual)

    def resolved_width(self, n_traits: int) -> int:
        return 1

    def prepare(self, phenotype, device=None, columns=None):
        """Factorise R once before the scan, keeping the traits the rounding target allows.

        The fast path is an FP64 Gram in sample blocks, FP64 Cholesky and
        triangular solve, and the check. The phenotype, correlation, factor,
        identity and inverse coexist during the solve; during the Gram the
        phenotype, correlation and one FP64 block (jagwas_blocks.gram_rows).

        columns: this test's traits within a wider scanned panel (JagwasGroups).
        reduce() then takes the whole panel's t and gathers its kept traits.
        """
        matrix = torch.as_tensor(phenotype)
        if device is not None:
            matrix = matrix.to(device)
        self._kept = None
        if columns is None:
            return self._prepare(matrix)
        index = torch.as_tensor(np.asarray(columns, dtype=np.int64), device=matrix.device)
        self._prepare(matrix.index_select(1, index))
        self._kept = index if self._kept is None else index.index_select(0, self._kept)
        return self

    def _prepare(self, matrix):
        samples, traits = matrix.shape
        if traits < 1:
            raise ValueError("jagwas needs at least one trait")
        _ensure_cuda_linalg_initialized(matrix.device)
        correlation = torch.zeros((traits, traits), dtype=torch.float64, device=matrix.device)
        rows = gram_rows(samples, traits)
        for start in range(0, samples, rows):
            block = matrix[start:start + rows].double()
            correlation.addmm_(block.T, block)
        del block
        correlation /= float(samples)
        if self.rcond is not None:
            return self._eigen_factor(correlation, samples, traits, matrix.dtype)
        factor, info = torch.linalg.cholesky_ex(correlation)
        identity = torch.eye(traits, dtype=torch.float64, device=correlation.device)
        inverse = torch.linalg.solve_triangular(factor, identity, upper=False)
        # ||L^-1||_F^2 = tr R^-1, read with the Cholesky info in one transfer.
        # With min_residual, also the largest VIF given all other traits,
        # R_jj (R^-1)_jj = R_jj ||L^-1[:, j]||^2: a trait's share of its own
        # variance given any subset is at least 1 / that.
        checks = [info.double(), torch.linalg.vector_norm(inverse)]
        if self.min_residual is not None:
            checks.append((inverse.square().sum(dim=0) * torch.diagonal(correlation)).max())
        check = torch.stack(checks)
        if check.device.type == 'meta':
            return self._set_factor(inverse)
        failed, norm, *largest_vif = check.tolist()
        precision, limit = self._selection.trace_limit(samples, traits, matrix.dtype)
        trace = norm * norm
        if (failed == 0 and math.isfinite(trace) and trace <= limit
                and (not largest_vif or largest_vif[0] * self.min_residual <= 1.0)):
            record = dict(method='cholesky', traits=traits, rank=traits, kept=None, target=self._selection.target,
                          precision=precision, rounding_error=2.0 * precision * math.sqrt(trace), dropped=[])
            if self.min_residual is not None:
                record['min_residual'] = self.min_residual
        else:
            del factor, identity, inverse
            inverse = None
            record = self._pivoted_selection(correlation, limit, precision)
        agreed = self._selection.agree(record)
        if agreed is record and record['dropped']:
            label = "jagwas" if self.name is None else f"jagwas group {self.name}"
            reason = (f"less than {self.min_residual:g} of each one's variance is its own"
                      if record.get('bound') == 'min_residual' else
                      f"keeping them would put the rounding error of T above {self._selection.target:g}")
            warnings.warn(f"{label}: {len(record['dropped'])} of {traits} traits are collinear with the others and "
                          f"were dropped ({reason}); the joint test has {record['rank']} degrees of freedom",
                          stacklevel=3)
        if agreed['kept'] is None and inverse is not None:
            return self._set_factor(inverse)
        return self._subset_factor(correlation, agreed['kept'])

    def _pivoted_selection(self, correlation, limit, precision):
        """The longest greedy prefix within the trace limit, and each dropped trait's residual variance."""
        from scipy.linalg import lapack, solve_triangular
        matrix = correlation.double().cpu().numpy()
        if not np.isfinite(matrix).all():
            raise ValueError("the trait correlation matrix is not finite, so the jagwas quadratic form is undefined")
        traits = matrix.shape[0]
        diagonal = np.diag(matrix).copy()
        lower, pivots, numerical_rank, info = lapack.dpstrf(matrix, tol=-1.0, lower=1)
        if info < 0 or numerical_rank < 1:
            raise ValueError("no jagwas trait has nonzero variance")
        lower = np.tril(lower)
        inverse = solve_triangular(lower[:numerical_rank, :numerical_rank], np.eye(numerical_rank), lower=True)
        prefix = np.cumsum(np.square(inverse).sum(axis=1))
        rank = max(1, int(np.searchsorted(prefix, limit, side='right')))
        order = pivots - 1
        bound = 'rounding'
        if self.min_residual is not None:
            # Pivot i is trait order[i]'s variance left after the traits before
            # it; the greedy order stops at the first below the threshold.
            share = np.square(np.diag(lower)[:numerical_rank]) / diagonal[order[:numerical_rank]]
            below = np.flatnonzero(share < self.min_residual)
            by_share = int(below[0]) if len(below) else numerical_rank
            if by_share < rank:
                rank, bound = max(1, by_share), 'min_residual'
        # dpstrf completes columns 0..numerical_rank-1 of L, so the Schur
        # complement diagonal of a dropped trait is R_jj - ||L[j, :rank]||^2.
        residual = diagonal[order[rank:]] - np.square(lower[rank:, :rank]).sum(axis=1)
        dropped = sorted((int(index), float(value / diagonal[index])) for index, value in zip(order[rank:], residual))
        record = dict(method='pivoted_cholesky', traits=traits, rank=rank,
                      kept=sorted(int(index) for index in order[:rank]), target=self._selection.target,
                      precision=precision, rounding_error=2.0 * precision * math.sqrt(prefix[rank - 1]),
                      # Exactly collinear traits have residual 0 up to rounding (either sign): VIF None.
                      dropped=[dict(index=index, residual_variance=value,
                                    vif=1.0 / value if value > 0 else None) for index, value in dropped])
        if self.min_residual is not None:
            record.update(min_residual=self.min_residual, bound=bound)
        return record

    def _eigen_factor(self, correlation, samples, traits, dtype):
        """Keep R's eigen-directions above rcond x the largest eigenvalue (numpy pinv's rule).

        T = ||Lambda_k^-1/2 U_k' z||^2 over the k kept directions, chi-square on
        k: the same statistic whatever basis the traits are expressed in, and a
        kept subspace that rounding cannot reorder unless an eigenvalue sits
        within rounding of the cutoff. The directions are also kept only while
        the rounding target allows (the sum of 1/lambda over them plays tr R_S^-1).

        The projection is R from Lambda_k^-1/2 U_k' = Q R: T = ||R z||^2, and R
        (k x K) is upper trapezoidal, so row block [s, e) needs only columns
        [s, K) and costs what the triangular factor does, about half the dense
        k x K product. R'R = U_k Lambda_k^-1 U_k' whatever basis eigh chose
        inside a repeated eigenvalue, so R is fixed by the kept subspace (up to
        row signs, which T does not see).
        """
        values, vectors = torch.linalg.eigh(correlation)
        if values.device.type == 'meta':
            rank = traits  # the kept count is data; a meta trace takes the full panel
        else:
            rank = self._eigen_rank(values.cpu().numpy()[::-1], samples, traits, dtype)
        # eigh is ascending: the kept directions are the last `rank` columns.
        top_values, top_vectors = values[-rank:], vectors[:, -rank:]
        scaled = (top_vectors / top_values.sqrt()).T
        del values, vectors, top_values, top_vectors
        return self._set_upper_factor(torch.linalg.qr(scaled, mode='r')[1], traits)

    def _eigen_rank(self, spectrum, samples, traits, dtype):
        """The kept count from R's descending spectrum (one K-element FP64 read from the device)."""
        largest = float(spectrum[0])
        if not largest > 0:
            raise ValueError("no jagwas trait has nonzero variance")
        precision, limit = self._selection.trace_limit(samples, traits, dtype)
        above = int((spectrum > self.rcond * largest).sum())
        positive = spectrum[spectrum > 0]
        within = int(np.searchsorted(np.cumsum(1.0 / positive), limit, side='right'))
        rank = max(1, min(above, within))
        inverse_sum = float((1.0 / spectrum[:rank]).sum())
        record = dict(method='eigen', traits=traits, rank=rank, kept=None, target=self._selection.target,
                      precision=precision, rounding_error=2.0 * precision * math.sqrt(inverse_sum), dropped=[],
                      rcond=self.rcond, bound='rcond' if above <= within else 'rounding',
                      largest_eigenvalue=largest, smallest_kept_eigenvalue=float(spectrum[rank - 1]),
                      largest_dropped_eigenvalue=float(spectrum[rank]) if rank < traits else None)
        agreed = self._selection.agree(record)
        if agreed is record and record['rank'] < traits:
            label = "jagwas" if self.name is None else f"jagwas group {self.name}"
            warnings.warn(f"{label}: keeping {record['rank']} of {traits} eigen-directions of the trait "
                          f"correlation (eigenvalue above {self.rcond:g} x the largest); the joint test has "
                          f"{record['rank']} degrees of freedom", stacklevel=5)
        return agreed['rank']

    def _set_upper_factor(self, factor, traits):
        """Adopt an upper-trapezoidal (k x K) factor: row block [s, e) multiplies columns [s, K)."""
        self._inverse_cholesky = factor
        self._n_traits = int(factor.shape[0])
        self._blocks = [(start, end, start, traits, factor[start:end, start:])
                        for start, end in triangular_blocks(self._n_traits)]
        return self

    def _subset_factor(self, correlation, kept):
        """L^-1 of R[kept, kept] in input order, on this device."""
        index = torch.arange(correlation.shape[0], device=correlation.device) if kept is None else \
            torch.as_tensor(kept, dtype=torch.long, device=correlation.device)
        block = correlation.index_select(0, index).index_select(1, index).double()
        factor, info = torch.linalg.cholesky_ex(block)
        if int(info):
            raise ValueError("the kept jagwas traits are not positive definite on this device")
        identity = torch.eye(int(index.numel()), dtype=torch.float64, device=block.device)
        self._kept = None if kept is None else index
        return self._set_factor(torch.linalg.solve_triangular(factor, identity, upper=False))

    def _set_factor(self, inverse):
        """Adopt L^-1 and its row blocks, which are views: no second K x K copy on the device."""
        self._inverse_cholesky = inverse
        self._n_traits = int(inverse.shape[0])
        self._blocks = [(start, end, 0, end, inverse[start:end, :end])
                        for start, end in triangular_blocks(self._n_traits)]
        return self

    @property
    def degrees_of_freedom(self) -> int:
        """r, the number of kept traits: the chi-square df of T."""
        if self._n_traits is not None:
            return self._n_traits
        if self._selection.record is None:
            raise ValueError("jagwas reduction was not prepared")
        return self._selection.record['rank']

    def rank_report(self, trait_names=None):
        """The run's kept-trait decision, with dropped trait names, for the manifest and summary."""
        record = self._selection.record
        if record is None:
            return None
        dropped = [dict(item, trait=str(trait_names[item['index']])) if trait_names is not None else dict(item)
                   for item in record['dropped']]
        return dict({key: value for key, value in record.items() if key != 'kept'}, dropped=dropped)

    def host_buffers(self, chunk_size: int, width: int, pin_memory: bool = True):
        return (
            torch.empty((chunk_size, 1), dtype=torch.float32, pin_memory=pin_memory),
            torch.empty((chunk_size, 1), dtype=torch.float32, pin_memory=pin_memory),
            torch.empty((chunk_size, 1), dtype=torch.int32, pin_memory=pin_memory),
            torch.empty(chunk_size, dtype=torch.uint8, pin_memory=pin_memory),
            torch.empty(chunk_size, pin_memory=pin_memory),
        )

    def reduce(self, beta, t, status, variant_df, width):
        """(chunk, K) t in, (chunk, 1) T out in the t slot; beta has no joint analogue and is NaN.

        A variant with any non-finite kept statistic has no joint statistic and
        is invalidated from the values, not only from status (backends without
        a status word pass NaN): the zeroing below keeps one NaN from poisoning
        the GEMM row, and the row is masked afterwards.
        """
        if self._inverse_cholesky is None:
            raise ValueError("jagwas reduction was not prepared")
        if variant_df is None:
            raise ValueError("the jagwas score statistic needs the per-variant df")
        scores = t if self._kept is None else t.index_select(1, self._kept)
        degenerate = ~torch.isfinite(scores).all(dim=1)
        scores = torch.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)
        # Score form z = t / sqrt(1 + t^2 / df), in the scan precision (its
        # rounding is that of t), then FP64 for the projection. Invalid variants
        # carry t = 0 and a status, so a non-positive df never reaches T.
        scores.mul_(scores.square().div_(variant_df.unsqueeze(1)).add_(1.0).rsqrt_())
        scores = scores.double()
        # (F z')[start:end] = F[start:end, first:last] z'[first:last] for the
        # factor's nonzero columns (the eigen R: [start, K); the rounding
        # cutoff's L^-1: [0, end)), written into contiguous row blocks of one
        # k x chunk buffer: one GEMM per block, then one square and one column
        # sum (few launches per block).
        transposed = scores.T
        projected = torch.empty((self._n_traits, scores.shape[0]), dtype=scores.dtype, device=scores.device)
        for start, end, first, last, rows in self._blocks:
            torch.mm(rows, transposed[first:last], out=projected[start:end])
        statistic = projected.square_().sum(dim=0).unsqueeze(1)
        statistic = statistic.masked_fill((degenerate | (status != 0)).unsqueeze(1), float("nan"))
        return (torch.full_like(statistic, float("nan"), dtype=beta.dtype),
                statistic.to(t.dtype),
                torch.zeros_like(statistic, dtype=torch.int32),
                status, variant_df)

    def merge(self, running, incoming, trait_offset: int):
        raise ValueError(
            "jagwas cannot be trait-blocked: the statistic is a quadratic form "
            "over the whole trait correlation, so a block of traits does not "
            "carry enough information to be merged with another block")


class JagwasGroups:
    """One joint test per trait group, all from one scan of the concatenated panel.

    A trait's t depends only on its own phenotype, the covariates and the
    genotype, so each group's t, kept traits, T and df are exactly those of a
    scan of that group alone; only the genotype pass is shared. Each group is
    an independent JagwasReduction (its own correlation, rank cutoff and
    factor) over its columns of the scanned panel. Groups may overlap.
    The result has one column per group, in order.

    A group is (name, columns) or (name, columns, cutoff): cutoff is a dict of
    JagwasReduction options (rcond or min_residual) or a bare rcond, and
    replaces the defaults given here for that group.
    """

    mode = "jagwas"

    def __init__(self, groups, reductions=None, rcond=None, min_residual=None):
        groups = [tuple(group) for group in (groups.items() if isinstance(groups, dict) else groups)]
        if not groups:
            raise ValueError("jagwas groups need at least one group")
        self.names = [str(group[0]) for group in groups]
        if len(set(self.names)) != len(self.names):
            raise ValueError("jagwas group names must be unique")
        self.columns = [np.asarray(group[1], dtype=np.int64).reshape(-1) for group in groups]
        for name, columns in zip(self.names, self.columns):
            if len(columns) == 0 or len(np.unique(columns)) != len(columns):
                raise ValueError(f"jagwas group {name} needs distinct trait columns")
        defaults = dict(rcond=rcond, min_residual=min_residual)
        self.cutoffs = []
        for group in groups:
            cutoff = group[2] if len(group) > 2 else None
            if cutoff is None or cutoff == {}:
                cutoff = defaults
            elif not isinstance(cutoff, dict):
                cutoff = dict(rcond=cutoff)
            unknown = set(cutoff) - set(defaults)
            if unknown:
                raise ValueError(f"jagwas group {group[0]}: unknown cutoff option(s) {sorted(unknown)}")
            self.cutoffs.append(dict(dict.fromkeys(defaults), **cutoff))
        self.reductions = ([JagwasReduction(name=name, **cutoff) for name, cutoff in zip(self.names, self.cutoffs)]
                           if reductions is None else list(reductions))
        self.width = len(groups)

    @property
    def column_groups(self):
        """Residualise each group on its own (preprocess.residualize_and_standardize):
        a near-collinear group's kept set can turn on the rounding of its panel,
        so it must see the same bits a scan of that group alone would."""
        return self.columns

    def spawn(self):
        """Groups for another device that share each group's kept-trait decision."""
        return type(self)(list(zip(self.names, self.columns, self.cutoffs)),
                          [reduction.spawn() for reduction in self.reductions])

    def resolved_width(self, n_traits: int) -> int:
        return self.width

    def prepare(self, phenotype, device=None):
        matrix = torch.as_tensor(phenotype)
        if device is not None:
            matrix = matrix.to(device)
        traits = matrix.shape[1]
        for name, columns in zip(self.names, self.columns):
            if columns.min() < 0 or columns.max() >= traits:
                raise ValueError(f"jagwas group {name} refers outside the {traits} scanned traits")
        for reduction, columns in zip(self.reductions, self.columns):
            reduction.prepare(matrix, columns=columns)
        return self

    @property
    def degrees_of_freedom(self) -> list[int]:
        return [reduction.degrees_of_freedom for reduction in self.reductions]

    def rank_report(self, trait_names=None):
        """Each group's kept-trait report, with its name, in column order."""
        reports = []
        for name, columns, reduction in zip(self.names, self.columns, self.reductions):
            report = reduction.rank_report(None if trait_names is None else [trait_names[i] for i in columns])
            if report is None:
                return None
            reports.append(dict(group=name, **report))
        return reports

    def host_buffers(self, chunk_size: int, width: int, pin_memory: bool = True):
        return (
            torch.empty((chunk_size, width), dtype=torch.float32, pin_memory=pin_memory),
            torch.empty((chunk_size, width), dtype=torch.float32, pin_memory=pin_memory),
            torch.empty((chunk_size, width), dtype=torch.int32, pin_memory=pin_memory),
            torch.empty(chunk_size, dtype=torch.uint8, pin_memory=pin_memory),
            torch.empty(chunk_size, pin_memory=pin_memory),
        )

    def reduce(self, beta, t, status, variant_df, width):
        """(chunk, K) t of the whole panel in, (chunk, groups) T out in the t slot."""
        statistic = torch.cat([reduction.reduce(beta, t, status, variant_df, 1)[1]
                               for reduction in self.reductions], dim=1)
        return (torch.full_like(statistic, float("nan"), dtype=beta.dtype),
                statistic,
                torch.zeros_like(statistic, dtype=torch.int32),
                status, variant_df)

    def merge(self, running, incoming, trait_offset: int):
        return self.reductions[0].merge(running, incoming, trait_offset)
