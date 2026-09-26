"""JAGWAS: score statistic, FP64 correlation, rank cutoff from a rounding target, block-triangular projection.

One joint statistic per variant, T = z' R^-1 z over the kept traits, chi-square
on r (the kept count) degrees of freedom. reduce.JagwasReduction is this class.

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

Cutoff. With R exact, what a small pivot amplifies is only the rounding of z.
To first order, under the null the rms rounding error of T over a trait set S
is 2 eps_z sqrt(tr R_S^-1) (tr R_S^-1 = sum of the VIFs), with eps_z = u
sqrt(N) the random-rounding bound on each z (u the scan dtype's unit
roundoff; measured on H100, A100 and 2080 Ti at 0.1-0.2 of it). The kept set
is the longest prefix of the greedy pivoted-Cholesky order whose error stays
within a target, TORCHGWAS_JAGWAS_T_ROUNDING (default 0.01, about 0.002 in
-log10 p), and never beyond R's FP64 numerical rank (a pivot at or below
K eps64 max diag R is zero: with FP64 statistics an exactly collinear trait
would otherwise pass the rounding target and add a degree of freedom that
carries no chi-square). tr R_S^-1 does not depend on the order, so a panel whose full set
meets the target keeps its unpivoted factor: the check is ||L^-1||_F, one
reduction and one read of two device scalars (Cholesky info, the norm).
Otherwise host LAPACK dpstrf (its FP64 numerical-rank default) orders the
traits; the leading block of L^-1 is the inverse of the leading block of L, so
every prefix's trace is a cumulative row sum. On the 22 real 35k-sample panels
it keeps every trait of the 17 well-conditioned ones and 76-99 of the five
near-singular ones, each at an estimated error just under 0.01; dropping traits
there carried no measurable signal (jagwas_rounding_cutoff_panels_20260926.py).
Variant shards share one JagwasRankSelection, so every device keeps the same
traits.

Projection. L^-1 is lower triangular, so (L^-1 z')[b] = L^-1[b, :end_b] z'[:end_b]:
one GEMM per row block, about (B+1)/(2B) of the dense work.
"""
from __future__ import annotations

import math
import threading
import warnings

import numpy as np
import torch

from .jagwas_blocks import (gram_rows, projection_flops_per_variant,  # noqa: F401
                            projection_gemm_dimensions, rounding_target_setting, triangular_blocks)

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
    """The joint test over the traits the rounding target keeps."""

    mode = "jagwas"
    width = 1

    def __init__(self, selection=None):
        self._inverse_cholesky = None
        self._n_traits = None
        self._selection = JagwasRankSelection() if selection is None else selection
        self._kept = None
        self._blocks = None

    def spawn(self):
        """A reduction for another device that shares this run's kept-trait decision."""
        return type(self)(self._selection)

    def resolved_width(self, n_traits: int) -> int:
        return 1

    def prepare(self, phenotype, device=None):
        """Factorise R once before the scan, keeping the traits the rounding target allows.

        The fast path is an FP64 Gram in sample blocks, FP64 Cholesky and
        triangular solve, and the check. The phenotype, correlation, factor,
        identity and inverse coexist during the solve; during the Gram the
        phenotype, correlation and one FP64 block (jagwas_blocks.gram_rows).
        """
        matrix = torch.as_tensor(phenotype)
        if device is not None:
            matrix = matrix.to(device)
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
        factor, info = torch.linalg.cholesky_ex(correlation)
        identity = torch.eye(traits, dtype=torch.float64, device=correlation.device)
        inverse = torch.linalg.solve_triangular(factor, identity, upper=False)
        # ||L^-1||_F^2 = tr R^-1, read with the Cholesky info in one transfer.
        check = torch.stack((info.double(), torch.linalg.vector_norm(inverse)))
        if check.device.type == 'meta':
            return self._set_factor(inverse)
        failed, norm = check.tolist()
        precision, limit = self._selection.trace_limit(samples, traits, matrix.dtype)
        trace = norm * norm
        if failed == 0 and math.isfinite(trace) and trace <= limit:
            record = dict(method='cholesky', traits=traits, rank=traits, kept=None, target=self._selection.target,
                          precision=precision, rounding_error=2.0 * precision * math.sqrt(trace), dropped=[])
        else:
            del factor, identity, inverse
            inverse = None
            record = self._pivoted_selection(correlation, limit, precision)
        agreed = self._selection.agree(record)
        if agreed is record and record['dropped']:
            warnings.warn(f"jagwas: {len(record['dropped'])} of {traits} traits are collinear with the others and "
                          f"were dropped (keeping them would put the rounding error of T above "
                          f"{self._selection.target:g}); the joint test has {record['rank']} degrees of freedom",
                          stacklevel=2)
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
        # dpstrf completes columns 0..numerical_rank-1 of L, so the Schur
        # complement diagonal of a dropped trait is R_jj - ||L[j, :rank]||^2.
        residual = diagonal[order[rank:]] - np.square(lower[rank:, :rank]).sum(axis=1)
        dropped = sorted((int(index), float(value / diagonal[index])) for index, value in zip(order[rank:], residual))
        return dict(method='pivoted_cholesky', traits=traits, rank=rank,
                    kept=sorted(int(index) for index in order[:rank]), target=self._selection.target,
                    precision=precision, rounding_error=2.0 * precision * math.sqrt(prefix[rank - 1]),
                    # Exactly collinear traits have residual 0 up to rounding (either sign): VIF None.
                    dropped=[dict(index=index, residual_variance=value,
                                  vif=1.0 / value if value > 0 else None) for index, value in dropped])

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
        self._blocks = [(start, end, inverse[start:end, :end]) for start, end in triangular_blocks(self._n_traits)]
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
        return dict(method=record['method'], traits=record['traits'], rank=record['rank'],
                    target=record['target'], precision=record['precision'],
                    rounding_error=record['rounding_error'], dropped=dropped)

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
        # (L^-1 z')[start:end] = L^-1[start:end, :end] z'[:end], written into
        # contiguous row blocks of one K x chunk buffer: one GEMM per block,
        # then one square and one column sum (few launches per block).
        transposed = scores.T
        projected = torch.empty((self._n_traits, scores.shape[0]), dtype=scores.dtype, device=scores.device)
        for start, end, rows in self._blocks:
            torch.mm(rows, transposed[:end], out=projected[start:end])
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
