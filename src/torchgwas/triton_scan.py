"""Fused per-chunk statistics in Triton: prepare, finish, and finish with min-p's ranking.

The default (Torch) statistics make about ten elementwise passes over the
chunk's (variants x samples) calls and several over its (variants x traits)
products: ~40 kernel launches per chunk. At K = 512 on four GPUs the scan is
host-bound, and issuing those launches is most of the per-chunk host cost
that does not parallelize across GPUs (docs/autotune_design_20260924.md). The
same fusion already exists as CUDA (native/scan_statistics.cu, the
`native_fused` backend), but it is opt-in and built only for sm80/sm90.
Triton compiles for whatever GPU is present, so these kernels can be a
default, with the Torch statistics as the fallback.

The contract is scan_gpu's, kernel for kernel:

- prepare(raw, scale, missing_value, encoding, n_samples) -> centred float32
  calls (a missing call is exactly zero), their sum of squares, the raw
  dosage min and max, and the observed count. Mean and sum of squares
  accumulate in FP64 (FP32 lane partials, FP64 across lanes).
- finish(products, centred_ss, min, max, phenotype_ss, present, df_offset)
  -> beta, t, status (0 valid, 1 non-finite residual, 2 invariant, 3 dosage
  out of [0, 2] when asked), invalid beta and t zero.
- finish_min_p: finish, but only each variant's largest |t| leaves the
  kernel (VariantReduction('min-p'): NaN and invalid variants rank last, the
  first maximum wins), so the (variants x traits) beta and t are never
  written.
"""
from __future__ import annotations

import functools
import threading

import torch

try:
    import triton
    import triton.language as tl
except ImportError:  # the Torch statistics are the fallback
    triton = None
    tl = None

# One kernel launch at a time across threads (variant shards, phenotype
# tiles). Triton 3.1 binds a compiled kernel to its C launcher lazily, and on
# four H100 shards one launch of finish_min_p once went through prepare's
# launcher ("function takes exactly 19 arguments (23 given)": 9 launch fields
# plus prepare's 10 arguments, against finish_min_p's 14). The GIL already
# serialises most of a launch; the lock also covers where Triton releases it.
_LAUNCH_LOCK = threading.Lock()

# Raw encodings: one byte or float per sample, or two-bit rows.
INT8, UINT8, FLOAT32, PGEN_2BIT, PLINK_2BIT = range(5)
# Samples per program step. A power of two, as tl.arange needs; one variant
# row per program, so a chunk launches `variants` programs per kernel.
SAMPLE_BLOCK = 1024
TRAIT_BLOCK = 1024
# prepare's warps per row. At a 1024-sample block the warp count moves no lane
# boundary, so every output is bit-identical across it; only speed changes.
# Packed rows of 22,250 samples, 4,096 per chunk, 16 warps against 4: A100
# 0.82 -> 0.56 ms, H100 0.30 -> 0.28, RTX 2080 Ti 1.35 -> 1.12. With 512 rows of
# 400,000 the A100 went the other way (1.16 -> 1.44 ms), so short chunks keep 4.
PREPARE_WARPS, PREPARE_WARPS_FEW_ROWS, PREPARE_MANY_ROWS = 16, 4, 2048


if triton is not None:
    @triton.jit
    def _decode(raw_ptr, offsets, in_row, scale, sentinel,
                KIND: tl.constexpr, HAS_SENTINEL: tl.constexpr):
        """(value, observed) for the samples at `offsets` of one raw row."""
        if KIND >= 3:
            byte = tl.load(raw_ptr + offsets // 4, mask=in_row, other=0).to(tl.int32)
            code = (byte >> ((offsets % 4) * 2)) & 3
            if KIND == 3:  # PGEN: the code is the ALT dosage, 3 missing
                observed = in_row & (code != 3)
                value = code.to(tl.float32)
            else:  # PLINK1: 01 missing, A2 dosage (code + 1) >> 1
                observed = in_row & (code != 1)
                value = ((code + 1) >> 1).to(tl.float32)
        else:
            value = tl.load(raw_ptr + offsets, mask=in_row, other=0).to(tl.float32)
            observed = in_row & (value == value)
            if HAS_SENTINEL:
                observed = observed & (value != sentinel)
            if KIND == 1:
                # Round-to-nearest: Triton's `/` is approximate, and 127 / 127.0
                # must be exactly 1, as in the CUDA and Torch paths.
                value = tl.div_rn(value, scale)
        return value, observed

    @triton.jit
    def _prepare_kernel(raw_ptr, row_stride, samples, scale, sentinel,
                        centered_ptr, ss_ptr, min_ptr, max_ptr, present_ptr,
                        KIND: tl.constexpr, HAS_SENTINEL: tl.constexpr, BLOCK: tl.constexpr):
        row = tl.program_id(0).to(tl.int64)
        raw_row = raw_ptr + row * row_stride
        lanes = tl.arange(0, BLOCK)
        total = tl.zeros([BLOCK], tl.float32)
        low = tl.full([BLOCK], float('inf'), tl.float32)
        high = tl.full([BLOCK], float('-inf'), tl.float32)
        count = tl.zeros([BLOCK], tl.int32)
        for start in range(0, samples, BLOCK):
            offsets = start + lanes
            value, observed = _decode(raw_row, offsets, offsets < samples, scale, sentinel, KIND, HAS_SENTINEL)
            total += tl.where(observed, value, 0.0)
            low = tl.minimum(low, tl.where(observed, value, float('inf')))
            high = tl.maximum(high, tl.where(observed, value, float('-inf')))
            count += observed.to(tl.int32)
        present = tl.sum(count, axis=0)
        mean = tl.where(present > 0, (tl.sum(total.to(tl.float64), axis=0)
                                      / tl.maximum(present, 1).to(tl.float64)).to(tl.float32), 0.0)
        # Nothing observed leaves the range empty, so the variant is invariant.
        tl.store(min_ptr + row, tl.where(present > 0, tl.min(low, axis=0), float('nan')))
        tl.store(max_ptr + row, tl.where(present > 0, tl.max(high, axis=0), float('nan')))
        tl.store(present_ptr + row, present)
        squares = tl.zeros([BLOCK], tl.float32)
        out_row = centered_ptr + row * samples
        for start in range(0, samples, BLOCK):
            offsets = start + lanes
            in_row = offsets < samples
            value, observed = _decode(raw_row, offsets, in_row, scale, sentinel, KIND, HAS_SENTINEL)
            # Centred, and masked: a missing call is exactly zero, so it drops
            # out of the GEMM and the sum of squares.
            centered = tl.where(observed, value - mean, 0.0)
            tl.store(out_row + offsets, centered, mask=in_row)
            squares += centered * centered
        tl.store(ss_ptr + row, tl.sum(squares.to(tl.float64), axis=0).to(tl.float32))

    @triton.jit
    def _row_terms(prod_row, ss_ptr, min_ptr, max_ptr, present_ptr, row, df_offset, traits, covariates,
                   VALIDATE: tl.constexpr, COV_BLOCK: tl.constexpr):
        """The variant's residual sum, df, validity and status (finish_kernel's rules)."""
        c = tl.arange(0, COV_BLOCK)
        gc = tl.load(prod_row + traits + c, mask=c < covariates, other=0.0)
        residual = tl.load(ss_ptr + row) - tl.sum(gc * gc, axis=0)
        low = tl.load(min_ptr + row)
        high = tl.load(max_ptr + row)
        df = tl.load(present_ptr + row).to(tl.float32) + df_offset
        valid = (residual > 1e-12) & (high > low) & (df > 0.0)
        finite = (residual == residual) & (tl.abs(residual) < float('inf'))
        status = tl.where(finite, tl.where(valid, 0, 2), 1)
        if VALIDATE:
            status = tl.where((low < 0.0) | (high > 2.0), 3, status)
        safe = tl.where(residual != residual, residual, tl.maximum(residual, 1e-12))
        return safe, df, valid, status

    @triton.jit
    def _pair_terms(prod_row, yss_ptr, offsets, in_row, safe, df):
        # Correctly rounded division and square root, as finish_kernel's.
        gy = tl.load(prod_row + offsets, mask=in_row, other=0.0)
        beta = tl.div_rn(gy, safe)
        yss = tl.load(yss_ptr + offsets, mask=in_row, other=1.0) - tl.div_rn(gy * gy, safe)
        yss = tl.where(yss != yss, yss, tl.maximum(yss, 1e-12))
        return beta, tl.div_rn(beta, tl.sqrt_rn(tl.div_rn(tl.div_rn(yss, df), safe)))

    @triton.jit
    def _finish_kernel(prod_ptr, prod_stride, ss_ptr, min_ptr, max_ptr, yss_ptr, present_ptr, df_offset,
                       traits, covariates, beta_ptr, t_ptr, status_ptr,
                       VALIDATE: tl.constexpr, COV_BLOCK: tl.constexpr, BLOCK: tl.constexpr):
        row = tl.program_id(0).to(tl.int64)
        prod_row = prod_ptr + row * prod_stride
        safe, df, valid, status = _row_terms(prod_row, ss_ptr, min_ptr, max_ptr, present_ptr, row, df_offset,
                                             traits, covariates, VALIDATE, COV_BLOCK)
        tl.store(status_ptr + row, status.to(tl.uint8))
        lanes = tl.arange(0, BLOCK)
        for start in range(0, traits, BLOCK):
            offsets = start + lanes
            in_row = offsets < traits
            beta, t = _pair_terms(prod_row, yss_ptr, offsets, in_row, safe, df)
            tl.store(beta_ptr + row * traits + offsets, tl.where(valid, beta, 0.0), mask=in_row)
            tl.store(t_ptr + row * traits + offsets, tl.where(valid, t, 0.0), mask=in_row)

    @triton.jit
    def _finish_min_p_kernel(prod_ptr, prod_stride, ss_ptr, min_ptr, max_ptr, yss_ptr, present_ptr, df_offset,
                             traits, covariates, beta_ptr, t_ptr, index_ptr, status_ptr,
                             VALIDATE: tl.constexpr, COV_BLOCK: tl.constexpr, BLOCK: tl.constexpr):
        row = tl.program_id(0).to(tl.int64)
        prod_row = prod_ptr + row * prod_stride
        safe, df, valid, status = _row_terms(prod_row, ss_ptr, min_ptr, max_ptr, present_ptr, row, df_offset,
                                             traits, covariates, VALIDATE, COV_BLOCK)
        tl.store(status_ptr + row, status.to(tl.uint8))
        lanes = tl.arange(0, BLOCK)
        best = float('-inf')
        best_index = 0
        for start in range(0, traits, BLOCK):
            offsets = start + lanes
            in_row = offsets < traits
            _, t = _pair_terms(prod_row, yss_ptr, offsets, in_row, safe, df)
            score = tl.where(in_row & (t == t) & (status == 0), tl.abs(t), float('-inf'))
            block_best = tl.max(score, axis=0)
            block_index = tl.argmax(score, axis=0) + start
            # Strictly greater: the first maximum wins, as torch.max's does.
            take = block_best > best
            best_index = tl.where(take, block_index, best_index)
            best = tl.where(take, block_best, best)
        beta, t = _pair_terms(prod_row, yss_ptr, best_index, best_index < traits, safe, df)
        tl.store(beta_ptr + row, tl.where(valid, beta, 0.0))
        tl.store(t_ptr + row, tl.where(valid, t, 0.0))
        tl.store(index_ptr + row, best_index.to(tl.int32))


# Why the kernels could not run on a device, by index: the scan falls back to
# the Torch statistics, and this is where to look when it did.
UNAVAILABLE = {}


@functools.lru_cache(maxsize=None)
def _launchable(index):
    try:
        with torch.cuda.device(index):
            raw = torch.zeros((1, 64), dtype=torch.uint8, device=f'cuda:{index}')
            prepare(raw, encoding='pgen_2bit', n_samples=4)
            torch.cuda.synchronize(index)
        return True
    except Exception as error:  # noqa: BLE001 - any failure means the Torch statistics
        UNAVAILABLE[index] = f'{type(error).__name__}: {error}'[:2000]
        return False


def available(device=None):
    """True when Triton can compile and run these kernels on `device` (probed once per device)."""
    if triton is None or not torch.cuda.is_available():
        return False
    device = torch.device('cuda', torch.cuda.current_device()) if device is None else torch.device(device)
    if device.type != 'cuda':
        return False
    return _launchable(device.index if device.index is not None else torch.cuda.current_device())


def _check(tensor, name):
    if not tensor.is_cuda or not tensor.is_contiguous():
        raise ValueError(f'{name} must be contiguous and CUDA-resident')


def prepare(raw, scale=1.0, missing_value=None, *, encoding=None, n_samples=None):
    """scan_gpu.prepare's contract: (centred, centred_ss, minimum, maximum, present)."""
    _check(raw, 'raw')
    if raw.ndim != 2 or min(raw.shape) <= 0:
        raise ValueError('unsupported raw shape')
    if encoding in ('pgen_2bit', 'plink_2bit'):
        if raw.dtype != torch.uint8 or not isinstance(n_samples, int) or isinstance(n_samples, bool) \
                or n_samples <= 0 or raw.shape[1] < (n_samples + 3) // 4:
            raise ValueError('packed two-bit input requires uint8 rows and a valid logical sample count')
        kind, samples, sentinel = (PGEN_2BIT if encoding == 'pgen_2bit' else PLINK_2BIT), n_samples, None
    elif encoding is None:
        if n_samples is not None:
            raise ValueError('n_samples is only used with packed input')
        kinds = {torch.int8: INT8, torch.uint8: UINT8, torch.float32: FLOAT32}
        if raw.dtype not in kinds:
            raise ValueError('unsupported raw dtype')
        kind, samples = kinds[raw.dtype], raw.shape[1]
        sentinel = -9 if raw.dtype == torch.int8 and missing_value is None else missing_value
    else:
        raise ValueError('unsupported preparation encoding')
    rows = raw.shape[0]
    device = raw.device
    centered = torch.empty((rows, samples), dtype=torch.float32, device=device)
    ss = torch.empty(rows, dtype=torch.float32, device=device)
    minimum, maximum = torch.empty_like(ss), torch.empty_like(ss)
    present = torch.empty(rows, dtype=torch.int32, device=device)
    with _LAUNCH_LOCK, torch.cuda.device(device):
        _prepare_kernel[(rows,)](raw, raw.stride(0), samples, float(scale),
                                 float(0 if sentinel is None else sentinel),
                                 centered, ss, minimum, maximum, present,
                                 KIND=kind, HAS_SENTINEL=sentinel is not None, BLOCK=SAMPLE_BLOCK,
                                 num_warps=PREPARE_WARPS if rows >= PREPARE_MANY_ROWS else PREPARE_WARPS_FEW_ROWS)
    return centered, ss, minimum, maximum, present


def _finish_arguments(products, traits):
    _check(products, 'products')
    covariates = products.shape[1] - traits
    if traits <= 0 or covariates < 0:
        raise ValueError('products must hold the traits then the covariates')
    return covariates, max(16, triton.next_power_of_2(max(covariates, 1))), \
        min(TRAIT_BLOCK, max(16, triton.next_power_of_2(traits)))


def finish(products, centered_ss, minimum, maximum, phenotype_ss, present, df_offset, validate_range=False):
    """scan_gpu.finish's contract: beta, t (variants x traits), uint8 status."""
    rows, traits = products.shape[0], phenotype_ss.numel()
    covariates, cov_block, block = _finish_arguments(products, traits)
    beta = torch.empty((rows, traits), dtype=torch.float32, device=products.device)
    t = torch.empty_like(beta)
    status = torch.empty(rows, dtype=torch.uint8, device=products.device)
    with _LAUNCH_LOCK, torch.cuda.device(products.device):
        _finish_kernel[(rows,)](products, products.stride(0), centered_ss, minimum, maximum, phenotype_ss,
                                present, float(df_offset), traits, covariates, beta, t, status,
                                VALIDATE=bool(validate_range), COV_BLOCK=cov_block, BLOCK=block, num_warps=4)
    return beta, t, status


def finish_min_p(products, centered_ss, minimum, maximum, phenotype_ss, present, df_offset, validate_range=False):
    """finish, keeping each variant's largest |t| only: beta, t, index (variants x 1), status."""
    rows, traits = products.shape[0], phenotype_ss.numel()
    covariates, cov_block, block = _finish_arguments(products, traits)
    beta = torch.empty((rows, 1), dtype=torch.float32, device=products.device)
    t = torch.empty_like(beta)
    index = torch.empty((rows, 1), dtype=torch.int32, device=products.device)
    status = torch.empty(rows, dtype=torch.uint8, device=products.device)
    with _LAUNCH_LOCK, torch.cuda.device(products.device):
        _finish_min_p_kernel[(rows,)](products, products.stride(0), centered_ss, minimum, maximum, phenotype_ss,
                                      present, float(df_offset), traits, covariates, beta, t, index, status,
                                      VALIDATE=bool(validate_range), COV_BLOCK=cov_block, BLOCK=block,
                                      num_warps=4)
    return beta, t, index, status


# -- the -log10 P tail ---------------------------------------------------------
#
# tails.py's device tail -- P(|T| > t) = I_x(df/2, 1/2) in logs, one Lentz
# continued fraction per cell (the direct one, or the reflected one near t =
# 0), 40 fixed iterations -- as one kernel with everything in registers. The
# AOTInductor build of the same stages exists only for the architectures it
# was built for; elsewhere the stages ran as ~12 torch.compile calls per tail,
# each holding the GIL through Dynamo's guards: on four 2080 Ti shards that was
# 40% of all GIL time (benchmarks/synthetic_multigpu_20260928.py). One launch
# here, on any GPU Triton supports. FP64 throughout, as the stages.
#
# One cell per thread (block = 32 x warps): each cell's 40 FP64 iterations
# already fill a thread's registers. A narrow tail -- min-p's 4,096 winners per
# chunk -- needs small blocks to reach every SM: the former 256-cell block ran
# them as 16 programs, 0.28 ms on an A100 against 0.08 at 64 (H100 0.058 ->
# 0.041, RTX 2080 Ti 0.41 -> 0.20). A dense chunk (4,096 x 512) has programs
# to spare and takes 8 warps: A100 2.20 -> 1.78 ms, H100 0.68 -> 0.64, 2080 Ti
# unchanged (FP64-bound at 31 ms). The tail is elementwise, so the outputs are
# bit-identical across these.
TAIL_BLOCK, TAIL_WARPS = 256, 8
TAIL_BLOCK_NARROW, TAIL_WARPS_NARROW, TAIL_WIDE_CELLS = 64, 2, 1 << 20
TAIL_ITERATIONS = 40


def _tail_constants():
    import math
    from scipy import special
    return float(special.gammaln(0.5)), 1.0 / math.log(10.0)


_LGAMMA_HALF, _INV_LN10 = _tail_constants()


if triton is not None:
    try:
        from triton.language.extra import libdevice as _libdevice
    except ImportError:  # older layouts
        from triton.language.extra.cuda import libdevice as _libdevice

    @triton.jit
    def _tail_kernel(t_ptr, df_ptr, out_ptr, cells, cols, t_row, t_col, df_row, df_col,
                     TINY: tl.constexpr, LGAMMA_HALF: tl.constexpr, INV_LN10: tl.constexpr,
                     ITERATIONS: tl.constexpr, BLOCK: tl.constexpr):
        offsets = tl.program_id(0).to(tl.int64) * BLOCK + tl.arange(0, BLOCK)
        valid = offsets < cells
        row = offsets // cols
        col = offsets % cols
        # FP64 constants built as FP64: a Python float literal would be FP32,
        # and 1e-300 would flush to zero.
        tiny = tl.full([BLOCK], TINY, tl.float64)
        t = tl.abs(tl.load(t_ptr + row * t_row + col * t_col, mask=valid, other=0.0).to(tl.float64))
        df = tl.load(df_ptr + row * df_row + col * df_col, mask=valid, other=1.0).to(tl.float64)
        a = df * 0.5
        squared = t * t
        total = df + squared
        x = df / total
        y = squared / total
        reflect = x >= (a + 1.0) / (a + 2.5)
        first = tl.where(reflect, 0.5, a)
        second = tl.where(reflect, a, 0.5)
        argument = tl.where(reflect, y, x)
        d = 1.0 - (first + second) * argument / (first + 1.0)
        d = 1.0 / tl.where(tl.abs(d) < tiny, tiny, d)
        c = tl.full([BLOCK], 1.0, tl.float64)
        h = d
        qab = first + second
        qap = first + 1.0
        qam = first - 1.0
        for index in range(ITERATIONS):
            m = (index + 1).to(tl.float64)
            m2 = 2.0 * m
            step = m * (second - m) * argument / ((qam + m2) * (first + m2))
            d = 1.0 + step * d
            d = 1.0 / tl.where(tl.abs(d) < tiny, tiny, d)
            c = 1.0 + step / c
            c = tl.where(tl.abs(c) < tiny, tiny, c)
            h = h * d * c
            step = -(first + m) * (qab + m) * argument / ((first + m2) * (qap + m2))
            d = 1.0 + step * d
            d = 1.0 / tl.where(tl.abs(d) < tiny, tiny, d)
            c = 1.0 + step / c
            c = tl.where(tl.abs(c) < tiny, tiny, c)
            h = h * d * c
        log_beta = (_libdevice.lgamma(a) + tl.full([BLOCK], LGAMMA_HALF, tl.float64)
                    - _libdevice.lgamma(a + 0.5))
        log_y = _libdevice.log(tl.maximum(y, tiny))
        log_direct = (a * _libdevice.log(tl.maximum(x, tiny)) + 0.5 * log_y
                      - _libdevice.log(a) - log_beta + _libdevice.log(h))
        upper = 2.0 * _libdevice.exp(-log_beta + 0.5 * log_y + a * _libdevice.log1p(-y)) * h
        log_reflect = _libdevice.log(tl.minimum(tl.maximum(1.0 - upper, tiny), 1.0))
        result = -tl.where(reflect, log_reflect, log_direct) * tl.full([BLOCK], INV_LN10, tl.float64)
        tl.store(out_ptr + offsets, result.to(out_ptr.dtype.element_ty), mask=valid)


def neg_log10_p(t, df, out=None):
    """tails.neg_log10_p_device's result for (rows, cols) t at df broadcasting to it, one launch."""
    if t.ndim != 2:
        raise ValueError('t must be (rows, cols)')
    if not isinstance(df, torch.Tensor) or df.device != t.device:
        df = torch.as_tensor(df, device=t.device)
    if df.ndim == 1:
        df = df.reshape(-1, 1)
    df = df.expand(t.shape)
    if out is None:
        out = torch.empty(t.shape, dtype=torch.float32, device=t.device)
    if not out.is_contiguous() or out.shape != t.shape:
        raise ValueError('out must be contiguous and shaped like t')
    rows, cols = t.shape
    cells = rows * cols
    if cells == 0:
        return out
    block, warps = ((TAIL_BLOCK, TAIL_WARPS) if cells >= TAIL_WIDE_CELLS
                    else (TAIL_BLOCK_NARROW, TAIL_WARPS_NARROW))
    with _LAUNCH_LOCK, torch.cuda.device(t.device):
        _tail_kernel[(triton.cdiv(cells, block),)](
            t, df, out, cells, cols, t.stride(0), t.stride(1), df.stride(0), df.stride(1),
            TINY=1e-300, LGAMMA_HALF=_LGAMMA_HALF, INV_LN10=_INV_LN10,
            ITERATIONS=TAIL_ITERATIONS, BLOCK=block, num_warps=warps)
    return out


# -- complete-case downdates ---------------------------------------------------
#
# CompleteCasePlan.correct's per-pattern terms -- each pattern's residual sum
# g^T P_S g and df for every variant of the chunk -- as one launch. The Torch
# version gathers each group of patterns' missing calls into a (variants x
# cells) tensor and makes about 20 calls per group: with a pattern per voxel
# (250,000 patterns, ~1,000 groups at the default group sizes) that is 20,000
# launches per chunk, issued under the GIL, which is what made one thread per
# GPU contend. Here a program takes one pattern and a block of variants and
# walks the pattern's missing samples; the calls are read transposed, so a
# block of variants at one sample is one contiguous load, and programs of a
# variant block run together and share it in L2.
#
# FP32 products of at most CC_CELLS calls each (no TF32: input_precision
# 'ieee'), summed in FP64, as the Torch version's FP32 bmm per segment summed
# in FP64. The quadratic form u^T (Z_S^T Z_S)^+ u is FP64.
CC_VARIANTS, CC_CELLS, CC_WARPS = 64, 32, 4


if triton is not None:
    @triton.jit
    def _complete_case_kernel(calls_ptr, observed_ptr, variants, offsets_ptr, positions_ptr, samples_ptr,
                              basis_ptr, inverse_ptr, counts_ptr, subset_rank_ptr,
                              projections_ptr, ss_ptr, missing_ptr, rss_ptr, df_ptr,
                              HAS_CALLS: tl.constexpr, RANK: tl.constexpr, RANK_BLOCK: tl.constexpr,
                              BLOCK_V: tl.constexpr, BLOCK_W: tl.constexpr):
        pattern = tl.program_id(0).to(tl.int64)
        v = tl.program_id(1) * BLOCK_V + tl.arange(0, BLOCK_V)
        in_chunk = v < variants
        r = tl.arange(0, RANK_BLOCK)
        w = tl.arange(0, BLOCK_W)
        first = tl.load(offsets_ptr + pattern)
        last = tl.load(offsets_ptr + pattern + 1)
        s2 = tl.zeros([BLOCK_V], tl.float64)
        u = tl.zeros([BLOCK_V, RANK_BLOCK], tl.float64)
        absent = tl.zeros([BLOCK_V], tl.int32)
        for start in range(first, last, BLOCK_W):
            cell = start + w
            in_cells = cell < last
            position = tl.load(positions_ptr + cell, mask=in_cells, other=0).to(tl.int64)
            sample = tl.load(samples_ptr + cell, mask=in_cells, other=0).to(tl.int64)
            both = in_cells[:, None] & in_chunk[None, :]
            # (cells, variants): the transposed calls, one row per sample.
            g = tl.load(calls_ptr + position[:, None] * variants + v[None, :], mask=both, other=0.0)
            s2 += tl.sum(g * g, axis=0).to(tl.float64)
            z = tl.load(basis_ptr + sample[:, None] * RANK_BLOCK + r[None, :], mask=in_cells[:, None], other=0.0)
            u += tl.dot(tl.trans(g), z, input_precision='ieee').to(tl.float64)
            if HAS_CALLS:
                seen = tl.load(observed_ptr + position[:, None] * variants + v[None, :], mask=both, other=1)
                absent += tl.sum((seen == 0).to(tl.int32), axis=0)
        u = tl.load(projections_ptr + v[:, None] * RANK + r[None, :],
                    mask=in_chunk[:, None] & (r[None, :] < RANK), other=0.0) - u
        # u (Z_S^T Z_S)^+, a column of u at a time (the pattern's RANK x RANK inverse, FP64).
        product = tl.zeros([BLOCK_V, RANK_BLOCK], tl.float64)
        inverse = inverse_ptr + pattern * RANK * RANK
        for s in tl.static_range(RANK):
            column = tl.sum(tl.where(r[None, :] == s, u, 0.0), axis=1)
            row = tl.load(inverse + s * RANK + r, mask=r < RANK, other=0.0)
            product += column[:, None] * row[None, :]
        quadratic = tl.sum(u * product, axis=1)
        rss = tl.load(ss_ptr + v, mask=in_chunk, other=0.0) - s2 - quadratic
        samples = tl.load(counts_ptr + pattern) + tl.zeros([BLOCK_V], tl.float64)
        if HAS_CALLS:
            # Observed calls in S: the subset less the calls missing there.
            samples -= tl.load(missing_ptr + v, mask=in_chunk, other=0.0) - absent.to(tl.float64)
        df = samples - tl.load(subset_rank_ptr + pattern) - 1.0
        tl.store(rss_ptr + pattern * variants + v, rss, mask=in_chunk)
        tl.store(df_ptr + pattern * variants + v, df, mask=in_chunk)


def complete_case_terms(calls_t, observed_t, offsets, positions, samples, basis, inverse, counts, subset_rank,
                        projections, sum_squares, missing_total):
    """Every pattern's (residual sum, df) for a chunk: two (patterns, variants) FP64 tensors.

    calls_t (columns, variants) FP32, the chunk's centred calls transposed
    (a missing call zero); observed_t (columns, variants) uint8, 1 where the
    call is observed, or None when every call is. Pattern p's missing samples
    are cells offsets[p]:offsets[p + 1] of positions (their columns in
    calls_t) and samples (their rows of basis). basis (samples, rank_block)
    FP32, zero beyond the rank; inverse (patterns, rank, rank), counts and
    subset_rank (patterns,), FP64. projections (variants, rank), sum_squares
    and missing_total (variants,), FP64.
    """
    for tensor, name in ((calls_t, 'calls_t'), (offsets, 'offsets'), (positions, 'positions'),
                         (samples, 'samples'), (basis, 'basis'), (inverse, 'inverse'), (counts, 'counts'),
                         (subset_rank, 'subset_rank'), (projections, 'projections'),
                         (sum_squares, 'sum_squares'), (missing_total, 'missing_total')):
        _check(tensor, name)
    if observed_t is not None:
        _check(observed_t, 'observed_t')
    if calls_t.dtype != torch.float32 or basis.dtype != torch.float32:
        raise ValueError('complete-case terms take FP32 calls and basis')
    variants = calls_t.shape[1]
    patterns, rank = inverse.shape[0], inverse.shape[1]
    rss = torch.empty((patterns, variants), dtype=torch.float64, device=calls_t.device)
    df = torch.empty_like(rss)
    if patterns == 0 or variants == 0:
        return rss, df
    with _LAUNCH_LOCK, torch.cuda.device(calls_t.device):
        _complete_case_kernel[(patterns, triton.cdiv(variants, CC_VARIANTS))](
            calls_t, calls_t if observed_t is None else observed_t, variants, offsets, positions, samples,
            basis, inverse, counts, subset_rank, projections, sum_squares, missing_total, rss, df,
            HAS_CALLS=observed_t is not None, RANK=rank, RANK_BLOCK=basis.shape[1],
            BLOCK_V=CC_VARIANTS, BLOCK_W=CC_CELLS, num_warps=CC_WARPS)
    return rss, df
