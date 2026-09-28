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

import torch

try:
    import triton
    import triton.language as tl
except ImportError:  # the Torch statistics are the fallback
    triton = None
    tl = None

# Raw encodings: one byte or float per sample, or two-bit rows.
INT8, UINT8, FLOAT32, PGEN_2BIT, PLINK_2BIT = range(5)
# Samples per program step. A power of two, as tl.arange needs; one variant
# row per program, so a chunk launches `variants` programs per kernel.
SAMPLE_BLOCK = 1024
TRAIT_BLOCK = 1024


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
    with torch.cuda.device(device):
        _prepare_kernel[(rows,)](raw, raw.stride(0), samples, float(scale),
                                 float(0 if sentinel is None else sentinel),
                                 centered, ss, minimum, maximum, present,
                                 KIND=kind, HAS_SENTINEL=sentinel is not None, BLOCK=SAMPLE_BLOCK,
                                 num_warps=4)
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
    with torch.cuda.device(products.device):
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
    with torch.cuda.device(products.device):
        _finish_min_p_kernel[(rows,)](products, products.stride(0), centered_ss, minimum, maximum, phenotype_ss,
                                      present, float(df_offset), traits, covariates, beta, t, index, status,
                                      VALIDATE=bool(validate_range), COV_BLOCK=cov_block, BLOCK=block,
                                      num_warps=4)
    return beta, t, index, status
