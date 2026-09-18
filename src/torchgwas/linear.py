from __future__ import annotations

import math
import os
import threading
import time
from collections.abc import Iterator
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor

import numpy as np
import torch
from scipy import special

from .kernels import linear_chunk_kernel
from .preprocess import residualize_and_standardize
from .streaming import ChunkedGenotype
from .tails import upper_tail_log10_from_t_torch
from .utils import choose_device, chunk_bounds


def _resolve_compute_dtypes(compute_dtype: str) -> tuple[np.dtype, torch.dtype]:
    if compute_dtype == "float32":
        return np.float32, torch.float32
    if compute_dtype == "float64":
        return np.float64, torch.float64
    raise ValueError(f"unsupported compute dtype: {compute_dtype}")


def _two_sided_t_pvalue(t_stat: np.ndarray, df) -> np.ndarray:
    """`df` may be a scalar or one value per variant, shaped to broadcast."""
    if isinstance(df, np.ndarray) and df.ndim == 1:
        df = df[:, None]
    return 2.0 * special.stdtr(df, -np.abs(t_stat))


def _unpack_plink_a2_float(
    packed: torch.Tensor,
    n_samples: int,
    sample_byte_indices: torch.Tensor | None = None,
    sample_bit_shifts: torch.Tensor | None = None,
) -> torch.Tensor:
    if sample_byte_indices is None:
        calls = torch.stack(
            (
                packed & 3,
                (packed >> 2) & 3,
                (packed >> 4) & 3,
                (packed >> 6) & 3,
            ),
            dim=2,
        ).reshape(packed.shape[0], -1)[:, :n_samples]
    else:
        calls = (packed[:, sample_byte_indices] >> sample_bit_shifts[None, :]) & 3
    dosage = ((calls + 1) >> 1).to(torch.float32)
    return torch.where(calls == 1, torch.nan, dosage)


def _packed_bed_statistics(
    packed: torch.Tensor,
    design: torch.Tensor,
    phenotype_ss: torch.Tensor,
    n_samples: int,
    n_traits: int,
    df: int,
    sample_byte_indices: torch.Tensor | None,
    sample_bit_shifts: torch.Tensor | None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    genotype = _unpack_plink_a2_float(
        packed,
        n_samples,
        sample_byte_indices,
        sample_bit_shifts,
    )
    # df = n - rank - 2, so the rank is recoverable and the compiled path
    # uses the same per-variant df the native kernel does.
    return _dosage_statistics(genotype, design, phenotype_ss, n_traits, df,
                              covariate_rank=n_samples - df - 2)


def _significant_pairs_iterator(chunks, significance, n_traits, df):
    """Turn full chunks into the (variant, trait) pairs that clear threshold.

    **There is no k.** An earlier version kept the top `k` traits per variant
    on the device and thresholded afterwards, which is exact only while fewer
    than `k` of them clear. That assumption fails precisely where this mode is
    meant to be used: a real association in a correlated trait set can clear
    significance on hundreds of traits at once, and every one past the `k`th
    was discarded with a warning standing in for a guarantee.

    Selection is on `|t|` against a critical value rather than on a p-value per
    cell, which is the same set: within a variant every trait shares the
    residual df, so p is strictly decreasing in `|t|`. Computing p for all
    `M x K` cells instead would cost 0.156 s per million on the host -- about
    178 s at this cohort -- to answer a question two comparisons settle.

    The critical value uses the run's nominal df rather than each variant's
    own. At N = 22,250 with 27 covariates, a variant missing 1% of its calls
    has df smaller by 1%, which moves the critical `|t|` at p = 4e-10 by under
    0.01% -- orders below any threshold anyone sets, and the same
    simplification the existing `p_value_threshold` path already makes.
    """
    critical = float(significance.critical_abs_t(
        np.asarray([float(df)]), n_traits)[0])
    for start, end, beta_chunk, t_chunk, _p_chunk in chunks:
        t_values = np.asarray(t_chunk)
        keep = np.isfinite(t_values) & (np.abs(t_values) >= critical)
        if not keep.any():
            yield (start, end, np.empty(0, dtype=np.int64),
                   np.empty(0, dtype=np.int32),
                   np.empty(0, dtype=t_values.dtype),
                   np.empty(0, dtype=t_values.dtype), float(df))
            continue
        rows, columns = np.nonzero(keep)
        yield (start, end, rows.astype(np.int64) + start,
               columns.astype(np.int32),
               np.asarray(beta_chunk)[rows, columns],
               t_values[rows, columns], float(df))


def _dosage_statistics(genotype, design, phenotype_ss, n_traits, df,
                       validate_range=False, covariate_rank=None):
    """Shared variant-major OLS for native dosage and unpacked hard calls.

    Returns beta, t, status and the per-variant residual degrees of
    freedom. With `covariate_rank` given, df is the variant's own observed
    count minus the rank and two, so a masked call costs its own sample
    rather than being treated as an observation at the mean. Without it the
    scalar df is broadcast and nothing changes.
    """
    # Missing calls arrive as NaN and are masked, not dropped: they take no
    # part in the mean or the range, and their centred value is zero, so they
    # contribute nothing to any dot product below. That is the same estimand
    # as imputing the variant mean. Without a fused kernel this costs an extra
    # pass over the chunk, which the native path avoids.
    observed = ~torch.isnan(genotype)
    present = observed.sum(dim=1, keepdim=True)
    genotype = torch.where(observed, genotype, torch.zeros_like(genotype))
    huge = torch.finfo(genotype.dtype).max
    minimum = torch.where(observed, genotype, torch.full_like(genotype, huge)).amin(dim=1)
    maximum = torch.where(observed, genotype, torch.full_like(genotype, -huge)).amax(dim=1)
    # Center before the projection to avoid cancellation for nearly fixed alleles.
    # The design includes an intercept, so this preserves the exact OLS estimand.
    mean = genotype.sum(dim=1, keepdim=True) / present.clamp(min=1)
    genotype = torch.where(observed, genotype - mean, torch.zeros_like(genotype))
    products = genotype @ design
    gy = products[:, :n_traits]
    gc = products[:, n_traits:]
    residual_ss = torch.sum(genotype * genotype, dim=1) - torch.sum(gc * gc, dim=1)
    valid = (residual_ss > 1e-12) & (maximum > minimum)
    status = torch.where(
        torch.isfinite(residual_ss),
        torch.where(
            valid,
            torch.zeros_like(residual_ss, dtype=torch.uint8),
            torch.full_like(residual_ss, 2, dtype=torch.uint8),
        ),
        torch.ones_like(residual_ss, dtype=torch.uint8),
    )
    if validate_range:
        status = torch.where((minimum < 0) | (maximum > 2),
                             torch.full_like(status, 3), status)
    if covariate_rank is None:
        variant_df = torch.full_like(residual_ss, float(df))
    else:
        variant_df = present.squeeze(1).to(residual_ss.dtype) - float(covariate_rank) - 2.0
        valid = valid & (variant_df > 0)
    safe_df = torch.clamp(variant_df, min=1.0)
    safe_ss = torch.clamp(residual_ss, min=1e-12)
    beta = gy / safe_ss[:, None]
    explained_ss = gy * gy / safe_ss[:, None]
    residual_y_ss = torch.clamp(phenotype_ss[None, :] - explained_ss, min=1e-12)
    standard_error = torch.sqrt(residual_y_ss / safe_df[:, None] / safe_ss[:, None])
    t_stat = beta / standard_error
    beta = torch.where(valid[:, None], beta, torch.zeros_like(beta))
    t_stat = torch.where(valid[:, None], t_stat, torch.zeros_like(t_stat))
    return beta, t_stat, status, variant_df


def _packed_bed_statistics_native(
    packed: torch.Tensor,
    design: torch.Tensor,
    phenotype_ss: torch.Tensor,
    n_samples: int,
    n_traits: int,
    df: int,
    sample_byte_indices: torch.Tensor | None,
    sample_bit_shifts: torch.Tensor | None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """The contract of _packed_bed_statistics, served by the CUDA kernel.

    The kernel reads samples in stored order, so a subsetting request is not
    something it can answer; the caller keeps the Torch path for that.
    """
    if sample_byte_indices is not None or sample_bit_shifts is not None:
        raise ValueError("native packed BED does not support sample subsetting")
    from . import scan_gpu

    centered, centered_ss, minimum, maximum, present = scan_gpu.prepare(
        packed, encoding="plink_2bit", n_samples=n_samples)
    products = centered @ design
    beta, t_stat, status = scan_gpu.finish(
        products, centered_ss, minimum, maximum, phenotype_ss, present,
        float(df) - float(n_samples))
    return beta, t_stat, status, present.to(torch.float32) + float(df) - float(n_samples)


def _select_packed_bed_statistics(sample_byte_indices, sample_bit_shifts):
    """Prefer the kernel; fall back to the compiled Torch path.

    `TORCHGWAS_BED_NATIVE=0` forces the Torch path, which is the escape hatch
    if the kernel is ever suspected of disagreeing with it.
    """
    if (sample_byte_indices is None and sample_bit_shifts is None
            and os.environ.get("TORCHGWAS_BED_NATIVE", "1") == "1"):
        try:
            from . import scan_gpu

            if scan_gpu.available():
                return _packed_bed_statistics_native
        except Exception:  # noqa: BLE001 - availability must never fail a scan
            pass
    return _get_packed_bed_statistics()


_COMPILED_PACKED_BED_STATISTICS = None
_PACKED_BED_COMPILE_ATTEMPTED = False
_PACKED_BED_COMPILE_LOCK = threading.Lock()


def _get_packed_bed_statistics():
    """Initialize the BED compiler only when a BED scan actually needs it.

    Compilation is per chunk shape, so its cost scales with the sample count
    while the work it saves does not scale with it in the same way. Measured
    on an H100: at N = 1,000,000 the warm-up compile takes 82 s and the scan
    it prepares takes 9.7 s. `TORCHGWAS_BED_COMPILE=0` runs the same function
    eagerly instead, which is the right trade at large N and the wrong one at
    small N.
    """
    global _COMPILED_PACKED_BED_STATISTICS, _PACKED_BED_COMPILE_ATTEMPTED
    with _PACKED_BED_COMPILE_LOCK:
        if not _PACKED_BED_COMPILE_ATTEMPTED:
            _PACKED_BED_COMPILE_ATTEMPTED = True
            if os.environ.get("TORCHGWAS_BED_COMPILE", "1") == "0":
                return _packed_bed_statistics
            if hasattr(torch, "compile"):
                try:
                    _COMPILED_PACKED_BED_STATISTICS = torch.compile(
                        _packed_bed_statistics, fullgraph=True, dynamic=False)
                except Exception:
                    pass
    return _COMPILED_PACKED_BED_STATISTICS or _packed_bed_statistics


def _packed_bed_cuda_iterator(
    genotype: ChunkedGenotype,
    pheno_proc: np.ndarray,
    q_matrix: np.ndarray | None,
    chunk_size: int,
    torch_device: torch.device,
    reader_workers: int | None,
    compute_p_values: bool,
    compute_log10_p: bool = False,
    variant_range: tuple[int, int] | None = None,
    reduction=None,
    borrow_results: bool = False,
) -> Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray | None]]:
    """Overlap packed BED reads, H2D, fused GPU decode, OLS, and result copies.

    `borrow_results` yields views into the pinned result ring instead of
    copies. Measured, that copy is **69% of the wall at K=2048** — 4.92 GB at
    4.4 GB/s, larger than the GPU time and 26x the producer — so it is the
    single biggest cost in a high-trait BED scan. It is opt-in because the
    contract it relaxes is real: a borrowed chunk is valid only until the
    consumer asks for the next one.
    """

    # `torch.cuda.set_device` demands an INDEX; a bare `cuda` raises
    # "Expected a torch.device with a specified index or an integer".
    # `device="cuda"` is the natural thing for a caller to write and is what
    # `run_linear_gwas` documents, so resolving it here is the difference
    # between the library picking a card and the user being told to.
    if torch_device.type == "cuda" and torch_device.index is None:
        torch_device = torch.device("cuda", torch.cuda.current_device())
    torch.cuda.set_device(torch_device)
    n_samples, n_traits = pheno_proc.shape
    covariate_rank = 0 if q_matrix is None else q_matrix.shape[1]
    df = n_samples - covariate_rank - 2
    workers = int(reader_workers or getattr(genotype, "reader_workers", 1))
    depth = max(3, workers + 2)
    bytes_per_variant = int(getattr(genotype, "_bytes_per_variant"))

    intercept_t = torch.full(
        (n_samples, 1),
        1.0 / math.sqrt(n_samples),
        dtype=torch.float32,
        device=torch_device,
    )
    if q_matrix is None:
        covariate_t = intercept_t
    else:
        q_t = torch.as_tensor(q_matrix, dtype=torch.float32, device=torch_device)
        covariate_t = torch.cat((intercept_t, q_t), dim=1)

    # Build the design ONCE and treat its leading columns as the phenotype,
    # rather than uploading the phenotype and then concatenating a second
    # copy beside it. At voxel-scale K this is the dominant allocation on the
    # card: the old form held the phenotype (n x K), the concatenated design
    # (n x K+c) AND the full-size temporary behind `phenotype_t *
    # phenotype_t` all at once -- three times n*K*4, measured as 84 GB at
    # 33,417 samples and 150,000 voxels where the ring model predicted 21 GB.
    # Copying host-to-device straight into the slice never materialises the
    # separate phenotype tensor at all.
    covariate_width = covariate_t.shape[1]
    design_t = torch.empty((n_samples, n_traits + covariate_width),
                           dtype=torch.float32, device=torch_device)
    design_t[:, n_traits:].copy_(covariate_t)
    phenotype_t = design_t[:, :n_traits]

    # Upload and reduce in column blocks. Two reasons, both about bounding a
    # full-size temporary: copying a contiguous host array into a STRIDED
    # device slice can stage a contiguous copy of the whole thing, and
    # `(A * A).sum(0)` allocates another n*K*4. Blocked, each costs one
    # block, and the upload overlaps nothing else so the split is free.
    phenotype_ss_t = torch.empty(n_traits, dtype=torch.float32,
                                 device=torch_device)
    column_block = max(1, min(n_traits, (1 << 28) // max(n_samples * 4, 1)))
    for begin in range(0, n_traits, column_block):
        stop = min(begin + column_block, n_traits)
        columns = phenotype_t[:, begin:stop]
        columns.copy_(torch.as_tensor(
            np.ascontiguousarray(pheno_proc[:, begin:stop])))
        torch.sum(columns * columns, dim=0, out=phenotype_ss_t[begin:stop])
    selected_samples = getattr(genotype, "_sample_indices", None)
    if selected_samples is None:
        sample_byte_indices_t = None
        sample_bit_shifts_t = None
    else:
        selected_samples = np.asarray(selected_samples, dtype=np.int64)
        sample_byte_indices_t = torch.as_tensor(
            selected_samples // 4,
            dtype=torch.int64,
            device=torch_device,
        )
        sample_bit_shifts_t = torch.as_tensor(
            (selected_samples % 4) * 2,
            dtype=torch.uint8,
            device=torch_device,
        )

    packed_device = [
        torch.empty((chunk_size, bytes_per_variant), dtype=torch.uint8, device=torch_device)
        for _ in range(4)
    ]
    packed_device[0].zero_()
    # Profiling, matching the keys `dosage_cuda_iterator` records. Without this
    # a profiled BED run reported setup, fetch and GPU time as exactly 0.00,
    # because only the dosage path ever populated `_last_scan_profile` -- which
    # made the BED duty cycle, and so the question of whether BED has a
    # multi-GPU crossover at all, unmeasurable.
    profiling = os.environ.get("TORCHGWAS_SCAN_PROFILE", "0") != "0"
    setup_started = time.perf_counter()
    timings = dict(enabled=profiling, setup_seconds=0.0, fetch_seconds=0.0,
                   result_wait_seconds=0.0, copy_wait_seconds=0.0,
                   gpu_compute_milliseconds=0.0, gpu_result_milliseconds=0.0,
                   gpu_input_conversion_milliseconds=0.0, chunks=0,
                   result_copy_seconds=0.0, statistics_backend="packed_bed")
    genotype._last_scan_profile = timings
    compute_done = [torch.cuda.Event() for _ in packed_device]
    copy_done = [torch.cuda.Event() for _ in packed_device]
    # Reused per device rather than created per scan: the caching allocator
    # keys cached blocks by stream, so a fresh stream each scan means the
    # previous scan's device pool cannot be reused. This path reserves far less
    # than the PGEN one -- BED's packed rows are a quarter the width -- which is
    # why it grew 1.26 GB over four scans where PGEN grew 62.2, but the
    # mechanism is the same and so is the remedy.
    from .native_scan import _scan_streams
    copy_stream, result_stream = _scan_streams(torch_device)
    # BED reading and Student-t evaluation have different scaling curves.  A
    # local NVMe device is often saturated by 4 readers, while exact tails can
    # still use many CPU cores.  Size the result ring independently, but cap its
    # pinned-memory footprint because a high-trait scan has much larger slots.
    if hasattr(os, "sched_getaffinity"):
        available_cpus = len(os.sched_getaffinity(0))
    else:
        available_cpus = os.cpu_count() or 1
    pvalue_workers = (
        max(1, min(24, available_cpus))
        if compute_p_values
        else max(1, min(8, workers))
    )
    # A reduced scan stages k columns, not K, so both the per-slot footprint and
    # the cap it is measured against shrink with it.
    reduction_width = None if reduction is None else reduction.resolved_width(n_traits)
    staged_traits = n_traits if reduction is None else reduction_width
    logp_bytes = staged_traits * 8 if compute_log10_p else 0
    result_slot_bytes = chunk_size * (staged_traits * 2 * 4 + logp_bytes + 1)
    max_slots_by_memory = max(2, (1 << 30) // max(1, result_slot_bytes))
    result_depth = max(2, min(pvalue_workers, max_slots_by_memory))
    beta_host = [
        torch.empty((chunk_size, staged_traits), dtype=torch.float32, pin_memory=True)
        for _ in range(result_depth)
    ]
    t_host = [
        torch.empty((chunk_size, staged_traits), dtype=torch.float32, pin_memory=True)
        for _ in range(result_depth)
    ]
    logp_host = (
        [torch.empty((chunk_size, staged_traits), dtype=torch.float64, pin_memory=True)
         for _ in range(result_depth)]
        if compute_log10_p else None
    )
    index_host = (
        None if reduction is None else [
            torch.empty((chunk_size, staged_traits), dtype=torch.int32, pin_memory=True)
            for _ in range(result_depth)
        ]
    )
    status_host = [
        torch.empty(chunk_size, dtype=torch.uint8, pin_memory=True)
        for _ in range(result_depth)
    ]
    df_host = [
        torch.empty(chunk_size, pin_memory=True)
        for _ in range(result_depth)
    ]
    result_done = [torch.cuda.Event(enable_timing=profiling)
                   for _ in range(result_depth)]
    # Timing events are indexed by *result* slot, not device slot, and are kept
    # separate from the synchronisation events above. The device ring turns over
    # every `len(packed_device)` iterations while `finish_result` runs on the
    # `result_depth` cadence, so reading `compute_done[device_slot]` from there
    # can read an event a later iteration has already re-recorded -- which is
    # not a wrong number but a hard `CUDA error: device not ready`. Indexed by
    # result slot, a timing event's lifetime matches exactly the one future that
    # reads it.
    result_start = ([torch.cuda.Event(enable_timing=True)
                     for _ in range(result_depth)] if profiling else [])
    compute_start = ([torch.cuda.Event(enable_timing=True)
                      for _ in range(result_depth)] if profiling else [])
    compute_end = ([torch.cuda.Event(enable_timing=True)
                    for _ in range(result_depth)] if profiling else [])

    statistics = _select_packed_bed_statistics(
        sample_byte_indices_t, sample_bit_shifts_t)
    try:
        warm_beta, warm_t, warm_status, warm_df = statistics(
            packed_device[0],
            design_t,
            phenotype_ss_t,
            n_samples,
            n_traits,
            df,
            sample_byte_indices_t,
            sample_bit_shifts_t,
        )
        (warm_beta.sum() + warm_t.sum() + warm_status.sum()
         + warm_df.sum()).item()
    except Exception:
        statistics = _packed_bed_statistics
    torch.cuda.synchronize(torch_device)

    loader = genotype.iter_packed_chunks(
        chunk_size=chunk_size,
        reader_workers=workers,
        depth=depth,
        variant_range=variant_range,
    )
    # Record the geometry this path RESOLVED, the way `native_scan` does.
    # Without it the fused BED scan reports chunk/depth/transport as None,
    # so a predicted ring can only be compared against a configuration that
    # never ran -- which is exactly how the memory model came to be checked
    # against chunk 4096 / depth 32 when BED actually uses
    # `preferred_gpu_chunk_size` (5000) and `max(3, workers + 2)`.
    timings['chunk_variants'] = int(loader.chunk_size)
    timings['depth'] = int(loader.depth)
    timings['decode_on_gpu'] = True
    timings['transfer_bytes_per_variant'] = int(
        getattr(genotype, '_bytes_per_variant', 0))
    timings['native_encoding'] = 'plink_2bit'
    pvalue_pool = ThreadPoolExecutor(
        max_workers=result_depth,
        thread_name_prefix="torchgwas-pvalue",
    )
    pending: deque[Future] = deque()
    exclusion_lock = threading.Lock()
    genotype._last_scan_exclusion_counts = {"missing": 0, "invariant": 0}

    def finish_result(result_slot: int, start: int, end: int):
        waited = time.perf_counter() if profiling else 0.0
        result_done[result_slot].synchronize()
        if profiling:
            # Every event read here is indexed by this same result slot and was
            # recorded before `result_done`, on streams `result_stream` waits
            # on -- so completing `result_done` completes all of them, and none
            # can have been re-recorded, since the slot does not come round
            # again until this future has been consumed.
            timings["result_wait_seconds"] += time.perf_counter() - waited
            timings["chunks"] += 1
            timings["gpu_compute_milliseconds"] += (
                compute_start[result_slot].elapsed_time(compute_end[result_slot]))
            timings["gpu_result_milliseconds"] += (
                result_start[result_slot].elapsed_time(result_done[result_slot]))
        count = end - start
        finish_started = time.perf_counter() if profiling else 0.0
        # Copy: these are views into a ring the pipeline reuses, and callers
        # may retain a chunk past that point.
        status_chunk = status_host[result_slot][:count].numpy().copy()
        df_chunk = df_host[result_slot][:count].numpy().copy()
        n_missing = int(np.count_nonzero(status_chunk == 1))
        n_invariant = int(np.count_nonzero(status_chunk == 2))
        if n_missing or n_invariant:
            with exclusion_lock:
                genotype._last_scan_exclusion_counts["missing"] += n_missing
                genotype._last_scan_exclusion_counts["invariant"] += n_invariant
        if borrow_results:
            # Views into the ring, not copies. Safe by the loop's own ordering:
            # slot `i % result_depth` is overwritten at iteration
            # `i + result_depth`, and the main loop calls `.result()` on that
            # slot's future and *yields it* before overwriting -- so the
            # consumer has already had the chunk and, being a generator's
            # consumer, has finished with it before asking for the next.
            #
            # Only for consumers that do not retain: `_drain_linear_chunks` and
            # the two text writers. The binary writer hands blocks to a writer
            # thread as their owner, and the trait-blocked accumulator keeps
            # tensors across iterations, so both still get copies.
            beta_chunk = beta_host[result_slot][:count].numpy()
            t_chunk = t_host[result_slot][:count].numpy()
        else:
            beta_chunk = beta_host[result_slot][:count].numpy().copy()
            t_chunk = t_host[result_slot][:count].numpy().copy()
        if profiling:
            # Charged separately because it is the largest term in a BED scan
            # and was invisible until this was added: two `chunk x traits`
            # copies out of the pinned ring, which at K=2048 over 300,000
            # variants is 4.9 GB of host memcpy.
            timings["result_copy_seconds"] += time.perf_counter() - finish_started
        invalid = status_chunk != 0
        if np.any(invalid):
            beta_chunk[invalid, :] = np.nan
            t_chunk[invalid, :] = np.nan
        if reduction is None:
            if compute_log10_p:
                logp_chunk = logp_host[result_slot][:count].numpy().copy()
                if np.any(invalid):
                    logp_chunk[invalid, :] = np.nan
                p_chunk = None
                if compute_p_values:
                    with np.errstate(under="ignore"):
                        p_chunk = np.power(10.0, -logp_chunk)
                return start, end, beta_chunk, t_chunk, p_chunk, logp_chunk
            p_chunk = (_two_sided_t_pvalue(t_chunk, df=df_chunk)
                       if compute_p_values else None)
            return start, end, beta_chunk, t_chunk, p_chunk
        index_chunk = index_host[result_slot][:count].numpy().copy()
        # df is per variant; broadcast it across the k traits kept so the
        # reduced p-values equal the unreduced ones exactly.
        p_chunk = (_two_sided_t_pvalue(t_chunk, df=df_chunk[:, None])
                   if compute_p_values else None)
        return start, end, beta_chunk, t_chunk, p_chunk, index_chunk

    try:
        if profiling:
            timings["setup_seconds"] = time.perf_counter() - setup_started

        def timed_loader():
            """Charge the wait for the next chunk to `fetch_seconds`.

            The producer runs ahead, so this is the time the consumer actually
            spends blocked on it -- which is the quantity the crossover model
            needs, not the producer's own wall time.
            """
            iterator = iter(loader)
            while True:
                started = time.perf_counter()
                try:
                    item = next(iterator)
                except StopIteration:
                    timings["fetch_seconds"] += time.perf_counter() - started
                    return
                timings["fetch_seconds"] += time.perf_counter() - started
                yield item

        for iteration, (buffer_index, host_packed, start, end) in enumerate(
                timed_loader() if profiling else loader):
            if len(pending) >= result_depth:
                yield pending.popleft().result()

            count = end - start
            device_slot = iteration % len(packed_device)
            if iteration >= len(packed_device):
                copy_stream.wait_event(compute_done[device_slot])
            with torch.cuda.stream(copy_stream):
                packed_device[device_slot].copy_(host_packed, non_blocking=True)
                copy_done[device_slot].record(copy_stream)

            compute_stream = torch.cuda.current_stream(torch_device)
            compute_stream.wait_event(copy_done[device_slot])
            result_slot = iteration % result_depth
            if profiling:
                compute_start[result_slot].record(compute_stream)
            beta_t, t_t, status_t, df_t = statistics(
                packed_device[device_slot],
                design_t,
                phenotype_ss_t,
                n_samples,
                n_traits,
                df,
                sample_byte_indices_t,
                sample_bit_shifts_t,
            )
            index_t = None
            logp_t = None
            if reduction is not None:
                # On the compute stream and before compute_done, so the copy
                # stream waits on the narrow result and the wide tensors are
                # released here instead of crossing the bus.
                beta_t, t_t, index_t, status_t, df_t = reduction.reduce(
                    beta_t[:count], t_t[:count], status_t[:count], df_t[:count],
                    reduction_width)
            elif compute_log10_p:
                logp_t = upper_tail_log10_from_t_torch(
                    t_t[:count], df_t[:count, None])
            if profiling:
                compute_end[result_slot].record(compute_stream)
            compute_done[device_slot].record(compute_stream)

            with torch.cuda.stream(result_stream):
                result_stream.wait_event(compute_done[device_slot])
                if profiling:
                    result_start[result_slot].record(result_stream)
                if reduction is None:
                    staged = ((beta_host, beta_t), (t_host, t_t),
                              (status_host, status_t), (df_host, df_t))
                    if compute_log10_p:
                        staged += ((logp_host, logp_t),)
                else:
                    # Already sliced to `count` by the reduction above.
                    staged = ((beta_host, beta_t), (t_host, t_t),
                              (index_host, index_t), (status_host, status_t),
                              (df_host, df_t))
                for ring, value in staged:
                    ring[result_slot][:count].copy_(
                        value if reduction is not None else value[:count],
                        non_blocking=True)
                    value.record_stream(result_stream)
                result_done[result_slot].record(result_stream)

            # Free the packed host buffer as soon as DMA no longer reads it.
            copy_done[device_slot].synchronize()
            loader.release(buffer_index)
            pending.append(pvalue_pool.submit(finish_result, result_slot, start, end))

        while pending:
            yield pending.popleft().result()
    finally:
        pvalue_pool.shutdown(wait=True, cancel_futures=True)
        loader.close()


def linear_scan(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
    compute_dtype: str = "float64",
    return_log10_p: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pheno_proc, q_matrix, phenotype_observed_counts = residualize_and_standardize(
        phenotype, covariates, return_observed_counts=True)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    n_traits = pheno_proc.shape[1]
    torch_device = choose_device(device)
    np_dtype, torch_dtype = _resolve_compute_dtypes(compute_dtype)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch_dtype, device=torch_device)
    q_t = None if q_matrix is None else torch.as_tensor(q_matrix, dtype=torch_dtype, device=torch_device)
    covariate_rank = 0 if q_matrix is None else q_matrix.shape[1]
    df = n_samples - covariate_rank - 2
    if df <= 0:
        raise ValueError(f"non-positive residual degrees of freedom: N={n_samples}, covariate_rank={covariate_rank}")
    trait_df = phenotype_observed_counts.astype(np.float64) - covariate_rank - 2
    if np.any(trait_df <= 0):
        raise ValueError("non-positive phenotype-specific residual degrees of freedom")
    trait_scale = np.sqrt(trait_df / float(df))
    trait_scale_t = torch.as_tensor(
        trait_scale, dtype=torch_dtype, device=torch_device
    ) if return_log10_p else None
    trait_df_factor_t = torch.as_tensor(
        trait_df / float(df), dtype=torch.float64, device=torch_device
    ) if return_log10_p else None

    beta = np.empty((n_markers, n_traits), dtype=np_dtype)
    t_stat = np.empty((n_markers, n_traits), dtype=np_dtype)
    p_value = (
        None if return_log10_p
        else np.empty((n_markers, n_traits), dtype=np.float64)
    )
    log10_p = (
        np.empty((n_markers, n_traits), dtype=np.float64)
        if return_log10_p else None
    )

    for start, end in chunk_bounds(n_markers, chunk_size):
        geno_chunk = np.asarray(genotype[:, start:end], dtype=np_dtype)
        geno_t = torch.as_tensor(geno_chunk, dtype=torch_dtype, device=torch_device)
        beta_t, t_chunk_t, df_chunk_t = linear_chunk_kernel(
            geno_t, pheno_t, q_t, df, covariate_rank=covariate_rank)
        beta_chunk = beta_t.cpu().numpy()
        if return_log10_p:
            adjusted_t = t_chunk_t * trait_scale_t[None, :]
            pair_df_t = df_chunk_t.double()[:, None] * trait_df_factor_t[None, :]
            logp_chunk = upper_tail_log10_from_t_torch(adjusted_t, pair_df_t)
            t_chunk = adjusted_t.cpu().numpy()
            logp_chunk = logp_chunk.cpu().numpy()
        else:
            t_chunk = t_chunk_t.cpu().numpy()
            variant_df = df_chunk_t.cpu().numpy()
            t_chunk *= trait_scale[None, :]
            pair_df = variant_df[:, None] * (trait_df[None, :] / float(df))
            p_chunk = _two_sided_t_pvalue(t_chunk, df=pair_df)
        beta[start:end] = beta_chunk
        t_stat[start:end] = t_chunk
        if p_value is not None:
            p_value[start:end] = p_chunk
        if return_log10_p:
            log10_p[start:end] = logp_chunk
    if return_log10_p:
        return beta, t_stat, log10_p, q_matrix
    return beta, t_stat, p_value, q_matrix


def linear_scan_streaming(
    genotype: ChunkedGenotype,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
    compute_dtype: str = "float32",
    reader_workers: int | None = None,
    prefetch_chunks: int | None = None,
    return_log10_p: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    requested_device = choose_device(device)
    if (
        requested_device.type == "cuda"
        and compute_dtype == "float32"
        and (hasattr(genotype, "iter_packed_chunks") or getattr(genotype, "supports_fused_qc", False))
    ):
        chunk_iterator, q_matrix = linear_scan_streaming_chunks(
            genotype,
            phenotype,
            covariates,
            chunk_size=chunk_size,
            device=str(requested_device),
            compute_dtype=compute_dtype,
            reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
            compute_p_values=not return_log10_p,
            compute_log10_p=return_log10_p,
        )
        n_markers = genotype.shape[1]
        n_traits = phenotype.shape[1]
        beta = np.empty((n_markers, n_traits), dtype=np.float32)
        t_stat = np.empty((n_markers, n_traits), dtype=np.float32)
        p_value = (
            None if return_log10_p
            else np.empty((n_markers, n_traits), dtype=np.float64)
        )
        log10_p = (
            np.empty((n_markers, n_traits), dtype=np.float64)
            if return_log10_p else None
        )
        for chunk_result in chunk_iterator:
            if return_log10_p:
                start, end, beta_chunk, t_chunk, _p_chunk, logp_chunk = chunk_result
            else:
                start, end, beta_chunk, t_chunk, p_chunk = chunk_result
            beta[start:end] = beta_chunk
            t_stat[start:end] = t_chunk
            if p_value is not None:
                p_value[start:end] = p_chunk
            if return_log10_p:
                log10_p[start:end] = logp_chunk
        if return_log10_p:
            return beta, t_stat, log10_p, q_matrix
        return beta, t_stat, p_value, q_matrix

    pheno_proc, q_matrix, phenotype_observed_counts = residualize_and_standardize(
        phenotype, covariates, return_observed_counts=True)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    n_traits = pheno_proc.shape[1]
    torch_device = requested_device
    covariate_rank = 0 if q_matrix is None else q_matrix.shape[1]
    chunk = chunk_size or _default_chunk_variants(
        genotype, torch_device, n_samples, n_markers,
        pheno_proc.shape[1], covariate_rank, prefetch_chunks,
    )
    df = n_samples - covariate_rank - 2
    if df <= 0:
        raise ValueError(f"non-positive residual degrees of freedom: N={n_samples}, covariate_rank={covariate_rank}")
    trait_df = phenotype_observed_counts.astype(np.float64) - covariate_rank - 2
    if np.any(trait_df <= 0):
        raise ValueError("non-positive phenotype-specific residual degrees of freedom")
    trait_scale = np.sqrt(trait_df / float(df))

    np_dtype, torch_dtype = _resolve_compute_dtypes(compute_dtype)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch_dtype, device=torch_device)
    q_t = None if q_matrix is None else torch.as_tensor(q_matrix, dtype=torch_dtype, device=torch_device)

    beta = np.empty((n_markers, n_traits), dtype=np_dtype)
    t_stat = np.empty((n_markers, n_traits), dtype=np_dtype)
    p_value = (
        None if return_log10_p
        else np.empty((n_markers, n_traits), dtype=np.float64)
    )
    log10_p = (
        np.empty((n_markers, n_traits), dtype=np.float64)
        if return_log10_p else None
    )
    trait_scale_t = torch.as_tensor(
        trait_scale, dtype=torch_dtype, device=torch_device
    ) if return_log10_p else None
    trait_df_factor_t = torch.as_tensor(
        trait_df / float(df), dtype=torch.float64, device=torch_device
    ) if return_log10_p else None

    for start, end, geno_chunk in genotype.iter_chunks(
        chunk_size=chunk,
        dtype=np_dtype,
        prefetch_chunks=prefetch_chunks,
        reader_workers=reader_workers,
    ):
        geno_t = torch.as_tensor(geno_chunk, dtype=torch_dtype, device=torch_device)
        beta_t, t_chunk_t, df_chunk_t = linear_chunk_kernel(
            geno_t, pheno_t, q_t, df, covariate_rank=covariate_rank)
        beta_chunk = beta_t.cpu().numpy()
        if return_log10_p:
            adjusted_t = t_chunk_t * trait_scale_t[None, :]
            pair_df_t = df_chunk_t.double()[:, None] * trait_df_factor_t[None, :]
            logp_chunk = upper_tail_log10_from_t_torch(adjusted_t, pair_df_t)
            t_chunk = adjusted_t.cpu().numpy()
            logp_chunk = logp_chunk.cpu().numpy()
        else:
            t_chunk = t_chunk_t.cpu().numpy()
            variant_df = df_chunk_t.cpu().numpy()
            t_chunk *= trait_scale[None, :]
            pair_df = variant_df[:, None] * (trait_df[None, :] / float(df))
            p_chunk = _two_sided_t_pvalue(t_chunk, df=pair_df)
        beta[start:end] = beta_chunk
        t_stat[start:end] = t_chunk
        if p_value is not None:
            p_value[start:end] = p_chunk
        if return_log10_p:
            log10_p[start:end] = logp_chunk
    if return_log10_p:
        return beta, t_stat, log10_p, q_matrix
    return beta, t_stat, p_value, q_matrix


def linear_scan_multigpu(
    genotype,
    phenotype,
    covariates=None,
    *,
    devices=None,
    chunk_size=None,
    reader_workers=None,
    prefetch_chunks=None,
    compute_p_values=True,
    **kwargs,
):
    """Run one scan across several GPUs, sharded by variant range.

    `devices` defaults to every visible CUDA device. With one device this is
    the ordinary single-device scan, so callers need not special-case it.

    Yields the same `(start, end, beta, t, p)` tuples as the single-device
    scan, in the file's own variant numbering and in variant order.
    """
    import threading
    import queue as queue_module

    if devices is None:
        count = torch.cuda.device_count() if torch.cuda.is_available() else 0
        devices = [f"cuda:{index}" for index in range(count)] or ["cpu"]
    devices = [str(device) for device in devices]

    n_markers = int(genotype.shape[1])
    if len(devices) < 2 or n_markers < len(devices):
        iterator, q_matrix = linear_scan_streaming_chunks(
            genotype, phenotype, covariates, chunk_size=chunk_size,
            device=devices[0], reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
            compute_p_values=compute_p_values, **kwargs)
        return iterator, q_matrix

    # Divide the reader threads among the shards instead of giving each the
    # full count: the file and the disk are shared, so replicating readers
    # multiplies contention without adding bandwidth.
    if reader_workers is not None:
        reader_workers = max(1, int(reader_workers) // len(devices))

    # Residualise once. Every shard needs the same processed phenotype, and at
    # a large trait count recomputing it per shard is both slow and large.
    shared_pheno, shared_q = residualize_and_standardize(phenotype, covariates)

    # Contiguous shards whose boundaries land on chunk multiples, so the scan
    # decomposes into exactly the chunks a single-device run would use and the
    # results agree bit for bit however many GPUs are visible.
    grain = int(chunk_size) if chunk_size else 1
    per_shard = max(grain, ((n_markers // len(devices)) // grain) * grain)
    bounds = []
    cursor = 0
    for index in range(len(devices)):
        last = index == len(devices) - 1
        end = n_markers if last else min(n_markers, cursor + per_shard)
        if cursor >= n_markers:
            break
        bounds.append((cursor, end))
        cursor = end
    devices = devices[:len(bounds)]

    # One bounded queue per shard. A shard runs ahead into its own queue and
    # blocks only when that queue is full, which cannot stall another shard.
    shard_queues = [queue_module.Queue(maxsize=4) for _ in devices]
    errors: list[BaseException] = []
    q_holder: dict[str, object] = {}
    finished = object()

    def run_shard(index, device, span):
        outbox = shard_queues[index]
        try:
            iterator, q_matrix = linear_scan_streaming_chunks(
                genotype, shared_pheno, shared_q, chunk_size=chunk_size,
                device=device, reader_workers=reader_workers,
                prefetch_chunks=prefetch_chunks,
                compute_p_values=compute_p_values,
                variant_range=span, already_processed=True, **kwargs)
            q_holder.setdefault("q", q_matrix)
            for item in iterator:
                outbox.put(item)
        except BaseException as error:  # noqa: BLE001 - re-raised to the caller
            errors.append(error)
        finally:
            outbox.put(finished)

    threads = [threading.Thread(target=run_shard, args=(index, device, span),
                                name=f"torchgwas-shard-{index}", daemon=True)
               for index, (device, span) in enumerate(zip(devices, bounds))]
    for thread in threads:
        thread.start()

    def drain():
        try:
            # Shards cover ascending contiguous ranges and each shard's chunks
            # arrive in order, so draining them in turn is already sorted.
            for outbox in shard_queues:
                while True:
                    item = outbox.get()
                    if item is finished:
                        break
                    yield item
        finally:
            for thread in threads:
                thread.join()
        if errors:
            raise errors[0]

    return drain(), q_holder.get("q")


def _default_chunk_variants(genotype, device, n_samples, n_markers, n_traits,
                            covariate_rank, prefetch_chunks):
    """Pick a chunk size for this shape rather than a per-source constant.

    The stored constants (5000 for BED, 4096 otherwise) do not depend on the
    sample count, but every per-chunk buffer does. At 22,250 samples a
    5000-variant chunk is a few hundred megabytes; at a million samples the
    same chunk is tens of gigabytes per ring slot and the scan cannot start.

    Variants are the axis to cut: a regression needs every sample of its
    variant, so samples cannot be split, and cutting traits shrinks only the
    result and projection terms, which stay small until K is very large.

    On CUDA the chunk is derived from free device memory. Elsewhere, or if
    anything about the query fails, the previous constants are used, so a
    machine this cannot measure behaves exactly as it did before.
    """
    fallback = (
        (getattr(genotype, "preferred_gpu_chunk_size", None)
         if device.type == "cuda" else None)
        or getattr(genotype, "preferred_chunk_size", None)
        or min(n_markers, 4096)
        or 1
    )
    if device.type != "cuda":
        return fallback
    try:
        from .pipeline_model import auto_chunk_variants

        free_bytes, _total = torch.cuda.mem_get_info(device)
        # What actually crosses the bus, which is not the same as the decoded
        # width. A packed two-bit transport sends a quarter of a byte per
        # sample where float32 sends four -- a 16x difference in every ring
        # buffer -- so reading this wrong makes the chunk 16x too small or too
        # large.
        #
        # This used to special-case `pgen_2bit` alone, which left PLINK BED --
        # equally two-bit on the wire, via `tg_scan_prepare_bed2` -- sized as
        # though it were float32: chunk 2319 instead of 4342 at 35,365 samples,
        # a ring of 10.5 GB where 1.2 GB would do.
        width = getattr(genotype, "native_row_width", 0)
        packed = getattr(genotype, "native_encoding", None) in {"pgen_2bit",
                                                                "plink_2bit"}
        if packed and width:
            transfer_bytes = float(width)
        elif hasattr(genotype, "iter_packed_chunks"):
            # BED advertises no `native_row_width`, and `iter_packed_chunks` is
            # what actually gates its two-bit path (see the dispatch below and
            # `scan_gpu.prepare(..., encoding="plink_2bit")`). Its wire format
            # is one two-bit call per sample, padded to a byte.
            transfer_bytes = float((n_samples + 3) // 4)
        else:
            itemsize = np.dtype(
                getattr(genotype, "native_dtype", np.float32)
            ).itemsize
            transfer_bytes = float(n_samples) * float(itemsize)
        derived = auto_chunk_variants(
            n_samples=n_samples,
            n_traits=n_traits,
            covariate_rank=covariate_rank,
            transfer_bytes_per_variant=transfer_bytes,
            device_memory_bytes=int(free_bytes),
            depth=max(1, int(prefetch_chunks or 4)),
            # A source that decodes on the device holds an extra ring slot and
            # a decode tile alongside the chunk, so its ring is
            # `(depth + 1) * max(tile, chunk)` rather than `depth * chunk`.
            # Omitting this sized BGEN against the host-decode formula and
            # under-counted its device memory.
            decode_on_gpu=hasattr(genotype, "iter_device_chunks"),
            # A framed source decompresses whole frames, so a chunk that
            # straddles them pays for what it discards -- measured at 4.6x for
            # chunk 1,536 against a 2,048 frame. The bisection above returns
            # arbitrary integers, so without this the chunk lands in the
            # penalised class by luck of the memory arithmetic.
            frame_variants=getattr(genotype, "chunk_alignment_variants", None),
        )
    except Exception:  # noqa: BLE001 - never fail a scan over a size heuristic
        return fallback
    # `min` can re-break the alignment when the file is shorter than one
    # chunk, but then there is a single chunk covering everything and framing
    # is irrelevant -- the whole file is read either way.
    return max(1, min(derived, n_markers))


def linear_scan_streaming_chunks(
    genotype: ChunkedGenotype,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
    compute_dtype: str = "float32",
    reader_workers: int | None = None,
    prefetch_chunks: int | None = None,
    compute_p_values: bool = True,
    compute_log10_p: bool = False,
    variant_range: tuple[int, int] | None = None,
    already_processed: bool = False,
    reduction=None,
    borrow_results: bool = False,
) -> tuple[
    Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray | None]],
    np.ndarray | None,
]:
    if compute_log10_p and reduction is not None:
        raise ValueError("compute_log10_p is supported only for full result matrices")

    # With `reduction` set, every backend yields a six-tuple
    # `(start, end, beta, t, p, trait_index)` whose trait axis is k wide instead
    # of K, and the wide matrix is never copied off the device. The five-tuple
    # contract is untouched when it is None, so no existing consumer changes.
    # `already_processed` lets a caller that residualises once -- the multi-GPU
    # driver -- hand the same result to every shard instead of each shard
    # recomputing it.
    # The device is resolved before residualising, not after, so the projection
    # can run there. It is one pass through a 27-column basis -- milliseconds on
    # a GPU that is idle at this point in the run -- against a host version
    # measured between 0.37 s and 10.26 s for the same call at 512 traits, the
    # spread coming from torch's thread pool degrading numpy's threaded GEMM.
    torch_device = choose_device(device)
    if already_processed:
        pheno_proc, q_matrix = phenotype, covariates
    else:
        pheno_proc, q_matrix, phenotype_observed_counts = residualize_and_standardize(
            phenotype, covariates, device=torch_device,
            return_observed_counts=True)
    if already_processed:
        phenotype_observed_counts = np.full(
            pheno_proc.shape[1], pheno_proc.shape[0], dtype=np.int64)
    # A reduction that needs a whole-trait-set quantity -- jagwas needs the
    # K x K trait correlation and its Cholesky factor -- gets one chance to
    # compute it, here, from the *processed* design the scan will actually use.
    # Deriving it from the raw phenotype instead would use different columns.
    if reduction is not None and hasattr(reduction, "prepare"):
        reduction.prepare(pheno_proc, device=torch_device)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    covariate_rank = 0 if q_matrix is None else q_matrix.shape[1]
    chunk = chunk_size or _default_chunk_variants(
        genotype, torch_device, n_samples, n_markers,
        pheno_proc.shape[1], covariate_rank, prefetch_chunks,
    )
    df = n_samples - covariate_rank - 2
    if df <= 0:
        raise ValueError(f"non-positive residual degrees of freedom: N={n_samples}, covariate_rank={covariate_rank}")
    trait_df = phenotype_observed_counts.astype(np.float64) - covariate_rank - 2
    if np.any(trait_df <= 0):
        raise ValueError("non-positive phenotype-specific residual degrees of freedom")
    phenotype_has_missing = bool(np.any(phenotype_observed_counts != n_samples))
    if phenotype_has_missing and reduction is not None:
        raise ValueError(
            "phenotype missingness is supported in full mode; trait reductions "
            "need pair-specific df and are not yet supported")
    trait_scale = np.sqrt(trait_df / float(df)).astype(np.float32)

    def apply_phenotype_missingness(iterator):
        if not phenotype_has_missing:
            yield from iterator
            return
        for chunk_result in iterator:
            if compute_log10_p:
                start, end, beta, t_stat, _p, _logp = chunk_result
            else:
                start, end, beta, t_stat, _p = chunk_result
            t_stat *= trait_scale[None, :]
            if compute_log10_p:
                t_tail = torch.as_tensor(
                    t_stat, dtype=torch.float64, device=torch_device)
                df_tail = torch.as_tensor(
                    trait_df[None, :], dtype=torch.float64,
                    device=torch_device)
                logp = upper_tail_log10_from_t_torch(t_tail, df_tail).cpu().numpy()
                p_value = None
                if compute_p_values:
                    with np.errstate(under="ignore"):
                        p_value = np.power(10.0, -logp)
                yield start, end, beta, t_stat, p_value, logp
            else:
                p_value = (
                    _two_sided_t_pvalue(t_stat, df=trait_df[None, :])
                    if compute_p_values else None
                )
                yield start, end, beta, t_stat, p_value
    np_dtype, torch_dtype = _resolve_compute_dtypes(compute_dtype)

    # Both CUDA fast paths stage results through a preallocated `chunk x k`
    # pinned ring, so a significance threshold -- whose result length is a
    # property of the data -- does not fit here. It is not handled by a second
    # transport either: an earlier version routed it to the generic device path
    # and that measured **307x slower** (1362.75 s against 4.44 s on the same
    # 200,000-variant scan). `run_linear_gwas` instead turns a significance
    # request into an ordinary top-k reduction, which this path already runs at
    # full speed, and applies the threshold in the writer. Exact whenever fewer
    # than k traits at a variant clear it, and the writer reports the variants
    # where it cannot be sure.
    if (
        torch_device.type == "cuda"
        and compute_dtype == "float32"
        and hasattr(genotype, "iter_packed_chunks")
    ):
        return (
            apply_phenotype_missingness(_packed_bed_cuda_iterator(
                genotype,
                pheno_proc,
                q_matrix,
                chunk,
                torch_device,
                reader_workers,
                compute_p_values,
                compute_log10_p,
                variant_range,
                reduction=reduction,
                borrow_results=borrow_results,
            )),
            q_matrix,
        )

    if (torch_device.type == "cuda" and compute_dtype == "float32"
                and getattr(genotype, "supports_fused_qc", False)):
        from .native_scan import dosage_cuda_iterator
        return apply_phenotype_missingness(dosage_cuda_iterator(genotype, pheno_proc, q_matrix, chunk,
                                    torch_device, reader_workers, prefetch_chunks,
                                    compute_p_values,
                                    compute_log10_p,
                                    variant_range=variant_range,
                                    reduction=reduction)), q_matrix

    pheno_t = torch.as_tensor(pheno_proc, dtype=torch_dtype, device=torch_device)
    q_t = None if q_matrix is None else torch.as_tensor(q_matrix, dtype=torch_dtype, device=torch_device)

    def _iterator() -> Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray | None]]:
        # Only ask for a range when one was requested: not every source
        # implements the parameter, and an unranged scan must call as before.
        extra = {} if variant_range is None else {"variant_range": variant_range}
        for start, end, geno_chunk in genotype.iter_chunks(
            chunk_size=chunk,
            dtype=np_dtype,
            prefetch_chunks=prefetch_chunks,
            reader_workers=reader_workers,
            **extra,
        ):
            geno_t = torch.as_tensor(geno_chunk, dtype=torch_dtype, device=torch_device)
            beta_t, t_chunk_t, df_chunk_t = linear_chunk_kernel(
                geno_t, pheno_t, q_t, df, covariate_rank=covariate_rank)
            if reduction is not None:
                # This backend has no per-variant status word -- degenerate
                # variants arrive as NaN rather than as a flag -- so pass an
                # all-clear status and let the reducer's NaN guard handle them.
                width = reduction.resolved_width(beta_t.shape[1])
                clear = torch.zeros(beta_t.shape[0], dtype=torch.uint8,
                                    device=beta_t.device)
                beta_t, t_chunk_t, index_t, _, _ = reduction.reduce(
                    beta_t, t_chunk_t, clear, df_chunk_t, width)
                index_chunk = index_t.cpu().numpy()
            beta_chunk = beta_t.cpu().numpy()
            t_chunk = t_chunk_t.cpu().numpy()
            df_chunk = df_chunk_t.cpu().numpy()
            if reduction is not None:
                # df is per variant, so it broadcasts across the k kept traits.
                p_chunk = (_two_sided_t_pvalue(t_chunk, df=df_chunk[:, None])
                           if compute_p_values else None)
                yield start, end, beta_chunk, t_chunk, p_chunk, index_chunk
                continue
            if compute_log10_p:
                logp_chunk = upper_tail_log10_from_t_torch(
                    t_chunk_t, df_chunk_t[:, None]).cpu().numpy()
                p_chunk = None
                if compute_p_values:
                    with np.errstate(under="ignore"):
                        p_chunk = np.power(10.0, -logp_chunk)
                yield start, end, beta_chunk, t_chunk, p_chunk, logp_chunk
                continue
            p_chunk = (_two_sided_t_pvalue(t_chunk, df=df_chunk)
                       if compute_p_values else None)
            yield start, end, beta_chunk, t_chunk, p_chunk

    return apply_phenotype_missingness(_iterator()), q_matrix
