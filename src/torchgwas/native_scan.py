"""One CUDA association pipeline for decoded host or device-native inputs."""
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import math
import os
import sys
import time

import numpy as np
import torch

from .streaming import PinnedDosageLoader


# One pair of auxiliary streams per device, reused by every scan.
#
# `torch.cuda.Stream()` draws round-robin from a pool of 32 handles, so a fresh
# call per scan yields a *different* stream each time -- and the caching
# allocator keys cached blocks by stream, refusing to hand a block cached for
# one stream to another until a `cudaMalloc` actually fails. Each scan therefore
# reserved its own 16 x 890 MB staging pool that the next scan could not touch:
# measured as reserved memory climbing 19.5 -> 33.7 -> 48.0 -> 62.2 GB over four
# scans with allocated returning to 0.03 GB each time, 100% of the new segments
# inactive, `inactive_split` at 0.3% and `num_alloc_retries` zero -- growth by
# stream partitioning, not fragmentation. The fifth scan would reach 76.5 GB and
# the sixth would exhaust an 80 GB card.
#
# Keyed by device index because streams belong to a device.
_auxiliary_streams: dict[int, tuple] = {}


def _scan_streams(device):
    """The copy and result streams for `device`, created once and reused."""
    index = torch.device(device).index
    if index is None:
        index = torch.cuda.current_device()
    streams = _auxiliary_streams.get(index)
    if streams is None:
        streams = (torch.cuda.Stream(device=index), torch.cuda.Stream(device=index))
        _auxiliary_streams[index] = streams
    return streams


def dosage_cuda_iterator(source, phenotype, q_matrix, chunk_size, device,
                         reader_workers=None, prefetch_chunks=None,
                         compute_p_values=True, compute_log10_p=False,
                         variant_range=None,
                         reduction=None):
    from .linear import _dosage_statistics, _two_sided_t_pvalue
    from .tails import upper_tail_log10_from_t_torch

    profiling = os.environ.get('TORCHGWAS_SCAN_PROFILE', '0') != '0'
    blocking_events = os.environ.get('TORCHGWAS_BLOCKING_EVENTS', '0') != '0'
    from .scan_gpu import resolve_statistics_backend
    statistics_backend = resolve_statistics_backend()
    native_statistics = statistics_backend == 'native_fused'
    native_encoding = getattr(source, 'native_encoding', 'dosage')
    if native_encoding not in ('dosage', 'pgen_2bit'):
        raise ValueError(f"unsupported native transfer encoding: {native_encoding}")
    if native_encoding == 'pgen_2bit' and not native_statistics:
        raise ValueError("packed PGEN transfer requires TORCHGWAS_NATIVE_STATS=1")
    if native_statistics:
        from .scan_gpu import prepare as native_prepare, finish as native_finish
    setup_started = time.perf_counter() if profiling else 0
    timings = dict(enabled=profiling, blocking_completion_events=blocking_events, setup_seconds=0.0, fetch_seconds=0.0,
                   result_wait_seconds=0.0, copy_wait_seconds=0.0,
                   gpu_compute_milliseconds=0.0, gpu_input_conversion_milliseconds=0.0, gpu_result_milliseconds=0.0,
                   chunks=0)
    timings['statistics_backend'] = statistics_backend
    timings['native_encoding'] = native_encoding
    timings['conversion_fused_into_statistics'] = native_statistics
    source._last_scan_profile = timings

    device = torch.device(device)
    if device.index is None:
        device = torch.device("cuda", torch.cuda.current_device())
    torch.cuda.set_device(device)
    n, traits = phenotype.shape
    rank = 0 if q_matrix is None else q_matrix.shape[1]
    df = n - rank - 2
    # Offset rather than a df: each variant adds its own observed count.
    df_offset = float(-rank - 2)
    depth = max(2, int(prefetch_chunks or 3))
    reader_workers = getattr(source, "decode_workers", reader_workers)
    intercept = torch.full((n, 1), 1 / math.sqrt(n), device=device)
    covariates = intercept if q_matrix is None else torch.cat((intercept,
        torch.as_tensor(q_matrix, dtype=torch.float32, device=device)), dim=1)

    # Build the design ONCE and read the phenotype out of its leading
    # columns. The obvious form -- upload the phenotype, concatenate a design
    # beside it, then square it for the sums -- holds THREE full n x K
    # float32 matrices at the same moment. At 33,417 samples and 150,000
    # voxels that is 3 x 20.05 GB and it was the peak of the entire run:
    # 55.2 GB allocated where the scan's own rings need about 21 GB. At the
    # full 2,085,000 voxels one copy is 278 GB, so the difference between one
    # and three decides how many trait blocks a card can take, and therefore
    # how many times the whole genotype has to be re-read.
    design = torch.empty((n, traits + covariates.shape[1]),
                         dtype=torch.float32, device=device)
    design[:, traits:].copy_(covariates)
    phenotype_t = design[:, :traits]

    # Column-blocked: copying a contiguous host array into a STRIDED device
    # slice can stage a contiguous copy of the whole thing, and
    # `(A * A).sum(0)` allocates another n x K. Each is now one block wide.
    phenotype_ss = torch.empty(traits, dtype=torch.float32, device=device)
    column_block = max(1, min(traits, (1 << 28) // max(n * 4, 1)))
    for begin in range(0, traits, column_block):
        stop = min(begin + column_block, traits)
        columns = phenotype_t[:, begin:stop]
        columns.copy_(torch.as_tensor(
            np.ascontiguousarray(phenotype[:, begin:stop]), dtype=torch.float32))
        torch.sum(columns * columns, dim=0, out=phenotype_ss[begin:stop])
    copy_stream, result_stream = _scan_streams(device)
    compute_stream = torch.cuda.current_stream(device)
    selected_backend = (source.resolve_decode_backend(device)
                        if hasattr(source, "resolve_decode_backend") else
                        getattr(source, "decode_backend", "gpu"))
    device_source = hasattr(source, "iter_device_chunks") and selected_backend != "cpu"
    loader = None
    native_dtype = np.dtype(getattr(source, "native_transfer_dtype",
                                   getattr(source, "native_dtype", np.float32)))
    row_width = int(getattr(source, 'native_row_width', n))
    # Device-native decoders account for their compressed transfers themselves.
    timings['transfer_bytes_per_variant'] = 0 if device_source else row_width * native_dtype.itemsize
    # The geometry the scan actually ran with, not what the caller asked for.
    # Chunk and depth are usually chosen by the planner, so without recording
    # them here nothing downstream can check a predicted ring against the real
    # one -- which is the whole point of having a model that both predicts and
    # chooses.
    timings['chunk_variants'] = int(chunk_size)
    timings['depth'] = int(depth)
    timings['decode_on_gpu'] = bool(device_source)
    transfer_dtype = {np.dtype(np.int8): torch.int8, np.dtype(np.uint8): torch.uint8}.get(
        native_dtype, torch.float32)

    # These buffers are written by an H2D on copy_stream, so they are allocated
    # on copy_stream: the caching allocator only guarantees that reusing a freed
    # block is safe for work on the stream it was allocated on. Allocating them
    # on the compute stream instead let the allocator hand back the address of a
    # setup temporary whose producing kernel was still queued -- the (n, traits)
    # intermediate behind phenotype_ss is exactly such a temporary -- and that
    # kernel could then land on top of an already-copied chunk, silently
    # corrupting its first `traits` variants. Compute reads them only after
    # copy_done, and the final stream drain covers their release.
    with torch.cuda.stream(copy_stream):
        device_buffers = [] if device_source else [torch.empty(
            (chunk_size, row_width), dtype=transfer_dtype, device=device) for _ in range(depth)]
    # Belt and braces for the same hazard in any other setup allocation: no
    # transfer may overtake work already queued on the compute stream.
    copy_stream.wait_stream(compute_stream)
    result_stream.wait_stream(compute_stream)
    copy_done = [torch.cuda.Event(blocking=blocking_events) for _ in range(depth)]
    compute_done = [torch.cuda.Event(enable_timing=profiling) for _ in range(depth)]
    result_done = [torch.cuda.Event(enable_timing=profiling, blocking=blocking_events) for _ in range(depth)]
    compute_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if profiling else []
    conversion_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if profiling else []
    result_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if profiling else []
    # A reduced scan stages `chunk x k`, not `chunk x K`. That is the point of
    # the mode: at K = 10^5 each of these slots would otherwise be 205 MB of
    # pinned memory per buffer per slot, which is where a high-trait scan dies
    # before it ever reaches the writer.
    reduction_width = None if reduction is None else reduction.resolved_width(traits)
    if reduction is None:
        result_buffers = [(torch.empty((chunk_size, traits), pin_memory=True),
                           torch.empty((chunk_size, traits), pin_memory=True),
                           torch.empty(chunk_size, dtype=torch.uint8, pin_memory=True),
                           # One residual df per variant, not per test.
                           torch.empty(chunk_size, pin_memory=True),
                           *(tuple([torch.empty((chunk_size, traits),
                                                dtype=torch.float64,
                                                pin_memory=True)])
                             if compute_log10_p else tuple()))
                          for _ in range(depth)]
    else:
        result_buffers = [reduction.host_buffers(chunk_size, reduction_width)
                          for _ in range(depth)]
    pending = deque()
    release_pool = None if device_source else ThreadPoolExecutor(
        max_workers=1, thread_name_prefix='torchgwas-copy-release')
    release_futures = [None] * depth

    def release_after_copy(slot, buffer_index):
        copy_done[slot].synchronize()
        loader.release(buffer_index)

    pool = ThreadPoolExecutor(max_workers=depth, thread_name_prefix="torchgwas-result")
    source._last_scan_exclusion_counts = {"missing": 0, "invariant": 0}

    def finish(slot, start, end):
        result_done[slot].synchronize()
        count = end - start
        # Return owned arrays: callers may retain results beyond ring reuse.
        staged = [value[:count].numpy().copy() for value in result_buffers[slot]]
        logp = None
        if reduction is None:
            beta, t, status, variant_df = staged[:4]
            logp = staged[4] if compute_log10_p else None
            trait_index = None
        else:
            beta, t, trait_index, status, variant_df = staged
        malformed = np.flatnonzero(status == 3)
        if malformed.size:
            raise ValueError(f"invalid native ALT1 dosage at variant {start + int(malformed[0])}")
        invalid = status != 0
        beta[invalid] = np.nan
        t[invalid] = np.nan
        if logp is not None:
            logp[invalid] = np.nan
        # df is per variant either way; a reduced chunk broadcasts it across the
        # k traits it kept, so the p-values match the unreduced scan exactly.
        df_for_p = variant_df[:, None] if reduction is not None else variant_df
        if logp is not None and compute_p_values:
            with np.errstate(under="ignore"):
                p = np.power(10.0, -logp)
        else:
            p = (_two_sided_t_pvalue(t, df_for_p) if compute_p_values else None)
        gpu_times = (compute_start[slot].elapsed_time(compute_done[slot]),
                     result_start[slot].elapsed_time(result_done[slot]),
                     0.0 if device_source else conversion_start[slot].elapsed_time(compute_start[slot])) if profiling else None
        emitted = ((start, end, beta, t, p, logp)
                   if reduction is None and compute_log10_p else
                   (start, end, beta, t, p) if reduction is None
                   else (start, end, beta, t, p, trait_index))
        return emitted, int((status == 1).sum()), int((status == 2).sum()), gpu_times

    def resolve(future):
        started = time.perf_counter() if profiling else 0
        result, missing, invariant, gpu_times = future.result()
        if profiling:
            timings['result_wait_seconds'] += time.perf_counter() - started
            timings['gpu_compute_milliseconds'] += gpu_times[0]
            timings['gpu_result_milliseconds'] += gpu_times[1]
            timings['gpu_input_conversion_milliseconds'] += gpu_times[2]
            timings['chunks'] += 1
        source._last_scan_exclusion_counts["missing"] += missing
        source._last_scan_exclusion_counts["invariant"] += invariant
        return result

    try:
        # Both branches must carry the variant range. The pinned loader accepts
        # and honours it, but this call omitted it, so a ranged scan of a
        # non-device source silently read the **whole file** instead of the
        # range asked for -- no error, just the wrong extent, and a plausible
        # runtime. Same failure the multi-GPU duplicate-coverage check exists to
        # catch; here it would have had every shard scan every variant.
        loader = (source.iter_device_chunks(chunk_size=chunk_size, device=device,
                  reader_workers=reader_workers, prefetch_chunks=depth,
                  variant_range=variant_range)
                  if device_source else PinnedDosageLoader(source, chunk_size, depth,
                                                           reader_workers,
                                                           variant_range=variant_range))
        if profiling:
            timings['setup_seconds'] = time.perf_counter() - setup_started
        def timed_items():
            iterator = iter(loader)
            while True:
                started = time.perf_counter()
                try:
                    item = next(iterator)
                except StopIteration:
                    timings['fetch_seconds'] += time.perf_counter() - started
                    return
                timings['fetch_seconds'] += time.perf_counter() - started
                yield item
        for iteration, item in enumerate(timed_items() if profiling else loader):
            if len(pending) >= depth:
                yield resolve(pending.popleft())
            slot = iteration % depth
            if device_source:
                start, end, genotype_t = item
                if genotype_t.device != device or genotype_t.shape != (end - start, n):
                    raise ValueError("device decoder returned the wrong device or shape")
                genotype_t.record_stream(compute_stream)
            else:
                buffer_index, host, start, end = item
                # Complete the previous lease before recording this slot's event again.
                if release_futures[slot] is not None:
                    copy_wait_started = time.perf_counter() if profiling else 0
                    release_futures[slot].result()
                    if profiling:
                        timings['copy_wait_seconds'] += time.perf_counter() - copy_wait_started
                if iteration >= depth:
                    copy_stream.wait_event(compute_done[slot])
                with torch.cuda.stream(copy_stream):
                    device_buffers[slot][:end-start].copy_(host, non_blocking=True)
                    copy_done[slot].record(copy_stream)
                compute_stream.wait_event(copy_done[slot])
                if profiling:
                    conversion_start[slot].record(compute_stream)
                genotype_t = device_buffers[slot][:end-start]
                if not native_statistics and transfer_dtype == torch.int8:
                    genotype_t = torch.where(genotype_t == -9, torch.nan,
                                             genotype_t.to(torch.float32))
                elif not native_statistics and transfer_dtype == torch.uint8:
                    genotype_t = genotype_t.to(torch.float32) / float(source.native_scale)
                elif (not native_statistics and getattr(source, 'allows_direct_native_fill', False) and
                      getattr(source, 'native_missing_value', None) is not None):
                    genotype_t = torch.where(genotype_t == source.native_missing_value,
                                             torch.nan, genotype_t)
            if profiling:
                compute_start[slot].record(compute_stream)
            if native_statistics:
                scale = (float(source.native_scale)
                         if genotype_t.dtype == torch.uint8 and native_encoding == 'dosage' else 1.0)
                missing = (getattr(source, 'native_missing_value', None)
                           if getattr(source, 'allows_direct_native_fill', False) else None)
                if native_encoding == 'pgen_2bit':
                    centered, centered_ss, minimum, maximum, present = native_prepare(
                        genotype_t, encoding='pgen_2bit', n_samples=n)
                else:
                    centered, centered_ss, minimum, maximum, present = native_prepare(genotype_t, scale, missing)
                products = centered @ design
                beta, t, status = native_finish(
                    products, centered_ss, minimum, maximum, phenotype_ss,
                    present, df_offset,
                    getattr(source, 'validate_native_range', False))
                variant_df = present.to(torch.float32) + df_offset
            else:
                beta, t, status, variant_df = _dosage_statistics(
                    genotype_t, design, phenotype_ss, traits, df,
                    getattr(source, "validate_native_range", False),
                    covariate_rank=rank)
            if reduction is not None:
                # Reduce on the compute stream, before compute_done is recorded,
                # so the narrow result is what the copy stream waits for and the
                # wide (chunk x K) tensors are freed here rather than travelling.
                beta, t, index_t, status, variant_df = reduction.reduce(
                    beta, t, status, variant_df, reduction_width)
                staged_values = (beta, t, index_t, status, variant_df)
            else:
                logp = (upper_tail_log10_from_t_torch(
                    t, variant_df[:, None]) if compute_log10_p else None)
                staged_values = (beta, t, status, variant_df)
                if logp is not None:
                    staged_values += (logp,)
            compute_done[slot].record(compute_stream)
            with torch.cuda.stream(result_stream):
                result_stream.wait_event(compute_done[slot])
                if profiling:
                    result_start[slot].record(result_stream)
                for destination, value in zip(result_buffers[slot], staged_values):
                    destination[:end-start].copy_(value, non_blocking=True)
                    value.record_stream(result_stream)
                result_done[slot].record(result_stream)
            if not device_source:
                # Return pinned storage only after DMA, without serializing submission.
                release_futures[slot] = release_pool.submit(release_after_copy, slot, buffer_index)
            pending.append(pool.submit(finish, slot, start, end))
        while pending:
            yield resolve(pending.popleft())
    finally:
        # A CUDA error must not skip producer/reader shutdown. Preserve the primary
        # exception while still checking asynchronous errors in the final slots.
        primary_error = sys.exc_info()[0] is not None
        cleanup_errors = []
        def cleanup(operation):
            try:
                operation()
            except BaseException as error:
                cleanup_errors.append(error)
        for stream in (copy_stream, result_stream, compute_stream):
            cleanup(stream.synchronize)
        if release_pool is not None:
            cleanup(lambda: release_pool.shutdown(wait=True))
        if hasattr(loader, "close"):
            cleanup(loader.close)
        cleanup(lambda: pool.shutdown(wait=True, cancel_futures=True))
        for future in release_futures:
            if future is not None:
                cleanup(future.result)
        # Drop the large buffers here so their release is deterministic rather
        # than waiting on the garbage collector.
        #
        # Hygiene, not a fix for anything measured: `torch.cuda.memory_allocated`
        # already returns to 0.03 GB after a scan without this, so nothing was
        # being retained. Repeated scans in one process *do* degrade badly --
        # four identical ones took 48, 257 and 88 seconds while **reserved**
        # device memory climbed 33.7 -> 48.0 -> 62.2 GB -- but that is the
        # caching allocator fragmenting and requesting fresh segments, not this
        # frame holding references. See the repeat-scan item; the remedy is on
        # the allocation-pattern side, not here.
        cleanup(pending.clear)
        release_futures[:] = [None] * len(release_futures)
        del device_buffers[:]
        del result_buffers[:]
        if cleanup_errors and not primary_error:
            raise cleanup_errors[0]
