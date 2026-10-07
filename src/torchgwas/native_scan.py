"""One CUDA association pipeline for decoded host or device-native inputs."""
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import contextlib
import math
import os
import sys
import threading
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


# Chunk graphs (TORCHGWAS_CHUNK_GRAPHS): captured on one cached stream per
# device (a fresh stream per scan would partition the allocator, as above), one
# capture at a time across the shard threads, in thread-local mode so the other
# shards' reader threads may keep pinning and synchronizing meanwhile.
_capture_streams: dict[int, torch.cuda.Stream] = {}
_GRAPH_CAPTURE_LOCK = threading.Lock()
# Chunk graphs per GPU, used in turn: one computes while the other's results copy out.
GRAPH_LANES = 2


def _capture_stream(device):
    index = torch.device(device).index
    stream = _capture_streams.get(index)
    if stream is None:
        stream = _capture_streams[index] = torch.cuda.Stream(device=index)
    return stream


def resolve_reader_workers(source, requested=None, limit=None):
    """Preserve source preference while enforcing an aggregate driver's cap."""
    workers = getattr(source, "decode_workers", requested)
    if limit is not None:
        if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
            raise ValueError("reader worker limit must be a positive integer")
        workers = limit if workers is None else min(int(workers), limit)
    return workers

def dosage_cuda_iterator(source, phenotype, q_matrix, chunk_size, device,
                         reader_workers=None, prefetch_chunks=None,
                         compute_p_values=True, variant_range=None,
                         reduction=None, reader_worker_limit=None,
                         borrow_results=False, return_df=False, significance=None, significance_n_traits=None,
                         chunk_size_selector=None, chunk_observer=None, return_beta=True,
                         shared_loader=None, log10_p=None, complete_case=None):
    """shared_loader: a SharedDecodeHub subscriber used instead of this scan's
    own PinnedDosageLoader, so several tile scans share one decode pass.

    complete_case (complete_case.CompleteCasePlan, missing phenotypes): traits
    with missing values get complete-case OLS on the device; unreduced chunks
    stage the (chunk, K) pair df, which is their df."""
    from .linear import (_device_log10_p, _dosage_statistics, _log10_p_dtype, _two_sided_t_pvalue,
                         _unpack_pgen_2bit_float)
    logp_dtype = _log10_p_dtype(log10_p)

    if type(return_beta) is not bool or (not return_beta and reduction is not None):
        raise ValueError('Beta omission requires an unreduced scan and a boolean return_beta')

    measurement_window = None
    if chunk_size_selector is not None or chunk_observer is not None:
        from .adaptive_chunks import (validate_chunk_control, _ChunkDelivery,
                                      InitialChunkMeasurements, ChunkDeviceTiming)
        validate_chunk_control(source, torch.device(device), 'float32', chunk_size,
                               chunk_size_selector, chunk_observer)
        if isinstance(chunk_observer, InitialChunkMeasurements):
            measurement_window = chunk_observer
    measure_cuda = measurement_window is not None and measurement_window.record_cuda
    profiling = os.environ.get('TORCHGWAS_SCAN_PROFILE', '0') != '0'
    timed_gpu = profiling or measure_cuda
    blocking_events = os.environ.get('TORCHGWAS_BLOCKING_EVENTS', '0') != '0'
    from .scan_gpu import fused_module, resolve_statistics_backend
    statistics_backend = resolve_statistics_backend(device)
    # The CUDA kernels or their Triton port: the same prepare/finish contract.
    fused = fused_module(statistics_backend)
    native_statistics = fused is not None
    native_encoding = getattr(source, 'native_encoding', 'dosage')
    if native_encoding not in ('dosage', 'pgen_2bit'):
        raise ValueError(f"unsupported native transfer encoding: {native_encoding}")
    # Packed rows on a device without fused statistics (Triton failed there):
    # unpacked by Torch ops, then the Torch statistics.
    unpack_packed = native_encoding == 'pgen_2bit' and not native_statistics
    if native_statistics:
        native_prepare, native_finish = fused.prepare, fused.finish
    # min-p on a complete panel ranks by |t| inside the Triton finish, so the
    # chunk's (variants x traits) beta and t are never written.
    fused_min_p = (getattr(fused, 'finish_min_p', None) if log10_p is not None and complete_case is None
                   and getattr(reduction, 'mode', None) == 'min-p' and hasattr(reduction, 'from_winners')
                   and reduction.resolved_width(phenotype.shape[1]) == 1
                   else None)
    setup_started = time.perf_counter() if profiling else 0
    timings = dict(enabled=profiling, blocking_completion_events=blocking_events, setup_seconds=0.0, fetch_seconds=0.0,
                   result_wait_seconds=0.0, copy_wait_seconds=0.0,
                   gpu_compute_milliseconds=0.0, gpu_input_conversion_milliseconds=0.0, gpu_result_milliseconds=0.0,
                   chunks=0)
    timings['statistics_backend'] = statistics_backend
    timings['native_encoding'] = native_encoding
    timings['return_beta'] = return_beta
    timings['result_payload_bytes'] = 0
    timings['conversion_fused_into_statistics'] = native_statistics
    timings['packed_unpacked_by_torch'] = unpack_packed
    source._last_scan_profile = timings

    device = torch.device(device)
    if device.index is None:
        device = torch.device("cuda", torch.cuda.current_device())
    torch.cuda.set_device(device)
    if measurement_window is not None and str(device) not in measurement_window.devices:
        raise ValueError('Measurement device is outside the configured window')
    n, traits = phenotype.shape
    rank = 0 if q_matrix is None else q_matrix.shape[1]
    df = n - rank - 2
    # Offset rather than a df: each variant adds its own observed count.
    df_offset = float(-rank - 2)
    if significance is not None and reduction is not None:
        raise ValueError('Significance and joint reduction are mutually exclusive')
    if significance is not None:
        from .reduce import device_significance_critical, device_significant_pairs
        critical=device_significance_critical(significance,n,significance_n_traits or traits,device)
        timings['result_selection']='device_significant'
        timings['result_payload_bytes']=0
    depth = max(2, int(prefetch_chunks or 3))
    reader_workers = resolve_reader_workers(source, reader_workers, reader_worker_limit)
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
    #
    # A row that carries every sample in file order (a sample selection on
    # PGEN, e.g. subjects dropped for a missing phenotype; pgen.py): the
    # unselected samples are marked missing on the device, so they leave every
    # genotype sum and the observed count, and each design row sits at its
    # sample's file position. Rows at unselected positions stay zero.
    physical = getattr(source, 'native_physical_samples', None)
    positions = packed_missing = unselected = None
    if physical is not None:
        positions = torch.as_tensor(source.native_sample_positions, dtype=torch.long, device=device)
        if native_encoding == 'pgen_2bit':
            packed_missing = torch.as_tensor(source.native_packed_missing_mask, dtype=torch.uint8, device=device)
        else:
            # One-byte-per-sample rows: the unselected columns take the
            # transport's missing sentinel, which every path already masks.
            unselected = torch.as_tensor(source.native_unselected_samples, dtype=torch.long, device=device)
            unselected_value = source.native_missing_value
    rows = n if physical is None else int(physical)
    if physical is not None and complete_case is not None and getattr(complete_case, 'needs_calls', True):
        complete_case = complete_case.at_positions(source.native_sample_positions)
    timings['physical_samples'] = rows
    design = (torch.empty if positions is None else torch.zeros)(
        (rows, traits + covariates.shape[1]), dtype=torch.float32, device=device)
    if positions is None:
        design[:, traits:].copy_(covariates)
    else:
        design[positions, traits:] = covariates
    phenotype_t = design[:, :traits]

    # Column-blocked: copying a contiguous host array into a STRIDED device
    # slice can stage a contiguous copy of the whole thing, and
    # `(A * A).sum(0)` allocates another n x K. Each is now one block wide.
    phenotype_ss = torch.empty(traits, dtype=torch.float32, device=device)
    column_block = design_column_block(n, traits)
    for begin in range(0, traits, column_block):
        stop = min(begin + column_block, traits)
        if positions is not None:
            columns = (phenotype[:, begin:stop].to(device=device, dtype=torch.float32)
                       if isinstance(phenotype, torch.Tensor) else torch.as_tensor(
                           np.ascontiguousarray(phenotype[:, begin:stop]), dtype=torch.float32, device=device))
            design[positions, begin:stop] = columns
            torch.sum(columns * columns, dim=0, out=phenotype_ss[begin:stop])
            continue
        columns = phenotype_t[:, begin:stop]
        if isinstance(phenotype, torch.Tensor):
            columns.copy_(phenotype[:, begin:stop])  # device to device (variant shards)
        else:
            columns.copy_(torch.as_tensor(
                np.ascontiguousarray(phenotype[:, begin:stop]), dtype=torch.float32))
        torch.sum(columns * columns, dim=0, out=phenotype_ss[begin:stop])
    copy_stream, result_stream = _scan_streams(device)
    compute_stream = torch.cuda.current_stream(device)
    selected_backend = (source.resolve_decode_backend(device)
                        if hasattr(source, "resolve_decode_backend") else
                        getattr(source, "decode_backend", "gpu"))
    device_source = hasattr(source, "iter_device_chunks") and selected_backend != "cpu"
    if shared_loader is not None and (device_source or shared_loader.capacity != int(chunk_size)):
        raise ValueError('A shared decode loader needs host decoding and the same chunk capacity')
    timings['shared_decode'] = shared_loader is not None
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
    timings['chunk_variants'] = int(chunk_size)  # Allocated capacity if adaptive.
    timings['adaptive_chunk_sizes'] = chunk_size_selector is not None
    timings['chunk_observer_enabled'] = chunk_observer is not None
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
    copy_done = [torch.cuda.Event(enable_timing=measure_cuda, blocking=blocking_events) for _ in range(depth)]
    compute_done = [torch.cuda.Event(enable_timing=timed_gpu) for _ in range(depth)]
    result_done = [torch.cuda.Event(enable_timing=timed_gpu, blocking=blocking_events) for _ in range(depth)]
    copy_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if measure_cuda else []
    compute_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if timed_gpu else []
    conversion_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if timed_gpu else []
    result_start = [torch.cuda.Event(enable_timing=True) for _ in range(depth)] if timed_gpu else []
    # A reduced scan stages `chunk x k`, not `chunk x K`. That is the point of
    # the mode: at K = 10^5 each of these slots would otherwise be 205 MB of
    # pinned memory per buffer per slot, which is where a high-trait scan dies
    # before it ever reaches the writer.
    reduction_width = None if reduction is None else reduction.resolved_width(traits)
    # The complete-case pair df crosses to the host only when read there.
    stage_pair_df = complete_case is not None and (return_df or (compute_p_values and log10_p is None))
    # The staged results, as (per-row shape, dtype), None where absent.
    if significance is not None:
        result_specs = None
    elif reduction is None:
        result_specs = (((traits,), torch.float32) if return_beta else None,
                        ((traits,), torch.float32),
                        ((), torch.uint8),
                        # One residual df per variant, not per test.
                        ((), torch.float32),
                        # float32 -log10 P, computed on the device (log10_p).
                        *((((traits,), logp_dtype),) if log10_p is not None else ()),
                        # Complete-case pair df (missing phenotypes).
                        *((((traits,), torch.float32),) if stage_pair_df else ()))
    else:
        templates = (reduction.host_buffers(1, reduction_width, pin_memory=False) if log10_p is None else
                     reduction.host_buffers(1, reduction_width, pin_memory=False, log10_p_dtype=logp_dtype))
        result_specs = tuple((tuple(template.shape[1:]), template.dtype) for template in templates)

    def result_layout(rows):
        """Each staged array's byte offset for a chunk of `rows`, 64-byte aligned, and the total."""
        offsets, total = [], 0
        for spec in result_specs:
            if spec is None:
                offsets.append(None)
                continue
            offsets.append(total)
            total += -(-rows * math.prod(spec[0]) * torch.empty((), dtype=spec[1]).element_size() // 64) * 64
        return offsets, total

    def carve(block, rows):
        """The staged arrays of a `rows`-variant chunk, laid out back to back in `block` (bytes)."""
        offsets, _ = result_layout(rows)
        views = []
        for spec, offset in zip(result_specs, offsets):
            if spec is None:
                views.append(None)
                continue
            size = rows * math.prod(spec[0]) * torch.empty((), dtype=spec[1]).element_size()
            views.append(block[offset:offset + size].view(spec[1]).view((rows,) + spec[0]))
        return tuple(views)

    # One pinned block per slot, pinned on the slot's first chunk and grown
    # once to the capacity, like the input ring (streaming.PinnedDosageLoader);
    # dense output at K = 512 and depth 32 pins 537 MB per GPU at chunk 4096.
    # A chunk's arrays are carved for its own row count, so they are contiguous
    # and a chunk graph moves them all with one copy (chunk_graph).
    result_blocks = [] if result_specs is None else [None] * depth
    result_rows = [0] * depth
    result_views = [dict() for _ in range(depth)]

    def result_slot(slot, rows):
        """This slot's pinned results for a `rows`-variant chunk. Its previous chunk was resolved
        before this iteration, so they are free."""
        if result_blocks[slot] is None or result_rows[slot] < rows:
            size = rows if result_blocks[slot] is None else chunk_size
            result_blocks[slot] = torch.empty(result_layout(size)[1], dtype=torch.uint8, pin_memory=True)
            result_rows[slot] = size
            result_views[slot].clear()
        views = result_views[slot].get(rows)
        if views is None:
            views = result_views[slot][rows] = carve(result_blocks[slot], rows)
        return views

    # Chunk graphs: a chunk's device work -- prepare, GEMM, finish (or min-p's
    # ranking and tail), packing the results into one buffer -- replays as one
    # CUDA graph, and its results cross to the host in one copy: three host
    # calls with the input's device copy, instead of ~25 from Python. With
    # variant shards those calls queue on the GIL: at K = 512 on four GPUs
    # each GPU was busy a third of the scan. Two graphs per GPU, used in turn
    # (GRAPH_LANES), each on its own input and output buffers, so a chunk can
    # compute while the one before is still copying out; a ring slot is copied
    # into a lane's input on the device. A graph per slot instead was captured
    # depth times per GPU, one capture at a time across shards while the rest
    # of the process held the GIL: four A100 shards at depth 4 spent 1.5 s,
    # summed over the shards, capturing their 16 slots.
    # Captured lazily, after one eager chunk has compiled and warmed
    # everything, for the steady chunk size only: one the chunk before also
    # had, and under the chunk-size tuner only once it has settled (its trials
    # stay eager, so they compare like with like). A short last chunk, and
    # every mode outside native statistics with dense output or fused min-p,
    # stay eager.
    graph_subset = (native_statistics and not device_source and significance is None
                    and complete_case is None and not timed_gpu
                    and (reduction is None or fused_min_p is not None))
    graphs_enabled = graph_subset and os.environ.get('TORCHGWAS_CHUNK_GRAPHS', '0') != '0'
    timings['chunk_graphs'] = graphs_enabled
    graphs = {}
    graph_inputs, graph_outputs = [None] * GRAPH_LANES, [None] * GRAPH_LANES
    lane_done = [torch.cuda.Event() for _ in range(GRAPH_LANES)] if graphs_enabled else []
    graph_state = dict(enabled=graphs_enabled, warmed=False, previous=None, next_lane=0,
                       pool=torch.cuda.graph_pool_handle() if graphs_enabled else None)
    validate_range = getattr(source, 'validate_native_range', False)
    graph_missing = (getattr(source, 'native_missing_value', None)
                     if getattr(source, 'allows_direct_native_fill', False) else None)

    def native_compute(genotype_t):
        """The eager loop's native-statistics work for dense output or fused min-p, as one function."""
        if packed_missing is not None:
            genotype_t.bitwise_or_(packed_missing)
        elif unselected is not None:
            genotype_t.index_fill_(1, unselected, unselected_value)
        if native_encoding == 'pgen_2bit':
            prepared = native_prepare(genotype_t, encoding='pgen_2bit', n_samples=rows)
        else:
            scale = (float(source.native_scale)
                     if genotype_t.dtype == torch.uint8 and native_encoding == 'dosage' else 1.0)
            prepared = native_prepare(genotype_t, scale, graph_missing)
        centered, centered_ss, minimum, maximum, present = prepared
        products = centered @ design
        variant_df = present.to(torch.float32) + df_offset
        if fused_min_p is not None:
            winners = fused_min_p(products, centered_ss, minimum, maximum, phenotype_ss, present,
                                  df_offset, validate_range)
            return reduction.from_winners(*winners, variant_df, logp_dtype)
        beta, t, status = native_finish(products, centered_ss, minimum, maximum, phenotype_ss,
                                        present, df_offset, validate_range)
        staged = (beta if return_beta else None, t, status, variant_df)
        return staged + ((_device_log10_p(t, variant_df, logp_dtype),) if log10_p is not None else ())

    def chunk_graph(count):
        """(lane, graph, payload bytes) for a steady-size chunk, or None to run it eagerly."""
        steady, graph_state['previous'] = count == graph_state['previous'], count
        if not (steady and graph_state['enabled'] and graph_state['warmed']):
            return None
        # Under the chunk-size tuner only once it has settled; any other
        # observer (a measurement window) keeps the eager chunks it measures.
        state = 'committed' if chunk_observer is None else getattr(chunk_observer, 'state', None)
        if state not in ('committed', 'fixed', 'skipped'):
            return None
        lane = graph_state['next_lane']
        captured = graphs.get((lane, count))
        if captured is None:
            captured = capture(lane, count)
            if captured is None:
                return None
        graph_state['next_lane'] = (lane + 1) % GRAPH_LANES
        return (lane,) + captured

    def capture(lane, count):
        capture_started = time.perf_counter()
        try:
            if graph_outputs[lane] is None:
                # Beside one chunk's intermediates, which the graph pool holds
                # again (the eager chunk's stay cached): leave room for both.
                need = chunk_size * (rows * 4 + design.shape[1] * 4 + traits * 16 + 64)
                if 2 * need > torch.cuda.mem_get_info(device)[0]:
                    raise MemoryError(f'{2 * need / 2**30:.1f} GiB free needed for chunk graphs')
                # Read and written on the compute stream only, so allocated there.
                graph_inputs[lane] = torch.empty((chunk_size, row_width), dtype=transfer_dtype, device=device)
                graph_outputs[lane] = torch.empty(result_layout(chunk_size)[1], dtype=torch.uint8, device=device)
            side = _capture_stream(device)
            side.wait_stream(compute_stream)
            graph = torch.cuda.CUDAGraph()
            with _GRAPH_CAPTURE_LOCK, torch.cuda.stream(side):
                timings['chunk_graph_lock_wait_seconds'] = (timings.get('chunk_graph_lock_wait_seconds', 0.0)
                                                            + time.perf_counter() - capture_started)
                # cuBLAS keeps a workspace per (handle, stream): create this
                # stream's outside the capture, not inside the graph's pool.
                torch.mm(design[:16], design[:16].T)
                graph.capture_begin(pool=graph_state['pool'], capture_error_mode='thread_local')
                try:
                    outputs = native_compute(graph_inputs[lane][:count])
                    for packed, value in zip(carve(graph_outputs[lane], count), outputs):
                        if packed is not None:
                            packed.copy_(value)
                finally:
                    graph.capture_end()
            compute_stream.wait_stream(side)
        except Exception as error:  # noqa: BLE001 - the eager path remains
            timings['chunk_graphs'] = graph_state['enabled'] = False
            timings['chunk_graph_error'] = f'{type(error).__name__}: {error}'[:500]
            graphs.clear()
            return None
        # The outputs' memory stays the graph's (its private pool); the
        # packed buffer is what the host copies.
        payload = result_layout(count)[1]
        graphs[lane, count] = (graph, payload)
        timings['chunk_graphs_captured'] = len(graphs)
        timings['chunk_graph_capture_seconds'] = (timings.get('chunk_graph_capture_seconds', 0.0)
                                                  + time.perf_counter() - capture_started)
        return graphs[lane, count]
    pending = deque()
    release_pool = None if device_source else ThreadPoolExecutor(
        max_workers=1, thread_name_prefix='torchgwas-copy-release')
    release_futures = [None] * depth

    def release_after_copy(slot, buffer_index):
        copy_done[slot].synchronize()
        loader.release(buffer_index)

    pool = ThreadPoolExecutor(max_workers=depth, thread_name_prefix="torchgwas-result")
    source._last_scan_exclusion_counts = {"missing": 0, "invariant": 0}
    # Significant pairs: select on a separate stream, one chunk behind, so the
    # host collects chunk i's pairs while chunk i+1 already runs.
    select_stream = torch.cuda.Stream(device) if significance is not None else None
    pending_selection = deque()
    selection_lag = 0
    if significance is not None:
        lag_setting = os.environ.get('TORCHGWAS_SELECTION_LAG', 'auto')
        if lag_setting == 'auto':
            free_bytes, _ = torch.cuda.mem_get_info(device)
            # beta, t and the product of one extra chunk stay resident.
            selection_lag = int(3 * chunk_size * traits * 4 < 0.2 * free_bytes)
        else:
            selection_lag = max(0, min(1, int(lag_setting)))
        timings['selection_lag'] = selection_lag

    def select(entry):
        """Collect one chunk's significant pairs (host side), then yield them."""
        slot, start, end, beta, t, status, variant_df, delivery = entry
        with torch.cuda.stream(select_stream):
            select_stream.wait_event(compute_done[slot])
            status_host = status.cpu().numpy()
            malformed = np.flatnonzero(status_host == 3)
            if malformed.size:
                raise ValueError(f"invalid native ALT1 dosage at variant {start+int(malformed[0])}")
            source._last_scan_exclusion_counts['missing'] += int((status_host == 1).sum())
            source._last_scan_exclusion_counts['invariant'] += int((status_host == 2).sum())
            timings['result_payload_bytes'] += status_host.nbytes
            # Gather before yielding: the stream context must not stay active
            # while the consumer runs.
            blocks = list(device_significant_pairs(beta, t, status, variant_df, critical, start=start))
        del beta, t, status, variant_df
        for selected in blocks:
            timings['result_payload_bytes'] += sum(a.nbytes for a in selected[2:])
            if delivery is None:
                yield selected
            else:
                yield from delivery.deliver(selected)
        if delivery is not None:
            if delivery.record_cuda:
                delivery.cuda = device_timing(slot, result_transfer=False)
            delivery.complete()
        if profiling:
            compute_done[slot].synchronize()
            timings['gpu_compute_milliseconds'] += compute_start[slot].elapsed_time(compute_done[slot])
            timings['chunks'] += 1

    def device_timing(slot, *, result_transfer):
        return ChunkDeviceTiming(
            copy_start[slot].elapsed_time(copy_done[slot]) / 1000.,
            conversion_start[slot].elapsed_time(compute_start[slot]) / 1000.,
            compute_start[slot].elapsed_time(compute_done[slot]) / 1000.,
            result_start[slot].elapsed_time(result_done[slot]) / 1000. if result_transfer else None)

    def finish(slot, start, end, delivery=None):
        result_done[slot].synchronize()
        if delivery is not None and delivery.record_cuda:
            delivery.cuda = device_timing(slot, result_transfer=True)
        count = end - start
        # The loop yields this slot before reusing it. Borrowed views are valid
        # until the caller requests another chunk; owned results may be kept.
        views = result_views[slot][count]  # carved for `count` rows (result_slot)
        if borrow_results:
            staged = [None if value is None else value[:count].numpy() for value in views]
        else:
            staged = [None if value is None else value[:count].numpy().copy() for value in views]
        logp = None
        pair_df = None
        if reduction is None:
            beta, t, status, variant_df, *extra = staged
            logp = extra.pop(0) if log10_p is not None else None
            pair_df = extra.pop(0) if stage_pair_df else None
            trait_index = None
        else:
            beta, t, trait_index, status, variant_df, *kept_pairs = staged
            if kept_pairs:
                # The kept pairs' -log10 P and df, from the device (VariantReduction.reduce).
                logp, pair_df = kept_pairs
        malformed = np.flatnonzero(status == 3)
        if malformed.size:
            raise ValueError(f"invalid native ALT1 dosage at variant {start + int(malformed[0])}")
        invalid = status != 0
        if beta is not None:beta[invalid] = np.nan
        t[invalid] = np.nan
        if logp is not None:logp[invalid] = np.nan
        # df is per variant, or per pair for complete-case traits; a reduced
        # chunk broadcasts it across the k traits it kept, so the p-values
        # match the unreduced scan exactly.
        df_for_p = (variant_df[:, None] if reduction is not None
                    else variant_df if pair_df is None else pair_df)
        if logp is not None:
            from .linear import _p_from_log10
            p = _p_from_log10(logp) if compute_p_values else None
        else:
            p = (_two_sided_t_pvalue(t, df_for_p) if compute_p_values else None)
        gpu_times = (compute_start[slot].elapsed_time(compute_done[slot]),
                     result_start[slot].elapsed_time(result_done[slot]),
                     0.0 if device_source else conversion_start[slot].elapsed_time(compute_start[slot])) if profiling else None
        emitted = ((start, end, beta, t, p, *((logp,) if logp is not None else ())) if reduction is None
                   else (start, end, beta, t, p, trait_index, *((logp, pair_df) if logp is not None else ())))
        if return_df:
            emitted = (*emitted, variant_df[:, None] if pair_df is None or reduction is not None else pair_df)
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

    def emit_observed(item):
        future, delivery = item
        result = resolve(future)
        if delivery is None:
            yield result
        else:
            yield from delivery.deliver(result)
            delivery.complete()

    try:
        # Both branches must carry the variant range. The pinned loader accepts
        # and honours it, but this call omitted it, so a ranged scan of a
        # non-device source silently read the **whole file** instead of the
        # range asked for -- no error, just the wrong extent, and a plausible
        # runtime. Same failure the multi-GPU duplicate-coverage check exists to
        # catch; here it would have had every shard scan every variant.
        loader = (shared_loader if shared_loader is not None else
                  source.iter_device_chunks(chunk_size=chunk_size, device=device,
                  reader_workers=reader_workers, prefetch_chunks=depth,
                  variant_range=variant_range,
                  **({} if chunk_size_selector is None else {'chunk_size_selector': chunk_size_selector}))
                  if device_source else PinnedDosageLoader(source, chunk_size, depth,
                                                           reader_workers,
                                                           variant_range=variant_range,
                                                           chunk_size_selector=chunk_size_selector,
                                                           record_timing=(lambda start, end: measurement_window.reserve_read(start, end, str(device)))
                                                               if measurement_window is not None else chunk_observer is not None))
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
                if chunk_observer is None:
                    yield resolve(pending.popleft())
                else:
                    yield from emit_observed(pending.popleft())
            slot = iteration % depth
            delivery = None
            if device_source:
                start, end, genotype_t = item
                if genotype_t.device != device or genotype_t.shape != (end - start, n):
                    raise ValueError("device decoder returned the wrong device or shape")
                genotype_t.record_stream(compute_stream)
                if chunk_observer is not None:
                    from .adaptive_chunks import MinimalDelivery
                    delivery = MinimalDelivery(start, end, chunk_size, device, chunk_observer)
            else:
                buffer_index, host, start, end = item[:4]
                if chunk_observer is not None and len(item) > 4 and item[4] is not None:
                    delivery = _ChunkDelivery(start, end, chunk_size, device,
                                              item[4], chunk_observer, record_cuda=measure_cuda)
                # Complete the previous lease before recording this slot's event again.
                if release_futures[slot] is not None:
                    copy_wait_started = time.perf_counter() if profiling else 0
                    release_futures[slot].result()
                    if profiling:
                        timings['copy_wait_seconds'] += time.perf_counter() - copy_wait_started
                if iteration >= depth:
                    copy_stream.wait_event(compute_done[slot])
                # Shared decode with GPU fan-out: `host` is the chunk already in
                # the root GPU's ring. Wait for that copy, then pull it GPU to GPU
                # (NVLink peer copy). PyTorch runs a cross-device copy on the
                # source GPU's current stream, so give it a dedicated one there.
                fanout_ready = item[5] if len(item) > 5 else None
                if fanout_ready is not None:
                    copy_stream.wait_event(fanout_ready)
                peer = (loader.peer_stream(device) if fanout_ready is not None and host.device != device
                        else None)
                with (torch.cuda.stream(peer) if peer is not None else contextlib.nullcontext()):
                    with torch.cuda.stream(copy_stream):
                        if delivery is not None and delivery.record_cuda:
                            copy_start[slot].record(copy_stream)
                        device_buffers[slot][:end-start].copy_(host, non_blocking=True)
                        copy_done[slot].record(copy_stream)
                compute_stream.wait_event(copy_done[slot])
                graph = (chunk_graph(end - start)
                         if graph_state['enabled'] and (delivery is None or not delivery.record_cuda) else None)
                if graph is not None:
                    lane, compute_graph, payload = graph
                    # The lane's packed results are rewritten: the copy-out of
                    # its previous chunk must be done first.
                    compute_stream.wait_event(lane_done[lane])
                    graph_inputs[lane][:end - start].copy_(device_buffers[slot][:end - start])
                    compute_graph.replay()
                    compute_done[slot].record(compute_stream)
                    result_slot(slot, end - start)  # sizes the block, carves the arrays finish reads
                    block = result_blocks[slot]
                    with torch.cuda.stream(result_stream):
                        result_stream.wait_event(compute_done[slot])
                        block[:payload].copy_(graph_outputs[lane][:payload], non_blocking=True)
                        result_done[slot].record(result_stream)
                        lane_done[lane].record(result_stream)
                    timings['result_payload_bytes'] += payload
                    timings['chunk_graph_replays'] = timings.get('chunk_graph_replays', 0) + 1
                    release_futures[slot] = release_pool.submit(release_after_copy, slot, buffer_index)
                    future = pool.submit(finish, slot, start, end, delivery)
                    pending.append(future if chunk_observer is None else (future, delivery))
                    continue
                if profiling or (delivery is not None and delivery.record_cuda):
                    conversion_start[slot].record(compute_stream)
                genotype_t = device_buffers[slot][:end-start]
                if packed_missing is not None:
                    genotype_t.bitwise_or_(packed_missing)
                elif unselected is not None:
                    genotype_t.index_fill_(1, unselected, unselected_value)
                if unpack_packed:
                    genotype_t = _unpack_pgen_2bit_float(genotype_t, rows)
                elif not native_statistics and transfer_dtype == torch.int8:
                    genotype_t = torch.where(genotype_t == -9, torch.nan,
                                             genotype_t.to(torch.float32))
                elif not native_statistics and transfer_dtype == torch.uint8:
                    raw = genotype_t
                    genotype_t = raw.to(torch.float32) / float(source.native_scale)
                    # PGEN's uint8 dosage keeps 255 free as its missing code
                    # (pgen.py); whole rows put it at every unselected sample.
                    sentinel = getattr(source, 'native_missing_value', None)
                    if sentinel is not None:
                        genotype_t = torch.where(raw == sentinel, torch.nan, genotype_t)
                elif (not native_statistics and getattr(source, 'allows_direct_native_fill', False) and
                      getattr(source, 'native_missing_value', None) is not None):
                    genotype_t = torch.where(genotype_t == source.native_missing_value,
                                             torch.nan, genotype_t)
            if profiling or (delivery is not None and delivery.record_cuda):
                compute_start[slot].record(compute_stream)
            if native_statistics:
                scale = (float(source.native_scale)
                         if genotype_t.dtype == torch.uint8 and native_encoding == 'dosage' else 1.0)
                missing = (getattr(source, 'native_missing_value', None)
                           if getattr(source, 'allows_direct_native_fill', False) else None)
                if native_encoding == 'pgen_2bit':
                    centered, centered_ss, minimum, maximum, present = native_prepare(
                        genotype_t, encoding='pgen_2bit', n_samples=rows)
                else:
                    centered, centered_ss, minimum, maximum, present = native_prepare(genotype_t, scale, missing)
                products = centered @ design
                variant_df = present.to(torch.float32) + df_offset
                if fused_min_p is not None:
                    winners = fused_min_p(products, centered_ss, minimum, maximum, phenotype_ss, present,
                                          df_offset, getattr(source, 'validate_native_range', False))
                    beta = t = status = None
                else:
                    beta, t, status = native_finish(
                        products, centered_ss, minimum, maximum, phenotype_ss,
                        present, df_offset,
                        getattr(source, 'validate_native_range', False))
                pair_df = None
                if complete_case is not None:
                    from .complete_case import call_mask
                    pair_df = complete_case.correct(
                        centered, centered_ss, products[:, traits:], products[:, :traits], phenotype_ss,
                        variant_df, beta, t,
                        calls_observed=call_mask(genotype_t, missing_value=missing,
                                                 encoding='pgen_2bit' if native_encoding == 'pgen_2bit' else None))
            else:
                beta, t, status, variant_df, *extra = _dosage_statistics(
                    genotype_t, design, phenotype_ss, traits, df,
                    getattr(source, "validate_native_range", False),
                    covariate_rank=rank, complete_case=complete_case)
                pair_df = extra[0] if extra else None
            if significance is not None:
                # Release the input lease after its own H2D event, independently
                # of how long selection and the downstream writer take.
                if not device_source:
                    release_futures[slot] = release_pool.submit(release_after_copy, slot, buffer_index)
                compute_done[slot].record(compute_stream)
                for tensor in (beta, t, status, variant_df):
                    tensor.record_stream(select_stream)  # read there, allocated here
                pending_selection.append((slot, start, end, beta, t, status, variant_df, delivery))
                del beta, t, status, variant_df
                while len(pending_selection) > selection_lag:
                    yield from select(pending_selection.popleft())
                continue
            if reduction is not None:
                # Reduce on the compute stream, before compute_done is recorded,
                # so the narrow result is what the copy stream waits for and the
                # wide (chunk x K) tensors are freed here rather than travelling.
                if fused_min_p is not None:
                    staged_values = reduction.from_winners(*winners, variant_df, logp_dtype)
                    del winners
                else:
                    from .linear import _joint_pair_df
                    staged_values = reduction.reduce(
                        beta, t, status, variant_df, reduction_width,
                        **({} if log10_p is None else dict(log10_p=(pair_df, logp_dtype))),
                        **_joint_pair_df(reduction, pair_df))
                del beta, t, pair_df
            else:
                staged_values = (beta if return_beta else None, t, status, variant_df)
                if log10_p is not None:
                    staged_values += (_device_log10_p(t, variant_df if pair_df is None else pair_df, logp_dtype),)
                if stage_pair_df:
                    staged_values += (pair_df,)
            timings['result_payload_bytes'] += sum(value.numel()*value.element_size() for value in staged_values if value is not None)
            compute_done[slot].record(compute_stream)
            with torch.cuda.stream(result_stream):
                result_stream.wait_event(compute_done[slot])
                if profiling or (delivery is not None and delivery.record_cuda):
                    result_start[slot].record(result_stream)
                for destination, value in zip(result_slot(slot, end - start), staged_values):
                    if destination is None:continue
                    destination[:end-start].copy_(value, non_blocking=True)
                    value.record_stream(result_stream)
                result_done[slot].record(result_stream)
            if not device_source:
                # Return pinned storage only after DMA, without serializing submission.
                release_futures[slot] = release_pool.submit(release_after_copy, slot, buffer_index)
            future = pool.submit(finish, slot, start, end, delivery)
            pending.append(future if chunk_observer is None else (future, delivery))
            # One eager chunk has compiled the kernels and probed the tail.
            graph_state['warmed'] = graph_state['enabled']
        while pending_selection:
            yield from select(pending_selection.popleft())
        while pending:
            if chunk_observer is None:
                yield resolve(pending.popleft())
            else:
                yield from emit_observed(pending.popleft())
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
        pending_selection.clear()
        for stream in (copy_stream, result_stream, compute_stream) + ((select_stream,) if select_stream is not None else ()):
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
        graphs.clear()
        graph_inputs[:] = graph_outputs[:] = [None] * GRAPH_LANES
        release_futures[:] = [None] * len(release_futures)
        del device_buffers[:]
        for views in result_views:
            views.clear()
        del result_blocks[:]
        if cleanup_errors and not primary_error:
            raise cleanup_errors[0]


def design_column_block(n_samples, traits):
    """Source design upload/square workspace cap: 256 MiB per trait block."""
    if n_samples < 1 or traits < 1:
        raise ValueError('Positive design dimensions required')
    return max(1, min(traits, (1 << 28) // max(n_samples * 4, 1)))
