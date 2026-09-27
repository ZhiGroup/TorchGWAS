"""Bounded changes to future chunks, with explicit delivery timing boundaries.

This is an execution mechanism, not a memory admission or tuning policy. The
caller must admit the fixed ring capacity plus the workspace and transient
allocation envelope for every permitted chunk shape and tail first.
"""
from dataclasses import dataclass
from numbers import Integral
import math
import threading
import time


def _positive_size(value, name):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return int(value)


class ChunkSizeControl:
    """Thread-safe choice among explicitly supplied chunk sizes.

    Changing size affects only reads not yet issued. Prefetched ranges retain
    their existing size and results; the ring allocation never changes.
    """
    def __init__(self, sizes, *, initial):
        values = tuple(_positive_size(v, 'chunk size') for v in sizes)
        if not values or len(set(values)) != len(values):
            raise ValueError('Nonempty unique chunk sizes required')
        self.sizes = tuple(sorted(values))
        self.capacity = max(self.sizes)
        self._lock = threading.Lock()
        self.set_size(initial)

    def set_size(self, size):
        size = _positive_size(size, 'chunk size')
        if size not in self.sizes:
            raise ValueError('Chunk size is outside the supplied choices')
        with self._lock:
            self._size = size

    def __call__(self, start, stop, capacity):
        if capacity < self.capacity:
            raise ValueError('Chunk ring is smaller than the supplied choices')
        with self._lock:
            return self._size


def aligned_chunk_sizes(sizes):
    """An explicit finite grid whose smallest member divides every choice."""
    values=tuple(_positive_size(value,'chunk size') for value in sizes)
    if not values or len(set(values))!=len(values):
        raise ValueError('Nonempty unique chunk sizes required')
    values=tuple(sorted(values))
    if any(value%values[0] for value in values):
        raise ValueError('Adaptive candidate sizes must be multiples of the smallest choice')
    return values


def aligned_chunk_shapes(sizes, markers):
    """All full sizes and the sole possible tail for a fixed source interval."""
    sizes=aligned_chunk_sizes(sizes)
    markers=_positive_size(markers,'markers')
    shapes={size for size in sizes if size<=markers}
    if markers%sizes[0]:shapes.add(markers%sizes[0])
    return tuple(sorted(shapes))


class AlignedChunkSizeControl(ChunkSizeControl):
    """Change future reads without introducing new intermediate tail shapes.

    Once fewer than the requested rows remain, issue the largest admitted full
    size that fits. Only the final remainder below the smallest choice is
    clipped. Every start remains on the smallest-size grid relative to its
    fixed source interval; all emitted shapes are known before execution.
    """
    def __init__(self,sizes,*,initial):
        super().__init__(aligned_chunk_sizes(sizes),initial=initial)

    def __call__(self,start,stop,capacity):
        requested=super().__call__(start,stop,capacity)
        remaining=stop-start
        if remaining<=0:raise ValueError('Nonempty source interval required')
        if requested<=remaining:return requested
        return next((size for size in reversed(self.sizes) if size<=remaining),remaining)


def selected_chunk_size(selector, start, stop, capacity):
    size = capacity if selector is None else selector(start, stop, capacity)
    size = _positive_size(size, 'selected chunk size')
    if size > capacity:
        raise ValueError('Selected chunk size exceeds the fixed ring capacity')
    return min(size, stop - start)


@dataclass(frozen=True)
class MinimalChunkObservation:
    """Chunk delivered and the consumer resumed; no read or CUDA timing.

    Emitted by paths without per-read instrumentation (packed BED, device
    BGEN decode, generic range reads). Only observers that declare
    `accepts_minimal_observations` (the empirical tuner) are given these.
    """
    start: int
    end: int
    capacity: int
    device: str
    completed: float
    boundary: str = 'chunk delivered and consumer resumed; no read or CUDA timing'


class MinimalDelivery:
    """Delivery wrapper for paths without read timing; reports on completion."""
    record_cuda = False

    def __init__(self, start, end, capacity, device, observer):
        self.start, self.end, self.capacity, self.device, self.observer = start, end, capacity, str(device), observer
        self.cuda = None

    def deliver(self, result):
        yield result

    def complete(self):
        self.observer(MinimalChunkObservation(self.start, self.end, self.capacity, self.device, time.perf_counter()))


def chunk_control_path(source, device):
    """Which chunk-switching mechanism a source offers on this device.

    'native' (direct native fill: PGEN, zstd store), 'packed' (BED packed
    reads), 'device' (BGEN GPU decoder slicing), 'range' (read_chunk: BGEN
    CPU decode, disk-backed matrices) or None.
    """
    backend = (source.resolve_decode_backend(device) if hasattr(source, 'resolve_decode_backend')
               else getattr(source, 'decode_backend', 'gpu'))
    if hasattr(source, 'iter_device_chunks') and backend != 'cpu':
        return 'device'
    if hasattr(source, 'iter_packed_chunks'):
        # Only with the native BED kernel, which accepts any chunk shape; the
        # compiled Torch fallback recompiles per shape (82 s at N=1M on H100).
        import os
        native = (getattr(source, '_sample_indices', None) is None
                  and os.environ.get('TORCHGWAS_BED_NATIVE', '1') == '1')
        if native:
            try:
                from . import scan_gpu
                native = bool(scan_gpu.available())
            except Exception:  # noqa: BLE001 - availability must never fail a scan
                native = False
        return 'packed' if native and hasattr(source, 'read_packed_into') else None
    if getattr(source, 'supports_fused_qc', False) and getattr(source, 'allows_direct_native_fill', False):
        return 'native'
    if hasattr(source, 'read_chunk'):
        return 'range'
    return None


def validate_chunk_control(source, device, compute_dtype, capacity, selector, observer):
    """Reject unsupported paths before phenotype preprocessing or allocation.

    Calibration/productive observers need the native dosage path's per-read
    instrumentation. The empirical tuner (accepts_minimal_observations) also
    works on packed BED, device BGEN and generic range readers.
    """
    if selector is None and observer is None:
        return
    for name, callback in [('chunk size selector', selector), ('chunk observer', observer)]:
        if callback is not None and not callable(callback):
            raise ValueError(f'{name} must be callable')
    capacity = _positive_size(capacity, 'explicit chunk capacity')
    if isinstance(selector, ChunkSizeControl) and selector.capacity > capacity:
        raise ValueError('Chunk ring is smaller than the supplied choices')
    if device.type != 'cuda' or compute_dtype != 'float32':
        raise ValueError('Chunk control requires direct native CPU fill and CUDA FP32 dosage statistics')
    empirical = getattr(observer, 'accepts_minimal_observations', False)
    if empirical and chunk_control_path(source, device) is not None:
        return
    if (not getattr(source, 'supports_fused_qc', False)
            or not getattr(source, 'allows_direct_native_fill', False)
            or hasattr(source, 'iter_packed_chunks')
            or getattr(source, 'native_encoding', 'dosage') != 'dosage'):
        raise ValueError('Chunk control requires direct native CPU fill and CUDA FP32 dosage statistics')
    backend = (source.resolve_decode_backend(device) if hasattr(source, 'resolve_decode_backend')
               else getattr(source, 'decode_backend', 'gpu'))
    if hasattr(source, 'iter_device_chunks') and backend != 'cpu':
        raise ValueError('Chunk control does not support device-native decoding')


@dataclass(frozen=True)
class ChunkObservation:
    """One completely delivered source chunk; times use perf_counter seconds.

    read_finished-read_started measures one reader's wall service, possibly
    concurrent with other reads. submitted is after fetch and before transfer
    submission. first_result includes dependency/queue waits, not just GPU work.
    consumer_seconds sums time suspended at yields. In an asynchronous driver
    that measures queue acceptance, not disk completion. completed excludes
    this observer's own cost and is emitted only after all yields resume.
    These overlapping intervals must not be summed as independent stage costs.
    CUDA event spans are optional and include stream scheduling gaps; they are
    not kernel-only service times. No profiler or extra stream drain is used.
    Reader CPU and runnable-wait counters are sampled only for reserved chunks.
    They bracket the read with small probe overhead and are diagnostic loaded
    observations, not independent service rates or a CPU-capacity estimate.
    Consumer CPU and runnable wait cover the downstream work before the next
    iterator request. Depending on output mode this can include selection,
    queue acceptance or synchronous writing; it is not durable-write service
    or a selector-only price.
    """
    start: int
    end: int
    capacity: int
    device: str
    read_started: float
    read_finished: float
    submitted: float
    first_result: float | None
    completed: float
    consumer_seconds: float
    result_blocks: int
    result_bytes: int
    cuda: 'ChunkDeviceTiming | None' = None
    boundary: str = 'source chunk through downstream iterator resumption; no durable-write guarantee'
    read_cpu_seconds: float | None = None
    read_runnable_wait_seconds: float | None = None
    read_probe_wall_seconds: float | None = None
    consumer_cpu_seconds: float | None = None
    consumer_runnable_wait_seconds: float | None = None
    consumer_probe_wall_seconds: float | None = None


class _ChunkDelivery:
    """Per-slot scalar bookkeeping, discarded immediately after observation."""
    def __init__(self, start, end, capacity, device, read_times, observer, *, record_cuda=False):
        self.start, self.end, self.capacity = start, end, capacity
        self.device = str(device)
        if len(read_times) not in (2, 5):
            raise ValueError('Reader timing must have two or five values')
        self.read_started, self.read_finished = read_times[:2]
        (self.read_cpu_seconds, self.read_runnable_wait_seconds,
         self.read_probe_wall_seconds) = (read_times[2:] if len(read_times) == 5
                                          else (None, None, None))
        self.submitted = time.perf_counter()
        self.first_result = None
        self.consumer_seconds = 0.
        self.consumer_cpu_seconds = 0.
        self.consumer_runnable_wait_seconds = 0.
        self.consumer_probe_wall_seconds = 0.
        self._consumer_cpu_valid = self._consumer_wait_valid = self._consumer_probe_valid = True
        self._resumed_blocks = 0
        self.result_blocks = self.result_bytes = 0
        self.observer = observer
        self.record_cuda = record_cuda
        self.cuda = None

    def deliver(self, result):
        from .streaming import _thread_runnable_wait_sample
        probe_started = time.perf_counter()
        thread_before = threading.get_ident()
        wait_before, enabled_before = _thread_runnable_wait_sample()
        cpu_before = time.thread_time()
        ready = time.perf_counter()
        if self.first_result is None:
            self.first_result = ready
        self.result_blocks += 1
        self.result_bytes += sum(value.nbytes for value in result[2:] if value is not None)
        yield result
        # Generator close/throw must not claim that the consumer accepted data.
        finished = time.perf_counter()
        cpu_after = time.thread_time()
        wait_after, enabled_after = _thread_runnable_wait_sample()
        probe_finished = time.perf_counter()
        self._resumed_blocks += 1
        self.consumer_seconds += finished - ready
        if threading.get_ident() != thread_before:
            self._consumer_cpu_valid = self._consumer_wait_valid = self._consumer_probe_valid = False
        else:
            self.consumer_probe_wall_seconds += max(0.,
                (ready - probe_started) + (probe_finished - finished))
            self.consumer_cpu_seconds += max(0., cpu_after - cpu_before)
            ambiguous = (wait_before == wait_after and
                         enabled_before is not True and enabled_after is not True)
            if (wait_before is None or wait_after is None or
                    enabled_before != enabled_after or ambiguous):
                self._consumer_wait_valid = False
            else:
                self.consumer_runnable_wait_seconds += max(0.,
                    (wait_after - wait_before) / 1e9)

    def complete(self):
        self.observer(ChunkObservation(
            self.start, self.end, self.capacity, self.device,
            self.read_started, self.read_finished, self.submitted,
            self.first_result, time.perf_counter(), self.consumer_seconds,
            self.result_blocks, self.result_bytes, cuda=self.cuda,
            read_cpu_seconds=self.read_cpu_seconds,
            read_runnable_wait_seconds=self.read_runnable_wait_seconds,
            read_probe_wall_seconds=self.read_probe_wall_seconds,
            consumer_cpu_seconds=(self.consumer_cpu_seconds if self._resumed_blocks == self.result_blocks
                                  and self.result_blocks
                                  and self._consumer_cpu_valid else None),
            consumer_runnable_wait_seconds=(self.consumer_runnable_wait_seconds
                if self._resumed_blocks == self.result_blocks and self.result_blocks
                and self._consumer_wait_valid else None),
            consumer_probe_wall_seconds=(self.consumer_probe_wall_seconds
                if self._resumed_blocks == self.result_blocks and self.result_blocks
                and self._consumer_probe_valid else None)))


@dataclass(frozen=True)
class ChunkDeviceTiming:
    """Seconds between CUDA events, not additive active-kernel times.

    Statistics includes joint reduction when requested. Device significance
    has no single result-stream transfer span; its value is None. Consumer
    time and host submission gaps can overlap these stream intervals.
    """
    h2d: float
    conversion: float
    statistics_and_reduction: float
    result_transfer: float | None


class InitialChunkMeasurements:
    """Sample a bounded early window of real work, separately on each GPU.

    Pass this as _chunk_observer. Measurements are spread over production
    chunks after warmup, retaining all GWAS output. Reserving at read issue
    bounds in-flight measurements as well as completed ones. Once the count
    or window limit is reached no additional timing events are recorded.
    The wall limit controls new reservations, not cancellation of issued work.
    This collects component observations; it neither fits GWAS runtime tables
    nor authorizes a configuration or claims an available hardware capacity.
    """
    def __init__(self, devices, *, max_chunks_per_device=8, warmup_chunks=1,
                 stride=4, max_window_seconds=10., cuda_events=True):
        self.devices = tuple(map(str, devices))
        if (not self.devices or len(set(self.devices)) != len(self.devices)
                or any(not d.startswith('cuda:') or not d[5:].isdigit() for d in self.devices)):
            raise ValueError('Unique explicit CUDA devices required')
        self.max_chunks = _positive_size(max_chunks_per_device, 'max_chunks_per_device')
        self.stride = _positive_size(stride, 'stride')
        if isinstance(warmup_chunks, bool) or not isinstance(warmup_chunks, Integral) or warmup_chunks < 0:
            raise ValueError('warmup_chunks must be a nonnegative integer')
        if isinstance(max_window_seconds, bool) or not math.isfinite(max_window_seconds) or max_window_seconds <= 0:
            raise ValueError('Positive finite max_window_seconds required')
        if type(cuda_events) is not bool:
            raise ValueError('cuda_events must be boolean')
        self.warmup_chunks = int(warmup_chunks)
        self.max_window_seconds = float(max_window_seconds)
        self.record_cuda = cuda_events
        self._lock = threading.Lock()
        self._started = None
        self._started_unix = None
        self._stopped = False
        self._seen = dict.fromkeys(self.devices, 0)
        self._reserved = dict.fromkeys(self.devices, 0)
        self._limits = dict.fromkeys(self.devices, self.max_chunks)
        self._pending = set()
        self._observations = []

    def reserve_read(self, start, end, device):
        with self._lock:
            if device not in self._seen:
                raise ValueError('Measurement device is outside the configured window')
            now = time.perf_counter()
            if self._started is None:
                self._started = now
                self._started_unix = time.time()
            index = self._seen[device]
            self._seen[device] += 1
            if (self._stopped or now - self._started >= self.max_window_seconds
                    or self._reserved[device] >= self._limits[device]
                    or index < self.warmup_chunks
                    or (index - self.warmup_chunks) % self.stride):
                return False
            key = (device, start, end)
            if key in self._pending:
                raise ValueError('Repeated source range in measurement window')
            self._pending.add(key)
            self._reserved[device] += 1
            return True

    def __call__(self, observation):
        with self._lock:
            key = (observation.device, observation.start, observation.end)
            if key not in self._pending:
                raise ValueError('Observation has no outstanding measurement reservation')
            self._pending.remove(key)
            self._observations.append(observation)

    def stop(self):
        """Stop new measurements; already-issued results are still accepted."""
        with self._lock:
            self._stopped = True

    def known_probe_wall_seconds(self):
        """Completed read/consumer probe brackets, excluding unobserved work.

        This is an additive conservative planning charge under overlap, not a
        hardware service price or an upper bound on future instrumentation.
        """
        with self._lock:
            return math.fsum(value for row in self._observations
                for value in (row.read_probe_wall_seconds,row.consumer_probe_wall_seconds)
                if value is not None)

    def snapshot(self):
        from dataclasses import asdict
        with self._lock:
            elapsed = 0. if self._started is None else time.perf_counter() - self._started
            exhausted = (self._stopped or elapsed >= self.max_window_seconds
                         or all(n >= self._limits[d] for d, n in self._reserved.items()))
            return dict(devices=list(self.devices), max_chunks_per_device=self.max_chunks,
                warmup_chunks=self.warmup_chunks, stride=self.stride,
                max_window_seconds=self.max_window_seconds, elapsed_seconds=elapsed,
                cuda_events=self.record_cuda, new_measurements_stopped=exhausted,
                source_chunks_seen=dict(self._seen), reserved=dict(self._reserved),
                current_limits=dict(self._limits), started_unix_seconds=self._started_unix,
                started_perf_counter=self._started,
                pending=[dict(device=d, start=s, end=e) for d, s, e in sorted(self._pending)],
                observations=[asdict(row) for row in self._observations],
                reader_probe_protocol='thread_cpu_schedstat_v1',
                reader_probe_scope='Reader thread only; scheduler wait is null if absent or ambiguous; probe bracketing wall time can include preemption and lies outside read wall service.',
                consumer_probe_protocol='thread_cpu_schedstat_probe_wall_v2',
                consumer_probe_scope='Downstream iterator suspension on its caller thread, including mode-dependent selection, queue acceptance or synchronous output. CPU and wait require normally resumed yields; scheduler wait is null if absent, ambiguous or resumed on a different thread. Probe bracketing wall time lies outside consumer_seconds; wait can exceed that narrower wall span. No durable-write or hardware-capacity inference.',
                scope='Bounded early production-chunk observations; overlapping host/CUDA spans, not independent hardware rates or durable output timings.')

    def publish(self, cache, *, dependencies, provenance, max_age_seconds):
        """Save completed observations for later jobs, with explicit expiry.

        Bind relevant source/library/device identities, workload dimensions,
        output mode, thread/math/storage context and the event protocol in
        dependencies. An incomplete window remains labeled incomplete; a cache
        hit never turns these loaded stream spans into hardware capacities.
        """
        value = self.snapshot()
        if not value['observations']:
            raise ValueError('No completed chunk measurements to publish')
        return cache.store('stage_observations', 'initial_chunks', value,
                           dependencies=dependencies, provenance=provenance,
                           max_age_seconds=max_age_seconds,
                           observed_unix_seconds=value["started_unix_seconds"])
