from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from collections.abc import Callable, Iterator
from typing import Protocol, runtime_checkable

import numpy as np


def _thread_runnable_wait_sample():
    """Best-effort (scheduler-wait ns, schedstats switch) for this thread."""
    try:
        with open('/proc/sys/kernel/sched_schedstats', encoding='ascii') as setting:
            switch = setting.read().strip() == '1'
    except OSError:
        switch = None
    try:
        with open('/proc/thread-self/schedstat', encoding='ascii') as stat:
            return int(stat.read().split()[1]), switch
    except (OSError, ValueError, IndexError):
        return None, switch


def _resolve_variant_range(variant_range, n_markers):
    """Half-open (start, end) in the source's own variant numbering.

    Keeping shard results in the file's numbering is what lets them be
    reassembled without a translation table.
    """
    if variant_range is None:
        return 0, int(n_markers)
    start, end = variant_range
    start = int(start)
    end = int(n_markers) if end is None else int(end)
    if not 0 <= start <= end <= int(n_markers):
        raise ValueError(
            f"variant range ({start}, {end}) is outside 0..{n_markers}")
    return start, end


from contextlib import contextmanager


@contextmanager
def _range_reader_session(source, dtype):
    """A native-fill style reader over source.read_chunk(start, end, dtype).

    read_chunk returns (samples, count); pinned rows are variant-major, so
    each fill copies the transpose into the destination slot.
    """
    def read_into(start, end, destination):
        np.copyto(destination, source.read_chunk(start, end, dtype).T)
    yield read_into


class PinnedDosageLoader:
    """Bounded CPU decode/staging with explicit DMA buffer ownership."""

    def __init__(self, source, chunk_size, depth=3, reader_workers=None,
                 variant_range=None, chunk_size_selector=None, record_timing=False):
        import queue
        import warnings
        import threading
        import torch

        if chunk_size_selector is not None and not callable(chunk_size_selector):
            raise ValueError('chunk size selector must be callable')
        if type(record_timing) is not bool and not callable(record_timing):
            raise ValueError('record_timing must be boolean or a reservation callback')
        direct = getattr(source, 'allows_direct_native_fill', False)
        controlled = chunk_size_selector is not None or bool(record_timing)
        # Variable-size reads need a range reader: a native fill (PGEN, zstd)
        # or read_chunk (BGEN CPU decode, disk-backed matrices) via an adapter.
        if controlled and not direct and not hasattr(source, 'read_chunk'):
            raise ValueError('Chunk control requires a direct native reader')
        self._queue_module = queue
        self.variant_start, self.variant_end = _resolve_variant_range(
            variant_range, int(source.shape[1]))
        self.stop = threading.Event()
        self.free = queue.Queue()
        self.ready = queue.Queue(maxsize=depth)
        native_dtype = np.dtype(getattr(source, "native_transfer_dtype",
                                       getattr(source, "native_dtype", np.float32)))
        row_width = int(getattr(source, "native_row_width", source.shape[0]))
        if row_width != source.shape[0] and not getattr(source, "allows_direct_native_fill", False):
            raise ValueError("packed transfer requires a direct native reader")
        torch_dtype = {np.dtype(np.int8): torch.int8, np.dtype(np.uint8): torch.uint8}.get(
            native_dtype, torch.float32)
        numpy_dtype = {torch.int8: np.int8, torch.uint8: np.uint8}.get(torch_dtype, np.float32)
        self.capacity, self.depth = int(chunk_size), int(depth)
        self.slot_shape, self.slot_dtype = (int(chunk_size), row_width), torch_dtype
        # Slots are pinned when first filled (_slot), not all at the capacity
        # here: a tuned scan starts at 1024 and may never use its 4096-row
        # capacity. H100 JAGWAS, 4 shards: pinning the capacity up front cost
        # 2.9 s of cudaHostAlloc against 0.8 s at the start size.
        self.buffers = [None] * depth
        for index in range(depth):
            self.free.put(index)

        def slot(index, rows):
            """Slot `index`, holding at least `rows` rows; taken only after its DMA completed."""
            buffer = self.buffers[index]
            if buffer is None or buffer.shape[0] < rows:
                # Sized for its first chunk when sizes vary; grows once, to the capacity.
                size = rows if buffer is None and chunk_size_selector is not None else self.capacity
                self.buffers[index] = buffer = torch.empty((size, row_width), dtype=torch_dtype, pin_memory=True)
            return buffer
        self._slot = slot

        def publish(value):
            while not self.stop.is_set():
                try:
                    self.ready.put(value, timeout=0.1)
                    return
                except queue.Full:
                    pass

        def produce():
            iterator = None
            try:
                ranged = (self.variant_start, self.variant_end) != (0, int(source.shape[1]))
                extra = ({"variant_range": (self.variant_start, self.variant_end)}
                         if ranged else {})
                try:
                    iterator = (source.iter_native_chunks(chunk_size, reader_workers=reader_workers,
                                                          prefetch_chunks=depth, **extra)
                                if hasattr(source, "iter_native_chunks") else
                                source.iter_chunks(chunk_size, dtype=np.float32,
                                                   reader_workers=reader_workers,
                                                   prefetch_chunks=depth, **extra))
                except TypeError as error:
                    if not ranged:
                        raise
                    raise TypeError(
                        f"{type(source).__name__} cannot read a variant range, "
                        f"so it cannot be sharded across devices"
                    ) from error
                for start, end, array in iterator:
                    while not self.stop.is_set():
                        try:
                            index = self.free.get(timeout=0.1)
                            break
                        except queue.Empty:
                            pass
                    else:
                        return
                    count = end - start
                    np.copyto(slot(index, count)[:count].numpy(), array.T)
                    publish((index, self.buffers[index][:count], start, end))
                publish(None)
            except BaseException as error:
                publish(error)
            finally:
                if hasattr(iterator, "close"):
                    iterator.close()

        def produce_direct():
            from collections import deque
            requested = int(reader_workers or getattr(source, 'decode_workers', 1))
            # A fill occupies a ring slot for its duration, so concurrency
            # cannot exceed the depth however many workers were asked for.
            workers = min(depth, requested)
            self.decode_workers_requested = requested
            self.decode_workers_effective = workers
            if workers < requested:
                warnings.warn(
                    f"decode concurrency limited to {workers} by prefetch depth "
                    f"{depth}, not the requested {requested} workers; raise "
                    f"prefetch_chunks to use more decoders, at the cost of one "
                    f"chunk buffer per added slot",
                    RuntimeWarning, stacklevel=2)
            pending = deque()
            def emit():
                index, start, end, future, measured = pending.popleft()
                result = future.result()
                read_times = result if measured else None
                item = (index, self.buffers[index][:end-start], start, end)
                publish((*item, read_times) if record_timing else item)
            try:
                # Executor exits first: no reader closes while a native fill owns a slot.
                with (source.native_reader_session() if direct else
                      _range_reader_session(source, np.dtype(numpy_dtype))) as read_into:
                    with ThreadPoolExecutor(max_workers=workers,
                                            thread_name_prefix='torchgwas-pinned-decode') as pool:
                        from .adaptive_chunks import selected_chunk_size
                        import time
                        def fill(start, end, destination):
                            probe_started = time.perf_counter()
                            wait_before, enabled_before = _thread_runnable_wait_sample()
                            cpu_before = time.thread_time()
                            started = time.perf_counter()
                            read_into(start, end, destination)
                            finished = time.perf_counter()
                            cpu_seconds = max(0., time.thread_time() - cpu_before)
                            wait_after, enabled_after = _thread_runnable_wait_sample()
                            probe_finished = time.perf_counter()
                            # Some hosts expose a useful per-thread counter
                            # even with the global switch off. An unchanged
                            # counter in that state cannot distinguish no wait
                            # from an inactive counter.
                            wait_ambiguous = (wait_before == wait_after
                                and enabled_before is not True and enabled_after is not True)
                            runnable_wait_seconds = (None if wait_before is None or wait_after is None
                                or enabled_before != enabled_after or wait_ambiguous
                                else max(0., (wait_after - wait_before) / 1e9))
                            probe_wall_seconds = max(0., (started - probe_started) + (probe_finished - finished))
                            return (started, finished, cpu_seconds, runnable_wait_seconds,
                                    probe_wall_seconds)
                        start = self.variant_start
                        while start < self.variant_end:
                            while not self.stop.is_set():
                                try:
                                    index = self.free.get(timeout=0.1)
                                    break
                                except queue.Empty:
                                    pass
                            else:
                                break
                            # Select only after acquiring a free slot: queued
                            # work keeps its range, and changes apply to the next
                            # undecoded rows without restarting the input pass.
                            count = selected_chunk_size(chunk_size_selector, start,
                                                        self.variant_end, chunk_size)
                            end = start + count
                            destination = slot(index, count)[:count].numpy()
                            measured = record_timing(start, end) if callable(record_timing) else record_timing
                            if type(measured) is not bool:
                                raise ValueError('Timing reservation must return a boolean')
                            pending.append((index, start, end,
                                            pool.submit(fill if measured else read_into, start, end, destination), measured))
                            start = end
                            if len(pending) >= workers:
                                emit()
                        while pending and not self.stop.is_set():
                            emit()
                publish(None)
            except BaseException as error:
                publish(error)

        target = produce_direct if direct or controlled else produce
        self.thread = threading.Thread(target=target, daemon=True,
                                       name="torchgwas-pinned-dosage")
        self.thread.start()

    def __iter__(self):
        while True:
            item = self.ready.get()
            if item is None:
                return
            if isinstance(item, BaseException):
                raise item
            yield item

    def release(self, index):
        self.free.put(index)

    def close(self):
        self.stop.set()
        self.thread.join()
        # Release the pinned staging buffers, not just the producer thread.
        #
        # These are `depth` page-locked host tensors of (chunk_size, row_width)
        # -- about 891 MB at depth 16, a 10,000-variant chunk and a PGEN row.
        # Pinned memory is expensive to allocate and is not handed back by the
        # caching allocator, so closing without dropping these left every
        # scan's staging resident and a caller scanning repeatedly in one
        # process accumulated them. The queues hold indices into this list, so
        # they are drained too rather than pinning it alive.
        self.buffers = []
        for pending in (self.free, self.ready):
            while True:
                try:
                    pending.get_nowait()
                except self._queue_module.Empty:
                    break


ChunkReader = Callable[[int, int, np.dtype], np.ndarray]


@runtime_checkable
class ChunkedGenotype(Protocol):
    sample_ids: np.ndarray
    marker_ids: np.ndarray

    @property
    def shape(self) -> tuple[int, int]: ...

    @property
    def genotype(self): ...

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
    ) -> Iterator[tuple[int, int, np.ndarray]]: ...


class OrderedChunkLoader:
    """Bounded, ordered multi-worker chunk prefetch.

    Workers may complete out of order, but chunks are yielded in marker order.
    Calling ``Future.result()`` on the consumer thread also guarantees that any
    read/decode exception is propagated instead of being mistaken for EOF.
    """

    def __init__(
        self,
        n_markers: int,
        read_chunk: ChunkReader,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int = 4,
        reader_workers: int = 4,
        variant_range: tuple[int, int] | None = None,
    ) -> None:
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if prefetch_chunks <= 0:
            raise ValueError("prefetch_chunks must be positive")
        if reader_workers <= 0:
            raise ValueError("reader_workers must be positive")
        self.n_markers = int(n_markers)
        self.read_chunk = read_chunk
        self.chunk_size = int(chunk_size)
        start, end = _resolve_variant_range(variant_range, self.n_markers)
        self.variant_start = start
        self.variant_end = end
        self.dtype = np.dtype(dtype)
        self.prefetch_chunks = int(prefetch_chunks)
        self.reader_workers = int(reader_workers)

    def __iter__(self) -> Iterator[tuple[int, int, np.ndarray]]:
        bounds = [
            (start, min(self.variant_end, start + self.chunk_size))
            for start in range(self.variant_start, self.variant_end,
                               self.chunk_size)
        ]
        if not bounds:
            return

        max_pending = min(self.prefetch_chunks, len(bounds))
        with ThreadPoolExecutor(max_workers=self.reader_workers, thread_name_prefix="torchgwas-reader") as pool:
            pending: dict[int, Future[np.ndarray]] = {}
            submit_index = 0

            def submit_one() -> None:
                nonlocal submit_index
                start, end = bounds[submit_index]
                pending[submit_index] = pool.submit(self.read_chunk, start, end, self.dtype)
                submit_index += 1

            while submit_index < max_pending:
                submit_one()

            for yield_index, (start, end) in enumerate(bounds):
                future = pending.pop(yield_index)
                chunk = future.result()
                expected_shape = (chunk.shape[0], end - start)
                if chunk.ndim != 2 or chunk.shape != expected_shape:
                    raise ValueError(
                        f"reader returned shape {chunk.shape} for markers [{start}, {end}); "
                        f"expected (*, {end - start})"
                    )
                yield start, end, chunk
                if submit_index < len(bounds):
                    submit_one()


class SelectedRangeLoader:
    """Ordered prefetch whose chunk bounds come from a selector at submission.

    The generic-path counterpart of OrderedChunkLoader for adaptive chunk
    sizes: reads go through read_chunk(start, end, dtype) in a small pool, in
    order, and a size change applies to the next unread rows only.
    """

    def __init__(self, read_chunk, capacity, selector, dtype, prefetch_chunks, reader_workers,
                 variant_range):
        self.read_chunk, self.capacity, self.selector = read_chunk, int(capacity), selector
        self.dtype = np.dtype(dtype)
        self.prefetch_chunks, self.reader_workers = max(1, int(prefetch_chunks)), max(1, int(reader_workers))
        self.variant_start, self.variant_end = variant_range

    def __iter__(self):
        from collections import deque
        from .adaptive_chunks import selected_chunk_size
        pending = deque()
        cursor = self.variant_start
        with ThreadPoolExecutor(max_workers=self.reader_workers, thread_name_prefix="torchgwas-range") as pool:
            def submit():
                nonlocal cursor
                if cursor >= self.variant_end:
                    return False
                count = selected_chunk_size(self.selector, cursor, self.variant_end, self.capacity)
                start, end = cursor, cursor + count
                pending.append((start, end, pool.submit(self.read_chunk, start, end, self.dtype)))
                cursor = end
                return True
            for _ in range(self.prefetch_chunks):
                if not submit():
                    break
            while pending:
                start, end, future = pending.popleft()
                chunk = future.result()
                if chunk.ndim != 2 or chunk.shape[1] != end - start:
                    raise ValueError(f"reader returned shape {chunk.shape} for markers [{start}, {end})")
                yield start, end, chunk
                submit()
