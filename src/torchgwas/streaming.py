from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from collections.abc import Callable, Iterator
from typing import Protocol, runtime_checkable

import numpy as np


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


class PinnedDosageLoader:
    """Bounded CPU decode/staging with explicit DMA buffer ownership."""

    def __init__(self, source, chunk_size, depth=3, reader_workers=None,
                 variant_range=None):
        import queue
        import warnings
        import threading
        import torch

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
        self.buffers = [torch.empty((chunk_size, row_width),
                                   dtype=torch_dtype, pin_memory=True)
                        for _ in range(depth)]
        for index in range(depth):
            self.free.put(index)

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
                    np.copyto(self.buffers[index][:count].numpy(), array.T)
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
                index, start, end, future = pending.popleft()
                future.result()
                publish((index, self.buffers[index][:end-start], start, end))
            try:
                # Executor exits first: no reader closes while a native fill owns a slot.
                with source.native_reader_session() as read_into:
                    with ThreadPoolExecutor(max_workers=workers,
                                            thread_name_prefix='torchgwas-pinned-decode') as pool:
                        for start in range(self.variant_start,
                                           self.variant_end, chunk_size):
                            while not self.stop.is_set():
                                try:
                                    index = self.free.get(timeout=0.1)
                                    break
                                except queue.Empty:
                                    pass
                            else:
                                break
                            end = min(self.variant_end, start + chunk_size)
                            destination = self.buffers[index][:end-start].numpy()
                            pending.append((index, start, end,
                                            pool.submit(read_into, start, end, destination)))
                            if len(pending) >= workers:
                                emit()
                        while pending and not self.stop.is_set():
                            emit()
                publish(None)
            except BaseException as error:
                publish(error)

        target = produce_direct if getattr(source, "allows_direct_native_fill", False) else produce
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
