from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from collections.abc import Callable, Iterator
from typing import Protocol, runtime_checkable

import numpy as np


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
        self.dtype = np.dtype(dtype)
        self.prefetch_chunks = int(prefetch_chunks)
        self.reader_workers = int(reader_workers)

    def __iter__(self) -> Iterator[tuple[int, int, np.ndarray]]:
        bounds = [
            (start, min(self.n_markers, start + self.chunk_size))
            for start in range(0, self.n_markers, self.chunk_size)
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
