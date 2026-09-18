"""Tiled binary summary-statistic output.

Text summary statistics do not scale to whole-cohort multi-trait scans. One
marker-by-trait row costs roughly eighty ASCII bytes plus a Python csv round
trip, so 8.1M variants by 128 traits is about 88 GB of text and over a billion
formatted rows. The default association result is three float32 values per
cell: effect size, t statistic and an underflow-safe significance measure.

This module writes that native representation directly:

``beta.f32``
    ``(n_variants, n_traits)`` little-endian float32, C order.
``tstat.f32``
    ``(n_variants, n_traits)`` little-endian float32, C order.
``neglog10p.f32``
    ``(n_variants, n_traits)`` little-endian float32 ``-log10(P)``, C order.
``manifest.json``
    shape, dtype, byte order, trait names, sample count, residual degrees of
    freedom, and the exclusion convention.

Excluded variants (missing or invariant) carry NaN in every array, matching
the scan contract; per-category counts stay in ``qc.json``. Raw P values are
not stored because they are redundant with ``-log10(P)`` and underflow for
strong associations.

Writes are decoupled from scan chunking. Chunks are appended into a staging
block and handed to the operating system only once a full block is ready, so
the write block size is a property of the storage rather than of the scan
geometry. A bounded ring of blocks is drained by one background thread per
array, which overlaps the write with the next chunk's GPU work and applies
backpressure when storage cannot keep up.
"""

from __future__ import annotations

import json
import os
import queue
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

MANIFEST_NAME = "manifest.json"
BETA_NAME = "beta.f32"
TSTAT_NAME = "tstat.f32"
LOGP_NAME = "neglog10p.f32"
FORMAT_VERSION = 2
DEFAULT_BLOCK_BYTES = 16 << 20
MIN_AUTO_BLOCK_BYTES = 1 << 20
MAX_AUTO_BLOCK_BYTES = 16 << 20
DEFAULT_INFLIGHT_BYTES = 512 << 20
DEFAULT_QUEUE_DEPTH = 3
DEFAULT_WRITEBACK_BYTES = 64 << 20

_SYNC_FILE_RANGE_WAIT_BEFORE = 1
_SYNC_FILE_RANGE_WRITE = 2
_SYNC_FILE_RANGE_WAIT_AFTER = 4


def _load_sync_file_range():
    """Return ``sync_file_range`` where the platform provides it, else None.

    Without it the writer still works; writeback scheduling is then left to
    the kernel, which defers the cost to the closing fsync.
    """
    try:
        import ctypes

        libc = ctypes.CDLL("libc.so.6", use_errno=True)
        function = libc.sync_file_range
        function.argtypes = [
            ctypes.c_int,
            ctypes.c_int64,
            ctypes.c_int64,
            ctypes.c_uint,
        ]
        function.restype = ctypes.c_int
        return function
    except (OSError, AttributeError):
        return None


_SYNC_FILE_RANGE = _load_sync_file_range()


@dataclass
class _StreamStats:
    payload_bytes: int = 0
    write_seconds: float = 0.0
    writeback_seconds: float = 0.0
    drain_seconds: float = 0.0
    fsync_seconds: float = 0.0
    blocks: int = 0
    stall_seconds: float = 0.0


class _BlockStream:
    """Append-only file fed fixed-size blocks by a background thread."""

    def __init__(
        self,
        path: Path,
        block_bytes: int,
        queue_depth: int,
        writeback_bytes: int = DEFAULT_WRITEBACK_BYTES,
        inflight_bytes: int = DEFAULT_INFLIGHT_BYTES,
    ) -> None:
        if block_bytes <= 0:
            raise ValueError("block_bytes must be positive")
        if queue_depth <= 0:
            raise ValueError("queue_depth must be positive")
        self.path = path
        self._block_bytes = block_bytes
        self._writeback_bytes = writeback_bytes if _SYNC_FILE_RANGE else 0
        self._offset = 0
        self._writeback_started = 0
        self._writeback_waited = 0
        self._fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        self._ready: queue.Queue = queue.Queue()
        self._free: queue.Queue = queue.Queue()
        for _ in range(queue_depth + 1):
            self._free.put(bytearray(block_bytes))
        self._staging = self._free.get()
        self._filled = 0
        self._credit = threading.Condition()
        # Memory budget for unwritten output. Deliberately not derived from
        # block_bytes: a small block is what makes pass-through engage, and
        # scaling the budget with it just trades the staging copy for a
        # backpressure stall of the same size.
        self._credit_bytes = max(inflight_bytes, 2 * block_bytes)
        self._inflight = 0
        self._error: BaseException | None = None
        self._closed = False
        self.stats = _StreamStats()
        self._thread = threading.Thread(
            target=self._drain, name=f"torchgwas-write-{path.name}", daemon=True
        )
        self._thread.start()

    def _drain(self) -> None:
        while True:
            item = self._ready.get()
            if item is None:
                return
            block, length, pooled = item
            try:
                if self._error is None:
                    started = time.perf_counter()
                    view = memoryview(block)[:length]
                    while view:
                        view = view[os.write(self._fd, view):]
                    self.stats.write_seconds += time.perf_counter() - started
                    self.stats.payload_bytes += length
                    self.stats.blocks += 1
                    self._offset += length
                    self._advance_writeback()
            except BaseException as exc:  # surfaced on the next append or on close
                if self._error is None:
                    self._error = exc
            finally:
                if pooled:
                    self._free.put(block)
                else:
                    self._release_credit(length)

    def _advance_writeback(self) -> None:
        """Push completed ranges to storage instead of banking dirty pages.

        Deferring everything to the closing fsync would expose the whole
        output as a tail after the scan, and the dirty pages would evict the
        read-ahead the scan still needs. Writeback of a finished range starts
        asynchronously; the range before it is waited on and dropped from the
        page cache, which bounds dirty data to roughly two intervals and keeps
        the final fsync short.
        """
        if not self._writeback_bytes:
            return
        while self._offset - self._writeback_started >= self._writeback_bytes:
            start = self._writeback_started
            length = self._writeback_bytes
            _SYNC_FILE_RANGE(self._fd, start, length, _SYNC_FILE_RANGE_WRITE)
            self._writeback_started += length
            lag = self._writeback_started - self._writeback_waited
            if lag >= 2 * self._writeback_bytes:
                waited = time.perf_counter()
                wait_start = self._writeback_waited
                _SYNC_FILE_RANGE(
                    self._fd,
                    wait_start,
                    self._writeback_bytes,
                    _SYNC_FILE_RANGE_WAIT_BEFORE
                    | _SYNC_FILE_RANGE_WRITE
                    | _SYNC_FILE_RANGE_WAIT_AFTER,
                )
                self.stats.writeback_seconds += time.perf_counter() - waited
                os.posix_fadvise(
                    self._fd, wait_start, self._writeback_bytes, os.POSIX_FADV_DONTNEED
                )
                self._writeback_waited += self._writeback_bytes

    def _raise_pending(self) -> None:
        if self._error is not None:
            raise RuntimeError(f"sumstats write to {self.path} failed") from self._error

    def _hand_off(self) -> None:
        self._ready.put((self._staging, self._filled, True))
        started = time.perf_counter()
        self._staging = self._free.get()
        self.stats.stall_seconds += time.perf_counter() - started
        self._filled = 0

    def _release_credit(self, length: int) -> None:
        with self._credit:
            self._inflight -= length
            self._credit.notify_all()

    def _acquire_credit(self, length: int) -> float:
        """Bound bytes queued for pass-through writes, as the pool bounds copies."""
        started = time.perf_counter()
        with self._credit:
            while self._inflight and self._inflight + length > self._credit_bytes:
                self._credit.wait()
            self._inflight += length
        return time.perf_counter() - started

    def append(self, payload: memoryview, owner=None) -> None:
        """Queue ``payload`` for the writer thread.

        A payload the caller has given up ownership of and that is already at
        least ``block_bytes`` is queued as-is. Copying it into a staging block
        first would cost a full memcpy of the result stream for no benefit:
        coalescing only earns anything when chunks are smaller than the write
        block. Smaller or borrowed payloads take the staging path.
        """
        self._raise_pending()
        if self._closed:
            raise RuntimeError(f"stream {self.path} is closed")
        total = len(payload)
        if owner is not None and self._filled == 0 and total >= self._block_bytes:
            self.stats.stall_seconds += self._acquire_credit(total)
            self._ready.put((payload, total, False))
            return
        offset = 0
        while offset < total:
            take = min(self._block_bytes - self._filled, total - offset)
            self._staging[self._filled:self._filled + take] = payload[offset:offset + take]
            self._filled += take
            offset += take
            if self._filled == self._block_bytes:
                self._hand_off()

    def close(self, fsync: bool = True) -> _StreamStats:
        if self._closed:
            return self.stats
        self._closed = True
        try:
            if self._filled:
                self._hand_off()
            # Time the tail drain. With a deep ring most of the payload can
            # still be buffered when the scan finishes, so charging only the
            # write calls and the fsync would understate the write.
            drain_started = time.perf_counter()
            self._ready.put(None)
            self._thread.join()
            self.stats.drain_seconds = time.perf_counter() - drain_started
            self._raise_pending()
            if fsync:
                started = time.perf_counter()
                os.fsync(self._fd)
                self.stats.fsync_seconds = time.perf_counter() - started
        finally:
            os.close(self._fd)
        return self.stats

    def abort(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            self._ready.put(None)
            self._thread.join(timeout=30.0)
        finally:
            os.close(self._fd)


@dataclass
class BinarySumstatsWriter:
    """Write beta, t_stat and -log10(P) tiles into a binary directory.

    Chunks must arrive in variant order and must together cover exactly
    ``n_variants`` rows; ``close`` verifies that.
    """

    directory: Path
    n_variants: int
    trait_names: list
    n_samples: int
    df: int | list[int]
    block_bytes: int | None = None
    """Write block size, or None to derive it from the first chunk.

    Pass-through only engages for a chunk that is at least one block, so a
    block larger than the per-array chunk payload silently forces every chunk
    through a staging copy. Deriving the default from the observed payload
    makes borrowing engage by construction, rather than depending on the
    caller matching two otherwise unrelated numbers.
    """

    queue_depth: int = DEFAULT_QUEUE_DEPTH
    writeback_bytes: int = DEFAULT_WRITEBACK_BYTES
    inflight_bytes: int = DEFAULT_INFLIGHT_BYTES
    fsync: bool = True
    borrow_chunks: bool = True
    """Queue caller arrays directly instead of copying them into staging.

    Safe only where the chunk iterator yields owned arrays, which is the
    documented contract of the scan backends. Set False for an iterator that
    reuses its result buffers.
    """

    store_beta: bool = True
    """Write the effect-size array alongside the statistic.

    Dropping it reduces the output from 12 to 8 bytes per cell. It is a
    screening-only choice: without beta there is no effect size and no standard
    error, so results cannot be meta-analysed or expressed on any effect scale.
    Keep the default unless the output is only ever used to rank or threshold.
    """

    extra_manifest: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.directory = Path(self.directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.n_traits = len(self.trait_names)
        self._written = 0
        self._append_seconds = 0.0
        self._closed = False
        self._summary: dict = {}
        self._beta = None
        self._tstat = None
        self._logp = None
        # Distinct from "_beta is None", which is also the steady state of a
        # t-only store; block sizing keys off whether streams exist yet.
        self._streams_open = False
        if self.block_bytes is not None:
            self._open_streams(self.block_bytes)

    def _open_streams(self, block_bytes: int) -> None:
        self.block_bytes = block_bytes
        self._streams_open = True
        if self.store_beta:
            self._beta = _BlockStream(
                self.directory / BETA_NAME,
                block_bytes,
                self.queue_depth,
                self.writeback_bytes,
                self.inflight_bytes,
            )
        try:
            self._tstat = _BlockStream(
                self.directory / TSTAT_NAME,
                block_bytes,
                self.queue_depth,
                self.writeback_bytes,
                self.inflight_bytes,
            )
            self._logp = _BlockStream(
                self.directory / LOGP_NAME,
                block_bytes,
                self.queue_depth,
                self.writeback_bytes,
                self.inflight_bytes,
            )
        except BaseException:
            if self._beta is not None:
                self._beta.abort()
                self._beta = None
            if self._tstat is not None:
                self._tstat.abort()
                self._tstat = None
            self._streams_open = False
            raise

    @property
    def rows_written(self) -> int:
        """Number of marker-by-trait cells appended so far."""
        return self._written * self.n_traits

    def write_chunk(self, start: int, end: int, beta, t_stat, neg_log10_p=None) -> None:
        if self._closed:
            raise RuntimeError("writer is closed")
        if start != self._written:
            raise ValueError(
                "sumstats chunks must arrive in order: expected start "
                f"{self._written}, got {start}"
            )
        count = end - start
        started = time.perf_counter()
        if neg_log10_p is None:
            from .tails import upper_tail_log10_from_t

            df = np.asarray(self.df, dtype=np.float64)
            if df.ndim == 1:
                df = df[None, :]
            neg_log10_p = upper_tail_log10_from_t(t_stat, df)
        if not self._streams_open:
            payload = count * self.n_traits * 4
            self._open_streams(
                min(MAX_AUTO_BLOCK_BYTES, max(MIN_AUTO_BLOCK_BYTES, payload))
            )
        pairs = [(t_stat, self._tstat), (neg_log10_p, self._logp)]
        if self._beta is not None:
            pairs.insert(0, (beta, self._beta))
        for array, stream in pairs:
            block = np.ascontiguousarray(array, dtype="<f4")
            if block.shape != (count, self.n_traits):
                raise ValueError(
                    f"expected chunk shape {(count, self.n_traits)}, got {block.shape}"
                )
            # The scan contract hands over owned result arrays, so a block
            # that ascontiguousarray already had to copy, or that the caller
            # will not touch again, can go to the writer without a second
            # copy through a staging buffer.
            owner = block if self.borrow_chunks else None
            stream.append(memoryview(block).cast("B"), owner=owner)
        self._append_seconds += time.perf_counter() - started
        self._written = end

    def close(self) -> dict:
        if self._closed:
            return self._summary
        self._closed = True
        if not self._streams_open:
            # A writer that never received a chunk still has to produce the
            # files its manifest promises, and only those.
            self._open_streams(self.block_bytes or DEFAULT_BLOCK_BYTES)
        beta_stats = _StreamStats() if self._beta is None else self._beta.close(fsync=self.fsync)
        tstat_stats = self._tstat.close(fsync=self.fsync)
        logp_stats = self._logp.close(fsync=self.fsync)
        if self._written != self.n_variants:
            raise ValueError(
                f"sumstats covered {self._written} variants, expected {self.n_variants}"
            )
        arrays = {"t_stat": TSTAT_NAME, "neg_log10_p": LOGP_NAME}
        if self.store_beta:
            arrays = {"beta": BETA_NAME, **arrays}
        manifest = {
            "format": "torchgwas-binary-sumstats",
            "version": FORMAT_VERSION,
            "byte_order": "little",
            "dtype": "float32",
            "shape": [int(self.n_variants), int(self.n_traits)],
            "order": "C",
            "arrays": arrays,
            "n_samples": int(self.n_samples),
            "df": (
                int(self.df)
                if np.ndim(self.df) == 0
                else [int(value) for value in self.df]
            ),
            "traits": list(self.trait_names),
            "excluded_convention": (
                "NaN marks a missing or invariant variant in every stored array"
            ),
            "significance": (
                "neg_log10_p is -log10 of the exact two-sided Student-t tail"
            ),
            **(
                {}
                if self.store_beta
                else {
                    "beta": (
                        "not stored; this store is screening-only and carries no "
                        "effect size or standard error"
                    )
                }
            ),
            **self.extra_manifest,
        }
        (self.directory / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2))
        payload = (beta_stats.payload_bytes + tstat_stats.payload_bytes
                   + logp_stats.payload_bytes)
        write_seconds = (beta_stats.write_seconds + tstat_stats.write_seconds
                         + logp_stats.write_seconds)
        # Writer-thread time actually spent facing storage: the write calls,
        # the blocking writeback waits and the closing fsync.
        service_seconds = (
            write_seconds
            + beta_stats.writeback_seconds
            + tstat_stats.writeback_seconds
            + logp_stats.writeback_seconds
            + beta_stats.fsync_seconds
            + tstat_stats.fsync_seconds
            + logp_stats.fsync_seconds
        )
        self._summary = {
            "directory": str(self.directory),
            "cells": self._written * self.n_traits,
            "payload_bytes": payload,
            "append_seconds": self._append_seconds,
            "stall_seconds": (beta_stats.stall_seconds + tstat_stats.stall_seconds
                              + logp_stats.stall_seconds),
            "write_seconds": write_seconds,
            "writeback_seconds": (
                beta_stats.writeback_seconds + tstat_stats.writeback_seconds
                + logp_stats.writeback_seconds
            ),
            "drain_seconds": max(beta_stats.drain_seconds, tstat_stats.drain_seconds,
                                 logp_stats.drain_seconds),
            "fsync_seconds": (beta_stats.fsync_seconds + tstat_stats.fsync_seconds
                              + logp_stats.fsync_seconds),
            "blocks": beta_stats.blocks + tstat_stats.blocks + logp_stats.blocks,
            "block_bytes": self.block_bytes,
            "queue_depth": self.queue_depth,
            "writeback_bytes": self.writeback_bytes,
            "inflight_bytes": self.inflight_bytes,
            "sync_file_range": _SYNC_FILE_RANGE is not None,
            "service_seconds": service_seconds,
            "write_gb_s": (
                (payload / service_seconds / 1e9) if service_seconds > 0 else None
            ),
        }
        return self._summary

    def abort(self) -> None:
        if self._closed:
            return
        self._closed = True
        if self._beta is not None:
            self._beta.abort()
        if self._tstat is not None:
            self._tstat.abort()
        if self._logp is not None:
            self._logp.abort()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.close()
        else:
            self.abort()


def read_manifest(directory) -> dict:
    return json.loads((Path(directory) / MANIFEST_NAME).read_text())


def open_binary_sumstats(directory):
    """Memory-map beta, t_stat and -log10(P) from a binary directory."""
    directory = Path(directory)
    manifest = read_manifest(directory)
    if manifest.get("format") != "torchgwas-binary-sumstats":
        raise ValueError(f"{directory} is not a torchgwas binary sumstats store")
    shape = tuple(manifest["shape"])
    arrays = manifest["arrays"]

    def _map(name):
        # mmap rejects a zero-length file, which is what a store covering no
        # variants legitimately contains.
        if 0 in shape:
            return np.empty(shape, dtype="<f4")
        return np.memmap(directory / name, dtype="<f4", mode="r", shape=shape)

    # A screening-only store has no beta; callers must handle None rather than
    # silently reading a zero array.
    beta = _map(arrays["beta"]) if "beta" in arrays else None
    tstat = _map(arrays["t_stat"])
    logp = _map(arrays["neg_log10_p"])
    return beta, tstat, logp, manifest
