"""Tiled binary summary-statistic output.

Text summary statistics do not scale to whole-cohort multi-trait scans. One
marker-by-trait row costs roughly eighty ASCII bytes plus a Python csv round
trip, so 8.1M variants by 128 traits is about 88 GB of text and over a billion
formatted rows. The association result itself is two float32 values per cell,
which is 8 bytes: the same information at a tenth of the bytes and none of the
formatting work.

This module writes that native representation directly:

``beta.f32``
    ``(n_variants, n_traits)`` little-endian float32, C order.
``tstat.f32``
    ``(n_variants, n_traits)`` little-endian float32, C order.
``neglog10p.f32``
    ``(n_variants, n_traits)`` little-endian float32 ``-log10(P)``, C order:
    the exact two-sided Student-t tail, computed in FP64 on the scan device.
``manifest.json``
    shape, dtype, byte order, trait names, sample count, residual degrees of
    freedom, and the exclusion convention.

Excluded variants (missing or invariant) carry NaN in both arrays, matching the
scan contract; per-category counts stay in ``qc.json``. P-values are not stored
because they are an exact function of ``t_stat`` and ``df``; recomputing them
costs far less than writing them.

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
import uuid
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

MANIFEST_NAME = "manifest.json"
BETA_NAME = "beta.f32"
TSTAT_NAME = "tstat.f32"
DF_NAME = "df.f32"
LOGP_NAME = "neglog10p.f32"
# Version 2: neglog10p.f32 is stored (the release's format). A df sidecar is
# declared by a df dict in the manifest, not by the version.
FORMAT_VERSION = 2
DEFAULT_BLOCK_BYTES = 16 << 20
MIN_AUTO_BLOCK_BYTES = 1 << 20
MAX_AUTO_BLOCK_BYTES = 16 << 20
DEFAULT_INFLIGHT_BYTES = 512 << 20
DEFAULT_QUEUE_DEPTH = 3
DEFAULT_WRITEBACK_BYTES = 64 << 20

@dataclass(frozen=True)
class DenseWriteProgress:
    """New rows written in every stored beta/t array, in store coordinates.

    Successful write(2) progress, not queue acceptance, final fsync or manifest
    publication. The df sidecar may still be coalescing; its written prefix is
    separate. A range can split or combine scan chunks.
    """
    start: int
    end: int
    trait_range: tuple[int, int]
    statistic_bytes: int
    variant_df_complete_to: int | None
    completed: float
    directory: str
    device: str | None = None
    writer_queue: dict | None = None

    @property
    def rows(self):
        return (self.end-self.start)*(self.trait_range[1]-self.trait_range[0])


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


def write_manifest(directory, manifest, *, fsync=True):
    """Publish completion only after the data files have closed successfully."""
    directory=Path(directory)
    temporary=directory/(MANIFEST_NAME+'.'+uuid.uuid4().hex+'.tmp')
    try:
        with temporary.open('x',encoding='utf-8') as handle:
            json.dump(manifest,handle,indent=2,allow_nan=False)
            handle.write('\n')
            handle.flush()
            if fsync:
                os.fsync(handle.fileno())
        os.replace(temporary,directory/MANIFEST_NAME)
        if fsync and os.name=='posix':
            fd=os.open(directory,os.O_RDONLY|os.O_DIRECTORY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
    finally:
        temporary.unlink(missing_ok=True)


@dataclass
class _StreamStats:
    payload_bytes: int = 0
    write_seconds: float = 0.0
    writeback_seconds: float = 0.0
    drain_seconds: float = 0.0
    fsync_seconds: float = 0.0
    blocks: int = 0
    stall_seconds: float = 0.0
    progress_callback_seconds: float = 0.0


class _BlockStream:
    """Append-only file fed fixed-size blocks by a background thread."""

    def __init__(
        self,
        path: Path,
        block_bytes: int,
        queue_depth: int,
        writeback_bytes: int = DEFAULT_WRITEBACK_BYTES,
        inflight_bytes: int = DEFAULT_INFLIGHT_BYTES,
        *, on_written=None,
    ) -> None:
        if block_bytes <= 0:
            raise ValueError("block_bytes must be positive")
        if queue_depth <= 0:
            raise ValueError("queue_depth must be positive")
        self.path = path
        self._block_bytes = block_bytes
        self._queue_depth = queue_depth
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
        self._on_written = on_written
        self._queue_state_lock = threading.Lock() if on_written is not None else None
        self._queue_state = dict(accepted_bytes=0, written_bytes=0,
                                 staging_bytes=0, queued_bytes=0,
                                 active_bytes=0, queued_blocks=0)
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
            if self._queue_state_lock is not None:
                with self._queue_state_lock:
                    self._queue_state['queued_bytes'] -= length
                    self._queue_state['queued_blocks'] -= 1
                    self._queue_state['active_bytes'] += length
            try:
                if self._error is None:
                    started = time.perf_counter()
                    view = memoryview(block)[:length]
                    while view:
                        written = os.write(self._fd, view)
                        if written <= 0:
                            raise OSError('sumstats write made no progress')
                        view = view[written:]
                    self.stats.write_seconds += time.perf_counter() - started
                    self.stats.payload_bytes += length
                    self.stats.blocks += 1
                    self._offset += length
                    if self._queue_state_lock is not None:
                        with self._queue_state_lock:
                            self._queue_state['active_bytes'] -= length
                            self._queue_state['written_bytes'] += length
                    self._advance_writeback()
                    if self._on_written is not None:
                        notified = time.perf_counter()
                        try:
                            self._on_written(self._offset)
                        finally:
                            self.stats.progress_callback_seconds += time.perf_counter()-notified
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
        if self._queue_state_lock is not None:
            with self._queue_state_lock:
                self._queue_state['staging_bytes'] -= self._filled
                self._queue_state['queued_bytes'] += self._filled
                self._queue_state['queued_blocks'] += 1
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
            if self._queue_state_lock is not None:
                with self._queue_state_lock:
                    self._queue_state['accepted_bytes'] += total
                    self._queue_state['queued_bytes'] += total
                    self._queue_state['queued_blocks'] += 1
            self._ready.put((payload, total, False))
            return
        offset = 0
        while offset < total:
            take = min(self._block_bytes - self._filled, total - offset)
            # Assign through the destination view: bytearray slice assignment
            # first materializes a temporary bytearray from the source view.
            # This keeps a single synchronous copy before the scan may reuse
            # its source buffer, without another payload-sized allocation.
            memoryview(self._staging)[self._filled:self._filled + take] = payload[offset:offset + take]
            if self._queue_state_lock is not None:
                with self._queue_state_lock:
                    self._queue_state['accepted_bytes'] += take
                    self._queue_state['staging_bytes'] += take
            self._filled += take
            offset += take
            if self._filled == self._block_bytes:
                self._hand_off()

    def _queue_snapshot_unlocked(self) -> dict:
        """Called while this stream's state lock is held."""
        state = dict(self._queue_state)
        state['pending_write_bytes_interval'] = [
            state['staging_bytes'] + state['queued_bytes'],
            state['accepted_bytes'] - state['written_bytes']]
        state.update(block_bytes=self._block_bytes,
                     queue_depth=self._queue_depth,
                     error=self._error is not None,
                     scope='Active bytes include an entire block until its final write; a partially completed os.write is not credited.')
        return state

    def queue_snapshot(self) -> dict | None:
        """Accepted payload not yet fully written, at one stream-local instant."""
        if self._queue_state_lock is None:
            return None
        with self._queue_state_lock:
            return self._queue_snapshot_unlocked()

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
            self._on_written = None
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
            self._on_written = None
            os.close(self._fd)


@dataclass
class BinarySumstatsWriter:
    """Write beta and t_stat tiles into a binary sumstats directory.

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

    Dropping it halves the output to 4 bytes per cell, which on storage where
    the write does not hide is a direct halving of the write cost. It is a
    screening-only choice: without beta there is no effect size and no standard
    error, so results cannot be meta-analysed or expressed on any effect scale.
    Keep the default unless the output is only ever used to rank or threshold.
    """

    extra_manifest: dict = field(default_factory=dict)
    store_variant_df: bool = False
    """Store the scan's one-column per-variant df sidecar (four bytes/variant)."""
    on_write_progress: object = None
    """Optional callback on writer threads; must not wait for the scan.

    Called after every stored beta/t array has written a new common row prefix.
    It neither flushes staging nor changes block or fsync policy.
    """

    def __post_init__(self) -> None:
        if self.on_write_progress is not None and (not callable(self.on_write_progress) or not self.trait_names):
            raise ValueError('Dense write progress requires a callback and nonempty traits')
        self.directory = Path(self.directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.n_traits = len(self.trait_names)
        self._written = 0
        self._append_seconds = 0.0
        self._closed = False
        self._summary: dict = {}
        self._beta = None
        self._tstat = None
        self._df = None
        self._logp = None
        self._progress_lock = threading.RLock() if self.on_write_progress is not None else None
        self._progress_bytes = {'t_stat': 0, 'neg_log10_p': 0, **({'beta': 0} if self.store_beta else {}),
                                **({'df': 0} if self.store_variant_df else {})}
        self._progress_end = 0
        self._progress_error = None
        self._progress_directory = str(self.directory.resolve()) if self.on_write_progress is not None else None
        # Existing payloads are about to be replaced. A failed replacement
        # must not leave an old completion manifest advertising those files.
        (self.directory / MANIFEST_NAME).unlink(missing_ok=True)
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
                **self._progress_hook('beta'),
            )
        try:
            self._tstat = _BlockStream(
                self.directory / TSTAT_NAME,
                block_bytes,
                self.queue_depth,
                self.writeback_bytes,
                self.inflight_bytes,
                **self._progress_hook('t_stat'),
            )
            self._logp = _BlockStream(
                self.directory / LOGP_NAME,
                block_bytes,
                self.queue_depth,
                self.writeback_bytes,
                self.inflight_bytes,
                **self._progress_hook('neg_log10_p'),
            )
            if self.store_variant_df:
                self._df = _BlockStream(self.directory / DF_NAME, min(block_bytes,MIN_AUTO_BLOCK_BYTES),
                    self.queue_depth, self.writeback_bytes, self.inflight_bytes, **self._progress_hook('df'))
        except BaseException:
            for stream in (self._beta,self._tstat,self._logp,self._df):
                if stream is not None:
                    stream.abort()
            self._beta = self._tstat = self._logp = self._df = None
            self._streams_open = False
            raise

    def _progress_hook(self, field):
        if self.on_write_progress is None:
            return {}
        return dict(on_written=lambda offset: self._stream_written(field, offset))

    def _stream_written(self, field, offset):
        # Serialize cross-array notifications and their callbacks. Retain only
        # byte prefixes, never a growing list of pending chunks.
        with self._progress_lock:
            if self._progress_error is not None:
                raise RuntimeError('Dense write progress callback failed') from self._progress_error
            self._progress_bytes[field] = offset
            matrix_fields = ('beta', 't_stat', 'neg_log10_p') if self.store_beta else ('t_stat', 'neg_log10_p')
            end = min(self._progress_bytes[name]//(4*self.n_traits) for name in matrix_fields)
            if end <= self._progress_end:
                return
            event = DenseWriteProgress(self._progress_end, end, (0, self.n_traits),
                (end-self._progress_end)*4*self.n_traits*len(matrix_fields),
                self._progress_bytes['df']//4 if self.store_variant_df else None,
                time.perf_counter(), self._progress_directory,
                writer_queue=self.queue_snapshot())
            try:
                self.on_write_progress(event)
            except BaseException as error:
                self._progress_error = error
                raise
            self._progress_end = end

    def queue_snapshot(self) -> dict:
        """One atomic accepted/written snapshot across this writer's streams."""
        began = time.perf_counter()
        active = [(name, stream) for name, stream in
                  (('beta', self._beta), ('t_stat', self._tstat), ('neg_log10_p', self._logp), ('df', self._df))
                  if stream is not None]
        locks = [stream._queue_state_lock for _, stream in active]
        if any(lock is None for lock in locks):
            streams = {name: stream.queue_snapshot() for name, stream in active}
            atomic = False
        else:
            for lock in locks:
                lock.acquire()
            try:
                streams = {name: stream._queue_snapshot_unlocked()
                           for name, stream in active}
            finally:
                for lock in reversed(locks):
                    lock.release()
            atomic = True
        finished = time.perf_counter()
        valid = all(row is not None and not row['error'] and
            row['accepted_bytes'] - row['written_bytes'] ==
            row['staging_bytes'] + row['queued_bytes'] + row['active_bytes']
            for row in streams.values())
        return dict(kind='torchgwas.dense_writer_queue_observation.v1',
                    capture_started_seconds=began,
                    capture_finished_seconds=finished,
                    streams=streams, valid=valid, atomic_writer_streams=atomic,
                    scope='Accepted/written counters are captured while all this writer\'s stream state locks are held. Active bytes conservatively include partial writes. Producer/GPU work, other writer directories and durable storage are not synchronized.')

    @property
    def rows_written(self) -> int:
        """Number of marker-by-trait cells appended so far."""
        return self._written * self.n_traits

    def write_chunk(self, start: int, end: int, beta, t_stat, neg_log10_p=None, *, variant_df=None) -> None:
        """Append rows [start, end). neg_log10_p is the scan's device result
        (compute_log10_p); without it the host computes it from t and df, which
        is exact but costs ~0.1 ms per thousand cells on one core."""
        if self._closed:
            raise RuntimeError("writer is closed")
        if start != self._written:
            raise ValueError(
                "sumstats chunks must arrive in order: expected start "
                f"{self._written}, got {start}"
            )
        count = end - start
        if count <= 0 or end > self.n_variants:
            raise ValueError('sumstats chunk must be nonempty and within the variant extent')
        if self.store_variant_df:
            df_block = np.ascontiguousarray(variant_df, dtype='<f4')
            if df_block.shape != (count,1) or not np.isfinite(df_block).all():
                raise ValueError('variant_df must be a finite (variant, 1) array')
        elif variant_df is not None:
            raise ValueError('variant_df requires store_variant_df=True')
        started = time.perf_counter()
        if neg_log10_p is None:
            from .tails import upper_tail_log10_from_t
            pair_df = (np.asarray(variant_df, dtype=np.float64) if variant_df is not None else
                       np.asarray(self.df, dtype=np.float64)[None, :] if np.ndim(self.df) else float(self.df))
            neg_log10_p = upper_tail_log10_from_t(t_stat, pair_df).astype(np.float32)
        if not self._streams_open:
            payload = count * self.n_traits * 4
            self._open_streams(
                min(MAX_AUTO_BLOCK_BYTES, max(MIN_AUTO_BLOCK_BYTES, payload))
            )
        pairs = [(t_stat, self._tstat), (neg_log10_p, self._logp)]
        if self._beta is not None:
            pairs.insert(0, (beta, self._beta))
        if self._df is not None:
            pairs.append((df_block,self._df))
        for array, stream in pairs:
            block = np.ascontiguousarray(array, dtype="<f4")
            expected = (count,1) if stream is self._df else (count,self.n_traits)
            if block.shape != expected:
                raise ValueError(
                    f"expected chunk shape {expected}, got {block.shape}"
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
        closed_stats, errors = [], []
        for stream in (self._beta,self._tstat,self._logp,self._df):
            try:
                closed_stats.append(_StreamStats() if stream is None else stream.close(fsync=self.fsync))
            except BaseException as error:
                closed_stats.append(_StreamStats())
                errors.append(error)
        if errors:
            raise errors[0]
        beta_stats, tstat_stats, logp_stats, df_stats = closed_stats
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
            "df": (dict(array=DF_NAME,axis='variant',shape=[int(self.n_variants),1],dtype='float32') if self.store_variant_df else
            (
                int(self.df)
                if np.ndim(self.df) == 0
                else [int(value) for value in self.df]
            )),
            "traits": list(self.trait_names),
            "excluded_convention": (
                "NaN marks a missing or invariant variant in every stored array"
            ),
            "significance": ("neg_log10_p is -log10 of the exact two-sided Student-t tail at "
                             + ("the per-variant df in df.f32" if self.store_variant_df else
                                "the per-trait df" if np.ndim(self.df) else "df")),
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
        write_manifest(self.directory,manifest,fsync=self.fsync)
        payload = sum(stats.payload_bytes for stats in closed_stats)
        write_seconds = sum(stats.write_seconds for stats in closed_stats)
        # Writer-thread time actually spent facing storage: the write calls,
        # the blocking writeback waits and the closing fsync.
        service_seconds = (
            write_seconds
            + beta_stats.writeback_seconds
            + tstat_stats.writeback_seconds
            + beta_stats.fsync_seconds
            + tstat_stats.fsync_seconds
            + df_stats.writeback_seconds + df_stats.fsync_seconds
            + logp_stats.writeback_seconds + logp_stats.fsync_seconds
        )
        self._summary = {
            "directory": str(self.directory),
            "cells": self._written * self.n_traits,
            "payload_bytes": payload,
            "append_seconds": self._append_seconds,
            "stall_seconds": sum(stats.stall_seconds for stats in closed_stats),
            "progress_callback_seconds": sum(stats.progress_callback_seconds for stats in closed_stats),
            "write_seconds": write_seconds,
            "writeback_seconds": (
                sum(stats.writeback_seconds for stats in closed_stats)
            ),
            "drain_seconds": max(stats.drain_seconds for stats in closed_stats),
            "fsync_seconds": sum(stats.fsync_seconds for stats in closed_stats),
            "blocks": sum(stats.blocks for stats in closed_stats),
            "df_payload_bytes": df_stats.payload_bytes,
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
        if self._df is not None:
            self._df.abort()
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


def open_binary_df(directory):
    """Return scalar, trait-row or mapped variant-column df for broadcasting."""
    directory=Path(directory)
    manifest=read_manifest(directory)
    if manifest.get('layout') in ('trait_tiles','variant_shards') and isinstance(manifest.get('df'),list):
        # Missing phenotypes: one df per trait, as in a single store.
        open_binary_sumstats(directory)
        return np.asarray(manifest['df'])[None,:]
    if manifest.get('layout')=='trait_tiles':
        from .sumstats_tiled import TiledSumstatsArray,open_trait_tiled_sumstats
        open_trait_tiled_sumstats(directory,manifest)
        return TiledSumstatsArray(directory,manifest,'df')
    if manifest.get('layout')=='variant_shards':
        from .sumstats_sharded import VariantShardedArray,open_variant_sharded_sumstats
        open_variant_sharded_sumstats(directory,manifest)
        return VariantShardedArray(directory,manifest,'df')
    df=manifest['df']
    if isinstance(df,dict):
        if df.get('axis')!='variant' or df.get('shape')!=[manifest['shape'][0],1] or df.get('dtype')!='float32':
            raise ValueError('Unsupported df sidecar layout')
        path=directory/df['array']
        if path.stat().st_size!=manifest['shape'][0]*4:
            raise ValueError('df sidecar length differs from manifest')
        if manifest['shape'][0]==0:
            return np.empty((0,1),dtype='<f4')
        return np.memmap(path,dtype='<f4',mode='r',shape=tuple(df['shape']))
    return np.asarray(df)[None,:] if np.ndim(df) else df


def open_binary_sumstats(directory):
    """Memory-map (beta, t_stat, neg_log10_p, manifest) from a binary sumstats directory.

    beta is None for a t-only store and neg_log10_p for a store written
    before it was stored (version 1).
    """
    directory = Path(directory)
    manifest = read_manifest(directory)
    if manifest.get("format") != "torchgwas-binary-sumstats":
        raise ValueError(f"{directory} is not a torchgwas binary sumstats store")
    if manifest.get('layout')=='trait_tiles':
        from .sumstats_tiled import open_trait_tiled_sumstats
        return open_trait_tiled_sumstats(directory,manifest)
    if manifest.get('layout')=='variant_shards':
        from .sumstats_sharded import open_variant_sharded_sumstats
        return open_variant_sharded_sumstats(directory,manifest)
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
    logp = _map(arrays["neg_log10_p"]) if "neg_log10_p" in arrays else None
    return beta, tstat, logp, manifest
