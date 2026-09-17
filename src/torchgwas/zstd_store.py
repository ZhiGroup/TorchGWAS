from __future__ import annotations

import contextlib
import csv
import json
import os
import sys
import time
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import zstandard as zstd

from . import metadata_cache

_VARIANT_CACHE_ARRAYS = ("marker_id", "chromosome", "position",
                         "effect_allele", "other_allele")




DEFAULT_ZSTD_LEVEL = 15
DEFAULT_ZSTD_CHUNK_SIZE = 2500


def _store_path(prefix: str | Path, suffix: str) -> Path:
    return Path(f"{Path(prefix)}{suffix}")


def encode_zstd_store(
    source,
    output_prefix: str | Path,
    *,
    chunk_size: int = DEFAULT_ZSTD_CHUNK_SIZE,
    level: int = DEFAULT_ZSTD_LEVEL,
    compression_workers: int = 4,
) -> dict[str, int | float | str]:
    """Encode sample-by-variant uint8 dosage codes into independent zstd frames.

    Each frame contains at most ``chunk_size`` complete variant rows in
    variant-major uint8 order. Frames are compressed concurrently and written
    in marker order, allowing bounded parallel decompression during scans.
    """

    if chunk_size <= 0 or compression_workers <= 0:
        raise ValueError("chunk_size and compression_workers must be positive")
    zst_path = _store_path(output_prefix, ".zst")
    idx_path = _store_path(output_prefix, ".idx.npz")
    samples_path = _store_path(output_prefix, ".samples.tsv")
    variants_path = _store_path(output_prefix, ".variants.tsv")
    manifest_path = _store_path(output_prefix, ".complete.json")
    zst_path.parent.mkdir(parents=True, exist_ok=True)

    offsets: list[int] = []
    sizes: list[int] = []
    rows: list[int] = []
    position = 0

    compression_worker_seconds = 0.0

    def compress(payload: bytes) -> tuple[bytes, float]:
        started = time.perf_counter()
        blob = zstd.ZstdCompressor(level=level).compress(payload)
        return blob, time.perf_counter() - started

    pending: deque[tuple[int, int, Future[tuple[bytes, float]]]] = deque()
    with zst_path.open("wb", buffering=4 << 20) as output, ThreadPoolExecutor(
        max_workers=compression_workers,
        thread_name_prefix="torchgwas-zstd-encode",
    ) as pool:

        def write_oldest() -> None:
            nonlocal compression_worker_seconds, position
            start, end, future = pending.popleft()
            blob, compression_seconds = future.result()
            compression_worker_seconds += compression_seconds
            output.write(blob)
            offsets.append(position)
            sizes.append(len(blob))
            rows.append(end - start)
            position += len(blob)

        # A cohort-scale conversion runs for the better part of an hour and
        # used to print nothing at all, so a run that died at minute fifty was
        # indistinguishable from one that died at minute two -- which is
        # exactly what happened, and cost an hour to notice. Progress goes to
        # stderr so it never contaminates anything parsing stdout.
        total_variants = int(source.shape[1])
        started_at = time.perf_counter()
        last_report = started_at
        report_every = float(os.environ.get("TORCHGWAS_ZSTD_PROGRESS_SECONDS", "30"))

        for start, end, chunk in source.iter_chunks(
            chunk_size=chunk_size,
            dtype=np.uint8,
            prefetch_chunks=max(2, compression_workers * 2),
        ):
            codes = np.asarray(chunk)
            if codes.dtype != np.uint8:
                raise ValueError("zstd dosage cache source must emit uint8 codes")
            payload = np.asarray(codes.T, dtype=np.uint8, order="C").tobytes()
            pending.append((start, end, pool.submit(compress, payload)))
            now = time.perf_counter()
            if report_every > 0 and now - last_report >= report_every:
                last_report = now
                done = max(end, 1)
                elapsed = now - started_at
                rate = done / elapsed
                remaining = (total_variants - done) / rate if rate > 0 else 0.0
                print(
                    "  zstd store: %d/%d variants (%.1f%%)  %.1f GB written  "
                    "%.0f variants/s  ~%.0f min left"
                    % (done, total_variants, 100.0 * done / total_variants,
                       position / 1e9, rate, remaining / 60),
                    file=sys.stderr, flush=True)
            # Match the proven BGEN converter: one in-flight frame per worker,
            # plus the next ordered result. A 2x queue only retains another
            # multi-gigabyte wave at cohort scale without adding concurrency.
            if len(pending) > compression_workers:
                write_oldest()
        while pending:
            write_oldest()

    np.savez(
        idx_path,
        # Keep the established TorchGWAS dosage-store index contract. Additional
        # keys are optional metadata and do not alter the core reader contract.
        offs=np.asarray(offsets, dtype=np.int64),
        sizes=np.asarray(sizes, dtype=np.int64),
        rows=np.asarray(rows, dtype=np.int32),
        nsamp=np.int64(source.shape[0]),
        nsnp=np.int64(source.shape[1]),
        chunk=np.int64(chunk_size),
        sub=np.int64(1),
        start=np.int64(0),
        level=np.int64(level),
    )
    with samples_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("FID", "IID"))
        writer.writerows(
            zip(getattr(source, "family_ids", source.sample_ids), source.sample_ids)
        )
    metadata = getattr(source, "variant_metadata", None)
    if metadata is None:
        raise ValueError("zstd cache source must provide chromosome/position/allele metadata")
    # Schema 2 writes the variant table as arrays, not as a TSV. On the
    # benchmark cohort the text table is 466 MB and parsing it cost **6.33 s of
    # a 6.53 s store open** -- 30% of a full scan's wall time, paid on every
    # open, for a format whose compressed genotypes are the smallest read of
    # any we support. A reader-side cache (added first) hides that after the
    # first open; writing arrays means a store never pays it even once, and
    # costs nothing at write time.
    #
    # The TSV is still written, because it is the only human-readable
    # description of what is in a store and it is cheap relative to the
    # genotypes. The `.variants` directory is what the reader loads.
    #
    # Separate `.npy` files in a directory rather than one `.npz`, because the
    # reader memory-maps them. A `.npz` is a zip: `np.load` reads and
    # decompresses every array eagerly, which measured **0.491 s** against
    # **0.018 s** for mmap'd `.npy` — so the first version of this change beat
    # parsing the text (4.474 s) but lost to the reader-side cache it was meant
    # to replace.
    variants_dir = _store_path(output_prefix, ".variants")
    variants_dir.mkdir(parents=True, exist_ok=True)
    for name, array in (
        ("marker_id", np.asarray(source.marker_ids, dtype=str)),
        ("chromosome", np.asarray(metadata["chromosome"], dtype=str)),
        ("position", np.asarray(metadata["position"], dtype=np.int64)),
        ("effect_allele", np.asarray(metadata["effect_allele"], dtype=str)),
        ("other_allele", np.asarray(metadata["other_allele"], dtype=str)),
    ):
        np.save(variants_dir / f"{name}.npy", array, allow_pickle=False)
    with variants_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("chromosome", "marker_id", "position", "effect_allele", "other_allele")
        )
        writer.writerows(
            zip(
                metadata["chromosome"],
                source.marker_ids,
                metadata["position"],
                metadata["effect_allele"],
                metadata["other_allele"],
            )
        )
    if source.shape[1] == 0:
        raise ValueError("genotype conversion excluded every variant")
    raw_bytes = int(source.shape[0] * source.shape[1])
    dosage_scale = float(getattr(source, "dosage_scale", 1.0))
    effect_allele = str(getattr(source, "effect_allele_convention", "unspecified"))
    manifest = {
        # 2 adds `.variants.npz`. Readers accept 1 and fall back to the TSV.
        "cache_schema": 2,
        "backend": "variant-major-uint8-zstd",
        "chunk_size": int(chunk_size),
        "zstd_level": int(level),
        "n_samples": int(source.shape[0]),
        "n_variants": int(source.shape[1]),
        "raw_uint8_bytes": raw_bytes,
        "compressed_bytes": int(position),
        "compression_ratio": float(raw_bytes / position) if position else 0.0,
        "dosage_encoding": "round(clip(expected_allele_count,0,2)*dosage_scale)",
        "dosage_scale": dosage_scale,
        "effect_allele": effect_allele,
        "missing_policy": str(getattr(source, "missing_policy", "unspecified")),
        "exclusion_counts": dict(getattr(source, "exclusion_counts", {})),
        "compression_worker_seconds": float(compression_worker_seconds),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


class ZstdGenotype:
    """Parallel reader for independently framed variant-major uint8 genotypes."""

    ndim = 2
    supports_fused_qc = True
    native_dtype = np.uint8

    def __init__(
        self,
        prefix: str | Path,
        *,
        reader_workers: int = 4,
        prefetch_chunks: int = 8,
        read_batch_bytes: int = 16 << 20,
        read_ahead_batches: int = 2,
        read_workers: int = 1,
        metadata_cache_dir: str | Path | None = None,
    ) -> None:
        from .zstd_native import ZstdInto
        self._decoder = ZstdInto()
        self.read_batch_bytes = int(read_batch_bytes)
        self.read_ahead_batches = int(read_ahead_batches)
        self.read_workers = int(read_workers)
        if min(self.read_batch_bytes, self.read_ahead_batches, self.read_workers) <= 0:
            raise ValueError("read-ahead dimensions must be positive")
        self.prefix = Path(prefix)
        self.zst_path = _store_path(prefix, ".zst")
        self.idx_path = _store_path(prefix, ".idx.npz")
        self.samples_path = _store_path(prefix, ".samples.tsv")
        self.variants_path = _store_path(prefix, ".variants.tsv")
        self.manifest_path = _store_path(prefix, ".complete.json")
        for path in (self.zst_path, self.idx_path, self.samples_path, self.variants_path, self.manifest_path):
            if not path.is_file():
                raise FileNotFoundError(path)

        with np.load(self.idx_path, allow_pickle=False) as index:
            self.offsets = np.asarray(index["offs"], dtype=np.int64)
            self.sizes = np.asarray(index["sizes"], dtype=np.int64)
            self._n_samples = int(index["nsamp"])
            self._n_variants = int(index["nsnp"])
            self.preferred_chunk_size = int(index["chunk"])
            if int(index["sub"]) != 1 or int(index["start"]) != 0:
                raise ValueError("TorchGWAS public zstd reader requires sub=1 and start=0")
            self.rows = (
                np.asarray(index["rows"], dtype=np.int64)
                if "rows" in index.files
                else np.asarray(
                    [
                        min(self.preferred_chunk_size, self._n_variants - frame * self.preferred_chunk_size)
                        for frame in range(len(self.offsets))
                    ],
                    dtype=np.int64,
                )
            )
            self.zstd_level = int(index["level"]) if "level" in index.files else None
        manifest = json.loads(self.manifest_path.read_text())
        self.dosage_scale = float(manifest.get("dosage_scale", 1.0))
        self.native_scale = self.dosage_scale
        if not np.isfinite(self.dosage_scale) or self.dosage_scale <= 0:
            raise ValueError("zstd manifest dosage_scale must be positive")
        expected_frames = (self._n_variants + self.preferred_chunk_size - 1) // self.preferred_chunk_size
        if not (len(self.offsets) == len(self.sizes) == len(self.rows) == expected_frames):
            raise ValueError("zstd index frame count does not match genotype dimensions")
        expected_rows = np.minimum(self.preferred_chunk_size,
                                   self._n_variants - np.arange(expected_frames) * self.preferred_chunk_size)
        if not np.array_equal(self.rows, expected_rows):
            raise ValueError('zstd frame row counts do not match chunk dimensions')
        if np.any(self.sizes <= 0) or (self.offsets.size and self.offsets[0] != 0):
            raise ValueError('zstd frame offsets and sizes are invalid')
        if self.offsets.size > 1 and not np.array_equal(self.offsets[1:], self.offsets[:-1] + self.sizes[:-1]):
            raise ValueError('zstd frame offsets must be contiguous')
        if self.offsets.size and int(self.offsets[-1] + self.sizes[-1]) != self.zst_path.stat().st_size:
            raise ValueError("zstd index does not cover the complete compressed file")

        samples = pd.read_table(self.samples_path, dtype=str)
        if len(samples) != self._n_samples:
            raise ValueError("zstd metadata row counts do not match the index")
        self.family_ids = samples["FID"].to_numpy(dtype=object)
        self.sample_ids = samples["IID"].to_numpy(dtype=object)

        # The variants table is essentially the whole cost of opening a store.
        # Measured on the 8.93M-variant benchmark store: **6.33 s of a 6.53 s
        # open**, a 466 MB TSV, against 0.02 s for the index and 0.05 s for the
        # samples. That is 30% of the store's wall time on a full scan, paid
        # again on every open, and it is why the format with the *smallest*
        # read of the five ranked third on wall time.
        def parse_variants():
            frame = pd.read_table(
                self.variants_path,
                dtype={"chromosome": str, "marker_id": str,
                       "effect_allele": str, "other_allele": str})
            return {
                "marker_id": frame["marker_id"].to_numpy(dtype=object),
                "chromosome": frame["chromosome"].to_numpy(dtype=object),
                "position": frame["position"].to_numpy(dtype=np.int64),
                "effect_allele": frame["effect_allele"].to_numpy(dtype=object),
                "other_allele": frame["other_allele"].to_numpy(dtype=object),
            }

        # Three layouts, newest first. A `.variants` directory of `.npy` files
        # is memory-mapped and costs nothing to open; a single `.variants.npz`
        # (the first cut of schema 2) has to be read and decompressed in full;
        # a schema 1 store has only the TSV and is parsed, with the
        # reader-side cache in front of it.
        variants_dir = _store_path(prefix, ".variants")
        legacy_npz = _store_path(prefix, ".variants.npz")
        if variants_dir.is_dir():
            variants = {
                name: np.load(variants_dir / f"{name}.npy", mmap_mode="r",
                              allow_pickle=False)
                for name in _VARIANT_CACHE_ARRAYS
            }
        elif legacy_npz.exists():
            with np.load(legacy_npz, allow_pickle=False) as handle:
                variants = {name: handle[name] for name in _VARIANT_CACHE_ARRAYS}
        else:
            variants = metadata_cache.cached_arrays(
                metadata_cache_dir, self.variants_path, "zstd-variants",
                _VARIANT_CACHE_ARRAYS, parse_variants)
        if len(variants["marker_id"]) != self._n_variants:
            raise ValueError("zstd metadata row counts do not match the index")
        self.marker_ids = variants["marker_id"]
        self.chromosomes = variants["chromosome"]
        self.positions = np.asarray(variants["position"], dtype=np.int64)
        self.effect_alleles = variants["effect_allele"]
        self.other_alleles = variants["other_allele"]
        self.reader_workers = int(reader_workers)
        self.prefetch_chunks = int(prefetch_chunks)
        if self.reader_workers <= 0 or self.prefetch_chunks <= 0:
            raise ValueError("reader_workers and prefetch_chunks must be positive")
        open_flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
        self._fd = os.open(self.zst_path, open_flags) if hasattr(os, "pread") else None
        self._closed = False

    def __del__(self) -> None:
        self.close()

    def close(self) -> None:
        if getattr(self, "_closed", True):
            return
        fd = getattr(self, "_fd", None)
        if fd is not None:
            try:
                os.close(fd)
            except OSError:
                pass
            self._fd = None
        self._closed = True

    def _read_at(self, size: int, offset: int) -> bytes:
        if self._closed:
            raise ValueError("zstd genotype store is closed")
        if hasattr(os, "pread"):
            return os.pread(self._fd, size, offset)
        # Windows has no os.pread(). Independent short-lived descriptors keep
        # parallel frame reads safe and avoid holding a non-delete-share handle.
        fd = os.open(self.zst_path, os.O_RDONLY | getattr(os, "O_BINARY", 0))
        try:
            os.lseek(fd, offset, os.SEEK_SET)
            chunks: list[bytes] = []
            remaining = size
            while remaining:
                chunk = os.read(fd, remaining)
                if not chunk:
                    break
                chunks.append(chunk)
                remaining -= len(chunk)
            return b"".join(chunks)
        finally:
            os.close(fd)

    @property
    def shape(self) -> tuple[int, int]:
        return self._n_samples, self._n_variants

    @property
    def genotype(self) -> "ZstdGenotype":
        return self

    @property
    def variant_metadata(self) -> dict[str, np.ndarray]:
        return {
            "chromosome": self.chromosomes,
            "position": self.positions,
            "effect_allele": self.effect_alleles,
            "other_allele": self.other_alleles,
        }

    def _decode_frame(self, frame: int) -> np.ndarray:
        blob = self._read_at(int(self.sizes[frame]), int(self.offsets[frame]))
        if len(blob) != int(self.sizes[frame]):
            raise OSError(f"short read for zstd frame {frame}")
        expected = int(self.rows[frame]) * self._n_samples
        raw = np.empty(expected, dtype=np.uint8)
        self._decoder.decompress_into(blob, raw)
        if len(raw) != expected:
            raise OSError(f"zstd frame {frame} decoded to {len(raw)} bytes; expected {expected}")
        return np.frombuffer(raw, dtype=np.uint8).reshape(int(self.rows[frame]), self._n_samples)

    def read_codes(self, start: int, end: int) -> np.ndarray:
        if not (0 <= start <= end <= self._n_variants):
            raise IndexError(f"invalid marker range [{start}, {end})")
        if start == end:
            return np.empty((self._n_samples, 0), dtype=np.uint8)
        first = start // self.preferred_chunk_size
        last = (end - 1) // self.preferred_chunk_size
        decoded = [self._decode_frame(frame) for frame in range(first, last + 1)]
        joined = np.concatenate(decoded, axis=0) if len(decoded) > 1 else decoded[0]
        frame_start = first * self.preferred_chunk_size
        selected = joined[start - frame_start : end - frame_start]
        return np.asarray(selected.T, dtype=np.uint8, order="C")

    def read_chunk(self, start: int, end: int, dtype: np.dtype = np.float32) -> np.ndarray:
        codes = self.read_codes(start, end)
        if np.dtype(dtype) == np.dtype(np.uint8):
            return codes
        output = np.asarray(codes, dtype=dtype, order="C")
        if self.dosage_scale != 1.0:
            output /= self.dosage_scale
        return output

    def __getitem__(self, key) -> np.ndarray:
        if not isinstance(key, tuple) or len(key) != 2:
            raise IndexError("zstd genotype access requires genotype[samples, variants]")
        sample_key, marker_key = key
        if not isinstance(marker_key, slice):
            raise IndexError("zstd marker access must be a contiguous slice")
        start, end, step = marker_key.indices(self._n_variants)
        if step != 1:
            raise IndexError("zstd marker slices must have step=1")
        return self.read_chunk(start, end, np.float32)[sample_key]

    def iter_native_chunks(self, chunk_size, prefetch_chunks=None,
                           reader_workers=None, variant_range=None):
        """Unscaled uint8 transport; consumer divides by native_scale for beta units.

        `variant_range` was added to `iter_chunks` and missed here, which made
        this a **sixth** loader the range had never reached -- and the one that
        matters most, because `PinnedDosageLoader` prefers `iter_native_chunks`
        whenever a source defines it, so the CUDA path took this route and no
        other. The symptom was a `TypeError` turned into "ZstdGenotype cannot
        read a variant range, so it cannot be sharded across devices", which
        reads like a deliberate limitation rather than a missing keyword.
        """
        return self.iter_chunks(chunk_size, dtype=np.uint8,
                                prefetch_chunks=prefetch_chunks,
                                reader_workers=reader_workers,
                                variant_range=variant_range)



    allows_direct_native_fill = True
    """Let `PinnedDosageLoader` fill its pinned ring without an intermediate.

    Without this the loader takes the general path, which is

        np.copyto(self.buffers[index][:count].numpy(), array.T)

    -- a full copy of every decoded byte into pinned memory, *transposed*, so
    the writes are strided as well. On this cohort that is **198.7 GB copied
    with a transpose** for one pass over the store, on top of the decode, and
    it is why the zstd reader sat near six effective cores while the same store
    read through a pinned ring reaches nearly twenty.

    The store is variant-major and so is the pinned buffer, so the direct fill
    needs no transpose at all: `read_into` inflates each frame straight into
    its slice of the destination.
    """

    @contextlib.contextmanager
    def native_reader_session(self):
        """Yield `read_into(start, end, destination)` for the pinned loader.

        `destination` is a `(end - start, n_samples)` uint8 view of a pinned
        buffer. Several of these run at once, one per loader slot, so the
        callable must be thread-safe: it opens its own descriptor and relies on
        `ZstdInto` keeping one native decompression context per calling thread.
        Reads are positional, so the descriptor carries no shared file offset.
        """
        if self._closed:
            raise ValueError('zstd genotype store is closed')
        fd = os.open(self.zst_path, os.O_RDONLY | getattr(os, 'O_BINARY', 0))

        def read_into(start, end, destination):
            rows_per_frame = int(self.preferred_chunk_size)
            first_frame = start // rows_per_frame
            last_frame = (end + rows_per_frame - 1) // rows_per_frame
            width = self._n_samples
            view = memoryview(destination).cast('B')
            for frame in range(first_frame, last_frame):
                frame_start = frame * rows_per_frame
                rows = int(self.rows[frame])
                offset = int(self.offsets[frame])
                size = int(self.sizes[frame])
                blob = bytearray(size)
                got = 0
                while got < size:
                    read = os.preadv(fd, [memoryview(blob)[got:]], offset + got)
                    if not read:
                        raise OSError(
                            f'short zstd read at byte {offset + got}')
                    got += read
                lo = max(start, frame_start)
                hi = min(end, frame_start + rows)
                if hi <= lo:
                    continue
                if lo == frame_start and hi == frame_start + rows:
                    # The common case: the frame lands whole inside this
                    # chunk, so it inflates straight into the destination and
                    # no scratch buffer exists at all.
                    target = view[(lo - start) * width:(hi - start) * width]
                    self._decoder.decompress_into(blob, target)
                else:
                    # A chunk boundary that does not fall on a frame boundary.
                    # Only the first and last frame of a range can do this, so
                    # the scratch costs one frame per chunk at worst, not one
                    # per frame.
                    scratch = np.empty((rows, width), dtype=np.uint8)
                    self._decoder.decompress_into(blob, scratch)
                    destination[lo - start:hi - start] = scratch[
                        lo - frame_start:hi - frame_start]

        try:
            yield read_into
        finally:
            os.close(fd)

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
        variant_range: tuple[int, int] | None = None,
    ):
        """Read contiguous compressed batches, decode each frame once, reblock.

        Decoder futures are bounded independently of read-ahead. Returned arrays
        own their storage and remain valid when subsequent chunks are consumed.

        `variant_range` reads only the frames the range touches. This was the
        fifth loader and the only one the range had never reached -- a ranged
        scan of a zstd store raised `TypeError` rather than producing anything,
        so sharding one across devices was impossible. Frames hold
        `preferred_chunk_size` variants, so a range that does not fall on a
        frame boundary is served by trimming the first and last frames rather
        than by reading whole ones and discarding.
        """
        from .streaming import _resolve_variant_range
        from .zstd_native import contiguous_frame_batches
        if self._closed:
            raise ValueError('zstd genotype store is closed')
        depth = self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks
        workers = self.reader_workers if reader_workers is None else reader_workers
        if min(chunk_size, depth, workers) <= 0:
            raise ValueError('chunk size, prefetch depth and workers must be positive')
        dtype = np.dtype(dtype)
        first, last = _resolve_variant_range(variant_range, self._n_variants)
        rows_per_frame = int(self.preferred_chunk_size)
        frame_first = first // rows_per_frame
        frame_last = (last + rows_per_frame - 1) // rows_per_frame
        wanted = last - first
        batches = contiguous_frame_batches(self.zst_path,
                                            self.offsets[frame_first:frame_last],
                                            self.sizes[frame_first:frame_last],
                                            target_bytes=self.read_batch_bytes,
                                            prefetch_batches=self.read_ahead_batches,
                                            read_workers=self.read_workers)
        pending = deque()

        # A ring of decompression buffers, used only when the consumer copies.
        #
        # The general path allocates a fresh (rows_per_frame x samples) array
        # per frame -- 55.6 MB here -- and first-touch faulting on per-chunk
        # allocation is what cost the BGEN decoder 1.3x. A ring of `depth`
        # buffers removes it, and the reuse is safe by the generator's own
        # ordering: buffer `n % depth` was last written by frame `n - depth`,
        # and by the time frame `n` is submitted below, frame `n - depth` has
        # already been yielded *and* the consumer has copied it, because a
        # generator does not resume until the consumer asks for the next value.
        #
        # It is NOT safe for the uint8 hand-off, where the consumer keeps the
        # array rather than copying it, so the ring is gated on the dtype that
        # forces a copy. `_scaled_copy_path` below is that same condition.
        ring = None
        if dtype != np.dtype(np.uint8):
            ring = [np.empty((rows_per_frame, self._n_samples), dtype=np.uint8)
                    for _ in range(depth)]

        def decode(frame, blob, slot):
            # `frame` indexes the slice handed to the batch reader, so the
            # store-wide index it names is offset by the first frame read.
            absolute = frame_first + frame
            rows = int(self.rows[absolute])
            if ring is None or rows > ring[slot].shape[0]:
                raw = np.empty((rows, self._n_samples), dtype=np.uint8)
            else:
                raw = ring[slot][:rows]
            self._decoder.decompress_into(blob, raw)
            start = absolute * rows_per_frame
            lo = max(0, first - start)
            hi = min(raw.shape[0], last - start)
            return raw if (lo == 0 and hi == raw.shape[0]) else raw[lo:hi]

        def decoded_frames():
            with ThreadPoolExecutor(max_workers=workers, thread_name_prefix='torchgwas-zstd-decode') as pool:
                try:
                    submitted = 0
                    for batch in batches:
                        for frame, blob in batch:
                            pending.append(pool.submit(decode, frame, blob,
                                                       submitted % depth))
                            submitted += 1
                            if len(pending) >= depth:
                                yield pending.popleft().result()
                    while pending:
                        yield pending.popleft().result()
                finally:
                    for future in pending:
                        future.cancel()
                    batches.close()

        # Decompress straight into the chunk that will be yielded, when the
        # geometry allows it.
        #
        # The general path below allocates a fresh (rows_per_frame x samples)
        # array per frame -- 55.6 MB here -- decompresses into it, and then
        # copies it again into the output block. That is one large allocation
        # and one full copy per frame that the aligned case does not need, and
        # first-touch page faulting on per-chunk allocation is what cost the
        # BGEN decoder 3.76x earlier. Tuning I/O and decode concurrency changed
        # nothing (18 configurations within 1.4%), which is what pointed here.
        #
        # Conditions: native uint8 transport, chunk a whole number of frames,
        # and a frame-aligned start, so each frame lands entirely inside one
        # output block. `output[a:b]` stays C-contiguous, which the decoder
        # requires. The general path still handles every other case.
        if (dtype == np.dtype(np.uint8)
                and chunk_size % rows_per_frame == 0
                and first % rows_per_frame == 0):
            def place(blob, whole_rows, target):
                """Decompress one frame into `target`, which may be shorter.

                zstd writes the *whole* frame or fails, so a destination
                clipped to the space left in the output block is refused with
                "Destination buffer is too small" -- which is what happens on
                the final frame whenever the requested range ends mid-frame.
                That case takes a scratch buffer and keeps the rows wanted; it
                is at most one frame per scan, so the allocation this path
                exists to avoid is still avoided for every other frame.
                """
                if target.shape[0] == whole_rows:
                    self._decoder.decompress_into(blob, target)
                    return
                scratch = np.empty((whole_rows, self._n_samples), dtype=np.uint8)
                self._decoder.decompress_into(blob, scratch)
                target[:] = scratch[:target.shape[0]]

            def aligned_chunks():
                with ThreadPoolExecutor(max_workers=workers,
                                        thread_name_prefix='torchgwas-zstd-decode') as pool:
                    inflight = deque()
                    output = None
                    filled = 0
                    produced = 0
                    try:
                        for batch in batches:
                            for frame, blob in batch:
                                if output is None:
                                    height = min(chunk_size, wanted - produced)
                                    output = np.empty((height, self._n_samples),
                                                      dtype=np.uint8)
                                    filled = 0
                                absolute = frame_first + frame
                                whole = int(self.rows[absolute])
                                rows = min(whole, output.shape[0] - filled)
                                target = output[filled:filled + rows]
                                inflight.append(pool.submit(
                                    place, blob, whole, target))
                                filled += rows
                                if filled == output.shape[0]:
                                    for future in inflight:
                                        future.result()
                                    inflight.clear()
                                    yield (first + produced,
                                           first + produced + filled, output.T)
                                    produced += filled
                                    output = None
                        for future in inflight:
                            future.result()
                        if output is not None and filled:
                            yield (first + produced, first + produced + filled,
                                   output[:filled].T)
                            produced += filled
                    finally:
                        for future in inflight:
                            future.cancel()
                        batches.close()
                if produced != wanted:
                    raise ValueError(
                        'zstd decoded rows do not match the requested variant range')

            yield from aligned_chunks()
            return

        frames = decoded_frames()
        emitted = 0
        filled = 0
        output = None
        # Hoisted: the condition depends only on the request, not on the frame,
        # and it decides which of two copy loops runs per block.
        scale_on_copy = (dtype != np.dtype(np.uint8) and self.dosage_scale != 1.0)
        try:
            for raw in frames:
                # Native transport can hand off an aligned, uniquely owned frame
                # directly. Its transpose is a view; pinned staging sees raw's
                # contiguous variant-major storage again via array.T.
                # Positions are absolute variant indices, so a ranged scan
                # reports the variants it actually read; `wanted` bounds the
                # remaining work instead of the whole file.
                if (filled == 0 and dtype == np.dtype(np.uint8)
                        and raw.shape[0] == min(chunk_size, wanted - emitted)):
                    yield first + emitted, first + emitted + raw.shape[0], raw.T
                    emitted += raw.shape[0]
                    continue
                consumed = 0
                while consumed < raw.shape[0]:
                    if output is None:
                        output = np.empty((min(chunk_size, wanted - emitted), self._n_samples), dtype=dtype)
                    count = min(output.shape[0] - filled, raw.shape[0] - consumed)
                    source = raw[consumed:consumed + count]
                    target = output[filled:filled + count]
                    if scale_on_copy:
                        # One pass, not two. Widening the codes into the block
                        # and scaling them used to be separate sweeps, so a
                        # float64 read touched every output element twice --
                        # and the block is 55.6 MB here. `np.divide` with
                        # `out=` widens and divides together.
                        #
                        # Deliberately a divide by `dosage_scale` and not a
                        # multiply by its reciprocal: the reciprocal of a scale
                        # like 254 is not exact in binary, so multiplying would
                        # move the last bits of every dosage and quietly break
                        # the byte-identity checks. Same arithmetic as before,
                        # half the traffic.
                        np.divide(source, self.dosage_scale, out=target)
                    else:
                        target[...] = source
                    consumed += count
                    filled += count
                    if filled == output.shape[0]:
                        yield first + emitted, first + emitted + filled, output.T
                        emitted += filled
                        filled = 0
                        output = None
            if emitted != wanted or filled:
                raise ValueError(
                    'zstd decoded rows do not match the requested variant range')
        finally:
            frames.close()
