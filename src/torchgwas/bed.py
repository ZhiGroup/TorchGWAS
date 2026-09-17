from __future__ import annotations

import os
import queue
import hashlib
import json
import tempfile
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .streaming import OrderedChunkLoader


_BED_MAGIC = b"\x6c\x1b\x01"
_BIM_CACHE_SCHEMA = 1


def _bim_cache_path(cache_dir: str | Path, bim_path: Path) -> Path:
    stat = bim_path.stat()
    identity = json.dumps(
        {
            "schema": _BIM_CACHE_SCHEMA,
            "path": str(bim_path.resolve()),
            "size": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
        },
        sort_keys=True,
    )
    digest = hashlib.sha1(identity.encode("utf-8")).hexdigest()[:12]
    directory = Path(cache_dir)
    directory.mkdir(parents=True, exist_ok=True)
    return directory / f"{bim_path.stem}_{digest}.bim.cache"


def write_plink_bim_cache(
    cache_dir: str | Path,
    bim_path: str | Path,
    *,
    chromosomes,
    marker_ids,
    positions,
    other_alleles,
    effect_alleles,
) -> Path:
    """Write a validated binary BIM cache for repeated large BED analyses."""

    bim_path = Path(bim_path)
    cache_path = _bim_cache_path(cache_dir, bim_path)
    arrays = {
        "chromosomes": np.asarray(chromosomes, dtype=str),
        "marker_ids": np.asarray(marker_ids, dtype=str),
        "positions": np.asarray(positions, dtype=np.int64),
        "other_alleles": np.asarray(other_alleles, dtype=str),
        "effect_alleles": np.asarray(effect_alleles, dtype=str),
    }
    lengths = {value.shape[0] for value in arrays.values()}
    if len(lengths) != 1:
        raise ValueError("BIM cache arrays have inconsistent lengths")
    if cache_path.is_dir():
        return cache_path
    stat = bim_path.stat()
    directory = cache_path.parent
    with tempfile.TemporaryDirectory(
        prefix=f".{cache_path.name}.",
        dir=directory,
    ) as temporary_name:
        temporary = Path(temporary_name)
        for name, array in arrays.items():
            np.save(temporary / f"{name}.npy", array, allow_pickle=False)
        (temporary / "manifest.json").write_text(
            json.dumps(
                {
                    "schema": _BIM_CACHE_SCHEMA,
                    "bim_path": str(bim_path.resolve()),
                    "bim_size": int(stat.st_size),
                    "bim_mtime_ns": int(stat.st_mtime_ns),
                    "n_variants": next(iter(lengths)),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        temporary.replace(cache_path)
    return cache_path


def _load_plink_bim_cache(cache_dir: str | Path, bim_path: Path):
    cache_path = _bim_cache_path(cache_dir, bim_path)
    if not cache_path.is_dir():
        return None, cache_path
    manifest_path = cache_path / "manifest.json"
    if not manifest_path.is_file():
        return None, cache_path
    manifest = json.loads(manifest_path.read_text())
    stat = bim_path.stat()
    valid = (
        int(manifest.get("schema", -1)) == _BIM_CACHE_SCHEMA
        and str(manifest.get("bim_path")) == str(bim_path.resolve())
        and int(manifest.get("bim_size", -1)) == stat.st_size
        and int(manifest.get("bim_mtime_ns", -1)) == stat.st_mtime_ns
    )
    if not valid:
        return None, cache_path
    arrays = tuple(
        np.load(cache_path / f"{key}.npy", mmap_mode="r", allow_pickle=False)
        for key in (
            "chromosomes",
            "marker_ids",
            "positions",
            "other_alleles",
            "effect_alleles",
        )
    )
    if len({array.shape[0] for array in arrays}) != 1:
        raise ValueError(f"BIM cache arrays have inconsistent lengths: {cache_path}")
    if arrays[0].shape[0] != int(manifest["n_variants"]):
        raise ValueError(f"BIM cache length does not match manifest: {cache_path}")
    return arrays, cache_path


_BIM_COLUMNS = 6
# chromosome, marker id, centimorgans (unused), position, A1, A2
_BIM_WANTED = (0, 1, 3, 4, 5)


def _bim_delimiter(path: Path) -> bytes | None:
    """The single character separating this `.bim`'s fields, if there is one.

    PLINK writes `.bim` with tabs, but space-separated files exist and
    column-aligned ones (runs of spaces) exist too. A fast single-character
    reader can only handle the first two, so this returns None for anything
    else and the caller falls back to the general parser rather than guessing.
    """
    with path.open("rb") as handle:
        head = handle.read(1 << 16)
    line, _, _ = head.partition(b"\n")
    if not line:
        return None
    for candidate in (b"\t", b" "):
        if len(line.split(candidate)) == _BIM_COLUMNS:
            return candidate
    return None


def _read_bim_fast(path: Path):
    """Parse a `.bim` with pandas' C engine and a literal delimiter, or None.

    `sep=r"\\s+"` is a REGEX separator, and pandas falls back to its Python
    engine for it. Passing the actual delimiter keeps it on the C engine for
    the same work and the same output types.

    **The gain is 1.1x to 7.6x, and which one depends on the page cache.**
    Measured on the real 253.9 MB / 8,931,083-row `.bim`, with all five columns
    verified identical against the regex parser across every row:

        condition                           sep=r"\\s+"   literal
        warm, same process                      6.54 s    5.81 s   1.1x
        cold, verified 0.0% resident           33.84 s    4.43 s   7.6x

    The informative part is not the ratio but the spread: the literal-delimiter
    arm sits at 4.4-6.0 s under every condition tried, while the regex arm
    ranges 6.5-33.8 s. pandas' Python engine is what degrades under a cold
    cache and a loaded host; the C engine barely notices. So this is worth
    taking for predictability as much as for speed.

    Do not trust a cold measurement here without checking residency:
    `posix_fadvise(DONTNEED)` is advisory and frequently does not evict.
    Earlier runs of `direct_metadata_parse_check.py` reported 2.2x and 1.8x
    for this same comparison while some arms were silently warm, which is why
    that script now measures residency with `mincore` and flags any arm that
    was not cold.

    **The real fix for this cost is the metadata cache**
    (`_load_plink_bim_cache`), which takes the parse to ~0.17 s on every run
    after the first. It matters because the parse is flat in M and serial, so
    it does not shrink when anything else does: once the hard-call store
    compresses the genotype read to ~1.15 s at full genome, the parse IS the
    first run.

    **pyarrow was tried here and is slower, despite parsing far faster.** A
    microbenchmark had it at 0.49 s against pandas' 13.08 s, which was
    measuring nothing: that arm only read `table.num_rows` and never
    materialised a value. Building the three string columns into numpy object
    arrays means ~27 million Python strings, and with that included the real
    figure was **16.40 s against pandas' 7.23 s** -- 2.3x SLOWER. Arrow holds
    strings in packed buffers and is fast for as long as they stay there; this
    interface hands back object arrays, so they cannot. Reviving it would mean
    changing what the metadata is represented as, all the way through.

    Returns None -- and the caller uses the regex separator -- when the file is
    not simple enough for a single-character delimiter. That fallback is not a
    formality: column-aligned `.bim` files are real.
    """
    delimiter = _bim_delimiter(path)
    if delimiter is None:
        return None
    try:
        table = pd.read_csv(
            path,
            sep=delimiter.decode("ascii"),
            header=None,
            usecols=list(_BIM_WANTED),
            # Declared, never inferred: a marker id of "00123" read as an
            # integer comes back as 123 with the leading zeros gone, and a
            # chromosome column holding both "1" and "X" is the same hazard.
            dtype={0: str, 1: str, 3: np.int64, 4: str, 5: str},
            engine="c",
            memory_map=True,
        )
    except Exception:  # noqa: BLE001 - any parse trouble falls back
        return None
    if table.shape[1] != len(_BIM_WANTED):
        return None
    return (
        table.iloc[:, 0].to_numpy(dtype=object),
        table.iloc[:, 1].to_numpy(dtype=object),
        table[3].to_numpy(dtype=np.int64),
        table[4].to_numpy(dtype=object),
        table[5].to_numpy(dtype=object),
    )


def _make_a2_dosage_lut() -> np.ndarray:
    """Decode PLINK 1 two-bit calls as BIM-column-6 (A2) dosage.

    **This counts A1, and it used to count A2.** A correctness fix, not a
    preference: the PGEN reader emits "the count of ALT1 (PVAR ALT)", and in a
    `.bim` written from a `.pvar` the A1 column *is* ALT, so the two readers
    were counting opposite alleles. Measured on the benchmark cohort before
    the fix, per-variant mean dosage summed to **exactly 2.0000** between the
    two paths on every variant -- `dosage_bed = 2 - dosage_pgen`. The same tool
    therefore reported betas of opposite sign for the same data depending on
    whether it was handed a `.bed` or a `.pgen`, with nothing in the output
    naming the allele.

    A1 is the target because it is what everything else already does: PGEN
    counts ALT, and plink2's `--glm` reports per copy of A1/ALT by default.
    Counting A2 made this reader the only dissenter.

    The GPU kernels carry the same table and are changed with it --
    `decode_two_bit` in `scan_statistics.cu` and `fused_scan.cu`. Changing one
    without the other would be worse than leaving both wrong, because the
    disagreement would then depend on which device the scan happened to pick.
    """

    # PLINK codes, least-significant pair first:
    # PLINK specification: 00=A1/A1, 01=missing, 10=A1/A2, 11=A2/A2.
    # A1 dosage is therefore 2, missing, 1, 0 -- the reverse of the A2 reading.
    code_to_a2 = np.asarray([0.0, np.nan, 1.0, 2.0], dtype=np.float32)
    lut = np.empty((256, 4), dtype=np.float32)
    for byte in range(256):
        for sample_offset in range(4):
            lut[byte, sample_offset] = code_to_a2[(byte >> (2 * sample_offset)) & 0b11]
    return lut


_A2_DOSAGE_LUT = _make_a2_dosage_lut()


def resolve_plink_triplet(
    genotype_path: str | Path,
    bim: str | Path | None = None,
    fam: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    path = Path(genotype_path)
    if path.suffix.lower() in {".bed", ".bim", ".fam"}:
        prefix = Path(str(path)[: -len(path.suffix)])
    else:
        prefix = path
    bed_path = path if path.suffix.lower() == ".bed" else Path(f"{prefix}.bed")
    bim_path = Path(bim) if bim is not None else Path(f"{prefix}.bim")
    fam_path = Path(fam) if fam is not None else Path(f"{prefix}.fam")
    return bed_path, bim_path, fam_path


class PlinkBedGenotype:
    """Bounded-memory, variant-major PLINK BED reader.

    The on-disk format is already variant-major, so each worker receives one
    contiguous variant range and performs one positional read. Decoded chunks
    are returned in the package's public sample-by-variant orientation.
    """

    ndim = 2
    preferred_gpu_chunk_size = 5000

    def __init__(
        self,
        genotype_path: str | Path,
        bim: str | Path | None = None,
        fam: str | Path | None = None,
        reader_workers: int = 4,
        prefetch_chunks: int = 4,
        metadata_cache_dir: str | Path | None = None,
        hardcall_store: str | Path | None = None,
    ) -> None:
        # A hard-call store replaces ONLY the genotype bytes. The `.bim` and
        # `.fam` are still read from the triplet, so the metadata path and its
        # cache are untouched and the store cannot introduce a metadata
        # difference. Frames hold bed-convention two-bit rows, so everything
        # downstream -- the lookup table, the GPU decode kernels, the pinned
        # ring -- is unchanged. Measured 7.636x smaller than the `.bed`, with
        # 414,106 variants verified byte-identical against it.
        #
        # Worth knowing before reaching for it: this pays only where the READ
        # binds. MEASURED at full genome, cold, two interleaved rounds --
        # **2.15x at K = 1** (7.58 s against 16.26 s) and **1.05x at K = 512**
        # (21.63 s against 22.61 s -- inside the round-to-round spread, so read
        # it as parity). At K = 512 the GEMM and the write set a floor, and the
        # store's zstd decode competes with the writer for host CPU where the
        # bed's decode runs on the GPU and does not.
        self._store = None
        if hardcall_store is not None:
            from .hardcall_store import HardcallStore
            self._store = HardcallStore(hardcall_store)

        self.bed_path, self.bim_path, self.fam_path = resolve_plink_triplet(genotype_path, bim=bim, fam=fam)
        # The `.bed` itself need not exist when a store stands in for it --
        # that is the point of having one -- but the metadata still must.
        required = ((self.bim_path, self.fam_path) if self._store is not None
                    else (self.bed_path, self.bim_path, self.fam_path))
        for path in required:
            if not path.is_file():
                raise FileNotFoundError(path)

        fam_table = pd.read_csv(
            self.fam_path,
            sep=r"\s+",
            header=None,
            dtype=str,
            usecols=[0, 1],
            memory_map=True,
        )
        cached_bim = None
        self.metadata_cache_path = None
        if metadata_cache_dir is not None:
            cached_bim, self.metadata_cache_path = _load_plink_bim_cache(
                metadata_cache_dir,
                self.bim_path,
            )
        bim_table = None
        parsed_bim = None
        if cached_bim is None:
            # C engine with a literal delimiter -- 5.81 s against 6.54 s in situ on
            # the 253.9 MB `.bim`, see `_read_bim_fast`. It declines rather
            # than guesses when the file is not simple enough for a
            # single-character separator, so the regex path stays below.
            parsed_bim = _read_bim_fast(self.bim_path)
            if parsed_bim is None:
                bim_table = pd.read_csv(
                    self.bim_path,
                    sep=r"\s+",
                    header=None,
                    usecols=[0, 1, 3, 4, 5],
                    dtype={0: str, 1: str, 3: np.int64, 4: str, 5: str},
                    memory_map=True,
                )
        if fam_table.shape[1] != 2:
            raise ValueError(f"invalid FAM file (expected at least 2 columns): {self.fam_path}")
        if bim_table is not None and bim_table.shape[1] != 5:
            raise ValueError(f"invalid BIM file (expected at least 6 columns): {self.bim_path}")

        self._stored_family_ids = fam_table.iloc[:, 0].to_numpy(dtype=object)
        self._stored_sample_ids = fam_table.iloc[:, 1].to_numpy(dtype=object)
        self.family_ids = self._stored_family_ids
        self.sample_ids = self._stored_sample_ids
        if cached_bim is None:
            if parsed_bim is not None:
                (self.chromosomes, self.marker_ids, self.positions,
                 self.other_alleles,        # A1
                 self.effect_alleles) = parsed_bim   # A2 dosage
            else:
                self.chromosomes = bim_table.iloc[:, 0].to_numpy(dtype=object)
                self.marker_ids = bim_table.iloc[:, 1].to_numpy(dtype=object)
                self.positions = bim_table[3].to_numpy(dtype=np.int64)
                self.other_alleles = bim_table[4].to_numpy(dtype=object)  # A1
                self.effect_alleles = bim_table[5].to_numpy(dtype=object)  # A2
            if metadata_cache_dir is not None:
                self.metadata_cache_path = write_plink_bim_cache(
                    metadata_cache_dir,
                    self.bim_path,
                    chromosomes=self.chromosomes,
                    marker_ids=self.marker_ids,
                    positions=self.positions,
                    other_alleles=self.other_alleles,
                    effect_alleles=self.effect_alleles,
                )
        else:
            (
                self.chromosomes,
                self.marker_ids,
                self.positions,
                self.other_alleles,
                self.effect_alleles,
            ) = cached_bim
        self.reader_workers = int(reader_workers)
        self.prefetch_chunks = int(prefetch_chunks)
        if self.reader_workers <= 0 or self.prefetch_chunks <= 0:
            raise ValueError("reader_workers and prefetch_chunks must be positive")

        self._stored_n_samples = int(self._stored_sample_ids.size)
        self._sample_indices: np.ndarray | None = None
        self._n_samples = self._stored_n_samples
        self._n_markers = int(self.marker_ids.size)
        self._bytes_per_variant = (self._stored_n_samples + 3) // 4
        if self._store is not None:
            # The store's own shape must agree with the metadata, which is the
            # same check the BED size test performs -- a store built from a
            # different cohort or a truncated one would otherwise read
            # plausible genotypes for the wrong samples.
            if self._store.shape != (self._stored_n_samples, self._n_markers):
                raise ValueError(
                    f"hard-call store {self._store.path} holds "
                    f"{self._store.shape} but {self.fam_path} and "
                    f"{self.bim_path} describe "
                    f"{(self._stored_n_samples, self._n_markers)}")
            if self._store.bytes_per_variant != self._bytes_per_variant:
                raise ValueError(
                    f"hard-call store stride {self._store.bytes_per_variant} "
                    f"!= {self._bytes_per_variant}")
            self._fd = None
        else:
            expected_size = 3 + self._n_markers * self._bytes_per_variant
            actual_size = self.bed_path.stat().st_size
            if actual_size != expected_size:
                raise ValueError(
                    f"BED size mismatch for {self.bed_path}: expected {expected_size} bytes "
                    f"for {self._n_samples} samples and {self._n_markers} variants, got {actual_size}"
                )
            with self.bed_path.open("rb") as handle:
                magic = handle.read(3)
            if magic != _BED_MAGIC:
                raise ValueError(
                    f"unsupported BED header {magic!r}; TorchGWAS requires PLINK 1 variant-major BED"
                )
            self._fd = os.open(self.bed_path, os.O_RDONLY)

    @property
    def chunk_alignment_variants(self) -> int | None:
        """Variant multiple a chunk should respect, or None if any size is equal.

        A plain `.bed` is a flat array of fixed-width rows: any range costs
        what it reads, so there is nothing to align to. A hard-call store is
        framed, and the frame is the unit of decompression -- a chunk that
        straddles frames decompresses every frame it touches and discards most
        of it. Measured on an idle H100 at frame 2,048, that costs 2.5x (chunk
        1,024), 4.6x (1,536) and 1.5x (3,072), while every multiple of 2,048
        lands within noise of the best. The planner reads this so it does not
        choose one of the bad sizes by accident.
        """
        store = getattr(self, "_store", None)
        return None if store is None else int(store.frame_variants)

    def read_packed_into(self, target, start: int, end: int) -> int:
        """Fill `target` with packed rows for [start, end), from whichever source.

        The single place the store and the plain `.bed` diverge. Both write
        straight into caller memory -- a pinned ring slot in the hot path --
        so routing through a store adds no per-chunk allocation.
        """
        if self._store is not None:
            return self._store.read_packed_into(target, start, end)
        byte_count = (end - start) * self._bytes_per_variant
        view = memoryview(target).cast("B")[:byte_count]
        offset = 3 + start * self._bytes_per_variant
        received = 0
        while received < byte_count:
            amount = os.preadv(self._fd, [view[received:]], offset + received)
            if amount == 0:
                raise OSError(
                    f"unexpected EOF in {self.bed_path} at byte {offset + received}")
            received += amount
        return received

    def __del__(self) -> None:
        fd = getattr(self, "_fd", None)
        if fd is not None:
            try:
                os.close(fd)
            except OSError:
                pass
            self._fd = None
        store = getattr(self, "_store", None)
        if store is not None:
            try:
                store.close()
            except Exception:  # noqa: BLE001 - interpreter teardown
                pass
            self._store = None

    @property
    def shape(self) -> tuple[int, int]:
        return self._n_samples, self._n_markers

    @property
    def genotype(self) -> "PlinkBedGenotype":
        return self

    @property
    def variant_metadata(self) -> dict[str, np.ndarray]:
        return {
            "chromosome": self.chromosomes,
            "position": self.positions,
            "effect_allele": self.effect_alleles,
            "other_allele": self.other_alleles,
        }

    def select_samples(self, sample_ids) -> "PlinkBedGenotype":
        """Restrict and reorder samples by FAM IID without rewriting the BED."""

        requested = np.asarray([str(value) for value in sample_ids], dtype=object)
        if requested.ndim != 1 or requested.size == 0:
            raise ValueError("sample_ids must be a non-empty one-dimensional sequence")
        if np.unique(requested).size != requested.size:
            raise ValueError("requested sample_ids contain duplicate IID values")
        stored = np.asarray([str(value) for value in self._stored_sample_ids], dtype=object)
        if np.unique(stored).size != stored.size:
            raise ValueError("BED FAM contains duplicate IID values; select samples by FID+IID")
        lookup = {sample_id: index for index, sample_id in enumerate(stored.tolist())}
        missing = [sample_id for sample_id in requested.tolist() if sample_id not in lookup]
        if missing:
            preview = ", ".join(missing[:5])
            raise ValueError(
                f"{len(missing)} requested sample IDs are absent from {self.fam_path} ({preview})"
            )
        indices = np.asarray([lookup[sample_id] for sample_id in requested], dtype=np.int64)
        self._sample_indices = indices
        self.sample_ids = stored[indices]
        self.family_ids = self._stored_family_ids[indices]
        self._n_samples = int(indices.size)
        return self

    def read_chunk(self, start: int, end: int, dtype: np.dtype = np.float32) -> np.ndarray:
        if not (0 <= start <= end <= self._n_markers):
            raise IndexError(f"invalid marker range [{start}, {end}) for {self._n_markers} variants")
        count = end - start
        packed = np.empty((count, self._bytes_per_variant), dtype=np.uint8)
        self.read_packed_into(packed, start, end)
        decoded = _A2_DOSAGE_LUT[packed].reshape(count, self._bytes_per_variant * 4)[
            :, : self._stored_n_samples
        ]
        if self._sample_indices is not None:
            decoded = decoded[:, self._sample_indices]
        return np.asarray(decoded.T, dtype=dtype, order="C")

    def __getitem__(self, key) -> np.ndarray:
        if not isinstance(key, tuple) or len(key) != 2:
            raise IndexError("PLINK BED access requires genotype[samples, variants]")
        sample_key, marker_key = key
        if not isinstance(marker_key, slice):
            raise IndexError("PLINK BED marker access must be a contiguous slice")
        start, end, step = marker_key.indices(self._n_markers)
        if step != 1:
            raise IndexError("PLINK BED marker slices must have step=1")
        chunk = self.read_chunk(start, end, np.float32)
        return chunk[sample_key]

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
        variant_range: tuple[int, int] | None = None,
    ) -> OrderedChunkLoader:
        return OrderedChunkLoader(
            self._n_markers,
            self.read_chunk,
            chunk_size=chunk_size,
            dtype=dtype,
            prefetch_chunks=self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks,
            reader_workers=self.reader_workers if reader_workers is None else reader_workers,
            variant_range=variant_range,
        )

    def iter_packed_chunks(
        self,
        chunk_size: int | None = None,
        reader_workers: int | None = None,
        depth: int | None = None,
        variant_range: tuple[int, int] | None = None,
    ) -> "PinnedPackedBedLoader":
        workers = self.reader_workers if reader_workers is None else int(reader_workers)
        chunk = self.preferred_gpu_chunk_size if chunk_size is None else int(chunk_size)
        ring_depth = max(3, workers + 2) if depth is None else int(depth)
        return PinnedPackedBedLoader(self, chunk, workers, ring_depth,
                                     variant_range=variant_range)


class PinnedPackedBedLoader:
    """Ordered BED reads into pinned buffers without CPU genotype decoding."""

    def __init__(
        self,
        genotype: PlinkBedGenotype,
        chunk_size: int,
        reader_workers: int,
        depth: int,
        variant_range: tuple[int, int] | None = None,
    ) -> None:
        import torch

        from .streaming import _resolve_variant_range

        if chunk_size <= 0 or reader_workers <= 0 or depth <= 0:
            raise ValueError("chunk_size, reader_workers, and depth must be positive")
        self.variant_start, self.variant_end = _resolve_variant_range(
            variant_range, int(genotype.shape[1]))
        self.genotype = genotype
        self.chunk_size = int(chunk_size)
        self.reader_workers = int(reader_workers)
        self.depth = int(depth)
        self.buffers = [
            torch.empty(
                (self.chunk_size, genotype._bytes_per_variant),
                dtype=torch.uint8,
                pin_memory=True,
            )
            for _ in range(self.depth)
        ]
        self.free: queue.Queue[int] = queue.Queue()
        for index in range(self.depth):
            self.free.put(index)
        self._pool = ThreadPoolExecutor(
            max_workers=self.reader_workers,
            thread_name_prefix="torchgwas-bed-packed",
        )

    def _fill(self, buffer_index: int, start: int, end: int):
        byte_count = (end - start) * self.genotype._bytes_per_variant
        target = memoryview(self.buffers[buffer_index].numpy()).cast("B")[:byte_count]
        self.genotype.read_packed_into(target, start, end)
        return buffer_index, start, end

    def __iter__(self):
        bounds = [
            (start, min(self.variant_end, start + self.chunk_size))
            for start in range(self.variant_start, self.variant_end,
                               self.chunk_size)
        ]
        pending: dict[int, Future] = {}
        submit_index = 0

        def submit_one() -> None:
            nonlocal submit_index
            buffer_index = self.free.get()
            start, end = bounds[submit_index]
            pending[submit_index] = self._pool.submit(
                self._fill,
                buffer_index,
                start,
                end,
            )
            submit_index += 1

        while submit_index < min(len(bounds), self.depth):
            submit_one()
        try:
            for output_index in range(len(bounds)):
                buffer_index, start, end = pending.pop(output_index).result()
                yield buffer_index, self.buffers[buffer_index], start, end
                if submit_index < len(bounds):
                    submit_one()
        finally:
            for future in pending.values():
                future.cancel()

    def release(self, buffer_index: int) -> None:
        self.free.put(buffer_index)

    def close(self) -> None:
        self._pool.shutdown(wait=True, cancel_futures=True)
