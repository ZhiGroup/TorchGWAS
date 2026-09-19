from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .bed import PlinkBedGenotype, resolve_plink_triplet
from .bgen import BgenDosageSource, BgenGenotype
from .pgen import (
    DEFAULT_PGEN_COMPRESSION_WORKERS,
    DEFAULT_PGEN_DECODE_BATCH_SIZE,
    DEFAULT_PGEN_DECODE_WORKERS,
    PgenDosageSource,
    PgenGenotype,
    pgenlib_version,
    resolve_pgen_triplet,
    sample_selection_identity,
)
from .streaming import ChunkedGenotype, OrderedChunkLoader
from .utils import TEXT_OUTPUT_COMPRESSLEVEL, mkdir
from .zstd_store import (
    DEFAULT_ZSTD_CHUNK_SIZE,
    DEFAULT_ZSTD_LEVEL,
    ZstdGenotype,
    encode_zstd_store,
)


class GenotypeChunkLoader(OrderedChunkLoader):
    def __init__(
        self,
        genotype: np.ndarray,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int = 4,
        reader_workers: int = 4,
    ) -> None:
        self.genotype = genotype

        def read_chunk(start: int, end: int, out_dtype: np.dtype) -> np.ndarray:
            return np.asarray(self.genotype[:, start:end], dtype=out_dtype).copy(order="C")

        super().__init__(
            n_markers=int(genotype.shape[1]),
            read_chunk=read_chunk,
            chunk_size=chunk_size,
            dtype=dtype,
            prefetch_chunks=prefetch_chunks,
            reader_workers=reader_workers,
        )


class DiskBackedGenotype:
    def __init__(
        self,
        memmap_path: str | Path,
        sample_ids: np.ndarray,
        marker_ids: np.ndarray,
        dtype: np.dtype = np.float32,
        storage_order: str = "sample-major",
        reader_workers: int = 4,
        prefetch_chunks: int = 4,
    ) -> None:
        self.memmap_path = Path(memmap_path)
        self.sample_ids = np.asarray(sample_ids, dtype=object)
        self.marker_ids = np.asarray(marker_ids, dtype=object)
        self.dtype = dtype
        if storage_order not in {"sample-major", "variant-major"}:
            raise ValueError("storage_order must be 'sample-major' or 'variant-major'")
        self.storage_order = storage_order
        self.reader_workers = int(reader_workers)
        self.prefetch_chunks = int(prefetch_chunks)
        self._genotype = np.load(self.memmap_path, mmap_mode="r")
        expected_shape = (
            (self.sample_ids.size, self.marker_ids.size)
            if storage_order == "sample-major"
            else (self.marker_ids.size, self.sample_ids.size)
        )
        if self._genotype.shape != expected_shape:
            raise ValueError(
                f"cache shape {self._genotype.shape} does not match metadata {expected_shape} "
                f"for storage_order={storage_order}"
            )

    @property
    def shape(self) -> tuple[int, int]:
        return int(self.sample_ids.size), int(self.marker_ids.size)

    @property
    def genotype(self):
        return self if self.storage_order == "variant-major" else self._genotype

    def __getitem__(self, key) -> np.ndarray:
        if self.storage_order == "sample-major":
            return self._genotype[key]
        if not isinstance(key, tuple) or len(key) != 2:
            raise IndexError("disk-backed access requires genotype[samples, variants]")
        sample_key, marker_key = key
        if not isinstance(marker_key, slice):
            raise IndexError("variant-major marker access must be a contiguous slice")
        return np.asarray(self._genotype[marker_key, sample_key]).T

    def read_chunk(self, start: int, end: int, dtype: np.dtype = np.float32) -> np.ndarray:
        if self.storage_order == "sample-major":
            source = self._genotype[:, start:end]
        else:
            source = self._genotype[start:end, :].T
        return np.asarray(source, dtype=dtype).copy(order="C")

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
    ) -> OrderedChunkLoader:
        return OrderedChunkLoader(
            n_markers=self.shape[1],
            read_chunk=self.read_chunk,
            chunk_size=chunk_size,
            dtype=dtype,
            prefetch_chunks=self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks,
            reader_workers=self.reader_workers if reader_workers is None else reader_workers,
        )


def load_array(path: str | Path) -> np.ndarray:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        return np.load(path, allow_pickle=False)
    if suffix in {".csv", ".tsv", ".txt"}:
        delimiter = "," if suffix == ".csv" else None
        return np.loadtxt(path, delimiter=delimiter)
    raise ValueError(f"unsupported array format for {path}")


def load_vector(path: str | Path | None) -> np.ndarray | None:
    if path is None:
        return None
    path = Path(path)
    values = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return np.asarray(values, dtype=object)


def load_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_table(path, sep=None, engine="python")


def _resolve_plink_triplet(genotype_path: str | Path, bim: str | Path | None = None, fam: str | Path | None = None) -> tuple[Path, Path, Path]:
    return resolve_plink_triplet(genotype_path, bim=bim, fam=fam)


def load_plink_genotype(
    genotype_path: str | Path,
    bim: str | Path | None = None,
    fam: str | Path | None = None,
    reader_workers: int = 4,
    prefetch_chunks: int = 4,
    metadata_cache_dir: str | Path | None = None,
    hardcall_store: str | Path | None = None,
) -> tuple[PlinkBedGenotype, np.ndarray, np.ndarray]:
    genotype = PlinkBedGenotype(
        genotype_path,
        bim=bim,
        fam=fam,
        reader_workers=reader_workers,
        prefetch_chunks=prefetch_chunks,
        metadata_cache_dir=metadata_cache_dir,
        hardcall_store=hardcall_store,
    )
    return genotype, genotype.sample_ids, genotype.marker_ids


def _find_plink2_binary(explicit_path: str | Path | None = None) -> str:
    candidates = []
    if explicit_path is not None:
        candidates.append(str(explicit_path))
    candidates.append("plink2")
    for candidate in candidates:
        if Path(candidate).exists():
            return candidate
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise FileNotFoundError("could not locate a plink2 binary")


def load_bgen_genotype(
    genotype_path: str | Path,
    sample_file: str | Path | None = None,
    plink2_binary: str | Path | None = None,
    cache_dir: str | Path | None = None,
    reader_workers: int = 4,
    prefetch_chunks: int = 4,
) -> tuple[ZstdGenotype, np.ndarray, np.ndarray]:
    """Convert BGEN probabilities directly to the native zstd dosage cache.

    ``plink2_binary`` is retained for API compatibility but is deliberately not
    used: BGEN is decoded directly, quantized to uint8 expected allele dosage,
    and compressed into independent 2,500-variant level-15 frames.
    """

    del plink2_binary
    genotype_path = Path(genotype_path)
    sample_path = None if sample_file is None else Path(sample_file)
    if not genotype_path.is_file():
        raise FileNotFoundError(genotype_path)
    if sample_path is not None and not sample_path.is_file():
        raise FileNotFoundError(sample_path)
    cache_prefix = _resolve_bgen_cache_prefix(genotype_path, sample_path, cache_dir=cache_dir)
    manifest_path = Path(f"{cache_prefix}.complete.json")
    required_suffixes = (".zst", ".idx.npz", ".samples.tsv", ".variants.tsv", ".complete.json")
    # Schema 2 also writes `.variants.npz`. It is published when present but is
    # NOT required, so a schema 1 cache written before the format change stays
    # valid rather than being silently rebuilt. Anything the writer emits that
    # is not listed here is written into the staging directory and then thrown
    # away with it -- which is exactly what happened to the first schema 2
    # store, after a 37-minute conversion.
    optional_suffixes = (".variants", ".variants.npz")
    if manifest_path.is_file() and all(Path(f"{cache_prefix}{suffix}").is_file() for suffix in required_suffixes):
        genotype = ZstdGenotype(
            cache_prefix,
            reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
        )
        return genotype, genotype.sample_ids, genotype.marker_ids

    with tempfile.TemporaryDirectory(prefix=f".{cache_prefix.name}.", dir=cache_prefix.parent) as tmpdir:
        staged_prefix = Path(tmpdir) / "decoded"
        source = BgenDosageSource(
            genotype_path,
            sample_path,
            metadata_path=Path(tmpdir) / "bgen.metadata2.mmm",
            reader_workers=reader_workers,
        )
        try:
            manifest = encode_zstd_store(
                source,
                staged_prefix,
                chunk_size=DEFAULT_ZSTD_CHUNK_SIZE,
                level=DEFAULT_ZSTD_LEVEL,
                compression_workers=reader_workers,
            )
        finally:
            source.close()
        manifest.update(
            {
                "source_bgen": _file_identity(genotype_path),
                "source_sample": None if sample_path is None else _file_identity(sample_path),
            }
        )
        Path(f"{staged_prefix}.complete.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        for suffix in required_suffixes[:-1]:
            Path(f"{staged_prefix}{suffix}").replace(Path(f"{cache_prefix}{suffix}"))
        for suffix in optional_suffixes:
            staged = Path(f"{staged_prefix}{suffix}")
            if staged.exists():
                staged.replace(Path(f"{cache_prefix}{suffix}"))
        Path(f"{staged_prefix}.complete.json").replace(manifest_path)

    genotype = ZstdGenotype(
        cache_prefix,
        reader_workers=reader_workers,
        prefetch_chunks=prefetch_chunks,
    )
    return genotype, genotype.sample_ids, genotype.marker_ids


def load_pgen_genotype(
    genotype_path: str | Path,
    *,
    pvar: str | Path | None = None,
    psam: str | Path | None = None,
    selected_sample_ids=None,
    cache_dir: str | Path | None = None,
    reader_workers: int = 4,
    prefetch_chunks: int = 4,
    pgen_mode: str = "auto",
    decode_workers: int = DEFAULT_PGEN_DECODE_WORKERS,
    decode_batch_size: int = DEFAULT_PGEN_DECODE_BATCH_SIZE,
    zstd_chunk_size: int = DEFAULT_ZSTD_CHUNK_SIZE,
    zstd_level: int = DEFAULT_ZSTD_LEVEL,
    compression_workers: int = DEFAULT_PGEN_COMPRESSION_WORKERS,
) -> tuple[ZstdGenotype, np.ndarray, np.ndarray]:
    """Convert a PLINK 2 PGEN fileset once, then scan the native zstd store."""

    pgen_path, pvar_path, psam_path = resolve_pgen_triplet(
        genotype_path, pvar=pvar, psam=psam
    )
    if pgen_mode not in {"auto", "hardcall", "dosage"}:
        raise ValueError("pgen_mode must be 'auto', 'hardcall', or 'dosage'")
    for name, value in (
        ("reader_workers", reader_workers),
        ("prefetch_chunks", prefetch_chunks),
        ("decode_workers", decode_workers),
        ("decode_batch_size", decode_batch_size),
        ("zstd_chunk_size", zstd_chunk_size),
        ("compression_workers", compression_workers),
    ):
        if int(value) <= 0:
            raise ValueError(f"{name} must be positive")

    selection = sample_selection_identity(selected_sample_ids)
    cache_prefix = _resolve_pgen_cache_prefix(
        pgen_path,
        pvar_path,
        psam_path,
        selection=selection,
        pgen_mode=pgen_mode,
        zstd_chunk_size=zstd_chunk_size,
        zstd_level=zstd_level,
        cache_dir=cache_dir,
    )
    manifest_path = Path(f"{cache_prefix}.complete.json")
    required_suffixes = (".zst", ".idx.npz", ".samples.tsv", ".variants.tsv", ".complete.json")
    # Published when present, never required -- see the note at the BGEN site.
    optional_suffixes = (".variants", ".variants.npz")
    expected_sources = {
        "source_format": "pgen",
        "source_pgen": _file_identity(pgen_path),
        "source_pvar": _file_identity(pvar_path),
        "source_psam": _file_identity(psam_path),
        "sample_selection": selection,
        "pgen_mode": pgen_mode,
    }
    if _valid_pgen_cache(cache_prefix, required_suffixes, expected_sources):
        genotype = ZstdGenotype(
            cache_prefix,
            reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
        )
        return genotype, genotype.sample_ids, genotype.marker_ids

    with tempfile.TemporaryDirectory(prefix=f".{cache_prefix.name}.", dir=cache_prefix.parent) as tmpdir:
        staged_prefix = Path(tmpdir) / "decoded"
        source = PgenDosageSource(
            pgen_path,
            pvar=pvar_path,
            psam=psam_path,
            selected_sample_ids=selected_sample_ids,
            mode=pgen_mode,
            reader_workers=decode_workers,
            decode_batch_size=decode_batch_size,
        )
        try:
            manifest = encode_zstd_store(
                source,
                staged_prefix,
                chunk_size=zstd_chunk_size,
                level=zstd_level,
                compression_workers=compression_workers,
            )
            manifest.update(
                {
                    **expected_sources,
                    "pgenlib_version": pgenlib_version(),
                    "resolved_pgen_mode": source.mode,
                    "pgen_decode_workers": int(decode_workers),
                    "pgen_decode_batch_size": int(decode_batch_size),
                    "pgen_compression_workers": int(compression_workers),
                    "decoded_logical_bytes": int(source.decoded_logical_bytes),
                    "decode_worker_seconds": float(source.decode_worker_seconds),
                    "pgen_read_worker_seconds": float(source.read_worker_seconds),
                    "pgen_transform_worker_seconds": float(
                        source.transform_worker_seconds
                    ),
                    "pgen_metadata_parse_seconds": float(
                        source.metadata_parse_seconds
                    ),
                    "pipeline": (
                        "parallel pgenlib range decode -> in-place uint8 transform -> "
                        "bounded parallel CPU zstd -> ordered write"
                    ),
                }
            )
        finally:
            source.close()
        Path(f"{staged_prefix}.complete.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        for suffix in required_suffixes[:-1]:
            Path(f"{staged_prefix}{suffix}").replace(Path(f"{cache_prefix}{suffix}"))
        for suffix in optional_suffixes:
            staged = Path(f"{staged_prefix}{suffix}")
            if staged.exists():
                staged.replace(Path(f"{cache_prefix}{suffix}"))
        Path(f"{staged_prefix}.complete.json").replace(manifest_path)

    genotype = ZstdGenotype(
        cache_prefix,
        reader_workers=reader_workers,
        prefetch_chunks=prefetch_chunks,
    )
    return genotype, genotype.sample_ids, genotype.marker_ids


def infer_genotype_format(genotype_path: str | Path, genotype_format: str = "auto") -> str:
    if genotype_format != "auto":
        return genotype_format
    path = Path(genotype_path)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        raise ValueError(
            "NumPy genotype files are not supported; use BED, PGEN, or BGEN input"
        )
    if suffix == ".bed":
        return "plink"
    if suffix == ".bgen":
        return "bgen"
    if suffix == ".pgen":
        return "pgen"
    if suffix == ".zst" or Path(f"{path}.zst").exists():
        return "zstd"
    if Path(f"{path}.bed").exists():
        return "plink"
    if Path(f"{path}.pgen").exists():
        return "pgen"
    raise ValueError(f"could not infer genotype format from {path}")


def load_genotype(
    genotype_path: str | Path,
    genotype_format: str = "auto",
    bim: str | Path | None = None,
    fam: str | Path | None = None,
    sample_file: str | Path | None = None,
    pvar: str | Path | None = None,
    psam: str | Path | None = None,
    selected_sample_ids=None,
    genotype_cache_dir: str | Path | None = None,
    plink2_binary: str | Path | None = None,
    reader_workers: int = 4,
    prefetch_chunks: int = 4,
    zstd_read_workers: int | None = None,
    pgen_mode: str = "auto",
    pgen_decode_workers: int | None = None,
    pgen_decode_batch_size: int = DEFAULT_PGEN_DECODE_BATCH_SIZE,
    pgen_compression_workers: int = DEFAULT_PGEN_COMPRESSION_WORKERS,
    bgen_decode_backend: str = "auto",
    hardcall_store: str | Path | None = None,
) -> tuple[np.ndarray | ChunkedGenotype, np.ndarray | None, np.ndarray | None, dict]:
    resolved_format = infer_genotype_format(genotype_path, genotype_format=genotype_format)
    if hardcall_store is not None and resolved_format != "plink":
        # REFUSE rather than ignore. Only the PLINK path can substitute a
        # hard-call store for its genotype bytes, and silently dropping the
        # argument for any other format would report a "store" measurement
        # taken without the store -- a benchmark that lies in the flattering
        # direction. This module has been bitten by exactly that before:
        # `selected_sample_ids` was accepted and dropped, and a sweep asking
        # for 4,000 of 35,365 samples got full-cohort scans at every size,
        # with identical times and identical peak memory to prove it.
        raise ValueError(
            f"hardcall_store is only supported for the PLINK format, not "
            f"{resolved_format!r}; the store holds bed-convention two-bit "
            f"rows and stands in for a .bed")
    if resolved_format == "zstd":
        prefix = str(genotype_path)
        if prefix.endswith(".zst"):
            prefix = prefix[:-4]
        # `read_workers` is the store's I/O lane count, a separate knob from
        # `reader_workers`, which counts decoder threads. It used not to be
        # reachable from here at all, so every store opened through the public
        # loader ran on the constructor default of **one** lane -- and that is
        # not a detail. Measured on this cohort, 2M variants at K=128 with 16
        # decoders: **25.38 s at one lane against 15.88 s at eight**, a 1.6x
        # left on the floor by a parameter the reader documents and the loader
        # did not forward. Default 8, where the sweep flattens.
        if selected_sample_ids is not None:
            # Refuse rather than ignore. The store holds variant-major rows in
            # file order and has no column-selection path, so a subset cannot
            # be honoured here -- but accepting one and scanning the whole
            # cohort anyway is worse than saying so: a sweep asking for 4,000
            # of 35,365 samples silently got all 35,365, at every size.
            raise ValueError(
                "the zstd store cannot subset samples: its rows are stored in "
                "file order with no column selection. Build a store for the "
                "cohort you want, or use a format that subsets (PGEN, BGEN, "
                "PLINK BED)")
        genotype = ZstdGenotype(
            prefix, reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
            read_workers=(8 if zstd_read_workers is None
                          else int(zstd_read_workers)),
            metadata_cache_dir=genotype_cache_dir)
        return genotype, genotype.sample_ids, genotype.marker_ids, {
            "genotype_format": "zstd", "genotype_backend": "native_zstd_store",
        }
    if resolved_format == "plink":
        genotype, sample_ids, marker_ids = load_plink_genotype(
            genotype_path,
            bim=bim,
            fam=fam,
            reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
            metadata_cache_dir=genotype_cache_dir,
            hardcall_store=hardcall_store,
        )
        if selected_sample_ids is not None:
            # This used to be accepted and dropped. A sweep asking for 4,000 of
            # 35,365 samples got a full-cohort scan at every requested size --
            # identical time and identical peak memory, which is what a
            # discarded subset looks like. `PlinkBedGenotype.select_samples`
            # existed the whole time; the loader simply never called it.
            genotype = genotype.select_samples(selected_sample_ids)
            sample_ids = genotype.sample_ids
        return genotype, sample_ids, marker_ids, {
            "genotype_format": "plink",
            "genotype_backend": "direct_variant_major_bed",
            "reader_workers": int(reader_workers),
            "prefetch_chunks": int(prefetch_chunks),
            "effect_allele": "BIM_A2",
            "metadata_cache_path": (
                None
                if genotype.metadata_cache_path is None
                else str(genotype.metadata_cache_path)
            ),
        }
    if resolved_format == "bgen":
        genotype = BgenGenotype(genotype_path, sample_file=sample_file,
            decode_backend=bgen_decode_backend, reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
            metadata_cache_dir=genotype_cache_dir)
        if selected_sample_ids is not None:
            genotype.select_samples(selected_sample_ids)
        return genotype, genotype.sample_ids, genotype.marker_ids, {
            "genotype_format": "bgen", "genotype_backend": "direct_bgen",
            "decode_backend_requested": bgen_decode_backend,
            "dosage_scale": 1.0, "effect_allele": "BGEN_ALLELE_2",
            "cache_conversion": False,
        }
    if resolved_format == "pgen":
        # PgenGenotype fixes its reader count at construction, and
        # native_scan prefers source.decode_workers over the argument it is
        # given, so a reader session built for one worker cannot be widened
        # later. Following reader_workers here is what makes --reader-workers
        # mean anything for PGEN; an explicit --pgen-decode-workers still wins.
        resolved_pgen_workers = (
            reader_workers if pgen_decode_workers is None else pgen_decode_workers
        )
        genotype = PgenGenotype(
            genotype_path, pvar=pvar, psam=psam if psam is not None else sample_file,
            selected_sample_ids=selected_sample_ids, mode=pgen_mode,
            reader_workers=resolved_pgen_workers,
            decode_batch_size=pgen_decode_batch_size,
            prefetch_chunks=prefetch_chunks,
            metadata_cache_dir=genotype_cache_dir,
        )
        return genotype, genotype.sample_ids, genotype.marker_ids, {
            "genotype_format": "pgen", "genotype_backend": "direct_pgen",
            "pgen_mode": pgen_mode, "resolved_pgen_mode": genotype.mode,
            "reader_workers": resolved_pgen_workers,
            "pgen_decode_workers_requested": pgen_decode_workers,
            "dosage_scale": 1.0,
            "effect_allele": "PVAR_ALT1", "cache_conversion": False,
            "missing_policy": genotype.missing_policy,
        }
    raise ValueError(f"unsupported genotype format: {resolved_format}")


def _resolve_bgen_cache_prefix(
    genotype_path: Path,
    sample_file: Path | None,
    cache_dir: str | Path | None = None,
) -> Path:
    source = json.dumps(
        {
            "cache_schema": 3,
            "bgen": _file_identity(genotype_path),
            "sample": None if sample_file is None else _file_identity(sample_file),
            "dosage_encoding": "round(expected_BGEN_allele_2_dosage*127.5)",
            "zstd_level": DEFAULT_ZSTD_LEVEL,
            "zstd_frame_variants": DEFAULT_ZSTD_CHUNK_SIZE,
        },
        sort_keys=True,
    )
    digest = hashlib.sha1(source.encode("utf-8")).hexdigest()[:12]
    base_dir = mkdir(Path(cache_dir) if cache_dir is not None else Path.cwd() / ".torchgwas_cache")
    return base_dir / f"{genotype_path.stem}_{digest}"


def _resolve_pgen_cache_prefix(
    pgen_path: Path,
    pvar_path: Path,
    psam_path: Path,
    *,
    selection: dict | None,
    pgen_mode: str,
    zstd_chunk_size: int,
    zstd_level: int,
    cache_dir: str | Path | None = None,
) -> Path:
    source = json.dumps(
        {
            "cache_schema": 1,
            "source_format": "pgen",
            "pgen": _file_identity(pgen_path),
            "pvar": _file_identity(pvar_path),
            "psam": _file_identity(psam_path),
            "sample_selection": selection,
            "pgen_mode": pgen_mode,
            "dosage_encoding": "round(PVAR_ALT1_count*127.5)",
            "missing_policy": "exclude_variant_if_any_selected_sample_is_missing",
            "zstd_level": int(zstd_level),
            "zstd_frame_variants": int(zstd_chunk_size),
        },
        sort_keys=True,
    )
    digest = hashlib.sha1(source.encode("utf-8")).hexdigest()[:12]
    base_dir = mkdir(Path(cache_dir) if cache_dir is not None else Path.cwd() / ".torchgwas_cache")
    return base_dir / f"{pgen_path.stem}_pgen_{digest}"


def _valid_pgen_cache(prefix: Path, required_suffixes: tuple[str, ...], expected: dict) -> bool:
    if not all(Path(f"{prefix}{suffix}").is_file() for suffix in required_suffixes):
        return False
    try:
        manifest = json.loads(Path(f"{prefix}.complete.json").read_text())
    except (OSError, ValueError, TypeError):
        return False
    return all(manifest.get(key) == value for key, value in expected.items())


def _file_identity(path: Path) -> dict[str, str | int]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def align_table_to_samples(
    table_path: str | Path,
    sample_ids: np.ndarray,
    value_columns: list[str] | None = None,
    sample_id_column: str = "IID",
    fallback_sample_id_columns: tuple[str, ...] = ("sampleid", "sample_id", "ID", "id"),
) -> tuple[np.ndarray, list[str]]:
    table = load_table(table_path)
    chosen_id_col = sample_id_column
    if chosen_id_col not in table.columns:
        for candidate in fallback_sample_id_columns:
            if candidate in table.columns:
                chosen_id_col = candidate
                break
        else:
            raise ValueError(f"could not find sample ID column in {table_path}; looked for {sample_id_column} and {fallback_sample_id_columns}")
    table[chosen_id_col] = table[chosen_id_col].astype(str)
    if value_columns is None:
        excluded = {"FID", "fid", chosen_id_col}
        value_columns = [col for col in table.columns if col not in excluded]
    missing_cols = [col for col in value_columns if col not in table.columns]
    if missing_cols:
        raise ValueError(f"missing requested columns in {table_path}: {missing_cols}")
    duplicate_ids = table.loc[table[chosen_id_col].duplicated(keep=False), chosen_id_col].unique()
    if duplicate_ids.size:
        preview = ", ".join(str(value) for value in duplicate_ids[:5])
        raise ValueError(
            f"{table_path} contains {duplicate_ids.size} duplicated {chosen_id_col} values "
            f"({preview}); use unique IID values or pre-align by FID+IID"
        )
    table = table.set_index(chosen_id_col)
    missing_ids = [sample_id for sample_id in sample_ids if sample_id not in table.index]
    if missing_ids:
        raise ValueError(f"{table_path} is missing {len(missing_ids)} samples present in genotype input")
    aligned = table.loc[list(sample_ids), value_columns].to_numpy(dtype=np.float64)
    return aligned, value_columns


def write_table(rows: list[dict], path: str | Path) -> None:
    path = Path(path)
    mkdir(path.parent)
    if not rows:
        raise ValueError("refusing to write an empty result table")
    with gzip.open(path, "wt", newline="", compresslevel=TEXT_OUTPUT_COMPRESSLEVEL) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
