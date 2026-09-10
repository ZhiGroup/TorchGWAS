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
from .bgen import BgenDosageSource
from .streaming import ChunkedGenotype, OrderedChunkLoader
from .utils import mkdir
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
) -> tuple[PlinkBedGenotype, np.ndarray, np.ndarray]:
    genotype = PlinkBedGenotype(
        genotype_path,
        bim=bim,
        fam=fam,
        reader_workers=reader_workers,
        prefetch_chunks=prefetch_chunks,
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
        return "npy"
    if suffix == ".bed":
        return "plink"
    if suffix == ".bgen":
        return "bgen"
    if Path(f"{path}.bed").exists():
        return "plink"
    raise ValueError(f"could not infer genotype format from {path}")


def load_genotype(
    genotype_path: str | Path,
    genotype_format: str = "auto",
    bim: str | Path | None = None,
    fam: str | Path | None = None,
    sample_file: str | Path | None = None,
    genotype_cache_dir: str | Path | None = None,
    plink2_binary: str | Path | None = None,
    reader_workers: int = 4,
    prefetch_chunks: int = 4,
) -> tuple[np.ndarray | ChunkedGenotype, np.ndarray | None, np.ndarray | None, dict]:
    resolved_format = infer_genotype_format(genotype_path, genotype_format=genotype_format)
    if resolved_format == "npy":
        genotype = load_array(genotype_path)
        return genotype, None, None, {"genotype_format": "npy"}
    if resolved_format == "plink":
        genotype, sample_ids, marker_ids = load_plink_genotype(
            genotype_path,
            bim=bim,
            fam=fam,
            reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
        )
        return genotype, sample_ids, marker_ids, {
            "genotype_format": "plink",
            "genotype_backend": "direct_variant_major_bed",
            "reader_workers": int(reader_workers),
            "prefetch_chunks": int(prefetch_chunks),
            "effect_allele": "BIM_A2",
        }
    if resolved_format == "bgen":
        genotype, sample_ids, marker_ids = load_bgen_genotype(
            genotype_path,
            sample_file=sample_file,
            plink2_binary=plink2_binary,
            cache_dir=genotype_cache_dir,
            reader_workers=reader_workers,
            prefetch_chunks=prefetch_chunks,
        )
        meta = {
            "genotype_format": "bgen",
            "genotype_backend": "direct_bgen_uint8_zstd_cache",
            "reader_workers": int(reader_workers),
            "prefetch_chunks": int(prefetch_chunks),
            "zstd_level": DEFAULT_ZSTD_LEVEL,
            "zstd_frame_variants": DEFAULT_ZSTD_CHUNK_SIZE,
            "dosage_scale": 127.5,
            "effect_allele": "BGEN_ALLELE_2",
        }
        if genotype_cache_dir is not None:
            meta["genotype_cache_dir"] = str(genotype_cache_dir)
        return genotype, sample_ids, marker_ids, meta
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
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
