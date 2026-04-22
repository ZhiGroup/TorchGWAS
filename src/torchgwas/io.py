from __future__ import annotations

import csv
import gzip
import hashlib
import queue
import shutil
import subprocess
import tempfile
import threading
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from numpy.lib.format import open_memmap
from pandas_plink import Chunk, read_plink1_bin

from .utils import mkdir


class GenotypeChunkLoader:
    def __init__(
        self,
        genotype: np.ndarray,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int = 2,
    ) -> None:
        self.genotype = genotype
        self.chunk_size = int(chunk_size)
        self.dtype = dtype
        self.prefetch_chunks = max(int(prefetch_chunks), 1)
        self.queue: queue.Queue[tuple[int, int, np.ndarray] | None] = queue.Queue(maxsize=self.prefetch_chunks)
        self._thread: threading.Thread | None = None

    def _producer(self) -> None:
        n_markers = self.genotype.shape[1]
        try:
            for start in range(0, n_markers, self.chunk_size):
                end = min(n_markers, start + self.chunk_size)
                chunk = np.asarray(self.genotype[:, start:end], dtype=self.dtype).copy()
                self.queue.put((start, end, chunk))
        finally:
            self.queue.put(None)

    def __iter__(self) -> Iterator[tuple[int, int, np.ndarray]]:
        self._thread = threading.Thread(target=self._producer, daemon=True)
        self._thread.start()
        return self

    def __next__(self) -> tuple[int, int, np.ndarray]:
        item = self.queue.get()
        if item is None:
            if self._thread is not None:
                self._thread.join()
            raise StopIteration
        return item


class DiskBackedGenotype:
    def __init__(
        self,
        memmap_path: str | Path,
        sample_ids: np.ndarray,
        marker_ids: np.ndarray,
        dtype: np.dtype = np.float32,
    ) -> None:
        self.memmap_path = Path(memmap_path)
        self.sample_ids = np.asarray(sample_ids, dtype=object)
        self.marker_ids = np.asarray(marker_ids, dtype=object)
        self.dtype = dtype
        self._genotype = np.load(self.memmap_path, mmap_mode="r")

    @property
    def shape(self) -> tuple[int, int]:
        return tuple(self._genotype.shape)

    @property
    def genotype(self) -> np.ndarray:
        return self._genotype

    def iter_chunks(self, chunk_size: int, dtype: np.dtype = np.float64, prefetch_chunks: int = 2) -> GenotypeChunkLoader:
        return GenotypeChunkLoader(self._genotype, chunk_size=chunk_size, dtype=dtype, prefetch_chunks=prefetch_chunks)


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
    genotype_path = Path(genotype_path)
    if genotype_path.suffix.lower() == ".bed":
        bed = genotype_path
        prefix = genotype_path.with_suffix("")
    else:
        bed = genotype_path.with_suffix(".bed")
        prefix = genotype_path
    bim_path = Path(bim) if bim is not None else prefix.with_suffix(".bim")
    fam_path = Path(fam) if fam is not None else prefix.with_suffix(".fam")
    return bed, bim_path, fam_path


def load_plink_genotype(
    genotype_path: str | Path,
    bim: str | Path | None = None,
    fam: str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bed, bim_path, fam_path = _resolve_plink_triplet(genotype_path, bim=bim, fam=fam)
    data = read_plink1_bin(bed=bed, bim=bim_path, fam=fam_path, verbose=False, chunk=Chunk())
    genotype = np.asarray(data.values, dtype=np.float64)
    marker_ids = np.asarray(data.snp.values, dtype=object)
    sample_ids = np.asarray(data.iid.values.astype(str), dtype=object)
    return genotype, sample_ids, marker_ids


def _find_plink2_binary(explicit_path: str | Path | None = None) -> str:
    candidates = []
    if explicit_path is not None:
        candidates.append(str(explicit_path))
    candidates.extend(
        [
            "plink2",
            "/data4012/zxie3/fast_GWAS/plink2",
            "/data4012/zxie3/Bgen-reader/thirdparty/plink-2.0/plink2",
        ]
    )
    for candidate in candidates:
        if Path(candidate).exists():
            return candidate
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise FileNotFoundError("could not locate a plink2 binary for BGEN fallback decoding")


def _load_bgen_via_plink2(
    genotype_path: str | Path,
    sample_file: str | Path,
    plink2_binary: str | Path | None = None,
    cache_dir: str | Path | None = None,
) -> tuple[DiskBackedGenotype, np.ndarray, np.ndarray]:
    plink2 = _find_plink2_binary(plink2_binary)
    genotype_path = Path(genotype_path)
    sample_file = Path(sample_file)
    cache_prefix = _resolve_bgen_cache_prefix(genotype_path, sample_file, cache_dir=cache_dir)
    cache_array = cache_prefix.with_suffix(".npy")
    sample_ids_path = cache_prefix.with_name(f"{cache_prefix.name}.sample_ids.txt")
    marker_ids_path = cache_prefix.with_name(f"{cache_prefix.name}.marker_ids.txt")
    if cache_array.exists() and sample_ids_path.exists() and marker_ids_path.exists():
        sample_ids = load_vector(sample_ids_path)
        marker_ids = load_vector(marker_ids_path)
        genotype = DiskBackedGenotype(cache_array, sample_ids=sample_ids, marker_ids=marker_ids, dtype=np.float32)
        return genotype, sample_ids, marker_ids
    with tempfile.TemporaryDirectory(prefix="torchgwas_bgen_") as tmpdir:
        out_prefix = Path(tmpdir) / "decoded"
        cmd = [
            plink2,
            "--bgen",
            str(genotype_path),
            "ref-first",
            "--sample",
            str(sample_file),
            "--export",
            "A",
            "--out",
            str(out_prefix),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        genotype, sample_ids, marker_ids = _convert_plink2_raw_to_memmap(
            out_prefix.with_suffix(".raw"),
            cache_array,
            sample_ids_path,
            marker_ids_path,
        )
        return genotype, sample_ids, marker_ids


def load_bgen_genotype(
    genotype_path: str | Path,
    sample_file: str | Path,
    plink2_binary: str | Path | None = None,
    cache_dir: str | Path | None = None,
) -> tuple[DiskBackedGenotype, np.ndarray, np.ndarray]:
    try:
        return _load_bgen_via_plink2(
            genotype_path,
            sample_file,
            plink2_binary=plink2_binary,
            cache_dir=cache_dir,
        )
    except Exception as exc:
        raise RuntimeError(
            "failed to prepare a disk-backed BGEN cache via plink2; "
            "TorchGWAS does not fall back to in-memory BGEN decoding because "
            "that path is unsafe for large cohorts"
        ) from exc


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
    if path.with_suffix(".bed").exists():
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
) -> tuple[np.ndarray | DiskBackedGenotype, np.ndarray | None, np.ndarray | None, dict]:
    resolved_format = infer_genotype_format(genotype_path, genotype_format=genotype_format)
    if resolved_format == "npy":
        genotype = load_array(genotype_path)
        return genotype, None, None, {"genotype_format": "npy"}
    if resolved_format == "plink":
        genotype, sample_ids, marker_ids = load_plink_genotype(genotype_path, bim=bim, fam=fam)
        return genotype, sample_ids, marker_ids, {"genotype_format": "plink"}
    if resolved_format == "bgen":
        if sample_file is None:
            raise ValueError("BGEN input requires --sample-file")
        genotype, sample_ids, marker_ids = load_bgen_genotype(
            genotype_path,
            sample_file=sample_file,
            plink2_binary=plink2_binary,
            cache_dir=genotype_cache_dir,
        )
        meta = {"genotype_format": "bgen", "genotype_backend": "disk_backed"}
        if genotype_cache_dir is not None:
            meta["genotype_cache_dir"] = str(genotype_cache_dir)
        return genotype, sample_ids, marker_ids, meta
    raise ValueError(f"unsupported genotype format: {resolved_format}")


def _resolve_bgen_cache_prefix(
    genotype_path: Path,
    sample_file: Path,
    cache_dir: str | Path | None = None,
) -> Path:
    source = f"{genotype_path.resolve()}::{sample_file.resolve()}::{genotype_path.stat().st_mtime_ns}"
    digest = hashlib.sha1(source.encode("utf-8")).hexdigest()[:12]
    base_dir = mkdir(Path(cache_dir) if cache_dir is not None else Path.cwd() / ".torchgwas_cache")
    return base_dir / f"{genotype_path.stem}_{digest}"


def _write_lines(path: Path, values: np.ndarray) -> None:
    path.write_text("\n".join(str(value) for value in values.tolist()) + "\n")


def _count_text_rows(path: Path) -> int:
    with path.open("r") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def _convert_plink2_raw_to_memmap(
    raw_path: Path,
    cache_array: Path,
    sample_ids_path: Path,
    marker_ids_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    header = pd.read_table(raw_path, sep=r"\s+", nrows=0)
    marker_ids = np.asarray(list(header.columns[6:]), dtype=object)
    n_samples = _count_text_rows(raw_path)
    genotype = open_memmap(cache_array, mode="w+", dtype=np.float32, shape=(n_samples, marker_ids.shape[0]))
    sample_ids = np.empty(n_samples, dtype=object)
    row_offset = 0
    for chunk in pd.read_table(raw_path, sep=r"\s+", chunksize=2048):
        current = len(chunk)
        sample_ids[row_offset : row_offset + current] = chunk["IID"].astype(str).to_numpy(dtype=object)
        genotype[row_offset : row_offset + current, :] = chunk.iloc[:, 6:].to_numpy(dtype=np.float32, copy=False)
        row_offset += current
    genotype.flush()
    _write_lines(sample_ids_path, sample_ids)
    _write_lines(marker_ids_path, marker_ids)
    return DiskBackedGenotype(cache_array, sample_ids=sample_ids, marker_ids=marker_ids, dtype=np.float32), sample_ids, marker_ids


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
    table = table.drop_duplicates(subset=[chosen_id_col]).set_index(chosen_id_col)
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
