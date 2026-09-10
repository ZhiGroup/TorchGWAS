from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from .streaming import OrderedChunkLoader


_BED_MAGIC = b"\x6c\x1b\x01"


def _make_a2_dosage_lut() -> np.ndarray:
    """Decode PLINK 1 two-bit calls as BIM-column-6 (A2) dosage."""

    # PLINK codes, least-significant pair first:
    # 00=A2/A2, 01=missing, 10=A1/A2, 11=A1/A1.
    code_to_a2 = np.asarray([2.0, np.nan, 1.0, 0.0], dtype=np.float32)
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

    def __init__(
        self,
        genotype_path: str | Path,
        bim: str | Path | None = None,
        fam: str | Path | None = None,
        reader_workers: int = 4,
        prefetch_chunks: int = 4,
    ) -> None:
        self.bed_path, self.bim_path, self.fam_path = resolve_plink_triplet(genotype_path, bim=bim, fam=fam)
        for path in (self.bed_path, self.bim_path, self.fam_path):
            if not path.is_file():
                raise FileNotFoundError(path)

        fam_table = pd.read_csv(self.fam_path, sep=r"\s+", header=None, dtype=str)
        bim_table = pd.read_csv(self.bim_path, sep=r"\s+", header=None, dtype=str)
        if fam_table.shape[1] < 2:
            raise ValueError(f"invalid FAM file (expected at least 2 columns): {self.fam_path}")
        if bim_table.shape[1] < 6:
            raise ValueError(f"invalid BIM file (expected at least 6 columns): {self.bim_path}")

        self.family_ids = fam_table.iloc[:, 0].to_numpy(dtype=object)
        self.sample_ids = fam_table.iloc[:, 1].to_numpy(dtype=object)
        self.chromosomes = bim_table.iloc[:, 0].to_numpy(dtype=object)
        self.marker_ids = bim_table.iloc[:, 1].to_numpy(dtype=object)
        self.positions = pd.to_numeric(bim_table.iloc[:, 3], errors="raise").to_numpy(dtype=np.int64)
        self.other_alleles = bim_table.iloc[:, 4].to_numpy(dtype=object)  # A1
        self.effect_alleles = bim_table.iloc[:, 5].to_numpy(dtype=object)  # A2 dosage
        self.reader_workers = int(reader_workers)
        self.prefetch_chunks = int(prefetch_chunks)
        if self.reader_workers <= 0 or self.prefetch_chunks <= 0:
            raise ValueError("reader_workers and prefetch_chunks must be positive")

        self._n_samples = int(self.sample_ids.size)
        self._n_markers = int(self.marker_ids.size)
        self._bytes_per_variant = (self._n_samples + 3) // 4
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

    def __del__(self) -> None:
        fd = getattr(self, "_fd", None)
        if fd is not None:
            try:
                os.close(fd)
            except OSError:
                pass
            self._fd = None

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

    def _read_exact(self, offset: int, length: int) -> bytes:
        chunks: list[bytes] = []
        received = 0
        while received < length:
            block = os.pread(self._fd, length - received, offset + received)
            if not block:
                raise OSError(f"unexpected EOF in {self.bed_path} at byte {offset + received}")
            chunks.append(block)
            received += len(block)
        return b"".join(chunks)

    def read_chunk(self, start: int, end: int, dtype: np.dtype = np.float32) -> np.ndarray:
        if not (0 <= start <= end <= self._n_markers):
            raise IndexError(f"invalid marker range [{start}, {end}) for {self._n_markers} variants")
        count = end - start
        raw = self._read_exact(3 + start * self._bytes_per_variant, count * self._bytes_per_variant)
        packed = np.frombuffer(raw, dtype=np.uint8).reshape(count, self._bytes_per_variant)
        decoded = _A2_DOSAGE_LUT[packed].reshape(count, self._bytes_per_variant * 4)[:, : self._n_samples]
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
    ) -> OrderedChunkLoader:
        return OrderedChunkLoader(
            self._n_markers,
            self.read_chunk,
            chunk_size=chunk_size,
            dtype=dtype,
            prefetch_chunks=self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks,
            reader_workers=self.reader_workers if reader_workers is None else reader_workers,
        )
