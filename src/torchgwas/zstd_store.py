from __future__ import annotations

import json
import os
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import zstandard as zstd

from .streaming import OrderedChunkLoader


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

    def compress(payload: bytes) -> bytes:
        return zstd.ZstdCompressor(level=level).compress(payload)

    pending: deque[tuple[int, int, Future[bytes]]] = deque()
    with zst_path.open("wb", buffering=4 << 20) as output, ThreadPoolExecutor(
        max_workers=compression_workers,
        thread_name_prefix="torchgwas-zstd-encode",
    ) as pool:

        def write_oldest() -> None:
            nonlocal position
            start, end, future = pending.popleft()
            blob = future.result()
            output.write(blob)
            offsets.append(position)
            sizes.append(len(blob))
            rows.append(end - start)
            position += len(blob)

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
            if len(pending) >= compression_workers * 2:
                write_oldest()
        while pending:
            write_oldest()

    np.savez(
        idx_path,
        # Keep the established torchGWAS1.1 zstore index contract. Additional
        # keys are metadata only; lowrank/zstore.py can read this cache directly.
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
    pd.DataFrame(
        {
            "FID": getattr(source, "family_ids", source.sample_ids),
            "IID": source.sample_ids,
        }
    ).to_csv(samples_path, sep="\t", index=False)
    metadata = getattr(source, "variant_metadata", None)
    if metadata is None:
        raise ValueError("zstd cache source must provide chromosome/position/allele metadata")
    pd.DataFrame(
        {
            "chromosome": metadata["chromosome"],
            "marker_id": source.marker_ids,
            "position": metadata["position"],
            "effect_allele": metadata["effect_allele"],
            "other_allele": metadata["other_allele"],
        }
    ).to_csv(variants_path, sep="\t", index=False)
    if source.shape[1] == 0:
        raise ValueError("BGEN conversion excluded every variant")
    raw_bytes = int(source.shape[0] * source.shape[1])
    dosage_scale = float(getattr(source, "dosage_scale", 1.0))
    effect_allele = str(getattr(source, "effect_allele_convention", "unspecified"))
    manifest = {
        "cache_schema": 1,
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
        "exclusion_counts": dict(getattr(source, "exclusion_counts", {})),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


class ZstdGenotype:
    """Parallel reader for independently framed variant-major uint8 genotypes."""

    ndim = 2

    def __init__(
        self,
        prefix: str | Path,
        *,
        reader_workers: int = 4,
        prefetch_chunks: int = 8,
    ) -> None:
        self.prefix = Path(prefix)
        self.zst_path = _store_path(prefix, ".zst")
        self.idx_path = _store_path(prefix, ".idx.npz")
        self.samples_path = _store_path(prefix, ".samples.tsv")
        self.variants_path = _store_path(prefix, ".variants.tsv")
        self.manifest_path = _store_path(prefix, ".complete.json")
        for path in (self.zst_path, self.idx_path, self.samples_path, self.variants_path, self.manifest_path):
            if not path.is_file():
                raise FileNotFoundError(path)

        index = np.load(self.idx_path, allow_pickle=False)
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
        if self.dosage_scale <= 0:
            raise ValueError("zstd manifest dosage_scale must be positive")
        expected_frames = (self._n_variants + self.preferred_chunk_size - 1) // self.preferred_chunk_size
        if not (len(self.offsets) == len(self.sizes) == len(self.rows) == expected_frames):
            raise ValueError("zstd index frame count does not match genotype dimensions")
        if self.offsets.size and int(self.offsets[-1] + self.sizes[-1]) != self.zst_path.stat().st_size:
            raise ValueError("zstd index does not cover the complete compressed file")

        samples = pd.read_table(self.samples_path, dtype=str)
        variants = pd.read_table(self.variants_path, dtype={"chromosome": str, "marker_id": str, "effect_allele": str, "other_allele": str})
        if len(samples) != self._n_samples or len(variants) != self._n_variants:
            raise ValueError("zstd metadata row counts do not match the index")
        self.family_ids = samples["FID"].to_numpy(dtype=object)
        self.sample_ids = samples["IID"].to_numpy(dtype=object)
        self.marker_ids = variants["marker_id"].to_numpy(dtype=object)
        self.chromosomes = variants["chromosome"].to_numpy(dtype=object)
        self.positions = variants["position"].to_numpy(dtype=np.int64)
        self.effect_alleles = variants["effect_allele"].to_numpy(dtype=object)
        self.other_alleles = variants["other_allele"].to_numpy(dtype=object)
        self.reader_workers = int(reader_workers)
        self.prefetch_chunks = int(prefetch_chunks)
        if self.reader_workers <= 0 or self.prefetch_chunks <= 0:
            raise ValueError("reader_workers and prefetch_chunks must be positive")
        self._fd = os.open(self.zst_path, os.O_RDONLY)

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
        blob = os.pread(self._fd, int(self.sizes[frame]), int(self.offsets[frame]))
        if len(blob) != int(self.sizes[frame]):
            raise OSError(f"short read for zstd frame {frame}")
        expected = int(self.rows[frame]) * self._n_samples
        raw = zstd.ZstdDecompressor().decompress(blob, max_output_size=expected)
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

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.float64,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
    ) -> OrderedChunkLoader:
        return OrderedChunkLoader(
            n_markers=self._n_variants,
            read_chunk=self.read_chunk,
            chunk_size=chunk_size,
            dtype=dtype,
            prefetch_chunks=self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks,
            reader_workers=self.reader_workers if reader_workers is None else reader_workers,
        )
