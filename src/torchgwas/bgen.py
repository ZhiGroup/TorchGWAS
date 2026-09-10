from __future__ import annotations

from collections import deque
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd


DOSAGE_SCALE = 127.5


def _open_bgen(*args, **kwargs):
    try:
        from bgen_reader import open_bgen
    except ImportError as exc:
        raise ImportError(
            "BGEN input requires the optional dependency bgen-reader; "
            "install TorchGWAS with `pip install -e '.[bgen]'`"
        ) from exc
    return open_bgen(*args, **kwargs)


def _sample_columns(sample_file: Path | None, reader_samples: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    reader_samples = np.asarray(reader_samples, dtype=str)
    if sample_file is None:
        ids = reader_samples.astype(object)
        return ids.copy(), ids

    table = pd.read_csv(sample_file, sep=r"\s+", dtype=str)
    if len(table) and table.iloc[0].astype(str).str.fullmatch(r"0").all():
        table = table.iloc[1:].reset_index(drop=True)
    if not {"ID_1", "ID_2"}.issubset(table.columns):
        raise ValueError(f"BGEN sample file must contain ID_1 and ID_2 columns: {sample_file}")
    family_ids = table["ID_1"].to_numpy(dtype=object)
    sample_ids = table["ID_2"].to_numpy(dtype=object)
    if sample_ids.size != reader_samples.size:
        raise ValueError(
            f"sample file has {sample_ids.size} samples but BGEN reports {reader_samples.size}"
        )
    if not np.array_equal(sample_ids.astype(str), reader_samples):
        raise ValueError("sample-file ID_2 order does not match the BGEN reader sample order")
    return family_ids, sample_ids


def _split_biallelic_alleles(value: object) -> tuple[str, str]:
    alleles = str(value).split(",")
    if len(alleles) != 2 or not all(alleles):
        raise ValueError(f"expected two comma-separated BGEN alleles, got {value!r}")
    return alleles[0], alleles[1]


class BgenDosageSource:
    """Stream BGEN probabilities as quantized expected allele-2 dosage.

    Retained values are ``round((P(A1/A2) + 2*P(A2/A2)) * 127.5)``.
    Unsupported or incomplete variants are skipped while decoding, and retained
    variants are repacked into full output chunks.
    """

    dosage_scale = DOSAGE_SCALE
    effect_allele_convention = "BGEN_ALLELE_2"

    def __init__(
        self,
        genotype_path: str | Path,
        sample_file: str | Path | None = None,
        *,
        metadata_path: str | Path | None = None,
        reader_workers: int = 4,
        decode_batch_size: int = 512,
    ) -> None:
        self.genotype_path = Path(genotype_path)
        self.sample_file = None if sample_file is None else Path(sample_file)
        if not self.genotype_path.is_file():
            raise FileNotFoundError(self.genotype_path)
        if self.sample_file is not None and not self.sample_file.is_file():
            raise FileNotFoundError(self.sample_file)
        if reader_workers <= 0 or decode_batch_size <= 0:
            raise ValueError("reader_workers and decode_batch_size must be positive")
        self.reader_workers = int(reader_workers)
        self.decode_batch_size = int(decode_batch_size)
        self._reader = _open_bgen(
            self.genotype_path,
            samples_filepath=self.sample_file,
            metadata_filepath=metadata_path,
            verbose=False,
        )
        self._raw_n_variants = int(self._reader.nvariants)
        self._n_samples = int(self._reader.nsamples)
        self.family_ids, self.sample_ids = _sample_columns(
            self.sample_file, np.asarray(self._reader.samples)
        )

        nalleles = np.asarray(self._reader.nalleles)
        phased = np.asarray(self._reader.phased, dtype=bool)
        if nalleles.shape != (self._raw_n_variants,) or phased.shape != (self._raw_n_variants,):
            raise ValueError("BGEN metadata arrays do not match the reported variant count")
        self._all_ids = np.asarray(self._reader.ids, dtype=object)
        self._all_rsids = np.asarray(self._reader.rsids, dtype=object)
        self._all_chromosomes = np.asarray(self._reader.chromosomes, dtype=object)
        self._all_positions = np.asarray(self._reader.positions, dtype=np.int64)
        self._all_allele_ids = np.asarray(self._reader.allele_ids, dtype=object)
        remaining = np.ones(self._raw_n_variants, dtype=bool)
        self.exclusion_counts = {
            "multiallelic": 0,
            "phased": 0,
            "non_diploid": 0,
            "missing": 0,
            "invalid_probability": 0,
        }
        rejected = remaining & (nalleles != 2)
        self.exclusion_counts["multiallelic"] = int(rejected.sum())
        remaining &= ~rejected
        rejected = remaining & phased
        self.exclusion_counts["phased"] = int(rejected.sum())
        remaining &= ~rejected
        self._metadata_eligible = remaining
        self._kept_indices = np.empty(0, dtype=np.int64)
        self.marker_ids = np.empty(0, dtype=object)
        self.chromosomes = np.empty(0, dtype=object)
        self.positions = np.empty(0, dtype=np.int64)
        self.effect_alleles = np.empty(0, dtype=object)
        self.other_alleles = np.empty(0, dtype=object)
        self._converted = False

    def close(self) -> None:
        reader = getattr(self, "_reader", None)
        if reader is not None:
            reader.close()
            self._reader = None

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    @property
    def shape(self) -> tuple[int, int]:
        n_variants = int(self._kept_indices.size) if self._converted else self._raw_n_variants
        return self._n_samples, n_variants

    @property
    def variant_metadata(self) -> dict[str, np.ndarray]:
        if not self._converted:
            raise RuntimeError("BGEN metadata is final only after the conversion stream is exhausted")
        return {
            "chromosome": self.chromosomes,
            "position": self.positions,
            "effect_allele": self.effect_alleles,
            "other_allele": self.other_alleles,
        }

    def _finalize_metadata(self, kept: list[np.ndarray]) -> None:
        self._kept_indices = (
            np.concatenate(kept).astype(np.int64, copy=False) if kept else np.empty(0, dtype=np.int64)
        )
        ids = self._all_ids[self._kept_indices]
        rsids = self._all_rsids[self._kept_indices]
        rsid_text = rsids.astype(str)
        usable_rsid = ~np.isin(rsid_text, ["", ".", "NA", "nan", "None"])
        self.marker_ids = np.where(usable_rsid, rsids, ids).astype(object)
        self.chromosomes = self._all_chromosomes[self._kept_indices]
        self.positions = self._all_positions[self._kept_indices]
        allele_values = self._all_allele_ids[self._kept_indices]
        split = [_split_biallelic_alleles(value) for value in allele_values]
        self.other_alleles = np.asarray([value[0] for value in split], dtype=object)
        self.effect_alleles = np.asarray([value[1] for value in split], dtype=object)
        self._converted = True

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.uint8,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
    ) -> Iterator[tuple[int, int, np.ndarray]]:
        del prefetch_chunks
        if self._converted:
            raise RuntimeError("a BgenDosageSource conversion stream can only be consumed once")
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if np.dtype(dtype) != np.dtype(np.uint8):
            raise ValueError("BgenDosageSource emits quantized uint8 dosage codes")
        nthreads = self.reader_workers if reader_workers is None else int(reader_workers)
        if nthreads <= 0:
            raise ValueError("reader_workers must be positive")

        code_parts: deque[np.ndarray] = deque()
        index_parts: deque[np.ndarray] = deque()
        buffered = 0
        emitted = 0
        kept: list[np.ndarray] = []

        def take(count: int) -> tuple[np.ndarray, np.ndarray]:
            nonlocal buffered
            codes_out: list[np.ndarray] = []
            indices_out: list[np.ndarray] = []
            remaining = count
            while remaining:
                codes = code_parts[0]
                indices = index_parts[0]
                use = min(remaining, codes.shape[1])
                codes_out.append(codes[:, :use])
                indices_out.append(indices[:use])
                if use == codes.shape[1]:
                    code_parts.popleft()
                    index_parts.popleft()
                else:
                    code_parts[0] = codes[:, use:]
                    index_parts[0] = indices[use:]
                remaining -= use
                buffered -= use
            out_codes = codes_out[0] if len(codes_out) == 1 else np.concatenate(codes_out, axis=1)
            out_indices = indices_out[0] if len(indices_out) == 1 else np.concatenate(indices_out)
            return out_codes, out_indices

        for raw_start in range(0, self._raw_n_variants, self.decode_batch_size):
            raw_end = min(self._raw_n_variants, raw_start + self.decode_batch_size)
            indices = np.flatnonzero(self._metadata_eligible[raw_start:raw_end]) + raw_start
            if not indices.size:
                continue
            probabilities, missing, ploidy = self._reader.read(
                index=(slice(None), indices),
                dtype=np.float16,
                order="C",
                max_combinations=3,
                return_probabilities=True,
                return_missings=True,
                return_ploidies=True,
                num_threads=nthreads,
            )
            probabilities = np.asarray(probabilities)
            missing = np.asarray(missing, dtype=bool)
            ploidy = np.asarray(ploidy)
            if probabilities.shape != (self._n_samples, indices.size, 3):
                raise ValueError(
                    f"unexpected BGEN probability shape {probabilities.shape}; "
                    f"expected {(self._n_samples, indices.size, 3)}"
                )

            has_missing = missing.any(axis=0)
            non_diploid = (~has_missing) & (ploidy != 2).any(axis=0)
            probability_sums = probabilities.astype(np.float32, copy=False).sum(axis=2)
            invalid = (
                (~has_missing)
                & (~non_diploid)
                & (
                    (~np.isfinite(probabilities)).any(axis=(0, 2))
                    | (probabilities < 0).any(axis=(0, 2))
                    | (probabilities > 1).any(axis=(0, 2))
                    | (~np.isclose(probability_sums, 1.0, atol=0.01, rtol=0.0)).any(axis=0)
                )
            )
            valid = ~(has_missing | non_diploid | invalid)
            self.exclusion_counts["missing"] += int(has_missing.sum())
            self.exclusion_counts["non_diploid"] += int(non_diploid.sum())
            self.exclusion_counts["invalid_probability"] += int(invalid.sum())
            if not valid.any():
                continue

            probs = probabilities[:, valid, :].astype(np.float32, copy=False)
            dosage = probs[:, :, 1] + 2.0 * probs[:, :, 2]
            codes = np.rint(np.clip(dosage, 0.0, 2.0) * DOSAGE_SCALE).astype(np.uint8)
            valid_indices = indices[valid]
            code_parts.append(codes)
            index_parts.append(valid_indices)
            buffered += codes.shape[1]

            while buffered >= chunk_size:
                output, output_indices = take(chunk_size)
                kept.append(output_indices)
                yield emitted, emitted + chunk_size, output
                emitted += chunk_size

        if buffered:
            output, output_indices = take(buffered)
            kept.append(output_indices)
            yield emitted, emitted + output.shape[1], output
            emitted += output.shape[1]
        self._finalize_metadata(kept)
