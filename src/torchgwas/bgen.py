from __future__ import annotations

from collections import deque
from collections.abc import Iterator
from pathlib import Path

import os
import threading

import numpy as np
import pandas as pd

from .bgen_cpu import (BgenDecodeError, CpuBgenDecoder, cpu_decoder_available,
                       transpose as _native_transpose)


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
            "masked_missing": 0,
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

            # A missing sample carries no meaningful ploidy or probabilities,
            # so every check below is made over observed samples only. The
            # variant stays; the missing sample is masked further down. The GPU
            # decoder already behaves this way, and the two disagreeing on the
            # same file was the defect.
            observed = ~missing
            has_missing = missing.any(axis=0)
            non_diploid = ((ploidy != 2) & observed).any(axis=0)
            probability_sums = probabilities.astype(np.float32, copy=False).sum(axis=2)
            bad_cell = (
                (~np.isfinite(probabilities)).any(axis=2)
                | (probabilities < 0).any(axis=2)
                | (probabilities > 1).any(axis=2)
                | (~np.isclose(probability_sums, 1.0, atol=0.01, rtol=0.0))
            )
            invalid = (~non_diploid) & (bad_cell & observed).any(axis=0)
            valid = ~(non_diploid | invalid | ~observed.any(axis=0))
            self.exclusion_counts["masked_missing"] += int(has_missing.sum())
            self.exclusion_counts["non_diploid"] += int(non_diploid.sum())
            self.exclusion_counts["invalid_probability"] += int(invalid.sum())
            if not valid.any():
                continue

            probs = probabilities[:, valid, :].astype(np.float32, copy=False)
            dosage = probs[:, :, 1] + 2.0 * probs[:, :, 2]
            seen = observed[:, valid]
            if not seen.all():
                # Mask: a missing sample takes the variant's observed mean, so
                # its centred contribution is zero, matching the GPU path.
                counts = seen.sum(axis=0)
                sums = np.where(seen, dosage, 0.0).sum(axis=0, dtype=np.float64)
                means = np.divide(sums, counts, out=np.zeros_like(sums),
                                  where=counts > 0)
                dosage = np.where(seen, dosage, means[None, :].astype(dosage.dtype))
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


_BGI_CACHE_SCHEMA = 1
_BGI_CACHE_ARRAYS = ("chromosomes", "positions", "marker_ids",
                     "other_alleles", "effect_alleles", "offsets", "lengths")


def _bgi_cache_path(cache_dir, bgi_path):
    import hashlib

    digest = hashlib.sha256(str(Path(bgi_path).resolve()).encode()).hexdigest()
    return Path(cache_dir) / f"bgi-{digest[:32]}"


def write_bgen_index_cache(cache_dir, bgi_path, **arrays):
    """Store the parsed `.bgi` so the next open does not re-read SQLite.

    Parsing the index is the single largest item in a BGEN run: **22.7 s of a
    59.6 s full-file scan, 38%**, and it is paid before a variant is read.
    Measured against the alternatives, it is also the *only* thing that can be
    fixed here -- reading the whole 391 MB index off cold disk takes **0.2 s**,
    so this is not I/O; SQLite's own `mmap_size` was tried and made no
    difference (39.96 s against 39.50 s); and converting the four `object`
    columns to fixed-width text made it **slower** (15.0 s against 13.9 s).
    What remains is that SQLite materialises 8.93M rows as Python tuples, which
    costs ~10 s of the ~14 s whatever is done with them afterwards. The only
    way not to pay it is not to do it twice.

    Same shape as the PLINK BIM cache next door: a directory of `.npy` files
    and a manifest keyed on the index's size and mtime, written through a
    temporary directory so a crash mid-write cannot leave a half-cache that a
    later run would trust.
    """
    import json
    import tempfile

    bgi_path = Path(bgi_path)
    cache_path = _bgi_cache_path(cache_dir, bgi_path)
    #  refuses object arrays without pickle, and that refusal is
    # correct: a cache directory that can execute code when it is loaded is a
    # hole, not a feature. The text columns arrive from SQLite as , so
    # they are stored as fixed-width unicode and restored on load.
    # Stored as numpy unicode, and the size that costs is deliberate. UCS-4
    # makes a  marker id 96 bytes for 15 ASCII characters, so the cache
    # is **928.8 MB against a 391 MB index** -- larger than the thing it
    # caches. Storing  bytes instead does shrink it to 393 MB, and was
    # tried: the load then needs  over 35.7M strings, and the
    # cache hit goes from **1.3 s to 11.7 s**, collapsing the speedup from
    # 14.9x to 1.6x. The whole point of this file is the open time, so it buys
    # back the space with the only currency that does not matter here: an
    # 8.5 TB volume holding 900 MB per BGEN.
    text = {"chromosomes", "marker_ids", "other_alleles", "effect_alleles"}
    stored = {name: np.asarray(arrays[name],
                               dtype=str if name in text else None)
              for name in _BGI_CACHE_ARRAYS}
    lengths = {value.shape[0] for value in stored.values()}
    if len(lengths) != 1:
        raise ValueError("BGI cache arrays have inconsistent lengths")
    if cache_path.is_dir():
        return cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    stat = bgi_path.stat()
    with tempfile.TemporaryDirectory(prefix=f".{cache_path.name}.",
                                     dir=cache_path.parent) as temporary_name:
        temporary = Path(temporary_name)
        for name, array in stored.items():
            np.save(temporary / f"{name}.npy", array, allow_pickle=False)
        (temporary / "manifest.json").write_text(json.dumps({
            "schema": _BGI_CACHE_SCHEMA,
            "bgi_path": str(bgi_path.resolve()),
            "bgi_size": stat.st_size,
            "bgi_mtime_ns": stat.st_mtime_ns,
            "variants": int(next(iter(lengths))),
        }, indent=2, sort_keys=True) + "\n")
        try:
            os.replace(temporary, cache_path)
        except OSError:
            # Another process won the race; its cache is as good as ours.
            return cache_path
    return cache_path


def _load_bgen_index_cache(cache_dir, bgi_path):
    """Return the cached arrays, or None if there is no trustworthy cache."""
    import json

    cache_path = _bgi_cache_path(cache_dir, bgi_path)
    manifest_path = cache_path / "manifest.json"
    if not manifest_path.is_file():
        return None, cache_path
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, ValueError):
        return None, cache_path
    stat = Path(bgi_path).stat()
    # An index that changed size or mtime is a different index. Trusting a
    # stale cache here would silently scan the wrong file offsets, which is
    # the one failure mode a metadata cache must not have.
    if not (int(manifest.get("schema", -1)) == _BGI_CACHE_SCHEMA
            and str(manifest.get("bgi_path")) == str(Path(bgi_path).resolve())
            and int(manifest.get("bgi_size", -1)) == stat.st_size
            and int(manifest.get("bgi_mtime_ns", -1)) == stat.st_mtime_ns):
        return None, cache_path
    try:
        arrays = {name: np.load(cache_path / f"{name}.npy", allow_pickle=False)
                  for name in _BGI_CACHE_ARRAYS}
        # Hand back exactly what the SQLite path produces, so nothing
        # downstream can tell which route the metadata took.
        for name in ("chromosomes", "marker_ids", "other_alleles",
                     "effect_alleles"):
            arrays[name] = arrays[name].astype(object)
    except (OSError, ValueError):
        return None, cache_path
    if len({value.shape[0] for value in arrays.values()}) != 1:
        return None, cache_path
    return arrays, cache_path


class BgenGenotype:
    """Direct, precision-preserving BGEN Layout-2 dosage source.

    **Supported scope**, decoded exactly and pinned by a test grid:

    - **Layout 2 only.** A Layout-1 file is refused at open, by name, with the
      advice to convert it. Layout 1 is a different record encoding, not a
      variation on this one.
    - **Compression 0 (none), 1 (zlib) and 2 (zstd).** Flag 3 is reserved by the
      specification and is refused at open, also by name. The two refusals are
      separate messages because they have different answers.
    - **Biallelic, unphased, diploid records, bit depth 1 to 32.** Anything else
      raises at the record rather than silently changing the marker axis.

    The three decoders do not all cover the same codecs, and that is deliberate
    rather than an oversight: the C and GPU decoders handle zlib (the C one also
    handles uncompressed), while zstd files fall back to the numpy reference
    path. A file decoding by a different route must not decode to different
    numbers, so the narrower backends decline rather than approximate.

    Physical file order is used for sequential I/O. BGI metadata is read once;
    its initialization time must be reported separately from the scan.
    """
    supports_fused_qc = True
    cpu_parse_worker_parameter = "reader_workers"
    decode_tile_multiple_of_chunk = True
    dosage_scale = 1.0
    effect_allele_convention = 'BGEN_ALLELE_2'

    # The two decode backends want different reader-worker counts, and the
    # difference is large enough that one default cannot serve both. Measured on
    # the 86 GB UK Biobank BGEN.
    #
    #   CPU decode  rises to a plateau at 12-16 workers and stays there. On a
    #               quiet host (`busy_percent` 17-42%): 15,264 variants/s at 8,
    #               21,187 at 12, 22,744 at 16, 22,887 at 32, 22,166 at 48 --
    #               monotone, 9.98x scaling, no cliff anywhere.
    #
    #   GPU decode  improves to 16 workers and is flat from there to 44 (233 /
    #               139 / 113 / 118 seconds at 4 / 8 / 16 / 44). The host
    #               threads only read and hand off.
    #
    # **This default was 8, chosen from a measurement that does not reproduce.**
    # That sweep reported a peak at 8 (6,606 variants/s) *halving* by 16, and
    # the quiet-host curve above shows no such collapse -- 16 is 49% faster than
    # 8, not half. Re-running the old sweep on its original host reproduced
    # nothing either: it peaked at 12, collapsed 7x across 24-48, and then
    # *recovered to its best at 96*, which no property of the decoder can
    # produce and which identifies the host rather than the code. Part of the
    # original effect may also have been real and since removed: the collapse
    # signature was more CPU for less output, which is allocation thrashing, and
    # the per-chunk staging buffer that caused it is now reused per thread.
    #
    # 12 rather than 16, deliberately. Both curves put 12 within a few percent
    # of peak, and 12 sits further from the region where the contended host fell
    # over, so it gives up ~7% on an idle machine to be robust on a busy one --
    # which is the machine most runs land on.
    CPU_DECODE_WORKERS = 12
    GPU_DECODE_WORKERS = 16

    def __init__(self, genotype_path, sample_file=None, *, decode_backend='auto',
                 reader_workers=None, prefetch_chunks=2, metadata_path=None,
                 bgi_path=None, max_variants=None, decode_batch_size=None,
                 metadata_cache_dir=None, **kwargs):
        import sqlite3
        import struct
        import threading
        self.genotype_path = Path(genotype_path)
        # Per-thread decoder and staging buffer, invalidated by generation so a
        # sample re-selection cannot leave a thread decoding the old cohort.
        self._thread_state = threading.local()
        self._decode_generation = 0
        self.decode_backend = decode_backend
        if decode_backend not in {'auto', 'cpu', 'gpu'}:
            raise ValueError('decode_backend must be auto, cpu, or gpu')
        if reader_workers is None:
            # 'auto' resolves to the GPU decoder when one is available, so it
            # takes the GPU figure; the CPU count applies only when CPU decode
            # was asked for outright.
            reader_workers = (self.CPU_DECODE_WORKERS if decode_backend == 'cpu'
                              else self.GPU_DECODE_WORKERS)
        if reader_workers <= 0 or prefetch_chunks <= 0:
            raise ValueError('reader_workers and prefetch_chunks must be positive')
        self.reader_workers, self.prefetch_chunks = int(reader_workers), int(prefetch_chunks)
        self.decode_batch_size=None if decode_batch_size is None else int(decode_batch_size)
        if self.decode_batch_size is not None and self.decode_batch_size<=0:
            raise ValueError("decode_batch_size must be positive")
        self.backend_reason = None
        with self.genotype_path.open('rb') as f:
            fixed = f.read(20)
            if len(fixed) != 20:
                raise ValueError('truncated BGEN header')
            offset, header_length, n_variants, n_samples = struct.unpack_from('<IIII', fixed)
            if header_length < 20 or fixed[16:20] not in (b'bgen', b'\x00'*4):
                raise ValueError('invalid BGEN header')
            f.seek(header_length)
            flags = struct.unpack('<I', f.read(4))[0]
            self._compression = flags & 3
            # Two different refusals, named separately. They used to share one
            # message, so a Layout-1 file and a file with a reserved codec were
            # reported identically and neither told the reader which it had --
            # and they have different answers: Layout 1 must be converted, while
            # a reserved compression value means the file is not standard.
            layout = (flags >> 2) & 15
            if layout != 2:
                raise ValueError(
                    'direct BGEN source requires Layout 2, this file declares '
                    f'Layout {layout}; convert it with qctool or plink2 first')
            if self._compression == 3:
                raise ValueError(
                    'BGEN compression flag 3 is reserved by the specification; '
                    'supported values are 0 (none), 1 (zlib) and 2 (zstd)')
            ids = None
            if flags & (1 << 31):
                block_length, count = struct.unpack('<II', f.read(8))
                if count != n_samples:
                    raise ValueError('BGEN sample count mismatch')
                ids = []
                for _ in range(n_samples):
                    length = struct.unpack('<H', f.read(2))[0]
                    value = f.read(length)
                    if len(value) != length:
                        raise ValueError('truncated BGEN sample ID')
                    ids.append(value.decode('utf-8'))
                if f.tell() > offset + 4:
                    raise ValueError('BGEN sample block overlaps genotype data')
        if sample_file is not None:
            table = pd.read_csv(sample_file, sep=r'\s+', dtype=str)
            if len(table) and table.iloc[0].astype(str).str.fullmatch(r'0').all():
                table = table.iloc[1:].reset_index(drop=True)
            if not {'ID_1','ID_2'}.issubset(table):
                raise ValueError('sample file must contain ID_1 and ID_2')
            supplied = table.ID_2.to_numpy(dtype=object)
            if len(supplied) != n_samples or (ids is not None and not np.array_equal(supplied, ids)):
                raise ValueError('sample file does not match BGEN sample order')
            ids = supplied
            families = table.ID_1.to_numpy(dtype=object)
        elif ids is None:
            raise ValueError('BGEN without embedded IDs requires a sample file')
        else:
            families = np.asarray(ids, dtype=object)
        self.sample_ids = np.asarray(ids, dtype=object)
        if len(set(self.sample_ids)) != n_samples:
            raise ValueError('BGEN sample IDs must be unique')
        self.family_ids = np.asarray(families, dtype=object)
        self._n_bgen_samples = int(n_samples)
        self._sample_indices = np.arange(n_samples, dtype=np.int32)
        index = Path(bgi_path) if bgi_path is not None else Path(str(self.genotype_path)+'.bgi')
        if not index.is_file():
            raise FileNotFoundError(f'BGEN index required: {index}; create with bgenix -g FILE -index')
        # Read bounded rows, then sort inexpensive integer offsets in memory.
        # SQLite ORDER BY file offset may spill millions of string rows to disk;
        # pandas fetchall also duplicates all tuples at peak initialization RAM.
        query = ('SELECT chromosome,position,rsid,number_of_alleles,allele1,allele2,'
                 'file_start_position,size_in_bytes FROM Variant')
        if max_variants is not None:
            if int(max_variants) <= 0:
                raise ValueError('max_variants must be positive')
            # Bounded diagnostic runs need the physically first records, so sort
            # in SQLite here. Full scans avoid the external metadata sort.
            query += f' ORDER BY file_start_position LIMIT {int(max_variants)}'
        count = n_variants if max_variants is None else min(n_variants,int(max_variants))
        # A cache is only valid for the whole index; a truncated diagnostic run
        # asks a different question and re-reads.
        cached = None
        self.metadata_cache_path = None
        if metadata_cache_dir is not None and max_variants is None:
            cached, self.metadata_cache_path = _load_bgen_index_cache(
                metadata_cache_dir, index)
        if cached is not None:
            self.chromosomes = cached["chromosomes"]
            self.positions = cached["positions"]
            self.marker_ids = cached["marker_ids"]
            self.other_alleles = cached["other_alleles"]
            self.effect_alleles = cached["effect_alleles"]
            self._offsets = cached["offsets"]
            self._lengths = cached["lengths"]
            if self._offsets.shape[0] != count:
                raise ValueError('cached BGI variant count does not match BGEN')
            self.input_bytes = int(self._lengths.sum())
            self.backend_used = 'cpu' if decode_backend == 'cpu' else None
            return
        columns = [np.empty(count,dtype=dtype) for dtype in
                   (object,np.int64,object,np.int32,object,object,np.uint64,np.uint64)]
        copied=0
        with sqlite3.connect(f'file:{index}?mode=ro', uri=True) as connection:
            cursor=connection.execute(query)
            while True:
                rows=cursor.fetchmany(65536)
                if not rows:break
                stop=copied+len(rows)
                if stop>count:raise ValueError('BGI has more variants than BGEN header')
                for target,values in zip(columns,zip(*rows)):
                    target[copied:stop]=values
                copied=stop
        if copied!=count:
            raise ValueError('BGI variant count does not match BGEN')
        if np.any(columns[3] != 2):
            raise ValueError('direct BGEN supports biallelic variants; explicitly filter multiallelic input first')
        if count>1 and np.any(columns[6][1:]<columns[6][:-1]):
            order=np.argsort(columns[6],kind='stable')
            columns=[column[order] for column in columns]
        (self.chromosomes,self.positions,self.marker_ids,_,self.other_alleles,
         self.effect_alleles,self._offsets,self._lengths)=columns
        size = self.genotype_path.stat().st_size
        if np.any(self._lengths < 16) or np.any(self._offsets > size) or np.any(self._lengths > size-self._offsets):
            raise ValueError('BGI contains invalid record extents')
        if len(self._offsets)>1 and np.any(self._offsets[1:] < self._offsets[:-1]+self._lengths[:-1]):
            raise ValueError('BGI contains overlapping records')
        self.input_bytes = int(self._lengths.sum())
        # Only after every validation above has passed: a cache of a file the
        # loader would have rejected is worse than no cache.
        if metadata_cache_dir is not None and max_variants is None:
            try:
                self.metadata_cache_path = write_bgen_index_cache(
                    metadata_cache_dir, index,
                    chromosomes=self.chromosomes, positions=self.positions,
                    marker_ids=self.marker_ids,
                    other_alleles=self.other_alleles,
                    effect_alleles=self.effect_alleles,
                    offsets=self._offsets, lengths=self._lengths)
            except OSError:
                # A cache that cannot be written is a missed speedup, never a
                # failed scan.
                pass
        self.backend_used = 'cpu' if decode_backend == 'cpu' else None

    @property
    def shape(self):
        return len(self.sample_ids), len(self.marker_ids)

    @property
    def genotype(self):
        return self

    @property
    def variant_metadata(self):
        return {'chromosome':self.chromosomes, 'position':self.positions,
                'effect_allele':self.effect_alleles, 'other_allele':self.other_alleles}

    def select_samples(self, sample_ids):
        wanted = np.asarray(sample_ids, dtype=object)
        lookup = {str(x):i for i,x in enumerate(self.sample_ids)}
        if len(set(map(str,wanted))) != len(wanted):
            raise ValueError('duplicate requested sample IDs')
        absent = [x for x in wanted if str(x) not in lookup]
        if absent:
            raise ValueError(f'samples absent from BGEN: {absent[:5]}')
        order = np.asarray([lookup[str(x)] for x in wanted], dtype=np.int64)
        self._sample_indices = np.ascontiguousarray(self._sample_indices[order], dtype=np.int32)
        self.sample_ids, self.family_ids = self.sample_ids[order], self.family_ids[order]
        # Cached per-thread decoders hold a copy of `_sample_indices`, and the
        # cached staging buffers are sized by the sample count. Both are stale
        # the moment the selection changes, and a stale decoder would keep
        # emitting the previous cohort's columns without any error, so the
        # generation is bumped and every thread rebuilds on its next chunk.
        self._decode_generation += 1
        return self

    def _records(self, start, end):
        import os
        fd = os.open(self.genotype_path, os.O_RDONLY)
        try:
            first = start
            while first < end:
                stop = first+1
                extent = int(self._offsets[first]+self._lengths[first])
                while stop < end and int(self._offsets[stop]) == extent:
                    extent += int(self._lengths[stop]); stop += 1
                begin = int(self._offsets[first])
                size = extent-begin
                parts = []
                received = 0
                while received < size:
                    part = os.pread(fd, size-received, begin+received)
                    if not part:
                        raise OSError('short BGEN read')
                    parts.append(part); received += len(part)
                data = memoryview(parts[0] if len(parts)==1 else b''.join(parts))
                for i in range(first,stop):
                    offset = int(self._offsets[i])-begin
                    yield i, data[offset:offset+int(self._lengths[i])]
                first = stop
        finally:
            os.close(fd)

    def _record_payload(self, record, index):
        """Parse a record's header and return its probability block.

        Both decoders share this, so the BGI cross-check below happens exactly
        once and cannot drift between them. The block is returned still
        compressed, with the inflated length the record declares.
        """
        import struct
        p=0
        def take(fmt):
            nonlocal p
            n=struct.calcsize(fmt)
            if p+n>len(record):raise ValueError('truncated BGEN record')
            value=struct.unpack_from(fmt,record,p)[0];p+=n
            return value
        def text(width):
            nonlocal p
            n=take(width)
            if p+n>len(record):raise ValueError('truncated BGEN text')
            value=bytes(record[p:p+n]).decode('utf-8');p+=n
            return value
        variant_id,rsid,chrom=text('<H'),text('<H'),text('<H')
        position,alleles=take('<I'),take('<H')
        names=[text('<I') for _ in range(alleles)]
        if (chrom != str(self.chromosomes[index]) or position != self.positions[index]
            or rsid != self.marker_ids[index] or names != [self.other_alleles[index],self.effect_alleles[index]]):
            raise ValueError(f'BGI/BGEN metadata mismatch at variant {index}')
        compressed_length=take('<I')
        if self._compression == 0:
            return record[p:],compressed_length
        raw_length=take('<I')
        if compressed_length != len(record)-p+4:
            raise ValueError('BGEN compressed length mismatch')
        return record[p:],raw_length

    def _decode_record(self, record, index, dtype):
        """Reference decoder. The C path must agree with this one exactly."""
        import struct
        import zlib
        payload,raw_length=self._record_payload(record,index)
        if self._compression == 0:
            raw=bytes(payload)
        elif self._compression == 1:
            raw=zlib.decompress(payload)
        else:
            import zstandard
            raw=zstandard.ZstdDecompressor().decompress(payload,max_output_size=raw_length)
        n=self._n_bgen_samples
        if len(raw)!=raw_length or len(raw)<10+n:
            raise ValueError('BGEN inflated length mismatch')
        samples,k=struct.unpack_from('<IH',raw)
        if samples!=n or k!=2 or raw[6]!=2 or raw[7]!=2 or raw[8+n]!=0:
            raise ValueError(f'unsupported BGEN dimensions/ploidy/phase at variant {index}')
        bits=raw[9+n]
        if bits<1 or bits>32 or len(raw)!=10+n+(n*2*bits+7)//8:
            raise ValueError('invalid BGEN packed probability length')
        pm=np.frombuffer(raw,dtype=np.uint8,count=n,offset=8)[self._sample_indices]
        if np.any((pm&63)!=2):
            raise ValueError('non-diploid selected sample')
        missing=(pm&128)!=0
        data=np.frombuffer(raw,dtype=np.uint8,offset=10+n)
        bit_offsets=self._sample_indices.astype(np.uint64)*(2*bits)
        def unpack(offsets):
            byte=(offsets>>3).astype(np.int64);shift=offsets&7
            padded=np.pad(data,(0,5))
            word=np.zeros(len(offsets),dtype=np.uint64)
            for j in range((bits+14)//8):
                word |= padded[byte+j].astype(np.uint64) << (8*j)
            return (word>>shift)&np.uint64((1<<bits)-1)
        p0,p1=unpack(bit_offsets),unpack(bit_offsets+bits)
        denominator=np.uint64((1<<bits)-1)
        if np.any((p0+p1>denominator)&~missing):
            raise ValueError('invalid BGEN probabilities')
        result=((2*denominator-2*p0-p1).astype(np.float64)/float(denominator)).astype(dtype)
        result[missing]=np.nan
        return result

    def read_chunk(self,start,end,dtype=np.float32):
        if not 0<=start<=end<=self.shape[1]:
            raise ValueError('invalid BGEN chunk bounds')
        dtype=np.dtype(dtype)
        if dtype.kind!='f':
            raise ValueError('direct BGEN emits floating point dosages')
        # One decoder per *thread*, not per chunk and certainly not per variant:
        # it owns the inflate scratch, so it must not be shared concurrently,
        # and rebuilding it every chunk threw away that scratch along with the
        # staging buffer.
        decoder,staging=(self._thread_decode_state(end-start,dtype)
                         if dtype==np.float32 else (None,None))
        if decoder is None:
            result=np.empty((self.shape[0],end-start),dtype=dtype)
            for i,record in self._records(start,end):
                result[:,i-start]=self._decode_record(record,i,dtype)
            return result
        # Decode variant-major, then transpose once.
        #
        # A chunk is (samples, variants), so writing a variant straight into
        # its column puts one float every `end-start` floats: at 35,365 samples
        # that is 35,365 separate cache lines scattered over a 72 MB buffer per
        # variant, and it cost 3.76x -- 1006 microseconds against 268 for the
        # identical decode into a contiguous vector. The threaded sweep read as
        # a concurrency ceiling (throughput flat past eight workers, under
        # three cores busy out of forty-eight) because the readers were stalled
        # on memory, not waiting on the GIL, which by then held only 1.2% of
        # the work.
        #
        # One blocked transpose per chunk moves 2 x 72 MB sequentially, which
        # costs far less than the scattered writes it replaces.
        for i,record in self._records(start,end):
            payload,raw_length=self._record_payload(record,i)
            try:
                decoder.decode_into(payload,raw_length,staging[i-start])
            except BgenDecodeError as exc:
                raise ValueError(f'{exc} at variant {i}') from exc
        return _native_transpose(staging)

    def _cpu_decoder(self):
        """A C decoder for this file, or None to use the numpy reference."""
        if not cpu_decoder_available(self._compression):
            return None
        return CpuBgenDecoder(self._n_bgen_samples,self._sample_indices,self._compression)

    def _thread_decode_state(self,rows,dtype):
        """A per-thread decoder and staging buffer for one chunk.

        Both used to be rebuilt on every chunk. Measured against a
        pre-faulted-buffer arm, that allocation is worth **1.10-1.26x** at every
        worker count -- a steady tax rather than the concurrency ceiling it was
        once suspected of being, but a real one: 512 x 35,365 float32 is 72 MB
        of fresh pages per chunk, and first-touch faulting is kernel-side.
        Reuse cut minor faults from 21,670 to 40 at four workers.

        **Reusing this buffer is safe, where the ring-buffer one was not.** It
        never escapes `read_chunk`: the transpose copies out of it into a freshly
        allocated array that the caller owns. The bug this project already paid
        for came from handing a reused buffer *out*.

        Binding is per thread, never per task index -- a `CpuBgenDecoder` holds
        an inflate scratch, and two concurrent decodes through one corrupt it
        into `zlib inflate failed`, which names the codec and not the aliasing.
        `threading.local()` gives that guarantee; `index % workers` only looks
        like it does.
        """
        state=self._thread_state
        if getattr(state,'generation',None)!=self._decode_generation:
            state.generation=self._decode_generation
            state.decoder=self._cpu_decoder()
            state.staging=None
        decoder=state.decoder
        if decoder is None:
            return None,None
        import os
        if os.environ.get('TORCHGWAS_BGEN_REUSE_STAGING','1')=='0':
            # Escape hatch, and the A/B switch the benefit was measured with.
            return decoder,np.empty((rows,self.shape[0]),dtype=dtype)
        staging=state.staging
        if (staging is None or staging.shape[0]<rows
                or staging.shape[1]!=self.shape[0] or staging.dtype!=dtype):
            # Grow to the largest chunk this thread has been asked for and keep
            # it. Sizing down again would reintroduce the faulting this exists
            # to avoid, and chunk sizes do not vary within a scan.
            staging=np.empty((rows,self.shape[0]),dtype=dtype)
            state.staging=staging
        return decoder,staging[:rows]

    def iter_chunks(self,chunk_size,dtype=np.float32,prefetch_chunks=None,reader_workers=None,
                    variant_range=None):
        from .streaming import OrderedChunkLoader
        self.backend_used='cpu'
        return OrderedChunkLoader(self.shape[1],self.read_chunk,chunk_size,dtype,
            self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks,
            self.reader_workers if reader_workers is None else reader_workers,
            variant_range=variant_range)

    def resolve_decode_backend(self, device):
        """Choose availability only; runtime decode/allocation errors propagate."""
        import torch
        from .bgen_gpu import load_decoder_library
        self.backend_reason=None
        if self.decode_backend=='cpu':
            self.backend_used='cpu'
            self.backend_reason='CPU decoder explicitly selected'
            return 'cpu'
        unavailable=None
        if torch.device(device).type!='cuda':
            unavailable='CPU scan device selected'
        elif self._compression!=1:
            unavailable='GPU BGEN decoder supports zlib only; using CPU for this compression'
        else:
            try:
                load_decoder_library()
            except ImportError as exc:
                unavailable=str(exc)
            else:
                # Loading is not enough. The library carries code only for the
                # architectures it was built for, and an unsupported card fails
                # at *launch* -- by which point this function has already
                # promised the GPU backend and the error surfaces mid-scan.
                from .bgen_gpu import decoder_launchable

                if not decoder_launchable(device):
                    unavailable=('GPU BGEN decoder has no kernel for this '
                                 'device architecture; rebuild with '
                                 'bash build_direct_bgen.sh')
        if unavailable is not None:
            if self.decode_backend=='gpu':
                raise ValueError(unavailable)
            self.backend_used='cpu'
            self.backend_reason=unavailable
        else:
            self.backend_used='gpu'
        return self.backend_used

    def iter_device_chunks(self,chunk_size,device,reader_workers=None,prefetch_chunks=None,
                           variant_range=None):
        import torch
        from concurrent.futures import ThreadPoolExecutor
        from collections import deque
        from .bgen_gpu import GpuBgenDecoder
        if chunk_size<=0:
            raise ValueError('chunk_size must be positive')
        depth=self.prefetch_chunks if prefetch_chunks is None else int(prefetch_chunks)
        if depth<=0:raise ValueError('prefetch_chunks must be positive')
        from .streaming import _resolve_variant_range
        span_start,span_end=_resolve_variant_range(variant_range,int(self.shape[1]))
        backend=self.resolve_decode_backend(device)
        import os as _os
        # Decode lanes. A single lane leaves the GPU idle for roughly 45% of
        # the producer wall: per tile the host does a 193 MB pageable->pinned
        # memcpy, then an H2D, then blocks on `cudaStreamSynchronize` before
        # tile k+1 begins, so nothing overlaps. Measured ~64 ms of GPU decode
        # work against ~119 ms of producer wall per 20,000-variant tile.
        #
        # Each lane needs its OWN decoder, because a decoder owns a handle, a
        # CUDA stream and its device buffers -- two threads sharing one would
        # race on all three. Batches are lane-independent (`tg_bgen_read` takes
        # no handle), so any lane can decode any batch.
        #
        # Cost is device memory: each lane holds its own staging buffers, so
        # lanes trade memory for overlap. Two is the useful step; set
        # TORCHGWAS_BGEN_DECODE_LANES=1 to restore the old behaviour.
        lanes=max(1,int(_os.environ.get('TORCHGWAS_BGEN_DECODE_LANES','2')))
        decoders=([GpuBgenDecoder(self._n_bgen_samples,self._sample_indices,device)
                   for _ in range(lanes)] if backend=='gpu' else [])
        decoder=decoders[0] if decoders else None
        if decoder is None:
            self.backend_used='cpu'
            for start,end,array in self.iter_chunks(chunk_size,np.float32,depth,reader_workers,
                                                   variant_range=(span_start,span_end)):
                yield start,end,torch.as_tensor(array.T.copy(),device=device)
            return
        self.backend_used='gpu'
        import os
        workers=self.reader_workers if reader_workers is None else int(reader_workers)
        if workers<=0:
            decoder.close()
            raise ValueError("reader_workers must be positive")
        tile=max(chunk_size,self.decode_batch_size or chunk_size)
        tile=((tile+chunk_size-1)//chunk_size)*chunk_size
        self.effective_decode_batch_size=tile
        from contextlib import ExitStack
        resources=ExitStack()
        def close_decoder():
            # Lane 0's profile is reported as before; the others are kept under
            # `lane_profiles` so a multi-lane run can still be accounted for
            # without changing the shape callers already read.
            try:
                self.last_decode_profile=decoder.profile()
                if len(decoders)>1:
                    self.last_decode_profile['lanes']=len(decoders)
                    self.last_decode_profile['lane_profiles']=[
                        lane.profile() for lane in decoders]
            finally:
                for lane in decoders:
                    lane.close()
        resources.callback(close_decoder)
        try:
            fd=os.open(self.genotype_path,os.O_RDONLY)
            resources.callback(os.close,fd)
            # Independent contiguous preads may be required to saturate storage;
            # GPU decode remains one ordered lane and outputs stay bounded.
            reads=resources.enter_context(ThreadPoolExecutor(
                max_workers=min(workers,depth),thread_name_prefix='bgen-read'))
            gpu=resources.enter_context(ThreadPoolExecutor(
                max_workers=len(decoders),thread_name_prefix='bgen-decode'))
        except BaseException:
            resources.close()
            raise
        pending=deque()
        def prepare(start,end):
            import time
            started=time.perf_counter();cpu_started=time.thread_time()
            metadata = ''.join('\x00'.join(row)+'\x00' for row in zip(self.chromosomes[start:end],self.marker_ids[start:end],self.other_alleles[start:end],self.effect_alleles[start:end])).encode('utf-8')
            if decoder.profile_enabled:
                decoder.host_timings.append(('metadata_pack',time.perf_counter()-started,time.thread_time()-cpu_started))
            return decoder.read(fd,self._offsets[start:end],self._lengths[start:end],metadata,self.positions[start:end])
        # One decoder per decode thread, claimed on that thread's first task and
        # kept for the life of the pool. A decoder is not thread-safe -- handle,
        # stream and device buffers are all per-instance -- so binding by thread
        # is what makes more than one lane sound.
        lane_lock=threading.Lock()
        unclaimed=list(decoders)
        lane_local=threading.local()
        def lane():
            mine=getattr(lane_local,'decoder',None)
            if mine is None:
                with lane_lock:
                    mine=unclaimed.pop()
                lane_local.decoder=mine
            return mine
        def decode(future,count):
            return lane().decode(future.result(),count)
        try:
            bounds=iter((s,min(s+tile,span_end)) for s in range(span_start,span_end,tile))
            def submit():
                bound=next(bounds,None)
                if bound is None:return False
                start,end=bound
                read=reads.submit(prepare,start,end)
                pending.append((start,end,gpu.submit(decode,read,end-start)))
                return True
            for _ in range(depth):
                if not submit():break
            while pending:
                start,end,future=pending.popleft()
                output=future.result()
                for offset in range(0,end-start,chunk_size):
                    stop=min(offset+chunk_size,end-start)
                    view=output[offset:stop]
                    view.record_stream(torch.cuda.current_stream(device))
                    yield start+offset,start+stop,view
                submit()
        finally:
            resources.close()

    def close(self):
        pass
