from __future__ import annotations

import hashlib
import json
import os
import threading
import time
from collections import deque
from contextlib import contextmanager
from collections.abc import Iterator, Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import metadata_cache


DOSAGE_SCALE = 127.5
# Preserve sequential access to large PGENs on shared storage. Additional
# readers stride through distant regions and can collapse source throughput.
DEFAULT_PGEN_DECODE_WORKERS = 1
DEFAULT_PGEN_DECODE_BATCH_SIZE = 5000
_PGEN_CPU_COUNT = max(1, os.cpu_count() or 1)
DEFAULT_PGEN_COMPRESSION_WORKERS = min(
    40,
    _PGEN_CPU_COUNT if _PGEN_CPU_COUNT < 8 else _PGEN_CPU_COUNT - 4,
)
PGEN_MISSING_POLICY = "mask_missing_calls_keep_variant"
_HARDCALL_TO_CODE = np.asarray([0, 128, 255], dtype=np.uint8)


def _fill_masked_codes(codes, gaps, raw):
    """Give each masked call its variant's observed mean, as an eight-bit code.

    The code space is full -- 0..255 already spans dosage 0..2 at scale 127.5 --
    so there is no sentinel to write and no way to defer the decision. Writing
    the observed mean is what masking reduces to once the value has to be a
    number: the centred contribution of such a call is zero either way.

    The mean is taken in dosage space when the caller still has the dosages
    (`raw`), because the hard-call code table 0/128/255 is not exactly linear
    in dosage -- 128 is the code for 1.0 but 127.5 is what the scale implies --
    so averaging codes and averaging dosages do not agree to the last bit.
    """
    observed = ~gaps
    counts = observed.sum(axis=1)
    source = codes if raw is None else raw
    sums = np.where(observed, source, 0).sum(axis=1, dtype=np.float64)
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
    if raw is not None:
        means = means * DOSAGE_SCALE
    return np.where(gaps, np.rint(means)[:, None].astype(np.uint8), codes)


def _open_pgen(*args, **kwargs):
    try:
        import pgenlib
    except ImportError as exc:
        # Reaching here means the native backend was declined or forced off --
        # the message must say *why* pgenlib is needed for this particular file,
        # because for a hard-call PGEN it is not needed at all any more and
        # "PGEN input requires pgenlib" would send the reader down a dead end.
        raise ImportError(
            "this PGEN needs the optional dependency pgenlib: dosage, phased "
            "and multiallelic records are outside the built-in decoder's scope, "
            "and TORCHGWAS_PGEN_BACKEND=pgenlib forces it as well. A hard-call "
            "PGEN needs no extra. Install with `pip install -e '.[pgen]'`"
        ) from exc
    return pgenlib.PgenReader(*args, **kwargs)


def resolve_pgen_mode_and_backend(path, requested_mode: str) -> tuple[str, str, str]:
    """Decide read mode and reader backend for a whole PGEN, and say why.

    Returns `(mode, backend, reason)` with `backend` in {"native", "pgenlib"}.

    Both decisions are made once per file rather than per variant, because a
    per-variant policy is not usable here: the real dosage PGEN has 8,930,997
    of 8,931,083 variants carrying a dosage track, so a reader refusing record
    by record would decode 86 of them and reject the rest.

    **`auto` now means what it says.** It used to resolve to `dosage`
    unconditionally, which made every PGEN -- including one holding nothing but
    hard calls -- require pgenlib, and fail outright when pgenlib was absent.
    It now asks the file: a PGEN whose records carry no dosage, phase or
    multiallelic track *is* a hard-call file, its dosages would be exactly
    0/1/2, and reading it as hard calls gives identical values with no optional
    C extension. An explicit `mode=` is still honoured exactly as given.

    `TORCHGWAS_PGEN_BACKEND` forces the backend; forcing `native` on a file
    outside its scope raises here rather than failing partway through a scan.
    """
    from . import pgen_native
    from .pgen_reader import file_scope

    forced = os.environ.get("TORCHGWAS_PGEN_BACKEND", "auto").strip().lower()
    if forced not in {"auto", "native", "pgenlib"}:
        raise ValueError(
            "TORCHGWAS_PGEN_BACKEND must be auto, native or pgenlib, "
            f"got {forced!r}")
    if requested_mode not in {"auto", "dosage", "hardcall"}:
        raise ValueError(
            f"PGEN mode must be auto, dosage or hardcall, got {requested_mode!r}")

    # Two independent questions, kept apart deliberately. What the *file* holds
    # is read from its index and does not depend on whether our shared library
    # happens to be built; conflating the two made `TORCHGWAS_PGEN_BACKEND=
    # pgenlib` change what `auto` believed about the data.
    try:
        scope = file_scope(path)
        scope_error = None
    except Exception as exc:  # noqa: BLE001 - any parse failure disqualifies it
        scope = None
        scope_error = f"the PGEN index did not parse: {type(exc).__name__}: {exc}"
    hardcall_only = scope is not None and scope.supported

    mode = ("hardcall" if hardcall_only else "dosage") if requested_mode == "auto" \
        else requested_mode

    if forced == "pgenlib":
        return mode, "pgenlib", "TORCHGWAS_PGEN_BACKEND=pgenlib"
    # A parse failure is reported ahead of the mode, because it is the root
    # cause of both: an index that did not parse is *why* `auto` fell back to
    # dosage, and "dosage records are read by pgenlib" would send the reader
    # looking in the wrong place.
    if scope_error is not None:
        blocked = scope_error
    elif mode != "hardcall":
        blocked = ("dosage records are read by pgenlib; the native decoder "
                   "reads hard calls only")
    elif not hardcall_only:
        blocked = f"the file is outside the native decoder's scope: {scope.reason}"
    elif not pgen_native.available():
        blocked = ("the native PGEN decoder is not built "
                   "(run build_pgen_decode.sh, or set TORCHGWAS_PGEN_LIBRARY)")
    else:
        return mode, "native", "hard-call records within the native decoder's scope"
    if forced == "native":
        raise ValueError(f"TORCHGWAS_PGEN_BACKEND=native refused: {blocked}")
    return mode, "pgenlib", blocked


def pgenlib_version() -> str:
    try:
        import pgenlib
    except ImportError:
        return "unavailable"
    return str(getattr(pgenlib, "__version__", "unknown"))


def resolve_pgen_triplet(
    genotype_path: str | Path,
    *,
    pvar: str | Path | None = None,
    psam: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    """Resolve a PLINK 2 prefix or .pgen path without truncating dotted prefixes."""

    path = Path(genotype_path)
    prefix = Path(str(path)[: -len(".pgen")]) if path.suffix.lower() == ".pgen" else path
    pgen_path = path if path.suffix.lower() == ".pgen" else Path(f"{prefix}.pgen")
    pvar_path = Path(pvar) if pvar is not None else Path(f"{prefix}.pvar")
    psam_path = Path(psam) if psam is not None else Path(f"{prefix}.psam")
    for companion in (pgen_path, pvar_path, psam_path):
        if not companion.is_file():
            raise FileNotFoundError(companion)
    return pgen_path, pvar_path, psam_path


def _split_fields(line: str) -> list[str]:
    return line.rstrip("\r\n").split()


def _read_psam(path: Path) -> tuple[np.ndarray, np.ndarray]:
    header: list[str] | None = None
    iid_col: int | None = None
    fid_col: int | None = None
    iids: list[str] = []
    fids: list[str] = []
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("##"):
                continue
            fields = _split_fields(line)
            if header is None:
                header = fields
                header[0] = header[0].lstrip("#")
                if "IID" not in header:
                    raise ValueError(f"PSAM file must contain an IID column: {path}")
                iid_col = header.index("IID")
                fid_col = header.index("FID") if "FID" in header else None
                continue
            if len(fields) != len(header):
                raise ValueError(
                    f"PSAM row {row_number} has {len(fields)} fields; expected {len(header)}"
                )
            iid = fields[iid_col]
            iids.append(iid)
            fids.append(fields[fid_col] if fid_col is not None else "0")
    if header is None:
        raise ValueError(f"PSAM file has no header: {path}")
    if not iids:
        raise ValueError(f"PSAM file contains no samples: {path}")
    if len(set(iids)) != len(iids):
        raise ValueError(f"PSAM IID values must be unique: {path}")
    return np.asarray(fids, dtype=object), np.asarray(iids, dtype=object)


def _read_pvar(path: Path) -> dict[str, np.ndarray]:
    required = ("CHROM", "POS", "ID", "REF", "ALT")
    skiprows = 0
    header_line = ""
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                skiprows += 1
                continue
            header_line = line
            break
        else:
            raise ValueError(f"PVAR file has no header: {path}")

    # Both tab and whitespace-separated PVAR are accepted. Keep chromosome,
    # allele, and ID text lexical, but parse POS directly as int64 to avoid an
    # additional full column of temporary Python strings. Ignore annotations
    # such as INFO which are not needed by this biallelic scan source.
    separator = "\t" if "\t" in header_line else r"\s+"
    columns = [name.lstrip("#") for name in header_line.strip().split()]
    missing = [name for name in required if name not in columns]
    if missing:
        raise ValueError(f"PVAR file is missing required columns {missing}: {path}")
    try:
        frame = pd.read_csv(
            path,
            sep=separator,
            skiprows=skiprows,
            header=0,
            names=columns,
            usecols=list(required),
            dtype={name: np.int64 if name == "POS" else str for name in required},
            keep_default_na=False,
            na_filter=False,
        )
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"PVAR file contains an invalid integer POS value or malformed row: {path}") from exc
    if frame.empty:
        raise ValueError(f"PVAR file contains no variants: {path}")
    positions = frame["POS"].to_numpy(dtype=np.int64, copy=False)
    chromosomes = frame["CHROM"].to_numpy(dtype=object, copy=False)
    refs = frame["REF"].to_numpy(dtype=object, copy=False)
    alts = frame["ALT"].to_numpy(dtype=object, copy=False)
    marker_ids = frame["ID"].to_numpy(dtype=object, copy=True)
    absent = frame["ID"].isin(("", ".", "NA", "nan", "None")).to_numpy(
        dtype=bool, copy=False
    )
    if bool(absent.any()):
        marker_ids[absent] = [
            f"{chrom}:{position}:{ref}:{alt}"
            for chrom, position, ref, alt in zip(
                chromosomes[absent], positions[absent], refs[absent], alts[absent]
            )
        ]
    multiallelic = np.flatnonzero(
        frame["ALT"].str.contains(",", regex=False).to_numpy(dtype=bool, copy=False)
    )
    if multiallelic.size:
        preview = ", ".join(str(index) for index in multiallelic[:5])
        raise ValueError(
            f"PGEN conversion currently supports biallelic variants only; "
            f"PVAR has {multiallelic.size} multiallelic rows (0-based indexes: {preview})"
        )
    return {
        "chromosome": np.asarray(chromosomes, dtype=object),
        "position": np.asarray(positions, dtype=np.int64),
        "marker_id": np.asarray(marker_ids, dtype=object),
        "other_allele": np.asarray(refs, dtype=object),
        "effect_allele": np.asarray(alts, dtype=object),
    }


_PVAR_CACHE_ARRAYS = ("chromosome", "position", "marker_id", "other_allele",
                      "effect_allele")


def _pvar_cache_path(cache_dir: str | Path, pvar_path: Path,
                     *, create: bool = False) -> Path:
    return metadata_cache.cache_path(cache_dir, Path(pvar_path), "pvar",
                                     create=create)


def write_pgen_pvar_cache(cache_dir: str | Path, pvar_path: str | Path,
                          metadata: dict[str, np.ndarray]) -> Path:
    """Store parsed PVAR columns so the next open does not re-parse them."""
    return metadata_cache.store_arrays(
        cache_dir, pvar_path, "pvar",
        {name: (np.asarray(metadata[name], dtype=np.int64) if name == "position"
                else np.asarray(metadata[name], dtype=str))
         for name in _PVAR_CACHE_ARRAYS})


def _load_pgen_pvar_cache(cache_dir: str | Path, pvar_path: Path):
    """Parsed PVAR columns, or None when there is no usable cache."""
    arrays = metadata_cache.load_arrays(cache_dir, pvar_path, "pvar",
                                        _PVAR_CACHE_ARRAYS)
    return arrays, _pvar_cache_path(cache_dir, Path(pvar_path))


def _read_pvar_cached(pvar_path: Path, cache_dir: str | Path | None):
    """`_read_pvar`, with the parse skipped when a valid cache exists.

    Worth caching because the parse is paid on every open regardless of how
    much of the file is then scanned: on the 8.93M-variant benchmark cohort it
    is 5.7-7.3 s, which on the fastest transport is 60% of the whole run.
    """
    return metadata_cache.cached_arrays(
        cache_dir, pvar_path, "pvar", _PVAR_CACHE_ARRAYS,
        lambda: _read_pvar(pvar_path))


def sample_selection_identity(sample_ids: Sequence[object] | np.ndarray | None) -> dict | None:
    if sample_ids is None:
        return None
    values = [str(value) for value in sample_ids]
    digest = hashlib.sha256()
    for value in values:
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
    return {"count": len(values), "ordered_iid_sha256": digest.hexdigest()}


class PgenDosageSource:
    """Bounded-memory PLINK 2 PGEN source for native zstd conversion.

    The output dosage is the count of ALT1 (PVAR ``ALT``), while PVAR ``REF``
    is the other allele. In ``hardcall`` mode calls 0/1/2 map exactly to
    uint8 codes 0/128/255. A variant is excluded when any selected sample is
    missing; invalid calls are also excluded and reported.

    The default ``auto`` mode uses contiguous dosage reads for both dosage and
    hardcall PGENs. ``auto`` and ``dosage`` require an installed pgenlib with
    ``read_dosages_range``; explicit ``hardcall`` mode is the compatibility
    path for older readers and must only be used on known hardcall-only input.
    """

    dosage_scale = DOSAGE_SCALE
    missing_policy = PGEN_MISSING_POLICY

    def __init__(
        self,
        genotype_path: str | Path,
        *,
        pvar: str | Path | None = None,
        psam: str | Path | None = None,
        selected_sample_ids: Sequence[object] | np.ndarray | None = None,
        mode: str = "auto",
        reader_workers: int = DEFAULT_PGEN_DECODE_WORKERS,
        decode_batch_size: int = DEFAULT_PGEN_DECODE_BATCH_SIZE,
        metadata_cache_dir: str | Path | None = None,
    ) -> None:
        if mode not in {"auto", "hardcall", "dosage"}:
            raise ValueError("PGEN mode must be 'auto', 'hardcall', or 'dosage'")
        if reader_workers <= 0 or decode_batch_size <= 0:
            raise ValueError("reader_workers and decode_batch_size must be positive")
        self.genotype_path, self.pvar_path, self.psam_path = resolve_pgen_triplet(
            genotype_path, pvar=pvar, psam=psam
        )
        self.requested_mode = mode
        self.reader_workers = int(reader_workers)
        self.decode_batch_size = int(decode_batch_size)
        metadata_started = time.perf_counter()
        all_family_ids, all_sample_ids = _read_psam(self.psam_path)
        metadata = _read_pvar_cached(self.pvar_path, metadata_cache_dir)
        self.metadata_parse_seconds = time.perf_counter() - metadata_started
        self._raw_n_samples = int(all_sample_ids.size)
        self._raw_n_variants = int(metadata["marker_id"].size)
        self._all_metadata = metadata

        if selected_sample_ids is None:
            requested_indices = np.arange(self._raw_n_samples, dtype=np.uint32)
        else:
            requested = np.asarray([str(value) for value in selected_sample_ids], dtype=object)
            if requested.size == 0:
                raise ValueError("PGEN sample selection cannot be empty")
            if len(set(requested.tolist())) != int(requested.size):
                raise ValueError("PGEN sample selection contains duplicate IID values")
            index_by_iid = {str(iid): index for index, iid in enumerate(all_sample_ids)}
            absent = [str(iid) for iid in requested if str(iid) not in index_by_iid]
            if absent:
                preview = ", ".join(absent[:5])
                raise ValueError(
                    f"PGEN sample selection contains {len(absent)} IID values absent from PSAM "
                    f"({preview})"
                )
            requested_indices = np.asarray(
                [index_by_iid[str(iid)] for iid in requested], dtype=np.uint32
            )
        sorted_indices = np.sort(requested_indices)
        self._reader_sample_subset = (
            None
            if np.array_equal(sorted_indices, np.arange(self._raw_n_samples, dtype=np.uint32))
            else np.ascontiguousarray(sorted_indices, dtype=np.uint32)
        )
        self._sample_reorder = np.searchsorted(sorted_indices, requested_indices)
        self._reorder_is_identity = np.array_equal(
            self._sample_reorder, np.arange(requested_indices.size)
        )
        self.family_ids = all_family_ids[requested_indices.astype(np.int64)]
        self.sample_ids = all_sample_ids[requested_indices.astype(np.int64)]
        self._n_samples = int(self.sample_ids.size)

        # Resolve mode and backend before the first reader exists: `_new_reader`
        # needs both, and they are properties of the file rather than of any
        # one reader.
        (self.mode, self.pgen_backend,
         self.pgen_backend_reason) = resolve_pgen_mode_and_backend(
            self.genotype_path, mode)
        self._probe_reader = self._new_reader()
        observed_samples = int(self._probe_reader.get_raw_sample_ct())
        observed_variants = int(self._probe_reader.get_variant_ct())
        if observed_samples != self._raw_n_samples or observed_variants != self._raw_n_variants:
            self.close()
            raise ValueError(
                "PGEN dimensions do not match PSAM/PVAR metadata: "
                f"PGEN={observed_samples}x{observed_variants}, "
                f"metadata={self._raw_n_samples}x{self._raw_n_variants}"
            )
        has_dosage_range = callable(getattr(self._probe_reader, "read_dosages_range", None))
        self._has_packed_range = callable(getattr(self._probe_reader, "read_packed_range_into", None))
        if self.mode == "dosage" and not has_dosage_range:
            self.close()
            raise RuntimeError(
                "PGEN dosage conversion requires pgenlib.read_dosages_range(); "
                "the installed pgenlib only exposes scalar dosage reads. "
                "Upgrade pgenlib, or use mode='hardcall' only when the PGEN is known "
                "to contain hardcalls without dosage precision."
            )
        self.effect_allele_convention = (
            "PVAR_ALT1_DOSAGE" if self.mode == "dosage" else "PVAR_ALT1_HARDCALL"
        )
        if self.reader_workers > 1:
            # Parallel decoding uses one independent reader per worker thread;
            # do not retain an otherwise-idle probe reader for the full conversion.
            self.close()

        # `masked_missing` counts variants that carried a missing call. They
        # are kept and analysed, so it is a record of how much masking the
        # file needed, not a count of anything discarded.
        self.exclusion_counts = {"masked_missing": 0, "invalid_call": 0}
        self.marker_ids = np.empty(0, dtype=object)
        self.chromosomes = np.empty(0, dtype=object)
        self.positions = np.empty(0, dtype=np.int64)
        self.effect_alleles = np.empty(0, dtype=object)
        self.other_alleles = np.empty(0, dtype=object)
        self._kept_indices = np.empty(0, dtype=np.int64)
        self._converted = False
        self.decode_worker_seconds = 0.0
        self.read_worker_seconds = 0.0
        self.transform_worker_seconds = 0.0
        self.decoded_logical_bytes = 0

    def _new_reader(self):
        # One reader per worker thread; the backend is a per-file property and
        # is resolved once, before the first reader is built.
        backend = self.pgen_backend
        if backend == "native":
            from .pgen_native_reader import NativePgenReader

            return NativePgenReader(
                self.genotype_path,
                raw_sample_ct=self._raw_n_samples,
                variant_ct=self._raw_n_variants,
                sample_subset=self._reader_sample_subset,
            )
        return _open_pgen(
            os.fsencode(self.genotype_path),
            raw_sample_ct=self._raw_n_samples,
            variant_ct=self._raw_n_variants,
            sample_subset=self._reader_sample_subset,
        )

    def close(self) -> None:
        reader = getattr(self, "_probe_reader", None)
        if reader is not None:
            reader.close()
            self._probe_reader = None

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
            raise RuntimeError("PGEN metadata is final only after the conversion stream is exhausted")
        return {
            "chromosome": self.chromosomes,
            "position": self.positions,
            "effect_allele": self.effect_alleles,
            "other_allele": self.other_alleles,
        }

    def _finalize_metadata(self, kept: list[np.ndarray]) -> None:
        self._kept_indices = (
            np.concatenate(kept).astype(np.int64, copy=False)
            if kept
            else np.empty(0, dtype=np.int64)
        )
        metadata = self._all_metadata
        self.marker_ids = metadata["marker_id"][self._kept_indices]
        self.chromosomes = metadata["chromosome"][self._kept_indices]
        self.positions = metadata["position"][self._kept_indices]
        self.effect_alleles = metadata["effect_allele"][self._kept_indices]
        self.other_alleles = metadata["other_allele"][self._kept_indices]
        self._converted = True

    def _decode_range(self, reader, start: int, end: int):
        row_count = end - start
        read_started = time.perf_counter()
        if self.mode == "hardcall":
            values = np.empty((row_count, self._n_samples), dtype=np.int8)
            reader.read_range(start, end, values, allele_idx=1, sample_maj=False)
            read_seconds = time.perf_counter() - read_started
            transform_started = time.perf_counter()
            # The range is taken over observed calls only: a missing call must
            # not be able to make a well-formed row look malformed, and a
            # genuinely malformed row must still be caught.
            observed = values != -9
            row_min = np.where(observed, values, 2).min(axis=1)
            row_max = np.where(observed, values, 0).max(axis=1)
            missing = ~observed.all(axis=1)
            invalid = (~observed.any(axis=1)) | (row_min < 0) | (row_max > 2)
            valid = ~invalid
            selected = values if bool(valid.all()) else values[valid]
            if not self._reorder_is_identity:
                selected = selected[:, self._sample_reorder]
            gaps = selected == -9
            codes = np.asarray(_HARDCALL_TO_CODE[np.where(gaps, 0, selected)],
                               dtype=np.uint8, order="C")
            if gaps.any():
                codes = _fill_masked_codes(codes, gaps, selected)
        else:
            values = np.empty((row_count, self._n_samples), dtype=np.float32)
            try:
                reader.read_dosages_range(
                    start, end, values, allele_idx=1, sample_maj=False
                )
            except TypeError as exc:
                raise RuntimeError(
                    "installed pgenlib.read_dosages_range() has an incompatible API; "
                    "TorchGWAS requires contiguous variant-major dosage range reads"
                ) from exc
            read_seconds = time.perf_counter() - read_started
            transform_started = time.perf_counter()
            # As above: judge the row on its observed calls. A missing call is
            # masked, not a reason to discard everything else in the row.
            observed = values != -9.0
            row_min = np.where(observed, values, 2.0).min(axis=1)
            row_max = np.where(observed, values, 0.0).max(axis=1)
            missing = ~observed.all(axis=1)
            invalid = (
                (~observed.any(axis=1))
                | ~np.isfinite(row_min)
                | ~np.isfinite(row_max)
                | (row_min < 0.0)
                | (row_max > 2.0)
            )
            valid = ~invalid
            selected = values if bool(valid.all()) else values[valid]
            if not self._reorder_is_identity:
                selected = selected[:, self._sample_reorder]
            gaps = selected == -9.0
            selected = np.where(gaps, 0.0, selected).astype(np.float32)
            np.clip(selected, 0.0, 2.0, out=selected)
            selected *= np.float32(DOSAGE_SCALE)
            np.rint(selected, out=selected)
            codes = selected.astype(np.uint8)
            if gaps.any():
                codes = _fill_masked_codes(codes, gaps, None)
        transform_seconds = time.perf_counter() - transform_started
        indices = np.flatnonzero(valid).astype(np.int64) + start
        return (
            codes,
            indices,
            int(missing.sum()),
            int(invalid.sum()),
            read_seconds,
            transform_seconds,
            int(values.size),
        )

    def iter_chunks(
        self,
        chunk_size: int,
        dtype: np.dtype = np.uint8,
        prefetch_chunks: int | None = None,
        reader_workers: int | None = None,
    ) -> Iterator[tuple[int, int, np.ndarray]]:
        if self._converted:
            raise RuntimeError("a PgenDosageSource conversion stream can only be consumed once")
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if np.dtype(dtype) != np.dtype(np.uint8):
            raise ValueError("PgenDosageSource emits quantized uint8 dosage codes")
        workers = self.reader_workers if reader_workers is None else int(reader_workers)
        if workers <= 0:
            raise ValueError("reader_workers must be positive")
        # One queued range per reader is sufficient to keep every independent
        # decoder busy; a deeper queue only multiplies large cohort buffers.
        max_pending = max(
            1,
            min(workers, workers if prefetch_chunks is None else int(prefetch_chunks)),
        )
        ranges = iter(
            (start, min(self._raw_n_variants, start + self.decode_batch_size))
            for start in range(0, self._raw_n_variants, self.decode_batch_size)
        )
        thread_local = threading.local()
        readers: list = []
        readers_lock = threading.Lock()

        def decode(start: int, end: int):
            reader = getattr(thread_local, "reader", None)
            if reader is None:
                reader = self._new_reader()
                thread_local.reader = reader
                with readers_lock:
                    readers.append(reader)
            return self._decode_range(reader, start, end)

        def decoded_results():
            if workers == 1:
                reader = self._probe_reader
                if reader is None:
                    reader = self._new_reader()
                    readers.append(reader)
                for start, end in ranges:
                    yield self._decode_range(reader, start, end)
                return
            pending: deque[Future] = deque()
            with ThreadPoolExecutor(
                max_workers=workers, thread_name_prefix="torchgwas-pgen-decode"
            ) as pool:
                for _ in range(max_pending):
                    try:
                        start, end = next(ranges)
                    except StopIteration:
                        break
                    pending.append(pool.submit(decode, start, end))
                while pending:
                    yield pending.popleft().result()
                    try:
                        start, end = next(ranges)
                    except StopIteration:
                        continue
                    pending.append(pool.submit(decode, start, end))

        code_parts: deque[np.ndarray] = deque()
        index_parts: deque[np.ndarray] = deque()
        buffered = 0
        emitted = 0
        kept: list[np.ndarray] = []

        def take(count: int) -> tuple[np.ndarray, np.ndarray]:
            nonlocal buffered
            code_rows: list[np.ndarray] = []
            index_rows: list[np.ndarray] = []
            remaining = count
            while remaining:
                codes = code_parts[0]
                indices = index_parts[0]
                use = min(remaining, codes.shape[0])
                code_rows.append(codes[:use])
                index_rows.append(indices[:use])
                if use == codes.shape[0]:
                    code_parts.popleft()
                    index_parts.popleft()
                else:
                    code_parts[0] = codes[use:]
                    index_parts[0] = indices[use:]
                remaining -= use
                buffered -= use
            output = code_rows[0] if len(code_rows) == 1 else np.concatenate(code_rows, axis=0)
            output_indices = (
                index_rows[0] if len(index_rows) == 1 else np.concatenate(index_rows)
            )
            return output, output_indices

        try:
            for (
                codes,
                indices,
                missing_count,
                invalid_count,
                read_seconds,
                transform_seconds,
                logical_bytes,
            ) in decoded_results():
                self.exclusion_counts["masked_missing"] += missing_count
                self.exclusion_counts["invalid_call"] += invalid_count
                self.read_worker_seconds += read_seconds
                self.transform_worker_seconds += transform_seconds
                self.decode_worker_seconds += read_seconds + transform_seconds
                self.decoded_logical_bytes += logical_bytes
                if codes.shape[0]:
                    code_parts.append(codes)
                    index_parts.append(indices)
                    buffered += codes.shape[0]
                while buffered >= chunk_size:
                    output, output_indices = take(chunk_size)
                    kept.append(output_indices)
                    yield emitted, emitted + chunk_size, output.T
                    emitted += chunk_size
            if buffered:
                output, output_indices = take(buffered)
                kept.append(output_indices)
                yield emitted, emitted + output.shape[0], output.T
            self._finalize_metadata(kept)
        finally:
            for reader in readers:
                reader.close()


class PgenGenotype(PgenDosageSource):
    """Direct, repeatable biallelic ALT1 PGEN scans without a converted store.

    Dosages retain pgenlib's float32 precision; hardcalls are exactly 0/1/2.
    Missing calls become NaN for the scan engine's existing missing-data policy.
    Rows are never removed, so variant metadata and output indexes stay aligned.
    Each iterator owns its readers and closes them even after early termination.
    """

    dosage_scale = 1.0
    dtype = np.dtype(np.float32)
    missing_policy = "nan_per_missing_call_preserve_variant_axis"
    supports_fused_qc = True
    allows_direct_native_fill = True
    validate_native_range = True
    decode_tile_matches_chunk = True

    def __init__(self, *args, prefetch_chunks: int | None = None, **kwargs):
        if prefetch_chunks is not None and prefetch_chunks <= 0:
            raise ValueError("prefetch_chunks must be positive")
        super().__init__(*args, **kwargs)
        self.prefetch_chunks = max(2, self.reader_workers) if prefetch_chunks is None else int(prefetch_chunks)
        self.decode_workers = self.reader_workers
        self.native_dtype = np.dtype(np.int8 if self.mode == "hardcall" else np.float32)
        # OPT-IN uint8 dosage transport. Dosage ships float32 today, which is
        # 141,460 B per variant at N=35,365 -- 1.26 TB over PCIe for the full
        # genome, and the measured scan (48.95 s) is exactly that at 25.8 GB/s,
        # i.e. the format is bus-bound rather than decode-bound. uint8 codes
        # carry the same information: the zstd store already stores dosage this
        # way, and the cross-format check measured max |delta| 0.003922 against
        # pgen float32 -- exactly the 1/255 step -- with median 0.000000.
        #
        # Scale 127.0 rather than zstd's 127.5 so codes span 0..254 and **255
        # is free as the missing sentinel**. At 127.5 the code space is full
        # (0..255 already spans dosage 0..2) and masked calls have to be filled
        # with the variant mean, which loses the per-variant observed count and
        # so perturbs the residual df. One spare code buys exact missingness for
        # a step of 2/254 against 2/255.
        #
        # Behind a flag because the win is NOT yet measured: the saving is ~35 s
        # of PCIe, but quantising on the host costs a pass over M*N elements
        # (3.16e11 at full scale), which is the same order. The right fix is to
        # quantise inside `read_dosages_range` in the native reader so there is
        # no second pass at all; this flag exists to measure whether the
        # transport saving is real before paying for that C change.
        self.dosage_uint8_transport = (
            self.mode != "hardcall"
            and os.environ.get("TORCHGWAS_PGEN_DOSAGE_UINT8") == "1")
        if self.dosage_uint8_transport:
            self.native_dtype = np.dtype(np.uint8)
            self.native_scale = 127.0
            self.native_missing_value = 255
        self.native_transfer_dtype = self.native_dtype
        self.native_row_width = self._n_samples
        self.native_encoding = "dosage"
        # Packed two-bit transport, used wherever it is possible rather than
        # only when asked for.
        #
        # This was opt-in for a good reason that has since expired. Packed
        # transport originally needed a separately built pgenlib bridge wheel
        # installed into an isolated environment (see
        # docs/direct_pgen_packed_transport.md), so defaulting to it would have
        # failed for everyone who had not built it. `NativePgenReader` now
        # implements `read_packed_range_into` itself, so the bridge is not a
        # dependency any more and the gate outlived its reason.
        #
        # What the gate cost while it stayed shut, measured on the benchmark
        # cohort: hard calls went across PCIe as int8, one byte per sample,
        # **22,250 bytes per variant = 198.7 GB** for the file, where the same
        # calls packed two bits wide are **5,568 bytes per variant = 49.7 GB**
        # -- and 49.7 GB is exactly what BED sends for the same genotypes. So
        # every PGEN measurement taken with the gate shut moved four times the
        # payload of the equivalent BED run, and the resulting "PGEN reads at
        # 0.95 GB/s, 3.6x off the disk rate" was a statement about the bus, not
        # about storage.
        #
        # Eligibility is a property of the file, the request and the
        # statistics backend -- not of a flag of its own. Hard calls only,
        # every sample in file order, a reader that can fill packed rows, and
        # the native fused statistics kernels, which are the only consumer that
        # can read a packed row: `native_scan.dosage_cuda_iterator` refuses
        # `pgen_2bit` outright without them. That coupling is why this cannot
        # simply default to on -- doing so turns a working torch-backend run
        # into a hard error at scan time, which is exactly what the first
        # attempt at this did.
        #
        # `TORCHGWAS_PGEN_PACKED=0` forces the old int8 transport, which is
        # what an A/B of the two needs; `=1` turns an ineligible source into an
        # error rather than a silent fallback.
        requested = os.environ.get("TORCHGWAS_PGEN_PACKED")
        if self.mode != "hardcall":
            refusal = ("packed PGEN transport requires explicit hardcall mode; "
                       "it carries two-bit calls, not dosages")
        elif self._reader_sample_subset is not None or not self._reorder_is_identity:
            refusal = "packed PGEN transport requires all samples in file order"
        elif not self._has_packed_range:
            refusal = ("packed PGEN transport requires a reader providing "
                       "read_packed_range_into()")
        else:
            refusal = None
        if requested == "1":
            # Explicit opt-in keeps its original contract, including raising
            # for dosage mode and for sample subsets. A caller asking for the
            # packed transport by name may be driving `native_reader_session`
            # directly and never running a scan at all, so this must not be
            # made conditional on a statistics backend it does not use.
            if refusal is not None:
                self.close()
                raise ValueError(f"TORCHGWAS_PGEN_PACKED=1 but {refusal}")
            enable_packed = True
        elif requested == "0":
            enable_packed = False
        else:
            # Automatic, and deliberately conservative: only when the native
            # fused statistics kernels are already in use, because they are the
            # only consumer that can read a packed row --
            # `native_scan.dosage_cuda_iterator` refuses `pgen_2bit` outright
            # without them. Enabling packed unconditionally would turn a
            # working torch-backend scan into a hard error at scan time, which
            # is exactly what the first attempt at this did.
            enable_packed = (refusal is None
                             and os.environ.get("TORCHGWAS_NATIVE_STATS", "0") == "1")
        if enable_packed:
            self.native_transfer_dtype = np.dtype(np.uint8)
            self.native_row_width = ((self._n_samples + 3) // 4 + 63) // 64 * 64
            self.native_encoding = "pgen_2bit"
        # -9 is the int8 sentinel. It must NOT be applied to the uint8
        # dosage transport: -9 as a uint8 is 247, which is a legitimate
        # dosage code (247/127 = 1.94), so every masked call would be read
        # as 1.94 and every genuine 1.94 would be read as missing. 255 is
        # free there precisely because the scale is 127.0 rather than 127.5.
        if not getattr(self, "dosage_uint8_transport", False):
            self.native_missing_value = -9
        self.native_host_staging_copies = 0 if self._reorder_is_identity else 1
        self.marker_ids = self._all_metadata["marker_id"]
        self.chromosomes = self._all_metadata["chromosome"]
        self.positions = self._all_metadata["position"]
        self.effect_alleles = self._all_metadata["effect_allele"]
        self.other_alleles = self._all_metadata["other_allele"]
        # The metadata probe is not an iterator reader and need not remain open.
        self.close()

    @property
    def shape(self):
        return self._n_samples, self._raw_n_variants

    @property
    def genotype(self):
        return self

    @property
    def variant_metadata(self):
        return {key: self._all_metadata[key] for key in
                ("chromosome", "position", "effect_allele", "other_allele")}

    def _read_direct(self, reader, start, end, dtype):
        dtype = np.dtype(dtype)
        native_calls = self.mode == "hardcall" and dtype == np.dtype(np.int8)
        if not native_calls and dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
            raise ValueError("direct PGEN scans require float32 or float64 to preserve dosages and NaN")
        if not 0 <= start <= end <= self.shape[1]:
            raise IndexError("PGEN variant range is out of bounds")
        # Native variant-major storage allows a zero-copy transpose for the
        # samples x variants protocol; avoid a full strided output copy.
        values = np.empty((end - start, self._n_samples), dtype=(
            np.int8 if self.mode == "hardcall" else np.float32))
        if self.mode == "hardcall":
            reader.read_range(start, end, values, allele_idx=1, sample_maj=False)
        else:
            reader.read_dosages_range(start, end, values, allele_idx=1, sample_maj=False)
        if native_calls:
            # pgenlib decodes a two-bit alphabet into 0/1/2/-9. The GPU scan
            # handles missingness; no host float cast or QC pass is required.
            if not self._reorder_is_identity:
                values = values[:, self._sample_reorder]
            return values.T
        # The common complete-call path requires two row reductions, without
        # allocating several sample x variant boolean arrays per chunk.
        if values.size:
            row_min = values.min(axis=1)
            row_max = values.max(axis=1)
            suspect = ((row_min < 0) | (row_max > 2) |
                       ~np.isfinite(row_min) | ~np.isfinite(row_max))
            if suspect.any():
                rows = np.flatnonzero(suspect)
                selected = values[rows]
                invalid = ((selected != -9) & ((selected < 0) |
                           (selected > 2) | ~np.isfinite(selected)))
                if invalid.any():
                    row, sample = np.argwhere(invalid)[0]
                    raise ValueError(f"invalid PGEN ALT1 value at variant {start + rows[row]}, selected sample {sample}")
            has_missing = bool((row_min == -9).any())
        else:
            has_missing = False
        values = values.astype(dtype, copy=False)
        if has_missing and not native_calls:
            values[values == -9] = np.nan
        if not self._reorder_is_identity:
            values = values[:, self._sample_reorder]
        return values.T

    @contextmanager
    def native_reader_session(self):
        """Yield a parallel-safe ``read_into(start, end, variant_major_out)``.

        Caller owns each writable C-contiguous native output until GPU transfer
        completes. Full/ordered samples decode directly into it. A reordered
        sample selection uses a reusable per-thread scratch and one reorder
        copy. Ordinary native missing calls remain -9, including float dosages.
        Opt-in pgen_2bit transport uses uint8 rows of native_row_width bytes,
        64-byte aligned, with little-endian two-bit ALT1 codes 0/1/2/3 where
        3 means missing. Public iterators remain unchanged. Each fill
        decodes exactly the requested compute tile; decode_batch_size applies
        only to the public iterators, which retain independent reblocking.

        Join the caller's fill pool before leaving this context. Exit also
        waits any active native calls before closing this session's readers.
        """
        local = threading.local()
        condition = threading.Condition()
        readers = []
        profiling = os.environ.get("TORCHGWAS_SCAN_PROFILE", "0") != "0"
        profile = dict(enabled=profiling, fill_calls=0, requested_variants=0,
                       fill_worker_seconds=0.0, fill_worker_cpu_seconds=0.0,
                       reader_count=0, reader_init_worker_seconds=0.0,
                       reader_init_worker_cpu_seconds=0.0,
                       fill_includes_reader_initialization=True)
        self.last_native_reader_profile = profile
        active = 0
        closing = False

        def read_into(start, end, out):
            nonlocal active
            if not 0 <= start <= end <= self.shape[1]:
                raise IndexError("PGEN variant range is out of bounds")
            if (not isinstance(out, np.ndarray) or out.dtype != self.native_transfer_dtype or
                    out.shape != (end - start, self.native_row_width) or
                    not out.flags.c_contiguous or not out.flags.writeable):
                raise ValueError("PGEN output must be writable C-contiguous variant-major native dtype with exact shape")
            if self.native_encoding == "pgen_2bit" and out.size and out.ctypes.data % 64:
                raise ValueError("packed PGEN output must be 64-byte aligned")
            with condition:
                if closing:
                    raise RuntimeError("PGEN native reader session is closed")
                active += 1
            fill_started = time.perf_counter() if profiling else 0.0
            fill_cpu_started = time.thread_time() if profiling else 0.0
            try:
                reader = getattr(local, "reader", None)
                if reader is None:
                    init_started = time.perf_counter() if profiling else 0.0
                    init_cpu_started = time.thread_time() if profiling else 0.0
                    reader = self._new_reader()
                    local.reader = reader
                    with condition:
                        readers.append(reader)
                        if profiling:
                            profile["reader_count"] += 1
                            profile["reader_init_worker_seconds"] += time.perf_counter() - init_started
                            profile["reader_init_worker_cpu_seconds"] += time.thread_time() - init_cpu_started
                if start == end:
                    return out
                if self.native_encoding == "pgen_2bit":
                    reader.read_packed_range_into(start, end, out)
                    return out
                target = out
                if not self._reorder_is_identity:
                    scratch = getattr(local, "scratch", None)
                    if scratch is None or scratch.shape[0] < end - start:
                        scratch = np.empty(out.shape, dtype=self.native_dtype)
                        local.scratch = scratch
                    target = scratch[:end - start]
                if self.mode == "hardcall":
                    reader.read_range(start, end, target, allele_idx=1, sample_maj=False)
                elif self.dosage_uint8_transport:
                    # pgenlib writes float32 and there is no uint8 entry
                    # point, so the quantisation happens here, once, into a
                    # reusable per-thread buffer. `out` is the uint8
                    # transport row; `floats` never leaves this worker.
                    floats = getattr(local, 'dosage_floats', None)
                    if floats is None or floats.shape[0] < end - start:
                        floats = np.empty((end - start, self._n_samples),
                                          dtype=np.float32)
                        local.dosage_floats = floats
                    values = floats[:end - start]
                    reader.read_dosages_range(start, end, values,
                                              allele_idx=1, sample_maj=False)
                    # 255 is the missing sentinel, which is why the scale is
                    # 127.0 and not zstd's 127.5: codes span 0..254 and the
                    # top code stays free, so a masked call keeps its identity
                    # and the per-variant residual df stays exact.
                    missing = np.isnan(values)
                    np.clip(values, 0.0, 2.0, out=values)
                    np.multiply(values, 127.0, out=values)
                    np.rint(values, out=values)
                    np.copyto(target, values, casting='unsafe')
                    target[missing] = 255
                else:
                    reader.read_dosages_range(start, end, target, allele_idx=1, sample_maj=False)
                if target is not out:
                    np.take(target, self._sample_reorder, axis=1, out=out, mode="clip")
                return out
            finally:
                with condition:
                    if profiling:
                        profile["fill_calls"] += 1
                        profile["requested_variants"] += int(end - start)
                        profile["fill_worker_seconds"] += time.perf_counter() - fill_started
                        profile["fill_worker_cpu_seconds"] += time.thread_time() - fill_cpu_started
                    active -= 1
                    condition.notify_all()
        try:
            yield read_into
        finally:
            with condition:
                closing = True
                condition.wait_for(lambda: active == 0)
            # Every reader belongs only to this session; none can still be in
            # native code or be used by a later pass when it is closed.
            first_error = None
            for reader in readers:
                try:
                    reader.close()
                except Exception as error:
                    if first_error is None:
                        first_error = error
            if first_error is not None:
                raise first_error

    def read_chunk(self, start, end, dtype=np.float32):
        reader = self._new_reader()
        try:
            return self._read_direct(reader, start, end, dtype)
        finally:
            reader.close()

    def iter_chunks(self, chunk_size, dtype=np.float32, prefetch_chunks=None,
                    reader_workers=None, variant_range=None):
        from .streaming import OrderedChunkLoader
        workers = self.reader_workers if reader_workers is None else reader_workers
        depth = self.prefetch_chunks if prefetch_chunks is None else prefetch_chunks
        readers = []
        lock = threading.Lock()
        local = threading.local()
        initialization_times = []

        def read(start, end, output_dtype):
            reader = getattr(local, "reader", None)
            if reader is None:
                reader_started = time.perf_counter()
                reader = self._new_reader()
                reader_seconds = time.perf_counter() - reader_started
                local.reader = reader
                with lock:
                    readers.append(reader)
                    initialization_times.append(reader_seconds)
            return self._read_direct(reader, start, end, output_dtype)

        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        loader = OrderedChunkLoader(self.shape[1], read, self.decode_batch_size, dtype, depth, workers)
        decoded = iter(loader)
        parts = deque()
        buffered = 0
        emitted = 0

        def take(count):
            nonlocal buffered
            blocks = []
            remaining = count
            while remaining:
                block = parts[0]
                use = min(remaining, block.shape[1])
                blocks.append(block[:, :use])
                if use == block.shape[1]:
                    parts.popleft()
                else:
                    parts[0] = block[:, use:]
                remaining -= use
                buffered -= use
            # Aligned decode/compute tiles remain zero-copy views. Only a
            # compute chunk crossing decode boundaries needs concatenation.
            return blocks[0] if len(blocks) == 1 else np.concatenate(blocks, axis=1)

        try:
            for start, end, block in decoded:
                parts.append(block)
                buffered += end - start
                while buffered >= chunk_size:
                    output = take(chunk_size)
                    yield emitted, emitted + chunk_size, output
                    emitted += chunk_size
            if buffered:
                output = take(buffered)
                yield emitted, emitted + output.shape[1], output
        finally:
            # Join decode workers before closing their independent readers.
            decoded.close()
            for reader in readers:
                reader.close()
            self.last_iteration_stats = {
                "reader_count": len(initialization_times),
                "reader_init_worker_seconds": sum(initialization_times),
                "reader_init_max_seconds": max(initialization_times, default=0.0),
            }




    def iter_native_chunks(self, chunk_size, prefetch_chunks=None, reader_workers=None,
                           variant_range=None):
        """Compact H2D payload; int8 hardcall missing -9 must become NaN on GPU."""
        return self.iter_chunks(chunk_size, self.native_dtype, prefetch_chunks,
                                reader_workers, variant_range=variant_range)
