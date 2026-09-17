"""A pgenlib-shaped reader backed by our own PGEN decoder.

`PgenDosageSource` talks to exactly five methods of a pgenlib reader --
`get_raw_sample_ct`, `get_variant_ct`, `read_range`, `read_packed_range_into`
and `close` -- so the cheapest way to stop requiring pgenlib at runtime is to
present those five from the decoder we already wrote and validated, rather than
to thread a second backend through the read path.

**Why this exists.** pgenlib is an optional C extension and it is installed in
none of the conda environments on the lab hosts, so PGEN input currently cannot
be opened at all, while the decoder that would open it ships in the package and
is reachable only from two benchmarks. The exactness check against pgenlib
stays valuable and stays in the tests; it is the *runtime* dependency that goes.

**Scope is decided once per file, not per variant.** The real dosage PGEN has
8,930,997 of 8,931,083 variants carrying a dosage track and 86 that do not, so a
reader that refused variant by variant would decode 86 records and reject the
rest -- a uselessly partial answer. :func:`torchgwas.pgen_reader.file_scope`
reads only the header and record index, so the decision is cheap even on a
whole-cohort file, and a file outside scope falls back to pgenlib with its
reason named.

**LD-compressed records force a replay.** Forms 2 and 3 are expressed against
the most recent record that was not itself LD-compressed, so a read entering the
file at an arbitrary variant backs up to that record, decodes forward, and
discards the rows before the one it wanted. Variant blocks never open with an
LD-compressed record, so the walk is bounded and in practice is a few variants.
"""

from __future__ import annotations

import os

import numpy as np

from . import pgen_native
from .pgen_reader import (
    _bytes_per_sample_id,
    file_scope,
    ld_safe_start,
    read_header,
)

MISSING_CODE = 3
MISSING_VALUE = -9

# Codes 0, 1, 2 are the ALT1 count and 3 is missing; the table is padded to 256
# so an out-of-range byte cannot index past its end. Applying it with
# `np.take(..., out=)` turns the remap into a single pass that writes straight
# into the caller's array. The obvious spelling -- view as int8, copy, then
# `putmask` where the code is 3 -- costs four passes over the block and a bool
# temporary the size of the genotypes.
_HARDCALL_TABLE = np.zeros(256, dtype=np.int8)
_HARDCALL_TABLE[:4] = (0, 1, 2, MISSING_VALUE)


class NativePgenReader:
    """Decode a biallelic hard-call PGEN without pgenlib.

    The constructor takes pgenlib's keyword names so it is a drop-in for
    `pgenlib.PgenReader` at the one call site that builds readers.
    """

    def __init__(self, path, raw_sample_ct=None, variant_ct=None,
                 sample_subset=None):
        self.path = os.fsdecode(path)
        self.header = read_header(self.path)
        self.sample_ct = int(self.header.sample_ct)
        self.variant_ct = int(self.header.variant_ct)
        if raw_sample_ct is not None and int(raw_sample_ct) != self.sample_ct:
            raise ValueError(
                f"PSAM declares {int(raw_sample_ct)} samples, "
                f"the PGEN header says {self.sample_ct}")
        if variant_ct is not None and int(variant_ct) != self.variant_ct:
            raise ValueError(
                f"PVAR declares {int(variant_ct)} variants, "
                f"the PGEN header says {self.variant_ct}")
        self.genovec_bytes = (self.sample_ct + 3) // 4
        self._id_bytes = _bytes_per_sample_id(self.sample_ct)
        self._vrtypes = np.ascontiguousarray(self.header.vrtypes, dtype=np.uint8)
        self._offsets = np.ascontiguousarray(self.header.record_offsets, dtype=np.uint64)
        self._lengths = np.ascontiguousarray(self.header.record_lengths, dtype=np.uint32)
        self._subset = (None if sample_subset is None
                        else np.ascontiguousarray(sample_subset, dtype=np.int64))
        # One descriptor per reader, read positionally: readers are created per
        # worker thread, and preadv carries its own offset so no seek is shared.
        self._fd = os.open(self.path, os.O_RDONLY)
        self._blob = None

    # -- pgenlib surface -------------------------------------------------

    def get_raw_sample_ct(self) -> int:
        return self.sample_ct

    def get_variant_ct(self) -> int:
        return self.variant_ct

    def close(self) -> None:
        fd = getattr(self, "_fd", None)
        if fd is not None:
            os.close(fd)
            self._fd = None

    def __enter__(self) -> "NativePgenReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # noqa: BLE001 - interpreter teardown
            pass

    def read_range(self, start, end, out, allele_idx=1, sample_maj=False):
        """Fill `out` with ALT1 hard calls: 0, 1, 2, or -9 for a missing call."""
        if allele_idx != 1:
            raise ValueError(
                "the native PGEN decoder reads biallelic records, so only "
                f"allele_idx=1 is meaningful; got {allele_idx}")
        if sample_maj:
            raise ValueError(
                "the native PGEN decoder emits variant-major rows; "
                "sample_maj=True is not supported")
        start, end = int(start), int(end)
        width = self.sample_ct if self._subset is None else int(self._subset.size)
        if out.shape != (end - start, width):
            raise ValueError(
                f"out must be {(end - start, width)}, got {out.shape}")
        if end == start:
            return out
        packed = self._decode_packed(start, end)
        if self._subset is None and pgen_native.expand_hardcall_available():
            # The whole-cohort case, which is the common one: the decoder emits
            # the final signed values straight into the caller's array, so there
            # is no wide intermediate and no second pass at all.
            pgen_native.expand_hardcall(packed, self.sample_ct, out,
                                        MISSING_VALUE)
            return out
        categories = pgen_native.expand(packed, self.sample_ct)
        if self._subset is not None:
            categories = categories[:, self._subset]
        # Fallback for a sample subset, or a library built before the one-pass
        # entry point existed. `mode='clip'` is not a safety net, it is the
        # point: the default `mode='raise'` bounds-checks every index and
        # measured *slower than the code this replaced*, while clip is not. The
        # check can never fire, since a uint8 index cannot leave a 256-entry
        # table, so clipping removes a test whose outcome is already known.
        np.take(_HARDCALL_TABLE, categories, out=out, mode='clip')
        return out

    def read_packed_range_into(self, start, end, out):
        """Fill `out` with raw two-bit rows: 0 ref, 1 het, 2 alt, 3 missing."""
        if self._subset is not None:
            raise ValueError(
                "packed PGEN transport requires every sample in file order; "
                "this reader was opened with a sample subset")
        start, end = int(start), int(end)
        if out.shape[0] != end - start:
            raise ValueError(f"out holds {out.shape[0]} rows, need {end - start}")
        if out.shape[1] < self.genovec_bytes:
            raise ValueError(
                f"out rows are {out.shape[1]} bytes, need at least "
                f"{self.genovec_bytes}")
        if end == start:
            return out
        # Straight into the caller's buffer. This used to decode into a fresh
        # array and then copy the rows across, which over one sweep of the
        # benchmark cohort is 49.7 GB of zero-fill plus 49.7 GB of copy --
        # against 22.6 GB of input actually read. The decoder already takes a
        # row stride, so the padded transport row costs it nothing, and it
        # zeroes the padding itself: leaving that undefined would make two runs
        # of the same scan differ byte for byte on identical genotypes.
        self._decode_into(start, end, out)
        return out

    # -- decode ----------------------------------------------------------

    def _records_for(self, safe: int, end: int):
        """The record bytes for `[safe, end)`, with offsets relative to them.

        The bytes land in a buffer owned by this reader and are overwritten by
        the next call, which is safe because every caller decodes before it
        returns. `os.pread` would allocate a fresh object per chunk instead --
        about 48 MB at the default chunk size, malloc'd, faulted in and freed
        446 times over a sweep of this cohort. Readers are per worker thread,
        so the buffer is not shared.
        """
        offsets = self._offsets[safe:end]
        lengths = self._lengths[safe:end]
        begin = int(offsets[0])
        span = int(offsets[-1]) + int(lengths[-1]) - begin
        if self._blob is None or self._blob.size < span:
            # Grown, never shrunk: chunk sizes are fixed for a run, so this
            # settles after the first chunk.
            self._blob = np.empty(span, dtype=np.uint8)
        blob = self._blob[:span]
        filled = 0
        while filled < span:
            # A single preadv is not promised to return everything asked for,
            # and a short read here would otherwise surface as "truncated" on a
            # file that is perfectly intact.
            read = os.preadv(self._fd, [blob[filled:]], begin + filled)
            if read <= 0:
                raise ValueError(
                    f"PGEN records [{safe}, {end}) are truncated: "
                    f"wanted {span} bytes, read {filled}")
            filled += read
        return (
            blob,
            np.ascontiguousarray(offsets - np.uint64(begin), dtype=np.uint64),
            np.ascontiguousarray(lengths, dtype=np.uint32),
            np.ascontiguousarray(self._vrtypes[safe:end], dtype=np.uint8),
        )

    def _decode_into(self, start: int, end: int, out: np.ndarray) -> np.ndarray:
        """Decode `[start, end)` into `out`, one row per variant.

        `out` rows may be wider than a genovec; the decoder zeroes the padding.

        The LD prefix is handled separately rather than by decoding it into the
        caller's buffer and slicing it off afterwards, because the caller's
        buffer has room for the requested variants only. Records before `start`
        exist solely to rebuild the LD base, so they go to a scratch row that is
        thrown away; `decode_range` is built to resume from a base it was
        handed, which is exactly what that needs.
        """
        if not 0 <= start <= end <= self.variant_ct:
            raise IndexError(
                f"PGEN variant range [{start}, {end}) is outside "
                f"[0, {self.variant_ct})")
        if end == start:
            return out
        safe = ld_safe_start(self._vrtypes, start)
        blob, offsets, lengths, vrtypes = self._records_for(safe, end)
        base = np.zeros(self.genovec_bytes, dtype=np.uint8)
        have_base = False
        prefix = start - safe
        if prefix:
            scratch = np.empty((prefix, self.genovec_bytes), dtype=np.uint8)
            have_base = pgen_native.decode_range(
                blob, offsets[:prefix], lengths[:prefix], vrtypes[:prefix],
                self.sample_ct, self._id_bytes, scratch, base,
                have_ld_base=False, first_variant=safe)
        pgen_native.decode_range(
            blob, offsets[prefix:], lengths[prefix:], vrtypes[prefix:],
            self.sample_ct, self._id_bytes, out, base,
            have_ld_base=have_base, first_variant=start)
        return out

    def _decode_packed(self, start: int, end: int) -> np.ndarray:
        """Packed rows for `[start, end)`, replaying any LD prefix first.

        `np.empty` rather than `np.zeros`: every record form writes the whole
        genovec, and a decode that fails raises rather than returning a
        partially filled array.
        """
        out = np.empty((end - start, self.genovec_bytes), dtype=np.uint8)
        return self._decode_into(start, end, out)


def native_backend_reason(path) -> str | None:
    """Why the native backend cannot read this file, or None if it can.

    Returns a sentence, not a flag, because the two reasons a caller can act on
    are different: an unbuilt shared library is fixed by building it, and a file
    carrying dosage or phase tracks is not fixable and needs pgenlib.
    """
    if not pgen_native.available():
        return ("the native PGEN decoder is not built "
                "(run build_pgen_decode.sh, or set TORCHGWAS_PGEN_LIBRARY)")
    try:
        scope = file_scope(path)
    except Exception as exc:  # noqa: BLE001 - any parse failure disqualifies it
        return f"the PGEN index did not parse: {type(exc).__name__}: {exc}"
    if not scope.supported:
        return f"the file is outside the native decoder's scope: {scope.reason}"
    return None
