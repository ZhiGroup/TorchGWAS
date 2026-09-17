"""ctypes binding for the native PGEN record decoder.

The Python reader in :mod:`torchgwas.pgen_reader` stays the reference
implementation and the oracle-validated definition of correct. This module
binds the C decoder that does the same work fast enough to ship, and keeps
the error codes distinguishable so an unsupported file is reported as
unsupported rather than as corruption.

ctypes releases the GIL for the duration of the call, so disjoint variant
runs decode concurrently from separate threads.
"""

from __future__ import annotations

import ctypes
import os
from pathlib import Path

import numpy as np

ABI_VERSION = 1
LIBRARY_NAME = "libtorchgwas_pgen.so"

PGEN_OK = 0
_ERRORS = {
    -1: ("truncated", "record ended before the decoder expected"),
    -2: ("unsupported", "record carries a multiallelic, phase or dosage track"),
    -3: ("unsupported", "record uses a reserved storage form"),
    -4: ("corrupt", "difflist names a sample outside the cohort"),
    -5: ("corrupt", "LD-compressed record with no base variant in hand"),
    -6: ("unsupported", "unreadable one-bit mode byte"),
    -7: ("corrupt", "variable-length integer overflowed"),
}


class PgenNativeUnavailable(RuntimeError):
    """The shared library is not built or not loadable."""


class PgenUnsupported(ValueError):
    """The file is valid PGEN but outside this decoder's scope."""


class PgenCorrupt(ValueError):
    """The file does not decode as PGEN."""


def _candidate_paths() -> list[Path]:
    override = os.environ.get("TORCHGWAS_PGEN_LIBRARY")
    if override:
        return [Path(override)]
    root = Path(__file__).resolve().parents[2]
    return [
        root / ".build-libs" / LIBRARY_NAME,
        root / LIBRARY_NAME,
        Path(LIBRARY_NAME),
    ]


_library = None


def load_library():
    """Load and memoize the decoder, verifying the ABI it was built with."""
    global _library
    if _library is not None:
        return _library
    errors = []
    for path in _candidate_paths():
        try:
            handle = ctypes.CDLL(str(path))
        except OSError as exc:
            errors.append(f"{path}: {exc}")
            continue
        handle.torchgwas_pgen_abi.restype = ctypes.c_int
        handle.torchgwas_pgen_abi.argtypes = []
        found = handle.torchgwas_pgen_abi()
        if found != ABI_VERSION:
            raise PgenNativeUnavailable(
                f"{path} reports ABI {found}, expected {ABI_VERSION}; rebuild with "
                "build_pgen_decode.sh"
            )
        handle.torchgwas_pgen_decode_range.restype = ctypes.c_int
        handle.torchgwas_pgen_decode_range.argtypes = [
            ctypes.POINTER(ctypes.c_uint8),   # records
            ctypes.POINTER(ctypes.c_uint64),  # offsets
            ctypes.POINTER(ctypes.c_uint32),  # lengths
            ctypes.POINTER(ctypes.c_uint8),   # vrtypes
            ctypes.c_uint64,                  # count
            ctypes.c_uint64,                  # sample_ct
            ctypes.c_uint32,                  # id_bytes
            ctypes.POINTER(ctypes.c_uint8),   # out
            ctypes.c_uint64,                  # out_stride
            ctypes.POINTER(ctypes.c_uint8),   # ld_base
            ctypes.POINTER(ctypes.c_int),     # have_ld_base
            ctypes.POINTER(ctypes.c_uint64),  # failed_index
        ]
        handle.torchgwas_pgen_expand.restype = ctypes.c_int
        handle.torchgwas_pgen_expand.argtypes = [
            ctypes.POINTER(ctypes.c_uint8),
            ctypes.c_uint64,
            ctypes.c_uint64,
            ctypes.c_uint64,
            ctypes.POINTER(ctypes.c_uint8),
        ]
        # Bound optionally, not through the ABI check: a library built before
        # this entry point existed is still correct, just slower, and failing
        # to load it outright would strand anyone with a stale build. Callers
        # ask `expand_hardcall_available()` and fall back.
        try:
            handle.torchgwas_pgen_expand_hardcall.restype = ctypes.c_int
            handle.torchgwas_pgen_expand_hardcall.argtypes = [
                ctypes.POINTER(ctypes.c_uint8),
                ctypes.c_uint64,
                ctypes.c_uint64,
                ctypes.c_uint64,
                ctypes.POINTER(ctypes.c_int8),
                ctypes.c_int8,
            ]
        except AttributeError:
            pass
        _library = handle
        return handle
    raise PgenNativeUnavailable(
        "could not load the PGEN decoder; run build_pgen_decode.sh. Tried:\n  "
        + "\n  ".join(errors)
    )


def available() -> bool:
    try:
        load_library()
    except PgenNativeUnavailable:
        return False
    return True


def _raise(status: int, variant_index: int) -> None:
    kind, message = _ERRORS.get(status, ("corrupt", f"decoder returned {status}"))
    text = f"variant {variant_index}: {message}"
    if kind == "unsupported":
        raise PgenUnsupported(
            text + "; use the pgenlib backend for this file"
        )
    raise PgenCorrupt(text)


def _u8(array: np.ndarray):
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))


def decode_range(
    records: np.ndarray,
    offsets: np.ndarray,
    lengths: np.ndarray,
    vrtypes: np.ndarray,
    sample_ct: int,
    id_bytes: int,
    out: np.ndarray,
    ld_base: np.ndarray,
    have_ld_base: bool,
    first_variant: int = 0,
) -> bool:
    """Decode consecutive records into packed rows of ``out``.

    ``offsets`` are byte positions within ``records``. ``out`` is
    ``(count, stride)`` uint8. ``ld_base`` is updated in place; the returned
    flag says whether it now holds a usable base, so a caller can decode a
    long run in pieces.
    """
    handle = load_library()
    count = len(offsets)
    if out.shape[0] < count:
        raise ValueError(f"out holds {out.shape[0]} rows, need {count}")
    for name, array, dtype in (
        ("records", records, np.uint8),
        ("offsets", offsets, np.uint64),
        ("lengths", lengths, np.uint32),
        ("vrtypes", vrtypes, np.uint8),
        ("out", out, np.uint8),
        ("ld_base", ld_base, np.uint8),
    ):
        if array.dtype != dtype:
            raise TypeError(f"{name} must be {np.dtype(dtype)}, got {array.dtype}")
        if not array.flags["C_CONTIGUOUS"]:
            raise ValueError(f"{name} must be C-contiguous")

    flag = ctypes.c_int(1 if have_ld_base else 0)
    failed = ctypes.c_uint64(0)
    status = handle.torchgwas_pgen_decode_range(
        _u8(records),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        lengths.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
        _u8(vrtypes),
        ctypes.c_uint64(count),
        ctypes.c_uint64(sample_ct),
        ctypes.c_uint32(id_bytes),
        _u8(out),
        ctypes.c_uint64(out.strides[0]),
        _u8(ld_base),
        ctypes.byref(flag),
        ctypes.byref(failed),
    )
    if status != PGEN_OK:
        _raise(status, first_variant + int(failed.value))
    return bool(flag.value)


def expand_hardcall_available() -> bool:
    """True when the built library can emit signed hard calls in one pass."""
    try:
        handle = load_library()
    except PgenNativeUnavailable:
        return False
    return hasattr(handle, "torchgwas_pgen_expand_hardcall")


def expand_hardcall(packed: np.ndarray, sample_ct: int, out: np.ndarray,
                    missing: int = -9) -> np.ndarray:
    """Expand packed rows directly into `out` as 0/1/2/`missing` int8.

    One pass, and no wide intermediate: the remap this replaces was 85% of a
    native read.
    """
    handle = load_library()
    packed = np.ascontiguousarray(packed, dtype=np.uint8)
    if out.dtype != np.int8 or not out.flags["C_CONTIGUOUS"]:
        raise TypeError("out must be a C-contiguous int8 array")
    if out.shape != (packed.shape[0], sample_ct):
        raise ValueError(
            f"out must be {(packed.shape[0], sample_ct)}, got {out.shape}")
    handle.torchgwas_pgen_expand_hardcall(
        _u8(packed),
        ctypes.c_uint64(packed.shape[0]),
        ctypes.c_uint64(packed.strides[0]),
        ctypes.c_uint64(sample_ct),
        out.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)),
        ctypes.c_int8(missing),
    )
    return out


def expand(packed: np.ndarray, sample_ct: int) -> np.ndarray:
    """Expand packed rows to one uint8 category per sample."""
    handle = load_library()
    packed = np.ascontiguousarray(packed, dtype=np.uint8)
    count = packed.shape[0]
    out = np.empty((count, sample_ct), dtype=np.uint8)
    handle.torchgwas_pgen_expand(
        _u8(packed),
        ctypes.c_uint64(count),
        ctypes.c_uint64(packed.strides[0]),
        ctypes.c_uint64(sample_ct),
        _u8(out),
    )
    return out
