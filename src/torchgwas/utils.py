from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import torch


def choose_device(device: str = "auto") -> torch.device:
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda", torch.cuda.current_device())
        return torch.device("cpu")
    return torch.device(device)


def ensure_2d(array: np.ndarray, name: str) -> np.ndarray:
    if array.ndim == 1:
        return array.reshape(-1, 1)
    if array.ndim != 2:
        raise ValueError(f"{name} must be 1D or 2D, got shape {array.shape}")
    return array


_FINITE_BLOCK_ROWS = 4096


def validate_no_missing(array: np.ndarray, name: str) -> None:
    """Refuse a matrix with missing or non-finite entries.

    Checked in row blocks because `np.isfinite(array)` materialises a bool
    array the full size of the input -- 20 GB on a 600,000-voxel phenotype
    and 70 GB at 2,085,000 -- purely to reduce it to one flag. Blocking also
    stops at the first offending block instead of scanning to the end.
    """
    rows = array.shape[0] if array.ndim else 0
    for begin in range(0, rows, _FINITE_BLOCK_ROWS):
        block = array[begin:begin + _FINITE_BLOCK_ROWS]
        if not np.isfinite(block).all():
            raise ValueError(
                f"{name} contains missing/non-finite values; "
                "v0.1 requires complete matrices")


def check_aligned_rows(*arrays: tuple[str, np.ndarray]) -> int:
    counts = {name: array.shape[0] for name, array in arrays if array is not None}
    if len(set(counts.values())) != 1:
        raise ValueError(f"row-count mismatch across inputs: {counts}")
    return next(iter(counts.values()))


_COLUMN_STD_BLOCK_ROWS = 256
_COLUMN_STD_FAST_PATH_BYTES = 64 << 20


def column_std_mask(array: np.ndarray) -> np.ndarray:
    """Which columns have non-zero standard deviation.

    `np.nanstd` over a large float64 matrix makes four full-size temporaries
    (the NaN mask, a copy of the data, and two in-place passes over that copy)
    and sums them single-threaded. On a 22,250 x 32,768 phenotype -- 5.8 GB --
    that is 37 s per call, and `run_linear_gwas` calls it on every scan, which
    at that trait count is fifty times the scan itself. The same decision is
    made here in two read passes with cache-sized temporaries.

    The mask, not the standard deviation, is the result, and the mask is
    order-independent: the sum of squared deviations is zero exactly when every
    deviation is zero, whatever the order the squares are added in. The mean is
    numpy's own `sum / n`, bit-identical to the one `nanvar` forms on a finite
    array, so each deviation is the same number `np.nanstd` would square.

    The fast path is taken only for a finite 2-D float64 array above a size
    where the temporaries matter; everything else -- other dtypes, small
    arrays, anything carrying a NaN -- goes through `np.nanstd` as before.
    """
    array = np.asarray(array)
    if (array.ndim != 2 or array.dtype != np.float64
            or array.nbytes < _COLUMN_STD_FAST_PATH_BYTES):
        return np.nanstd(array, axis=0) > 0
    n_rows = array.shape[0]
    mean = np.sum(array, axis=0, keepdims=True) / n_rows
    if not np.isfinite(mean).all():
        # A NaN or infinity somewhere in the data: nanstd's missing-value
        # semantics apply, so let it make the decision.
        return np.nanstd(array, axis=0) > 0
    sum_squares = np.zeros(array.shape[1], dtype=np.float64)
    for start in range(0, n_rows, _COLUMN_STD_BLOCK_ROWS):
        deviation = array[start:start + _COLUMN_STD_BLOCK_ROWS] - mean
        np.multiply(deviation, deviation, out=deviation)
        sum_squares += deviation.sum(axis=0)
    return np.sqrt(sum_squares / n_rows) > 0


def chunk_bounds(total: int, chunk_size: int | None) -> list[tuple[int, int]]:
    chunk = chunk_size or min(total, 4096) or 1
    return [(start, min(total, start + chunk)) for start in range(0, total, chunk)]


def timestamp() -> float:
    return time.perf_counter()


def elapsed(start: float) -> float:
    return round(time.perf_counter() - start, 6)


def mkdir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


TEXT_OUTPUT_COMPRESSLEVEL = 1
"""Compression level for every gzipped text output this package writes.

`gzip.open` defaults to **9**, and that default was costing more than the
whole rest of the writer.  Measured on an idle H100 host, 2,000 variants by
K=512 = 1.02M output rows, best of two:

    level 9 (the old default)   21.57 s    56.1 MB
    level 6                     15.22 s    56.3 MB
    level 1                      8.58 s    60.8 MB
    no compression               7.03 s   144.8 MB

So level 9 was **67% of the wall**, and buying it over level 6 costs 42% more
time for 0.4% fewer bytes -- indefensible at any file size.  Level 1 is
**2.51x faster than level 9 for 8.4% more bytes**, which is the right trade
for a tool whose argument is throughput.  The floor without compression is
7.03 s, so level 1 spends only 1.55 s compressing and still cuts the file to
42% of the plain text.
"""


def write_json(data: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def as_list(values: np.ndarray | None, prefix: str) -> list[str]:
    if values is None:
        return []
    if values.dtype.kind in {"U", "S", "O"}:
        return [str(v) for v in values.tolist()]
    return [f"{prefix}{i}" for i in range(values.shape[0])]


def upper_tail_log10(p_values: np.ndarray, t=None, df=None) -> np.ndarray:
    """`-log10 p`, falling back to a log-space tail where `p` has underflowed.

    Clipping at `np.finfo(float).tiny` caps this at **307.6527**, and
    `2 * stdtr(df, -|t|)` reaches exactly zero from |t| = 38.354 at df = 20,000
    -- so without a fallback every association past that point reports the same
    number, and the ranking a min-P scan exists to produce stops existing
    precisely where the answer is.

    Passing `t` and `df` enables the fallback. It is deliberately *only* a
    fallback: where `p` survived it carries that variant's own degrees of
    freedom and is used unchanged, so no existing row moves. Where `p` is zero
    there is nothing left to preserve, and `tails.upper_tail_log10_from_t`
    supplies a finite, strictly increasing value instead.
    """
    clipped = np.clip(p_values, np.finfo(float).tiny, 1.0)
    result = -np.log10(clipped)
    if t is None or df is None:
        return result
    underflowed = np.asarray(p_values) <= 0.0
    if not np.any(underflowed):
        return result
    from .tails import upper_tail_log10_from_t

    t_array, df_array = np.broadcast_arrays(
        np.asarray(t, dtype=np.float64), np.asarray(df, dtype=np.float64))
    result = np.array(result, dtype=np.float64, copy=True)
    result[underflowed] = upper_tail_log10_from_t(t_array[underflowed],
                                                  df_array[underflowed])
    return result
