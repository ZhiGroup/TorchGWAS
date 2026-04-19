from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import torch


def choose_device(device: str = "auto") -> torch.device:
    if device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device)


def ensure_2d(array: np.ndarray, name: str) -> np.ndarray:
    if array.ndim == 1:
        return array.reshape(-1, 1)
    if array.ndim != 2:
        raise ValueError(f"{name} must be 1D or 2D, got shape {array.shape}")
    return array


def validate_no_missing(array: np.ndarray, name: str) -> None:
    if np.isnan(array).any():
        raise ValueError(f"{name} contains missing values; v0.1 requires complete matrices")


def check_aligned_rows(*arrays: tuple[str, np.ndarray]) -> int:
    counts = {name: array.shape[0] for name, array in arrays if array is not None}
    if len(set(counts.values())) != 1:
        raise ValueError(f"row-count mismatch across inputs: {counts}")
    return next(iter(counts.values()))


def column_std_mask(array: np.ndarray) -> np.ndarray:
    return np.nanstd(array, axis=0) > 0


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


def write_json(data: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def as_list(values: np.ndarray | None, prefix: str) -> list[str]:
    if values is None:
        return []
    if values.dtype.kind in {"U", "S", "O"}:
        return [str(v) for v in values.tolist()]
    return [f"{prefix}{i}" for i in range(values.shape[0])]


def upper_tail_log10(p_values: np.ndarray) -> np.ndarray:
    clipped = np.clip(p_values, np.finfo(float).tiny, 1.0)
    return -np.log10(clipped)

