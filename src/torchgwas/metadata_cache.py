"""Binary sidecar caches for text metadata that is re-parsed on every open.

Every input format in this project carries a text metadata file whose parse is
paid in full on each open, regardless of how much of the genotype is then read,
and each one has been found the same way -- as a flat cost sitting in the
`open` column of a benchmark that nothing else explained:

| format | file | rows | parse |
| --- | --- | --- | --- |
| PLINK | `.bim` | 8.93M | cached since earlier work |
| BGEN | `.bgi` (SQLite) | 8.93M | 19.2 s, cached, now 1.3 s |
| PGEN | `.pvar` | 8.93M | 6.09 s, cached, now 0.17 s |
| zstd store | `.variants.tsv` | 8.93M | 6.33 s of a 6.53 s open |

This module is the shared mechanism, extracted after the third one so the
fourth does not become a fourth copy. It is deliberately small: identify the
source, store named arrays, refuse anything that does not match.

**Two decisions worth keeping.** Text columns are stored as fixed-width unicode
rather than pickled objects -- the BGEN cache was first written with bytes and
a warm open took 11.7 s against 1.3 s, which made it nearly worthless. And
computing a cache path never creates a directory: a caller merely looking for a
cache would otherwise bring it into existence, which a test asserting that an
unused cache directory stays unused caught.
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np

SCHEMA = 1


def cache_path(cache_dir: str | Path, source: Path, tag: str,
               *, create: bool = False) -> Path:
    """Where `source`'s cache lives, keyed on its identity and this schema.

    The key includes size and modification time, so a rewritten source simply
    gets a different path and the stale entry is never consulted again rather
    than being detected and repaired.
    """
    stat = source.stat()
    identity = json.dumps(
        {
            "schema": SCHEMA,
            "tag": tag,
            "path": str(source.resolve()),
            "size": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
        },
        sort_keys=True,
    )
    digest = hashlib.sha1(identity.encode("utf-8")).hexdigest()[:12]
    directory = Path(cache_dir)
    if create:
        directory.mkdir(parents=True, exist_ok=True)
    return directory / f"{source.stem}_{digest}.{tag}.cache"


def store_arrays(cache_dir: str | Path, source: str | Path, tag: str,
                 arrays: dict[str, np.ndarray]) -> Path:
    """Write named arrays for `source`, atomically."""
    source = Path(source)
    target = cache_path(cache_dir, source, tag, create=True)
    if target.is_dir():
        return target
    lengths = {value.shape[0] for value in arrays.values()}
    if len(lengths) != 1:
        raise ValueError(f"{tag} cache arrays have inconsistent lengths")
    stat = source.stat()
    # Built in a temporary directory and renamed into place, so a reader never
    # sees a half-written cache and two writers cannot interleave.
    with tempfile.TemporaryDirectory(prefix=f".{target.name}.",
                                     dir=target.parent) as temporary_name:
        temporary = Path(temporary_name)
        for name, array in arrays.items():
            np.save(temporary / f"{name}.npy", array, allow_pickle=False)
        (temporary / "manifest.json").write_text(
            json.dumps(
                {
                    "schema": SCHEMA,
                    "tag": tag,
                    "source_path": str(source.resolve()),
                    "source_size": int(stat.st_size),
                    "source_mtime_ns": int(stat.st_mtime_ns),
                    "rows": next(iter(lengths)),
                    "arrays": sorted(arrays),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        temporary.replace(target)
    return target


def load_arrays(cache_dir: str | Path, source: str | Path, tag: str,
                names):
    """Named arrays for `source`, or None when there is no usable cache.

    Arrays come back memory-mapped and read-only, so a large cache costs
    address space rather than resident memory.
    """
    source = Path(source)
    target = cache_path(cache_dir, source, tag)
    manifest_path = target / "manifest.json"
    if not target.is_dir() or not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
    except ValueError:
        return None
    stat = source.stat()
    valid = (
        int(manifest.get("schema", -1)) == SCHEMA
        and str(manifest.get("tag")) == tag
        and str(manifest.get("source_path")) == str(source.resolve())
        and int(manifest.get("source_size", -1)) == stat.st_size
        and int(manifest.get("source_mtime_ns", -1)) == stat.st_mtime_ns
    )
    if not valid:
        return None
    arrays = {}
    for name in names:
        array_path = target / f"{name}.npy"
        if not array_path.is_file():
            return None
        arrays[name] = np.load(array_path, mmap_mode="r", allow_pickle=False)
    lengths = {array.shape[0] for array in arrays.values()}
    if len(lengths) != 1 or next(iter(lengths)) != int(manifest["rows"]):
        raise ValueError(f"{tag} cache length does not match its manifest: {target}")
    return arrays


def storable(array) -> np.ndarray:
    """An array `np.save` will accept without pickling.

    Parsed text columns arrive as object arrays of Python strings, and
    `np.save(..., allow_pickle=False)` refuses those outright -- so a cache
    that stored them verbatim would raise on every write. Fixed-width unicode
    is also the fast choice on the way back in: the BGEN cache was first
    written with bytes and a warm open took 11.7 s against 1.3 s.
    """
    array = np.asarray(array)
    return np.asarray(array, dtype=str) if array.dtype.kind in "OUSV" else array


def cached_arrays(cache_dir: str | Path | None, source: str | Path, tag: str,
                  names, parse):
    """`parse()`'s result, from cache when possible and stored when not.

    A cache that cannot be written is a missed optimisation rather than a
    failed run: a read-only or full cache directory must not stop work that has
    already parsed everything it needs. That tolerance covers the filesystem,
    not the data -- an array that cannot be represented is a bug here and is
    left to raise.
    """
    if cache_dir is None:
        return parse()
    cached = load_arrays(cache_dir, source, tag, names)
    if cached is not None:
        return cached
    parsed = parse()
    try:
        store_arrays(cache_dir, source, tag,
                     {name: storable(parsed[name]) for name in names})
    except OSError:
        pass
    return parsed
