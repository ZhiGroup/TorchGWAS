"""NumPy's transparent-hugepage advice, turned off while a run is active.

NumPy (1.22+, Linux) madvises MADV_HUGEPAGE on every array of 4 MB or more.
Under THP defrag 'madvise' (the kernel default) a page fault in such a region
compacts memory synchronously, and on a host whose NUMA nodes are full of page
cache the compaction mostly fails after the wait. the H100 host, 2026-09-27: nodes
1-3 at 0-2 GB free, 811M compaction stalls of which 98% failed.
- A fresh 1 GiB array took 10.8-120 s to touch with the advice, and 0.36-0.45 s
  without it (pinned to the node with 166 GB free).
- Converting 8.09M variant IDs to fixed-width unicode took 75-154 s, against
  0.46 s: the size of the 38-160 s post-scan tail of indexed runs, which
  converted them while publishing.
- Resident hugepages were worth at most ~1.6x on a streaming add (16.7-17.9
  against 10.9-15.5 GB/s), and nothing on random gathers.
A run allocates host arrays per chunk, part and publication, so the fault cost
dominates. The advice is off for the length of the call and restored after;
TORCHGWAS_NUMPY_HUGEPAGE=1 keeps NumPy's own setting. PyTorch's CPU allocator
does not madvise hugepages by default, and pinned buffers are not THP-backed.
"""
from __future__ import annotations

import os
import threading
from functools import wraps

_LOCK = threading.Lock()
_DEPTH = 0
_SAVED = None


def _setter():
    import numpy as np
    core = getattr(np, '_core', None)
    if core is None:
        core = np.core
    return getattr(core.multiarray, '_set_madvise_hugepage', None)


def numpy_hugepage_advice():
    """True, False, or None when this NumPy has no switch."""
    setter = _setter()
    if setter is None:
        return None
    value = setter(False)
    setter(value)
    return bool(value)


def _enter():
    global _DEPTH, _SAVED
    setter = None if os.environ.get('TORCHGWAS_NUMPY_HUGEPAGE') == '1' else _setter()
    if setter is None:
        return None
    with _LOCK:
        if _DEPTH == 0:
            _SAVED = setter(False)
        _DEPTH += 1
    return setter


def _exit(setter):
    global _DEPTH
    if setter is None:
        return
    with _LOCK:
        _DEPTH -= 1
        if _DEPTH == 0:
            setter(_SAVED)


def without_numpy_hugepages(function):
    """Run `function` with NumPy's hugepage advice off (reentrant, thread-counted)."""
    @wraps(function)
    def run(*args, **kwargs):
        setter = _enter()
        try:
            return function(*args, **kwargs)
        finally:
            _exit(setter)
    return run
