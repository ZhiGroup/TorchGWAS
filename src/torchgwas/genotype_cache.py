"""Decode once, replay later phenotype-tile rounds from host memory.

When more phenotype tiles than GPUs are needed (voxel scale), each GPU scans
its tiles in turn and every round re-reads and re-decodes the whole genotype.
The data moved per round is small next to the phenotype panel
(docs/autotune_design_20260924.md, section 2); the repeated decode is not.
CachedFillSource wraps a native-fill source: the first time a variant range
is decoded, the transfer-form rows (int8 hard calls, uint8 dosage or packed
2-bit rows) are copied into a host array; later fills of cached rows are a
memory copy. Readers fill disjoint rows, and two tiles decoding the same
range in the first round write identical data, so no locking is needed on
the cache itself.
"""
from __future__ import annotations

import contextlib
import threading

import numpy as np


def available_host_bytes():
    try:
        with open('/proc/meminfo') as handle:
            for line in handle:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return None


class GenotypeFillCache:
    """Host rows for variants [first, last) in the source's transfer form."""

    def __init__(self, first, last, row_width, dtype):
        self.first, self.last = int(first), int(last)
        self.values = np.empty((self.last - self.first, int(row_width)), dtype=dtype)
        self.ready = np.zeros(self.last - self.first, dtype=bool)
        self.hit_rows = self.miss_rows = 0
        self._lock = threading.Lock()

    @property
    def nbytes(self):
        return self.values.nbytes

    def count(self, hit, rows):
        with self._lock:
            if hit:
                self.hit_rows += rows
            else:
                self.miss_rows += rows

    def audit(self):
        return dict(variants=[self.first, self.last], bytes=self.nbytes, cached_rows=int(self.ready.sum()),
                    hit_rows=self.hit_rows, miss_rows=self.miss_rows)


class CachedFillSource:
    """A source view whose native fills go through a shared GenotypeFillCache."""

    def __init__(self, source, cache):
        self._source, self._cache = source, cache

    def __getattr__(self, name):
        return getattr(self._source, name)

    @contextlib.contextmanager
    def native_reader_session(self):
        cache, source = self._cache, self._source
        state = dict(session=None, fill=None)
        opening = threading.Lock()

        def inner():
            with opening:
                if state['fill'] is None:  # the real reader only when something misses
                    state['session'] = source.native_reader_session()
                    state['fill'] = state['session'].__enter__()
            return state['fill']

        def fill(start, end, out):
            lo, hi = start - cache.first, end - cache.first
            if 0 <= lo and hi <= len(cache.ready) and cache.ready[lo:hi].all():
                np.copyto(out, cache.values[lo:hi])
                cache.count(True, end - start)
                return
            inner()(start, end, out)
            if 0 <= lo and hi <= len(cache.ready):
                cache.values[lo:hi] = out
                cache.ready[lo:hi] = True
            cache.count(False, end - start)

        try:
            yield fill
        finally:
            if state['session'] is not None:
                state['session'].__exit__(None, None, None)


def fill_cache_for(genotype, variant_range, *, share=0.5):
    """A cache for this source's native fills, or None (unsupported or too large)."""
    if not (getattr(genotype, 'allows_direct_native_fill', False) and hasattr(genotype, 'native_reader_session')):
        return None
    dtype = np.dtype(getattr(genotype, 'native_transfer_dtype', getattr(genotype, 'native_dtype', np.int8)))
    width = int(getattr(genotype, 'native_row_width', 0) or genotype.shape[0])
    first, last = variant_range
    need = (last - first) * width * dtype.itemsize
    free = available_host_bytes()
    if free is None or need > share * free:
        return None
    return GenotypeFillCache(first, last, width, dtype)
