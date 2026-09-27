"""Independent generic-buffer CPU prototype; no published calibration prices."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import statistics
import subprocess
import time
import numpy as np
from torchgwas.host_significance import select_host_pairs, ceil_float32, predicate_block_shape


def main():
    root = Path('results/host_predicate_prototype_20260922')
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).with_suffix('.cpp')
    flags = ['-std=c++17', '-O3', '-fPIC', '-shared', '-fno-fast-math', '-ffp-contract=off']
    library = root / 'predicate.so'
    subprocess.run(['c++', *flags, str(source), '-o', str(library)], check=True)
    native = ctypes.CDLL(str(library.resolve())).host_predicate
    native.argtypes = [ctypes.c_void_p] * 3 + [ctypes.c_size_t] * 2
    native.restype = None
    def predicate(values, limits, mask):
        native(values.ctypes.data, limits.ctypes.data, mask.ctypes.data, *values.shape)
    def numpy_predicate(values, limits, mask):
        height, width, _ = predicate_block_shape(*values.shape)
        for row in range(0, values.shape[0], height):
            for column in range(0, values.shape[1], max(1, width)):
                sl = np.s_[row:row+height, column:column+width]
                absolute = np.abs(values[sl])
                np.greater_equal(absolute, limits[row:row+height], out=mask[sl])
                mask[sl] &= np.isfinite(absolute)
    def select(values, df, critical):
        limits = ceil_float32(critical)
        mask = np.empty(values.shape, bool)
        predicate(values, limits, mask)
        rows = np.flatnonzero(mask).astype(np.int64, copy=False)
        columns = np.empty_like(rows)
        if values.shape[1]: np.divmod(rows, values.shape[1], out=(rows, columns))
        return rows, columns, values[rows, columns], values[rows, columns], df[rows, 0]
    rng = np.random.default_rng(922921)
    # Exact independent reference including float32 thresholds and nonfinite values.
    limits = np.asarray([0., -0., -np.inf, np.nan, np.inf, 1e-50, 1e-40,
        1., np.nextafter(1., 2.), 7.123456789, np.finfo(np.float32).max, 1e40])[:, None]
    with np.errstate(over='ignore', invalid='ignore'):
        center = limits.astype(np.float32)
        values = np.concatenate([center, -center, np.nextafter(center, np.float32(-np.inf)),
            np.nextafter(center, np.float32(np.inf)), np.full_like(center, np.nan),
            np.full_like(center, np.inf), np.full_like(center, -np.inf)], axis=1)
    df = np.arange(len(limits), dtype=np.float32)[:, None] + 10
    for got, expected in zip(select(values, df, limits), select_host_pairs(values, values, df, limits)):
        np.testing.assert_array_equal(got, expected)
    records = []
    order_rng = random.Random(92333)
    for shape in [(256, 4096), (1024, 8192), (1024, 8193)]:
        ring = [rng.standard_normal(shape, dtype=np.float32) for _ in range(4)]
        df = np.full((shape[0], 1), 100., np.float32)
        mask = np.empty(shape, bool)
        for limit, density in [(7., 'empty'), (2., 'sparse'), (0., 'dense')]:
            critical = np.full((shape[0], 1), limit, np.float64)
            critical32 = ceil_float32(critical)
            for value in ring:
                for got, expected in zip(select(value, df, critical), select_host_pairs(value, value, df, critical)):
                    np.testing.assert_array_equal(got, expected)
            for repeat in range(8):
                values = ring[repeat % len(ring)]
                cases = ['numpy_predicate', 'native_predicate', 'numpy_selector', 'native_selector']
                order_rng.shuffle(cases)
                for case in cases:
                    faults = resource.getrusage(resource.RUSAGE_THREAD)
                    began = time.perf_counter(); cpu = time.thread_time()
                    if case == 'numpy_predicate': numpy_predicate(values, critical32, mask)
                    elif case == 'native_predicate': predicate(values, critical32, mask)
                    elif case == 'numpy_selector': result = select_host_pairs(values, values, df, critical)
                    else: result = select(values, df, critical)
                    cpu = time.thread_time() - cpu; wall = time.perf_counter() - began
                    after = resource.getrusage(resource.RUSAGE_THREAD)
                    records.append(dict(shape=shape, density=density, case=case, repeat=repeat,
                        cpu_seconds=cpu, wall_seconds=wall, minor_faults=after.ru_minflt-faults.ru_minflt))
                    if case.endswith('selector'): del result
            print(shape, density, {case: statistics.median(r['cpu_seconds'] for r in records
                if r['shape']==shape and r['density']==density and r['case']==case) for case in cases}, flush=True)
    report = dict(records=records, correctness='exact arrays vs existing selector on every buffer and boundary case',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0], flags=flags,
        library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(), numpy=np.__version__,
        affinity=sorted(os.sched_getaffinity(0)), scope=__doc__)
    (root/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__=='__main__': main()
