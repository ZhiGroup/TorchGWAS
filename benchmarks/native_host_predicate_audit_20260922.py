"""Independent production-selector controls, including nonempty output.

Counterbalanced generic pinned buffers. Not a GWAS speedup or published price.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import statistics
import time
import numpy as np
import torch
from torchgwas.host_significance import fill_predicate_mask, select_host_pairs, ceil_float32
from torchgwas.native_host_predicate import context
from torchgwas.detailed_calibration import source_identity, _numpy_core_context


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='results/native_host_predicate_audit_20260922')
    parser.add_argument('--concurrent-control', default=None)
    parser.add_argument('--torch-threads', type=int)
    args = parser.parse_args()
    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=False)
    if args.torch_threads is not None: torch.set_num_threads(args.torch_threads)
    source = source_identity(); native_context = context(); numpy_core = _numpy_core_context()
    settings = dict(numpy_core=numpy_core,numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
        torch_threads=torch.get_num_threads(), environment={key:os.getenv(key) for key in
            ['NUMPY_MADVISE_HUGEPAGE','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'OMP_WAIT_POLICY','GOMP_SPINCOUNT']})
    records = []; rng = random.Random(924911)
    for shape in [(256, 4096), (1024, 8193)]:
        tensors = [torch.empty(shape, dtype=torch.float32, pin_memory=True) for _ in range(4)]
        ring = [tensor.numpy() for tensor in tensors]
        for values in ring:
            values.fill(.25); values.ravel()[::31] = 3.; values.ravel()[::71] = np.nan
            values.ravel()[1::127] = np.inf; values.ravel()[2::127] = -np.inf
        df = np.full((shape[0], 1), 35000., np.float32)
        mask = np.empty(shape, bool)
        for limit, density in [(7., 'empty'), (2., 'sparse'), (0., 'dense')]:
            critical = np.full((shape[0], 1), limit, np.float64)
            limits = np.broadcast_to(ceil_float32(critical), shape)
            os.environ['TORCHGWAS_HOST_PREDICATE'] = 'numpy'
            expected = select_host_pairs(None, ring[0], df, critical)
            os.environ['TORCHGWAS_HOST_PREDICATE'] = 'native'
            actual = select_host_pairs(None, ring[0], df, critical)
            for got, want in zip(actual, expected):
                if want is None: assert got is None
                else:
                    assert got.shape == want.shape and got.dtype == want.dtype
                    assert got.tobytes() == want.tobytes()
            retained = len(actual[0]); del actual, expected
            cases = [(backend, operation) for backend in ('numpy', 'native') for operation in ('predicate', 'selector')]
            for repeat in range(8):
                order = list(cases); rng.shuffle(order)
                values = ring[repeat % 4]
                for backend, operation in order:
                    os.environ['TORCHGWAS_HOST_PREDICATE'] = backend
                    faults = resource.getrusage(resource.RUSAGE_THREAD)
                    start = time.perf_counter(); cpu = time.thread_time()
                    if operation == 'predicate': fill_predicate_mask(values, limits, mask)
                    else: result = select_host_pairs(None, values, df, critical)
                    cpu = time.thread_time()-cpu; wall = time.perf_counter()-start
                    after = resource.getrusage(resource.RUSAGE_THREAD)
                    records.append(dict(shape=shape, density=density, retained=retained, backend=backend,
                        operation=operation, repeat=repeat, cpu_seconds=cpu, wall_seconds=wall,
                        minor_faults=after.ru_minflt-faults.ru_minflt, major_faults=after.ru_majflt-faults.ru_majflt))
                    if operation == 'selector': del result
            print(json.dumps(dict(shape=shape, density=density, retained=retained,
                medians={backend+'_'+operation:statistics.median(r['cpu_seconds'] for r in records
                    if r['shape']==shape and r['density']==density and r['backend']==backend and r['operation']==operation)
                    for backend,operation in cases})), flush=True)
        del tensors, ring, values, mask
    assert source_identity() == source and context() == native_context
    assert _numpy_core_context() == numpy_core
    (root/'report.json').write_text(json.dumps(dict(source_sha256=source,
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        native_context=native_context, records=records, numpy=np.__version__, torch=torch.__version__,
        settings=settings,
        affinity=sorted(os.sched_getaffinity(0)), concurrent_control=args.concurrent_control,
        scope=__doc__), indent=2)+'\n')


if __name__ == '__main__': main()
