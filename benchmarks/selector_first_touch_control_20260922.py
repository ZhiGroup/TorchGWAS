"""Fixed 64-MiB first/warm copy controls followed by held-out selectors.

No fitted selector/GWAS timings and no calibration publication. All raw signed
first-minus-warm observations remain available; no sample selection or clipping.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import statistics
import time
from datetime import datetime, timezone

from direct_bounded_host_selection_prices_20260921 import allocation_counters, measured_call
import direct_bounded_host_selection_prices_20260921 as meter_module
import numpy as np
import torch
from torchgwas.detailed_calibration import source_identity, _numpy_core_context
from torchgwas.host_significance import select_host_pairs
from torchgwas.native_host_predicate import context as native_context


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4)
    assert not np._core.multiarray._get_madvise_hugepage()
    source=source_identity();core=_numpy_core_context();native=native_context()
    harness={str(Path(p).name):hashlib.sha256(Path(p).read_bytes()).hexdigest()
        for p in [__file__,meter_module.__file__]}
    report=dict(observation_started_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=source,harness_sha256=harness,context=dict(numpy_core=core,native=native,
        numpy=np.__version__,torch=torch.__version__,python=platform.python_version(),libc=list(platform.libc_ver()),
        page_bytes=os.sysconf('SC_PAGE_SIZE'),numpy_madvise_hugepage=False,torch_threads=torch.get_num_threads(),
        affinity=sorted(os.sched_getaffinity(0)),environment={key:os.getenv(key) for key in
        ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_WAIT_POLICY','GOMP_SPINCOUNT',
         'NUMPY_MADVISE_HUGEPAGE','MALLOC_MMAP_THRESHOLD_','MALLOC_TRIM_THRESHOLD_','GLIBC_TUNABLES']}),
        controls=[],selectors=[],prediction_complete=False,published=False,scope=__doc__)
    def save():
        (root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    size=64<<20;pages=size//report['context']['page_bytes']
    resident=np.full(size,71,np.uint8)
    # Fixed workload independent of selected shapes and survivor counts.
    # First must precede warm on the same allocation. Record this ordered pair
    # in every repeat; do not claim that order/frequency effects are eliminated.
    for repeat in range(7):
        destination,allocation=measured_call(lambda:np.empty(size,np.uint8))
        _,first=measured_call(lambda:np.copyto(destination,resident))
        _,warm=measured_call(lambda:np.copyto(destination,resident))
        assert np.array_equal(destination,resident)
        before=allocation_counters();faults=resource.getrusage(resource.RUSAGE_THREAD)
        started=time.perf_counter();cpu=time.thread_time()
        del destination
        release=dict(cpu_seconds=time.thread_time()-cpu,wall_seconds=time.perf_counter()-started)
        after=resource.getrusage(resource.RUSAGE_THREAD)
        release.update(minor_faults=after.ru_minflt-faults.ru_minflt,
            major_faults=after.ru_majflt-faults.ru_majflt,allocator_before=before,allocator_after=allocation_counters())
        report['controls'].append(dict(repeat=repeat,bytes=size,pages=pages,allocation=allocation,
            first=first,warm=warm,release=release,
            signed_first_less_warm_cpu=first['cpu_seconds']-warm['cpu_seconds'],
            signed_first_less_warm_faults=first['minor_faults']-warm['minor_faults']))
        save()
    report['control_summary']=dict(
        first_touch_cpu_seconds_per_page=statistics.median(r['signed_first_less_warm_cpu']/pages for r in report['controls']),
        fresh_minor_faults=[r['first']['minor_faults'] for r in report['controls']],
        warm_minor_faults=[r['warm']['minor_faults'] for r in report['controls']],
        scope='Diagnostic signed paired median at one fixed extent. Not added to allocation-inclusive primitive prices.')
    print(json.dumps(report['control_summary']),flush=True)
    del resident
    rng=random.Random(9221347)
    for shape in [(256,4096),(1024,8193)]:
        tensors=[torch.empty(shape,dtype=torch.float32,pin_memory=True) for _ in range(4)]
        ring=[tensor.numpy() for tensor in tensors]
        df=np.full((shape[0],1),35000.,np.float32)
        critical=np.full((shape[0],1),2.,np.float64)
        for density in ['sparse','dense']:
            for values in ring:
                values.fill(3. if density=='dense' else .25)
                if density=='sparse':values.ravel()[::31]=3.
            expected=None
            for backend in ['numpy','native']:
                os.environ['TORCHGWAS_HOST_PREDICATE']=backend
                result=select_host_pairs(None,ring[0],df,critical)
                fingerprint=[None if x is None else (x.shape,x.dtype.str,hashlib.sha256(memoryview(x)).hexdigest()) for x in result]
                if expected is None:expected=fingerprint
                else:assert fingerprint==expected
                retained=len(result[0]);del result
            for repeat in range(7):
                order=['numpy','native'];rng.shuffle(order)
                for backend in order:
                    os.environ['TORCHGWAS_HOST_PREDICATE']=backend
                    result,observation=measured_call(lambda:select_host_pairs(None,ring[repeat%4],df,critical))
                    assert len(result[0])==retained
                    before=allocation_counters();cpu=time.thread_time();started=time.perf_counter()
                    del result
                    release=dict(cpu_seconds=time.thread_time()-cpu,wall_seconds=time.perf_counter()-started,
                        allocator_before=before,allocator_after=allocation_counters())
                    report['selectors'].append(dict(shape=shape,density=density,backend=backend,repeat=repeat,
                        retained=retained,selected_bytes=24*retained,observation=observation,release=release))
            save()
            print(json.dumps(dict(shape=shape,density=density,retained=retained)),flush=True)
        del tensors,ring,values,df,critical
    assert source==source_identity() and core==_numpy_core_context() and native==native_context()
    assert harness=={str(Path(p).name):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [__file__,meter_module.__file__]}
    report['observation_finished_at_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':main()
