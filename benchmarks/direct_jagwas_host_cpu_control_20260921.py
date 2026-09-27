"""Unhooked CPU control for the same fixed operation bank and timer boundary."""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import threading
from concurrent.futures import ThreadPoolExecutor
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
                  OMP_WAIT_POLICY='PASSIVE',GOMP_SPINCOUNT='0')
import torch
from direct_jagwas_host_primitives import build_bank


def worker(device,barrier,library,checkpoint):
    torch.cuda.set_device(device);torch.backends.cuda.matmul.allow_tf32=False
    meter=ctypes.PyDLL(library)
    begin,end=meter.gil_probe_begin_mode,meter.gil_probe_end
    begin.argtypes=[ctypes.c_int];begin.restype=None
    end.argtypes=[ctypes.POINTER(ctypes.c_double)];end.restype=None
    state=(ctypes.c_double*4)();previous=[None]
    def observe(fn):
        begin(0);previous[0]=fn();end(state)
        if any(state[i] for i in [1,2,3]):raise RuntimeError('Unhooked CPU timer reports detached intervals')
        return state[0]
    bank=build_bank('cuda:'+str(device));names=list(bank);random.Random(920119).shuffle(names)
    for _ in range(20):
        for name in names:bank[name]()
    torch.cuda.synchronize();rows=[]
    for repeat in range(5):
        barrier.wait(timeout=120)
        for phase in range(2):
            barrier.wait(timeout=120)
            previous[0]=None
            empty=[observe(lambda:None) for _ in range(100)]
            samples={name:[] for name in names}
            for _ in range(100):
                for name in names:samples[name].append(observe(bank[name]))
            torch.cuda.synchronize();overhead=statistics.fmean(empty)
            for name,values in samples.items():
                raw=statistics.fmean(values)
                rows.append(dict(device=device,repeat=repeat,phase=phase,primitive=name,sample_count=len(values),
                    raw_cpu_sum_seconds=sum(values),raw_cpu_seconds_per_call=raw,empty_cpu_seconds_per_call=overhead,
                    cpu_seconds_per_call=max(0.,raw-overhead)))
            Path(checkpoint).write_text(json.dumps(dict(device=device,rows=rows,complete=False),indent=2))
    result=dict(device=device,rows=rows,complete=True)
    Path(checkpoint).write_text(json.dumps(result,indent=2));return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--library',required=True)
    parser.add_argument('--out',required=True);args=parser.parse_args()
    if os.environ.get('LD_AUDIT'):raise ValueError('This control must run without LD_AUDIT')
    path=Path(args.out)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    library=str(Path(args.library).resolve())
    os.sched_setaffinity(0,range(12,20));torch.set_num_threads(4)
    results=[]
    for mode,devices in [('single0',[0]),('threads',[0,2])]:
        barrier=threading.Barrier(len(devices))
        with ThreadPoolExecutor(max_workers=len(devices)) as pool:
            futures=[pool.submit(worker,device,barrier,library,
                str(path.with_name(path.stem+'.'+mode+'.'+str(device)+'.json'))) for device in devices]
            results.append(dict(mode=mode,workers=[future.result() for future in futures]))
    report=dict(results=results,ld_audit=None,host=os.uname().nodename,torch_version=torch.__version__,
        affinity=sorted(os.sched_getaffinity(0)),primitive_bank='jagwas',fixed_shape=[32,32],
        devices={str(d):dict(name=torch.cuda.get_device_name(d),capability=list(torch.cuda.get_device_capability(d))) for d in [0,2]},
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),
            Path('benchmarks/direct_jagwas_host_primitives.py'),Path('benchmarks/direct_calculator_gil_audit.c'),Path(library)]},
        scope='Same fixed32 bank, seeded mixed order, C thread-CPU timer, empty control, and prior-result destruction as the audited probe, '
            'but without loader hooks. Total CPU only; no GIL partition. Separate processes and changing server/driver load remain explicit confounders.')
    with path.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(dict(contexts=[row['mode'] for row in results],ld_audit=None)),flush=True)


if __name__=='__main__':main()
