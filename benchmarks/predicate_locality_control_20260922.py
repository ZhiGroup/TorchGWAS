"""Fixed predicate locality controls, separating reusable from rotating buffers."""
import argparse
from datetime import datetime,timezone
import ctypes
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
from torchgwas.host_significance import fill_predicate_mask
from torchgwas.native_host_predicate import context
from torchgwas.detailed_calibration import source_identity,_numpy_core_context


def mapping(array):
    pointer=array.__array_interface__['data'][0]
    line=next(line for line in Path('/proc/self/maps').read_text().splitlines()
        if int(line.split()[0].split('-')[0],16)<=pointer<int(line.split()[0].split('-')[1],16))
    base=line.split('-')[0]
    numa=next((line for line in Path('/proc/self/numa_maps').read_text().splitlines()
        if line.split()[0]==base),None)
    return dict(bytes=array.nbytes,mapping=line,numa_mapping=numa,
        scope='Containing mapping only; it can include more than this array.')


def main(args):
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(4);torch.set_num_interop_threads(1)
    assert not np._core.multiarray._get_madvise_hugepage()
    source=source_identity();core=_numpy_core_context();native=context()
    started=datetime.now(timezone.utc).isoformat()
    b,k,slots=256,4096,16
    tensor=torch.empty((slots,b,k),dtype=torch.float32,pin_memory=True)
    arrays=dict(pageable=np.empty((slots,b,k),np.float32),pinned=tensor.numpy())
    masks={name:np.empty((slots,b,k),bool) for name in arrays}
    for values in arrays.values():
        values.fill(.25);values.ravel()[::31]=3.;values.ravel()[::71]=np.nan
    critical=np.broadcast_to(np.full((b,1),2.,np.float32),(b,k))
    layouts={name:mapping(values) for name,values in arrays.items()}
    cpu=ctypes.CDLL(None).sched_getcpu;cpu.restype=ctypes.c_int;cpu.argtypes=[]
    cache={}
    for logical_cpu in sorted(os.sched_getaffinity(0)):
        base=Path('/sys/devices/system/cpu')/f'cpu{logical_cpu}'/'cache'
        cache[str(logical_cpu)]=[{key:(entry/key).read_text().strip() for key in
            ['level','type','size','shared_cpu_list']} for entry in sorted(base.glob('index*'))]
    cases=[(backend,allocation,pattern) for backend in ['numpy','native']
        for allocation in arrays for pattern in ['reused','rotating']]
    rng=random.Random(922255);records=[]
    # First touch every output and verify identical masks before timing.
    for name,values in arrays.items():
        os.environ['TORCHGWAS_HOST_PREDICATE']='numpy'
        for slot in range(slots):fill_predicate_mask(values[slot],critical,masks[name][slot])
    mask_hashes={name:hashlib.sha256(memoryview(mask).cast('B')).hexdigest() for name,mask in masks.items()}
    for repeat in range(24):
        order=list(cases);rng.shuffle(order)
        for backend,allocation,pattern in order:
            os.environ['TORCHGWAS_HOST_PREDICATE']=backend
            slot=0 if pattern=='reused' else 1+repeat%(slots-1)
            values=arrays[allocation][slot];out=masks[allocation][slot]
            if pattern=='reused':
                for _ in range(2):fill_predicate_mask(values,critical,out)
            before=resource.getrusage(resource.RUSAGE_THREAD);first_cpu=cpu()
            wall=time.perf_counter();thread=time.thread_time()
            fill_predicate_mask(values,critical,out)
            thread=time.thread_time()-thread;wall=time.perf_counter()-wall
            last_cpu=cpu();after=resource.getrusage(resource.RUSAGE_THREAD)
            records.append(dict(backend=backend,allocation=allocation,pattern=pattern,slot=slot,
                repeat=repeat,cpu_seconds=thread,wall_seconds=wall,first_cpu=first_cpu,last_cpu=last_cpu,
                user_seconds=after.ru_utime-before.ru_utime,system_seconds=after.ru_stime-before.ru_stime,
                minor_faults=after.ru_minflt-before.ru_minflt,major_faults=after.ru_majflt-before.ru_majflt))
    summaries=[]
    for backend,allocation,pattern in cases:
        rows=[r for r in records if (r['backend'],r['allocation'],r['pattern'])==(backend,allocation,pattern)]
        item=dict(backend=backend,allocation=allocation,pattern=pattern,
            median_cpu_seconds=statistics.median(r['cpu_seconds'] for r in rows),
            min_cpu_seconds=min(r['cpu_seconds'] for r in rows),max_cpu_seconds=max(r['cpu_seconds'] for r in rows),
            cpu_changes_during_call=sum(r['first_cpu']!=r['last_cpu'] for r in rows),
            minor_faults=sum(r['minor_faults'] for r in rows),major_faults=sum(r['major_faults'] for r in rows))
        summaries.append(item);print(json.dumps(item),flush=True)
    assert source==source_identity() and core==_numpy_core_context() and native==context()
    assert mask_hashes=={name:hashlib.sha256(memoryview(mask).cast('B')).hexdigest() for name,mask in masks.items()}
    report=dict(summaries=summaries,observations=records,source_sha256=source,numpy_core=core,
        native_context=native,input_mapping=layouts,cache_geometry=cache,mask_hashes=mask_hashes,
        shape=[b,k],slots=slots,settings=dict(affinity=sorted(os.sched_getaffinity(0)),
            torch_threads=torch.get_num_threads(),torch_interop_threads=torch.get_num_interop_threads(),
            numpy_madvise_hugepage=False),
        observation_started_at_utc=started,observation_finished_at_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Randomized fixed 1M-cell predicates. Reused buffers receive two untimed warm calls; rotating buffers span 64MiB input and 16MiB masks per allocator, without forced eviction or a claim of cold cache. No GWAS timing fit or calibrated profile publication.')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    main(parser.parse_args())
