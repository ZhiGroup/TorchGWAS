"""Record caller-thread CPU inside CUDA copy/sync APIs on fixed primitives."""
import argparse
import ctypes
from datetime import datetime,timezone
import gc
import json
import os
from pathlib import Path
import random
import statistics
import time

import numpy as np
import torch
from direct_device_significance_primitives import build_bank
from direct_device_significance_primitive_measure_20260921 import bulk_bank, telemetry
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record

class Event(ctypes.Structure):
    _fields_ = [(name,ctypes.c_uint64) for name in ['kind','bytes','direction','begin_wall','end_wall','begin_cpu','end_cpu','result']]
class Record(ctypes.Structure):
    _fields_ = [(name,ctypes.c_uint64) for name in ['begin_wall','end_wall','begin_cpu','end_cpu','return_wall','return_cpu','count','overflow']]+[('events',Event*64)]

parser=argparse.ArgumentParser()
parser.add_argument('--out',required=True)
parser.add_argument('--library',required=True)
parser.add_argument('--device',default='cuda:0')
args=parser.parse_args()
root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
source=source_identity()
library=ctypes.CDLL(str(Path(args.library).resolve()))
library.selection_probe_begin.argtypes=[];library.selection_probe_begin.restype=None
library.selection_probe_begin_uninstrumented.argtypes=[];library.selection_probe_begin_uninstrumented.restype=None
library.selection_probe_return.argtypes=[];library.selection_probe_return.restype=None
library.selection_probe_end.argtypes=[ctypes.POINTER(Record)];library.selection_probe_end.restype=None
library.selection_probe_record_bytes.argtypes=[];library.selection_probe_record_bytes.restype=ctypes.c_uint64
assert library.selection_probe_record_bytes()==ctypes.sizeof(Record)
torch.set_num_threads(2);torch.set_num_interop_threads(1)
torch.cuda.set_device(args.device)
limit=512<<20
prop=torch.cuda.get_device_properties(args.device)
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,args.device)
if torch.cuda.mem_get_info(args.device)[0]<limit:raise ValueError('Insufficient free GPU memory')
banks={'tiny':build_bank(args.device),'bulk':bulk_bank(args.device)}
banks['tiny']={name:function for name,function in banks['tiny'].items() if name.startswith(('copy_','nonzero_','loop_','ready_'))}
for scale,bank in banks.items():
    # Independent uint8 QC vectors, with all successful native status values.
    # Invalid input exits instead of contributing to a sustained scan price.
    status=np.arange(32 if scale=='tiny' else 1<<20,dtype=np.uint8)%3
    bank['qc_malformed_empty']=lambda status=status:np.flatnonzero(status==3)
    bank['qc_count_status']=lambda status=status:int((status==1).sum())
for bank in banks.values():
    for function in bank.values():
        for _ in range(8):function()
        torch.cuda.synchronize(args.device)
plan=[];rng=random.Random(9211842)
for repeat in range(9):
    order=[(scale,name,mode) for scale,bank in banks.items() for name in bank for mode in ['plain','record']]
    rng.shuffle(order);plan.append(dict(repeat=repeat,order=order))
paths=[Path(__file__),Path(__file__).with_name('device_selection_runtime_probe.cpp'),
       Path(__file__).with_name('direct_device_significance_primitives.py'),
       Path(__file__).with_name('direct_device_significance_primitive_measure_20260921.py')]
for path in paths:(root/path.name).write_bytes(path.read_bytes())
hashes={p.name:sha256_file(p) for p in paths}
write_record(root/'protocol.json',dict(source_sha256=source,harness_sha256=hashes,library_sha256=sha256_file(args.library),
    device=dict(index=args.device,name=prop.name,uuid=str(prop.uuid)),torch_version=torch.__version__,cuda_runtime=torch.version.cuda,
    affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),torch_interop_threads=torch.get_num_interop_threads(),
    environment={k:v for k,v in sorted(os.environ.items()) if k.startswith(('CUDA_','TORCH','OMP_','MKL_','OPENBLAS_','GOMP_','LD_'))},
    memory_limit_bytes=limit,plan=plan,calls_per_batch=8,
    boundary='One Python API invocation with a separate return checkpoint before output destruction; release then ends at the final marker. CUDA runtime events record caller thread CPU and monotonic wall directly. External GPU drain lies outside every marked invocation. Both modes use the same return checkpoint.',
    scope='Independent fixed primitive decomposition. The host-only preload shim forwards all CUDA calls unchanged and records copy/synchronization intervals. It launches no kernels. Both plain and instrumented invocations use identical boundary markers; plain disables internal runtime clocks. All repetitions are retained.'))
rows=[];states=[];record=Record();gc.collect();gc.disable()
try:
    for round_ in plan:
        states.append(dict(repeat=round_['repeat'],phase='before',telemetry=telemetry()))
        round_rows=[]
        for scale,name,mode in round_['order']:
            function=banks[scale][name]
            traces=[]
            torch.cuda.synchronize(args.device)
            batch_wall=time.perf_counter_ns();batch_cpu=time.thread_time_ns()
            begin = library.selection_probe_begin if mode=='record' else library.selection_probe_begin_uninstrumented
            for sample in range(8):
                begin()
                result=function()
                library.selection_probe_return()
                del result
                library.selection_probe_end(ctypes.byref(record))
                assert not record.overflow and record.count<=64
                events=[{field:getattr(event,field) for field,_ in Event._fields_} for event in record.events[:record.count]]
                assert all(event['result']==0 for event in events)
                if mode=='plain': assert not events
                elif name.startswith('qc_'):assert not events
                elif name.startswith(('copy_','nonzero_')):
                    # Verified production CUDA primitive contract: one D2H
                    # runtime submission followed by exactly one stream wait.
                    assert [event['kind'] for event in events]==[1,2],(name,events)
                    expected_bytes=(4 if name.startswith('nonzero_') else
                        (32 if scale=='tiny' else 1<<20)*{'fp32':4,'int64':8,'uint8':1}[name.rsplit('_',1)[1]])
                    assert events[0]['bytes']==expected_bytes and events[0]['direction']==2
                before={field:getattr(record,field) for field,_ in Record._fields_ if field!='events'}
                traces.append(dict(sample=sample,**before,events=events))
            batch_cpu=(time.thread_time_ns()-batch_cpu)*1e-9/8
            batch_wall=(time.perf_counter_ns()-batch_wall)*1e-9/8
            torch.cuda.synchronize(args.device)
            row=dict(repeat=round_['repeat'],scale=scale,primitive=name,mode=mode,
                instrumented_batch_wall_seconds=batch_wall,instrumented_batch_thread_cpu_seconds=batch_cpu,traces=traces)
            rows.append(row);round_rows.append(row)
        write_record(root/f"repeat_{round_['repeat']}.json",dict(observations=round_rows))
        states.append(dict(repeat=round_['repeat'],phase='after',telemetry=telemetry()))
        print(json.dumps(dict(repeat=round_['repeat'],batches=len(round_rows),traces=sum(len(r['traces']) for r in round_rows))),flush=True)
finally:
    gc.enable();write_record(root/'telemetry.json',states)
assert source==source_identity() and hashes=={p.name:sha256_file(p) for p in paths}
write_record(root/'report.json',dict(protocol_sha256=sha256_file(root/'protocol.json'),observations=rows,
    instrumentation_qualified=False,
    limitations=['ctypes marker and Python bookkeeping costs are retained',
        'Both batch envelopes include trace validation/serialization; compare per-invocation marker intervals, not batch envelopes',
        'No per-phase GIL ownership attribution', 'Pageable cudaMemcpyAsync may block inside its runtime call',
        'Other shared-host workloads and instrumentation can affect timing'],
    scope='Measured caller-thread CPU within copy and synchronization runtime calls, not a fitted cost model or complete kernel service bank.'))
print(json.dumps(dict(complete=True,batches=len(rows),marked_invocations=sum(len(r['traces']) for r in rows),runtime_traced_invocations=sum(len(r['traces']) for r in rows if r['mode']=='record'))),flush=True)
