"""Independent fixed N65/K32 factor phases, no association or shape timing grid."""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import tempfile
import time
import numpy as np
import torch
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record,kernel_census
from torchgwas.jagwas_preparation import FACTOR_REFERENCE,factor_phases
from torchgwas.reduction_tensor_work import jagwas_cutoff_method,jagwas_factor_arithmetic_work

METHOD=jagwas_cutoff_method()
FACTOR_PHASES=factor_phases(METHOD)


def build_bank(device):
    n,k=FACTOR_REFERENCE
    y=np.random.default_rng(9219164).normal(size=(n,k)).astype(np.float32)
    matrix=torch.as_tensor(y).to(device)
    # The source's phases at the reference (one sample block: gram_rows(65, 32) = 65):
    # FP32 block cast to FP64, FP64 Gram accumulate and normalise, then per method.
    block=matrix.double()
    double=torch.zeros((k,k),dtype=torch.float64,device=device).addmm_(block.T,block).div_(float(n))
    common=dict(upload=lambda:torch.as_tensor(y).to(device),
        correlation=lambda:torch.zeros((k,k),dtype=torch.float64,device=device).addmm_(block.T,block).div_(float(n)),
        cast=lambda:matrix.double())
    if METHOD=='eigen':
        # eigh, the one read of the K-value spectrum, scaling of the kept
        # directions (all K at the reference) and QR (R only).
        values,vectors=torch.linalg.eigh(double)
        scaled=(vectors/values.sqrt()).T
        return dict(common,eigh=lambda:torch.linalg.eigh(double),
            spectrum=lambda:values.cpu().numpy(),
            scale=lambda:(vectors/values.sqrt()).T,
            qr=lambda:torch.linalg.qr(scaled,mode='r'))
    # Rounding cutoff: cholesky_ex (no info sync), solve, then the rank
    # check's one read of (info, ||L^-1||_F).
    factor,info=torch.linalg.cholesky_ex(double)
    identity=torch.eye(k,dtype=torch.float64,device=device)
    inverse=torch.linalg.solve_triangular(factor,identity,upper=False)
    return dict(common,cholesky=lambda:torch.linalg.cholesky_ex(double),
        identity=lambda:torch.eye(k,dtype=torch.float64,device=device),
        solve=lambda:torch.linalg.solve_triangular(factor,identity,upper=False),
        rank_check=lambda:torch.stack((info.double(),torch.linalg.vector_norm(inverse))).tolist())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--devices',nargs='+',default=['cuda:0','cuda:2']);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False);source=source_identity()
    os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    rows=[];banks={}
    for device in args.devices:
        torch.cuda.set_device(device);properties=torch.cuda.get_device_properties(device)
        torch.cuda.set_per_process_memory_fraction((256<<20)/properties.total_memory,device)
        if torch.cuda.mem_get_info(device)[0]<512<<20:raise ValueError('Insufficient free memory for bounded fixed reference')
        bank=build_bank(device)
        for fn in bank.values():
            for _ in range(8):value=fn();torch.cuda.synchronize(device)
        del value
        geometry={}
        for phase,fn in bank.items():
            with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,torch.profiler.ProfilerActivity.CUDA]) as profiler:
                value=fn();torch.cuda.synchronize(device)
            with tempfile.TemporaryDirectory(dir=root,prefix='phase-') as temporary:
                trace=Path(temporary)/'trace.json';profiler.export_chrome_trace(str(trace))
                events=json.loads(trace.read_text())['traceEvents']
                # upload (H2D) and the eigen factor's spectrum read (D2H) are transfers only.
                transfer_only=phase in ('upload','spectrum')
                kernels=([ ] if transfer_only else kernel_census(events))
                if transfer_only and any(e.get('cat')=='kernel' for e in events):
                    raise ValueError('Unexpected compute kernel in transfer-only primitive: '+phase)
                geometry[phase]=dict(kernels=kernels,
                    cuda_runtime=[e['name'] for e in events if e.get('cat')=='cuda_runtime'],
                    transfers=[e['name'] for e in events if e.get('cat') in ('gpu_memcpy','gpu_memset')])
            del value,profiler
        records=[];rng=random.Random(9219164)
        for repeat in range(5):
            order=list(FACTOR_PHASES);rng.shuffle(order)
            for phase in order:
                fn=bank[phase];value=None;torch.cuda.synchronize(device)
                start=time.perf_counter();cpu=time.thread_time()
                for _ in range(32):value=fn();torch.cuda.synchronize(device)
                cpu=(time.thread_time()-cpu)/32;wall=(time.perf_counter()-start)/32
                del value
                records.append(dict(repeat=repeat,phase=phase,cpu_seconds=cpu,wall_seconds=wall,
                    non_cpu_seconds=wall-cpu,calls=32))
            write_record(root/(device.replace(':','_')+'_repeat_'+str(repeat)+'.json'),dict(records=records[-len(FACTOR_PHASES):]))
        context=dict(device=device,device_name=properties.name,torch_version=torch.__version__,cuda_version=torch.version.cuda,
            affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),allow_tf32=False,
            boundary='single_worker_call_then_device_synchronize')
        phases={name:{key:statistics.median(r[key] for r in records if r['phase']==name)
                      for key in ['cpu_seconds','non_cpu_seconds']} for name in FACTOR_PHASES}
        valid=all(value>=0 for p in phases.values() for value in p.values())
        record=dict(device=device,context=context,records=records,geometry=geometry,prices_nonnegative=valid,
            source_sha256=jagwas_factor_arithmetic_work(*FACTOR_REFERENCE)['source_sha256'],reference_shape=FACTOR_REFERENCE,compute_dtype='float32',method=METHOD,
            boundary='call_then_device_synchronize',phases=phases,transfer_qualified=False,prediction_complete=False)
        rows.append(record);banks[device]=record
        write_record(root/(device.replace(':','_')+'.json'),record)
        print(json.dumps(dict(device=device,prices_nonnegative=valid,phases=phases)),flush=True)
        del bank,fn;gc.collect();torch.cuda.empty_cache()
    if source!=source_identity():raise ValueError('Source changed during fixed primitive measurement')
    write_record(root/'report.json',dict(banks=banks,source_sha256=source,
        benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Fixed independent N65/K32 phase services with call-plus-synchronize boundary, randomized phase order, five complete repetitions. Separate serial worker contexts, not concurrent-factor or runtime-rank validation.',
        prediction_complete=False,transfer_qualified=False))


if __name__=='__main__':main()
