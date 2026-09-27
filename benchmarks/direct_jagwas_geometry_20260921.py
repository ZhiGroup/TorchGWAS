"""Duration-free reduction census on bounded synthetic tensors.

No GWAS timings or timing coefficients enter the exported record. Profiler
trace durations are discarded, and temporary raw traces are removed.
"""
import argparse
import gc
import json
from pathlib import Path
import tempfile
import numpy as np
import torch
from torch.utils._python_dispatch import TorchDispatchMode
from torch.utils._pytree import tree_flatten
from torchgwas.jagwas_projection import JagwasReduction
from torchgwas.reduction_tensor_work import jagwas_tensor_work,jagwas_factor_memory_floor
from torchgwas.geometry_collection import kernel_census,write_record
from torchgwas.detailed_calibration import source_identity

parser=argparse.ArgumentParser()
parser.add_argument('--out',default='results/jagwas_geometry_20260921')
parser.add_argument('--device',default='cuda:0')
args=parser.parse_args()
ROOT=Path(args.out)
ROOT.mkdir(parents=True,exist_ok=False)
DEVICE=args.device
LIMIT=2<<30
SHAPES=[(257,13,7),(2049,128,512),(4097,512,2048)]
torch.cuda.set_device(DEVICE)
torch.cuda.set_per_process_memory_fraction(LIMIT/torch.cuda.get_device_properties(DEVICE).total_memory,DEVICE)
torch.set_num_threads(2);torch.set_num_interop_threads(1)
torch.backends.cuda.matmul.allow_tf32=False
source=source_identity()


def describe(value):
    return dict(shape=list(value.shape),dtype=str(value.dtype),device_type=value.device.type,
                bytes=value.numel()*value.element_size(),stride=list(value.stride()))


class Observe(TorchDispatchMode):
    def __init__(self):self.steps=[]
    def __torch_dispatch__(self,func,types,args=(),kwargs=None):
        kwargs=kwargs or {}
        inputs=[describe(v) for v in tree_flatten((args,kwargs))[0] if isinstance(v,torch.Tensor)]
        result=func(*args,**kwargs)
        outputs=[describe(v) for v in tree_flatten(result)[0] if isinstance(v,torch.Tensor)]
        self.steps.append(dict(op=str(func),inputs=inputs,outputs=outputs))
        return result


rows=[]
for n,b,k in SHAPES:
    for precision in ['float32','float64']:
        dtype=getattr(torch,precision)
        ledgers={phase:jagwas_tensor_work(n,b,k,phase=phase,compute_dtype=precision) for phase in ['prepare','reduce']}
        floor=jagwas_factor_memory_floor(n,k,compute_dtype=precision)
        explicit=sum(v['distinct_temporary_bytes'] for v in ledgers.values())+sum(v['bytes'] for v in ledgers['reduce']['initial_storages'])+floor['arrays']['phenotype']
        if explicit>LIMIT//2:raise ValueError('Synthetic census exceeds explicit probe storage budget')
        free=torch.cuda.mem_get_info(DEVICE)[0]
        if free<LIMIT:raise ValueError('Insufficient free memory for bounded census')
        # Avoid real phenotype data; the synthetic matrix is full rank with a
        # broad spectral margin. Host upload is part of real factor preparation.
        rng=np.random.default_rng(9219122)
        y=rng.normal(size=(n,k)).astype(precision)
        beta=torch.empty((b,k),dtype=dtype,device=DEVICE)
        stat=torch.randn((b,k),dtype=dtype,device=DEVICE)
        status=torch.zeros(b,dtype=torch.uint8,device=DEVICE)
        df=torch.full((b,),float(n-2),dtype=torch.float32,device=DEVICE)
        factor=JagwasReduction().prepare(y,device=DEVICE)
        for phase in ['prepare','reduce']:
            def run():
                return (JagwasReduction().prepare(y,device=DEVICE) if phase=='prepare'
                        else factor.reduce(beta,stat,status,df,1))
            value=run();torch.cuda.synchronize(DEVICE);del value
            observer=Observe()
            with observer:value=run()
            torch.cuda.synchronize(DEVICE);del value
            gc.collect();torch.cuda.synchronize(DEVICE)
            baseline=torch.cuda.memory_allocated(DEVICE)
            torch.cuda.reset_peak_memory_stats(DEVICE)
            with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,
                    torch.profiler.ProfilerActivity.CUDA]) as profiler:
                value=run();torch.cuda.synchronize(DEVICE)
            peak=torch.cuda.max_memory_allocated(DEVICE)
            del value
            with tempfile.TemporaryDirectory(prefix='joint-geometry-') as temporary:
                trace=Path(temporary)/'trace.json';profiler.export_chrome_trace(str(trace))
                events=json.loads(trace.read_text())['traceEvents']
                kernels=kernel_census(events)
                runtimes=[e['name'] for e in events if e.get('cat')=='cuda_runtime']
                # Retain names/counts only, never timestamps or durations.
                transfers=[e['name'] for e in events if e.get('cat') in ('gpu_memcpy','gpu_memset')]
            record=dict(N=n,B=b,K=k,compute_dtype=precision,phase=phase,
                tensor_work=ledgers[phase],factor_capacity_floor=floor,
                observed_tensor_steps=observer.steps,kernels=kernels,
                cuda_runtime_calls=runtimes,transfer_events=transfers,
                observed_extra_allocated_bytes=peak-baseline,
                source_explicit_probe_bytes=explicit,durations_recorded=False)
            rows.append(record)
            write_record(ROOT/f'{n}_{b}_{k}_{precision}_{phase}.json',record)
            print(json.dumps(dict(N=n,B=b,K=k,compute_dtype=precision,phase=phase,
                tensor_dispatches=len(observer.steps),kernels=len(kernels),extra_allocated_bytes=peak-baseline)),flush=True)
        del factor,y,beta,stat,status,df,run,profiler,observer,record
        gc.collect();torch.cuda.empty_cache()
if source!=source_identity():raise ValueError('Source changed during census')
properties=torch.cuda.get_device_properties(DEVICE)
write_record(ROOT/'census.json',dict(source_sha256=source,torch_version=torch.__version__,
    device=properties.name,compute_capability=[properties.major,properties.minor],
    memory_probe_limit_bytes=LIMIT,allow_tf32=False,durations_recorded=False,rows=rows,
    scope='Synthetic source-operation, CUDA launch, API synchronization and allocation census. No timing coefficients. Observed allocator peaks do not establish an admission bound.'))
