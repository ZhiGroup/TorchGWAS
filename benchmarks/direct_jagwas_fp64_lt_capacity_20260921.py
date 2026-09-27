"""Fixed-size FP64 capacities selected by explicit cuBLASLt instruction flags."""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import time
import numpy as np
import torch
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record
from direct_cublaslt_fp64 import Gemm
from direct_jagwas_fp64_capacity_20260921 import census,N


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--devices',nargs='+',default=['cuda:0','cuda:2']);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False);source=source_identity()
    os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    def load():return subprocess.check_output(['nvidia-smi','--query-gpu=index,name,utilization.gpu,memory.used,memory.total,clocks.sm','--format=csv,noheader,nounits'],text=True)
    before=load();results=[]
    for device in args.devices:
        torch.cuda.set_device(device);properties=torch.cuda.get_device_properties(device)
        if torch.cuda.mem_get_info(device)[0]<2<<30:raise ValueError('Insufficient free memory')
        torch.cuda.set_per_process_memory_fraction((1<<30)/properties.total_memory,device);torch.manual_seed(9219167)
        a=torch.randn((N,N),device=device,dtype=torch.float64);b=torch.randn_like(a);out=torch.empty_like(a)
        workspace=torch.empty(64<<20,device=device,dtype=torch.uint8)
        torch.mm(a[:32,:32],b[:32,:32]);torch.cuda.synchronize(device)
        gemm=Gemm(a,b,out,workspace)
        try:
            indices=np.random.default_rng(9167).integers(0,N,size=(2,32))
            left=a[indices[0]].cpu().numpy();right=b[:,indices[1]].cpu().numpy();expected=np.einsum('ij,ji->i',left,right)
            geometries={};choices={};errors={};algorithms={}
            for kind in ['scalar','tensor']:
                choice=gemm.select(kind)
                for _ in range(3):gemm.run()
                torch.cuda.synchronize(device)
                got=out[indices[0],indices[1]].cpu().numpy();np.testing.assert_allclose(got,expected,rtol=1e-11,atol=1e-9)
                errors[kind]=float(np.max(np.abs(got-expected)))
                geometry=census(gemm.run,root);geometries[kind]=geometry;choices[kind]=choice
                write_record(root/(device.replace(':','_')+'_'+kind+'_geometry.json'),dict(selection=choice,kernels=geometry,correctness_max_abs=errors[kind],durations_recorded=False))
                names=[r['name'] for r in geometry];tensor=any('tensorop_d884gemm' in s or 'tensorop_d1688gemm' in s for s in names)
                scalar=any('dgemm' in s and 'tensorop' not in s for s in names)
                if kind=='scalar' and (tensor or not scalar):raise ValueError('Unexpected scalar kernel census '+str(names))
                if kind=='tensor' and not tensor:raise ValueError('Unexpected Tensor Core census '+str(names))
                algorithms[kind]=(gemm.algo,gemm.work_bytes)
            rows=[];rng=random.Random(9219167);begin=torch.cuda.Event(enable_timing=True);end=torch.cuda.Event(enable_timing=True)
            for repeat in range(5):
                order=['scalar','tensor'];rng.shuffle(order)
                for kind in order:
                    gemm.algo,gemm.work_bytes=algorithms[kind];torch.cuda.synchronize(device)
                    cpu=time.thread_time();wall=time.perf_counter();begin.record()
                    for _ in range(16):gemm.run()
                    end.record();submit=time.thread_time()-cpu;end.synchronize();seconds=begin.elapsed_time(end)/1000/16
                    assert seconds>0
                    row=dict(repeat=repeat,arithmetic=kind,loops=16,dimension=N,useful_flops=2*N**3,gpu_seconds=seconds,
                        flops_per_second=2*N**3/seconds,submit_cpu_seconds=submit/16,wall_seconds=(time.perf_counter()-wall)/16)
                    rows.append(row);write_record(root/(device.replace(':','_')+'_'+kind+'_'+str(repeat)+'.json'),row)
            result=dict(device=device,device_name=properties.name,torch_version=torch.__version__,cuda_version=torch.version.cuda,
                affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),library=gemm.identity,rows=rows,
                kernels=geometries,selections=choices,max_abs_cpu_reference_error=errors,
                resources={field:statistics.median(r['flops_per_second'] for r in rows if r['arithmetic']==kind)
                    for field,kind in [('fp64_flops_per_second','scalar'),('fp64_tensor_flops_per_second','tensor')]},
                calibration_transfer_qualified=False)
            results.append(result);write_record(root/(device.replace(':','_')+'.json'),result)
            print(json.dumps(dict(device=device,resources=result['resources'],errors=errors)),flush=True)
        finally:gemm.close()
        del gemm,a,b,out,workspace,left,right;gc.collect();torch.cuda.empty_cache()
    assert source==source_identity()
    names=['direct_cublaslt_fp64.py','direct_jagwas_fp64_capacity_20260921.py',Path(__file__).name]
    write_record(root/'report.json',dict(records=results,source_sha256=source,
        benchmark_sha256={n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
        nvidia_before=before,nvidia_after=load(),prediction_complete=False,
        documentation='https://docs.nvidia.com/cuda/archive/12.4.1/cublas/index.html#cublasltnumericalimplflags-t',
        scope='Fixed 4096-square, independently constrained FMA/DMMA vendor heuristics, verified by numerical flags, untimed kernel census and sampled CPU dot products. Five paired randomized repeats, median observed rate. No association data, grid interpolation, custom CUDA, peak guarantee or loaded-context qualification.'))


if __name__=='__main__':main()
