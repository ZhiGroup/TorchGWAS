"""Independent fixed-square FP64 cuBLAS GemmEx capacities with observed kernel families.

Uses vendor DGEMM with an owned handle; no custom CUDA or association kernel.
Default and pedantic modes are independently checked by untimed kernel census.
"""
import argparse
import ctypes as C
import gc
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import time
import numpy as np
import torch
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import kernel_census,write_record

DOC='https://docs.nvidia.com/cuda/archive/12.4.1/cublas/index.html#tensor-core-usage'
N=4096


def check(status):
    if status:raise RuntimeError('cuBLAS status '+str(status))


def library():
    paths={line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
           if '/libcublas.so.' in line and line.split()[-1].startswith('/')}
    if len(paths)!=1:raise ValueError('Require the single cuBLAS library loaded by Torch')
    path=Path(paths.pop());lib=C.CDLL(str(path))
    signatures={'cublasCreate_v2':[C.POINTER(C.c_void_p)],'cublasDestroy_v2':[C.c_void_p],
        'cublasSetStream_v2':[C.c_void_p,C.c_void_p],
        'cublasSetMathMode':[C.c_void_p,C.c_int],'cublasGetMathMode':[C.c_void_p,C.POINTER(C.c_int)],
        'cublasGetVersion_v2':[C.c_void_p,C.POINTER(C.c_int)],
        'cublasGemmEx':[C.c_void_p,C.c_int,C.c_int,C.c_int,C.c_int,C.c_int,C.c_void_p,
            C.c_void_p,C.c_int,C.c_int,C.c_void_p,C.c_int,C.c_int,C.c_void_p,C.c_void_p,C.c_int,C.c_int,C.c_int,C.c_int],
        'cublasDgemm_v2':[C.c_void_p,C.c_int,C.c_int,C.c_int,C.c_int,C.c_int,C.POINTER(C.c_double),
            C.c_void_p,C.c_int,C.c_void_p,C.c_int,C.POINTER(C.c_double),C.c_void_p,C.c_int]}
    for name,args in signatures.items():getattr(lib,name).argtypes=args;getattr(lib,name).restype=C.c_int
    return lib,dict(path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def census(fn,root):
    with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,torch.profiler.ProfilerActivity.CUDA]) as profiler:
        fn();torch.cuda.synchronize()
    with tempfile.TemporaryDirectory(dir=root,prefix='fp64-trace-') as temp:
        path=Path(temp)/'trace.json';profiler.export_chrome_trace(str(path))
        kernels=kernel_census(json.loads(path.read_text())['traceEvents'])
    return kernels


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--devices',nargs='+',default=['cuda:0','cuda:2']);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False;source=source_identity();records=[]
    def load():
        return subprocess.check_output(['nvidia-smi','--query-gpu=index,name,utilization.gpu,memory.used,memory.total,clocks.sm',
            '--format=csv,noheader,nounits'],text=True)
    before=load()
    for device in args.devices:
        torch.cuda.set_device(device);properties=torch.cuda.get_device_properties(device)
        if torch.cuda.mem_get_info(device)[0]<2<<30:raise ValueError('Insufficient free memory for bounded FP64 primitive')
        torch.cuda.set_per_process_memory_fraction((1<<30)/properties.total_memory,device)
        torch.manual_seed(9219167)
        a=torch.randn((N,N),device=device,dtype=torch.float64);b=torch.randn_like(a);out=torch.empty_like(a)
        # Initialize Torch's vendor library before binding its exact loaded path.
        torch.mm(a[:32,:32],b[:32,:32]);torch.cuda.synchronize(device)
        lib,identity=library();handle=C.c_void_p();check(lib.cublasCreate_v2(C.byref(handle)))
        try:
            stream=torch.cuda.current_stream(device);check(lib.cublasSetStream_v2(handle,C.c_void_p(stream.cuda_stream)))
            version=C.c_int();check(lib.cublasGetVersion_v2(handle,C.byref(version)))
            alpha,beta=C.c_double(1.),C.c_double(0.)
            # Row-major C=A@B corresponds to column-major C.T=B.T@A.T.
            def run():
                check(lib.cublasGemmEx(handle,0,0,N,N,N,C.byref(alpha),C.c_void_p(b.data_ptr()),1,N,
                    C.c_void_p(a.data_ptr()),1,N,C.byref(beta),C.c_void_p(out.data_ptr()),1,N,
                    71 if mode==2 else 70,-1))
            indices=np.random.default_rng(9167).integers(0,N,size=(2,32))
            left=a[indices[0]].cpu().numpy();right=b[:,indices[1]].cpu().numpy()
            expected=np.einsum('ij,ji->i',left,right)
            geometries={};errors={}
            for name,mode in [('scalar',2),('tensor',0)]:
                check(lib.cublasSetMathMode(handle,mode));actual=C.c_int()
                check(lib.cublasGetMathMode(handle,C.byref(actual)));assert actual.value==mode
                for _ in range(3):run()
                torch.cuda.synchronize(device)
                got=out[indices[0],indices[1]].cpu().numpy()
                np.testing.assert_allclose(got,expected,rtol=1e-11,atol=1e-9)
                errors[name]=float(np.max(np.abs(got-expected)))
                geometry=census(run,root);names=[k['name'] for k in geometry]
                tensor=any('tensorop_d884gemm' in n or 'tensorop_d1688gemm' in n for n in names)
                scalar=any('dgemm' in n and 'tensorop' not in n for n in names)
                if name=='scalar' and (tensor or not scalar):raise ValueError('Pedantic census did not identify scalar DGEMM: '+str(names))
                if name=='tensor' and not tensor:raise ValueError('Default census did not identify FP64 Tensor Core GEMM: '+str(names))
                geometries[name]=geometry
            rows=[];rng=random.Random(9219167);begin=torch.cuda.Event(enable_timing=True);end=torch.cuda.Event(enable_timing=True)
            for repeat in range(5):
                order=[('scalar',2),('tensor',0)];rng.shuffle(order)
                for name,mode in order:
                    check(lib.cublasSetMathMode(handle,mode));torch.cuda.synchronize(device)
                    cpu=time.thread_time();wall=time.perf_counter();begin.record(stream)
                    for _ in range(16):run()
                    end.record(stream);submit=time.thread_time()-cpu;end.synchronize()
                    seconds=begin.elapsed_time(end)/1000/16
                    if seconds<=0:raise ValueError('Nonpositive GPU interval')
                    row=dict(repeat=repeat,arithmetic=name,math_mode=mode,loops=16,dimension=N,
                        useful_flops=2*N**3,gpu_seconds=seconds,flops_per_second=2*N**3/seconds,
                        submit_cpu_seconds=submit/16,wall_seconds=(time.perf_counter()-wall)/16)
                    rows.append(row)
                    write_record(root/(device.replace(':','_')+'_'+name+'_'+str(repeat)+'.json'),row)
            result=dict(device=device,device_name=properties.name,torch_version=torch.__version__,cuda_version=torch.version.cuda,
                affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),cublas_version=version.value,
                library=identity,rows=rows,kernels=geometries,max_abs_cpu_reference_error=errors,
                resources={field:statistics.median(r['flops_per_second'] for r in rows if r['arithmetic']==kind)
                    for field,kind in [('fp64_flops_per_second','scalar'),('fp64_tensor_flops_per_second','tensor')]},
                calibration_transfer_qualified=False)
            records.append(result);write_record(root/(device.replace(':','_')+'.json'),result)
            print(json.dumps(dict(device=device,resources=result['resources'],errors=errors)),flush=True)
        finally:check(lib.cublasDestroy_v2(handle))
        del a,b,out,left,right,run;gc.collect();torch.cuda.empty_cache()
    if source!=source_identity():raise ValueError('Source changed during independent FP64 capacity probe')
    write_record(root/'report.json',dict(records=records,source_sha256=source,
        benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        nvidia_before=before,nvidia_after=load(),documentation=DOC,prediction_complete=False,
        scope='Fixed 4096-square GemmEx FP64 resource observations, explicit scalar/pedantic and Tensor Core/default modes verified by untimed kernel census and independent sampled CPU dots. Median of five complete randomized paired repeats; no association durations, shape grid, peak-throughput guarantee or loaded-context qualification.'))


if __name__=='__main__':main()
