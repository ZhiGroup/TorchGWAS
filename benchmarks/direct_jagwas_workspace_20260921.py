"""Untimed cuSOLVER workspace queries and small PyTorch allocator checks.

No custom GPU kernels, association timings, or large input matrices. The query
matches PyTorch 2.5.1's single-matrix FP64 lower xpotrf/default-params call.
"""
import argparse
import ctypes as ct
import gc
import hashlib
import json
import os
from pathlib import Path
import torch
from torchgwas.detailed_calibration import source_identity

SOURCES={
    'pytorch_dispatch':'https://raw.githubusercontent.com/pytorch/pytorch/v2.5.1/aten/src/ATen/native/cuda/linalg/BatchLinearAlgebra.cpp',
    'pytorch_workspace':'https://raw.githubusercontent.com/pytorch/pytorch/v2.5.1/aten/src/ATen/native/cuda/linalg/BatchLinearAlgebraLib.cpp',
    'pytorch_api_selection':'https://raw.githubusercontent.com/pytorch/pytorch/v2.5.1/aten/src/ATen/native/cuda/linalg/CUDASolver.h',
    'cusolver_api':'https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdnxpotrf',
}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--device',type=int,default=0)
    args=parser.parse_args();out=Path(args.out)
    if out.exists():raise FileExistsError(out)
    if torch.__version__.split('+')[0]!='2.5.1' or torch.version.cuda!='12.4':
        raise ValueError('This source contract requires PyTorch 2.5.1 CUDA 12.4')
    preferred=str(torch.backends.cuda.preferred_linalg_library())
    if preferred not in ('_LinalgBackend.Default','_LinalgBackend.Cusolver'):
        raise ValueError('Single-matrix cuSOLVER dispatch required: '+preferred)
    torch.cuda.set_device(args.device)
    value=torch.linalg.cholesky(torch.ones((1,1),dtype=torch.float64,device='cuda'))
    torch.cuda.synchronize();del value;gc.collect();torch.cuda.empty_cache()
    paths={line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines() if '/libcusolver.so' in line}
    if len(paths)!=1:raise ValueError('One actual loaded cuSOLVER library required: '+str(paths))
    library_path=Path(next(iter(paths)));library=ct.CDLL(str(library_path))
    def function(name,args):
        fn=getattr(library,name);fn.argtypes=args;fn.restype=ct.c_int
        return fn
    def call(fn,*values):
        status=fn(*values)
        if status:raise RuntimeError(fn.__name__+' returned '+str(status))
    create=function('cusolverDnCreate',[ct.POINTER(ct.c_void_p)])
    destroy=function('cusolverDnDestroy',[ct.c_void_p])
    create_params=function('cusolverDnCreateParams',[ct.POINTER(ct.c_void_p)])
    destroy_params=function('cusolverDnDestroyParams',[ct.c_void_p])
    query=function('cusolverDnXpotrf_bufferSize',[ct.c_void_p,ct.c_void_p,ct.c_int,ct.c_int64,
        ct.c_int,ct.c_void_p,ct.c_int64,ct.c_int,ct.POINTER(ct.c_size_t),ct.POINTER(ct.c_size_t)])
    prop=function('cusolverGetProperty',[ct.c_int,ct.POINTER(ct.c_int)])
    version=[]
    for kind in range(3):
        part=ct.c_int();call(prop,kind,ct.byref(part));version.append(part.value)
    if version[0]<11:raise ValueError('PyTorch 64-bit cuSOLVER API is unavailable')
    handle=ct.c_void_p();params=ct.c_void_p();rows=[]
    call(create,ct.byref(handle))
    try:
        call(create_params,ct.byref(params))
        try:
            for k in [1,7,512,2048,8192,16384]:
                device_bytes=ct.c_size_t();host_bytes=ct.c_size_t()
                # CUDA 12.4 enums: LOWER=0, CUDA_R_64F=1. PyTorch passes a
                # null matrix pointer for this size query, with lda=max(1,K).
                call(query,handle,params,0,k,1,None,max(1,k),1,
                    ct.byref(device_bytes),ct.byref(host_bytes))
                rows.append(dict(traits=k,device_workspace_bytes=device_bytes.value,
                    host_workspace_bytes=host_bytes.value,matrix_allocated=False))
        finally:call(destroy_params,params)
    finally:call(destroy,handle)
    observed=[];rounded=lambda size:((size+511)//512)*512
    for row in rows:
        k=row['traits']
        if k not in [7,512,2048]:continue
        matrix=torch.eye(k,dtype=torch.float64,device='cuda')
        warm=torch.linalg.cholesky(matrix);torch.cuda.synchronize();del warm
        gc.collect();torch.cuda.empty_cache()
        baseline=torch.cuda.memory_allocated();torch.cuda.reset_peak_memory_stats()
        result=torch.linalg.cholesky(matrix);torch.cuda.synchronize()
        extra=torch.cuda.max_memory_allocated()-baseline
        explicit=rounded(8*k*k)+rounded(row['device_workspace_bytes'])
        observed.append(dict(traits=k,extra_peak_allocated_bytes=extra,
            output_plus_queried_workspace_bytes=explicit,remaining_peak_bytes=extra-explicit,
            # Diagnostic only: info/check-error tensors and allocator rounding
            # are still separate from a full model of prepare/scan lifetimes.
            remaining_within_4096_bytes=0<=extra-explicit<=4096))
        del matrix,result;gc.collect();torch.cuda.empty_cache()
    stat=library_path.stat()
    record=dict(rows=rows,observations=observed,source_sha256=source_identity(),sources=SOURCES,
        benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        torch_version=torch.__version__,cuda_version=torch.version.cuda,preferred_linalg=preferred,
        device=args.device,gpu=torch.cuda.get_device_name(),compute_capability=list(torch.cuda.get_device_capability()),
        library=dict(path=str(library_path),version=version,bytes=stat.st_size,mtime_ns=stat.st_mtime_ns),
        host=os.uname().nodename,durations_recorded=False,prediction_complete=False,
        scope='Exact workspace-size query for the installed single-matrix FP64 lower/default cuSOLVER factorization. '
            'Small allocator controls exclude input allocation. This does not establish factorization runtime, '
            'triangular-solve workspace, allocator reservation, host lifetime or complete scan admission.')
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:json.dump(record,stream,indent=2)
    print(json.dumps(dict(rows=rows,observations=observed,library=record['library'])),flush=True)


if __name__=='__main__':main()
