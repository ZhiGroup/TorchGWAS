"""Untimed cuSOLVER workspace queries for the eigen JAGWAS factor, with small allocator checks.

The default factor (rcond cutoff, JagwasReduction._eigen_factor) is
torch.linalg.eigh of the K x K FP64 correlation and torch.linalg.qr(mode='r')
of the kept k x K scaled eigenvectors. PyTorch 2.5.1 runs them as one
cusolverDnXsyevd (vectors, lower, n = K, lda = K) and one cusolverDnXgeqrf
(m = k, n = K, lda = k) in place on R, both with NULL params (BatchLinearAlgebra
Lib.cpp apply_syevd / apply_geqrf; a single FP64 matrix never takes the syevj
or batched-cuBLAS branches). k is only known once the spectrum is read, so
geqrf is queried for every k in [1, K] and the largest request is recorded.

No custom GPU kernels, association timings, or large input matrices.
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
    'pytorch_qr_buffers':'https://raw.githubusercontent.com/pytorch/pytorch/v2.5.1/aten/src/ATen/native/BatchLinearAlgebra.cpp',
    'cusolver_syevd':'https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdnxsyevd',
    'cusolver_geqrf':'https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdnxgeqrf',
}
# CUDA 12.4 enums.
VECTOR=1;LOWER=0;R_64F=1


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--device',type=int,default=0)
    parser.add_argument('--traits',type=int,nargs='+',default=[1,7,512,2048,8192,16384])
    args=parser.parse_args();out=Path(args.out)
    if out.exists():raise FileExistsError(out)
    if torch.__version__.split('+')[0]!='2.5.1' or torch.version.cuda!='12.4':
        raise ValueError('This source contract requires PyTorch 2.5.1 CUDA 12.4')
    preferred=str(torch.backends.cuda.preferred_linalg_library())
    if preferred not in ('_LinalgBackend.Default','_LinalgBackend.Cusolver'):
        raise ValueError('Single-matrix cuSOLVER dispatch required: '+preferred)
    torch.cuda.set_device(args.device)
    value=torch.linalg.eigh(torch.ones((1,1),dtype=torch.float64,device='cuda'))
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
    syevd=function('cusolverDnXsyevd_bufferSize',[ct.c_void_p,ct.c_void_p,ct.c_int,ct.c_int,ct.c_int64,
        ct.c_int,ct.c_void_p,ct.c_int64,ct.c_int,ct.c_void_p,ct.c_int,ct.POINTER(ct.c_size_t),ct.POINTER(ct.c_size_t)])
    geqrf=function('cusolverDnXgeqrf_bufferSize',[ct.c_void_p,ct.c_void_p,ct.c_int64,ct.c_int64,
        ct.c_int,ct.c_void_p,ct.c_int64,ct.c_int,ct.c_void_p,ct.c_int,ct.POINTER(ct.c_size_t),ct.POINTER(ct.c_size_t)])
    prop=function('cusolverGetProperty',[ct.c_int,ct.POINTER(ct.c_int)])
    version=[]
    for kind in range(3):
        part=ct.c_int();call(prop,kind,ct.byref(part));version.append(part.value)
    if version[0]<11:raise ValueError('PyTorch 64-bit cuSOLVER API is unavailable')
    handle=ct.c_void_p();rows=[]
    def query_syevd(k,a=None,w=None):
        device=ct.c_size_t();host=ct.c_size_t()
        call(syevd,handle,None,VECTOR,LOWER,k,R_64F,a,max(1,k),R_64F,w,R_64F,ct.byref(device),ct.byref(host))
        return device.value,host.value
    def query_geqrf(m,k,a=None,tau=None):
        device=ct.c_size_t();host=ct.c_size_t()
        call(geqrf,handle,None,m,k,R_64F,a,max(1,m),R_64F,tau,R_64F,ct.byref(device),ct.byref(host))
        return device.value,host.value
    call(create,ct.byref(handle))
    try:
        for k in sorted(set(args.traits)):
            if k<1:raise ValueError('Positive trait counts required')
            syevd_device,syevd_host=query_syevd(k)
            requests=[query_geqrf(m,k) for m in range(1,k+1)]
            device=[request[0] for request in requests];host=[request[1] for request in requests]
            largest=max(range(k),key=lambda i:(device[i],host[i]))
            rows.append(dict(traits=k,syevd_device_workspace_bytes=syevd_device,syevd_host_workspace_bytes=syevd_host,
                geqrf_rows_queried=k,geqrf_full_rank_device_workspace_bytes=device[-1],
                geqrf_full_rank_host_workspace_bytes=host[-1],
                geqrf_max_device_workspace_bytes=max(device),geqrf_max_host_workspace_bytes=max(host),
                geqrf_max_device_rows=largest+1,
                geqrf_device_nondecreasing_in_rows=all(a<=b for a,b in zip(device,device[1:])),
                matrix_allocated=False))
            print(json.dumps(rows[-1]),flush=True)
        # PyTorch passes its matrix pointers; the sizes must not depend on them.
        pointer_checks=[]
        for row in rows:
            k=row['traits']
            if k not in (7,512,2048):continue
            matrix=torch.eye(k,dtype=torch.float64,device='cuda');values=torch.empty(k,dtype=torch.float64,device='cuda')
            tau=torch.empty(k,dtype=torch.float64,device='cuda')
            same=(query_syevd(k,matrix.data_ptr(),values.data_ptr())==(row['syevd_device_workspace_bytes'],row['syevd_host_workspace_bytes'])
                  and query_geqrf(k,k,matrix.data_ptr(),tau.data_ptr())==(row['geqrf_full_rank_device_workspace_bytes'],row['geqrf_full_rank_host_workspace_bytes']))
            pointer_checks.append(dict(traits=k,same_with_matrix_pointers=same))
            del matrix,values,tau
    finally:call(destroy,handle)
    observed=[];rounded=lambda size:((size+511)//512)*512
    by_traits={row['traits']:row for row in rows}
    for k in (7,512,2048):
        if k not in by_traits:continue
        row=by_traits[k]
        correlation=torch.eye(k,dtype=torch.float64,device='cuda')+0.1
        warm=torch.linalg.eigh(correlation);torch.cuda.synchronize();del warm
        gc.collect();torch.cuda.empty_cache()
        baseline=torch.cuda.memory_allocated();torch.cuda.reset_peak_memory_stats()
        values,vectors=torch.linalg.eigh(correlation);torch.cuda.synchronize()
        extra=torch.cuda.max_memory_allocated()-baseline
        # Outputs: eigenvalues and the column-major eigenvector copy of the input.
        explicit=rounded(8*k)+rounded(8*k*k)+rounded(row['syevd_device_workspace_bytes'])
        observed.append(dict(operation='eigh',traits=k,extra_peak_allocated_bytes=extra,
            output_plus_queried_workspace_bytes=explicit,remaining_peak_bytes=extra-explicit,
            remaining_within_4096_bytes=0<=extra-explicit<=4096))
        for m in sorted({k,max(1,k//2)}):
            # The production operand: (vectors[:, -m:] / sqrt(values)).T, a column-major m x K view.
            scaled=(vectors[:,-m:]/values[-m:].abs().sqrt()).T
            warm=torch.linalg.qr(scaled,mode='r');torch.cuda.synchronize();del warm
            gc.collect();torch.cuda.empty_cache()
            baseline=torch.cuda.memory_allocated();torch.cuda.reset_peak_memory_stats()
            result=torch.linalg.qr(scaled,mode='r')[1];torch.cuda.synchronize()
            extra=torch.cuda.max_memory_allocated()-baseline
            device_request=None
            handle=ct.c_void_p();call(create,ct.byref(handle))
            try:device_request=query_geqrf(m,k)[0]
            finally:call(destroy,handle)
            # Outputs: R (m x K, the in-place geqrf buffer) and tau (min(m, K)).
            explicit=rounded(8*m*k)+rounded(8*min(m,k))+rounded(device_request)
            observed.append(dict(operation='qr_r',traits=k,rows=m,extra_peak_allocated_bytes=extra,
                output_plus_queried_workspace_bytes=explicit,remaining_peak_bytes=extra-explicit,
                remaining_within_4096_bytes=0<=extra-explicit<=4096))
            del scaled,result;gc.collect();torch.cuda.empty_cache()
        del correlation,values,vectors;gc.collect();torch.cuda.empty_cache()
    stat=library_path.stat()
    record=dict(method='eigen',rows=rows,pointer_checks=pointer_checks,observations=observed,
        source_sha256=source_identity(),sources=SOURCES,
        benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        torch_version=torch.__version__,cuda_version=torch.version.cuda,preferred_linalg=preferred,
        device=args.device,gpu=torch.cuda.get_device_name(),compute_capability=list(torch.cuda.get_device_capability()),
        library=dict(path=str(library_path),version=version,bytes=stat.st_size,mtime_ns=stat.st_mtime_ns),
        host=os.uname().nodename,durations_recorded=False,prediction_complete=False,
        scope='Exact workspace-size queries for the installed single-matrix FP64 eigen factor: Xsyevd (vectors, lower) '
            'and Xgeqrf over every kept row count. Small allocator controls exclude input allocation. This does not '
            'establish factorization runtime, allocator reservation, host lifetime or complete scan admission.')
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:json.dump(record,stream,indent=2)
    print(json.dumps(dict(pointer_checks=pointer_checks,observations=observed,library=record['library'])),flush=True)


if __name__=='__main__':main()
