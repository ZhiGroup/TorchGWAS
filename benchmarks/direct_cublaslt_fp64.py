"""Owned cuBLASLt handle and explicit FP64 instruction-family preferences."""
import ctypes as C
import hashlib
from pathlib import Path
import torch


class Algo(C.Structure):
    _fields_=[('data',C.c_uint64*8)]


class Heuristic(C.Structure):
    _fields_=[('algo',Algo),('workspaceSize',C.c_size_t),('state',C.c_int),('wavesCount',C.c_float),('reserved',C.c_int*4)]


def check(code):
    if code:raise RuntimeError('cuBLASLt status '+str(code))


class Gemm:
    def __init__(self,a,b,out,workspace):
        self.a,self.b,self.out,self.workspace=a,b,out,workspace
        assert a.dtype==b.dtype==out.dtype==torch.float64 and a.is_contiguous() and b.is_contiguous() and out.is_contiguous()
        assert a.shape==b.shape==out.shape and a.shape[0]==a.shape[1]
        paths={line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines() if '/libcublasLt.so.' in line}
        if len(paths)!=1:raise ValueError('Single Torch-loaded cuBLASLt library required')
        path=Path(paths.pop());self.lib=C.CDLL(str(path));self.owned=[]
        self.identity=dict(path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        pointer=C.c_void_p;pp=C.POINTER(pointer)
        signatures={'cublasLtCreate':[pp],'cublasLtDestroy':[pointer],
            'cublasLtMatmulDescCreate':[pp,C.c_int,C.c_int],'cublasLtMatmulDescDestroy':[pointer],
            'cublasLtMatrixLayoutCreate':[pp,C.c_int,C.c_uint64,C.c_uint64,C.c_int64],'cublasLtMatrixLayoutDestroy':[pointer],
            'cublasLtMatmulPreferenceCreate':[pp],'cublasLtMatmulPreferenceDestroy':[pointer],
            'cublasLtMatmulPreferenceSetAttribute':[pointer,C.c_int,pointer,C.c_size_t],
            'cublasLtMatmulAlgoGetHeuristic':[pointer,pointer,pointer,pointer,pointer,pointer,pointer,C.c_int,C.POINTER(Heuristic),C.POINTER(C.c_int)],
            'cublasLtMatmulAlgoCapGetAttribute':[C.POINTER(Algo),C.c_int,pointer,C.c_size_t,C.POINTER(C.c_size_t)],
            'cublasLtMatmul':[pointer,pointer,pointer,pointer,pointer,pointer,pointer,pointer,pointer,pointer,pointer,pointer,C.POINTER(Algo),pointer,C.c_size_t,pointer]}
        for name,args in signatures.items():getattr(self.lib,name).argtypes=args;getattr(self.lib,name).restype=C.c_int
        try:
            self.handle=self.create('cublasLtCreate','cublasLtDestroy')
            self.desc=self.create('cublasLtMatmulDescCreate','cublasLtMatmulDescDestroy',70,1)
            n=a.shape[0];self.layout=self.create('cublasLtMatrixLayoutCreate','cublasLtMatrixLayoutDestroy',1,n,n,n)
            self.pref=self.create('cublasLtMatmulPreferenceCreate','cublasLtMatmulPreferenceDestroy')
            amount=C.c_uint64(workspace.numel());check(self.lib.cublasLtMatmulPreferenceSetAttribute(self.pref,1,C.byref(amount),C.sizeof(amount)))
            self.alpha,self.beta=C.c_double(1.),C.c_double(0.)
            self.stream=pointer(torch.cuda.current_stream(a.device).cuda_stream)
        except BaseException:
            self.close();raise

    def create(self,name,destroy,*args):
        value=C.c_void_p();check(getattr(self.lib,name)(C.byref(value),*args));self.owned.append((destroy,value));return value

    def select(self,kind):
        # CUDA 12.4 header: PREF_IMPL_MASK=12; all precision flags allowed,
        # only FMA bit 1 or double Tensor Core DMMA bit 8 permitted.
        instruction=1 if kind=='scalar' else 8 if kind=='tensor' else None
        if instruction is None:raise ValueError('Explicit instruction family required')
        mask=C.c_uint64(((1<<64)-1)^0xff|instruction)
        check(self.lib.cublasLtMatmulPreferenceSetAttribute(self.pref,12,C.byref(mask),C.sizeof(mask)))
        result=Heuristic();count=C.c_int()
        check(self.lib.cublasLtMatmulAlgoGetHeuristic(self.handle,self.desc,*([self.layout]*4),self.pref,1,C.byref(result),C.byref(count)))
        if count.value!=1 or result.state:raise ValueError('No supported '+kind+' algorithm under bounded workspace')
        flags=C.c_uint64();written=C.c_size_t()
        check(self.lib.cublasLtMatmulAlgoCapGetAttribute(C.byref(result.algo),15,C.byref(flags),C.sizeof(flags),C.byref(written)))
        if written.value!=8 or flags.value&0xff!=instruction:raise ValueError('Heuristic instruction flags differ from requested family')
        if result.workspaceSize>self.workspace.numel():raise ValueError('Heuristic exceeds workspace cap')
        self.algo=result.algo;self.work_bytes=result.workspaceSize
        return dict(instruction_flags=flags.value,preference_mask=mask.value,workspace_bytes=result.workspaceSize,
            waves_count=result.wavesCount,algorithm_data=list(result.algo.data),selection='First vendor heuristic; no algorithm timing search')

    def run(self):
        check(self.lib.cublasLtMatmul(self.handle,self.desc,C.byref(self.alpha),C.c_void_p(self.b.data_ptr()),self.layout,
            C.c_void_p(self.a.data_ptr()),self.layout,C.byref(self.beta),C.c_void_p(self.out.data_ptr()),self.layout,
            C.c_void_p(self.out.data_ptr()),self.layout,C.byref(self.algo),C.c_void_p(self.workspace.data_ptr()),self.work_bytes,self.stream))

    def close(self):
        while self.owned:
            name,value=self.owned.pop();check(getattr(self.lib,name)(value))
