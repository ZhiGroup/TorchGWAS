"""Optional shared native preparation/final statistics; Torch owns the GEMM.

Calls enqueue work on the caller's current CUDA stream and return owned tensors.
They do not synchronize or change global CUDA math/scheduling settings. Inputs
must be contiguous variant-major tensors. Availability errors permit fallback;
launch/runtime errors must propagate. The caller opts in explicitly.
"""
import ctypes as ct
import math
import os
from functools import lru_cache
from pathlib import Path

@lru_cache(maxsize=1)
def load_library():
    path=Path(__file__).with_name('native')/'libtorchgwas_scan.so'
    try:library=ct.CDLL(str(path))
    except OSError as exc:raise ImportError('native statistics unavailable; run bash build_direct_scan.sh') from exc
    try:
        library.tg_scan_abi_version.restype=ct.c_int
        if library.tg_scan_abi_version()!=4:raise ImportError('native statistics ABI mismatch')
    except AttributeError as exc:raise ImportError('native statistics ABI missing') from exc
    pointer=ct.c_void_p;integer=ct.c_int64
    library.tg_scan_error.restype=ct.c_char_p
    library.tg_scan_prepare.argtypes=[pointer,ct.c_int,integer,integer,ct.c_float,ct.c_int,ct.c_float,pointer,pointer,pointer,pointer,pointer,pointer]
    library.tg_scan_prepare.restype=ct.c_int
    library.tg_scan_finish.argtypes=[pointer,pointer,pointer,pointer,pointer,pointer,ct.c_float,integer,integer,integer,ct.c_int,pointer,pointer,pointer,pointer]
    library.tg_scan_finish.restype=ct.c_int
    return library

def resolve_statistics_backend():
    """Resolve explicit opt-in; unavailable native kernels raise without fallback."""
    if os.environ.get('TORCHGWAS_NATIVE_STATS', '0') == '1':
        load_library()
        return 'native_fused'
    return 'torch'


_LAUNCHABLE={}

def available(device=None):
    """True when the native kernels can actually run on this device.

    Loading the shared library is not enough to answer this. The build emits
    `sm_80` and `sm_90` and no PTX, so on any other architecture the library
    loads cleanly, ctypes resolves every symbol, and the failure arrives only at
    launch: `no kernel image is available for execution on the device`. Found on
    an RTX 2080 Ti (Turing, sm_75), where this returned True and nineteen tests
    then failed inside a kernel launch.

    The probe is therefore a real launch on a four-sample variant, cached per
    device -- the answer is a property of the architecture and the library, and
    neither changes within a run.
    """
    try:load_library()
    except ImportError:return False
    try:
        import torch
        if not torch.cuda.is_available():return False
        index=torch.cuda.current_device() if device is None else torch.device(device).index
        if index is None:index=torch.cuda.current_device()
    except Exception:return False
    cached=_LAUNCHABLE.get(index)
    if cached is not None:return cached
    try:
        import torch
        with torch.cuda.device(index):
            raw=torch.zeros((1,64),dtype=torch.uint8,device='cuda:%d'%index)
            prepare(raw,encoding='pgen_2bit',n_samples=4)
            torch.cuda.synchronize(index)
        ok=True
    except Exception:
        # Any failure here means the caller must take the Torch path. The reason
        # is not actionable by them, and raising would turn a working fallback
        # into a failed scan.
        ok=False
    _LAUNCHABLE[index]=ok
    return ok

def _check(tensor,name,device=None,dtype=None):
    if not tensor.is_cuda or not tensor.is_contiguous():
        raise ValueError(name+' must be contiguous and CUDA-resident')
    if device is not None and tensor.device!=device:raise ValueError(name+' is on the wrong CUDA device')
    if dtype is not None and tensor.dtype!=dtype:raise ValueError(name+' has the wrong dtype')

def _launched(library,result):
    if result:raise RuntimeError(library.tg_scan_error().decode('utf8',errors='replace'))

def prepare(raw,scale=1.0,missing_value=None,*,encoding=None,n_samples=None):
    """Return centered float32 dosage, centered SS, raw dosage min, raw dosage max.

    int8 defaults to the PGEN -9 missing sentinel. uint8 is divided by scale;
    float32 may provide an explicit sentinel. NaN min/max propagation and Inf
    range information match Torch. Mean and centered SS use double accumulation,
    while centered storage and all returned arrays are float32.
    """
    import torch
    _check(raw,'raw')
    kinds={torch.int8:0,torch.uint8:1,torch.float32:2}
    if raw.dtype not in kinds or raw.ndim!=2 or min(raw.shape)<=0:raise ValueError('unsupported raw shape/dtype')
    scale=float(scale)
    if not math.isfinite(scale) or scale<=0:raise ValueError('positive finite scale required')
    if raw.dtype==torch.int8 and missing_value is None:missing_value=-9
    packed=encoding in ('pgen_2bit','plink_2bit')
    if encoding not in (None,'pgen_2bit','plink_2bit'):raise ValueError('unsupported preparation encoding')
    if packed:
        if raw.dtype!=torch.uint8 or isinstance(n_samples,bool) or not isinstance(n_samples,int) or n_samples<=0 or raw.shape[1]<(n_samples+3)//4:
            raise ValueError('packed two-bit input requires uint8 rows and a valid logical sample count')
        if scale!=1 or missing_value is not None:raise ValueError('packed two-bit input fixes dosage scale and missing code')
    elif n_samples is not None:raise ValueError('n_samples is only used with packed PGEN')
    library=load_library()
    if packed:
        name=('tg_scan_prepare_bed2' if encoding=='plink_2bit'
              else 'tg_scan_prepare_pgen2')
        try:packed_prepare=getattr(library,name)
        except AttributeError as exc:raise ImportError('rebuild native statistics for packed two-bit support') from exc
        packed_prepare.argtypes=[ct.c_void_p,ct.c_int64,ct.c_int64,ct.c_int64]+[ct.c_void_p]*6
        packed_prepare.restype=ct.c_int
    with torch.cuda.device(raw.device):
        stream=torch.cuda.current_stream(raw.device)
        centered=torch.empty((raw.shape[0],n_samples) if packed else raw.shape,dtype=torch.float32,device=raw.device)
        ss=torch.empty(raw.shape[0],dtype=torch.float32,device=raw.device)
        minimum=torch.empty_like(ss);maximum=torch.empty_like(ss)
        # Observed calls per variant; the variant's residual df follows.
        present=torch.empty(raw.shape[0],dtype=torch.int32,device=raw.device)
        if packed:
            result=packed_prepare(raw.data_ptr(),raw.shape[0],n_samples,raw.shape[1],centered.data_ptr(),ss.data_ptr(),minimum.data_ptr(),maximum.data_ptr(),present.data_ptr(),stream.cuda_stream)
        else:
            result=library.tg_scan_prepare(raw.data_ptr(),kinds[raw.dtype],*raw.shape,scale,
            int(missing_value is not None),float(missing_value or 0),centered.data_ptr(),
            ss.data_ptr(),minimum.data_ptr(),maximum.data_ptr(),present.data_ptr(),stream.cuda_stream)
        raw.record_stream(stream)
        _launched(library,result)
        return centered,ss,minimum,maximum,present

def finish(products,centered_ss,minimum,maximum,phenotype_ss,present,df_offset,validate_range=False):
    """Return beta, t, uint8 status from genotype @ [phenotype, nuisance] products.

    Status codes: 0 valid, 1 nonfinite residual SS, 2 invariant, 3 invalid dosage
    range when requested. Invalid-variant beta/t masking matches shared Torch OLS.
    """
    import torch
    _check(products,'products',dtype=torch.float32)
    if products.ndim!=2:raise ValueError('products must be a matrix')
    variants=products.shape[0]
    tensors=(centered_ss,minimum,maximum,phenotype_ss)
    for name,tensor in zip(('centered_ss','minimum','maximum','phenotype_ss'),tensors):
        _check(tensor,name,products.device,torch.float32)
        if tensor.ndim!=1:raise ValueError(name+' must be a vector')
    if any(t.numel()!=variants for t in tensors[:3]):raise ValueError('variant dimensions differ')
    traits=phenotype_ss.numel();covariates=products.shape[1]-traits
    _check(present,'present',products.device,torch.int32)
    if present.ndim!=1 or present.numel()!=variants:raise ValueError('present must be one count per variant')
    if variants<=0 or traits<=0 or covariates<0 or not math.isfinite(float(df_offset)):raise ValueError('invalid statistics dimensions/df')
    library=load_library()
    with torch.cuda.device(products.device):
        stream=torch.cuda.current_stream(products.device)
        beta=torch.empty((variants,traits),dtype=torch.float32,device=products.device)
        statistic=torch.empty_like(beta)
        status=torch.empty(variants,dtype=torch.uint8,device=products.device)
        result=library.tg_scan_finish(products.data_ptr(),centered_ss.data_ptr(),minimum.data_ptr(),maximum.data_ptr(),
            phenotype_ss.data_ptr(),present.data_ptr(),ct.c_float(float(df_offset)),variants,traits,covariates,int(validate_range),
            beta.data_ptr(),statistic.data_ptr(),status.data_ptr(),stream.cuda_stream)
        for tensor in (products,*tensors,present):tensor.record_stream(stream)
        _launched(library,result)
        return beta,statistic,status
