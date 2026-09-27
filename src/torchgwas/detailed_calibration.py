"""Bind independent calculator profiles to their actual execution context.

This packages supplied component prices; it neither collects association
timings nor certifies price accuracy. Workload geometry and search bounds stay
outside the hardware profile so an input cannot silently carry an old plan.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

SCHEMA='torchgwas.detailed_calibration.v1'
PRICE_SCHEMA='torchgwas.detailed_calibration.v2'
PROFILE_FIELDS={'schema','bound_at_utc','contexts','execution_context','source_sha256',
    'component_artifacts','limitations','selection_validated','runtime_prediction_validated','scope'}
ENVIRONMENT_PREFIXES=('CUDA_','CUBLAS_','NVIDIA_','PYTORCH_','TORCH_','TORCHGWAS_',
    'OMP_','GOMP_','KMP_','MKL_','OPENBLAS_','MALLOC_','NPY_','NUMPY_')
ENVIRONMENT_KEYS=(
    'CUDA_VISIBLE_DEVICES','CUDA_DEVICE_ORDER','CUBLAS_WORKSPACE_CONFIG',
    'CUDA_LAUNCH_BLOCKING','CUDA_MODULE_LOADING','PYTORCH_NO_CUDA_MEMORY_CACHING',
    'PYTORCH_CUDA_ALLOC_CONF','PYTORCH_ALLOC_CONF','OMP_NUM_THREADS',
    'NVIDIA_TF32_OVERRIDE','TORCH_ALLOW_TF32_CUBLAS_OVERRIDE',
    'OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_WAIT_POLICY','GOMP_SPINCOUNT',
    'NUMPY_MADVISE_HUGEPAGE','MALLOC_MMAP_THRESHOLD_','MALLOC_TRIM_THRESHOLD_',
    'MALLOC_ARENA_MAX','MALLOC_MMAP_MAX_','GLIBC_TUNABLES','TORCHGWAS_NATIVE_STATS',
    'TORCHGWAS_PGEN_PACKED','TORCHGWAS_PGEN_BACKEND','TORCHGWAS_BLOCKING_EVENTS',
    'TORCHGWAS_SCAN_PROFILE')


def sha256_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8<<20),b''):h.update(block)
    return h.hexdigest()


def source_identity(*,digest_cache=None):
    paths=[p for p in sorted(Path(__file__).parent.iterdir()) if p.is_file() and p.suffix in ('.py','.cpp')]
    digests=[sha256_file(p) for p in paths] if digest_cache is None else digest_cache.digests(paths)
    return {p.name:digest for p,digest in zip(paths,digests)}


def _optional_text(path):
    try:return Path(path).read_text().strip()
    except FileNotFoundError:return None


def _numpy_core_context(*,digest_cache=None):
    """Bind installed core bytes and freshly read CPU dispatch capabilities.

    The version alone cannot distinguish builds or import-time CPU feature
    controls. Only the file digest may be reused; runtime metadata is reread.
    These NumPy runtime attributes must be available for a qualified profile.
    """
    from numpy._core import _multiarray_umath as core
    path=getattr(core,'__file__',None)
    features=getattr(core,'__cpu_features__',None)
    baseline=getattr(core,'__cpu_baseline__',None)
    dispatch=getattr(core,'__cpu_dispatch__',None)
    if not isinstance(path,str) or not path:
        raise ValueError('NumPy core binary identity unavailable')
    if (not isinstance(features,dict) or not features or
            any(not isinstance(k,str) or not k or type(v) is not bool for k,v in features.items())):
        raise ValueError('NumPy runtime CPU features unavailable')
    for name,values in [('baseline',baseline),('dispatch',dispatch)]:
        if (not isinstance(values,(list,tuple)) or
                any(not isinstance(v,str) or not v for v in values)):
            raise ValueError('NumPy CPU '+name+' unavailable')
    digest=sha256_file(path) if digest_cache is None else digest_cache.digests([path])[0]
    return dict(library_sha256=digest,cpu_features=dict(sorted(features.items())),
                cpu_baseline=list(baseline),cpu_dispatch=list(dispatch))


def storage_identity(path):
    """Identify the mounted backing store, allowing a not-yet-created output."""
    target=Path(path).expanduser().absolute()
    existing=target
    while not existing.exists():
        if existing==existing.parent:raise ValueError('No existing output ancestor')
        existing=existing.parent
    existing=existing.resolve(strict=True)
    data=json.loads(subprocess.check_output(['findmnt','-J','-T',str(existing),'-o','SOURCE,FSTYPE,TARGET'],text=True))
    filesystems=data.get('filesystems',[])
    if len(filesystems)!=1:raise ValueError('Ambiguous storage mount')
    mount=filesystems[0];stat=existing.stat()
    return dict(source=mount['source'],fstype=mount['fstype'],target=mount['target'],
                device_major=os.major(stat.st_dev),device_minor=os.minor(stat.st_dev))


def execution_context(devices, *, input_path, output_path, digest_cache=None):
    """Read identities/settings without changing affinity, flags or threads.

    Capturing device properties initializes CUDA. Call outside timed scans.
    This is a context binding, not a measurement of current available capacity.
    """
    import numpy as np
    import torch
    from .gpu_identity import canonical_cuda_device,physical_gpu_identity
    from .pgen_native import load_library
    from .numa_context import memory_policy_context
    if sys.platform!='linux' or not hasattr(os,'sched_getaffinity'):
        raise ValueError('Detailed calibration context currently requires Linux')
    devices=list(devices)
    if not devices or len(set(devices))!=len(devices):raise ValueError('Unique explicit devices required')
    for device in devices:canonical_cuda_device(device)
    properties={device:torch.cuda.get_device_properties(device) for device in devices}
    physical=physical_gpu_identity([str(getattr(p,'uuid','')).lower() for p in properties.values()])
    gpu={}
    for device in devices:
        p=properties[device]
        uuid=str(getattr(p,'uuid','')).lower()
        if uuid not in physical:raise ValueError('CUDA UUID must identify a measured physical GPU: '+device)
        gpu[device]=dict(uuid=uuid,name=p.name,total_memory=p.total_memory,
            sm_count=p.multi_processor_count,max_threads_per_sm=p.max_threads_per_multi_processor,
            compute_capability=[p.major,p.minor],**physical[uuid])
        domain,bus,slot=physical[uuid]['pci_bus_id'].split(':')
        pci=Path('/sys/bus/pci/devices')/f'{int(domain,16):04x}:{bus.lower()}:{slot.lower()}'
        gpu[device]['numa_node']=_optional_text(pci/'numa_node')
        gpu[device]['max_link_speed']=_optional_text(pci/'max_link_speed')
        gpu[device]['max_link_width']=_optional_text(pci/'max_link_width')
    machine=Path('/etc/machine-id').read_text().strip()
    if not machine:raise ValueError('Machine identity unavailable')
    # Bind loaded BLAS/OpenMP thread counts as well as environment declarations:
    # a threadpoolctl change after import must invalidate the old price context.
    try:
        from threadpoolctl import threadpool_info
    except ImportError as e:
        raise ValueError('threadpoolctl is required to verify numerical CPU pools') from e
    current_pools=threadpool_info()
    pool_paths=list(dict.fromkeys(p['filepath'] for p in current_pools))
    pool_hashes=([sha256_file(path) for path in pool_paths] if digest_cache is None or not pool_paths else digest_cache.digests(pool_paths))
    pool_hashes=dict(zip(pool_paths,pool_hashes))
    pools=sorted([dict(internal_api=p.get('internal_api'),prefix=p.get('prefix'),
        version=p.get('version'),num_threads=p.get('num_threads'),
        library_sha256=pool_hashes[p['filepath']]) for p in current_pools],
        key=lambda p:(str(p['internal_api']),str(p['prefix']),p['library_sha256']))
    affinity=sorted(os.sched_getaffinity(0));topology={}
    for cpu in affinity:
        path=Path('/sys/devices/system/cpu')/f'cpu{cpu}'
        topology[str(cpu)]={key:_optional_text(path/'topology'/key) for key in ['core_id','physical_package_id','thread_siblings_list']}
        topology[str(cpu)]['numa_nodes']=sorted(p.name for p in path.glob('node[0-9]*'))
    native_path=load_library()._name
    native_hash=sha256_file(native_path) if digest_cache is None else digest_cache.digests([native_path])[0]
    from .host_significance import host_selector,NATIVE_HOST_SELECTOR
    predicate_context={}
    if host_selector()==NATIVE_HOST_SELECTOR:
        from .native_host_predicate import context
        predicate_context['host_predicate']=context(digest_cache=digest_cache)
    return dict(machine_sha256=hashlib.sha256(machine.encode()).hexdigest(),
        kernel=platform.release(),architecture=platform.machine(),libc=list(platform.libc_ver()),
        python_version=sys.version,numpy_version=np.__version__,
        numpy_core=_numpy_core_context(digest_cache=digest_cache),torch_version=torch.__version__,
        cuda_runtime=torch.version.cuda,affinity=affinity,cpu_topology=topology,numa_policy=memory_policy_context(),
        page_bytes=os.sysconf('SC_PAGE_SIZE'),transparent_hugepage_policy=_optional_text('/sys/kernel/mm/transparent_hugepage/enabled'),
        torch_threads=torch.get_num_threads(),torch_interop_threads=torch.get_num_interop_threads(),
        numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
        torch_default_dtype=str(torch.get_default_dtype()),
        allow_tf32=bool(torch.backends.cuda.matmul.allow_tf32),float32_matmul_precision=torch.get_float32_matmul_precision(),cpu_pools=pools,
        environment={key:os.getenv(key) for key in sorted(set(ENVIRONMENT_KEYS)|{k for k in os.environ if k.startswith(ENVIRONMENT_PREFIXES)})},devices=gpu,
        storage=dict(input=storage_identity(input_path),output=storage_identity(output_path)),
        native_library_sha256=native_hash,**predicate_context)


def _context_devices(contexts):
    if not isinstance(contexts,list) or not contexts:raise ValueError('Independent contexts required')
    names=set();devices=[]
    for context in contexts:
        name=context.get('name')
        if not isinstance(name,str) or not name or name in names:raise ValueError('Unique context names required')
        names.add(name)
        transfer=context.get('shared_transfer_capacities')
        if transfer is not None and (not isinstance(transfer,dict) or
                set(transfer)!={'h2d','d2h'} or any(
                    isinstance(value,bool) or not isinstance(value,(int,float)) or
                    not math.isfinite(value) or value<=0
                    for value in transfer.values())):
            raise ValueError('Positive shared H2D/D2H context capacities required')
        group=context.get('devices',[])
        if not group or len(set(group))!=len(group) or set(context.get('profiles',{}))!=set(group):
            raise ValueError('One profile per active context device required')
        for device in group:
            if device not in devices:devices.append(device)
    return devices


def _differences(expected,actual,prefix='context'):
    if isinstance(expected,dict) and isinstance(actual,dict):
        rows=[]
        for key in sorted(set(expected)|set(actual)):
            name=prefix+'.'+str(key)
            if key not in expected or key not in actual:rows.append(name)
            else:rows.extend(_differences(expected[key],actual[key],name))
        return rows
    return [] if expected==actual else [prefix]


def bind_detailed_profile(contexts, execution, *, component_artifacts, limitations, sources=None, price_bindings=None):
    """Package caller-validated independent prices and an explicit live binding.

    Binding time is not measurement time. Artifact checks protect provenance;
    each producer still has to validate that its prices describe this context.
    """
    devices=_context_devices(contexts)
    if set(devices)!=set(execution.get('devices',{})):raise ValueError('Execution device coverage differs')
    for context in contexts:
        for device,profile in context['profiles'].items():
            gpu=execution['devices'][device]
            if 'sm_count' in profile.get('gpu_resources',{}) and profile['gpu_resources']['sm_count']!=gpu.get('sm_count'):
                raise ValueError('Priced GPU properties differ: '+device)
            for key,value in profile.get('reduction_gpu_properties',{}).items():
                if gpu.get(key)!=value:raise ValueError('Priced reduction properties differ: '+device+'.'+key)
            if profile.get('pageable_host_service') is not None:
                prices=profile['pageable_host_service']['prices'];measured=prices['context']
                if prices['devices']!=context['devices'] or prices['device']!=device:
                    raise ValueError('Pageable worker prices differ from context')
                for key in ['numpy_version','torch_version','python_version','libc','affinity','torch_threads','numpy_madvise_hugepage']:
                    if measured[key]!=execution.get(key):raise ValueError('Measured pageable context differs: '+key)
                for key,value in measured['allocator_environment'].items():
                    if execution.get('environment',{}).get(key)!=value:
                        raise ValueError('Measured allocator environment differs: '+key)
    if not isinstance(component_artifacts,dict) or not component_artifacts:
        raise ValueError('Component artifact hashes required')
    artifacts={}
    for path,digest in component_artifacts.items():
        resolved=Path(path).expanduser().resolve(strict=True)
        if not isinstance(digest,str) or len(digest)!=64 or sha256_file(resolved)!=digest:
            raise ValueError('Component artifact changed: '+str(path))
        if str(resolved) in artifacts and artifacts[str(resolved)]!=digest:
            raise ValueError('Conflicting component identity: '+str(resolved))
        artifacts[str(resolved)]=digest
    if not isinstance(limitations,list) or any(not isinstance(x,str) or not x for x in limitations):
        raise ValueError('Explicit limitation list required')
    result=dict(schema=SCHEMA,bound_at_utc=datetime.now(timezone.utc).isoformat(),
        contexts=copy.deepcopy(contexts),execution_context=copy.deepcopy(execution),
        source_sha256=source_identity() if sources is None else copy.deepcopy(sources),
        component_artifacts=artifacts,limitations=list(limitations),
        selection_validated=False,runtime_prediction_validated=False,
        scope='Independent component profiles bound to source, hardware, CPU/allocator/backend settings and mounted storage. No workload, automatic capacity rescaling or association-time fit. Binding does not certify current load, externally imposed allocator caps or prediction accuracy.')
    if price_bindings is not None:
        from .price_binding import canonical_price_bindings,validate_price_bindings
        result.update(schema=PRICE_SCHEMA,price_bindings=canonical_price_bindings(price_bindings))
        validate_price_bindings(result)
    return result


def validate_detailed_profile(profile, current, *, sources=None, verify_artifacts=True, deferred_bindings=()):
    _profile_schema(profile)
    if type(verify_artifacts) is not bool:raise ValueError('verify_artifacts must be boolean')
    devices=_context_devices(profile['contexts'])
    if set(devices)!=set(profile['execution_context'].get('devices',{})):
        raise ValueError('Execution device coverage differs')
    expected_sources=source_identity() if sources is None else sources
    changed=_differences(profile['source_sha256'],expected_sources,'source')
    changed+=_differences(profile['execution_context'],current)
    if changed:raise ValueError('Detailed calibration context changed: '+', '.join(changed))
    if verify_artifacts:
        for path,digest in profile['component_artifacts'].items():
            if sha256_file(path)!=digest:raise ValueError('Component artifact changed: '+path)
    from .price_binding import validate_price_bindings
    prices=validate_price_bindings(profile,deferred_bindings=deferred_bindings)
    return dict(devices=devices,context_matches=True,artifacts_verified=verify_artifacts,price_evidence=prices,
        limitations=list(profile['limitations'])+['Live CPU/GPU/storage contention and allocator history are not certified by identity checks.'])


def _profile_schema(profile):
    if not isinstance(profile,dict) or profile.get('schema') not in (SCHEMA,PRICE_SCHEMA):
        raise ValueError('Unknown detailed calibration schema')
    fields=PROFILE_FIELDS|({'price_bindings'} if profile['schema']==PRICE_SCHEMA else set())
    if set(profile)!=fields:raise ValueError('Detailed calibration fields differ; workload and plans belong outside the profile')
    if profile['selection_validated'] is not False or profile['runtime_prediction_validated'] is not False:
        raise ValueError('Context binding cannot certify selection or prediction accuracy')
    for key in ['source_sha256','component_artifacts']:
        if not isinstance(profile[key],dict) or not profile[key]:raise ValueError('Nonempty '+key+' required')


def write_detailed_profile(profile,path):
    """Publish a complete JSON profile atomically; never replace another one."""
    _profile_schema(profile)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    data=json.dumps(profile,indent=2,allow_nan=False)+'\n'
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=path.parent,prefix='.'+path.name+'.',delete=False) as stream:
            temporary=Path(stream.name);stream.write(data);stream.flush();os.fsync(stream.fileno())
        os.link(temporary,path)  # Atomic no-overwrite publication on one filesystem.
        return path
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)


def read_detailed_profile(path):
    def reject_constant(value):raise ValueError('Non-finite profile value: '+value)
    def finite_float(value):
        result=float(value)
        if not math.isfinite(result):reject_constant(value)
        return result
    def unique_object(pairs):
        result={}
        for key,value in pairs:
            if key in result:raise ValueError('Duplicate profile key: '+key)
            result[key]=value
        return result
    with Path(path).open(encoding='utf-8') as stream:
        profile=json.load(stream,parse_constant=reject_constant,parse_float=finite_float,object_pairs_hook=unique_object)
    _profile_schema(profile)
    return profile
