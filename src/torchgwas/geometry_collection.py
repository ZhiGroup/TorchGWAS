"""Complete missing duration-free geometry for the existing detailed calculator.

Component prices remain unchanged. Probe requests come only from feasible
candidates, run in fresh child processes, and publish no measured durations.
"""
from __future__ import annotations

import argparse
import copy
import ctypes
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

from .detailed_calibration import (bind_detailed_profile,execution_context,
    sha256_file,source_identity,validate_detailed_profile,write_detailed_profile)
from .mechanistic_plan import _integer
from .pageable_host_service import setup_host_allocations,validate_pageable_geometry
from .setup_work import setup_work
from .trait_candidate_space import prepare_trait_candidates,bounded_trait_plan
from .trait_tiling_model import trait_tiled_shape,trait_tiled_memory,reuse_trait_work

SCHEMA='torchgwas.duration_free_geometry.v1'
GEOMETRY_FIELDS={'grid','block','registers per thread','shared memory'}


def write_record(path,value):
    """Atomic no-overwrite publication; partial worker records are never valid."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=path.parent,delete=False) as stream:
            temporary=Path(stream.name);json.dump(value,stream,indent=2,allow_nan=False)
            stream.write('\n');stream.flush();os.fsync(stream.fileno())
        os.link(temporary,path)
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)


def kernel_census(events):
    """Whitelist launch work; raw profiler timing fields cannot enter a profile."""
    rows=[]
    for event in events:
        if event.get('cat')!='kernel':continue
        geometry={key:value for key,value in event.get('args',{}).items() if key in GEOMETRY_FIELDS}
        for name in ['grid','block']:
            value=geometry.get(name)
            if not isinstance(value,list) or len(value)!=3 or any(type(v) is not int or v<1 for v in value):
                raise ValueError('Kernel census requires a positive three-dimensional '+name)
        for name in GEOMETRY_FIELDS-{'grid','block'}:
            if name in geometry and (type(geometry[name]) is not int or geometry[name]<0):
                raise ValueError('Invalid kernel '+name)
        if not isinstance(event.get('name'),str) or not event['name']:raise ValueError('Kernel name required')
        rows.append(dict(name=event['name'],geometry=geometry))
    if not rows:raise ValueError('Profiler returned no CUDA kernels')
    return rows


@reuse_trait_work
def geometry_requirements(workload,contexts,*,bounds,joint,output,max_kernel_shapes=128,
                          max_host_extents=128,max_device_probe_bytes=2<<30,
                          max_host_probe_bytes=1<<30,max_host_touch_bytes=8<<30):
    """Use the same source memory ledger before proposing any real allocation."""
    limits=dict(max_kernel_shapes=max_kernel_shapes,max_host_extents=max_host_extents,
        max_device_probe_bytes=max_device_probe_bytes,max_host_probe_bytes=max_host_probe_bytes,
        max_host_touch_bytes=max_host_touch_bytes)
    for key,value in limits.items():_integer(key,value)
    space=prepare_trait_candidates(workload,contexts,output=output,**bounds)
    by_name={context['name']:context for context in contexts}
    gpu={};host={};rejected=[];feasible=[];ambiguities=[]
    checked_host=set()
    for index,candidate in enumerate(space['candidates']):
        shape=trait_tiled_shape(candidate)
        if not set(shape['devices'])<=set(joint['device_memory_bytes']):raise ValueError('Missing device memory budget')
        if shape['reader_workers']>joint['cpu_workers']:
            rejected.append(dict(candidate_index=index,reason='shared_reader_budget'));continue
        memory=trait_tiled_memory(candidate,host_reserve_bytes=joint.get('host_reserve_bytes',0),
            device_reserve_bytes=joint.get('device_reserve_bytes',0),device_memory_profiles=joint.get('device_memory_profiles'))
        reason=('host_memory' if memory['host_bytes']>joint['host_memory_bytes'] else
                'device_memory' if any(memory['device_bytes'][d]>joint['device_memory_bytes'][d] for d in shape['devices']) else None)
        if reason:
            rejected.append(dict(candidate_index=index,reason=reason));continue
        feasible.append(index)
        context_name=space['assignments'][index]['context']
        for tile in candidate['tiles']:
            data=tile['data'];device=tile['device'];profile=tile['profile']
            n,m,k,c=[data[key] for key in ['samples','markers','traits_analyzed','covariates']]
            chunk=profile['chunk_markers'];sizes={min(chunk,m)}
            if m%chunk:sizes.add(m%chunk)
            for b in sorted(sizes):
                dimensions=[n,b,k,c]
                matches=[row for row in profile.get('kernel_geometry',[]) if
                    [row.get(key) for key in ['N','B','K','C']]==dimensions and row.get('validate_range') is True]
                if len(matches)>1:
                    ambiguities.append(dict(context=context_name,device=device,shape=dimensions));continue
                if matches:continue
                if memory['device_bytes'][device]>max_device_probe_bytes:
                    raise ValueError('Required geometry exceeds max_device_probe_bytes before collection')
                key=(device,*dimensions)
                row=gpu.setdefault(key,dict(device=device,shape=dimensions,validate_range=True,
                    contexts=[],candidate_indices=[],device_budget_bytes=min(max_device_probe_bytes,joint['device_memory_bytes'][device]),
                    source_memory_bound_bytes=memory['device_bytes'][device]))
                row['source_memory_bound_bytes']=max(row['source_memory_bound_bytes'],memory['device_bytes'][device])
                if context_name not in row['contexts']:row['contexts'].append(context_name)
                if index not in row['candidate_indices']:row['candidate_indices'].append(index)
            host_inputs=profile.get('pageable_host_service')
            if host_inputs is None:continue
            identity=(context_name,device)
            if identity not in checked_host:
                validate_pageable_geometry(host_inputs['geometry'],host_inputs['prices']);checked_host.add(identity)
            work=setup_work(n,k,c,reuse_observed_counts=True,input_contiguous=data.get('phenotype_c_contiguous'))
            for array in setup_host_allocations(work):
                kind,size=array['kind'],array['bytes']
                if str(size) in host_inputs['geometry']['arrays'].get(kind,{}):continue
                if size>max_host_probe_bytes:raise ValueError('Required allocator extent exceeds max_host_probe_bytes before collection')
                row=host.setdefault((kind,size),dict(kind=kind,bytes=size,targets=[]))
                if list(identity) not in row['targets']:row['targets'].append(list(identity))
    if ambiguities:raise ValueError('Ambiguous existing kernel geometry: '+str(ambiguities))
    if not feasible:raise ValueError('No memory-feasible candidate for geometry collection: '+str(rejected))
    if len(gpu)>max_kernel_shapes:raise ValueError('Missing geometry exceeds max_kernel_shapes before collection')
    if len(host)>max_host_extents:raise ValueError('Missing geometry exceeds max_host_extents before collection')
    touched=4*sum(row['bytes'] for row in host.values())
    if touched>max_host_touch_bytes:raise ValueError('Allocator geometry exceeds max_host_touch_bytes before collection')
    return dict(schema=SCHEMA,kernel_requests=list(gpu.values()),host_requests=list(host.values()),
        memory_feasible_candidates=feasible,memory_rejected=rejected,limits=limits,host_touch_bytes=touched,
        input_file_identity=space['input_file_identity'],source_sha256=source_identity(),
        scope='Only missing geometry for candidates admitted by the existing memory/reader budgets. GPU requests deduplicate by device and exact shape; host requests deduplicate by allocator and extent. No association timings or new component prices.')


def _gpu_geometry(requests,directory):
    import torch
    from .linear import _dosage_statistics
    rows=[]
    for request in requests:
        device=request['device'];n,b,k,c=request['shape']
        torch.cuda.set_device(device)
        available=torch.cuda.mem_get_info(device)[0]+torch.cuda.memory_reserved(device)-torch.cuda.memory_allocated(device)
        if request['source_memory_bound_bytes']>available:raise ValueError('Insufficient live memory for geometry capture')
        total=torch.cuda.get_device_properties(device).total_memory
        torch.cuda.set_per_process_memory_fraction(min(request['device_budget_bytes'],total)/total,device)
        x=torch.randint(0,3,(b,n),dtype=torch.int8,device=device)
        design=torch.randn((n,k+c+1),dtype=torch.float32,device=device)
        ss=torch.full((k,),float(n),dtype=torch.float32,device=device)
        def compute():
            g=torch.where(x==-9,torch.nan,x.to(torch.float32))
            return _dosage_statistics(g,design,ss,k,n-c-2,True,covariate_rank=c)
        value=compute();del value;torch.cuda.synchronize(device)
        with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,torch.profiler.ProfilerActivity.CUDA]) as profiler:
            value=compute();torch.cuda.synchronize(device)
        del value
        with tempfile.TemporaryDirectory(dir=directory,prefix='trace-') as temporary:
            trace=Path(temporary)/'trace.json';profiler.export_chrome_trace(str(trace))
            kernels=kernel_census(json.loads(trace.read_text())['traceEvents'])
        rows.append(dict(device=device,N=n,B=b,K=k,C=c,validate_range=True,kernels=kernels))
        del x,design,ss,compute,profiler
        torch.cuda.empty_cache()
    return dict(rows=rows)


class _Mallinfo(ctypes.Structure):
    _fields_=[(name,ctypes.c_size_t) for name in
        ['arena','ordblks','smblks','hblks','hblkhd','usmblks','fsmblks','uordblks','fordblks','keepcost']]


def validate_host_capture(geometry):
    """Reject duration fields at every level of a fresh allocator artifact."""
    fields={'arrays','durations_recorded','page_bytes','numpy_version','torch_version',
        'python_version','libc','affinity','torch_threads','numpy_madvise_hugepage',
        'allocator_environment','source_sha256','scope'}
    if set(geometry)!=fields:raise ValueError('Unexpected allocator geometry fields')
    for bank in geometry['arrays'].values():
        for row in bank.values():
            if set(row)!={'route','mapped_bytes','observations'}:
                raise ValueError('Unexpected allocator extent fields')
            for observation in row['observations']:
                if set(observation)!={'repeat','route','allocate_count','allocate_bytes','release_count','release_bytes'}:
                    raise ValueError('Unexpected allocator observation fields')

def _host_geometry(requests):
    import numpy as np
    import torch
    from .api import _available_host_bytes
    if requests and max(row['bytes'] for row in requests)>_available_host_bytes():
        raise ValueError('Insufficient live host memory for allocator capture')
    sizes={kind:sorted(row['bytes'] for row in requests if row['kind']==kind) for kind in ['numpy','torch']}
    def collect():
        libc=ctypes.CDLL(None);libc.mallinfo2.restype=_Mallinfo
        arrays={}
        for kind,extents in sizes.items():
            allocate=(lambda n:np.empty(n,dtype=np.uint8)) if kind=='numpy' else (lambda n:torch.empty(n,dtype=torch.uint8,device='cpu'))
            value=allocate(128);del value
            records={size:[] for size in extents}
            for repeat in range(4):
                for size in sorted(extents,reverse=bool(repeat%2)):
                    before=libc.mallinfo2();value=allocate(size);allocated=libc.mallinfo2()
                    if kind=='numpy':value.fill(0)
                    else:value.fill_(0)
                    del value;freed=libc.mallinfo2()
                    delta=dict(allocate_count=allocated.hblks-before.hblks,allocate_bytes=allocated.hblkhd-before.hblkhd,
                        release_count=freed.hblks-allocated.hblks,release_bytes=freed.hblkhd-allocated.hblkhd)
                    if delta['allocate_count']==1 and delta['release_count']==-1 and delta['allocate_bytes']>=size and delta['release_bytes']==-delta['allocate_bytes']:
                        route='mmap'
                    elif not any(delta.values()):route='arena'
                    else:raise ValueError('Unresolved allocator geometry: '+str((kind,size,delta)))
                    records[size].append(dict(repeat=repeat,route=route,**delta))
            arrays[kind]={str(size):dict(route=rows[0]['route'] if len({r['route'] for r in rows})==1 else 'variable',
                mapped_bytes=rows[0]['allocate_bytes'] if all(r['route']=='mmap' and r['allocate_bytes']==rows[0]['allocate_bytes'] for r in rows) else None,
                observations=rows) for size,rows in records.items()}
        return arrays
    with ThreadPoolExecutor(max_workers=1) as pool:arrays=pool.submit(collect).result()
    return dict(arrays=arrays,durations_recorded=False,page_bytes=os.sysconf('SC_PAGE_SIZE'),
        numpy_version=np.__version__,torch_version=torch.__version__,python_version=sys.version,
        libc=list(platform.libc_ver()),affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),
        numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
        allocator_environment={key:os.getenv(key) for key in ['MALLOC_MMAP_THRESHOLD_','MALLOC_TRIM_THRESHOLD_','GLIBC_TUNABLES']},
        source_sha256={str(Path(__file__).resolve()):sha256_file(__file__)},
        scope='Duration-free four-pass allocation-route census in one worker. Future allocator history and concurrent-route behavior remain conditional.')


def _worker(request_path,result_path):
    import torch
    request=json.loads(Path(request_path).read_text())
    if request['schema']!=SCHEMA or request['source_sha256']!=source_identity():raise ValueError('Worker source changed')
    context=request['execution_context']
    torch.set_num_threads(context['torch_threads']);torch.set_num_interop_threads(context['torch_interop_threads'])
    torch.backends.cuda.matmul.allow_tf32=context['allow_tf32']
    actual=execution_context(list(context['devices']),input_path=request['input_path'],output_path=request['output_path'])
    if actual!=context:raise ValueError('Collector execution context differs from calibrated parent')
    result=(_gpu_geometry(request['requests'],Path(result_path).parent) if request['kind']=='gpu' else
            _host_geometry(request['requests']) if request['kind']=='host' else None)
    if result is None:raise ValueError('Unknown geometry worker')
    if request['source_sha256']!=source_identity():raise ValueError('Source changed during collection')
    write_record(result_path,dict(schema=SCHEMA,kind=request['kind'],source_sha256=request['source_sha256'],
        execution_context=context,durations_recorded=False,request_sha256=sha256_file(request_path),result=result))


def _run_worker(kind,requests,root,controller,timeout):
    request_path=root/(kind+'_request.json');result_path=root/(kind+'_geometry.json')
    write_record(request_path,dict(schema=SCHEMA,kind=kind,requests=requests,execution_context=controller.context,
        source_sha256=source_identity(),input_path=controller.input_path,output_path=controller.output_path))
    with (root/(kind+'_worker.log')).open('x') as stream:
        subprocess.run([sys.executable,'-m','torchgwas.geometry_collection','worker',
            '--request',str(request_path),'--result',str(result_path)],stdout=stream,stderr=subprocess.STDOUT,
            check=True,timeout=timeout)
    record=json.loads(result_path.read_text())
    if (record.get('schema')!=SCHEMA or record.get('kind')!=kind or record.get('durations_recorded') is not False
        or record.get('request_sha256')!=sha256_file(request_path) or record.get('source_sha256')!=source_identity()
        or record.get('execution_context')!=controller.context):raise ValueError('Invalid geometry worker artifact')
    return record['result'],result_path


def _price_contexts(contexts):
    result=copy.deepcopy(contexts)
    for context in result:
        for profile in context['profiles'].values():
            profile.pop('kernel_geometry',None)
            if 'pageable_host_service' in profile:profile['pageable_host_service'].pop('geometry',None)
    return result


def complete_profile_geometry(profile,workload,config,*,output,output_path,collection_dir,
                              collect=True,worker_timeout=300,**limits):
    from .detailed_autotune import DetailedAutotune
    from .analytical_plan_cache import input_identity
    _integer('worker_timeout',worker_timeout)
    if type(collect) is not bool:raise ValueError('collect must be boolean')
    controller=DetailedAutotune(profile,config,input_path=workload['genotype'],output_path=output_path)
    before_input=input_identity(controller.input_path)
    coverage=geometry_requirements(workload,controller.profile['contexts'],bounds=controller.config['bounds'],
                                   joint=controller.config['joint'],output=output,**limits)
    if not collect:return coverage
    root=Path(collection_dir);root.mkdir(parents=True,exist_ok=False)
    write_record(root/'coverage.json',coverage)
    contexts=copy.deepcopy(controller.profile['contexts']);by_name={c['name']:c for c in contexts}
    artifacts=dict(controller.profile['component_artifacts'])
    if coverage['kernel_requests']:
        result,path=_run_worker('gpu',coverage['kernel_requests'],root,controller,worker_timeout)
        expected={(r['device'],*r['shape']) for r in coverage['kernel_requests']}
        actual={(r['device'],*[r[key] for key in ['N','B','K','C']]) for r in result['rows']}
        if actual!=expected or len(result['rows'])!=len(expected):raise ValueError('Incomplete or duplicate kernel capture')
        captured={(r['device'],r['N'],r['B'],r['K'],r['C']):r for r in result['rows']}
        for request in coverage['kernel_requests']:
            row=copy.deepcopy(captured[request['device'],*request['shape']]);row.pop('device')
            if set(row)!={'N','B','K','C','validate_range','kernels'} or row['validate_range'] is not True:
                raise ValueError('Unexpected kernel geometry fields')
            sanitized=kernel_census([dict(cat='kernel',name=r['name'],args=r['geometry']) for r in row['kernels']])
            if row['kernels']!=sanitized:raise ValueError('Only duration-free kernel fields may be merged')
            for name in request['contexts']:
                by_name[name]['profiles'][request['device']].setdefault('kernel_geometry',[]).append(copy.deepcopy(row))
        artifacts[str(path.resolve())]=sha256_file(path)
    if coverage['host_requests']:
        result,path=_run_worker('host',coverage['host_requests'],root,controller,worker_timeout)
        validate_host_capture(result)
        actual={(kind,int(size)) for kind,bank in result['arrays'].items() for size in bank}
        if actual!={(r['kind'],r['bytes']) for r in coverage['host_requests']}:raise ValueError('Incomplete allocator capture')
        for request in coverage['host_requests']:
            for name,device in request['targets']:
                inputs=by_name[name]['profiles'][device]['pageable_host_service']
                validate_pageable_geometry(result,inputs['prices'])
                geometry=inputs['geometry'];geometry['arrays'].setdefault(request['kind'],{})[str(request['bytes'])]=copy.deepcopy(result['arrays'][request['kind']][str(request['bytes'])])
                geometry['source_sha256'].update(result['source_sha256'])
                geometry['scope']='Merged duration-free allocator-route censuses; existing extents retained, newly required extents collected in a fresh worker. Future allocator history remains conditional.'
        artifacts[str(path.resolve())]=sha256_file(path)
    if _price_contexts(contexts)!=_price_contexts(controller.profile['contexts']):raise ValueError('Geometry collection changed component prices')
    # Require full coverage through the same planner. A captured but unrecognized
    # compiled path is an explicit failure, never a partial-subset recommendation.
    plan=bounded_trait_plan(workload,contexts,bounds=controller.config['bounds'],joint=controller.config['joint'],output=output)
    if input_identity(controller.input_path)!=before_input:raise ValueError('Input changed during geometry collection')
    current=execution_context(controller.devices,input_path=controller.input_path,output_path=controller.output_path)
    validate_detailed_profile(controller.profile,current)
    write_record(root/'plan.json',plan)
    write_record(root/'request.json',dict(workload=workload,config=controller.config,output=output,output_path=str(output_path)))
    artifacts[str((root/'coverage.json').resolve())]=sha256_file(root/'coverage.json')
    complete=bind_detailed_profile(contexts,current,component_artifacts=artifacts,
        price_bindings=controller.profile.get('price_bindings'),
        limitations=controller.profile['limitations']+[
            'Missing feasible kernel/allocator geometry collected without retaining durations; component prices unchanged.',
            'Fresh-worker capture avoids changing the scan process allocator history. Concurrent allocator routes remain conditional.',
            'Complete geometry coverage does not certify absolute runtime or selection accuracy for this new workload.'])
    write_detailed_profile(complete,root/'profile.json')
    manifest=dict(schema=SCHEMA,profile_sha256=sha256_file(root/'profile.json'),plan_sha256=sha256_file(root/'plan.json'),
        kernel_shapes_collected=len(coverage['kernel_requests']),host_extents_collected=len(coverage['host_requests']),
        component_prices_unchanged=True,durations_retained=False,input_file_identity=before_input,
        memory_feasible_candidates=len(coverage['memory_feasible_candidates']),selected=plan['selected']['candidate_index'],
        source_sha256=source_identity())
    write_record(root/'manifest.json',manifest)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['inspect','complete','worker'])
    p.add_argument('--profile');p.add_argument('--workload');p.add_argument('--config');p.add_argument('--output-settings')
    p.add_argument('--output-path');p.add_argument('--collection-dir');p.add_argument('--request');p.add_argument('--result')
    for name,default in [('max-kernel-shapes',128),('max-host-extents',128),('max-device-probe-bytes',2<<30),
                         ('max-host-probe-bytes',1<<30),('max-host-touch-bytes',8<<30),('worker-timeout',300)]:
        p.add_argument('--'+name,type=int,default=default)
    a=p.parse_args()
    if a.mode=='worker':_worker(a.request,a.result);return
    if not all([a.profile,a.workload,a.config,a.output_settings,a.output_path,a.collection_dir]):
        p.error('profile, workload, config, output-settings, output-path and collection-dir are required')
    result=complete_profile_geometry(a.profile,json.loads(Path(a.workload).read_text()),json.loads(Path(a.config).read_text()),
        output=json.loads(Path(a.output_settings).read_text()),output_path=a.output_path,collection_dir=a.collection_dir,
        collect=a.mode=='complete',worker_timeout=a.worker_timeout,**{name:getattr(a,name) for name in
        ['max_kernel_shapes','max_host_extents','max_device_probe_bytes','max_host_probe_bytes','max_host_touch_bytes']})
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
