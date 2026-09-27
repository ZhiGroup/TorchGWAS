"""Independent tiny setup-phase prices for explicit covariate ranks.

No genotype or association execution. One trait at N=max(32,rank+3), seven
complete repeats of fixed generic operations; exact references, no interpolation.
"""
from __future__ import annotations
import argparse
import copy
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

from .detailed_calibration import (bind_detailed_profile,execution_context,read_detailed_profile,
    sha256_file,source_identity,validate_detailed_profile,write_detailed_profile)
from .geometry_collection import write_record
from .setup_work import setup_reference_shape,setup_primitive_bank

SCHEMA='torchgwas.setup_primitives.v1'
PHASES=['residual_common','residual_block','design_common','design_block']


def primitive_requests(contexts,ranks,*,max_ranks=8,max_covariates=256):
    if any(type(v) is not int or v<1 for v in [max_ranks,max_covariates]):raise ValueError('Positive calibration limits required')
    if not isinstance(ranks,(list,tuple)) or not ranks or len(ranks)>max_ranks or len(set(ranks))!=len(ranks):
        raise ValueError('Explicit unique bounded covariate ranks required')
    for rank in ranks:
        if type(rank) is not int or not 0<=rank<=max_covariates:raise ValueError('Covariate rank exceeds collection bound')
    requests={}
    for context in contexts:
        for device,profile in context['profiles'].items():
            for rank in ranks:
                explicit=profile.get('setup_primitives_by_rank',{}).get(str(rank))
                legacy=profile.get('setup_primitives',{})
                present=explicit is not None or any(row.get('reference_shape')==setup_reference_shape(rank) for row in legacy.values())
                if present:
                    setup_primitive_bank(profile,rank);continue
                row=requests.setdefault((device,rank),dict(device=device,rank=rank,contexts=[]))
                row['contexts'].append(context['name'])
    return list(requests.values())


def _measure(device,rank):
    import numpy as np
    import torch
    from .preprocess import _covariate_basis
    torch.cuda.set_device(device)
    # Fixed tiny reference and allocator cap, independent of workload N/K/B.
    if torch.cuda.mem_get_info(device)[0]<128<<20:raise ValueError('Insufficient live memory for tiny setup probe')
    torch.cuda.set_per_process_memory_fraction((128<<20)/torch.cuda.get_device_properties(device).total_memory,device)
    n,t,c=setup_reference_shape(rank)
    rng=np.random.default_rng(816);ph=rng.standard_normal((n,1),dtype=np.float32)
    basis=None if not c else _covariate_basis(rng.standard_normal((n,c),dtype=np.float32))
    if c and (basis is None or basis.shape[1]!=c):raise ValueError('Tiny probe basis rank differs')
    def residual_common():
        return None if basis is None else torch.as_tensor(np.ascontiguousarray(basis),dtype=torch.float32,device=device)
    q=residual_common()
    def residual_block():
        values=torch.as_tensor(np.ascontiguousarray(ph),dtype=torch.float32,device=device)
        values-=values.mean(dim=0,keepdim=True)
        if q is not None:values-=q@(q.T@values)
        sd=values.std(dim=0,keepdim=True,correction=0)
        sd=torch.where(sd==0,torch.ones_like(sd),sd);values/=sd
        return values.cpu().numpy()
    def design_common():
        intercept=torch.full((n,1),1/np.sqrt(n),device=device)
        cv=intercept if basis is None else torch.cat((intercept,torch.as_tensor(basis,dtype=torch.float32,device=device)),dim=1)
        matrix=torch.empty((n,1+c+1),dtype=torch.float32,device=device)
        matrix[:,1:].copy_(cv)
        return matrix,torch.empty(1,dtype=torch.float32,device=device)
    matrix,ss=design_common();columns=matrix[:,:1]
    def design_block():
        columns.copy_(torch.as_tensor(np.ascontiguousarray(ph),dtype=torch.float32))
        return torch.sum(columns*columns,dim=0,out=ss)
    rows=[]
    for name,fn in zip(PHASES,[residual_common,residual_block,design_common,design_block]):
        if name=='residual_common' and not c:
            rows.extend(dict(phase=name,repeat=i,loops=100,reference_shape=[n,1,c],measured=False,
                cpu_seconds=0.,non_cpu_seconds=0.,wall_seconds=0.) for i in range(7))
            continue
        for _ in range(5):value=fn()
        torch.cuda.synchronize(device)
        for repeat in range(7):
            loops=100;torch.cuda.synchronize(device);wall=time.perf_counter();cpu=time.thread_time()
            for _ in range(loops):value=fn()
            cpu=time.thread_time()-cpu;torch.cuda.synchronize(device);wall=time.perf_counter()-wall
            rows.append(dict(phase=name,repeat=repeat,loops=loops,reference_shape=[n,1,c],measured=True,
                cpu_seconds=cpu/loops,non_cpu_seconds=max(0.,wall-cpu)/loops,wall_seconds=wall/loops))
    phases={name:dict(reference_shape=[n,1,c],cpu_seconds=statistics.median(r['cpu_seconds'] for r in rows if r['phase']==name),
        non_cpu_seconds=statistics.median(r['non_cpu_seconds'] for r in rows if r['phase']==name)) for name in PHASES}
    return dict(device=device,rank=rank,rows=rows,phases=phases,
        scope='Isolated worker, warm tiny generic setup operations. CPU and non-CPU components retained separately; no genotype, association or workload-shape timing. Concurrent ownership and cold first use remain separate model terms.')


def validate_measurement(record):
    from .first_principles import positive
    rank=record['rank'];expected=setup_reference_shape(rank);rows=record['rows']
    if len(rows)!=28 or {(r['phase'],r['repeat']) for r in rows}!={(p,i) for p in PHASES for i in range(7)}:
        raise ValueError('Seven complete setup repetitions required')
    for row in rows:
        if row['loops']!=100 or row['reference_shape']!=expected:raise ValueError('Tiny setup reference differs')
        if row.get('measured') is not (rank!=0 or row['phase']!='residual_common'):raise ValueError('Incorrect structural-zero declaration')
        for key in ['cpu_seconds','non_cpu_seconds','wall_seconds']:positive(key,row[key],True)
        if abs(row['non_cpu_seconds']-max(0.,row['wall_seconds']-row['cpu_seconds']))>1e-12:
            raise ValueError('Setup wall/CPU accounting differs')
    bank=setup_primitive_bank({'setup_primitives':record['phases']},rank)
    for name,value in bank.items():
        for key in ['cpu_seconds','non_cpu_seconds']:
            if value[key]!=statistics.median(r[key] for r in rows if r['phase']==name):raise ValueError('Setup price does not retain full median')
    if rank==0 and any(bank['residual_common'][key] for key in ['cpu_seconds','non_cpu_seconds']):
        raise ValueError('Rank zero must not charge a basis upload')
    return bank


def worker(request_path,result_path):
    import torch
    request=json.loads(Path(request_path).read_text());context=request['execution_context']
    if request['schema']!=SCHEMA or request['source_sha256']!=source_identity():raise ValueError('Setup worker source differs')
    torch.set_num_threads(context['torch_threads']);torch.set_num_interop_threads(context['torch_interop_threads'])
    torch.backends.cuda.matmul.allow_tf32=context['allow_tf32']
    current=execution_context(list(context['devices']),input_path=request['input_path'],output_path=request['output_path'])
    if current!=context:raise ValueError('Setup worker execution context differs')
    with ThreadPoolExecutor(max_workers=1) as pool:
        rows=[pool.submit(_measure,r['device'],r['rank']).result() for r in request['requests']]
    if request['source_sha256']!=source_identity():raise ValueError('Source changed during setup measurement')
    for row in rows:validate_measurement(row)
    write_record(result_path,dict(schema=SCHEMA,source_sha256=source_identity(),execution_context=current,
        request_sha256=sha256_file(request_path),measurements=rows))


def complete_setup_primitives(profile,ranks,*,input_path,output_path,collection_dir,
                              max_ranks=8,max_covariates=256,worker_timeout=300):
    if type(worker_timeout) is not int or worker_timeout<1:raise ValueError('Positive worker timeout required')
    original=read_detailed_profile(profile) if isinstance(profile,(str,Path)) else copy.deepcopy(profile)
    devices=list(dict.fromkeys(d for c in original['contexts'] for d in c['devices']))
    current=execution_context(devices,input_path=input_path,output_path=output_path)
    validate_detailed_profile(original,current)
    requests=primitive_requests(original['contexts'],ranks,max_ranks=max_ranks,max_covariates=max_covariates)
    root=Path(collection_dir);root.mkdir(parents=True,exist_ok=False)
    request_path=root/'request.json';result_path=root/'primitives.json'
    write_record(request_path,dict(schema=SCHEMA,requests=requests,source_sha256=source_identity(),
        execution_context=current,input_path=str(input_path),output_path=str(output_path)))
    contexts=copy.deepcopy(original['contexts']);by_name={c['name']:c for c in contexts}
    artifacts=dict(original['component_artifacts'])
    if requests:
        with (root/'worker.log').open('x') as log:
            subprocess.run([sys.executable,'-m','torchgwas.setup_calibration','worker','--request',str(request_path),'--result',str(result_path)],
                stdout=log,stderr=subprocess.STDOUT,check=True,timeout=worker_timeout)
        result=json.loads(result_path.read_text())
        if (result['schema']!=SCHEMA or result['source_sha256']!=source_identity() or result['execution_context']!=current
            or result['request_sha256']!=sha256_file(request_path)):raise ValueError('Setup primitive provenance differs')
        rows=result['measurements'];identities={(r['device'],r['rank']) for r in rows}
        if len(rows)!=len(requests) or identities!={(r['device'],r['rank']) for r in requests}:
            raise ValueError('Incomplete setup reference coverage')
        prices={(r['device'],r['rank']):validate_measurement(r) for r in rows}
        for request in requests:
            for name in request['contexts']:
                by_name[name]['profiles'][request['device']].setdefault('setup_primitives_by_rank',{})[str(request['rank'])]=copy.deepcopy(prices[request['device'],request['rank']])
        artifacts[str(result_path.resolve())]=sha256_file(result_path)
    for context in contexts:
        for p in context['profiles'].values():
            for rank in ranks:setup_primitive_bank(p,rank)
    # Assert the only change is insertion of the requested setup references.
    restored=copy.deepcopy(contexts)
    for request in requests:
        for name in request['contexts']:
            p=next(c for c in restored if c['name']==name)['profiles'][request['device']]
            del p['setup_primitives_by_rank'][str(request['rank'])]
            old=next(c for c in original['contexts'] if c['name']==name)['profiles'][request['device']]
            if 'setup_primitives_by_rank' not in old and not p['setup_primitives_by_rank']:del p['setup_primitives_by_rank']
    if restored!=original['contexts']:raise ValueError('Unrelated component prices changed')
    current=execution_context(devices,input_path=input_path,output_path=output_path);validate_detailed_profile(original,current)
    artifacts[str(request_path.resolve())]=sha256_file(request_path)
    completed=bind_detailed_profile(contexts,current,component_artifacts=artifacts,
        price_bindings=original.get('price_bindings'),limitations=original['limitations']+[
        'Additional exact-rank tiny setup primitives collected in an isolated worker. No GWAS timing, fitted correction or interpolation.',
        'Setup bulk GEMM efficiency, cold library initialization and concurrent allocation behavior remain unresolved.'])
    write_detailed_profile(completed,root/'profile.json')
    manifest=dict(schema=SCHEMA,profile_sha256=sha256_file(root/'profile.json'),ranks=list(ranks),
        device_rank_references_collected=len(requests),other_component_prices_unchanged=True,source_sha256=source_identity())
    write_record(root/'manifest.json',manifest)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['complete','worker'])
    for key in ['profile','input-path','output-path','collection-dir','request','result']:p.add_argument('--'+key)
    p.add_argument('--ranks',nargs='+',type=int);p.add_argument('--max-ranks',type=int,default=8)
    p.add_argument('--max-covariates',type=int,default=256);p.add_argument('--worker-timeout',type=int,default=300)
    a=p.parse_args()
    if a.mode=='worker':worker(a.request,a.result);return
    if not all([a.profile,a.input_path,a.output_path,a.collection_dir,a.ranks]):p.error('profile, paths, collection-dir and ranks required')
    result=complete_setup_primitives(a.profile,a.ranks,input_path=a.input_path,output_path=a.output_path,
        collection_dir=a.collection_dir,max_ranks=a.max_ranks,max_covariates=a.max_covariates,worker_timeout=a.worker_timeout)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
