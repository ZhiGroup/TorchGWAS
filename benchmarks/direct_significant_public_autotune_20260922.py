"""Public significant autotune across fresh processes with immutable prices.

Synthetic resource prices exercise binding/admission/cache/execution only.
This is not throughput, ranking, or independent capacity qualification.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
import torch
from direct_significant_bounded_execution_20260922 import synthetic_profile,reference,read_result
from test_significant_host_model import bank
from torchgwas.api import run_linear_gwas
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import (bind_detailed_profile,execution_context,source_identity,
                                          sha256_file,write_detailed_profile)
from torchgwas.geometry_collection import _gpu_geometry,write_record
from torchgwas.significant_candidate_space import prepare_significant_host_candidates
from torchgwas.significant_host_work import significant_host_memory
from torchgwas.tensor_service import host_primitive_name
from torchgwas.tensor_work import eager_statistics_work


def configure():
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False


def prepare(fixture,out):
    out.mkdir(parents=True,exist_ok=False)
    n,m,k,c=2049,1025,4097,2;devices=['cuda:0','cuda:1'];contexts=[]
    for count in (1,2):
        active=devices[:count]
        contexts.append(dict(name=f'gpus{count}',devices=active,
            profiles={d:synthetic_profile(n,c,2//count,d) for d in active},
            shared_capacities=dict(cpu=2.,dram=1e9,input=1e8,output=1e8)))
    memory_profiles={}
    for device in devices:
        p=torch.cuda.get_device_properties(device)
        memory_profiles[device]=dict(torch_version=torch.__version__,sm_count=p.multi_processor_count,
            max_threads_per_sm=p.max_threads_per_multi_processor,compute_capability=[p.major,p.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
    workload=dict(genotype=str(fixture/'input.pgen'),samples=n,markers=m,traits=k,covariates=c,
        matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
    output=dict(block_bytes=None,queue_depth=1,store_beta=True,fsync=True)
    bounds=dict(chunks=[128,256],trait_blocks=[512,1024,k],max_candidates=12)
    space=prepare_significant_host_candidates(workload,contexts,output=output,**bounds)
    reserve=64<<20
    memories=[significant_host_memory(row,device_reserve_bytes=reserve,device_memory_profiles=memory_profiles)
              for row in space['candidates']]
    small=max(max(memory['device_bytes'].values()) for row,memory in zip(space['candidates'],memories) if row['trait_block']!=k)
    full=min(max(memory['device_bytes'].values()) for row,memory in zip(space['candidates'],memories) if row['trait_block']==k)
    assert small<full;limit=(small+full)//2;requests={}
    for candidate,memory in zip(space['candidates'],memories):
        if max(memory['device_bytes'].values())>limit:continue
        for tile in candidate['tiles']:
            width=tile['trait_range'][1]-tile['trait_range'][0];device=tile['device'];chunk=tile['profile']['chunk_markers']
            for b in sorted({chunk,m%chunk}):
                requests[(device,n,b,width,c)]=dict(device=device,shape=[n,b,width,c],
                    source_memory_bound_bytes=memory['device_bytes'][device],device_budget_bytes=limit)
    print(json.dumps(dict(phase='geometry',shapes=len(requests),budget_bytes=limit)),flush=True)
    captured=_gpu_geometry(list(requests.values()),out);write_record(out/'geometry.json',captured)
    for context in contexts:
        for device,profile in context['profiles'].items():
            rows=[{key:value for key,value in row.items() if key!='device'} for row in captured['rows'] if row['device']==device]
            profile['kernel_geometry']=copy.deepcopy(rows)
            for row in rows:
                work=eager_statistics_work(n,row['B'],row['K'],c,True)
                profile['host_primitives'].update({host_primitive_name(call):1e-6 for call in work['host_calls']})
    print(json.dumps(dict(phase='reference')),flush=True)
    calls=np.random.default_rng(9222207).integers(0,3,size=(m,n),dtype=np.uint8)
    y=np.load(fixture/'phenotype.npy',mmap_mode='r');cov=np.load(fixture/'covariates.npy')
    truth,margin=reference(calls,y,cov,.001);np.savez(out/'reference.npz',**truth)
    assert len(truth['variant_index'])==4220
    current=execution_context(devices,input_path=fixture/'input.pgen',output_path=out)
    sources=source_identity();cache=CalibrationParameterCache(out/'measurements')
    record=cache.store('cpu_capacity','significant_host_components',bank(),
        dependencies=dict(source_sha256=sources,execution_context=current),
        provenance={'scope':'Synthetic controls for public API correctness only; not measured service rates.'},
        observed_unix_seconds=time.time(),max_age_seconds=3600.)
    profile=bind_detailed_profile(contexts,current,component_artifacts={record['path']:sha256_file(record['path']),
        str((out/'geometry.json').resolve()):sha256_file(out/'geometry.json')},
        limitations=['Synthetic component prices: this binding validates identity and execution wiring only.'],sources=sources)
    write_detailed_profile(profile,out/'profile.json')
    config=dict(bounds=bounds,qc_trait_block=128,significant_host_prices=record['path'],plan_cache_dir=str(out/'plan_cache'),
        joint=dict(occupancy_scenarios={'none':'empty','all':'dense'},host_scenarios={'fluid':dict(host_serial_fraction=0.)},
            cpu_workers=2,host_memory_bytes=4<<30,host_reserve_bytes=512<<20,device_reserve_bytes=reserve,
            device_memory_bytes={d:limit for d in devices},device_memory_profiles=memory_profiles))
    write_record(out/'config.json',config)
    write_record(out/'prepared.json',dict(source_sha256=sources,calibration_record_sha256=record['record_sha256'],
        artifact_sha256=sha256_file(record['path']),budget_bytes=limit,largest_tiled_bound_bytes=small,
        smallest_full_panel_bound_bytes=full,reference_threshold_margin=margin,
        inputs={str(fixture/name):sha256_file(fixture/name) for name in ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}))
    print(json.dumps(dict(phase='prepared',calibration_record=record['record_sha256'])),flush=True)


def execute(fixture,out,case):
    prepared=json.loads((out/'prepared.json').read_text());config=json.loads((out/'config.json').read_text())
    assert source_identity()==prepared['source_sha256']
    for device,limit in config['joint']['device_memory_bytes'].items():
        assert torch.cuda.mem_get_info(device)[0]>limit
        torch.cuda.set_per_process_memory_fraction(limit/torch.cuda.get_device_properties(device).total_memory,device)
    profile=json.loads((out/'profile.json').read_text())
    directory=out/case
    if case=='expired':
        # A separately published expired control must fail before output exists.
        record=CalibrationParameterCache(out/'expired_measurements').store('cpu_capacity','significant_host_components',bank(),
            dependencies=dict(source_sha256=profile['source_sha256'],execution_context=profile['execution_context']),
            provenance={'scope':'Deliberately expired synthetic test control'},observed_unix_seconds=time.time()-100.,max_age_seconds=1.)
        profile['component_artifacts'][str(Path(record['path']).resolve())]=sha256_file(record['path'])
        config['significant_host_prices']=record['path']
        try:
            run_linear_gwas(fixture/'input.pgen',fixture/'phenotype.npy',fixture/'covariates.npy',
                pgen_mode='hardcall',reduce='significant',significance_threshold=.001,
                autotune_profile=profile,autotune_config=config,output_dir=directory)
        except ValueError as error:
            assert 'Expired bound calibration' in str(error)
        else:raise AssertionError('Expired evidence was accepted')
        assert not directory.exists()
        write_record(out/'expired_check.json',dict(rejected_before_output=True));return
    result=run_linear_gwas(fixture/'input.pgen',fixture/'phenotype.npy',fixture/'covariates.npy',
        pgen_mode='hardcall',reduce='significant',significance_threshold=.001,sumstats_fields='beta+t',
        sumstats_queue_depth=1,autotune_profile=out/'profile.json',autotune_config=out/'config.json',output_dir=directory)
    audit=result.run_metadata['autotune']
    assert audit['analytical_cache']['status']==('miss' if case=='first' else 'hit')
    assert audit['reduction']=='significant' and audit['significance_threshold']==.001
    assert audit['reduction_calibration']['record_sha256']==prepared['calibration_record_sha256']
    assert sha256_file(config['significant_host_prices'])==prepared['artifact_sha256']
    manifest,values=read_result(directory,1025,4097)
    with np.load(out/'reference.npz') as truth:
        for field in ['variant_index','trait_index','df']:np.testing.assert_array_equal(values[field],truth[field])
        for field in ['beta','t_stat']:np.testing.assert_allclose(values[field],truth[field],rtol=3e-5,atol=3e-5)
        max_error=float(np.max(np.abs(values['t_stat']-truth['t_stat'])))
    assert manifest['rows']==4220 and audit['selected']['trait_block']<4097
    assert audit['candidates_feasible']==8 and len(audit['rejected'])==2
    check=dict(case=case,pid=os.getpid(),cache_status=audit['analytical_cache']['status'],cache_key=audit['analytical_cache']['key'],
        selected_candidate=audit['selected']['candidate_index'],chunk_size=audit['selected']['chunk_size'],
        trait_block=audit['selected']['trait_block'],devices=audit['selected']['devices'],retained_pairs=manifest['rows'],
        calibration_record_sha256=audit['reduction_calibration']['record_sha256'],
        observed_unix_seconds=audit['reduction_calibration']['observed_unix_seconds'],
        calibration_age_seconds=audit['reduction_calibration']['age_seconds'],max_absolute_t_error=max_error,
        pair_sha256=hashlib.sha256(np.column_stack([values['variant_index'],values['trait_index']]).tobytes()).hexdigest())
    write_record(out/(case+'_check.json'),check);print(json.dumps(check),flush=True)


def finish(fixture,out):
    prepared=json.loads((out/'prepared.json').read_text())
    for path,digest in prepared['inputs'].items():assert sha256_file(path)==digest
    first=json.loads((out/'first_check.json').read_text());second=json.loads((out/'second_check.json').read_text())
    assert first['pid']!=second['pid'] and first['cache_key']==second['cache_key']
    for field in ['selected_candidate','pair_sha256','calibration_record_sha256','observed_unix_seconds']:
        assert first[field]==second[field]
    assert second['calibration_age_seconds']>=first['calibration_age_seconds']
    expired=json.loads((out/'expired_check.json').read_text());assert expired['rejected_before_output']
    assert source_identity()==prepared['source_sha256']
    revision=json.loads((out/'harness_revision.json').read_text()) if (out/'harness_revision.json').exists() else None
    write_record(out/'report.json',dict(preparation=prepared,runs=[first,second],expired_control=expired,harness_revision=revision,
        benchmark_sha256=sha256_file(__file__),scope='Fresh-process public significant autotune, indexed nonempty output, bounded tiles, immutable calibration and disk-plan reuse. Synthetic resource prices; no ranking or throughput qualification.'))
    print(json.dumps(dict(phase='complete',fresh_processes=2,retained_pairs=4220,expired_rejected=True)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--phase',choices=['prepare','first','second','expired','finish'],required=True)
    args=parser.parse_args();configure();fixture=Path(args.fixture);out=Path(args.out)
    if args.phase=='prepare':prepare(fixture,out)
    elif args.phase=='finish':finish(fixture,out)
    else:execute(fixture,out,args.phase)
