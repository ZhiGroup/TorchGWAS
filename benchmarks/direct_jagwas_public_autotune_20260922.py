"""Public JAGWAS binding/cache/output audit in fresh processes.

Synthetic service prices and retained duration-free kernel captures exercise
wiring only. This does not qualify throughput, ranking or empirical capacities.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
import torch
from test_jagwas_candidate_space import specification
from direct_jagwas_bounded_execution_20260922 import reference,read_result
from torchgwas.api import run_linear_gwas
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import (bind_detailed_profile,execution_context,source_identity,
                                          sha256_file,write_detailed_profile)
from torchgwas.geometry_collection import write_record


def configure():
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False


def prepare(fixture,out):
    out.mkdir(parents=True,exist_ok=False)
    spec=specification(fixture/'input.pgen');sources=source_identity()
    devices=['cuda:0','cuda:1'];limit=2<<30;memory={}
    for device in devices:
        p=torch.cuda.get_device_properties(device)
        memory[device]=dict(torch_version=torch.__version__,sm_count=p.multi_processor_count,
            max_threads_per_sm=p.max_threads_per_multi_processor,compute_capability=[p.major,p.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
    prepared={}
    for case,count in [('all',None),('single',1),('dual',2)]:
        directory=out/case;directory.mkdir()
        contexts=copy.deepcopy([c for c in spec['contexts'] if count is None or len(c['devices'])==count])
        active=list(dict.fromkeys(d for c in contexts for d in c['devices']))
        services={c['name']:copy.deepcopy(spec['preparation_services'][c['name']]) for c in contexts}
        current=execution_context(active,input_path=fixture/'input.pgen',output_path=directory)
        prices=dict(writer_prices=spec['prices'],preparation_services=services)
        record=CalibrationParameterCache(out/'measurements').store('cpu_capacity','jagwas_components',prices,
            dependencies=dict(source_sha256=sources,execution_context=current),
            provenance={'scope':'Synthetic resource controls, not measured capacities.'},
            observed_unix_seconds=time.time(),max_age_seconds=3600.)
        capture=Path('tests/fixtures/jagwas_chunk_geometry.json').resolve()
        profile=bind_detailed_profile(contexts,current,component_artifacts={
            str(Path(record['path']).resolve()):sha256_file(record['path']),str(capture):sha256_file(capture)},
            limitations=['Synthetic component prices and retained compiled geometry: identity/execution audit only.'],
            sources=sources)
        write_detailed_profile(profile,directory/'profile.json')
        config=dict(bounds=spec['bounds'],qc_trait_block=128,jagwas_services=record['path'],
            plan_cache_dir=str(directory/'plan_cache'),joint=dict(spec['joint'],
                host_memory_bytes=4<<30,host_reserve_bytes=512<<20,device_reserve_bytes=256<<20,
                device_memory_bytes={d:limit for d in active},device_memory_profiles={d:memory[d] for d in active}))
        write_record(directory/'config.json',config)
        prepared[case]=dict(calibration_record_sha256=record['record_sha256'],
            artifact_sha256=sha256_file(record['path']),devices=active)
    rng=np.random.default_rng(9221703)
    calls=rng.integers(0,3,size=(1025,2049),dtype=np.uint8)
    y=np.load(fixture/'phenotype.npy');cov=np.load(fixture/'covariates.npy')
    # Separate FP64 OLS with an intercept for sampled variants, including tails.
    indices=np.unique(np.r_[np.arange(4),np.random.default_rng(935).choice(1025,29,replace=False),1024])
    truth=reference(calls,y,cov,indices);np.savez(out/'reference.npz',indices=indices,chi2=truth)
    write_record(out/'prepared.json',dict(source_sha256=sources,cases=prepared,budget_bytes=limit,
        inputs={str(fixture/name):sha256_file(fixture/name)
                for name in ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']},
        source_mount=json.loads(subprocess.check_output(
            ['findmnt','-J','-T',str(fixture/'input.pgen'),'-o','SOURCE,FSTYPE,TARGET'],text=True))))
    print(json.dumps(dict(phase='prepared',dimensions=[2049,1025,512,2])),flush=True)


def execute(fixture,out,case):
    prepared=json.loads((out/'prepared.json').read_text())
    binding=case if case in ('single','dual') else 'all'
    root=out/binding
    config=json.loads((root/'config.json').read_text())
    profile=json.loads((root/'profile.json').read_text())
    assert source_identity()==prepared['source_sha256']
    for device,limit in config['joint']['device_memory_bytes'].items():
        assert torch.cuda.mem_get_info(device)[0]>limit
        torch.cuda.set_per_process_memory_fraction(limit/torch.cuda.get_device_properties(device).total_memory,device)
    directory=root/('execution_'+case)
    if case=='expired':
        old=json.loads(Path(config['jagwas_services']).read_text())
        record=CalibrationParameterCache(out/'expired_measurements').store('cpu_capacity','jagwas_components',old['value'],
            dependencies=old['dependencies'],provenance={'scope':'Expired synthetic control'},
            observed_unix_seconds=time.time()-100.,max_age_seconds=1.)
        config['jagwas_services']=record['path']
        profile['component_artifacts'][str(Path(record['path']).resolve())]=sha256_file(record['path'])
        try:
            run_linear_gwas(fixture/'input.pgen',fixture/'phenotype.npy',fixture/'covariates.npy',
                pgen_mode='hardcall',reduce='jagwas',autotune_profile=profile,autotune_config=config,output_dir=directory)
        except ValueError as error:assert 'Expired bound calibration' in str(error)
        else:raise AssertionError('Expired calibration accepted')
        assert not directory.exists()
        write_record(out/'expired_check.json',dict(rejected_before_output=True));return
    result=run_linear_gwas(fixture/'input.pgen',fixture/'phenotype.npy',fixture/'covariates.npy',
        pgen_mode='hardcall',reduce='jagwas',sumstats_queue_depth=2,autotune_profile=profile,
        autotune_config=config,output_dir=directory)
    audit=result.run_metadata['autotune']
    assert audit['analytical_cache']['status']==('hit' if case=='second' else 'miss')
    assert audit['reduction']=='jagwas' and result.run_metadata['trait_block'] is None
    assert audit['output']['store_beta'] is False
    assert audit['reduction_calibration']['record_sha256']==prepared['cases'][binding]['calibration_record_sha256']
    assert sha256_file(config['jagwas_services'])==prepared['cases'][binding]['artifact_sha256']
    if case in ('single','dual'):
        assert len(audit['selected']['devices'])==(1 if case=='single' else 2)
    manifest,values=read_result(directory,1025)
    assert manifest['shape']==[1025,512] and manifest['df']==512 and manifest['rows']==1025
    with np.load(out/'reference.npz') as truth:
        actual=values[truth['indices']]
        np.testing.assert_allclose(actual,truth['chi2'],rtol=3e-4,atol=3e-4)
        max_error=float(np.max(np.abs(actual-truth['chi2'])))
    np.save(out/(case+'_values.npy'),values)
    check=dict(case=case,pid=os.getpid(),cache_status=audit['analytical_cache']['status'],
        cache_key=audit['analytical_cache']['key'],selected_candidate=audit['selected']['candidate_index'],
        chunk_size=audit['selected']['chunk_size'],devices=audit['selected']['devices'],retained_variants=1025,
        calibration_record_sha256=audit['reduction_calibration']['record_sha256'],
        observed_unix_seconds=audit['reduction_calibration']['observed_unix_seconds'],
        calibration_age_seconds=audit['reduction_calibration']['age_seconds'],
        max_absolute_reference_error=max_error,output_sha256=hashlib.sha256(values.tobytes()).hexdigest())
    write_record(out/(case+'_check.json'),check);print(json.dumps(check),flush=True)


def finish(fixture,out):
    prepared=json.loads((out/'prepared.json').read_text())
    for path,digest in prepared['inputs'].items():assert sha256_file(path)==digest
    runs=[json.loads((out/(case+'_check.json')).read_text()) for case in ['first','second','single','dual']]
    assert len({row['pid'] for row in runs})==4
    for key in ['cache_key','selected_candidate','output_sha256','calibration_record_sha256','observed_unix_seconds']:
        assert runs[0][key]==runs[1][key]
    assert runs[1]['calibration_age_seconds']>=runs[0]['calibration_age_seconds']
    baseline=np.load(out/'first_values.npy')
    for case in ['second','single','dual']:
        np.testing.assert_allclose(np.load(out/(case+'_values.npy')),baseline,rtol=6e-5,atol=3e-4)
    expired=json.loads((out/'expired_check.json').read_text());assert expired['rejected_before_output']
    assert source_identity()==prepared['source_sha256']
    write_record(out/'report.json',dict(preparation=prepared,runs=runs,expired_control=expired,
        benchmark_sha256=sha256_file(__file__),
        scope='Four fresh-process public full-panel JAGWAS scans on one/two GPUs; immutable component and plan reuse. Synthetic prices, numerical and binding validation only; no throughput or ranking qualification.'))
    print(json.dumps(dict(phase='complete',fresh_processes=4,retained_variants=1025,expired_rejected=True)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--phase',choices=['prepare','first','second','single','dual','expired','finish'],required=True)
    args=parser.parse_args();configure();fixture=Path(args.fixture);out=Path(args.out)
    if args.phase=='prepare':prepare(fixture,out)
    elif args.phase=='finish':finish(fixture,out)
    else:execute(fixture,out,args.phase)
