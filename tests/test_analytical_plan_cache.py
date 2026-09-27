import copy
import hashlib
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from torchgwas.analytical_plan_cache import AnalyticalPlanCache,canonical,input_identity
from torchgwas.detailed_autotune import DetailedAutotune
from test_detailed_autotune import fixture,source


def cached_fixture(tmp_path,monkeypatch):
    import torch
    import torchgwas.detailed_autotune as module
    profile,config,path,context=fixture(tmp_path,monkeypatch)
    config['plan_cache_dir']=str(tmp_path/'cache')
    monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(1<<30,1<<30))
    monkeypatch.setattr(torch.cuda,'memory_reserved',lambda d:0)
    monkeypatch.setattr(torch.cuda,'memory_allocated',lambda d:0)
    monkeypatch.setattr(torch.backends.cuda.matmul,'allow_tf32',False)
    monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:1<<30)
    calls=[]
    def build(workload,contexts,*,bounds,joint,output):
        calls.append(copy.deepcopy((workload,contexts,bounds,joint,output)))
        identity=input_identity(path)
        selected=dict(candidate_index=0,worst_supplied_scenario_seconds=10.,
            devices=['cuda:1'],required_environment={},
            memory=dict(device_bytes={'cuda:1':128},host_bytes=256),
            api_kwargs=dict(chunk_size=bounds['chunks'][0],trait_block=bounds['trait_blocks'][0],
                trait_devices=['cuda:1'],reader_workers=4,prefetch_chunks=4))
        alternate=copy.deepcopy(selected);alternate['candidate_index']=1
        alternate['worst_supplied_scenario_seconds']=20.
        alternate['memory']=dict(device_bytes={'cuda:1':64},host_bytes=128)
        alternate['api_kwargs']['chunk_size']=bounds['chunks'][-1]
        return dict(selected=selected,shortlist=[selected],admission_candidates=[selected,alternate],
            max_slowdown_fraction=0.,rejected=[],search_space=dict(input_file_identity=identity),
            candidates_evaluated=4,candidates_feasible=2,objective='analytical test',
            unresolved_timing_terms=[],scope='test')
    monkeypatch.setattr(module,'bounded_trait_plan',build)
    rng=np.random.default_rng(490)
    y=rng.normal(size=(64,9)).astype(np.float32);cov=rng.normal(size=(64,8)).astype(np.float32)
    qc=dict(phenotype_missing_cells=0,dropped_phenotype_columns=0,dropped_covariate_columns=0)
    def select(output=None):
        tuner=DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'output')
        return tuner.select(source(path),y,cov,qc,output={} if output is None else output)
    # Existing production inputs are stable before reuse. Keep this real-file
    # check outside the filesystem timestamp grace interval.
    time.sleep(1.05)
    return profile,config,path,calls,select


def test_second_select_reuses_analytical_plan_but_checks_live_memory(tmp_path,monkeypatch):
    import torch
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    first,audit,basis=select()
    assert audit['analytical_cache']['status']=='miss'
    second,hit,_=select()
    assert second==first and len(calls)==1 and hit['analytical_cache']['status']=='hit'
    assert hit['analytical_cache']['write_status']=='not_attempted'
    monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(0,1<<30))
    with pytest.raises(ValueError,match='current device capacity'):select()
    assert len(calls)==1
    assert not (tmp_path/'output').exists()


@pytest.mark.parametrize('resource',['device','host'])
def test_dense_cache_reselects_all_admitted_layouts_when_live_memory_shrinks(tmp_path,monkeypatch,resource):
    import torch
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    first,_,_=select();assert first['chunk_size']==128
    if resource=='device':monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(100,1<<30))
    else:monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:200)
    fallback,audit,_=select()
    assert fallback['chunk_size']==512 and len(calls)==1
    assert audit['analytical_cache']['status']=='hit'
    assert audit['planned_candidate_index']==0 and audit['selected']['candidate_index']==1
    assert audit['live_admission_rejections']==[dict(candidate_index=0,reason='current '+resource+' capacity')]
    monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(1<<30,1<<30))
    monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:1<<30)
    restored,audit,_=select()
    assert restored==first and len(calls)==1 and not audit['live_admission_rejections']


def test_dense_cache_without_full_admission_list_is_rebuilt(tmp_path,monkeypatch):
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    first,_,_=select()
    entry=next((tmp_path/'cache'/'torchgwas-analytical-plans-v1').glob('*.json'))
    record=json.loads(entry.read_text())
    record['plan'].pop('admission_candidates')
    record['sha256']=hashlib.sha256(canonical(record['plan'])).hexdigest()
    entry.write_text(json.dumps(record))
    restored,audit,_=select()
    assert restored==first and len(calls)==2
    assert audit['analytical_cache']['status']=='schema_upgrade'
    assert 'admission_candidates' in json.loads(entry.read_text())['plan']


@pytest.mark.parametrize('change',['source','price','bounds','budget','input','output'])
def test_material_request_change_invalidates_plan(tmp_path,monkeypatch,change):
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    select();stat=path.stat()
    if change=='source':profile['source_sha256']['test.py']='new'
    elif change=='price':profile['contexts'][0]['independent_price']=123
    elif change=='bounds':config['bounds']['chunks']=[32,128]
    elif change=='budget':config['joint']['device_memory_bytes']['cuda:1']-=1
    elif change=='input':
        path.write_bytes(b'replacement');os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        assert path.stat().st_size==stat.st_size and path.stat().st_mtime_ns==stat.st_mtime_ns
        assert path.stat().st_ctime_ns!=stat.st_ctime_ns
    _,audit,_=select(output=dict(store_beta=False) if change=='output' else None)
    assert len(calls)==2 and audit['analytical_cache']['status']==('recent_input' if change=='input' else 'miss')
    if change=='input':assert audit['analytical_cache']['write_status']=='recent_input'


def test_corrupt_cache_recomputes_and_unwritable_cache_does_not_change_answer(tmp_path,monkeypatch):
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    first,_,_=select()
    entry=next((tmp_path/'cache'/'torchgwas-analytical-plans-v1').glob('*.json'))
    record=json.loads(entry.read_text());record['plan']['selected']['api_kwargs']['chunk_size']=999
    entry.write_text(json.dumps(record))
    second,audit,_=select()
    assert first==second and len(calls)==2 and audit['analytical_cache']['status']=='invalid'
    obstacle=tmp_path/'not_a_directory';obstacle.write_text('preserve')
    config['plan_cache_dir']=str(obstacle)
    third,audit,_=select()
    assert first==third and len(calls)==3
    assert audit['analytical_cache']['write_status']=='unavailable'
    assert obstacle.read_text()=='preserve'


def test_hit_does_not_bypass_live_context_revalidation(tmp_path,monkeypatch):
    import torchgwas.detailed_autotune as module
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    select();checked=[]
    def validate(*args,**kwargs):
        checked.append(1)
        if len(checked)==2:raise ValueError('source changed after cache lookup')
        return dict(context_matches=True)
    monkeypatch.setattr(module,'validate_detailed_profile',validate)
    with pytest.raises(ValueError,match='after cache lookup'):select()
    assert len(calls)==1 and len(checked)==2


def test_cache_is_bounded_and_atomic_under_concurrent_access(tmp_path,monkeypatch):
    import torchgwas.analytical_plan_cache as module
    monkeypatch.setattr(module,'MAX_ENTRIES',3)
    plan=dict(selected=dict(chunk_size=128),search_space=dict(input_file_identity={}))
    keep=tmp_path/'torchgwas-analytical-plans-v1';keep.mkdir()
    (keep/'unrelated.json').write_text('preserve')
    def store_and_read(i):
        cache=AnalyticalPlanCache(tmp_path,dict(request=i))
        cache.store(plan)
        result=cache.load()
        assert result is None or result==plan
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(store_and_read,range(16)))
    final=AnalyticalPlanCache(tmp_path,dict(request=99));final.store(plan)
    entries=[p for p in keep.iterdir() if module._NAME.fullmatch(p.name)]
    assert len(entries)<=3 and final.load()==plan
    assert (keep/'unrelated.json').read_text()=='preserve'
    assert not list(keep.glob('.pending-*'))
    monkeypatch.setattr(module,'MAX_ENTRY_BYTES',32)
    assert final.load() is None and final.state=='oversized'
    final.store(plan);assert final.write_status=='oversized'


@pytest.mark.parametrize('directory',['',123,[]])
def test_bad_cache_path_fails_before_loading_or_creating_output(tmp_path,monkeypatch,directory):
    profile,config,path,_=fixture(tmp_path,monkeypatch)
    config['plan_cache_dir']=directory
    with pytest.raises(ValueError,match='plan_cache_dir'):
        DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_input_change_during_cache_hit_is_rejected_even_with_restored_mtime(tmp_path,monkeypatch):
    import torch
    profile,config,path,calls,select=cached_fixture(tmp_path,monkeypatch)
    select();stat=path.stat()
    def mutate(device):
        path.write_bytes(b'replacement');os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        return (1<<30,1<<30)
    monkeypatch.setattr(torch.cuda,'mem_get_info',mutate)
    with pytest.raises(ValueError,match='during analytical plan reuse'):select()
    assert len(calls)==1 and not (tmp_path/'output').exists()
