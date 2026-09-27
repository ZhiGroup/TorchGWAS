"""Mode-correct public planning, immutable selector binding and cache reuse."""
import copy
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
from test_detailed_autotune import fixture, source, request
from test_significant_host_model import bank
from torchgwas.api import run_linear_gwas
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import sha256_file
from torchgwas.detailed_autotune import DetailedAutotune, validate_autotune_request
from torchgwas.analytical_plan_cache import input_identity


def bound_fixture(tmp_path,monkeypatch):
    profile,config,path,context=fixture(tmp_path,monkeypatch)
    dependencies=dict(source_sha256=profile['source_sha256'],execution_context=context)
    with patch('torchgwas.calibration_cache.time.time',return_value=100.):
        record=CalibrationParameterCache(tmp_path/'measurements').store('cpu_capacity',
            'significant_host_components',bank(),dependencies=dependencies,
            provenance={'scope':'synthetic test controls'},observed_unix_seconds=90.,max_age_seconds=20.)
    config.update(significant_host_prices=record['path'],plan_cache_dir=str(tmp_path/'plans'))
    config['joint']['occupancy_scenarios']={'none':'empty','all':'dense'}
    profile['component_artifacts']={record['path']:sha256_file(record['path'])}
    return profile,config,path,context


def controller(profile,config,path,**kwargs):
    return DetailedAutotune(profile,config,input_path=path,output_path=path.parent/'output',
        reduction='significant',**kwargs)


def test_entry_contract_allows_significant_but_rejects_dense_coalescing():
    options=request();options['reduce']='significant'
    validate_autotune_request(genotype=options['genotype'],phenotype=options['phenotype'],
                             output_dir=options['output_dir'],options=options)
    options['sumstats_block_bytes']=1024
    with pytest.raises(ValueError,match='indexed output'):
        validate_autotune_request(genotype=options['genotype'],phenotype=options['phenotype'],
                                 output_dir=options['output_dir'],options=options)


@pytest.mark.parametrize('fault',['missing','unbound','changed','dependencies','selector','expired','partition','threshold'])
def test_selector_binding_refuses_wrong_identity_age_or_mode(tmp_path,monkeypatch,fault):
    profile,config,path,context=bound_fixture(tmp_path,monkeypatch)
    kwargs={}
    if fault=='missing':config.pop('significant_host_prices')
    if fault=='unbound':profile['component_artifacts']={}
    if fault=='changed':Path(config['significant_host_prices']).write_text('{}')
    if fault=='dependencies':profile['source_sha256']={'test.py':'changed'}
    if fault=='selector':
        # A content-valid record with the wrong selector is still inapplicable.
        prices=bank();prices['host_selector']='obsolete'
        with patch('torchgwas.calibration_cache.time.time',return_value=100.):
            record=CalibrationParameterCache(tmp_path/'measurements').store('cpu_capacity','significant_host_components',
                prices,dependencies=dict(source_sha256=profile['source_sha256'],execution_context=context),
                provenance={'test':'wrong selector'},observed_unix_seconds=90.,max_age_seconds=20.)
        config['significant_host_prices']=record['path'];profile['component_artifacts']={record['path']:sha256_file(record['path'])}
    if fault=='partition':config['bounds']['partition_axes']=['variant']
    if fault=='threshold':kwargs['significance_threshold']=True
    with patch('torchgwas.calibration_cache.time.time',return_value=110. if fault=='expired' else 101.):
        with pytest.raises(ValueError):controller(profile,config,path,**kwargs)
    assert not (tmp_path/'output').exists()


def planning_fixture(tmp_path,monkeypatch):
    import torch
    import torchgwas.detailed_autotune as module
    profile,config,path,context=bound_fixture(tmp_path,monkeypatch)
    monkeypatch.setattr('torchgwas.analytical_plan_cache.input_is_stable',lambda identity:True)
    monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(1<<30,1<<30))
    monkeypatch.setattr(torch.cuda,'memory_reserved',lambda d:0)
    monkeypatch.setattr(torch.cuda,'memory_allocated',lambda d:0)
    monkeypatch.setattr(torch.backends.cuda.matmul,'allow_tf32',False)
    monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:1<<30)
    calls=[]
    def build(workload,contexts,*,bounds,joint,output,prices,significance_threshold):
        calls.append(dict(threshold=significance_threshold,output=output,prices=copy.deepcopy(prices)))
        chosen=dict(devices=['cuda:1'],required_environment={},memory=dict(device_bytes={'cuda:1':128},host_bytes=256),
            api_kwargs=dict(reduce='significant',significance_threshold=significance_threshold,chunk_size=128,
                trait_block=4,trait_devices=['cuda:1'],device='cuda:1',reader_workers=4,prefetch_chunks=4))
        return dict(selected=chosen,candidates=[chosen],rejected=[],search_space=dict(input_file_identity=input_identity(path)),
            candidates_evaluated=4,candidates_feasible=1,objective='test minimax',unpriced_terms=['test gap'],scope='test')
    monkeypatch.setattr(module,'bounded_significant_host_plan',build)
    monkeypatch.setattr(module,'bounded_trait_plan',lambda *a,**kw:pytest.fail('Reduced mode reached dense planner'))
    rng=np.random.default_rng(982);y=rng.normal(size=(64,9)).astype(np.float32)
    cov=rng.normal(size=(64,8)).astype(np.float32)
    qc=dict(phenotype_missing_cells=0,dropped_phenotype_columns=0,dropped_covariate_columns=0)
    output=dict(block_bytes=None,queue_depth=1,store_beta=True,fsync=True)
    def select(alpha=.01):
        return controller(profile,config,path,significance_threshold=alpha).select(source(path),y,cov,qc,output=output)
    return profile,config,path,calls,select,(y,cov,qc,output)


def test_threshold_and_mode_survive_dispatch_and_cached_reuse(tmp_path,monkeypatch):
    profile,config,path,calls,select,_=planning_fixture(tmp_path,monkeypatch)
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        settings,first,_=select();second,hit,_=select()
        assert settings==second and len(calls)==1
        assert settings['reduce']=='significant' and settings['significance_threshold']==.01
        assert first['analytical_cache']['status']=='miss' and hit['analytical_cache']['status']=='hit'
        assert first['reduction']=='significant' and first['significance_threshold']==.01
        assert first['unresolved_timing_terms']==['test gap'] and len(first['shortlist'])==1
        assert hit['reduction_calibration']['observed_unix_seconds']==90.
        _,changed,_=select(.02)
        assert len(calls)==2 and changed['analytical_cache']['key']!=hit['analytical_cache']['key']
        _,default,_=select(None)
        assert len(calls)==3 and default['significance_threshold'] is None
    with patch('torchgwas.calibration_cache.time.time',return_value=110.):
        with pytest.raises(ValueError,match='Expired'):select()
    assert len(calls)==3


def test_cache_hit_still_checks_live_memory(tmp_path,monkeypatch):
    import torch
    _,_,_,calls,select,_=planning_fixture(tmp_path,monkeypatch)
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        select();monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(0,1<<30))
        with pytest.raises(ValueError,match='current device capacity'):select()
    assert len(calls)==1


@pytest.mark.parametrize('after_qc',[True,False])
def test_expiry_during_qc_or_planning_prevents_output(tmp_path,monkeypatch,after_qc):
    profile,config,path,calls,_,values=planning_fixture(tmp_path,monkeypatch)
    y,cov,qc,output=values
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        tuner=controller(profile,config,path)
    times=[110.] if after_qc else [101.,110.]
    with patch('torchgwas.calibration_cache.time.time',side_effect=times):
        with pytest.raises(ValueError,match='Expired'):
            tuner.select(source(path),y,cov,qc,output=output)
    assert len(calls)==(0 if after_qc else 1) and not (tmp_path/'output').exists()


@pytest.mark.parametrize('fields',['t','beta+t'])
def test_api_routes_selected_tiles_to_indexed_significant_writer(tmp_path,monkeypatch,fields):
    # Mock only selection/context: execute actual CPU statistics and both writers.
    import torchgwas.detailed_autotune as module
    from torchgwas.bed import PlinkBedGenotype
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    from test_statistics import _write_bed
    rng=np.random.default_rng(715)
    bed=_write_bed(tmp_path/'input',rng.integers(0,3,size=(64,17)).astype(float))
    y=rng.normal(size=(64,9)).astype(np.float32);cov=rng.normal(size=(64,2)).astype(np.float32)
    np.save(tmp_path/'y.npy',y)
    class Controller:
        devices=['cpu'];qc_trait_block=2
        def __init__(self,*a,reduction,significance_threshold,**kwargs):
            assert reduction=='significant' and significance_threshold==1.
        def validate_inputs(self,*args):pass
        def select(self,*args,output):
            assert output['store_beta']==(fields=='beta+t')
            return dict(chunk_size=5,trait_block=4,trait_devices=['cpu'],device='cpu',
                reader_workers=2,prefetch_chunks=2,reduce='significant'),{'mode':'test'},None
    monkeypatch.setattr(module,'DetailedAutotune',Controller)
    kwargs=dict(reduce='significant',significance_threshold=1.,sumstats_fields=fields,compute_dtype='float32',sumstats_queue_depth=1)
    with patch('torchgwas.api.load_genotype',return_value=(PlinkBedGenotype(bed),None,None,{})):
        result=run_linear_gwas(tmp_path/'fake.pgen',tmp_path/'y.npy',cov,pgen_mode='hardcall',
            autotune_profile={},autotune_config={},output_dir=tmp_path/'auto',**kwargs)
    run_linear_gwas(PlinkBedGenotype(bed),y,cov,device='cpu',chunk_size=5,trait_block=4,trait_devices=['cpu'],
        reader_workers=2,prefetch_chunks=2,output_dir=tmp_path/'explicit',**kwargs)
    def rows(directory):
        manifest,parts=open_indexed_sumstats(directory/'sumstats')
        chunks=list(parts)
        assert manifest['format']=='torchgwas-indexed-sumstats' and manifest['rows']==17*9
        values={key:np.concatenate([chunk[key] for chunk in chunks]) for key in chunks[0]}
        order=np.argsort(values['variant_index']*9+values['trait_index'])
        return {key:value[order] for key,value in values.items()}
    actual,expected=rows(tmp_path/'auto'),rows(tmp_path/'explicit')
    assert set(actual)==set(expected) and ('beta' in actual)==(fields=='beta+t')
    for key in actual:np.testing.assert_array_equal(actual[key],expected[key])
    assert result.run_metadata['trait_block']==4 and result.run_metadata['autotune']=={'mode':'test'}


@pytest.mark.parametrize('resource',['device','host'])
def test_fresh_admission_can_choose_alternate_without_mutating_cached_ranking(tmp_path,monkeypatch,resource):
    import torch
    import torchgwas.detailed_autotune as module
    _,_,_,calls,select,_=planning_fixture(tmp_path,monkeypatch)
    original=module.bounded_significant_host_plan
    def ranked(*args,**kwargs):
        plan=original(*args,**kwargs)
        plan['selected']['candidate_index']=0
        alternate=copy.deepcopy(plan['selected']);alternate['candidate_index']=1
        alternate['memory']=dict(device_bytes={'cuda:1':64},host_bytes=128)
        alternate['api_kwargs']['chunk_size']=32
        plan['candidates'].append(alternate)
        return plan
    monkeypatch.setattr(module,'bounded_significant_host_plan',ranked)
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        first,_,_=select();assert first['chunk_size']==128
        if resource=='device':monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(100,1<<30))
        else:monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:200)
        fallback,audit,_=select()
        assert fallback['chunk_size']==32 and len(calls)==1
        assert audit['analytical_cache']['status']=='hit'
        assert audit['planned_candidate_index']==0 and audit['selected']['candidate_index']==1
        assert audit['live_admission_rejections']==[dict(candidate_index=0,reason='current '+resource+' capacity')]
        monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(1<<30,1<<30))
        monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:1<<30)
        restored,audit,_=select()
        assert restored==first and len(calls)==1 and not audit['live_admission_rejections']
