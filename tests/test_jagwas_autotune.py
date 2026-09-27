"""Full-panel public JAGWAS dispatch and immutable service reuse."""
import copy
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
from test_detailed_autotune import fixture, source, request
from test_jagwas_actual_candidate import writer_prices
from torchgwas.api import run_linear_gwas
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import sha256_file
from torchgwas.detailed_autotune import DetailedAutotune, validate_autotune_request
from torchgwas.analytical_plan_cache import input_identity


def bound_fixture(tmp_path,monkeypatch):
    profile,config,path,context=fixture(tmp_path,monkeypatch)
    values=dict(writer_prices=writer_prices(),preparation_services={'test':{'fluid':dict(
        library_arithmetic={'cuda:1':'scalar'},shared_cpu_steps=[dict(seconds=.001,resources={'cpu':1.})],
        finalize=[dict(seconds=.001,resources={'cpu':1.})])}})
    with patch('torchgwas.calibration_cache.time.time',return_value=100.):
        record=CalibrationParameterCache(tmp_path/'measurements').store('cpu_capacity','jagwas_components',values,
            dependencies=dict(source_sha256=profile['source_sha256'],execution_context=context),
            provenance={'scope':'synthetic test controls'},observed_unix_seconds=90.,max_age_seconds=20.)
    config.update(jagwas_services=record['path'],plan_cache_dir=str(tmp_path/'plans'))
    config['bounds'].pop('trait_blocks')
    config['joint']['occupancy_scenarios']={'all':'dense'}
    profile['component_artifacts']={record['path']:sha256_file(record['path'])}
    return profile,config,path,context


def controller(profile,config,path,**kwargs):
    return DetailedAutotune(profile,config,input_path=path,output_path=path.parent/'output',
        reduction='jagwas',**kwargs)


def test_entry_contract_allows_jagwas_without_dense_coalescing():
    options=request();options['reduce']='jagwas'
    validate_autotune_request(genotype=options['genotype'],phenotype=options['phenotype'],
                             output_dir=options['output_dir'],options=options)
    options['sumstats_block_bytes']=1024
    with pytest.raises(ValueError,match='indexed output'):
        validate_autotune_request(genotype=options['genotype'],phenotype=options['phenotype'],
                                 output_dir=options['output_dir'],options=options)


@pytest.mark.parametrize('fault',['missing','unbound','changed','dependencies','expired','traits','partition','mode','threshold'])
def test_binding_and_full_panel_contract(tmp_path,monkeypatch,fault):
    profile,config,path,_=bound_fixture(tmp_path,monkeypatch);kwargs={}
    if fault=='missing':config.pop('jagwas_services')
    if fault=='unbound':profile['component_artifacts']={}
    if fault=='changed':Path(config['jagwas_services']).write_text('{}')
    if fault=='dependencies':profile['source_sha256']={'changed.py':'hash'}
    if fault=='traits':config['bounds']['trait_blocks']=[9]
    if fault=='partition':config['bounds']['partition_axes']=['variant']
    if fault=='mode':config['significant_host_prices']=config['jagwas_services']
    if fault=='threshold':kwargs['significance_threshold']=1.
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
    def build(workload,contexts,*,bounds,joint,output,prices,preparation_services):
        assert set(bounds)=={'chunks'} and output['store_beta'] is False
        calls.append(copy.deepcopy(dict(prices=prices,preparation_services=preparation_services)))
        rows=[]
        for index,chunk in enumerate((128,32)):
            rows.append(dict(candidate_index=index,devices=['cuda:1'],chunk_size=chunk,required_environment={},
                memory=dict(device_bytes={'cuda:1':chunk},host_bytes=2*chunk),
                api_kwargs=dict(reduce='jagwas',chunk_size=chunk,variant_devices=['cuda:1'],
                    device='cuda:1',reader_workers=4,prefetch_chunks=4)))
        return dict(selected=rows[0],candidates=rows,rejected=[],search_space=dict(input_file_identity=input_identity(path)),
            candidates_evaluated=2,candidates_feasible=2,objective='test minimax',unpriced_terms=['test gap'],scope='test')
    monkeypatch.setattr(module,'bounded_jagwas_plan',build)
    for name in ('bounded_significant_host_plan','bounded_trait_plan'):
        monkeypatch.setattr(module,name,lambda *a,**kw:pytest.fail('JAGWAS reached another mode planner'))
    rng=np.random.default_rng(912);y=rng.normal(size=(64,9)).astype(np.float32)
    cov=rng.normal(size=(64,8)).astype(np.float32)
    qc=dict(phenotype_missing_cells=0,dropped_phenotype_columns=0,dropped_covariate_columns=0)
    output=dict(block_bytes=None,queue_depth=1,store_beta=False,fsync=True)
    def select():return controller(profile,config,path).select(source(path),y,cov,qc,output=output)
    return profile,config,path,calls,select,(y,cov,qc,output)


def test_cache_reuse_preserves_original_record_and_checks_live_memory(tmp_path,monkeypatch):
    import torch
    _,_,_,calls,select,_=planning_fixture(tmp_path,monkeypatch)
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        first,audit,_=select()
        assert audit['analytical_cache']['status']=='miss' and audit['reduction']=='jagwas'
        assert 'trait_block' not in first and first['variant_devices']==['cuda:1']
        monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(64,1<<30))
        fallback,hit,_=select()
        assert fallback['chunk_size']==32 and hit['analytical_cache']['status']=='hit'
        assert hit['planned_candidate_index']==0 and hit['selected']['candidate_index']==1
        assert hit['reduction_calibration']['observed_unix_seconds']==90.
        assert hit['reduction_calibration']['record_sha256']==audit['reduction_calibration']['record_sha256']
        monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(1<<30,1<<30))
        restored,_,_=select();assert restored==first and len(calls)==1
    with patch('torchgwas.calibration_cache.time.time',return_value=110.):
        with pytest.raises(ValueError,match='Expired'):select()
    assert len(calls)==1


@pytest.mark.parametrize('after_qc',[True,False])
def test_expiry_during_qc_or_planning_prevents_output(tmp_path,monkeypatch,after_qc):
    profile,config,path,calls,_,values=planning_fixture(tmp_path,monkeypatch)
    y,cov,qc,output=values
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        tuner=controller(profile,config,path)
    with patch('torchgwas.calibration_cache.time.time',side_effect=[110.] if after_qc else [101.,110.]):
        with pytest.raises(ValueError,match='Expired'):tuner.select(source(path),y,cov,qc,output=output)
    assert len(calls)==(0 if after_qc else 1) and not (tmp_path/'output').exists()


def test_residual_rank_rejected_before_planning(tmp_path,monkeypatch):
    profile,config,path,calls,_,values=planning_fixture(tmp_path,monkeypatch)
    _,cov,qc,output=values
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        tuner=controller(profile,config,path)
        with pytest.raises(ValueError,match='residual phenotype rank'):
            tuner.select(source(path),np.ones((64,56),np.float32),cov,qc,output=output)
    assert not calls


@pytest.mark.parametrize('fields',['t','beta+t'])
def test_api_separates_qc_blocks_and_routes_full_panel_joint_output(tmp_path,monkeypatch,fields):
    import torchgwas.detailed_autotune as module
    import torchgwas.preprocess as prep
    import torchgwas.reduction_tensor_work as factor
    from torchgwas.bed import PlinkBedGenotype
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    from test_statistics import _write_bed
    rng=np.random.default_rng(901)
    bed=_write_bed(tmp_path/'input',rng.integers(0,3,size=(64,17)).astype(float))
    y=rng.normal(size=(64,9)).astype(np.float32);cov=rng.normal(size=(64,2)).astype(np.float32)
    np.save(tmp_path/'y.npy',y)
    widths=[];mask=prep._phenotype_column_mask;events=[]
    def checked(values):widths.append(values.shape[1]);return mask(values)
    monkeypatch.setattr(prep,'_phenotype_column_mask',checked)
    capacity=factor.require_jagwas_factor_capacity
    def capacity_check(n,k,devices,**kw):
        events.append(('capacity',devices));assert events[0]=='selection'
        return capacity(n,k,devices,**kw)
    monkeypatch.setattr(factor,'require_jagwas_factor_capacity',capacity_check)
    class Controller:
        devices=['cpu'];qc_trait_block=2
        def __init__(self,*a,reduction,**kw):assert reduction=='jagwas'
        def validate_inputs(self,*args):pass
        def select(self,genotype,phenotype,covariates,qc,*,output):
            events.append('selection');assert max(widths)<=2 and phenotype.shape[1]==9
            assert output['store_beta'] is False
            return dict(chunk_size=5,variant_devices=['cpu'],device='cpu',
                reader_workers=2,prefetch_chunks=2,reduce='jagwas'),{'mode':'test'},None
    monkeypatch.setattr(module,'DetailedAutotune',Controller)
    kwargs=dict(reduce='jagwas',sumstats_fields=fields,compute_dtype='float32',sumstats_queue_depth=1)
    with patch('torchgwas.api.load_genotype',return_value=(PlinkBedGenotype(bed),None,None,{})):
        result=run_linear_gwas(tmp_path/'fake.pgen',tmp_path/'y.npy',cov,pgen_mode='hardcall',
            autotune_profile={},autotune_config={},output_dir=tmp_path/'auto',**kwargs)
    assert events==['selection',('capacity',['cpu'])]
    monkeypatch.setattr(factor,'require_jagwas_factor_capacity',capacity)
    run_linear_gwas(PlinkBedGenotype(bed),y,cov,device='cpu',chunk_size=5,variant_devices=['cpu'],
        reader_workers=2,prefetch_chunks=2,output_dir=tmp_path/'explicit',**kwargs)
    def rows(directory):
        manifest,parts=open_indexed_sumstats(directory/'sumstats');chunks=list(parts)
        assert manifest['shape']==[17,9] and manifest['df']==9
        values={key:np.concatenate([chunk[key] for chunk in chunks]) for key in chunks[0]}
        assert set(values)=={'variant_index','chi2'} and len(values['variant_index'])==17
        order=np.argsort(values['variant_index'])
        return {key:value[order] for key,value in values.items()}
    actual,expected=rows(tmp_path/'auto'),rows(tmp_path/'explicit')
    for key in actual:np.testing.assert_array_equal(actual[key],expected[key])
    assert result.run_metadata['trait_block'] is None and result.run_metadata['variant_devices']==['cpu']
    assert result.qc_summary['phenotype_qc_trait_block']==2


def test_public_selection_updates_native_reader_preference(tmp_path,monkeypatch):
    import torch
    import torchgwas.detailed_autotune as module
    import torchgwas.native_scan as native
    from torchgwas.pgen import PgenGenotype
    from test_jagwas_variant_devices import fixture as data_fixture, rows, reference
    if not torch.cuda.is_available():pytest.skip('CUDA device required')
    for key,value in dict(TORCHGWAS_PGEN_BACKEND='native',TORCHGWAS_PGEN_PACKED='0',TORCHGWAS_NATIVE_STATS='0').items():
        monkeypatch.setenv(key,value)
    path,calls,y,cov=data_fixture(tmp_path,'pgen',missing=False)
    source=PgenGenotype(path,mode='hardcall',reader_workers=1)
    observed=[];loader=native.PinnedDosageLoader
    def load(*args,**kwargs):
        value=loader(*args,**kwargs)
        observed.append((value.decode_workers_requested,value.decode_workers_effective))
        return value
    monkeypatch.setattr(native,'PinnedDosageLoader',load)
    class Controller:
        devices=['cuda:0'];qc_trait_block=2
        def __init__(self,*args,**kwargs):pass
        def validate_inputs(self,*args):pass
        def select(self,*args,**kwargs):
            return dict(chunk_size=5,variant_devices=['cuda:0'],device='cuda:0',
                reader_workers=2,prefetch_chunks=4,reduce='jagwas'),{},None
    monkeypatch.setattr(module,'DetailedAutotune',Controller)
    with patch('torchgwas.api.load_genotype',return_value=(source,None,None,{})):
        result=run_linear_gwas(path,y,cov,pgen_mode='hardcall',reduce='jagwas',
            autotune_profile={},autotune_config={},output_dir=tmp_path/'auto')
    assert observed==[(2,2)] and result.run_metadata['reader_workers']==2
    _,actual=rows(tmp_path/'auto');expected=reference(calls,y,cov,(0,23))
    assert actual.keys()==expected.keys()
    np.testing.assert_allclose(list(actual.values()),[expected[i] for i in actual],rtol=3e-4,atol=3e-4)
