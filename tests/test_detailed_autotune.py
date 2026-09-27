import copy
import inspect
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.api import run_linear_gwas
from torchgwas.detailed_autotune import DetailedAutotune,validate_autotune_request
from torchgwas.pgen import PgenGenotype


def request():
    defaults={key:p.default for key,p in inspect.signature(run_linear_gwas).parameters.items()}
    defaults.update(genotype='input.pgen',phenotype='y.npy',output_dir='output',pgen_mode='hardcall')
    return defaults


@pytest.mark.parametrize('field,value',[
    ('genotype',np.zeros((32,4))),('genotype','input.bed'),('phenotype','y.csv'),
    ('output_dir',None),('genotype_format','bgen'),('pgen_mode','dosage'),
    ('compute_dtype','float64'),('device','cuda:1'),('sumstats_format','none'),
    ('sumstats_fsync',False),('pipeline_profile',{}),('chunk_size',128),
    ('trait_block',3),('trait_devices',['cuda:1']),('reader_workers',4),
    ('prefetch_chunks',4),('pgen_decode_workers',4),('reduce','top-k'),
    ('_internal_reduction',object()),('variant_range',(0,2)),
    ('p_value_threshold',0.01),('phenotype_table','y.tsv'),('covariates_table','c.tsv'),
    ('hardcall_store','store'),('zstd_read_workers',1)])
def test_conflicts_fail_before_any_load_or_writer(tmp_path,field,value):
    kwargs=request();kwargs.update(output_dir=tmp_path/'out',autotune_profile={},autotune_config={})
    kwargs[field]=value
    with patch('torchgwas.api.load_genotype',side_effect=AssertionError('must fail before load')):
        with pytest.raises(ValueError,match='Detailed autotune'):run_linear_gwas(**kwargs)
    assert not (tmp_path/'out').exists()


def test_profile_and_search_are_required_together():
    for key in ['autotune_profile','autotune_config']:
        with pytest.raises(ValueError,match='together'):
            run_linear_gwas('input.pgen','y.npy',**{key:{}})


def fixture(tmp_path,monkeypatch):
    import torchgwas.detailed_autotune as module
    path=tmp_path/'input.pgen';path.write_bytes(b'placeholder')
    gpu=dict(sm_count=12,max_threads_per_sm=2048,compute_capability=[9,0])
    context=dict(devices={'cuda:1':gpu},torch_version='test',torch_default_dtype='torch.float32',environment={'CUBLAS_WORKSPACE_CONFIG':None})
    bank={name:dict(reference_shape=[32,1,8],cpu_seconds=0.,non_cpu_seconds=0.) for name in ['residual_common','residual_block','design_common','design_block']}
    profile=dict(contexts=[dict(devices=['cuda:1'],profiles={'cuda:1':dict(setup_primitives=bank)})],source_sha256={'test.py':'hash'})
    memory=dict(**gpu,torch_version='test',cublas_workspace_config=None,cublas_handle_stream_pairs=2)
    config=dict(bounds=dict(chunks=[128,512],trait_blocks=[4,9]),qc_trait_block=4,
        joint=dict(cpu_workers=4,host_memory_bytes=1<<30,host_reserve_bytes=1<<20,
            device_memory_bytes={'cuda:1':1<<30},device_reserve_bytes=1<<20,
            device_memory_profiles={'cuda:1':memory},host_scenarios={}))
    monkeypatch.setattr(module,'execution_context',lambda *a,**k:context)
    monkeypatch.setattr(module,'validate_detailed_profile',lambda *a,**k:dict(context_matches=True))
    return profile,config,path,context


def source(path):
    value=object.__new__(PgenGenotype)
    value.__dict__.update(genotype_path=path,mode='hardcall',pgen_backend='native',
        _reader_sample_subset=None,_reorder_is_identity=True,native_dtype=np.dtype('int8'),
        native_encoding='dosage',_n_samples=64,_raw_n_variants=17)
    return value


@pytest.mark.parametrize('change,match',[
    (lambda c:c.update(selected={}), 'exactly'),
    (lambda c:c['bounds'].update(chunks=[128,128]),'duplicates'),
    (lambda c:c.update(qc_trait_block=5),'narrowest'),
    (lambda c:c['joint'].pop('host_reserve_bytes'),'Explicit'),
    (lambda c:c['joint']['device_memory_bytes'].update({'cuda:2':123}),'Exactly'),
    (lambda c:c['joint']['device_memory_profiles']['cuda:1'].update(sm_count=13),'properties'),
    (lambda c:c['joint']['device_memory_profiles']['cuda:1'].update(cublas_handle_stream_pairs=1),'two cuBLAS')])
def test_config_cannot_smuggle_a_plan_or_mismatch_resources(tmp_path,monkeypatch,change,match):
    profile,config,path,_=fixture(tmp_path,monkeypatch);change(config)
    with pytest.raises(ValueError,match=match):
        DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'out')
    assert not (tmp_path/'out').exists()


@pytest.mark.parametrize('change,match',[
    (lambda s,y,c:(setattr(s,'_reader_sample_subset',[0]) or y,c),'every sample'),
    (lambda s,y,c:(setattr(s,'_reorder_is_identity',False) or y,c),'every sample'),
    (lambda s,y,c:(setattr(s,'native_encoding','pgen_2bit') or y,c),'native int8'),
    (lambda s,y,c:(y.astype(np.float64),c),'float32'),
    (lambda s,y,c:(np.asfortranarray(y),c),'C-contiguous'),
    (lambda s,y,c:(y,c[:-1]),'aligned')])
def test_input_contract_is_checked_without_materializing_panel(tmp_path,monkeypatch,change,match):
    profile,config,path,_=fixture(tmp_path,monkeypatch)
    tuner=DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'out')
    src=source(path);y=np.ones((64,9),np.float32);c=np.ones((64,8),np.float32)
    y,c=change(src,y,c)
    with pytest.raises(ValueError,match=match):tuner.validate_inputs(src,y,c)


def test_rank_and_filtered_columns_fail_before_census(tmp_path,monkeypatch):
    profile,config,path,_=fixture(tmp_path,monkeypatch)
    tuner=DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'out')
    rng=np.random.default_rng(824);y=rng.normal(size=(64,9)).astype(np.float32)
    c=rng.normal(size=(64,8)).astype(np.float32);c[:,7]=c[:,0]
    qc=dict(phenotype_missing_cells=0,dropped_phenotype_columns=0,dropped_covariate_columns=0)
    with patch('torchgwas.detailed_autotune.bounded_trait_plan',side_effect=AssertionError('no census')):
        with pytest.raises(ValueError,match='setup primitives for covariate rank 7'):
            tuner.select(source(path),y,c,qc,output={})
        for field in ['phenotype_missing_cells','dropped_phenotype_columns']:
            with pytest.raises(ValueError,match='complete retained'):
                tuner.select(source(path),y,c,dict(qc,**{field:1}),output={})


@pytest.mark.parametrize('partition_axis',['trait','variant'])
def test_api_applies_selection_after_bounded_mmap_qc(tmp_path,monkeypatch,partition_axis):
    # Isolate the bridge from pricing; execute the real CPU scan and writer.
    import torchgwas.detailed_autotune as module
    import torchgwas.preprocess as prep
    from torchgwas.bed import PlinkBedGenotype
    from torchgwas.sumstats import open_binary_sumstats,open_binary_df
    from test_statistics import _write_bed
    rng=np.random.default_rng(881)
    bed=_write_bed(tmp_path/'input',rng.integers(0,3,size=(64,17)).astype(float))
    y=rng.normal(size=(64,9)).astype(np.float32);c=rng.normal(size=(64,8)).astype(np.float32)
    np.save(tmp_path/'y.npy',y)
    src=PlinkBedGenotype(bed)
    widths=[];mask=prep._phenotype_column_mask
    def checked(values):widths.append(values.shape[1]);return mask(values)
    monkeypatch.setattr(prep,'_phenotype_column_mask',checked)
    class Controller:
        devices=['cpu'];qc_trait_block=2
        def __init__(self,*a,**kw):pass
        def validate_inputs(self,source,phenotype,covariates):
            assert isinstance(phenotype,np.memmap) and phenotype.mode=='r'
        def select(self,source,phenotype,covariates,qc,*,output):
            assert max(widths)<=2 and qc['phenotype_columns_kept']==9
            assert output==dict(block_bytes=None,queue_depth=1,store_beta=True,fsync=True)
            selected=dict(chunk_size=5,device='cpu',reader_workers=2,prefetch_chunks=2)
            selected.update(dict(trait_block=4,trait_devices=['cpu']) if partition_axis=='trait' else dict(variant_devices=['cpu']))
            return selected,dict(selected_for_test=True),prep._covariate_basis(covariates)
    monkeypatch.setattr(module,'DetailedAutotune',Controller)
    with patch('torchgwas.api.load_genotype',return_value=(src,None,None,{})):
        result=run_linear_gwas(tmp_path/'fake.pgen',tmp_path/'y.npy',c,pgen_mode='hardcall',
            autotune_profile={},autotune_config={},sumstats_queue_depth=1,output_dir=tmp_path/'auto')
    meta=result.run_metadata
    assert (meta['chunk_size'],meta['trait_block'],meta['reader_workers'],meta['prefetch_chunks'])==(5,4 if partition_axis=='trait' else None,2,2)
    assert meta['autotune']==dict(selected_for_test=True)
    assert result.qc_summary['phenotype_qc_trait_block']==2
    run_linear_gwas(PlinkBedGenotype(bed),y,c,device='cpu',compute_dtype='float32',
        chunk_size=5,reader_workers=2,prefetch_chunks=2,
        sumstats_queue_depth=1,output_dir=tmp_path/'explicit',
        **(dict(trait_block=4,trait_devices=['cpu']) if partition_axis=='trait' else {}))
    actual=open_binary_sumstats(tmp_path/'auto'/'sumstats')
    expected=open_binary_sumstats(tmp_path/'explicit'/'sumstats')
    for i in (0,1,2):np.testing.assert_array_equal(np.asarray(actual[i]),np.asarray(expected[i]))
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path/'auto'/'sumstats')),
                                  np.asarray(open_binary_df(tmp_path/'explicit'/'sumstats')))


def test_cli_passes_profile_and_config_without_geometry_overrides(monkeypatch):
    import torchgwas.cli as cli
    received=[];monkeypatch.setattr(cli,'run_linear_gwas',lambda **k:received.append(k))
    args=cli._build_parser().parse_args(['linear','--genotype','x.pgen','--phenotype','y.npy',
        '--pgen-mode','hardcall','--autotune-profile','prices.json','--autotune-config','bounds.json',
        '--output-dir','out'])
    assert cli._run_linear(args)==0
    assert received[0]['autotune_profile']=='prices.json' and received[0]['autotune_config']=='bounds.json'
    assert all(received[0][k] is None for k in ['trait_block','chunk_size','reader_workers','prefetch_chunks'])

def test_unpriced_default_dtype_is_refused_before_planning(tmp_path,monkeypatch):
    profile,config,path,context=fixture(tmp_path,monkeypatch)
    context['torch_default_dtype']='torch.float64'
    with pytest.raises(ValueError,match='default tensor dtype'):
        DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'out')
    assert not (tmp_path/'out').exists()
