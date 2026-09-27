"""Effective-rank, source work and independent setup reference contracts."""
import copy
from unittest.mock import patch
import numpy as np
import pytest
from torchgwas.setup_work import setup_reference_shape,setup_primitive_bank,setup_work,setup_service,phase_work
from torchgwas.setup_calibration import primitive_requests,validate_measurement
from torchgwas.trait_candidate_space import prepare_trait_candidates,bounded_trait_plan
from torchgwas.trait_tiling_model import trait_tiled_memory,_prepare_graph
from test_trait_candidate_space import spec
from test_trait_tiling_model import component
from test_detailed_autotune import fixture,source
from torchgwas.detailed_autotune import DetailedAutotune

PHASES=['residual_common','residual_block','design_common','design_block']


def bank(rank):
    return {phase:dict(reference_shape=setup_reference_shape(rank),cpu_seconds=0. if not rank and phase=='residual_common' else .001,
                       non_cpu_seconds=0. if not rank and phase=='residual_common' else .002) for phase in PHASES}


@pytest.mark.parametrize('rank',[0,1,7,8,27,64])
def test_exact_references_and_no_rank_interpolation(rank):
    p=dict(setup_primitives=bank(8),setup_primitives_by_rank={str(rank):bank(rank)})
    assert setup_primitive_bank(p,rank)==bank(rank)
    assert setup_reference_shape(rank)==[max(32,rank+3),1,rank]
    with pytest.raises(ValueError,match='Missing independent'):setup_primitive_bank(p,rank+2)


def test_zero_rank_work_has_no_projection_or_covariate_cat():
    zero=setup_work(4096,512,0);positive=setup_work(4096,512,8)
    assert zero['h2d_bytes']==2*4*4096*512
    residual=next(p for p in zero['phases'] if p['phase']=='residual_block')
    assert residual['gemm_flops']==residual['gemm_bytes']==0
    assert residual['vector_ops']==8*4096*512
    assert phase_work('design_common',4096,1,0)['vector_bytes']==12*4096
    assert positive['phases'][1]['gemm_flops']>0
    for rank,columns in [(-1,None),(1,0),(1,4094),(True,None)]:
        with pytest.raises(ValueError):setup_work(4096,512,rank,covariate_columns=columns)


def test_rank_reference_requests_deduplicate_devices_and_preserve_existing(tmp_path):
    value=spec(tmp_path);contexts=value['contexts'];old=copy.deepcopy(contexts)
    requests=primitive_requests(contexts,[0,8,27])
    assert {(r['device'],r['rank']) for r in requests}=={(d,r) for d in ['cuda:0','cuda:1'] for r in [0,27]}
    assert contexts==old
    for ranks,limits in [([1,1],{}),([257],{}),([0,1],dict(max_ranks=1)),([True],{})]:
        with pytest.raises(ValueError):primitive_requests(contexts,ranks,**limits)


def test_complete_raw_repeats_required_and_zero_is_explicit():
    for rank in [0,27,64]:
        prices=bank(rank);rows=[]
        for phase,p in prices.items():
            for repeat in range(7):rows.append(dict(phase=phase,repeat=repeat,loops=100,
                reference_shape=setup_reference_shape(rank),measured=not(rank==0 and phase=='residual_common'),
                cpu_seconds=p['cpu_seconds'],non_cpu_seconds=p['non_cpu_seconds'],wall_seconds=p['cpu_seconds']+p['non_cpu_seconds']))
        record=dict(rank=rank,phases=prices,rows=rows)
        assert validate_measurement(record)==prices
        for mutate in [lambda r:r['rows'].pop(),lambda r:r['phases']['design_block'].update(cpu_seconds=1.),
                       lambda r:r['rows'][0].update(reference_shape=[32,2,rank])]:
            bad=copy.deepcopy(record);mutate(bad)
            with pytest.raises(ValueError):validate_measurement(bad)


def test_columns_bound_host_memory_but_rank_prices_gpu_and_df(tmp_path):
    value=spec(tmp_path);value['workload'].update(covariates=7,covariate_columns=9)
    for context in value['contexts']:
        for p in context['profiles'].values():
            p['setup_primitives_by_rank']={'7':bank(7)}
            for row in p['kernel_geometry']:row['C']=7
    space=prepare_trait_candidates(value['workload'],value['contexts'],output=value['output'],**value['bounds'])
    assert all(t['data']['covariates']==7 and t['data']['covariate_columns']==9 for c in space['candidates'] for t in c['tiles'])
    candidate=space['candidates'][0];larger=trait_tiled_memory(candidate)
    compact=copy.deepcopy(candidate)
    for t in compact['tiles']:t['data'].pop('covariate_columns')
    smaller=trait_tiled_memory(compact)
    assert larger['device_bytes']==smaller['device_bytes'] and larger['host_bytes']>smaller['host_bytes']
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        plan=bounded_trait_plan(**value)
    assert plan['selected']['candidate_index']>=0
    p=candidate['tiles'][0]['profile']
    a,_=_prepare_graph(setup_work(32,2,7,covariate_columns=9),p,.5,0,0)
    b,_=_prepare_graph(setup_work(32,2,7),p,.5,0,0)
    # Full-output workers reuse the API's shared basis. Retained columns still
    # affect the conservative host bound, but there is no per-tile SVD service.
    assert a.nodes['covariate_basis'][0]==b.nodes['covariate_basis'][0]==0.


@pytest.mark.parametrize('rank,columns',[ (0,0),(1,1),(8,9),(27,27) ])
def test_api_uses_effective_rank_after_qc_before_any_census(tmp_path,monkeypatch,rank,columns):
    profile,config,path,_=fixture(tmp_path,monkeypatch)
    profile['contexts'][0]['profiles']['cuda:1']['setup_primitives_by_rank']={str(rank):bank(rank)}
    rng=np.random.default_rng(34);y=rng.standard_normal((64,9),dtype=np.float32)
    c=None if not columns else rng.standard_normal((64,columns),dtype=np.float32)
    if columns>rank:c[:,-1]=c[:,0]
    tuner=DetailedAutotune(profile,config,input_path=path,output_path=tmp_path/'out')
    qc=dict(phenotype_missing_cells=0,dropped_phenotype_columns=0,dropped_covariate_columns=1)
    def capture(workload,*args,**kwargs):
        assert workload['covariates']==rank
        assert workload.get('covariate_columns',rank)==columns
        raise RuntimeError('captured rank')
    with patch('torchgwas.detailed_autotune.bounded_trait_plan',side_effect=capture):
        with pytest.raises(RuntimeError,match='captured rank'):tuner.select(source(path),y,c,qc,output={})


def test_tiny_rank_zero_and_positive_services_use_exact_reference():
    for rank in [0,1,27,64]:
        n,_,_=setup_reference_shape(rank)
        p=dict(setup_primitives_by_rank={str(rank):bank(rank)},cpu_fraction=1.,
            gpu_resources=dict(gpu_fraction=1.,hbm_bytes_per_second=1e9,fp32_flops_per_second=1e10),
            h2d_bytes_per_second=1e9,d2h_bytes_per_second=1e9,process_units=dict(numpy_copy_bytes=1e-10))
        result=setup_service(setup_work(n,1,rank,input_contiguous=True),p)
        assert result['seconds']==pytest.approx(.009 if rank==0 else .012)


def test_low_level_planner_rejects_different_covariate_column_counts(tmp_path):
    from torchgwas.trait_tiling_plan import detailed_trait_plan
    value=spec(tmp_path)
    a=prepare_trait_candidates(value['workload'],value['contexts'],output=value['output'],**value['bounds'])['candidates'][0]
    b=copy.deepcopy(a)
    for tile in b['tiles']:tile['data']['covariate_columns']=9
    with patch('torchgwas.trait_tiling_plan.torch_trait_tiled_runtime',side_effect=AssertionError('cost evaluation began')):
        with pytest.raises(ValueError,match='covariate column'):
            detailed_trait_plan([a,b],**value['joint'])
