import copy
import inspect
from unittest.mock import patch

import numpy as np
import pytest

from test_mechanistic_shapes import fixture, component
from test_pgen_native_reader import write_pgen
from torchgwas.api import run_linear_gwas
from torchgwas.autotune import detailed_trait_plan
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.native_control_work import native_control_work
from torchgwas.pgen_work_census import census
from torchgwas.trait_tiling_model import trait_tiled_shape, trait_tiled_memory, torch_trait_tiled_runtime


def candidate(path, width=2, count=2, traits=5, block_bytes=48):
    devices = [f'cuda:{i}' for i in range(count)]
    tiles = []
    encoded = census(path,4,include_chunks=True)
    for index,start in enumerate(range(0,traits,width)):
        k = min(width,traits-start)
        data,profile = fixture(k=k,c=8)
        data['encoded'] = copy.deepcopy(encoded)
        profile.update(chunk_markers=4,decode_workers=2//count,validate_range=True,event_wait_cpu_fraction=1.,
            result_ownership='borrowed',write_bytes_per_second=1e8,fsync_seconds=1e-6,
            pin_cpu_seconds_per_page=1e-6,pin_driver_seconds_per_page=1e-6,pin_cached_cpu_seconds_per_call=1e-7,
            writeback_service=dict(pagecache_seconds_per_byte=1e-9,storage_seconds_per_byte=1e-8,
                submit_seconds=1e-6,wait_seconds=1e-6,fadvise_seconds=1e-6),
            control_primitives={name:1e-6 for counts in native_control_work().values() for name in counts},
            result_finish_service=dict(cpu_seconds=1e-6,serial_cpu_seconds=0.,baseline_copy_bytes=0,
                replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False),
            setup_primitives={name:dict(reference_shape=[32,1,8],cpu_seconds=1e-6,non_cpu_seconds=1e-6)
                for name in ['residual_common','residual_block','design_common','design_block']})
        profile['gpu_resources']['gpu_fraction'] = 1.
        profile['process_units'].update(bytearray_zero_bytes=1e-10,covariate_basis_work=1e-10)
        profile['kernel_geometry'] = [dict(N=32,B=b,K=k,C=8,validate_range=True,kernels=[]) for b in (4,2)]
        tiles.append(dict(trait_range=[start,start+k],device=devices[index%count],data=data,profile=profile))
    return dict(tiles=tiles,trait_block=width,devices=devices,
        shared_capacities=dict(cpu=2.,dram=1e9,input=1e8,output=1e8),
        output=dict(block_bytes=block_bytes,queue_depth=1,store_beta=True,fsync=True))


@pytest.fixture
def input_path(tmp_path):
    path = tmp_path/'input.pgen'
    write_pgen(path,np.arange(10*32,dtype=np.uint8).reshape(10,32)%4)
    return path


def runtime(c,**kwargs):
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return torch_trait_tiled_runtime(c,host_serial_fraction=.5,**kwargs)


@pytest.mark.parametrize('block_bytes',[None,48])
@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
def test_tile_graph_preserves_lifetimes_and_complete_output(input_path,block_bytes,policy):
    c = candidate(input_path,block_bytes=block_bytes)
    original = copy.deepcopy(c)
    graph = runtime(c,host_serial_policy=policy,return_graph=True)
    result = graph.solve()
    assert c==original
    assert result['start']['tile2:prepare:covariate_basis']>=result['end']['tile0:complete']
    assert result['start']['tile2:submit_decode:0']>=result['end']['tile0:writer:df:fsync']
    for tile in range(3):
        for chunk in range(3):
            # No writer copy may use a ring slot before its result is yielded.
            for name in graph.nodes:
                if name.startswith(f'tile{tile}:writer:beta:copy:{chunk}:'):
                    assert result['start'][name]>=result['end'][f'tile{tile}:resolve:{chunk}']
        if block_bytes is not None:
            assert result['start'][f'tile{tile}:prepare:covariate_basis']>=result['end'][f'tile{tile}:writer:open']
        else:
            assert result['start'][f'tile{tile}:writer:open']>=result['end'][f'tile{tile}:resolve:0']
    # Resource work conservation supplies a lower bound independent of DAG topology.
    for resource,capacity in graph.capacities.items():
        work = sum(graph.nodes[name][0]*demands.get(resource,0.) for name,demands in graph.demands.items())
        assert result['seconds']+1e-10>=work/capacity
    report = runtime(c,host_serial_policy=policy)
    assert report['genotype_passes']==3
    assert report['binary_payload_bytes']==12*10*5+4*10*3
    assert report['df_payload_bytes']==4*10*3
    assert sum(tile['writer_copy_bytes'] for tile in report['tiles'])==report['binary_payload_bytes']
    assert report['prediction_complete'] is False


def test_cached_pins_follow_size_classes_and_writer_close(input_path):
    c = candidate(input_path,traits=7,count=1)
    report = runtime(c)
    assert report['tiles'][0]['pin_fresh_pages']>0
    assert report['tiles'][1]['pin_fresh_pages']==0
    assert report['tiles'][1]['pin_cached_calls']==10
    memory = trait_tiled_memory(c)
    assert memory['host_bytes']>memory['pinned_cache_bytes_by_device']['cuda:0']
    assert memory['device_bytes']['cuda:0']>0
    assert memory['unresolved_memory_terms']


def test_shared_storage_and_pcie_capacities_are_not_multiplied(input_path):
    c = candidate(input_path)
    c['shared_storage_bytes_per_second'] = 1.
    c['shared_links'] = [dict(devices=c['devices'],h2d_bytes_per_second=1.,d2h_bytes_per_second=1.)]
    graph = runtime(c,return_graph=True)
    result = graph.solve()
    for resource in ['storage','link:0:h2d','link:0:d2h']:
        work = sum(graph.nodes[name][0]*demands.get(resource,0.) for name,demands in graph.demands.items())
        assert work>0
        assert result['seconds']+1e-8>=work
    assert runtime(c)['estimated_tile_seconds']>runtime(candidate(input_path))['estimated_tile_seconds']


@pytest.mark.parametrize('mutate,match',[
    (lambda c:c.update(observed_seconds=1),'Unknown tiled'),
    (lambda c:c['tiles'][1].update(device='cuda:0'),'round-robin'),
    (lambda c:c['tiles'][1].update(trait_range=[3,5]),'consecutive'),
    (lambda c:c['tiles'][0]['profile'].update(result_ownership='owned'),'borrowed'),
    (lambda c:c['tiles'][0]['data']['encoded'].update(variant_range=[0,5]),'full genotype'),
    (lambda c:c['output'].update(fsync=False),'durable'),
    (lambda c:c['shared_capacities'].update(input=c['shared_capacities']['input']*1.01),'Shared capacity'),
])
def test_refuses_different_executor_contracts(input_path,mutate,match):
    c = candidate(input_path);mutate(c)
    with pytest.raises(ValueError,match=match):trait_tiled_shape(c)


def options():
    return dict(host_scenarios={'fluid':dict(host_serial_fraction=0.,host_serial_policy='fluid'),
                               'held':dict(host_serial_fraction=1.,host_serial_policy='held-first')},
        cpu_workers=2,host_memory_bytes=2**30,device_memory_bytes={'cuda:0':2**30,'cuda:1':2**30})


def plan(candidates,**kwargs):
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return detailed_trait_plan(candidates,**(options()|kwargs))


def test_planner_produces_valid_api_arguments_and_enforces_bounds(input_path):
    candidates = [candidate(input_path,count=1),candidate(input_path)]
    result = plan(candidates)
    assert result['candidates_feasible']==2
    assert len(result['admission_candidates'])==2
    assert all('scenarios' not in row for row in result['admission_candidates'])
    assert result['selection_validated'] is False
    assert result['predicted_selection_penalty_fraction']==0
    for row in result['shortlist']:
        kwargs = row['api_kwargs']|dict(phenotype='Y.npy',covariates='C.npy',output_dir='out')
        inspect.signature(run_linear_gwas).bind(**kwargs)
        assert kwargs['sumstats_fields']=='beta+t'
        assert set(row['scenarios'])=={'fluid','held'}
    assert plan(candidates,device_memory_bytes={'cuda:0':2**30,'cuda:1':1})['candidates_feasible']==1
    with pytest.raises(ValueError,match='max_candidates'):plan(candidates,max_candidates=1)
    with pytest.raises(ValueError,match='max_scenario_evaluations'):plan(candidates,max_scenario_evaluations=1)
    with pytest.raises(ValueError,match='max_tiles'):plan(candidates,max_tiles=1)
    with pytest.raises(ValueError,match='max_chunk_evaluations'):plan(candidates,max_chunk_evaluations=1)
    with pytest.raises(ValueError,match='shared_reader_budget'):plan(candidates,cpu_workers=1)


def test_missing_component_prices_never_use_coarse_fallback(input_path):
    candidates = [candidate(input_path,count=1),candidate(input_path)]
    candidates[0]['tiles'][0]['profile'].pop('setup_primitives')
    result = plan(candidates)
    assert result['candidates_feasible']==1
    assert result['rejected'][0]['reason']=='unsupported_model_context'
    assert 'setup primitives' in result['rejected'][0]['detail']


def test_graph_composition_preserves_private_tokens_and_waits():
    local = ExecutionGraph();local.token_capacities['slot']=1
    local.add('ready',2.);local.add('attempt',1.)
    local.add_wait_delay('wake',3.,'attempt','ready')
    local.add('use',1.,['wake'])
    local.token_actions['use']={'acquire':{'slot':1},'release_finish':{'slot':1}}
    graph = ExecutionGraph();end=graph.compose(local,'first:')
    graph.compose(local,'second:',[end])
    result=graph.solve()
    assert result['seconds']==12.
    assert result['conditional_delays']['first:wake']['blocked']
    assert set(graph.token_capacities)=={'first:slot','second:slot'}
    with pytest.raises(ValueError,match='prefix'):graph.compose(local,'first:')


def test_bulk_numeric_setup_copy_uses_cpu_and_dram_without_held_gil(input_path):
    from torchgwas.trait_tiling_model import _prepare_graph
    from torchgwas.setup_work import setup_work
    profile=candidate(input_path)['tiles'][0]['profile']
    for n,k in [(32,1),(32,15),(32,16),(4096,4096)]:
        work=setup_work(n,k,reuse_observed_counts=True)
        graph,report=_prepare_graph(work,profile,1.,0,0)
        for i,(phase,cost) in enumerate(zip(work['phases'],report['phases'])):
            seconds=cost['seconds'];r=graph.demands['phase:'+str(i)];q=profile['cpu_fraction']
            held=cost['fixed_cpu_seconds']+(cost['host_copy_seconds'] if n*phase['traits']<=500 else 0)
            assert r['host_serial']*seconds==pytest.approx(q*held)
            assert r['cpu']*seconds==pytest.approx(q*(cost['fixed_cpu_seconds']+cost['host_copy_seconds']))
            assert r['dram']*seconds==pytest.approx(2*phase['host_copy_bytes'])

def test_scoped_scan_work_cache_keys_all_data_and_prices_and_expires():
    from torchgwas.trait_tiling_model import reuse_trait_work,_cached_scan_work,_SCAN_WORK_CACHE
    @reuse_trait_work
    def plan():
        data={'samples':32};profile={'capacity':1.}
        one=_cached_scan_work(data,profile)
        assert _cached_scan_work(copy.deepcopy(data),copy.deepcopy(profile)) is one
        profile['capacity']=2.
        assert _cached_scan_work(data,profile) is not one
        data['samples']=64
        _cached_scan_work(data,profile)
    with patch('torchgwas.trait_tiling_model.torch_scan_work',side_effect=lambda *a:{}) as scan:
        plan();assert scan.call_count==3
        plan();assert scan.call_count==6
    assert _SCAN_WORK_CACHE.get() is None


def test_scoped_scan_cache_preserves_every_graph_field_and_returned_work(input_path):
    from torchgwas.trait_tiling_model import reuse_trait_work
    c=candidate(input_path,traits=7,count=1)
    expected=runtime(c,host_serial_policy='held-last',return_graph=True)
    @reuse_trait_work
    def plan():
        first=runtime(c,host_serial_policy='held-last',return_graph=True)
        first.nodes.clear()
        return runtime(c,host_serial_policy='held-last',return_graph=True)
    actual=plan()
    assert actual.__dict__==expected.__dict__
    assert actual.solve()==expected.solve()
