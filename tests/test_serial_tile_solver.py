"""Sequential cuts must preserve source work, sharing, waits and pinned reuse."""
import copy
from unittest.mock import patch

import pytest

from test_trait_tiling_model import candidate,input_path,component,runtime,plan
from torchgwas import trait_tiling_model as model
from torchgwas.mechanistic_torch import handoff_summary


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
@pytest.mark.parametrize('block',[None,48])
@pytest.mark.parametrize('spin',[0.,1.])
def test_decomposition_matches_full_contention_graph(input_path,policy,block,spin):
    c=candidate(input_path,count=1,traits=13,block_bytes=block)
    c['shared_storage_bytes_per_second']=3e7
    c['shared_links']=[dict(devices=c['devices'],h2d_bytes_per_second=2e7,d2h_bytes_per_second=1e7)]
    for tile in c['tiles']:tile['profile']['event_wait_cpu_fraction']=spin
    before=copy.deepcopy(c)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        expected=model._trait_tiled_runtime(c,host_serial_fraction=.5,host_serial_policy=policy)
    actual=runtime(c,host_serial_policy=policy)
    assert c==before
    assert actual['estimated_tile_seconds']==pytest.approx(expected['estimated_tile_seconds'],rel=2e-12,abs=1e-12)
    assert actual['handoffs']==pytest.approx(expected['handoffs'])
    assert {k:v for k,v in actual.items() if k not in ('handoffs','estimated_tile_seconds')}=={
        k:v for k,v in expected.items() if k not in ('handoffs','estimated_tile_seconds')}
    solved=runtime(c,host_serial_policy=policy,return_graph=True).solve()
    assert solved['seconds']==expected['estimated_tile_seconds']
    assert handoff_summary(solved)==expected['handoffs']


def test_reuses_only_equal_work_and_pinned_state_and_expires(input_path):
    c=candidate(input_path,count=1,traits=13)
    with patch.object(model,'_trait_tiled_runtime',wraps=model._trait_tiled_runtime) as expand:
        baseline=runtime(c)
        assert expand.call_count==3  # fresh full, cached full, tail
        assert model.trait_tiled_graph_chunks(c)==9
        c['tiles'][3]['profile']['executor_cpu_seconds']*=10
        changed=runtime(c)
        assert expand.call_count==7
        assert model.trait_tiled_graph_chunks(c)==12
        assert changed['estimated_tile_seconds']!=baseline['estimated_tile_seconds']
        # All source data, not just shape and prices, enter the key.
        c['tiles'][4]['data']['phenotype_c_contiguous']=False
        runtime(c)
        assert expand.call_count==12
        assert model.trait_tiled_graph_chunks(c)==15
    baseline['tiles'][1]['setup'].clear()
    assert runtime(candidate(input_path,count=1,traits=13))['tiles'][1]['setup']


def test_shared_census_serialization_is_scoped_and_value_based(input_path):
    c=candidate(input_path,count=1,traits=13)
    encoded=c['tiles'][0]['data']['encoded']
    for tile in c['tiles']:tile['data']['encoded']=encoded
    with patch.object(model.json,'dumps',wraps=model.json.dumps) as dumps:
        first=list(model._serial_tile_specs(c))
        assert sum(call.args[0] is encoded for call in dumps.call_args_list)==1
    encoded['new_census_field']='changed'
    second=list(model._serial_tile_specs(c))
    assert all(a[2]!=b[2] for a,b in zip(first,second))
    # Unsupported serialization only disables reuse; original evaluation works.
    encoded['new_census_field']=object()
    assert model.trait_tiled_graph_chunks(c)==21
    assert runtime(c)['genotype_passes']==7


def test_multiple_devices_keep_one_joint_graph(input_path):
    c=candidate(input_path,count=2,traits=13)
    with patch.object(model,'_trait_tiled_runtime',wraps=model._trait_tiled_runtime) as expand:
        runtime(c)
        assert expand.call_count==1 and expand.call_args.args[0] is c
    assert model.trait_tiled_graph_chunks(c)==21


def test_changing_transfer_capacity_keeps_composed_semantics(input_path):
    c=candidate(input_path,count=1,traits=13)
    c['tiles'][3]['profile']['h2d_bytes_per_second']*=.5
    assert model.trait_tiled_graph_chunks(c)==21
    with patch.object(model,'_trait_tiled_runtime',wraps=model._trait_tiled_runtime) as expand:
        runtime(c)
        assert expand.call_count==1 and expand.call_args.args[0] is c


def test_budget_counts_unique_graphs_but_keeps_logical_tile_limit(input_path):
    c=candidate(input_path,count=1,traits=13)
    # Two scenarios, three unique graphs each, three chunks per graph.
    assert plan([c],max_chunk_evaluations=18)['candidates_feasible']==1
    with patch('torchgwas.trait_tiling_plan.torch_trait_tiled_runtime',side_effect=AssertionError('graph started')):
        with pytest.raises(ValueError,match='max_chunk_evaluations'):
            plan([c],max_chunk_evaluations=17)
        with pytest.raises(ValueError,match='max_tiles'):
            plan([c],max_tiles=13)
        # A multi-device candidate still costs all of its overlapping tiles.
        with pytest.raises(ValueError,match='max_chunk_evaluations'):
            plan([c,candidate(input_path,count=2,traits=13)],max_chunk_evaluations=18+41)
