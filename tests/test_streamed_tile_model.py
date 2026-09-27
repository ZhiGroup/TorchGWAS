"""Real scan/writer graphs retain exact events when tiles are admitted lazily."""
import copy
from unittest.mock import patch

import pytest

from test_trait_tiling_model import candidate,input_path,component,runtime
from torchgwas import trait_tiling_model as model
from torchgwas.execution_graph import ExecutionGraph


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
@pytest.mark.parametrize('block',[None,48])
@pytest.mark.parametrize('spin',[0.,1.])
def test_full_multigpu_event_trace_and_all_report_fields_exact(input_path,policy,block,spin):
    c=candidate(input_path,count=2,traits=13,block_bytes=block)
    c['shared_storage_bytes_per_second']=2e7
    c['shared_links']=[dict(devices=c['devices'],h2d_bytes_per_second=2e7,d2h_bytes_per_second=1e7)]
    for tile in c['tiles']:
        tile['profile']['event_wait_cpu_fraction']=spin
        if tile['device']=='cuda:0':tile['profile']['fsync_seconds']*=100
    before=copy.deepcopy(c)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        expanded=model._trait_tiled_runtime(c,host_serial_fraction=.5,host_serial_policy=policy)
    gold=runtime(c,host_serial_policy=policy,return_graph=True).solve()
    original=ExecutionGraph.solve_chains;captured=[]
    def inspect(g,chains,**kwargs):
        result=original(g,chains,trace=True,**kwargs)
        captured.append(result);return result
    with patch.object(ExecutionGraph,'solve_chains',inspect):
        with patch.object(model,'_output_work',wraps=model._output_work) as writer:
            actual=runtime(c,host_serial_policy=policy)
            assert writer.call_count==5  # first/cached full for two devices, tail
    assert actual==expanded
    for key,value in gold.items():assert captured[0][key]==value,key
    assert captured[0]['peak_active_nodes']<captured[0]['scheduled_nodes']
    assert c==before


def test_changed_prices_in_later_request_cannot_reuse_previous_templates(input_path):
    c=candidate(input_path,count=2,traits=13)
    original=runtime(c)
    c['tiles'][4]['profile']['executor_cpu_seconds']*=50
    changed=runtime(c)
    assert changed['estimated_tile_seconds']!=original['estimated_tile_seconds']
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        expected=model._trait_tiled_runtime(c,host_serial_fraction=.5)
    assert changed==expected


def test_changing_per_device_capacity_keeps_original_full_graph(input_path):
    c=candidate(input_path,count=2,traits=13)
    c['tiles'][4]['profile']['h2d_bytes_per_second']*=.5
    with patch.object(ExecutionGraph,'solve_chains',side_effect=AssertionError('changed capacity')):
        actual=runtime(c)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        assert actual==model._trait_tiled_runtime(c,host_serial_fraction=.5)
