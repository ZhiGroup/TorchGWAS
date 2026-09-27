"""Finite significant-pairs search, scenario identity, and within-plan reuse."""
import copy
import json
from unittest.mock import patch

import pytest

from test_trait_candidate_space import spec as dense_spec
from test_trait_tiling_model import candidate, input_path
from test_mechanistic_shapes import component
from test_significant_host_model import bank, options
from torchgwas.mechanistic_torch import torch_scan_work
from torchgwas.significant_candidate_space import prepare_significant_host_candidates, bounded_significant_host_plan
from torchgwas.significant_host_plan import detailed_significant_host_plan
from torchgwas.significant_host_work import significant_host_memory


def specification(tmp_path):
    value=dense_spec(tmp_path)
    value['joint']=options();threshold=value['joint'].pop('significance_threshold')
    value.update(significance_threshold=threshold,prices=bank())
    value['output']['block_bytes']=None
    return value


def proposed(value):
    return prepare_significant_host_candidates(value['workload'],value['contexts'],output=value['output'],**value['bounds'])


def test_bounded_space_preserves_tiles_and_matches_explicit_plan(tmp_path):
    value=specification(tmp_path);before=copy.deepcopy(value);space=proposed(value)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        expected=detailed_significant_host_plan(space['candidates'],value['prices'],
            significance_threshold=value['significance_threshold'],**value['joint'])
        actual=bounded_significant_host_plan(**value)
    assert actual['selected']==expected['selected'] and actual['candidates']==expected['candidates']
    assert actual['search_space']['census_passes']==1 and value==before
    assert actual['candidates_evaluated']==6
    assert actual['selected']['api_kwargs']['reduce']=='significant'
    assert actual['selected']['api_kwargs']['significance_threshold']==.05
    assert actual['selected']['api_kwargs']['trait_block'] in (2,5)
    assert 'variant_devices' not in actual['selected']['api_kwargs']


def test_colliding_names_cannot_hide_worst_scenario_or_change_selected_candidate(input_path):
    choices=[candidate(input_path,width=2,count=1,block_bytes=None),
             candidate(input_path,width=5,count=1,block_bytes=None)]
    args=options();args['occupancy_scenarios']={'dense:x':'dense','dense':'empty'}
    args['host_scenarios']={'y':dict(host_serial_fraction=0.),'x:y':dict(host_serial_fraction=1.)}
    def costs(choice,prices,*,occupancy,host_serial_fraction,**kwargs):
        # Old string concatenation overwrote the first 100-second combination
        # with the last 1-second combination, incorrectly selecting candidate 0.
        seconds=(100. if occupancy=='dense' and host_serial_fraction==0. else 1.) if choice['trait_block']==2 else 50.
        return dict(estimated_tile_seconds=seconds,unpriced_terms=[])
    with patch('torchgwas.significant_host_plan.significant_host_runtime',side_effect=costs) as run:
        result=detailed_significant_host_plan(choices,bank(),**args)
    assert run.call_count==8 and result['selected']['candidate_index']==1
    bad=next(row for row in result['candidates'] if row['candidate_index']==0)
    assert bad['worst_supplied_scenario_seconds']==100.
    assert len(bad['scenarios'])==4
    assert {(row['occupancy'],row['host']) for row in bad['scenarios']}=={
        ('dense:x','y'),('dense:x','x:y'),('dense','y'),('dense','x:y')}
    assert json.loads(json.dumps(result))['selected']['candidate_index']==1


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),-1.,0.,True,None])
def test_every_scenario_score_must_be_finite_positive(input_path,bad):
    choice=candidate(input_path,width=2,count=1,block_bytes=None)
    def costs(*args,occupancy,**kwargs):
        return dict(estimated_tile_seconds=1. if occupancy=='empty' else bad,unpriced_terms=[])
    with patch('torchgwas.significant_host_plan.significant_host_runtime',side_effect=costs):
        with pytest.raises(ValueError,match='finite positive'):
            detailed_significant_host_plan([choice],bank(),**options())


def test_repeated_tile_work_is_reused_only_within_each_plan_without_changing_scores(input_path):
    choice=candidate(input_path,width=2,count=1,block_bytes=None)
    args=options();args['host_scenarios']['serial']=dict(host_serial_fraction=1.)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        with patch('torchgwas.trait_tiling_model.torch_scan_work',wraps=torch_scan_work) as collect:
            cached=detailed_significant_host_plan([choice],bank(),**args)
            assert collect.call_count==2  # K=2 and the K=1 tail, shared across four scenarios.
            repeated=detailed_significant_host_plan([choice],bank(),**args)
            assert collect.call_count==4  # No reuse of mutable inputs across jobs.
        with patch('torchgwas.significant_host_model._cached_scan_work',wraps=torch_scan_work) as uncached:
            fresh=detailed_significant_host_plan([choice],bank(),**args)
            assert uncached.call_count==12
    assert cached==repeated==fresh


@pytest.mark.parametrize('fault',['occupancy_name','host_name','host_fraction','host_policy','capacity','threshold','reserve'])
def test_invalid_scenario_or_capacity_is_refused_before_graph_work(input_path,fault):
    args=options();choices=[candidate(input_path,width=2,count=1,block_bytes=None)]
    if fault=='occupancy_name':args['occupancy_scenarios']={'':'dense'}
    if fault=='host_name':args['host_scenarios']={1:dict(host_serial_fraction=0.)}
    if fault=='host_fraction':args['host_scenarios']['fluid']['host_serial_fraction']=True
    if fault=='host_policy':args['host_scenarios']['fluid']['host_serial_policy']='unpriced'
    if fault=='capacity':args['device_memory_bytes']['cuda:0']=-1
    if fault=='threshold':args['significance_threshold']='invalid'
    if fault=='reserve':args['host_reserve_bytes']=-1
    with patch('torchgwas.significant_host_plan.significant_host_runtime') as run:
        with pytest.raises(ValueError):detailed_significant_host_plan(choices,bank(),**args)
        run.assert_not_called()


def test_infeasible_whole_panel_does_not_need_geometry_capture(tmp_path):
    value=specification(tmp_path);value['workload']['traits']=33
    value['contexts']=value['contexts'][:1];value['bounds']=dict(chunks=[4],trait_blocks=[2,33])
    space=proposed(value);small,large=space['candidates']
    lo=significant_host_memory(small)['device_bytes']['cuda:0']
    hi=significant_host_memory(large)['device_bytes']['cuda:0']
    assert lo<hi
    value['joint']['device_memory_bytes']['cuda:0']=(lo+hi)//2
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        plan=bounded_significant_host_plan(**value)
    assert plan['selected']['trait_block']==2 and plan['candidates_feasible']==1
    assert plan['rejected'][0]['reason']=='device_memory'


@pytest.mark.parametrize('fault',['partition','coalescing','budget'])
def test_unsupported_output_and_search_bounds_fail_before_census(tmp_path,fault):
    value=specification(tmp_path)
    if fault=='partition':value['bounds']['partition_axes']=['variant']
    if fault=='coalescing':value['output']['block_bytes']=1024
    if fault=='budget':value['joint']['max_candidates']=1
    with patch('torchgwas.trait_candidate_space.read_header') as header:
        with pytest.raises(ValueError):bounded_significant_host_plan(**value)
        header.assert_not_called()


def test_json_cli_keeps_reduction_and_scenario_names(tmp_path,capsys):
    from torchgwas.autotune import main
    value=specification(tmp_path);value['model']='detailed_significant_host_space'
    path=tmp_path/'space.json';path.write_text(json.dumps(value))
    with patch('sys.argv',['autotune',str(path)]),patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        main()
    result=json.loads(capsys.readouterr().out)
    assert result['selected']['api_kwargs']['reduce']=='significant'
    assert len(result['selected']['scenarios'])==2
