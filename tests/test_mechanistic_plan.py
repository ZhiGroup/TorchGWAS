"""Selection policy, executable geometry and hard resource constraints."""
import copy
from unittest.mock import patch

import numpy as np
import pytest

from test_pgen_native_reader import write_pgen
from test_mechanistic_shapes import fixture, component
from torchgwas.autotune import detailed_joint_plan
from torchgwas.linear import multigpu_variant_ranges
from torchgwas.pgen_work_census import census
from torchgwas.native_control_work import native_control_work


@pytest.fixture
def search(tmp_path):
    path = tmp_path/'input.pgen'
    write_pgen(path, np.arange(10*32, dtype=np.uint8).reshape(10,32)%4)
    candidates = []
    for count in (1, 2):
        shards = []
        for index, span in enumerate(multigpu_variant_ranges(10, 4, count)):
            data, profile = fixture()
            data['markers'] = span[1]-span[0]
            data['encoded'] = census(path, 4, span, include_chunks=True)
            profile.update(chunk_markers=4, decode_workers=2//count, validate_range=True,event_wait_cpu_fraction=1.,
                result_ownership='borrowed',
                control_primitives={name:1e-6 for counts in native_control_work().values() for name in counts},
                result_finish_service=dict(cpu_seconds=1e-6, serial_cpu_seconds=0., baseline_copy_bytes=0,
                    replaces_fixed_finish_and_tensor_conversion=True, includes_ready_cuda_event=False),
                acknowledgement_service=dict(create_cpu_seconds=1e-6,publish_cpu_seconds=1e-6,
                    receive_cpu_seconds=1e-6,wakeup_seconds=0.))
            profile['kernel_geometry'] = [dict(N=32,B=b,K=3,C=2,validate_range=True,kernels=[]) for b in (4,2)]
            shards.append(dict(device=f'cuda:{index}',data=data,profile=profile))
        candidates.append(dict(shards=shards, shared_capacities=dict(cpu=2.,dram=1e9,input=1e8),
            queue_service=dict(put_cpu_seconds=1e-6,get_cpu_seconds=1e-6,cpu_fraction=1.)))
    options = dict(host_scenarios={'low':dict(host_serial_fraction=0.,host_serial_policy='fluid'),
                                  'high':dict(host_serial_fraction=1.,host_serial_policy='held-first')},
        cpu_workers=2, host_memory_bytes=2**30,device_memory_bytes={'cuda:0':2**30,'cuda:1':2**30})
    return candidates, options


def run_search(search, **options):
    candidates, defaults = search
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return detailed_joint_plan(candidates, **(defaults | options))


def test_real_graph_path_preserves_inputs_and_executor_geometry(search):
    before = copy.deepcopy(search)
    result = run_search(search)
    assert search == before
    assert result['candidates_feasible'] == 2
    assert result['runtime_prediction_validated'] is False
    for row in result['shortlist']:
        assert row['scan_kwargs']['borrow_results'] is True
        assert row['scan_kwargs']['reader_workers'] == 2
        assert row['scan_kwargs']['chunk_size'] == 4
        assert row['scan_kwargs']['ordered'] is False
        assert row['source_kwargs']['pgen_decode_workers'] == 2
        assert row['required_environment']['TORCHGWAS_PGEN_PACKED'] == '0'
        assert set(row['scenarios']) == {'low','high'}
        assert row['memory']['host_bytes'] > 0
    assert result['predicted_selection_penalty_fraction'] == 0


def test_throughput_default_does_not_discard_small_predicted_gain(search):
    def predict(shards, *args, host_serial_fraction, **kwargs):
        seconds = (1.01 if len(shards)==1 else 1.) + .1*host_serial_fraction
        return dict(estimated_scan_seconds=seconds, unpriced_terms=[])
    candidates, options = search
    with patch('torchgwas.mechanistic_plan.torch_multigpu_scan_runtime',side_effect=predict):
        fast = detailed_joint_plan(candidates, **options)
        frugal = detailed_joint_plan(candidates, **options, max_slowdown_fraction=.02)
    assert len(fast['selected']['devices']) == 2
    assert len(frugal['selected']['devices']) == 1
    assert frugal['predicted_selection_penalty_fraction'] == pytest.approx(1.11/1.1-1)


def test_minimax_uses_every_scenario(search):
    def predict(shards, *args, host_serial_fraction, **kwargs):
        seconds = 1.1 if len(shards)==1 else (1. if host_serial_fraction==0 else 1.3)
        return dict(estimated_scan_seconds=seconds, unpriced_terms=['test unknown'])
    candidates, options = search
    with patch('torchgwas.mechanistic_plan.torch_multigpu_scan_runtime',side_effect=predict):
        result = detailed_joint_plan(candidates, **options)
    assert len(result['selected']['devices']) == 1
    assert result['unresolved_timing_terms'] == ['test unknown']


def test_host_device_and_reader_budgets_exclude_candidates_before_pricing(search):
    result = run_search(search)
    one = next(row for row in result['shortlist'] if len(row['devices'])==1)
    limited = run_search(search, host_memory_bytes=one['memory']['host_bytes'])
    assert limited['candidates_feasible'] == 1
    assert limited['rejected'][0]['reason'] == 'host_memory'
    limited = run_search(search,device_memory_bytes={'cuda:0':2**30,'cuda:1':1})
    assert limited['rejected'][0]['reason'] == 'device_memory'
    with pytest.raises(ValueError, match='shared_reader_budget'):
        run_search(search,cpu_workers=1)
    with pytest.raises(ValueError, match='No feasible'):
        run_search(search,device_reserve_bytes=2**30)


@pytest.mark.parametrize('mutation,match', [
    (lambda c: c.update(observed_seconds=1.), 'Unknown candidate'),
    (lambda c: c['shared_capacities'].update(input=1.), 'Shared capacity'),
    (lambda c: c['shards'][0]['data']['encoded'].update(variant_range=[1,10]), 'Census ranges'),
    (lambda c: c['shards'][0]['data']['encoded'].pop('chunks'), 'per-chunk'),
    (lambda c: c['shards'][0]['profile'].update(result_ownership='owned'), 'borrowed'),
    (lambda c: c['shards'][0]['profile'].update(validate_range=False), 'range validation'),
    (lambda c: c['shards'][0]['profile'].update(event_wait_cpu_fraction=.5), 'completion events'),
    (lambda c: c['shards'][0].update(device='cuda:00'), 'canonical'),
])
def test_refuses_nonexecutable_or_mismatched_inputs(search, mutation, match):
    mutation(search[0][0])
    with pytest.raises(ValueError,match=match):
        run_search(search)


def test_uneven_worker_budget_must_follow_executor_assignment(search):
    shards=search[0][1]['shards']
    shards[1]['profile']['decode_workers']=2
    with pytest.raises(ValueError,match='Reader allocation'):
        run_search(search,cpu_workers=3)


def test_missing_geometry_is_reported_without_coarse_fallback(search):
    search[0][0]['shards'][0]['profile']['kernel_geometry']=[]
    result=run_search(search)
    assert result['candidates_feasible']==1
    assert result['rejected'][0]['reason']=='unsupported_model_context'
    assert 'geometry' in result['rejected'][0]['detail']


def test_explicit_queue_depth_reaches_detailed_schedule(search):
    search[0][1]['result_queue_depth']=1
    result=run_search(search)
    two=next(row for row in result['shortlist'] if len(row['devices'])==2)
    assert two['scan_kwargs']['result_queue_depth']==1
    assert all(r['resource_policy']['result_queue_capacity']==2 for r in two['scenarios'].values())


@pytest.mark.parametrize('kwargs,match', [
    ({'max_candidates':1},'max_candidates'),
    ({'max_scenario_evaluations':3},'max_scenario_evaluations'),
    ({'host_scenarios':{}},'scenarios'),
    ({'max_slowdown_fraction':True},'max_slowdown'),
    ({'host_memory_bytes':True},'integer'),
])
def test_search_is_bounded_and_requires_explicit_context(search, kwargs, match):
    with pytest.raises(ValueError,match=match):
        run_search(search,**kwargs)
