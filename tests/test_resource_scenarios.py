import pytest
from torchgwas.resource_scenarios import primitive_cpu_scenarios


def test_cpu_scenarios_use_independent_repetitions_without_workload_data():
    rows=[dict(primitive=name,cpu_seconds_per_call=value) for name,values in
          [('copy',[3.,1.,2.]),('view',[.5,.7,.6])] for value in values]
    result=primitive_cpu_scenarios(rows)
    assert result['nominal']=={'copy':2.,'view':.6}
    assert result['low_dispatch']=={'copy':1.,'view':.5}
    assert result['high_dispatch']=={'copy':3.,'view':.7}
    assert rows[0]['cpu_seconds_per_call']==3.
    with pytest.raises(ValueError):primitive_cpu_scenarios([])
    for value in [-1.,float('nan'),float('inf')]:
        with pytest.raises(ValueError):primitive_cpu_scenarios([dict(primitive='x',cpu_seconds_per_call=value)])


def test_scenario_composition_keeps_resources_separate_and_reports_worst():
    from unittest.mock import patch
    from torchgwas.resource_scenarios import torch_runtime_scenarios
    profile={'host_primitives':{'copy':99.},'chunk_markers':8}
    def runtime(data,candidate):
        return {'estimated_seconds':candidate['host_primitives']['copy']+data['fixed'],
                'unpriced_terms':['allocation']}
    with patch('torchgwas.mechanistic_torch.torch_runtime',side_effect=runtime):
        result=torch_runtime_scenarios({'fixed':5.},profile,
            {'nominal':{'copy':2.},'low_dispatch':{'copy':1.},'high_dispatch':{'copy':3.}})
    assert result['nominal_seconds']==7.
    assert result['scenario_min_seconds']==6.
    assert result['scenario_max_seconds']==8.
    assert result['worst_scenario']=='high_dispatch'
    assert result['unpriced_terms']==['allocation']
    assert profile['host_primitives']['copy']==99.
    assert not result['prediction_complete']
    with pytest.raises(ValueError,match='nominal'):
        torch_runtime_scenarios({},profile,{})


def test_comparison_preserves_matched_scenarios_and_distinguishes_objectives():
    from torchgwas.resource_scenarios import compare_runtime_scenarios
    def report(a, b):
        return dict(prediction_complete=False, unpriced_terms=['allocation'],
                    estimates={s: dict(estimated_seconds=t) for s,t in [('idle',a),('busy',b)]})
    # Absolute worst time chooses A, but relative worst regret chooses B.
    result=compare_runtime_scenarios({'A':report(2,10),'B':report(1,11)})
    assert result['minimax_time_candidates']==['A']
    assert result['minimax_regret_candidates']==['B']
    assert result['common_winners']==[]
    assert result['worst_scenario_regret']['A']==1
    assert result['worst_scenario_regret']['B']==pytest.approx(.1)
    assert not result['automatic_selection_ready']
    assert not result['prediction_complete']
    assert result['unpriced_terms']==['allocation']


def test_comparison_keeps_ties_and_refuses_unmatched_or_invalid_times():
    from torchgwas.resource_scenarios import compare_runtime_scenarios
    def report(value, scenario='idle'):
        return dict(estimates={scenario:dict(estimated_seconds=value)})
    result=compare_runtime_scenarios({'A':report(1),'B':report(1)})
    assert result['common_winners']==['A','B']
    assert result['minimax_regret_candidates']==['A','B']
    with pytest.raises(ValueError):compare_runtime_scenarios({})
    with pytest.raises(ValueError,match='identical'):
        compare_runtime_scenarios({'A':report(1),'B':report(1,'busy')})
    for value in [None,0,-1,float('nan'),float('inf'),True]:
        with pytest.raises(ValueError,match='finite and positive'):
            compare_runtime_scenarios({'A':report(value)})
