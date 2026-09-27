"""Conservation under concurrency, dependency stalls, and conditional waits."""
import pytest
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.resource_balance import resource_balance


def test_two_gpu_workers_share_one_cpu_budget():
    g = ExecutionGraph(); g.capacities = {'cpu': 1., 'dram': 100.}
    for name in ['gpu0_host', 'gpu1_host']:
        g.add(name, 4., resources={'cpu': 1., 'dram': 20.})
    result = resource_balance(g, g.solve())
    assert result['scheduled_seconds'] == 8.
    assert result['resources']['cpu']['work'] == 8.
    assert result['resources']['dram']['work'] == 160.
    assert result['resources']['dram']['mean_capacity_fraction'] == .2
    assert result['largest_resource_bounds'] == ['cpu']
    assert result['resource_lower_bound_seconds'] == 8.


def test_dependency_time_is_not_mistaken_for_consumed_capacity():
    g = ExecutionGraph(); g.capacities = {'cpu': 1., 'transfer': 1.}
    g.add('prepare', 4., resources={'cpu': 1.})
    g.add('copy', 4., ['prepare'], resources={'transfer': 1.})
    result = resource_balance(g, g.solve())
    assert result['scheduled_seconds'] == 8.
    assert result['resource_lower_bound_seconds'] == 4.
    assert set(result['largest_resource_bounds']) == {'cpu', 'transfer'}
    assert all(row['mean_capacity_fraction'] == .5 for row in result['resources'].values())


@pytest.mark.parametrize('blocked', [False, True])
def test_conditional_handoff_only_charged_when_executed(blocked):
    g = ExecutionGraph(); g.capacities = {'cpu': 1.}
    g.add('attempt')
    g.add('ready', 5. if blocked else 0.)
    g.add('handoff', 2., ['attempt', 'ready'], resources={'cpu': 1.})
    g.conditional_delays['handoff'] = {'attempt': 'attempt', 'ready': 'ready'}
    result = resource_balance(g, g.solve())
    assert result['resources']['cpu']['work'] == (2. if blocked else 0.)
    assert result['scheduled_seconds'] == (7. if blocked else 0.)


def test_active_cpu_wait_is_counted_at_its_shared_rate():
    g = ExecutionGraph(); g.capacities = {'cpu': 1., 'gpu': 1.}
    g.add('start')
    g.add('kernel', 5., ['start'], resources={'gpu': 1.})
    g.add('host', 5., ['start'], resources={'cpu': 1.})
    g.resource_waits['spin'] = {'after': 'start', 'until': 'kernel', 'resources': {'cpu': 1.}}
    solved = g.solve()
    assert solved['wait_resource_seconds']['cpu'] == 2.5
    result = resource_balance(g, solved)
    assert result['resources']['cpu']['work'] == 7.5
    assert result['resource_lower_bound_seconds'] == result['scheduled_seconds'] == 7.5


def test_plain_dependency_graph_has_no_invented_resource_pressure():
    g = ExecutionGraph(); g.add('unpriced_resource', 3.)
    result = resource_balance(g, g.solve())
    assert result['resources'] == {} and result['largest_resource_bounds'] == []
    assert result['resource_lower_bound_seconds'] == 0.
    assert result['scheduled_seconds'] == 3.


def test_partial_or_wrong_solution_is_refused():
    g = ExecutionGraph(); g.capacities = {'cpu': 1.}
    g.add('work', 3., resources={'cpu': 1.})
    solved = g.solve()
    with pytest.raises(ValueError, match='complete'):
        resource_balance(g, dict(solved, start={}))
    with pytest.raises(ValueError, match='exceeds capacity'):
        resource_balance(g, dict(solved, seconds=1.))
