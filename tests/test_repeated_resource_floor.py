"""Repeated mandatory work scales without replicating chunks or fixed waits."""
import pytest
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.resource_balance import repeated_resource_floor


def test_billions_of_chunks_share_one_cpu_and_storage_budget():
    full=ExecutionGraph();full.capacities={'cpu':2.,'input':100.}
    full.add('decode',3.,resources={'cpu':1.,'input':50.})
    tail=ExecutionGraph();tail.capacities=dict(full.capacities)
    tail.add('decode',1.,resources={'cpu':1.,'input':50.})
    result=repeated_resource_floor([dict(graph=full,copies=10**9),dict(graph=tail,copies=1)])
    assert result['template_nodes']==2 and result['copies']==10**9+1
    assert result['resources']['cpu']['work']==3*10**9+1
    assert result['resources']['input']['work']==150*10**9+50
    assert result['lower_bound_seconds']==(3*10**9+1)/2
    assert set(result['largest_resource_bounds'])=={'cpu','input'}
    assert len(full.nodes)==len(tail.nodes)==1


def test_contention_dependent_waits_are_not_multiplied_as_mandatory_work():
    g=ExecutionGraph();g.capacities={'cpu':1.,'gpu':1.}
    g.add('start');g.add('kernel',5.,['start'],resources={'gpu':1.})
    g.add_wait_delay('resume',10.,'start','kernel');g.demands['resume']={'cpu':1.}
    g.resource_waits['spin']=dict(after='start',until='kernel',resources={'cpu':1.})
    r=repeated_resource_floor([dict(graph=g,copies=100)])
    assert r['resources']['cpu']['work']==0.
    assert r['resources']['gpu']['work']==500.
    assert r['single_copy_dependency_lower_bound_seconds']==5.
    assert r['templates'][0]['omitted_conditional_nodes']==r['templates'][0]['omitted_active_waits']==1


def test_independent_copy_critical_paths_are_not_assumed_serial():
    g=ExecutionGraph();g.add('a',2.);g.add('b',3.,['a'])
    r=repeated_resource_floor([dict(graph=g,copies=1000)])
    assert r['lower_bound_seconds']==5. and r['resource_lower_bound_seconds']==0.


def test_floor_is_below_complete_shared_resource_schedule():
    part=ExecutionGraph();part.capacities={'cpu':1.,'gpu':1.}
    part.add('decode',3.,resources={'cpu':1.});part.add('kernel',5.,['decode'],resources={'gpu':1.})
    whole=ExecutionGraph()
    for i in range(7):whole.compose(part,f'{i}:')
    r=repeated_resource_floor([dict(graph=part,copies=7)])
    assert r['lower_bound_seconds']<=whole.solve()['seconds']
    assert r['lower_bound_seconds']==35.


def test_individual_nominal_service_also_respects_capacity():
    g=ExecutionGraph();g.capacities={'cpu':1.}
    g.add('a',2.,resources={'cpu':4.});g.add('b',3.,['a'])
    assert repeated_resource_floor([dict(graph=g,copies=1)])['lower_bound_seconds']==11.


@pytest.mark.parametrize('copies',[0,-1,True,1.5])
def test_noninteger_multiplicities_are_refused(copies):
    with pytest.raises(ValueError):repeated_resource_floor([dict(graph=ExecutionGraph(),copies=copies)])


def test_conflicts_missing_capacities_and_expansion_budgets_are_refused():
    a=ExecutionGraph();a.capacities={'cpu':1.};a.add('x',1.,resources={'cpu':1.})
    b=ExecutionGraph();b.capacities={'cpu':2.}
    with pytest.raises(ValueError,match='Conflicting'):repeated_resource_floor([dict(graph=a,copies=1),dict(graph=b,copies=1)])
    with pytest.raises(ValueError,match='bounded'):repeated_resource_floor([dict(graph=a,copies=1)]*2,max_templates=1)
    a.add('y',1.)
    with pytest.raises(ValueError,match='max_template_nodes'):repeated_resource_floor([dict(graph=a,copies=1)],max_template_nodes=1)
    a.demands['x']={'missing':1.}
    with pytest.raises(ValueError,match='Missing'):repeated_resource_floor([dict(graph=a,copies=1)])
