import copy
import random

import pytest

from torchgwas.execution_graph import ExecutionGraph
from torchgwas.remaining_work import remaining_work_bounds


def chain_graph():
    graph=ExecutionGraph();graph.capacities={'cpu':.5,'io':2.}
    for i in range(5):
        deps=[] if not i else [f'done:{i-1}']
        start=graph.add(f'start:{i}',after=deps)
        ready=graph.add(f'io:{i}',.2,[start],{'io':3.})
        graph.resource_waits[f'wait:{i}']=dict(after=start,until=ready,resources={'cpu':.7})
        graph.add(f'done:{i}',.1,[ready],{'cpu':1.})
    return graph


def test_wait_chain_proof_tightens_ceiling_and_every_continuation_is_inside():
    graph=chain_graph();solved=graph.solve()
    groups=[[f'wait:{i}' for i in range(5)]]
    bound=remaining_work_bounds(graph,unissued_nodes=['start:0'],wait_chains=groups)
    loose=remaining_work_bounds(graph,unissued_nodes=['start:0'])
    assert bound['upper_seconds']<loose['upper_seconds']
    assert bound['wait_groups']==1 and bound['unchained_waits']==0
    for fraction in [0.,.01,.2,.49,.8,.99]:
        at=solved['seconds']*fraction;state=graph.checkpoint(at)
        unissued=[name for name in graph.nodes if name.startswith('start:') and name not in state['start']]
        result=remaining_work_bounds(graph,unissued_nodes=unissued,wait_chains=groups)
        remaining=graph.resume(state)['remaining_seconds']
        assert result['lower_seconds']<=remaining+1e-10<=result['upper_seconds']+1e-10


@pytest.mark.parametrize('seed',range(25))
def test_random_shared_resource_continuations_respect_bounds(seed):
    rng=random.Random(seed);graph=ExecutionGraph();graph.capacities={'cpu':.3+rng.random(),'dram':.5+rng.random()}
    graph.add('root')
    for i in range(25):
        names=list(graph.nodes)
        deps=rng.sample(names,min(len(names),rng.randint(0,4)))
        graph.add(str(i),rng.random(),deps,{key:2*rng.random() for key in graph.capacities})
    graph.resource_waits['external_wait']=dict(after='root',until='24',resources={'cpu':.7})
    original=copy.deepcopy(graph.__dict__);solved=graph.solve()
    for fraction in [0.,.2,.8]:
        state=graph.checkpoint(fraction*solved['seconds'])
        nodes=[name for name in graph.nodes if name not in state['start']]
        result=remaining_work_bounds(graph,unissued_nodes=nodes[::3])
        remaining=graph.resume(state)['remaining_seconds']
        assert result['lower_seconds']<=remaining+1e-8
        assert remaining<=result['upper_seconds']+1e-8
    assert graph.__dict__==original


def test_optional_delay_is_excluded_from_floor_and_included_in_ceiling():
    graph=ExecutionGraph();graph.add('a',1.);graph.add('b',2.)
    graph.add_wait_delay('maybe',100.,'a','b')
    graph.add('finish',1.,['maybe'])
    bound=remaining_work_bounds(graph,unissued_nodes=['a','b'])
    assert bound['lower_seconds']==3. and bound['upper_seconds']==104.
    assert bound['upper_seconds']>=graph.solve()['seconds']


@pytest.mark.parametrize('damage',['order','duplicate','unknown','budget','node_budget','cycle','capacity','demand','unissued'])
def test_inconsistent_or_unbounded_requests_fail(damage):
    graph=chain_graph();args=dict(unissued_nodes=['start:0'],wait_chains=[[f'wait:{i}' for i in range(5)]])
    if damage=='order':args['wait_chains'][0].reverse()
    elif damage=='duplicate':args['wait_chains'].append(['wait:0'])
    elif damage=='unknown':args['wait_chains']=[['other']]
    elif damage=='budget':args['max_reachability_visits']=1
    elif damage=='node_budget':args['max_nodes']=1
    elif damage=='cycle':graph.nodes['start:0']=(0.,('done:4',))
    elif damage=='capacity':graph.capacities['cpu']=0.
    elif damage=='demand':graph.demands['done:0']['cpu']=float('nan')
    else:args['unissued_nodes']=['missing']
    with pytest.raises(ValueError):remaining_work_bounds(graph,**args)
