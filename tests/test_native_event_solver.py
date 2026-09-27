"""Native scheduling must reproduce the independent Python event loop exactly."""
import copy
import random
from unittest.mock import patch

import pytest

from torchgwas.execution_graph import ExecutionGraph
from torchgwas import native_event_solver as native
from torchgwas.detailed_calibration import source_identity
from test_streamed_graphs import graph,full

pytestmark=pytest.mark.skipif(native.library() is None,reason='Build native CPU event solver first')


def compare(g):
    before=copy.deepcopy(g.__dict__)
    expected=g._solve_shared_python()
    actual=g._solve_shared()
    assert actual==expected
    assert g.__dict__==before
    return actual


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
@pytest.mark.parametrize('seed',range(24))
def test_complete_and_streamed_event_traces_match_python(seed,policy):
    chains=[[],[],[]]
    for index in range(15):
        g=graph(seed*37+index,index%3)
        if policy!='fluid':g=g.with_serial_sections(policy)
        chains[index%3].append((index,g))
    expanded=full(chains);expected=expanded._solve_shared_python()
    assert expanded._solve_shared()==expected
    actual=ExecutionGraph().solve_chains(chains,shared_tokens={'exclusive:host_serial'},trace=True)
    for key,value in expected.items():assert actual[key]==value,key
    reference=ExecutionGraph()._solve_shared_python(_chains=chains,_shared_tokens={'exclusive:host_serial'},_trace=True)
    assert actual==reference
    compact=ExecutionGraph().solve_chains(chains,shared_tokens={'exclusive:host_serial'})
    reference=ExecutionGraph()._solve_shared_python(_chains=chains,_shared_tokens={'exclusive:host_serial'})
    assert compact==reference


@pytest.mark.parametrize('seed',range(36))
def test_varied_dags_zero_demands_ties_and_token_sections(seed):
    rng=random.Random(seed);g=ExecutionGraph();g.capacities={'cpu':1.3,'dram':3.,'gpu':.9}
    g.token_capacities={'critical':2}
    for index in range(45):
        deps=rng.sample(list(g.nodes),min(len(g.nodes),rng.randrange(4)))
        resources={name:rng.choice([0.,.01,.3,1.,7.]) for name in g.capacities if rng.random()<.6}
        if rng.random()<.2:resources['zero-without-capacity']=0.
        name=str(index);g.add(name,rng.choice([0.,1e-16,.0001,.002,.2]),deps,resources)
        if rng.random()<.5:g.token_actions[name]={'acquire':{'critical':1},'release_finish':{'critical':1}}
    compare(g)


def test_future_shared_token_and_reused_templates_keep_global_state():
    a=graph(9).with_serial_sections('held-last');b=graph(3).with_serial_sections('held-last')
    chains=[[(i,a) for i in range(0,101,2)],[(i,b) for i in range(1,101,2)]]
    for trace in [False,True]:
        expected=ExecutionGraph()._solve_shared_python(_chains=chains,_shared_tokens={'exclusive:host_serial'},_trace=trace)
        actual=ExecutionGraph().solve_chains(chains,shared_tokens={'exclusive:host_serial'},trace=trace)
        assert actual==expected
    assert actual['peak_active_nodes']==len(a.nodes)+len(b.nodes)+2


@pytest.mark.parametrize('case',['cycle','missing_dependency','negative_demand','missing_capacity','invalid_token',
    'release','condition_endpoint','wait_endpoint','fifo_deadlock'])
def test_invalid_graphs_fail_with_same_exception_class(case):
    g=graph(1)
    if case=='cycle':g.nodes['read']=(1.,('fsync',))
    elif case=='missing_dependency':g.nodes['read']=(1.,('absent',))
    elif case=='negative_demand':g.demands['read']['cpu']=-1.
    elif case=='missing_capacity':del g.capacities['storage']
    elif case=='invalid_token':g.token_actions['read']={'acquire':{'missing':1}}
    elif case=='release':g.token_actions['read']={'release_finish':{'writer_pool':1}}
    elif case=='condition_endpoint':g.conditional_delays['resume']['ready']='absent'
    elif case=='wait_endpoint':g.resource_waits['spin']['until']='absent'
    else:g.fifo_enqueues.clear()
    with pytest.raises(Exception) as expected:g._solve_shared_python()
    with pytest.raises(type(expected.value)):g._solve_shared()


def test_shared_token_may_be_declared_in_another_part():
    # This unusual cross-template declaration is intentionally handled by the
    # Python fallback, preserving the accepted graph API beyond model inputs.
    a=ExecutionGraph();a.add('a',1.);a.token_capacities['shared']=1
    b=ExecutionGraph();b.add('b',2.);b.token_actions['b']={'acquire':{'shared':1},'release_finish':{'shared':1}}
    chains=[[(0,a),(1,b)]]
    expected=ExecutionGraph()._solve_shared_python(_chains=chains,_shared_tokens={'shared'})
    actual=ExecutionGraph().solve_chains(chains,shared_tokens={'shared'})
    assert actual==expected


def test_unusual_keys_keep_python_semantics():
    g=ExecutionGraph();g.capacities={'cpu':1.};g.add(1,1.,resources={'cpu':.5});g.add(2,1.,[1])
    compare(g)
    g=ExecutionGraph();g.capacities={'cpu':1.};g.add('a\0b',1.,resources={'cpu':.5});compare(g)
    chains=[[(0,graph(1))]]
    expected=ExecutionGraph()._solve_shared_python(_chains=chains,_prefix='x\0')
    assert ExecutionGraph().solve_chains(chains,prefix='x\0')==expected


def test_empty_schedules_and_empty_parts():
    compare(ExecutionGraph())
    chains=[[],[(0,ExecutionGraph()),(1,ExecutionGraph())]]
    assert ExecutionGraph().solve_chains(chains)==ExecutionGraph()._solve_shared_python(_chains=chains)
    assert ExecutionGraph().solve_chains([])==ExecutionGraph()._solve_shared_python(_chains=[])


def test_missing_binary_preserves_reference_and_cpp_is_bound():
    g=graph(3)
    with patch.object(native,'library',return_value=None):assert g._solve_shared()==g._solve_shared_python()
    assert '_event_solver.cpp' in source_identity()
    assert native.library().BUILD_KEY==native._identity()[2]

def test_streamed_nonstring_wait_key_preserves_composition_error():
    g=graph(1)
    g.resource_waits[7]=g.resource_waits.pop('spin')
    chains=[[(0,g)]]
    with pytest.raises(TypeError):ExecutionGraph()._solve_shared_python(_chains=chains)
    with pytest.raises(TypeError):ExecutionGraph().solve_chains(chains)
