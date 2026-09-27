"""Lazy graph admission must preserve every event in the expanded schedule."""
import copy
import random

import pytest

from torchgwas.execution_graph import ExecutionGraph
from torchgwas.mechanistic_torch import handoff_summary


def graph(seed,device=0):
    rng=random.Random(seed)
    g=ExecutionGraph();g.capacities={'cpu':1.3,'dram':2.,'host_serial':1.,'storage':1.}
    def duration():return rng.choice([.001,.003,.007,.013])
    g.add('read',duration(),resources={'storage':.8,'cpu':.2})
    g.add('api',duration(),resources={'cpu':.7,'host_serial':.3})
    g.add('gpu',duration(),['read','api'])
    g.resource_waits['spin']=dict(after='api',until='gpu',resources={'cpu':.4})
    g.add_wait_delay('resume',.002,'read','gpu')
    g.add('copy',duration(),['resume'],resources={'cpu':.6,'dram':1.5,'host_serial':.1})
    g.token_capacities['writer_pool']=1
    for i in range(2):
        enqueue=g.add('enqueue:'+str(i),0.,['copy'] if i==0 else ['enqueue:0'])
        write=g.add('write:'+str(i),duration(),['copy'],resources={'cpu':.2,'storage':.7})
        g.token_actions[write]={'acquire':{'writer_pool':1},'release_finish':{'writer_pool':1}}
        g.fifo_enqueues[enqueue]=('queue',write);g.fifo_dequeues[write]='queue'
    g.add('fsync',duration(),['write:0','write:1'])
    return g


def full(chains,caps=None):
    g=ExecutionGraph();g.capacities=dict(caps or {})
    rows=sorted((order,stream,tile) for stream,chain in enumerate(chains) for order,tile in chain)
    done={}
    for order,stream,tile in rows:
        done[stream]=g.compose(tile,f'part{order}:',[done[stream]] if stream in done else [],
                               shared_tokens={'exclusive:host_serial'})
    return g


def check(chains):
    expected=full(chains)._solve_shared()
    actual=ExecutionGraph().solve_chains(chains,shared_tokens={'exclusive:host_serial'},trace=True)
    for key,value in expected.items():assert actual[key]==value,key
    compact=ExecutionGraph().solve_chains(chains,shared_tokens={'exclusive:host_serial'})
    for key in ('seconds','resource_event_steps','wait_resource_seconds','resource_policy'):
        assert compact[key]==expected[key]
    summary=handoff_summary(expected);summary.pop('scope')
    assert compact['conditional_summary']==summary
    assert compact['start']==compact['end']==compact['conditional_delays']=={}
    assert compact['scheduled_nodes']==sum(len(g.nodes)+1 for c in chains for _,g in c)
    assert compact['peak_active_nodes']<=sum(max((len(g.nodes)+1 for _,g in c),default=0) for c in chains)
    return actual


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
def test_exact_full_event_traces_with_waits_fifo_and_shared_serial_token(policy):
    for seed in range(12):
        chains=[[],[],[]]
        for i in range(13):
            g=graph(seed*37+i,i%3)
            if policy!='fluid':g=g.with_serial_sections(policy)
            chains[i%3].append((i,g))
        before=copy.deepcopy([[g.__dict__ for _,g in c] for c in chains])
        check(chains)
        assert [[g.__dict__ for _,g in c] for c in chains]==before


def test_repeated_templates_do_not_accumulate_live_nodes_or_timestamps():
    a,b=graph(3).with_serial_sections('held-first'),graph(7).with_serial_sections('held-first')
    chains=[[(i,a) for i in range(0,101,2)],[(i,b) for i in range(1,101,2)]]
    result=check(chains)
    assert result['peak_active_nodes']==len(a.nodes)+len(b.nodes)+2


def test_later_global_wait_order_can_be_admitted_before_earlier_order():
    slow=graph(9);slow.nodes['fsync']=(100.,slow.nodes['fsync'][1])
    chains=[[(0,slow),(2,graph(1))],[(1,graph(3)),(3,graph(4)),(4,graph(5))]]
    result=check(chains)
    assert result['start']['part3:read']<result['start']['part2:read']


def test_shared_token_held_by_other_chain_survives_admission():
    first=ExecutionGraph();first.capacities={'host_serial':1.}
    first.add('wait',1.)
    second=ExecutionGraph();second.capacities={'host_serial':1.}
    second.add('critical',10.)
    second.token_capacities['exclusive:host_serial']=1
    second.token_actions['critical']={'acquire':{'exclusive:host_serial':1},'release_finish':{'exclusive:host_serial':1}}
    chains=[[(0,first),(2,second)],[(1,second)]]
    result=check(chains)
    assert result['start']['part2:critical']==10.
    assert result['seconds']==20.


@pytest.mark.parametrize('chains',[
    [[(0,graph(1))],[(0,graph(2))]],
    [[(1,graph(1)),(0,graph(2))]],
    [[(-1,graph(1))]],
    [[(True,graph(1))]],
])
def test_rejects_ambiguous_order(chains):
    with pytest.raises(ValueError,match='order'):ExecutionGraph().solve_chains(chains)


def test_rejects_conflicting_future_capacities_and_cycles():
    a,b=graph(1),graph(2);b.capacities['cpu']=9.
    with pytest.raises(ValueError,match='capacity'):ExecutionGraph().solve_chains([[(0,a),(1,b)]])
    b=ExecutionGraph();b.add('x',1.,['y']);b.add('y',1.,['x'])
    with pytest.raises(ValueError,match='Cyclic'):ExecutionGraph().solve_chains([[(0,a),(1,b)]])


def test_empty_graph_parts_and_empty_chains():
    check([[],[(0,ExecutionGraph()),(1,ExecutionGraph())]])
    result=ExecutionGraph().solve_chains([])
    assert result['seconds']==result['scheduled_nodes']==result['peak_active_nodes']==0
