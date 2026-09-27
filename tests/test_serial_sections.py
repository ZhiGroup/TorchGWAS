import copy
import pytest
from torchgwas.execution_graph import ExecutionGraph,torch_multigpu_schedule


def work(g,resource):
    return sum(seconds*g.demands.get(name,{}).get(resource,0.) for name,(seconds,_) in g.nodes.items())


def test_long_held_call_blocks_short_host_call_while_gpu_keeps_running():
    g=ExecutionGraph();g.capacities={'cpu':3.,'host_serial':1.}
    g.add('array_free',2.,resources={'cpu':1.,'host_serial':1.})
    g.add('gpu',3.)
    g.add('api_ready',.5)
    g.add('api',.1,['api_ready'],resources={'cpu':1.,'host_serial':1.})
    fluid=g.solve()
    serial=g.with_serial_sections('held-first').solve()
    assert fluid['end']['api']==pytest.approx(.7)
    assert serial['start']['api']==2. and serial['end']['api']==pytest.approx(2.1)
    assert serial['end']['gpu']==3.
    assert g.token_capacities=={}


@pytest.mark.parametrize('position',['held-first','held-last'])
def test_mixed_calls_preserve_cpu_dram_and_all_completion_dependencies(position):
    g=ExecutionGraph();g.capacities={'cpu':2.,'dram':4.,'host_serial':1.}
    g.add('call',10.,resources={'cpu':.5,'host_serial':.1,'dram':2.})
    g.add('done',1.,['call'])
    old=copy.deepcopy(g.__dict__);s=g.with_serial_sections(position)
    assert work(s,'cpu')==work(g,'cpu')==5.
    assert work(s,'dram')==work(g,'dram')==20.
    assert work(s,'host_serial')==0.
    held=[name for name,actions in s.token_actions.items() if 'exclusive:host_serial' in actions.get('acquire',{})]
    assert len(held)==1 and s.nodes[held[0]][0]==2.
    result=s.solve()
    assert result['end']['call']==10. and result['start']['done']==10.
    assert g.__dict__==old


def test_cpu_contention_extends_critical_section_ownership():
    g=ExecutionGraph();g.capacities={'cpu':1.,'host_serial':1.}
    g.add('holder',1.,resources={'cpu':1.,'host_serial':1.})
    g.add('other_cpu',1.,resources={'cpu':1.})
    g.add('delay',.5)
    g.add('waiter',.1,['delay'],resources={'cpu':1.,'host_serial':1.})
    result=g.with_serial_sections('held-first').solve()
    assert result['end']['holder']==2.
    assert result['start']['waiter']==2.


@pytest.mark.parametrize('position',['held-first','held-last'])
def test_existing_queue_tokens_and_fifo_survive_split(position):
    g=ExecutionGraph();g.capacities={'cpu':2.,'host_serial':1.}
    g.token_capacities={'queue':1,'consumer':1}
    for i in range(2):
        put=g.add(f'put{i}')
        get=g.add(f'get{i}',2.,[put],resources={'cpu':1.,'host_serial':.5})
        g.token_actions[put]={'acquire':{'queue':1}}
        g.token_actions[get]={'acquire':{'consumer':1},'release_start':{'queue':1},'release_finish':{'consumer':1}}
        g.fifo_enqueues[put]=('queue',get);g.fifo_dequeues[get]='queue'
    s=g.with_serial_sections(position);result=s.solve()
    assert result['seconds']==4.
    assert result['start']['get1:serial_section_start']==2.
    assert result['end']['get0']==2. and result['end']['get1']==4.


def test_unknown_section_orders_and_invalid_ledgers_rejected():
    g=ExecutionGraph();g.capacities={'cpu':1.,'host_serial':1.}
    g.add('call',1.,resources={'cpu':1.,'host_serial':2.})
    with pytest.raises(ValueError,match='within'):g.with_serial_sections('held-first')
    with pytest.raises(ValueError,match='order'):g.with_serial_sections('invented')
    g.demands['call']['host_serial']=1.
    g.resource_waits['wait']={'after':'call','until':'call','resources':{'host_serial':1.}}
    with pytest.raises(ValueError,match='waits'):g.with_serial_sections('held-first')


@pytest.mark.parametrize('position',['held-first','held-last'])
def test_two_gpu_api_partitions_keep_one_shared_critical_section(position):
    b=dict(decode_seconds=0.,h2d_seconds=0.,d2h_seconds=0.,finish_seconds=0.,consumer_seconds=0.,
        host_resources={'cpu':1.},host_submit_seconds=10.,host_submit_serial_cpu_seconds=2.,
        operations=[dict(host_submit_finish=10.,host_serial_cpu_finish=2.,kernel_service_seconds=0.)])
    shards=[dict(device=str(i),blocks=[b],depth=2,decode_workers=1) for i in range(2)]
    result=torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=1.,host_serial_policy=position)
    assert result['seconds']==pytest.approx(12.)
    assert result['resource_policy']['host_serial_policy']==position


def test_owned_discard_is_charged_once_on_shared_consumer():
    b=dict(decode_seconds=0.,h2d_seconds=0.,d2h_seconds=0.,finish_seconds=1.,consumer_seconds=0.,discard_seconds=2.,
        host_resources={'cpu':1.},host_submit_seconds=0.,operations=[],
        finish_operations=[dict(seconds=1.,host_serial_fraction=0.)])
    shards=[dict(device=str(i),blocks=[b],depth=2,decode_workers=1) for i in range(2)]
    one=torch_multigpu_schedule(shards[:1],{'cpu':2.},host_serial_fraction=1.,host_serial_policy='held-first')
    two=torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=1.,host_serial_policy='held-first')
    assert one['seconds']==3. and two['seconds']==5.
    assert all(two['end'][f'shard{i}:finish:0']==1. for i in range(2))
    assert sorted(two['end'][f'shard{i}:deliver:0'] for i in range(2))==[3.,5.]
