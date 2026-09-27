import pytest
from torchgwas.execution_graph import ExecutionGraph,torch_multigpu_schedule


def graph(cpu_capacity):
    g=ExecutionGraph();g.capacities={'cpu':cpu_capacity}
    g.add('begin');g.add('device',1.)
    g.add('cpu',1.,resources={'cpu':1.})
    g.resource_waits['spin']=dict(after='begin',until='device',resources={'cpu':1.})
    return g


def test_wait_cpu_stops_when_device_completes_not_after_fixed_cpu_work():
    result=graph(1.).solve()
    assert result['end']['device']==1.
    assert result['seconds']==pytest.approx(1.5)
    assert result['wait_resource_seconds']['cpu']==pytest.approx(.5)
    uncongested=graph(2.).solve()
    assert uncongested['seconds']==1.
    assert uncongested['wait_resource_seconds']['cpu']==1.


def test_already_ready_event_does_not_consume_wait_cpu():
    g=graph(1.)
    g.nodes['begin']=(0.,('device',))
    result=g.solve()
    assert result['seconds']==1.
    assert result['wait_resource_seconds'].get('cpu',0.)==0.


@pytest.mark.parametrize('zero_duration', [False, True])
def test_wait_with_same_endpoint_never_becomes_active(zero_duration):
    g=ExecutionGraph();g.capacities={'cpu':1.}
    g.add('endpoint',0. if zero_duration else .5)
    g.add('work',1.,resources={'cpu':1.})
    g.resource_waits['spin']=dict(after='endpoint',until='endpoint',resources={'cpu':1.})
    result=g.solve()
    assert result['seconds']==1.
    assert result['wait_resource_seconds'].get('cpu',0.)==0.


@pytest.mark.parametrize('reverse_order', [False, True])
def test_simultaneous_wait_endpoints_do_not_leave_stale_demand(reverse_order):
    g=ExecutionGraph();g.capacities={'cpu':1.}
    for name in (('until','after') if reverse_order else ('after','until')):
        g.add(name,.5)
    g.add('work',1.,resources={'cpu':1.})
    g.resource_waits['spin']=dict(after='after',until='until',resources={'cpu':1.})
    result=g.solve()
    assert result['seconds']==1.
    assert result['wait_resource_seconds'].get('cpu',0.)==0.


def test_wait_requires_known_endpoints_and_capacities():
    g=graph(1.);g.resource_waits['spin']['until']='missing'
    with pytest.raises(ValueError,match='endpoint'):g.solve()
    g=graph(1.);g.resource_waits['spin']['resources']={'missing':1.}
    with pytest.raises(ValueError,match='capacity'):g.solve()


def test_multigpu_composes_waits_as_shared_cpu_not_python_serial_work():
    b=dict(decode_seconds=0.,h2d_seconds=0.,h2d_bytes=0.,d2h_seconds=0.,d2h_bytes=0.,
           operations=[dict(host_submit_finish=0.,kernel_service_seconds=1.)],
           host_submit_seconds=0.,finish_seconds=0.,consumer_seconds=0.,
           host_resources={'cpu':1.},event_wait_resources={'cpu':1.})
    shards=[dict(device=str(i),blocks=[b],depth=2,decode_workers=1) for i in range(2)]
    result=torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=1.)
    assert result['seconds']==1.
    assert result['wait_resource_seconds']['cpu']==2.
    for i in range(2):
        assert result['end'][f'shard{i}:begin_finish_wait:0']==0.
