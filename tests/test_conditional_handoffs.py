import pytest
from torchgwas.execution_graph import ExecutionGraph,torch_scan_schedule,torch_multigpu_schedule


def test_resume_delay_depends_on_endpoint_order_not_configured_duration_alone():
    g=ExecutionGraph();g.add('attempt',.1);g.add('ready',.5)
    wake=g.add_wait_delay('wake',.05,'attempt','ready');g.add('work',.1,[wake])
    result=g.solve()
    assert result['seconds']==pytest.approx(.65)
    assert result['conditional_delays']['wake']['blocked']
    assert result['conditional_delays']['wake']['dependency_wait_seconds']==pytest.approx(.4)
    g.nodes['attempt']=(.6,())
    result=g.solve()
    assert result['seconds']==pytest.approx(.7)
    assert not result['conditional_delays']['wake']['blocked']
    assert result['conditional_delays']['wake']['extra_elapsed_service_seconds']==0


def test_resource_contention_can_change_whether_the_caller_blocks():
    g=ExecutionGraph();g.capacities={'cpu':1.}
    g.add('attempt',1.,resources={'cpu':1.});g.add('competing',1.,resources={'cpu':1.})
    g.add('ready',1.5);g.add_wait_delay('wake',.5,'attempt','ready')
    result=g.solve()
    assert result['seconds']==pytest.approx(2.)
    assert not result['conditional_delays']['wake']['blocked']


def test_equal_endpoint_times_do_not_create_a_wakeup():
    g=ExecutionGraph();g.add('attempt',1.);g.add('ready',1.)
    g.add_wait_delay('wake',5.,'attempt','ready')
    assert g.solve()['seconds']==1.


def test_conditional_delays_validate_endpoints_and_dependencies():
    g=ExecutionGraph();g.add('attempt');g.add('wake',1.)
    g.conditional_delays['wake']=dict(attempt='attempt',ready='missing')
    with pytest.raises(ValueError,match='endpoint'):g.solve()
    g.add('missing')
    with pytest.raises(ValueError,match='dependencies'):g.solve()


def blocks():
    return [dict(decode_seconds=1.,h2d_seconds=0.,h2d_bytes=0.,d2h_seconds=0.,d2h_bytes=0.,
                 operations=[],host_submit_seconds=0.,finish_seconds=0.,consumer_seconds=0.,
                 handoff_wakeup_seconds=dict(queue=.01,future=.02)) for _ in range(3)]


def test_full_ring_handoffs_follow_the_source_publish_order():
    result=torch_scan_schedule(blocks(),depth=2,decode_workers=2)
    assert result['seconds']==pytest.approx(2.07)
    # The second decode is complete before the producer can submit the next
    # job and emit that completed future. No future wake belongs on that path.
    assert not result['conditional_delays']['wake_publish:1']['blocked']
    assert result['conditional_delays']['wake_free:2']['blocked']
    assert result['start']['host_start:2']>=result['end']['fetch:2']
    assert result['start']['resolve:0']>=result['end']['fetch:2']


def test_multigpu_prefixes_conditional_endpoints_without_cross_shard_dependencies():
    shards=[dict(device=str(i),blocks=blocks(),depth=2,decode_workers=2) for i in range(2)]
    result=torch_multigpu_schedule(shards,{})
    assert result['seconds']==pytest.approx(2.07)
    for i in range(2):
        assert result['conditional_delays'][f'shard{i}:wake_free:2']['blocked']
        assert not result['conditional_delays'][f'shard{i}:wake_publish:1']['blocked']
    single=torch_multigpu_schedule(shards[:1],{})
    assert single['seconds']==pytest.approx(2.07)
    assert 'wake_free:2' in single['conditional_delays']


def test_partial_or_invalid_block_handoff_inputs_fail():
    data=blocks();data[-1].pop('handoff_wakeup_seconds')
    with pytest.raises(ValueError,match='Every block'):torch_scan_schedule(data,depth=2,decode_workers=2)
    data=blocks();data[0]['handoff_wakeup_seconds']['queue']=-1
    with pytest.raises(ValueError,match='Invalid handoff'):torch_scan_schedule(data,depth=2,decode_workers=2)
