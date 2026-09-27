import copy
import pytest
from torchgwas.execution_graph import torch_scan_schedule,torch_multigpu_schedule


def block(read=0.,compute=1.,transfer=0.):
    return dict(decode_seconds=0.,decode_read_seconds=read,read_resources={'disk':1} if read else {},
        h2d_seconds=transfer,h2d_bytes=transfer,d2h_seconds=0.,d2h_bytes=0.,
        operations=[dict(host_submit_finish=0.,kernel_service_seconds=compute)],
        host_submit_seconds=0.,finish_seconds=0.,consumer_seconds=0.)


def shard(name,blocks):
    return dict(device=name,blocks=blocks,depth=2,decode_workers=1)


def test_independent_gpus_overlap_but_common_disk_does_not_double():
    b=[block(compute=1.) for _ in range(4)]
    one=torch_scan_schedule(b,depth=2,decode_workers=1)['seconds']
    two=torch_multigpu_schedule([shard('0',b),shard('1',b)],{})['seconds']
    assert one==two==4.
    read=[block(read=1.,compute=0.) for _ in range(4)]
    shared=torch_multigpu_schedule([shard('0',read),shard('1',read)],{'disk':1})
    assert shared['seconds']==8.


def test_shared_pcie_link_is_not_two_independent_links():
    b=[block(compute=0.,transfer=1.) for _ in range(3)]
    shards=[shard('0',b),shard('1',b)]
    independent=torch_multigpu_schedule(shards,{})['seconds']
    common=torch_multigpu_schedule(shards,{},[dict(devices=['0','1'],h2d_bytes_per_second=1.,d2h_bytes_per_second=1.)])['seconds']
    assert independent==3. and common==6.


def test_ordered_delivery_can_stall_later_shard():
    shards=[shard('0',[block(compute=1.) for _ in range(12)]),
            shard('1',[block(compute=1.) for _ in range(12)])]
    concurrent=torch_multigpu_schedule(shards,{},ordered=False,result_queue_depth=1)
    ordered=torch_multigpu_schedule(shards,{},ordered=True,result_queue_depth=1)
    assert concurrent['seconds']==12.
    assert ordered['seconds']>concurrent['seconds']
    assert ordered['end']['shard1:deliver:0']>=ordered['end']['shard0:deliver:11']


def test_composition_keeps_original_blocks_and_counts():
    b=[block() for _ in range(3)];snapshot=copy.deepcopy(b)
    composed=torch_multigpu_schedule([shard('0',b)],{})
    assert b==snapshot
    assert composed['seconds']==torch_scan_schedule(b,depth=2,decode_workers=1)['seconds']
    with pytest.raises(ValueError,match='One active'):
        torch_multigpu_schedule([shard('0',b),shard('0',b)],{})


@pytest.mark.parametrize('ordered',[False,True])
@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
def test_borrowed_acknowledgement_protects_ring_and_next_result(ordered,policy):
    blocks=[dict(block(compute=0.),host_resources={'cpu':1.},handoff_wakeup_seconds={'queue':.001,'future':.001}) for _ in range(5)]
    shards=[shard('0',blocks),shard('1',blocks)]
    original=copy.deepcopy(shards)
    ack=dict(create_cpu_seconds=.25,publish_cpu_seconds=2.,receive_cpu_seconds=.5,wakeup_seconds=.1)
    result=torch_multigpu_schedule(shards,{'cpu':8.},ordered=ordered,borrow_results=True,
        acknowledgement_service=ack,queue_service=dict(put_cpu_seconds=0.,get_cpu_seconds=1.,cpu_fraction=1.),
        host_serial_fraction=1.,host_serial_policy=policy)
    for device in range(2):
        for i in range(5):
            resume=result['end'][f'shard{device}:ack_resume:{i}']
            assert resume>=result['end'][f'shard{device}:deliver:{i}']+.5
            if i+1<5:assert result['start'][f'shard{device}:resolve:{i+1}']>=resume
            if i+1<5:assert result['start'][f'shard{device}:attempt_resolve:{i+1}']>=resume
            if i+2<5:assert result['start'][f'shard{device}:host_start:{i+2}']>=resume
    assert result['seconds']>=30. and shards==original


def test_single_device_borrowing_has_no_cross_thread_acknowledgement():
    blocks=[block() for _ in range(4)]
    result=torch_multigpu_schedule([shard('0',blocks)],{},borrow_results=True)
    assert result['seconds']==torch_scan_schedule(blocks,depth=2,decode_workers=1)['seconds']
    assert not any('ack_' in name for name in result['end'])
    with pytest.raises(ValueError,match='acknowledgement'):
        torch_multigpu_schedule([shard('0',blocks),shard('1',blocks)],{},borrow_results=True)


def test_shared_host_serial_resource_limits_cpu_dispatch_without_serializing_gpu():
    b=block(compute=0.)
    b.update(host_resources={'cpu':1.},operations=[dict(host_submit_finish=1.,kernel_service_seconds=0.)],host_submit_seconds=1.)
    shards=[shard('0',[b]),shard('1',[b])]
    original=copy.deepcopy(shards)
    independent=torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=0.)
    serialized=torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=1.)
    assert independent['seconds']==1.
    assert serialized['seconds']==2.
    assert shards==original
    b['operations'][0]['kernel_service_seconds']=3.
    # Shared Python dispatch still permits each device to compute concurrently.
    concurrent=torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=1.)
    assert concurrent['seconds']==5.
    assert concurrent['resource_policy']['host_serial_fraction']==1.


def test_host_serial_scenario_validates_fraction_and_cpu_demand():
    for value in [-.1,1.1,float('nan'),True]:
        with pytest.raises(ValueError,match='fraction'):
            torch_multigpu_schedule([shard('0',[block()])],{},host_serial_fraction=value)
    with pytest.raises(ValueError,match='host CPU demand'):
        torch_multigpu_schedule([shard('0',[block()])],{},host_serial_fraction=.5)


def test_shared_result_consumer_serializes_delivery_across_devices():
    b=[block(compute=0.) for _ in range(4)]
    result=torch_multigpu_schedule([shard('0',b),shard('1',b)],{'cpu':8.},result_queue_depth=1,
        queue_service=dict(put_cpu_seconds=0.,get_cpu_seconds=1.,cpu_fraction=1.))
    assert result['seconds']==8.
    deliveries=sorted((result['start'][name],result['end'][name]) for name in result['start'] if ':deliver:' in name)
    assert all(left[1]<=right[0] for left,right in zip(deliveries,deliveries[1:]))
    assert result['resource_policy']['result_queue_capacity']==2


def test_queue_tokens_block_producers_until_consumer_takes_item():
    from torchgwas.execution_graph import ExecutionGraph
    g=ExecutionGraph();g.token_capacities={'queue':2,'consumer':1}
    for i in range(5):
        put=g.add(f'put{i}')
        get=g.add(f'get{i}',1.,[put]+([f'get{i-1}'] if i else []))
        g.token_actions[put]={'acquire':{'queue':1}}
        g.token_actions[get]={'acquire':{'consumer':1},'release_start':{'queue':1},'release_finish':{'consumer':1}}
    result=g.solve()
    assert result['seconds']==5.
    assert result['start']['put4']==2.
    assert result['start']['put3']==1.


def test_impossible_token_dependency_fails_instead_of_spinning():
    from torchgwas.execution_graph import ExecutionGraph
    g=ExecutionGraph();g.token_capacities={'slot':1}
    g.add('holder');g.add('waiter',1.,['holder'])
    g.token_actions={'holder':{'acquire':{'slot':1}},'waiter':{'acquire':{'slot':1}}}
    with pytest.raises(ValueError,match='deadlock'):g.solve()


def test_fifo_delivery_follows_enqueue_completion_not_get_readiness():
    from torchgwas.execution_graph import ExecutionGraph
    g=ExecutionGraph();g.token_capacities={'queue':3,'consumer':1}
    # A second item from producer A arrives before B, but its get depends
    # on A's first get. B must not jump ahead just because its get is ready.
    for name,duration,after in [('a0',0.,[]),('a1',.1,[]),('b0',.2,[])]:
        put=g.add('put:'+name,duration,after)
        get=g.add('get:'+name,1.,[put]+(['get:a0'] if name=='a1' else []))
        g.token_actions[put]={'acquire':{'queue':1}}
        g.token_actions[get]={'acquire':{'consumer':1},'release_start':{'queue':1},'release_finish':{'consumer':1}}
        g.fifo_enqueues[put]=('queue',get);g.fifo_dequeues[get]='queue'
    result=g.solve()
    assert result['start']['get:a0']==0.
    assert result['start']['get:a1']==1.
    assert result['start']['get:b0']==2.


def test_single_device_bypasses_multigpu_result_queue_costs():
    blocks=[block(compute=1.) for _ in range(4)]
    direct=torch_scan_schedule(blocks,depth=2,decode_workers=1)
    result=torch_multigpu_schedule([shard('0',blocks)],{'cpu':1.},result_queue_depth=1,
        queue_service=dict(put_cpu_seconds=100.,get_cpu_seconds=100.,cpu_fraction=1.))
    assert result['seconds']==direct['seconds']
    assert result['resource_policy']['result_queue_capacity']==0
    assert result['resource_policy']['queue_service'] is None
    assert not any('deliver:' in name for name in result['start'])


def test_gil_free_finish_copy_keeps_cpu_and_memory_limits():
    b=block(compute=0.)
    b.update(host_resources={'cpu':1.},finish_seconds=1.1,
             finish_operations=[dict(seconds=.1),dict(seconds=1.,host_serial_fraction=0.,resources={'dram':1.})])
    shards=[shard('0',[b]),shard('1',[b])]
    parallel=torch_multigpu_schedule(shards,{'cpu':2.,'dram':2.},host_serial_fraction=1.)
    assert parallel['seconds']==pytest.approx(1.2)
    assert parallel['end']['shard0:finish:0']>=1.1
    cpu_limited=torch_multigpu_schedule(shards,{'cpu':1.,'dram':2.},host_serial_fraction=1.)
    assert cpu_limited['seconds']==pytest.approx(2.2)
    memory_limited=torch_multigpu_schedule(shards,{'cpu':2.,'dram':1.},host_serial_fraction=1.)
    assert memory_limited['seconds']>=2.1
    assert shards[0]['blocks'][0]['finish_operations'][1]['host_serial_fraction']==0.


def test_finish_phase_service_must_be_conserved():
    b=block(compute=0.)
    b.update(finish_seconds=1.,finish_operations=[dict(seconds=2.)])
    with pytest.raises(ValueError,match='conserve'):
        torch_scan_schedule([b],depth=2,decode_workers=1)
