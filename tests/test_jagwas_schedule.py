"""Source-topology and hand-calculated resource controls for joint reductions."""
import copy
import pytest
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.indexed_schedule import jagwas_variant_schedule
from test_significant_schedule import tile,step


def prepare(seconds=7.):
    graph=ExecutionGraph();graph.add('preprocess',seconds)
    return graph


def shard(device='cuda:1',*,factor=3.,**kwargs):
    result=tile(device=device,backend='host',**kwargs)
    result['prepare']=prepare(factor)
    return result


def run(shards,**kwargs):
    multiple=len(shards)>1
    opts=dict(shared_prepare=prepare(),queue_depth=1 if multiple else 0,
        shared_capacities={},queue_service=dict(put=step(0.),get=step(0.)) if multiple else None,
        finalize=step(2.))
    opts.update(kwargs)
    return jagwas_variant_schedule(shards,**opts)


def test_shared_preprocessing_runs_once_and_gates_every_private_factor():
    shards=[shard(chunks=2),shard('cuda:2',chunks=2,factor=5.)]
    result=run(shards)
    zero=run(shards,shared_prepare=prepare(0.))
    assert result['seconds']==zero['seconds']+7.
    assert result['end']['shared_prepare:complete']==7.
    assert len([key for key in result['start'] if key=='shared_prepare:preprocess'])==1
    for i,factor in enumerate([3.,5.]):
        assert result['start'][f'tile:{i}:prepare:preprocess']==7.
        assert result['end'][f'tile:{i}:prepare:complete']==7.+factor
        assert result['start'][f'tile:{i}:submit_decode:0']>=7.+factor
    assert result['shared_preprocessing_passes']==1


@pytest.mark.parametrize('retained',[0,1])
def test_host_selection_is_under_single_consumer_after_dequeue_even_for_empty_chunks(retained):
    result=run([shard(chunks=3,selection=20.,retained=retained),
                shard('cuda:2',chunks=3,selection=20.,retained=retained)])
    start,end=result['start'],result['end']
    intervals=sorted((start[name],end[name]) for name in start
                     if name.startswith('tile:') and name.endswith((':select:0',':write:0')))
    assert all(a[1]<=b[0] for a,b in zip(intervals,intervals[1:]))
    assert len(intervals)==6*(1+bool(retained))
    for i in range(2):
        for j in range(3):
            base=f'tile:{i}:jagwas:{j}:0'
            assert start[base+':select:0']>=end[base+':get:done']
            assert end[base+':put:done']<=start[base+':select:0']
    assert result['retained_variants']==6*retained
    assert result['parts']==6*bool(retained)
    assert result['seconds']>=6*20+7


def test_slow_writer_blocks_queue_not_all_gpu_work_and_preserves_arrival_order():
    result=run([shard(chunks=4,kernel=10.,writer=40.),shard('cuda:2',chunks=4,kernel=1.,writer=40.)])
    start,end=result['start'],result['end']
    fast='tile:1:jagwas:0:0';slow='tile:0:jagwas:0:0'
    assert start[fast+':select:0']<start[slow+':select:0']
    assert start['tile:1:host_start:1']<end[fast+':write:done']
    events=[]
    for name in start:
        if name.endswith(':acquire'):events.append((start[name],1))
        if name.endswith(':get_start'):events.append((start[name],-1))
    count=0
    for when in sorted({time for time,_ in events}):
        count+=sum(delta for time,delta in events if time==when)
        assert 0<=count<=1
    assert count==0
    assert len([name for name in start if name.startswith('terminal:') and name.endswith(':get_start')])==2
    assert result['seconds']==end['jagwas:finalize:done']


def test_shared_storage_cpu_and_dram_work_are_conserved():
    shards=[shard(chunks=2),shard('cuda:2',chunks=2)]
    for s in shards:
        for outputs in s['outputs']:
            outputs[0]['selection']=[dict(seconds=2.,resources=dict(cpu=1.,dram=5.))]
            outputs[0]['writer']=[dict(seconds=4.,resources=dict(cpu=.5,output=10.))]
    graph=run(shards,shared_capacities=dict(cpu=1.,dram=3.,output=2.),return_graph=True)
    result=graph.solve()
    for resource,capacity in graph.capacities.items():
        demand=sum(graph.nodes[name][0]*values.get(resource,0.) for name,values in graph.demands.items())
        assert result['seconds']+1e-10>=demand/capacity
    assert result==graph._solve_shared_python()


def test_single_device_has_no_queue_or_sentinel_and_keeps_ring_overlap():
    result=run([shard(chunks=3,writer=20.)])
    assert result['queue_depth']==0
    assert not any(name.startswith('terminal:') or name.endswith(':get_start') for name in result['start'])
    assert result['start']['tile:0:host_start:1']<result['end']['tile:0:jagwas:0:0:write:done']


@pytest.mark.parametrize('bad',['shared','private','duplicate','device','multiple_outputs'])
def test_joint_graph_refuses_unsupported_or_missing_execution_stages(bad):
    shards=[shard(),shard('cuda:2')];opts={}
    if bad=='shared':opts['shared_prepare']=None
    if bad=='private':shards[0].pop('prepare')
    if bad=='duplicate':shards[1]['device']=shards[0]['device']
    if bad=='device':shards[0]['backend']='device'
    if bad=='multiple_outputs':shards[0]['outputs'][0]*=2
    with pytest.raises(ValueError):run(shards,**opts)


def test_supplied_model_inputs_are_unchanged_and_expansion_is_bounded():
    shards=[shard(chunks=2),shard('cuda:2',chunks=2)]
    before=copy.deepcopy(shards)
    run(shards)
    for left,right in zip(shards,before):
        assert left['prepare'].nodes==right['prepare'].nodes
        assert {k:v for k,v in left.items() if k!='prepare'}=={k:v for k,v in right.items() if k!='prepare'}
    with pytest.raises(ValueError,match='max_source_chunks'):
        run(shards,max_source_chunks=3)
