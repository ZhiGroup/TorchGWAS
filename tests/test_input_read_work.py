import pytest
from torchgwas.input_read_work import buffered_read_service


def price(size=100,calls=1,**kwargs):
    rates=dict(cpu_seconds_per_byte=.01,cpu_seconds_per_call=.1);rates.update(kwargs.pop('prices',{}))
    resources=dict(storage_bytes_per_second=100,dram_bytes_per_second=1000,cpu_fraction=1.);resources.update(kwargs)
    return buffered_read_service(size,calls,rates,**resources)


def test_input_read_conserves_each_resource_service_without_adding_to_decode():
    result=price()
    assert result['cpu_seconds']==pytest.approx(1.1)
    assert result['seconds']==pytest.approx(1.1)
    assert result['seconds']*result['resources']['cpu']==pytest.approx(1.1)
    assert result['seconds']*result['resources']['input']==pytest.approx(100)
    assert result['seconds']*result['resources']['dram']==pytest.approx(200)


def test_storage_cpu_and_memory_limits_are_independent():
    assert price(storage_bytes_per_second=10)['seconds']==pytest.approx(10)
    assert price(cpu_fraction=.5)['seconds']==pytest.approx(2.2)
    assert price(dram_bytes_per_second=10)['seconds']==pytest.approx(20)


def test_empty_read_can_still_have_syscall_cpu():
    result=price(size=0)
    assert result['seconds']==pytest.approx(.1)
    assert result['resources']['input']==0
    assert price(size=0,calls=0)['seconds']==0


@pytest.mark.parametrize('kwargs',[{'calls':0},{'calls':True},{'cpu_fraction':0},{'storage_bytes_per_second':0},{'prices':{'cpu_seconds_per_byte':float('nan')}}])
def test_invalid_service_is_refused(kwargs):
    with pytest.raises(ValueError):price(**kwargs)


def direct_fixture():
    amount=64*1024*1024;block=amount//4
    rows=[]
    for repeat,seconds in enumerate([.004,.005,.006]):
        calls=[dict(part=part,worker=part,bytes=block,offset=part*block,start=0.,end=seconds,wall_seconds=seconds) for part in range(4)]
        rows.append(dict(mode='direct',cache='bypass',repeat=repeat,workers=4,bytes=amount,
                         seconds=seconds,calls=calls,resident_at_launch=0,resident_after=0,
                         process_io_delta={'read_bytes':amount}))
    return dict(affinity=[12,13],input=dict(block_bytes=block,blocks=4,bytes=amount,filesystem='private ext4'),rows=rows)


def test_direct_capacity_uses_aggregate_io_without_replacing_buffered_cpu_cost():
    from torchgwas.input_read_work import direct_read_capacity
    capacity=direct_read_capacity(direct_fixture(),workers=4,cpu_affinity=[12,13],filesystem='private ext4')
    assert capacity['bytes_per_second']==pytest.approx((64*1024*1024)/.005)
    before=price(storage_bytes_per_second=100)
    after=price(storage_bytes_per_second=capacity['bytes_per_second'])
    assert before['cpu_seconds']==after['cpu_seconds']
    assert after['seconds']==pytest.approx(after['cpu_seconds'])
    assert capacity['unpriced_terms']


@pytest.mark.parametrize('bad',['cache','resident','read_bytes','coverage','interval','boundary','repeat','filesystem','affinity','worker'])
def test_direct_capacity_rejects_unverified_io_or_context(bad):
    from torchgwas.input_read_work import direct_read_capacity
    probe=direct_fixture();row=probe['rows'][0]
    kwargs=dict(workers=4,cpu_affinity=[12,13],filesystem='private ext4')
    if bad=='cache':row['cache']='warm'
    elif bad=='resident':row['resident_after']=1
    elif bad=='read_bytes':row['process_io_delta']['read_bytes']=0
    elif bad=='coverage':row['calls'][0]['part']=1
    elif bad=='interval':row['calls'][0]['wall_seconds']=.1
    elif bad=='boundary':row['seconds']=.1
    elif bad=='repeat':probe['rows'][1]['repeat']=0
    elif bad=='filesystem':kwargs['filesystem']='unrelated nfs'
    elif bad=='affinity':kwargs['cpu_affinity']=[0,1]
    else:kwargs['workers']=2
    with pytest.raises(ValueError):direct_read_capacity(probe,**kwargs)