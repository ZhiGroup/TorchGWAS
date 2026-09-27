import copy

import pytest

from torchgwas.output_write_work import direct_write_capacity, with_direct_write_capacity


def fixture():
    block, count = 16 << 20, 64
    rows = []
    for workers in (1, 4):
        for repeat in range(5):
            duration = .001*(repeat+1)
            calls = [dict(worker=worker, part=part, bytes=block, offset=part*block,
                          start=part*duration, end=(part+1)*duration,
                          wall_seconds=duration, cpu_seconds=duration/2)
                     for worker in range(workers) for part in range(count//workers)]
            rows.append(dict(mode='direct', repeat=repeat, workers=workers, bytes=block*count,
                seconds=count//workers*duration, fsync_seconds=9., calls=calls,
                sum_cpu_seconds=sum(call['cpu_seconds'] for call in calls),
                resident_at_launch=[0]*workers, resident_after=[0]*workers, verified_payload=True,
                write_io_delta={'write_bytes': block*count}, durable_io_delta={'write_bytes': block*count}))
    return dict(affinity=[12, 13], input=dict(filesystem='private ext4', block_bytes=block,
        blocks=count, bytes=block*count, preallocated=True), rows=rows)


def read(probe):
    return direct_write_capacity(probe, workers=4, cpu_affinity=[12, 13], filesystem='private ext4')


def test_output_service_separates_stream_shared_cpu_and_commit_costs():
    probe = fixture()
    probe['rows'].append(dict(mode='buffered', workers=4, seconds=.000001))
    price = read(probe)
    assert price['bytes_per_second'] == pytest.approx((1 << 30)/(.003*16))
    assert price['unpriced_terms']
    original = dict(write_bytes_per_second=1., writeback_service=dict(
        pagecache_seconds_per_byte=.4, storage_seconds_per_byte=1., wait_seconds=.2))
    snapshot = copy.deepcopy(original)
    updated = with_direct_write_capacity(original, probe, cpu_affinity=[12, 13], filesystem='private ext4')
    assert original == snapshot
    assert updated['writeback_service']['pagecache_seconds_per_byte'] == .4
    assert updated['writeback_service']['wait_seconds'] == .2
    assert updated['writeback_service']['storage_seconds_per_byte'] == pytest.approx((.003*64)/(1 << 30))
    assert updated['write_bytes_per_second'] == price['bytes_per_second']


@pytest.mark.parametrize('bad', ['buffered', 'pages', 'write_bytes', 'durable_bytes', 'payload', 'repeat',
    'coverage', 'offset', 'interval', 'boundary', 'sequential', 'nan', 'cpu', 'filesystem', 'affinity', 'extent'])
def test_unverified_capacity_is_refused(bad):
    probe = fixture()
    row = next(row for row in probe['rows'] if row['workers'] == 4)
    if bad == 'buffered': row['mode'] = 'buffered'
    elif bad == 'pages': row['resident_after'][0] = 1
    elif bad == 'write_bytes': row['write_io_delta']['write_bytes'] = 0
    elif bad == 'durable_bytes': row['durable_io_delta']['write_bytes'] = 0
    elif bad == 'payload': row['verified_payload'] = False
    elif bad == 'repeat': row['repeat'] = 1
    elif bad == 'coverage': row['calls'][0]['part'] = 1
    elif bad == 'offset': row['calls'][0]['offset'] = 1
    elif bad == 'interval': row['calls'][0]['wall_seconds'] = .4
    elif bad == 'boundary': row['seconds'] += .001
    elif bad == 'sequential':
        row['calls'][1].update(start=0., end=.001)
    elif bad == 'nan': row['seconds'] = float('nan')
    elif bad == 'cpu': row['sum_cpu_seconds'] += .1
    elif bad == 'filesystem': probe['input']['filesystem'] = 'unrelated nfs'
    elif bad == 'affinity': probe['affinity'] = [0, 1]
    else: probe['input']['preallocated'] = False
    with pytest.raises(ValueError): read(probe)


def copy_fixture():
    rows=[]
    for repeat in range(5):
        for size,loops,cost in [(0,10000,1e-6),(64<<20,16,1e-6+(64<<20)*1e-10)]:
            cpu=loops*cost
            rows.append(dict(mechanism='numpy_copyto',workers=1,payload_bytes=size,loops=loops,
                repeat=repeat,cpu_seconds=cpu,seconds=cpu*2,threads=[dict(worker=0,cpu_seconds=cpu)]))
    return dict(affinity=[12,13],python_version='test-python',numpy_version='test-numpy',rows=rows)


@pytest.mark.parametrize('mechanism',['numpy_copyto','memoryview_slice'])
def test_staging_service_separates_view_dispatch_and_bulk_without_shape_fit(mechanism):
    from torchgwas.output_write_work import staging_copy_prices,writer_copy_cost
    probe=copy_fixture()
    for row in probe['rows']:row['mechanism']=mechanism
    probe['rows'].append(dict(mechanism='bytearray_slice',workers=1,payload_bytes=64<<20))
    price=staging_copy_prices(probe,cpu_affinity=[12,13],python_version='test-python',numpy_version='test-numpy',mechanism=mechanism)
    assert price['cpu_seconds_per_call']==pytest.approx(1e-6)
    assert price['cpu_seconds_per_byte']==pytest.approx(1e-10)
    assert writer_copy_cost(dict(writer_copy_service=price))==dict(cpu_seconds_per_call=price['cpu_seconds_per_call'],cpu_seconds_per_byte=price['cpu_seconds_per_byte'])
    probe['rows'][0]['repeat']=1
    with pytest.raises(ValueError,match='Incomplete'):
        staging_copy_prices(probe,cpu_affinity=[12,13],python_version='test-python',numpy_version='test-numpy',mechanism=mechanism)


@pytest.mark.parametrize('bad',['affinity','python','numpy','loops','worker','cpu','missing_bulk','nonpositive_bulk'])
def test_staging_controls_require_complete_matching_context(bad):
    from torchgwas.output_write_work import staging_copy_prices
    probe=copy_fixture()
    if bad=='affinity':probe['affinity']=[0,1]
    elif bad=='python':probe['python_version']='other'
    elif bad=='numpy':probe['numpy_version']='other'
    elif bad=='loops':probe['rows'][0]['loops']=1
    elif bad=='worker':probe['rows'][0]['threads'][0]['worker']=1
    elif bad=='cpu':probe['rows'][0]['threads'][0]['cpu_seconds']=0.
    elif bad=='missing_bulk':probe['rows']=[row for row in probe['rows'] if row['payload_bytes']==0]
    else:
        for row in probe['rows']:
            if row['payload_bytes']:
                row['cpu_seconds']=row['loops']*1e-7
                row['threads'][0]['cpu_seconds']=row['cpu_seconds']
    with pytest.raises(ValueError):
        staging_copy_prices(probe,cpu_affinity=[12,13],python_version='test-python',numpy_version='test-numpy')


@pytest.mark.parametrize('key,value',[('cpu_seconds_per_call',-1),('cpu_seconds_per_byte',float('nan'))])
def test_invalid_staging_service_is_refused(key,value):
    from torchgwas.output_write_work import writer_copy_cost
    price=dict(cpu_seconds_per_call=0.,cpu_seconds_per_byte=1.)
    price[key]=value
    with pytest.raises(ValueError):writer_copy_cost(dict(writer_copy_service=price))
