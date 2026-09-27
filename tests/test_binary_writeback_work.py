import pytest
from torchgwas.binary_writeback_work import binary_writeback_work
from torchgwas.binary_output_work import binary_output_work


@pytest.mark.parametrize('borrow',[False,True])
@pytest.mark.parametrize('store_beta',[False,True])
def test_variant_df_stream_matches_actual_writer_and_graph(tmp_path,borrow,store_beta):
    import numpy as np
    from torchgwas.sumstats import BinarySumstatsWriter
    from torchgwas.binary_schedule import BinaryWriterSchedule
    from torchgwas.execution_graph import ExecutionGraph
    m,k,b=11,7,4
    work=binary_output_work(m,k,b,block_bytes=48,queue_depth=1,borrow_chunks=borrow,
                            store_beta=store_beta,store_variant_df=True,fsync=False,writeback_bytes=0)
    actual=BinarySumstatsWriter(tmp_path,m,list(range(k)),40,38,block_bytes=48,
        queue_depth=1,borrow_chunks=borrow,store_beta=store_beta,store_variant_df=True,fsync=False,writeback_bytes=0)
    model=BinaryWriterSchedule(work,copy_seconds_per_byte=1e-6,zero_seconds_per_byte=1e-6,
        write_seconds_per_byte=1e-6,fsync_seconds_per_array=0,writeback_bytes=0)
    g=ExecutionGraph();last=g.add('start')
    for i,start in enumerate(range(0,m,b)):
        end=min(start+b,m)
        actual.write_chunk(start,end,np.ones((end-start,k)),np.ones((end-start,k)),np.full((end-start,1),38))
        last=model.append(g,i,[last])
    summary=actual.close();model.close(g,[last]);g.solve()
    assert model.payload==work['binary_payload_bytes']==summary['payload_bytes']
    assert model.write_calls==work['binary_write_calls_minimum']==summary['blocks']
    assert work['df_payload_bytes']==summary['df_payload_bytes']==4*m
    assert model.copied==work['staging_copy_bytes']
    assert model.copy_calls==work['staging_copy_calls']
    assert work['staging_view_creations']==2*model.copy_calls


def test_staging_call_overhead_counts_partial_block_splits():
    from torchgwas.binary_schedule import BinaryWriterSchedule
    from torchgwas.execution_graph import ExecutionGraph
    work=binary_output_work(5,3,2,block_bytes=17,store_beta=False,
                            borrow_chunks=False,fsync=False,writeback_bytes=0)
    # Chunks of24,24,12 bytes cross17-byte blocks in2,2,2 pieces.
    assert work['staging_copy_calls']==6
    writer=BinaryWriterSchedule(work,copy_seconds_per_byte=0.,copy_seconds_per_call=.5,
        zero_seconds_per_byte=0.,write_seconds_per_byte=0.,fsync_seconds_per_array=0.,writeback_bytes=0)
    graph=ExecutionGraph();last=graph.add('start')
    for index in range(3):last=writer.append(graph,index,[last])
    end=writer.close(graph,[last])
    assert graph.solve()['end'][end]==pytest.approx(3.)


def test_write_crossing_multiple_intervals_waits_previous_and_retains_tail():
    work=binary_writeback_work([{'bytes':35}],10)
    assert [(a['kind'],a['offset']) for a in work['events'][0]['actions']]==[
        ('submit',0),('submit',10),('wait',0),('drop_cache',0),
        ('submit',20),('wait',10),('drop_cache',10)]
    assert work['submitted_bytes']==30
    assert work['waited_bytes']==20
    assert work['unsubmitted_tail_bytes']==5
    assert work['not_explicitly_waited_before_fsync_bytes']==15


@pytest.mark.parametrize('sizes',[[9],[10],[19],[20],[7,7,7],[35],[10,10,5]])
def test_ledger_matches_actual_writer_advance_method(sizes,monkeypatch):
    # Execute the production method with mocked syscalls, without storage timing.
    from torchgwas import sumstats
    calls=[]
    monkeypatch.setattr(sumstats,'_SYNC_FILE_RANGE',lambda fd,start,length,flags:
                        calls.append(('submit' if flags==sumstats._SYNC_FILE_RANGE_WRITE else 'wait',start,length)))
    monkeypatch.setattr(sumstats.os,'posix_fadvise',lambda fd,start,length,flags:
                        calls.append(('drop_cache',start,length)))
    stream=sumstats._BlockStream.__new__(sumstats._BlockStream)
    stream._fd=-1;stream._offset=0;stream._writeback_started=0
    stream._writeback_waited=0;stream._writeback_bytes=10
    stream.stats=sumstats._StreamStats()
    for length in sizes:
        stream._offset+=length
        stream._advance_writeback()
    work=binary_writeback_work([{'bytes':n} for n in sizes],10)
    assert calls==[(a['kind'],a['offset'],a['bytes']) for row in work['events'] for a in row['actions']]
    assert stream._writeback_started==work['submitted_bytes']
    assert stream._writeback_waited==work['waited_bytes']


def test_disabled_and_close_only_output():
    work=binary_output_work(9,traits=1,chunk_markers=3,block_bytes=40,writeback_bytes=10)
    assert work['events_per_array']==[dict(chunk=3,bytes=36,pooled=True,at_close=True)]
    assert work['writeback_per_array']['submit_calls']==3
    disabled=binary_writeback_work([{'bytes':36}],10,False)
    assert disabled['submit_calls']==0
    assert disabled['not_explicitly_waited_before_fsync_bytes']==36
    assert binary_writeback_work([],0)['submit_calls']==0
    with pytest.raises(ValueError):binary_writeback_work([], -1)


def test_schedule_waits_before_releasing_buffer_and_counts_storage_once():
    from torchgwas.binary_schedule import BinaryWriterSchedule
    from torchgwas.execution_graph import ExecutionGraph
    work=binary_output_work(9,traits=1,chunk_markers=3,block_bytes=12,store_beta=False,writeback_bytes=10)
    service=dict(pagecache_seconds_per_byte=0.,storage_seconds_per_byte=1.,
                 submit_seconds=0.,wait_seconds=0.,fadvise_seconds=0.)
    writer=BinaryWriterSchedule(work,copy_seconds_per_byte=0.,zero_seconds_per_byte=0.,
        write_seconds_per_byte=99.,fsync_seconds_per_array=2.,writeback_bytes=10,writeback_service=service)
    g=ExecutionGraph();g.capacities={'output':1.}
    last=g.add('start')
    for i in range(3):last=writer.append(g,i,[last])
    end=writer.close(g,[last]);result=g.solve()
    assert result['end'][end]==pytest.approx(38.)
    assert writer.range_schedules['t'].bytes==36
    wait='writer:t:range:1:1:wait'
    first_storage='writer:t:range:0:0:storage'
    assert result['start'][wait]>=result['end'][first_storage]
    # The next os.write must follow the prior range wait, not just os.write.
    assert result['start']['writer:t:write:2']>=result['end'][wait]


def test_unpriced_writeback_remains_refused():
    from torchgwas.binary_schedule import BinaryWriterSchedule
    work=binary_output_work(9,chunk_markers=3,block_bytes=12,writeback_bytes=10)
    with pytest.raises(ValueError,match='explicit storage-writeback'):
        BinaryWriterSchedule(work,copy_seconds_per_byte=0.,zero_seconds_per_byte=0.,
            write_seconds_per_byte=1.,fsync_seconds_per_array=0.,writeback_bytes=10)


def test_two_array_writeback_shares_storage_and_conserves_payload():
    from torchgwas.binary_schedule import BinaryWriterSchedule
    from torchgwas.execution_graph import ExecutionGraph
    work=binary_output_work(9,traits=1,chunk_markers=3,block_bytes=12,writeback_bytes=10)
    writer=BinaryWriterSchedule(work,copy_seconds_per_byte=0.,zero_seconds_per_byte=0.,
        write_seconds_per_byte=0.,fsync_seconds_per_array=0.,writeback_bytes=10,
        writeback_service=dict(pagecache_seconds_per_byte=0.,storage_seconds_per_byte=1.,
                              submit_seconds=0.,wait_seconds=0.,fadvise_seconds=0.))
    g=ExecutionGraph();g.capacities={'output':1.}
    last=g.add('start')
    for i in range(3):last=writer.append(g,i,[last])
    end=writer.close(g,[last])
    assert g.solve()['end'][end]==pytest.approx(72.)
    assert sum(s.bytes for s in writer.range_schedules.values())==work['binary_payload_bytes']


def test_populated_cache_eviction_is_priced_only_for_waited_ranges():
    from torchgwas.binary_writeback_schedule import RangeWritebackSchedule
    from torchgwas.execution_graph import ExecutionGraph
    ledger=binary_writeback_work([{'bytes':35}],10)
    service=dict(pagecache_seconds_per_byte=0.,storage_seconds_per_byte=1.,
                 submit_seconds=0.,wait_seconds=0.,fadvise_seconds=2.,
                 fadvise_eviction_seconds_per_byte=.3)
    schedule=RangeWritebackSchedule(ledger,service)
    g=ExecutionGraph();start=g.add('start')
    schedule.after_write(g,'t',0,start)
    drop=[seconds for name,(seconds,deps) in g.nodes.items() if name.endswith(':drop')]
    assert drop==[5.,5.]  # two waited ranges; final range and tail are not advised
    service['fadvise_eviction_seconds_per_byte']=float('nan')
    with pytest.raises(ValueError,match='Invalid writeback service'):
        RangeWritebackSchedule(ledger,service)
