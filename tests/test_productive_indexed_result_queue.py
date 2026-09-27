"""Bounded read-only JAGWAS result-queue observations."""
import queue
import threading
import time
from types import SimpleNamespace

import numpy as np
import pytest

from torchgwas.initial_chunk_autotune import PublicInitialChunkTuning
from torchgwas.indexed_writer_progress import IndexedWriterProgress
from torchgwas.productive_indexed_result_queue import snapshot_indexed_result_queue
from torchgwas.productive_indexed_queue_join import bind_jagwas_result_queue
from torchgwas.sumstats_indexed import (IndexedChunkWrite,
    IndexedOutputPartition, write_indexed_sumstats)


def test_queue_snapshot_counts_owned_payload_and_finished_sentinel():
    source=queue.Queue(maxsize=2);finished=object()
    source.put((0,4,None,np.arange(4,dtype=np.float64),None,None))
    source.put(finished)
    observed=snapshot_indexed_result_queue(source,[(0,4),(4,8)],
        ['cuda:0','cuda:1'],finished)
    assert observed['capacity']==2 and observed['queued_items']==2
    assert observed['queued_results']==1 and observed['queued_sentinels']==1
    assert observed['resident_array_bytes']==32
    assert observed['results']==[dict(device='cuda:0',variant_range=[0,4],
                                      resident_array_bytes=32)]
    assert observed['capture_started_seconds']<=observed['capture_anchor_seconds']<=observed['capture_finished_seconds']
    with pytest.raises(ValueError,match='partition'):
        snapshot_indexed_result_queue(source,[(0,2),(2,8)],
            ['cuda:0','cuda:1'],finished)


def test_productive_jagwas_registers_and_observes_live_shared_queue():
    parts=[dict(id='a',device='cuda:0',variant_range=[0,4],trait_range=[0,2]),
           dict(id='b',device='cuda:1',variant_range=[4,8],trait_range=[0,2])]
    cfg=dict(chunk_size=4,window_markers=[4,8,12],
             budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction='jagwas',
                          refresh=None)
    tuner=PublicInitialChunkTuning(owner,dict(partitions=parts,
        chunk_sizes=[4],initial_size=4,memory={}),{}, {})
    source=queue.Queue(maxsize=2);finished=object()
    tuner.register_indexed_result_queue(source,[(0,4),(4,8)],
                                        ['cuda:0','cuda:1'],finished)
    source.put((4,8,None,np.arange(4,dtype=np.float64),None,None))
    observed=tuner.capture_indexed_result_queue()
    assert observed['observation_valid'] and observed['stable_revision']
    assert not observed['checkpoint_valid']  # No useful output or bound issue prefix yet.
    assert observed['observation']['results'][0]['device']=='cuda:1'
    staged=tuner._source_writer_observation()
    assert staged['observation_valid'] and staged['observation']['resident_array_bytes']==32
    tuner._planning_done=True
    assert tuner.for_partition('cuda:0',(0,4),(0,2))(0,4,4)==4
    assert tuner.for_partition('cuda:1',(4,8),(0,2))(4,8,4)==4
    now=time.perf_counter()
    tuner.output_written(IndexedChunkWrite(0,4,'jagwas',0,0,None,now,now,
        False,IndexedOutputPartition('cuda:0',(0,4),(0,2)),(0,4)))
    bound=tuner.capture_indexed_result_queue()
    assert bound['checkpoint_valid'] and bound['issued_output_boundary']['valid']
    assert bound['joined_issued_queue']['queued_chunks']==1
    assert bound['joined_issued_queue']['queued_source_chunks'][0]['partition_id']=='b'
    assert bound['joined_issued_queue']['issued_not_part_written_outside_queue']==[]
    bad=dict(bound['observation'],results=[dict(bound['observation']['results'][0],device='cuda:0')])
    with pytest.raises(ValueError,match='unfinished issued'):
        bind_jagwas_result_queue(bound['issued_output_boundary'],bad)
    tuner.finish(successful=False)


def test_queue_anchor_binds_stably_active_indexed_writer_chunk():
    parts=[dict(id='left',device='cuda:0',variant_range=[0,8],trait_range=[0,2]),
           dict(id='right',device='cuda:1',variant_range=[8,16],trait_range=[0,2])]
    cfg=dict(chunk_size=4,window_markers=[4,8,12],
             budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction='jagwas',
                          refresh=None)
    tuner=PublicInitialChunkTuning(owner,dict(partitions=parts,
        chunk_sizes=[4],initial_size=4,memory={}),{}, {})
    source=queue.Queue(maxsize=2);finished=object()
    tuner.register_indexed_result_queue(source,[(0,8),(8,16)],
                                        ['cuda:0','cuda:1'],finished)
    writer=IndexedWriterProgress('jagwas')
    tuner.register_indexed_writer(writer)
    left=tuner.for_partition('cuda:0',(0,8),(0,2))
    right=tuner.for_partition('cuda:1',(8,16),(0,2))
    assert left(0,8,4)==4 and left(4,8,4)==4
    assert right(8,16,4)==4
    source.put((8,12,None,np.arange(4,dtype=np.float64),None,None))
    tuner._planning_done=True
    now=time.perf_counter()
    tuner.output_written(IndexedChunkWrite(0,4,'jagwas',0,0,None,now,now,
        False,IndexedOutputPartition('cuda:0',(0,8),(0,2)),(0,4)))
    writer.begin_chunk((4,8),IndexedOutputPartition('cuda:0',(0,8),(0,2)))
    observed=tuner.capture_indexed_result_queue()
    assert observed['checkpoint_valid']
    assert observed['indexed_writer_bracket']['active_at_queue_anchor']['variant_range']==[4,8]
    joined=observed['joined_issued_queue']
    assert joined['queued_chunks']==1 and joined['completed_producer_chunks']==2
    assert joined['active_writer_source_chunk']['variant_range']==[4,8]
    assert joined['issued_not_part_written_upstream_or_unresolved']==[]
    writer.begin_emit();writer.complete_chunk(0,0,False)
    tuner.finish(successful=False)


def test_queue_anchor_spans_real_active_indexed_part_write(tmp_path, monkeypatch):
    import torchgwas.sumstats_indexed as indexed
    parts=[dict(id='left',device='cuda:0',variant_range=[0,8],trait_range=[0,2]),
           dict(id='right',device='cuda:1',variant_range=[8,16],trait_range=[0,2])]
    cfg=dict(chunk_size=4,window_markers=[4,8,12],
             budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction='jagwas',
                          refresh=None)
    tuner=PublicInitialChunkTuning(owner,dict(partitions=parts,
        chunk_sizes=[4],initial_size=4,memory={}),{}, {})
    source=queue.Queue(maxsize=2);finished=object()
    tuner.register_indexed_result_queue(source,[(0,8),(8,16)],
                                        ['cuda:0','cuda:1'],finished)
    writer=IndexedWriterProgress('jagwas')
    tuner.register_indexed_writer(writer)
    left=tuner.for_partition('cuda:0',(0,8),(0,2))
    right=tuner.for_partition('cuda:1',(8,16),(0,2))
    assert left(0,8,4)==4 and left(4,8,4)==4
    assert right(8,16,4)==4
    tuner._planning_done=True
    now=time.perf_counter()
    tuner.output_written(IndexedChunkWrite(0,4,'jagwas',0,0,None,now,now,
        False,IndexedOutputPartition('cuda:0',(0,8),(0,2)),(0,4)))
    source.put((8,12,None,np.arange(4,dtype=np.float64),None,None))
    entered=threading.Event();release=threading.Event()
    original=indexed.np.savez
    def slow_savez(*args,**kwargs):
        entered.set()
        assert release.wait(5.)
        return original(*args,**kwargs)
    monkeypatch.setattr(indexed.np,'savez',slow_savez)
    failures=[]
    def write():
        try:
            write_indexed_sumstats(tmp_path,[str(i) for i in range(16)],
                ['a','b'],129,[(4,8,None,np.arange(4,dtype=np.float64),None)],
                kind='jagwas',df=100,chi2_df=2,fsync=False,
                on_chunk_written=tuner.output_written,
                partition_for_range=lambda *_:IndexedOutputPartition(
                    'cuda:0',(0,8),(0,2)),live_progress=writer)
        except BaseException as error:failures.append(error)
    thread=threading.Thread(target=write)
    thread.start()
    try:
        assert entered.wait(5.)
        captured=tuner.capture_indexed_result_queue()
        assert captured['checkpoint_valid']
        assert captured['joined_issued_queue']['completed_producer_chunks']==2
        assert captured['joined_issued_queue']['active_writer_source_chunk'][
            'variant_range']==[4,8]
    finally:
        release.set();thread.join(5.)
    assert not thread.is_alive() and not failures
    assert writer.snapshot()['phase']=='published'
    tuner.finish(successful=False)


def test_significant_tile_writer_capture_preserves_producer_identity():
    parts=[dict(id='low',device='cuda:0',variant_range=[0,4],trait_range=[0,2]),
           dict(id='high',device='cuda:1',variant_range=[0,4],trait_range=[2,4])]
    cfg=dict(chunk_size=2,window_markers=[2,4,6],
             budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction='significant',
                          refresh=None)
    tuner=PublicInitialChunkTuning(owner,dict(partitions=parts,
        chunk_sizes=[2],initial_size=2,memory={}),{}, {})
    writer=IndexedWriterProgress('significant')
    tuner.register_indexed_writer(writer)
    assert tuner.for_partition('cuda:1',(0,4),(2,4))(0,4,2)==2
    writer.begin_chunk((0,2),IndexedOutputPartition('cuda:1',(0,4),(2,4)))
    observed=tuner.capture_significant_indexed_writer()
    assert observed['observation_valid']
    assert observed['writer']['active']==dict(device='cuda:1',
        variant_range=[0,2],trait_range=[2,4])
    assert observed['status']=='indexed_writer_observed_host_selector_inflight_unobserved'
    writer.begin_emit();writer.complete_chunk(0,0,False)
    tuner.finish(successful=False)
