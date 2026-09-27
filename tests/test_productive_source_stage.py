"""Useful output may stage exact PGEN work without holding writer callbacks."""

import queue
import threading
import time
from types import SimpleNamespace

import pytest

from test_pgen_work_bounds import fixture, write_records
from torchgwas.initial_chunk_autotune import PublicInitialChunkTuning
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.productive_source_stage import (ProductiveSourceStage,
    validate_source_stage_config)
from torchgwas.sumstats import DenseWriteProgress, BinarySumstatsWriter
from torchgwas.sumstats_indexed import IndexedChunkWrite, IndexedOutputPartition
import numpy as np


def settings():
    return dict(records_per_step=4, max_steps=4, max_cpu_seconds=10.,
                max_window_seconds=10., max_retained_bytes=1 << 20,
                extra_host_reserve_bytes=2 << 20)


def test_each_written_event_earns_one_background_segment(tmp_path, monkeypatch):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    direct = header.schedule_bounds(0, 15, 2)
    original = header.schedule_bounds
    entered = threading.Event()
    release = threading.Event()
    calls = []

    def slow_first(*args, **kwargs):
        calls.append(args[:3])
        if len(calls) == 1:
            entered.set()
            assert release.wait(5.)
        return original(*args, **kwargs)

    monkeypatch.setattr(header, 'schedule_bounds', slow_first)
    stage = ProductiveSourceStage(lambda: header, header.input_identity,
                                  2, settings())
    try:
        stage.output_written()
        assert entered.wait(5.)
        # The first metadata step is held; later writer notifications return.
        for _ in range(3):
            stage.output_written()
        assert stage.snapshot()['pending_events'] == 3
    finally:
        release.set()
    deadline = time.monotonic() + 5.
    while not stage.snapshot()['complete'] and time.monotonic() < deadline:
        time.sleep(.01)
    audit = stage.finish()
    assert audit['complete'] and audit['stop_reason'] == 'complete'
    assert len(audit['steps']) == 4
    assert len(calls) == 4
    assert [row['cursor'] for row in audit['steps']] == [4, 8, 12, 15]
    assert audit['cpu_seconds'] > 0
    for key in direct:
        if key != 'scope':
            assert stage.result[key] == direct[key], key


def test_close_discards_queued_steps_after_current_segment(tmp_path, monkeypatch):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    original = header.schedule_bounds
    entered = threading.Event()
    release = threading.Event()
    calls = []

    def slow_first(*args, **kwargs):
        calls.append(args[:3])
        entered.set()
        assert release.wait(5.)
        return original(*args, **kwargs)

    monkeypatch.setattr(header, 'schedule_bounds', slow_first)
    stage = ProductiveSourceStage(lambda: header, header.input_identity,
                                  2, settings())
    stage.output_written()
    assert entered.wait(5.)
    stage.output_written()
    result = []
    closing = threading.Thread(target=lambda: result.append(stage.finish()))
    closing.start()
    deadline = time.monotonic() + 5.
    while not stage.closed and time.monotonic() < deadline:
        time.sleep(.01)
    assert stage.closed
    release.set()
    closing.join(5.)
    assert not closing.is_alive()
    audit = result[0]
    assert audit['closed'] and not audit['complete']
    assert audit['stop_reason'] == 'job_finished'
    assert len(calls) == 1
    stage.output_written()
    assert len(calls) == 1


def test_source_stage_limits_are_explicit():
    valid = settings()
    validate_source_stage_config(valid, 2)
    enlarged=dict(valid,max_cached_signatures=2048,max_cached_bounds=128,
                  extra_host_reserve_bytes=20<<20)
    validate_source_stage_config(enlarged,2)
    large=dict(valid,max_cached_signatures=32768,max_cached_bounds=256,
               extra_host_reserve_bytes=500<<20)
    validate_source_stage_config(large,2)
    for changed in (dict(valid, max_steps=33),
                    dict(valid, records_per_step=1),
                    dict(valid, max_cpu_seconds=float('nan')),
                    dict(valid, max_retained_bytes=3 << 20),
                    dict(valid, extra_host_reserve_bytes=0),
                    dict(enlarged,extra_host_reserve_bytes=2<<20),
                    dict(enlarged,max_cached_signatures=32769),
                    dict(large,extra_host_reserve_bytes=(500<<20)-1),
                    dict(enlarged,max_cached_bounds=257)):
        with pytest.raises(ValueError):
            validate_source_stage_config(changed, 2)


def test_worker_start_failure_does_not_escape_writer_callback(monkeypatch):
    class CannotStart:
        def __init__(self, **kwargs):
            pass
        def start(self):
            raise RuntimeError('worker unavailable')
    monkeypatch.setattr('torchgwas.productive_source_stage.threading.Thread',
                        CannotStart)
    stage = ProductiveSourceStage(lambda: None, {}, 2, settings())
    stage.output_written()
    audit = stage.finish()
    assert audit['stop_reason'] == 'worker_start_error'
    assert audit['start_error'] == 'RuntimeError: worker unavailable'
    assert not audit['complete']


def test_source_completion_callback_reads_ledger_once_and_contains_errors(tmp_path):
    path=tmp_path/'plain.pgen'
    write_records(path,8,[0]*7,[bytes([0,0])]*7)
    header=PgenHeaderWork(path)
    calls=[]
    done=threading.Event()
    def completed():
        calls.append(stage.completed_ledger()[1]['steps'])
        done.set()
        raise RuntimeError('optional evidence failed')
    stage=ProductiveSourceStage(lambda:header,header.input_identity,2,
        dict(settings(),records_per_step=4,max_steps=2),on_complete=completed)
    stage.output_written();stage.output_written()
    assert done.wait(5.)
    audit=stage.finish()
    assert audit['complete'] and audit['stop_reason']=='complete'
    assert calls==[2]
    assert audit['completion_callback_error']=='RuntimeError: optional evidence failed'


@pytest.mark.parametrize('change', ['none', 'frontier', 'binding',
                                    'rebase', 'budget'])
def test_live_staged_screen_runs_off_writer_and_rejects_stale_frontier(
        tmp_path, monkeypatch, change):
    path=tmp_path/'plain.pgen'
    write_records(path,8,[0]*31,[bytes([0,0])]*31)
    header=PgenHeaderWork(path)
    cfg=dict(chunk_size=4,window_markers=[8,16,24],
        budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.),
        source_staging=dict(settings(),records_per_step=8))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path=path,
        reduction=None,refresh=None,profile=dict(source_sha256='synthetic'))
    part=dict(id='only',device='cuda:0',variant_range=[0,31],trait_range=[0,3])
    start=dict(partitions=[part],chunk_sizes=[4,8],initial_size=4,
        context='fixture',
        input_file_identity=header.input_identity,
        memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    monkeypatch.setattr(tuner,'_check',lambda:None)
    entered=threading.Event();release=threading.Event();seen=[]
    def screen(frontier,stage,candidates,**options):
        assert stage.completed_ledger()[1]['steps']==4
        assert frontier['rectangles'][0]['variant_range']==(
            [20,31] if change=='rebase' and seen else [16,31])
        assert candidates==['priced-baseline']
        assert options['output_boundary']['valid']
        seen.append(threading.current_thread().name)
        entered.set()
        if len(seen)==1:
            assert release.wait(5.)
            if change=='budget':time.sleep(.02)
        return dict(stop_reason='complete',prediction_complete=False,
                    selection_validated=False)
    monkeypatch.setattr('torchgwas.productive_staged_screen.productive_staged_partial_screen',screen)
    monkeypatch.setattr('torchgwas.productive_staged_source_binding.audit_staged_source_price_binding',
        lambda *args,**kwargs:dict(status='declared_source_prices_verified'))
    monkeypatch.setattr('torchgwas.productive_staged_work_binding.audit_staged_work_price_binding',
        lambda *args,**kwargs:dict(status='declared_work_prices_verified'))
    options=dict(shared_source_capacities=dict(cpu=1.,dram=1e9,input=1e8),
        occupancy_scenario=None,max_candidates=2,max_partitions=2,
        max_unique_records=100,max_chunks_per_partition=100,
        max_cpu_seconds=2.,max_wall_seconds=5.)
    if change in ('rebase','budget'):
        options['max_rebases']=1
    if change=='budget':options['max_wall_seconds']=.001
    tuner.register_staged_screen(lambda frontier:['priced-baseline'],options)
    with pytest.raises(ValueError,match='once before first output'):
        tuner.register_staged_screen(lambda frontier:[],options)
    control=tuner.for_partition('cuda:0',[0,31],[0,3])
    for first in (0,4,8,12):
        assert control(first,31,8)==4
        tuner.output_written(DenseWriteProgress(first,first+4,(0,3),48,
            first+4,time.perf_counter(),'test','cuda:0'))
    assert entered.wait(5.)
    if change in ('frontier','rebase','budget'):assert control(16,31,8)==4
    if change=='binding':owner.profile['source_sha256']='changed'
    release.set()
    if change in ('rebase','budget'):
        tuner._staged_screen_worker.join(5.)
        assert not tuner._staged_screen_worker.is_alive()
    audit=tuner.finish(successful=False)
    evidence=audit['staged_screen_evidence']
    assert evidence['status']==dict(none='current_evidence',
        frontier='stale_frontier',binding='stale_binding',
        rebase='current_evidence',budget='stale_frontier')[change]
    assert evidence['screen']['stop_reason']=='complete'
    assert seen==['torchgwas-staged-screen']*(2 if change=='rebase' else 1)
    assert len(evidence['rebase_attempts'])==len(seen)
    if change=='budget':assert evidence['rebase_stop_reason']=='budget'
    assert not evidence['prediction_complete'] and not evidence['selection_validated']
    assert audit['source_staging']['complete']


def test_live_controller_runs_real_staged_partial_screen(tmp_path, monkeypatch):
    from test_productive_staged_screen import candidate
    path=tmp_path/'plain.pgen'
    write_records(path,8,[0]*31,[bytes([0,0])]*31)
    header=PgenHeaderWork(path)
    cfg=dict(chunk_size=4,window_markers=[8,16,24],
        budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.),
        source_staging=dict(settings(),records_per_step=8))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path=path,
        reduction=None,refresh=None,profile=dict(source_sha256='synthetic'))
    part=dict(id='only',device='cuda:0',variant_range=[0,31],trait_range=[0,3])
    start=dict(partitions=[part],chunk_sizes=[4,8],initial_size=4,
        context='fixture',
        input_file_identity=header.input_identity,
        memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    monkeypatch.setattr(tuner,'_check',lambda:None)
    schedule=header.schedule_bounds(16,31,4)
    source_profile=dict(decode_units={key:1e-8 for key in schedule['source_units']},
        cpu_fraction=.5,depth=2,decode_workers=2,cpu_available_cores=2.,
        shared_dram_bytes_per_second=1e8,read_bytes_per_second=1e7,
        input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-8,
                                   cpu_seconds_per_call=1e-5))
    owner.profile=dict(source_sha256='synthetic',contexts=[dict(
        name='fixture',profiles={'cuda:0':source_profile},
        shared_capacities=dict(cpu=2.,dram=1e8,input=1e7))])
    def factory(frontier):
        row=candidate(None,source_profile,4,False)
        row['partitions'][0]['variant_range']=frontier['rectangles'][0]['variant_range']
        row['partitions'][0]['trait_range']=[0,3]
        return [row]
    tuner.register_staged_screen(factory,dict(
        shared_source_capacities=dict(cpu=1.,dram=5e7,input=5e6),
        occupancy_scenario=None,max_candidates=2,max_partitions=2,
        max_unique_records=100,max_chunks_per_partition=100,
        max_cpu_seconds=5.,max_wall_seconds=10.))
    control=tuner.for_partition('cuda:0',[0,31],[0,3])
    for first in (0,4,8,12):
        assert control(first,31,8)==4
        tuner.output_written(DenseWriteProgress(first,first+4,(0,3),48,
            first+4,time.perf_counter(),'test','cuda:0'))
    deadline=time.monotonic()+5.
    while tuner._staged_screen_worker is None and time.monotonic()<deadline:
        time.sleep(.01)
    assert tuner._staged_screen_worker is not None
    evidence=tuner.finish(successful=False)['staged_screen_evidence']
    assert evidence['status']=='source_prices_unbound'
    assert evidence['source_price_binding']['status']=='source_prices_match_but_unbound'
    assert evidence['work_price_binding']['status']=='work_prices_match_but_unbound'
    screen=evidence['screen']
    assert screen['evaluated_candidates']==1
    assert screen['candidates'][0]['partial']['envelope']['coverage']['exact_coverage']
    assert screen['candidates'][0]['source_schedule_method']=='staged_primary_rebase'
    assert screen['prior_stage_cost_once']['steps']==4


def test_configured_native_chunk_screen_waits_for_output_and_uses_live_writer(
        tmp_path, monkeypatch):
    from test_productive_staged_work_binding import case as work_case
    path=tmp_path/'plain.pgen'
    write_records(path,8,[0]*31,[bytes([0,0])]*31)
    header=PgenHeaderWork(path)
    schedule=header.schedule_bounds(16,31,4)
    source=dict(decode_units={key:1e-8 for key in schedule['source_units']},
        cpu_fraction=.5,depth=2,decode_workers=2,cpu_available_cores=2.,
        shared_dram_bytes_per_second=1e8,read_bytes_per_second=1e7)
    profile,_=work_case()
    context=profile['contexts'][0]
    context['name']='fixture'
    context['shared_capacities'].update(cpu=2.,dram=1e8,input=1e7)
    active=context['profiles']['cuda:0']
    active.update(source)
    active['host_primitives']=dict(synthetic=.001)
    active['kernel_geometry']=[dict(N=8,B=b,K=3,C=3,kernels=[])
                               for b in (4,8)]
    cfg=dict(chunk_size=4,partition_axis='trait',window_markers=[8,16,24],
        budget=dict(max_steps=1,max_cpu_seconds=1.,max_window_seconds=10.),
        source_staging=dict(settings(),records_per_step=8),
        staged_screen=dict(chunk_sizes=[4,8],occupancy_scenario=None,
            max_partitions=2,max_unique_records=100,
            max_chunks_per_partition=100,max_cpu_seconds=5.,
            max_wall_seconds=10.,max_rebases=0))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path=path,
        reduction=None,refresh=None,profile=profile)
    part=dict(id='only',device='cuda:0',variant_range=[0,31],trait_range=[0,3])
    start=dict(partitions=[part],chunk_sizes=[4,8],initial_size=4,
        context='fixture',input_file_identity=header.input_identity,
        candidate=dict(tiles=[dict(data=dict(covariates=3))],
            output=dict(block_bytes=48,queue_depth=2,store_beta=True,fsync=True)),
        memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    monkeypatch.setattr(tuner,'_check',lambda:None)
    assert tuner._staged_screen_request is not None
    assert tuner._staged_screen_worker is None
    captured=[]
    def screen(frontier,stage,candidates,**options):
        assert frontier['rectangles'][0]['variant_range']==[16,31]
        assert [row['chunk_markers'] for row in candidates]==[4,8]
        assert all(row['partitions']==frontier['rectangles'] for row in candidates)
        assert candidates[0]['output_options']['dense_writer_options'][
            'store_variant_df'] is False
        captured.append(threading.current_thread().name)
        return dict(stop_reason='complete',prediction_complete=False,
                    selection_validated=False)
    monkeypatch.setattr('torchgwas.productive_staged_screen.productive_staged_partial_screen',screen)
    writer=BinarySumstatsWriter(tmp_path/'output',31,['a','b','c'],8,3,
        block_bytes=48,queue_depth=2,on_write_progress=tuner.output_written)
    tuner.register_writer(writer,'cuda:0')
    try:
        control=tuner.for_partition('cuda:0',[0,31],[0,3])
        values=np.ones((4,3),np.float32)
        for first in (0,4,8,12):
            assert control(first,31,8)==4
            writer.write_chunk(first,first+4,values,values)
        deadline=time.monotonic()+5.
        while tuner._staged_screen_worker is None and time.monotonic()<deadline:
            time.sleep(.01)
        with tuner._lock:worker=tuner._staged_screen_worker
        assert worker is not None
        worker.join(5.)
        assert not worker.is_alive()
        evidence=tuner.finish(successful=False)['staged_screen_evidence']
        assert captured==['torchgwas-staged-screen'],evidence
        assert evidence['status']=='source_prices_unbound'
        assert evidence['screen']['stop_reason']=='complete'
    finally:
        writer.abort()
        tuner.finish(successful=False)


def test_public_writer_events_complete_real_staged_source(tmp_path):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    cfg = dict(chunk_size=4, window_markers=[8, 16, 24],
        budget=dict(max_steps=1, max_cpu_seconds=1., max_window_seconds=10.),
        source_staging=dict(settings(),max_cached_signatures=2048,
            max_cached_bounds=128,extra_host_reserve_bytes=20<<20))
    owner = SimpleNamespace(config=dict(initial_chunks=cfg), input_path=path,
                            reduction=None, refresh=None)
    partition = dict(id='0', device='cuda:0', variant_range=[0, 15],
                     trait_range=[0, 3])
    start = dict(partitions=[partition], chunk_sizes=[4, 8],
                 initial_size=4, input_file_identity=header.input_identity,
                 memory=dict(retained_index_bases_bytes=0))
    tuner = PublicInitialChunkTuning(owner, start, {}, {})
    bound = tuner.for_partition('cuda:0', [0, 15], [0, 3])
    for number, first in enumerate((0, 4, 8, 12), 1):
        width = bound(first, 15, 8)
        last = first + width
        tuner.output_written(DenseWriteProgress(first, last, (0, 3),
            width * 3 * 4, last, time.perf_counter(), 'test', 'cuda:0'))
        deadline = time.monotonic() + 5.
        while len(tuner._source_stage.snapshot()['steps']) < number and time.monotonic() < deadline:
            time.sleep(.01)
        assert len(tuner._source_stage.snapshot()['steps']) == number
    assert tuner.run.snapshot()['current_chunk_size'] == 4
    assert tuner._source_stage.result is not None
    staged_header=tuner._source_stage.stage.header
    assert staged_header.cache_info()['maxsize']==2048
    assert staged_header.bounds_cache_info()['max_entries']==128
    assert tuner._source_stage.snapshot()['header_cache']['signatures']['maxsize']==2048
    direct = header.schedule_bounds(0, 15, 4)
    for key in direct:
        if key != 'scope':
            assert tuner._source_stage.result[key] == direct[key], key
    audit = tuner.finish(successful=False)
    assert audit['source_staging']['complete']
    assert audit['source_staging']['stop_reason'] == 'complete'


def test_staged_step_captures_and_prices_live_writer_without_callback_planning(tmp_path):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    cfg = dict(chunk_size=4, window_markers=[8, 16, 24],
        budget=dict(max_steps=1, max_cpu_seconds=1., max_window_seconds=10.),
        source_staging=settings())
    profile = dict(cpu_fraction=.5,
        writer_copy_service=dict(cpu_seconds_per_byte=0.,cpu_seconds_per_call=0.),
        process_units=dict(bytearray_zero_bytes=0.),executor_cpu_seconds=0.,
        writeback_service=dict(pagecache_seconds_per_byte=.1,
            storage_seconds_per_byte=.2,submit_seconds=0.,wait_seconds=0.,
            fadvise_seconds=0.),fsync_seconds=0.)
    owner = SimpleNamespace(config=dict(initial_chunks=cfg), input_path=path,
        reduction=None, refresh=None,
        profile=dict(contexts=[dict(name='test',profiles={'cuda:0':profile})]))
    part = dict(id='only',device='cuda:0',variant_range=[0,15],trait_range=[0,2])
    start = dict(partitions=[part],chunk_sizes=[4,8],initial_size=4,
                 context='test',input_file_identity=header.input_identity,
                 memory=dict(retained_index_bases_bytes=0))
    tuner = PublicInitialChunkTuning(owner,start,{}, {})
    writer = BinarySumstatsWriter(tmp_path/'output',15,['a','b'],129,123,
        block_bytes=32,writeback_bytes=0,on_write_progress=tuner.output_written)
    tuner.register_writer(writer,'cuda:0')
    try:
        control=tuner.for_partition('cuda:0',[0,15],[0,2])
        assert control(0,15,8)==4
        values=np.ones((4,2),np.float32)
        writer.write_chunk(0,4,values,values)
        deadline=time.monotonic()+5.
        while not tuner._source_stage.snapshot()['steps'] and time.monotonic()<deadline:
            time.sleep(.01)
        steps=tuner._source_stage.snapshot()['steps']
        assert len(steps)==1
        observed=steps[0]['writer_observation']
        assert observed['observation_valid']
        assert observed['bracket']['os_write_bytes_upper']==0
        assert observed['priced_service']['kind']=='torchgwas.priced_bracketed_dense_writer_queues.v1'
        assert tuner.run.revision_token()['current_chunk_size']==4
    finally:
        writer.abort()
        tuner.finish(successful=False)


def test_oversized_writer_observation_is_dropped_without_stopping_source(tmp_path):
    path=tmp_path/'mixed.pgen'
    fixture(path,129)
    header=PgenHeaderWork(path)
    stage=ProductiveSourceStage(lambda:header,header.input_identity,2,settings(),
        on_step_observation=lambda:dict(payload=bytearray(2<<20)))
    stage.output_written()
    deadline=time.monotonic()+5.
    while not stage.snapshot()['steps'] and time.monotonic()<deadline:
        time.sleep(.01)
    audit=stage.finish()
    assert len(audit['steps'])==1
    assert audit['steps'][0]['writer_observation']==dict(status='omitted_retained_budget')
    assert audit['steps'][0]['error'] is None
    assert audit['stop_reason']=='job_finished'


def test_staged_jagwas_step_binds_live_queue_to_retained_first_chunks(tmp_path):
    path = tmp_path / 'mixed.pgen'
    fixture(path, 129)
    header = PgenHeaderWork(path)
    cfg = dict(chunk_size=3, window_markers=[6, 12, 18],
        budget=dict(max_steps=1, max_cpu_seconds=1., max_window_seconds=10.),
        source_staging=dict(settings(), max_steps=5))
    owner = SimpleNamespace(config=dict(initial_chunks=cfg), input_path=path,
        reduction='jagwas', refresh=None)
    parts = [dict(id='left', device='cuda:0', variant_range=[0, 6],
                  trait_range=[0, 3]),
             dict(id='right', device='cuda:1', variant_range=[6, 15],
                  trait_range=[0, 3])]
    start = dict(partitions=parts, chunk_sizes=[3, 6], initial_size=3,
                 input_file_identity=header.input_identity,
                 memory=dict(retained_index_bases_bytes=0))
    tuner = PublicInitialChunkTuning(owner, start, {}, {})
    source = queue.Queue(maxsize=2)
    tuner.register_indexed_result_queue(source, [(0, 6), (6, 15)],
                                        ['cuda:0', 'cuda:1'], object())
    left = tuner.for_partition('cuda:0', (0, 6), (0, 3))
    right = tuner.for_partition('cuda:1', (6, 15), (0, 3))
    assert left(0, 6, 6) == 3
    assert right(6, 15, 6) == 3
    source.put((6, 9, None, np.arange(3, dtype=np.float64), None, None))
    now = time.perf_counter()
    tuner.output_written(IndexedChunkWrite(0, 3, 'jagwas', 0, 0, None,
        now, now, False, IndexedOutputPartition('cuda:0', (0, 6), (0, 3)),
        (0, 3)))
    deadline = time.monotonic() + 5.
    while not tuner._source_stage.snapshot()['steps'] and time.monotonic() < deadline:
        time.sleep(.01)
    steps = tuner._source_stage.snapshot()['steps']
    assert len(steps) == 1
    observed = steps[0]['writer_observation']
    assert observed['checkpoint_valid']
    assert observed['joined_issued_queue']['queued_chunks'] == 1
    assert observed['joined_issued_queue']['queued_source_chunks'][0]['partition_id'] == 'right'
    assert tuner.run.snapshot()['prefix_complete']
    assert tuner.run.snapshot()['stop_reason'] is None
    tuner.finish(successful=False)
