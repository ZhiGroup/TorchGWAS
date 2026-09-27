from concurrent.futures import ThreadPoolExecutor
from dataclasses import FrozenInstanceError
import threading
import time
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.sumstats_indexed import IndexedChunkWrite, write_indexed_sumstats


def run_state(**kw):
    return ProductiveTuningRun([
        dict(id='a', device='cuda:0', variant_range=[0,64], trait_range=[0,3]),
        dict(id='b', device='cuda:1', variant_range=[64,129], trait_range=[0,3])],
        chunk_sizes=[4,8,16], initial=4, **kw)


def event(*, rows=4, fsync=True):
    now=time.perf_counter()
    return IndexedChunkWrite(0,4,'jagwas',rows,100 if rows else 0,
        'part_000000.npz' if rows else None,now,now,bool(rows and fsync))


def proposal(size=16, before=10., after=5.):
    return dict(chunk_size=size, baseline_seconds=before, candidate_seconds=after)


def test_released_planner_never_applies_a_stale_source_frontier():
    entered=threading.Event();resume=threading.Event();result=[]
    run=run_state(budget=IncrementalPlanningBudget(max_cpu_seconds=10.,max_window_seconds=30.))
    run.output_written(event())
    def build(snapshot):
        assert snapshot['issued_revision']==0
        entered.set()
        assert resume.wait(5.)
        return proposal()
    worker=threading.Thread(target=lambda:result.append(run.planning_step(
        build,release_issue_frontier=True,**FORECAST)))
    worker.start()
    try:
        assert entered.wait(5.)
        # The reader can issue while the analytical build is still running.
        assert run.for_partition('a')(0,64,run.capacity)==4
    finally:
        resume.set();worker.join(5.)
    assert not worker.is_alive() and len(result)==1
    assert result[0]['evaluated'] and result[0]['stale_frontier']
    assert not result[0]['applied'] and result[0]['issued_revision']==1
    assert result[0]['planned_issued_revision']==0
    assert run.snapshot()['current_chunk_size']==4
    assert run.snapshot()['stop_reason'] is None
    run.finish(successful=False)


def test_released_planner_applies_only_at_unchanged_frontier():
    run=run_state(budget=IncrementalPlanningBudget(max_cpu_seconds=10.,max_window_seconds=30.))
    run.output_written(event())
    result=run.planning_step(lambda snapshot: proposal(),release_issue_frontier=True,**FORECAST)
    assert result['evaluated'] and result['applied']
    assert result['planned_issued_revision']==result['issued_revision']==0
    assert run.for_partition('a')(0,64,run.capacity)==16
    run.finish(successful=False)


FORECAST=dict(remaining_seconds=10., expected_cpu_seconds=.001, expected_wall_seconds=.002,
              expected_gain_seconds=5., switching_seconds=.001,publication_seconds=.06)


def test_no_planner_before_written_work_and_empty_output_starts_window():
    run=run_state()
    result=run.planning_step(lambda state: pytest.fail('premature planning'), **FORECAST)
    assert not result['evaluated'] and result['reason']=='no_useful_output_yet'
    assert run.snapshot()['first_written'] is None
    run.output_written(event(rows=0))
    snap=run.snapshot()
    assert snap['first_written'] is not None and snap['first_material_part'] is None
    assert snap['first_fsynced_part'] is None and snap['planning']['started']


def finish_source(run):
    for row in run.snapshot()['partitions']:
        reserve=run.for_partition(row['id']);cursor=row['cursor'];end=row['variant_range'][1]
        while cursor<end:cursor+=reserve(cursor,end,run.capacity)


def test_structural_reuse_is_lazy_charged_and_published_only_after_success(tmp_path):
    from torchgwas.tensor_work import eager_statistics_work
    from torchgwas.structural_tensor_cache import StructuralTensorWorkCache
    # This tests persistence/lifecycle, not cold PyTorch trace latency. Budget
    # overruns are checked separately with a controlled clock below.
    options=dict(structural_cache_dir=tmp_path,budget=IncrementalPlanningBudget(max_cpu_seconds=60.,max_window_seconds=120.))
    forecasts=dict(FORECAST,remaining_seconds=100.,expected_gain_seconds=50.)
    run=run_state(**options)
    with patch('torchgwas.productive_run.StructuralTensorWorkCache',wraps=StructuralTensorWorkCache) as constructor:
        run.planning_step(lambda state:pytest.fail('before output'),**forecasts)
        assert constructor.call_count==0 and run.snapshot()['structural_cache']['state'] is None
        run.output_written(event())
        def build(state):
            assert state['structural_cache']['state'] is not None
            eager_statistics_work(32,4,3,8,True)
            return proposal(before=100.,after=50.)
        result=run.planning_step(build,**forecasts)
        assert constructor.call_count==1 and result['applied'],result
        assert result['cpu_seconds']>0 and not list(tmp_path.glob('*/*.json'))
    assert run.snapshot()['structural_cache']['state']['pending']==1
    finish_source(run);first=run.finish(successful=True)
    assert first['structural_cache']['state']['closed']
    assert len(first['structural_cache']['publication']['stored'])==1
    before={str(p):p.read_bytes() for p in tmp_path.glob('*/*.json')}
    later=run_state(structural_cache_dir=tmp_path,budget=IncrementalPlanningBudget(max_cpu_seconds=60.,max_window_seconds=120.))
    later.output_written(event());later.planning_step(build,**forecasts)
    assert later.snapshot()['structural_cache']['state']['disk_hits']==1
    finish_source(later);final=later.finish(successful=True)
    assert final['structural_cache']['publication']['stored']==[]
    assert {str(p):p.read_bytes() for p in tmp_path.glob('*/*.json')}==before
    assert later.finish(successful=True)==final


def test_failed_run_discards_structural_work_and_publication_error_preserves_output(tmp_path):
    from torchgwas.tensor_work import eager_statistics_work
    def build(state):eager_statistics_work(32,4,3,8,True);return proposal()
    failed=run_state(structural_cache_dir=tmp_path,budget=IncrementalPlanningBudget(max_cpu_seconds=1.))
    failed.output_written(event());failed.planning_step(build,**FORECAST)
    report=failed.finish(successful=False)
    assert report['structural_cache']['publication']['status']=='unsuccessful'
    assert report['structural_cache']['state']['closed'] and not list(tmp_path.glob('*/*.json'))
    complete=run_state(structural_cache_dir=tmp_path,budget=IncrementalPlanningBudget(max_cpu_seconds=1.))
    complete.output_written(event());complete.planning_step(build,**FORECAST);finish_source(complete)
    with patch.object(complete._structural_cache,'publish',side_effect=OSError('unavailable')):
        report=complete.finish(successful=True)
    assert report['finished']['successful'] and report['structural_cache']['state']['closed']
    assert report['structural_cache']['publication']['status']=='cache_error'


def test_structural_initialization_overrun_is_charged_and_cannot_apply(tmp_path):
    from torchgwas.structural_tensor_cache import StructuralTensorWorkCache
    from test_planning_session import Clock
    from torchgwas import planning_session
    clock=Clock()
    def slow_create(directory):
        clock.advance(.2,.2)
        return StructuralTensorWorkCache(directory)
    with patch.object(planning_session.time,'perf_counter',lambda:clock.wall),patch.object(planning_session.time,'thread_time',lambda:clock.cpu):
        run=run_state(structural_cache_dir=tmp_path)
        run.output_written(event())
        with patch('torchgwas.productive_run.StructuralTensorWorkCache',side_effect=slow_create):
            result=run.planning_step(lambda state:proposal(),**FORECAST)
        assert result['cpu_seconds']==.2 and not result['applied']
        assert run.snapshot()['current_chunk_size']==4
        run.finish(successful=False)


def test_structural_persistence_requires_publication_cost_before_any_cache_io(tmp_path):
    run=run_state(structural_cache_dir=tmp_path);run.output_written(event())
    options={key:value for key,value in FORECAST.items() if key!='publication_seconds'}
    result=run.planning_step(lambda state:pytest.fail('missing publication estimate'),**options)
    assert not result['evaluated'] and result['reason']=='publication_cost_required'
    assert run.snapshot()['structural_cache']['state'] is None
    assert not list(tmp_path.glob('*/*.json'))
    assert run.for_partition('a')(0,64,16)==4
    run.finish(successful=False)


def test_publication_and_prior_planning_are_charged_when_applying_a_proposal():
    from test_planning_session import Clock
    from torchgwas import planning_session
    clock=Clock()
    with patch.object(planning_session.time,'perf_counter',lambda:clock.wall),patch.object(planning_session.time,'thread_time',lambda:clock.cpu):
        run=run_state(budget=IncrementalPlanningBudget(max_cpu_seconds=1.));run.output_written(event())
        options=dict(FORECAST,expected_gain_seconds=None,publication_seconds=.15,reserve_seconds=.05)
        def build(state):
            clock.advance(.1,.001)
            return proposal(before=1.,after=.75)
        first=run.planning_step(build,**options)
        assert not first['applied'] and first['total_tuning_cost_seconds']==pytest.approx(.301)
        second=run.planning_step(build,**dict(options,publication_seconds=0.,reserve_seconds=.05))
        # The latest 100 ms step alone would appear to repay a 250 ms gain.
        # Both early steps and the reserve cost exhaust that gain.
        assert not second['applied'] and second['total_tuning_cost_seconds']==pytest.approx(.251)
        assert second['cumulative_planning_wall_seconds']==pytest.approx(.2)
        assert run.for_partition('a')(0,64,16)==4
        run.finish(successful=False)


def test_exact_prefix_includes_unsampled_reserved_reads_and_preserves_old_sizes():
    run=run_state();a=run.for_scan('cuda:0',(0,64),3);b=run.for_scan('cuda:1',(64,129),3)
    assert [a(0,64,16),a(4,64,16),b(64,129,16)]==[4,4,4]
    run.output_written(event())
    original=run.snapshot()
    seen=[]
    result=run.planning_step(lambda state: seen.append(state) or proposal(), **FORECAST)
    assert result['applied'] and seen[0]['issued_revision']==3
    assert seen[0]['partitions'][0]['ranges']==[[0,4],[4,8]]
    assert a(8,64,16)==16 and b(68,129,16)==16
    assert run.snapshot()['partitions'][0]['ranges']==[[0,4],[4,8],[8,24]]
    assert original['partitions'][0]['ranges']==[[0,4],[4,8]]
    seen[0]['partitions'][0]['ranges'][0][0]=1000
    assert run.snapshot()['partitions'][0]['ranges'][0][0]==0


def test_step_holds_issue_frontier_and_does_not_drop_a_waiting_read():
    run=run_state(budget=IncrementalPlanningBudget(max_cpu_seconds=1., max_window_seconds=5.))
    a=run.for_partition('a');a(0,64,16);run.output_written(event())
    attempting=threading.Event();issued=threading.Event()
    def read():
        attempting.set();value=a(4,64,16);issued.set();return value
    with ThreadPoolExecutor(max_workers=1) as pool:
        def build(state):
            future=pool.submit(read)
            assert attempting.wait(1.) and not issued.is_set()
            assert state['partitions'][0]['cursor']==4
            return proposal()
        result=run.planning_step(build, **FORECAST)
        assert result['applied'] and issued.wait(1.)
    assert run.snapshot()['partitions'][0]['ranges']==[[0,4],[4,20]]


@pytest.mark.parametrize('value', [proposal(8,10,10),proposal(8,10,11),proposal(8,10,9.9999)])
def test_unprofitable_proposal_never_changes_future_work(value):
    run=run_state();run.output_written(event())
    result=run.planning_step(lambda state:value, **FORECAST)
    assert result['evaluated'] and not result['applied']
    assert run.for_partition('a')(0,64,16)==4


@pytest.mark.parametrize('value', [proposal(7),proposal(before=11),proposal(after=float('nan')),
                                  dict(chunk_size=8),proposal(size=True)])
def test_invalid_planning_stops_tuning_but_does_not_discard_scan(value):
    run=run_state();run.output_written(event())
    result=run.planning_step(lambda state:value, **FORECAST)
    assert result['error']=='ValueError' and not result['applied']
    assert run.snapshot()['stop_reason']=='planning_error'
    assert run.for_partition('a')(0,64,16)==4
    for _ in range(100):
        assert not run.planning_step(lambda state:pytest.fail('closed'), **FORECAST)['evaluated']
    assert len(run.snapshot()['decisions'])==1


def test_full_source_prefix_and_prefix_budget_skip_planning_without_losing_work():
    run=run_state(max_issued_chunks=1);a=run.for_partition('a')
    assert a(0,64,16)==a(4,64,16)==4
    run.output_written(event())
    snap=run.snapshot()
    assert not snap['prefix_complete'] and snap['stop_reason']=='issued_prefix_budget'
    assert snap['partitions'][0]['ranges']==[[0,4]] and snap['partitions'][0]['cursor']==8
    assert not run.planning_step(lambda state:pytest.fail('truncated'), **FORECAST)['evaluated']
    run=run_state()
    for row in run.snapshot()['partitions']:
        control=run.for_partition(row['id']);cursor,stop=row['variant_range']
        while cursor<stop:cursor+=control(cursor,stop,16)
    run.output_written(event())
    assert run.planning_step(lambda state:pytest.fail('nothing unissued'), **FORECAST)['reason']=='source_fully_issued'
    final=run.finish(successful=True)
    assert final['finished']['successful'] and final['cache']['closed']


def test_deadline_stops_tuning_and_repeated_writes_do_not_extend_it():
    run=run_state(budget=IncrementalPlanningBudget(max_window_seconds=.1))
    run.output_written(event());first=run.snapshot()['first_written']
    with patch('torchgwas.productive_run.time.perf_counter', return_value=first+1):
        run.output_written(event())
        assert not run.planning_step(lambda state:pytest.fail('late'), **FORECAST)['evaluated']
    assert run.snapshot()['first_written']==first
    assert run.for_partition('a')(0,64,16)==4


def test_partition_binding_and_contiguous_coverage_are_enforced():
    run=run_state();a=run.for_partition('a')
    for args in [('cuda:0',(0,64),4),('cuda:0',(64,129),3)]:
        with pytest.raises(ValueError):run.for_scan(*args)
    with pytest.raises(ValueError):a(1,64,16)
    with pytest.raises(ValueError):a(0,64,8)
    with pytest.raises(ValueError):run.finish(successful=True)
    run.finish(successful=False)
    with pytest.raises(ValueError):a(0,64,16)


def test_late_planner_result_cannot_change_the_size():
    clock=[10.]
    with patch('torchgwas.productive_run.time.perf_counter', side_effect=lambda:clock[0]):
        run=run_state(budget=IncrementalPlanningBudget(max_window_seconds=.1))
        run.output_written(event())
        def build(state):
            clock[0]+=.2
            return proposal()
        result=run.planning_step(build,**FORECAST)
        assert result['evaluated'] and not result['usable_for_decision'] and not result['applied']
        assert run.snapshot()['current_chunk_size']==4
        assert run.snapshot()['stop_reason']=='planning_budget_or_horizon'


def test_invalid_multigpu_partition_is_rejected_before_shared_preprocessing():
    import torch
    from torchgwas.linear import linear_scan_multigpu
    from test_adaptive_chunks import Direct
    source=Direct(np.zeros((129,129),np.int8));phenotype=np.zeros((129,3),np.float32)
    run=run_state()
    # Driver partitions on the allocated capacity, not our 64/65 boundary.
    # Force a different configured interval to test fail-before-work behavior.
    run._partitions['a']['variant_range']=[0,60]
    with patch('torchgwas.linear.residualize_and_standardize',side_effect=AssertionError('preprocessed')):
        with pytest.raises(ValueError,match='partition'):
            linear_scan_multigpu(source,phenotype,devices=['cuda:0','cuda:1'],chunk_size=16,
                reader_workers=2,compute_dtype='float32',_chunk_size_selector=run)


@pytest.mark.parametrize('kind', ['jagwas','significant'])
@pytest.mark.parametrize('fsync', [True,False])
def test_writer_reports_empty_and_material_chunks_after_requested_fsync(tmp_path,kind,fsync):
    completed=[];syncs=[]
    if kind=='jagwas':
        chunks=[(0,2,None,np.array([np.nan,np.nan]),None), (2,4,None,np.array([1.,2.]),None)]
    else:
        empty=np.empty(0)
        chunks=[(0,2,empty,empty,empty,empty,123),
                (2,4,np.array([2,3]),np.array([0,1]),np.ones(2),np.ones(2),123)]
    def observe(row):
        assert row.completed>=row.started
        if row.rows:
            assert row.part_bytes==(tmp_path/row.part_file).stat().st_size
            if fsync:assert syncs
        else:assert row.part_file is None and not row.part_file_fsynced
        completed.append(row)
    with patch('torchgwas.sumstats_indexed.os.fsync',side_effect=lambda fd:syncs.append(fd)):
        rows,_=write_indexed_sumstats(tmp_path,['a','b','c','d'],['x','y'],129,iter(chunks),
            kind=kind,df=123,chi2_df=2,fsync=fsync,on_chunk_written=observe)
    assert rows==2 and [x.rows for x in completed]==[0,2]
    assert completed[1].part_file_fsynced is fsync
    with pytest.raises(FrozenInstanceError):completed[1].rows=10


def test_failed_write_never_reports_completion_and_closes_iterator(tmp_path):
    closed=[];completed=[]
    def chunks():
        try:yield (0,1,None,np.ones(1),None)
        finally:closed.append(True)
    with patch('torchgwas.sumstats_indexed.os.fsync',side_effect=OSError('disk failed')):
        with pytest.raises(OSError):
            write_indexed_sumstats(tmp_path,['a'],['x'],129,chunks(),kind='jagwas',
                df=123,chi2_df=1,on_chunk_written=completed.append)
    assert closed==[True] and completed==[] and not (tmp_path/'manifest.json').exists()


def test_tiled_partition_identity_needs_explicit_binding():
    specs=[dict(id=str(i),device='cuda:0',variant_range=[0,10],trait_range=[i*3,(i+1)*3]) for i in range(2)]
    run=ProductiveTuningRun(specs,chunk_sizes=[2,4],initial=2)
    with pytest.raises(ValueError):run.for_scan('cuda:0',(0,10),3)
    assert run.for_partition('0')(0,10,4)==2
    assert run.for_partition('1')(0,10,4)==2
