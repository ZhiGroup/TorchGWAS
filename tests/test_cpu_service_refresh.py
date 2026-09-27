"""Fresh checks must never relabel old evidence as a new measurement."""
from copy import deepcopy
from pathlib import Path
import pytest

from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.cpu_service_refresh import CpuServiceRefresh


DEPS = dict(source_sha256={'cpu.py':'a'}, execution_context={'host':'h'},
            measurement_protocol={'operation':'controlled_cpu_work','units':100,'loops':4})


@pytest.fixture
def clock(monkeypatch):
    now = [100.]
    monkeypatch.setattr('torchgwas.cpu_service_refresh.time.time', lambda: now[0])
    return now


def controller(tmp_path, **kw):
    options = dict(dependencies=DEPS, work_units=100, max_age_seconds=100.)
    options.update(kw)
    return CpuServiceRefresh(tmp_path/'cache', 'work.v1', **options)


def sample(clock, cpu=.01, elapsed=.1):
    observed = clock[0]; clock[0] += elapsed
    return dict(cpu_seconds=cpu, wall_seconds=elapsed, work_units=100, observed_unix_seconds=observed)


def complete(value, clock, cpu=.01):
    for _ in range(8):
        row = value.advance(lambda: sample(clock, cpu))
        if row['state'] not in ('checking','measuring'):
            return row
    pytest.fail('Unbounded refresh')


def records(tmp_path):
    return {str(p):p.read_bytes() for p in (tmp_path/'cache').glob('*/*.json')}


def test_lazy_bounded_batches_and_original_observation_time(tmp_path, clock, monkeypatch):
    calls = []
    original = CalibrationParameterCache.lookup
    def lookup(*args, **kw):
        calls.append(True); return original(*args, **kw)
    monkeypatch.setattr(CalibrationParameterCache, 'lookup', lookup)
    value = controller(tmp_path)
    assert not calls and not records(tmp_path)
    states = []
    for index in range(4):
        states.append(value.advance(lambda: sample(clock)))
        assert len(records(tmp_path)) == (1 if index==3 else 0)
        clock[0] += 1
    assert [len(row['samples']) for row in states] == [2,4,6,7]
    result = states[-1]['result']
    assert result['record']['observed_unix_seconds']==100.
    assert result['record']['created_unix_seconds']>103.
    assert result['record']['value']['cpu_seconds_per_unit']==.0001
    assert len(calls)==1
    assert value.advance(lambda: pytest.fail('Repeated terminal measurement')) == states[-1]
    states[-1]['result']['record']['value']['samples'].clear()
    assert len(value.snapshot()['result']['record']['value']['samples'])==7


def test_matching_check_reuses_bytes_and_age_without_publication(tmp_path, clock):
    first = complete(controller(tmp_path), clock); before = records(tmp_path)
    clock[0] = 120.
    later = controller(tmp_path).advance(lambda: sample(clock, .011))
    assert later['state']=='ready' and len(later['samples'])==2
    assert later['check']['status']=='consistent'
    assert later['result']['status']=='reused_original'
    assert later['result']['record_sha256']==first['result']['record_sha256']
    assert later['result']['record']['observed_unix_seconds']==100.
    assert later['result']['age_seconds']==pytest.approx(20.2)
    assert records(tmp_path)==before


@pytest.mark.parametrize('rate', [.04,.0025])
def test_drift_refreshes_across_later_chunks_and_keeps_check_samples(tmp_path, clock, rate):
    old = complete(controller(tmp_path), clock); before=records(tmp_path); clock[0]=120.
    value = controller(tmp_path); first=value.advance(lambda: sample(clock,rate))
    assert first['state']=='measuring' and first['check']['status']=='drift'
    assert first['result'] is None and records(tmp_path)==before
    result = complete(value,clock,rate)
    assert result['steps']==4 and len(result['samples'])==7
    assert result['samples'][:2]==first['samples']
    assert result['result']['record']['observed_unix_seconds']==120.
    assert result['result']['record_sha256']!=old['result']['record_sha256']
    assert result['result']['record']['value']['cpu_seconds_per_unit']==rate/100
    assert len(records(tmp_path))==2
    assert all(records(tmp_path)[path]==body for path,body in before.items())


def test_record_expiring_during_check_is_not_reused(tmp_path, clock):
    complete(controller(tmp_path,max_age_seconds=10.),clock)
    clock[0]=109.
    value=controller(tmp_path,max_age_seconds=10.)
    check=value.advance(lambda:sample(clock,elapsed=.6))
    assert check['state']=='measuring'
    assert check['check']['status']=='expired_or_invalid_during_check'
    fresh=complete(value,clock)
    assert fresh['result']['record']['observed_unix_seconds']==109.
    assert fresh['result']['status']=='measured_and_published'


def test_newer_concurrent_record_cannot_replace_checked_record(tmp_path, clock):
    first=complete(controller(tmp_path),clock)
    value=controller(tmp_path); count=[]
    def probe():
        if not count:
            count.append(True)
            other=controller(tmp_path,max_age_seconds=.01)
            # A separately completed compatible record with a different value
            # cannot substitute for the specific baseline already being checked.
            fresh=value._value([sample(clock,.02) for _ in range(7)])
            other.cache.store('cpu_capacity','work.v1',fresh,dependencies=DEPS,
                provenance={'job':'other'},max_age_seconds=100.,
                observed_unix_seconds=fresh['samples'][0]['observed_unix_seconds'])
        return sample(clock)
    result=value.advance(probe)
    assert result['result']['record_sha256']==first['result']['record_sha256']
    assert len(records(tmp_path))==2


@pytest.mark.parametrize('change', ['source','protocol','shorter_age'])
def test_changed_dependencies_or_shorter_lifetime_trigger_fresh_measurement(tmp_path,clock,change):
    first=complete(controller(tmp_path),clock);clock[0]=110.
    deps=deepcopy(DEPS);age=100.
    if change=='source':deps['source_sha256']['cpu.py']='b'
    if change=='protocol':deps['measurement_protocol']['loops']=8
    if change=='shorter_age':age=5.
    result=complete(controller(tmp_path,dependencies=deps,max_age_seconds=age),clock)
    assert result['result']['status']=='measured_and_published'
    assert result['result']['record_sha256']!=first['result']['record_sha256']
    assert result['check'] is None


def test_publication_cannot_renew_expired_new_samples(tmp_path,clock):
    value=controller(tmp_path,max_age_seconds=1.)
    value.advance(lambda:sample(clock));clock[0]+=2.
    result=complete(value,clock)
    assert result['state']=='expired_measurement' and result['result'] is None
    assert not records(tmp_path)


def test_publication_latency_does_not_make_expired_result_usable(tmp_path,clock,monkeypatch):
    value=controller(tmp_path,max_age_seconds=1.); original=value.cache.store
    def delayed(*args,**kw):
        result=original(*args,**kw);clock[0]+=2.;return result
    monkeypatch.setattr(value.cache,'store',delayed)
    result=complete(value,clock)
    assert result['state']=='expired_measurement' and result['result'] is None
    assert len(records(tmp_path))==1  # Complete evidence, unusable at return.


def test_unstable_complete_window_is_not_published(tmp_path,clock):
    value=controller(tmp_path);rates=iter([.001,.1,.001,.1,.001,.1,.001])
    for _ in range(4):result=value.advance(lambda:sample(clock,next(rates)))
    assert result['state']=='unstable_measurement' and result['result'] is None
    assert not records(tmp_path)


@pytest.mark.parametrize('reverse',[False,True])
def test_drift_within_permitted_spread_cannot_publish_a_median(tmp_path,clock,reverse):
    values=[.01,.01,.01,.02,.03,.03,.03]
    if reverse:values.reverse()
    rates=iter(values);value=controller(tmp_path)
    for _ in range(4):result=value.advance(lambda:sample(clock,next(rates)))
    stability=result['window_stability']
    assert stability['spread']['stable'] and stability['spread']['ratio']==3.
    assert stability['temporal_drift'] and not stability['stable']
    assert result['state']=='unstable_measurement' and result['result'] is None
    assert not records(tmp_path)


def test_inconsistent_short_check_cannot_hide_behind_its_matching_median(tmp_path,clock):
    first=complete(controller(tmp_path),clock);before=records(tmp_path);clock[0]=120.
    value=controller(tmp_path);rates=iter([.005,.015])
    checked=value.advance(lambda:sample(clock,next(rates)))
    assert checked['check']['ratio']==pytest.approx(1.) and not checked['check']['drift']
    assert checked['check']['spread']['stable']
    assert checked['check']['window_stability']['temporal_drift']
    assert checked['check']['status']=='unstable_check'
    assert checked['state']=='measuring' and checked['result'] is None
    assert records(tmp_path)==before
    # A later stable window can retain those samples and publish new evidence.
    refreshed=complete(value,clock,.015)
    assert refreshed['result']['status']=='measured_and_published'
    assert refreshed['result']['record_sha256']!=first['result']['record_sha256']
    assert refreshed['samples'][:2]==checked['samples']
    assert refreshed['result']['record']['provenance']['window_stability']['stable']
    assert all(records(tmp_path)[path]==body for path,body in before.items())


def test_previously_saved_drifting_window_is_rejected_without_rewriting_it(tmp_path,clock):
    value=controller(tmp_path)
    rows=[sample(clock,rate) for rate in [.01,.01,.01,.02,.03,.03,.03]]
    value.cache.store('cpu_capacity','work.v1',value._value(rows),dependencies=DEPS,
        provenance={'scope':'Prior producer allowed spread without a temporal check'},
        max_age_seconds=100.,observed_unix_seconds=rows[0]['observed_unix_seconds'])
    before=records(tmp_path)
    result=complete(value,clock)
    assert result['lookup_reason']=='invalid_or_expired_cpu_window'
    assert result['result']['status']=='measured_and_published'
    assert all(records(tmp_path)[path]==body for path,body in before.items())


@pytest.mark.parametrize('values',[
    [.01,.02,.01,.03,.02,.01,.02],
    [1e-7,1e-7,1e-7,2e-7,3e-7,3e-7,3e-7]])
def test_window_check_keeps_declared_ratio_and_absolute_tolerances(tmp_path,clock,values):
    rates=iter(values);value=controller(tmp_path)
    for _ in range(4):result=value.advance(lambda:sample(clock,next(rates)))
    assert result['state']=='ready' and result['window_stability']['stable']
    assert not result['window_stability']['temporal_drift']


@pytest.mark.parametrize('fault', ['missing','nan','zero','wrong_work','old_timestamp','future_timestamp','probe_error'])
def test_failed_or_malformed_partial_window_cannot_be_published(tmp_path,clock,fault):
    value=controller(tmp_path);value.advance(lambda:sample(clock))
    def bad():
        row=sample(clock)
        if fault=='missing':row.pop('wall_seconds')
        if fault=='nan':row['cpu_seconds']=float('nan')
        if fault=='zero':row['cpu_seconds']=0
        if fault=='wrong_work':row['work_units']=200
        if fault=='old_timestamp':row['observed_unix_seconds']-=1
        if fault=='future_timestamp':row['observed_unix_seconds']+=1
        if fault=='probe_error':raise RuntimeError('probe failed')
        return row
    with pytest.raises((ValueError,RuntimeError)):value.advance(bad)
    assert value.snapshot()['state']=='failed' and not records(tmp_path)
    value.advance(lambda:pytest.fail('Failed sampler was retried'))


def test_invalid_cached_window_triggers_measurement_not_price_use(tmp_path,clock):
    value=controller(tmp_path)
    value.cache.store('cpu_capacity','work.v1',{'samples':[]},dependencies=DEPS,
        provenance={'job':'incomplete'},max_age_seconds=100.,observed_unix_seconds=clock[0])
    result=complete(value,clock)
    assert result['lookup_reason']=='invalid_or_expired_cpu_window'
    assert result['result']['status']=='measured_and_published'


@pytest.mark.parametrize('kw', [dict(work_units=0),dict(max_age_seconds=True),dict(samples_per_step=129),
    dict(check_samples=1),dict(refresh_samples=2),dict(drift_ratio=1),dict(max_sample_ratio=float('inf')),
    dict(dependencies={'source':'a'})])
def test_invalid_policy_is_rejected(tmp_path,clock,kw):
    with pytest.raises(ValueError):controller(tmp_path,**kw)


def test_productive_budget_charges_each_batch_and_stops_without_losing_source(tmp_path,clock,monkeypatch):
    from test_planning_session import Clock
    from test_productive_run import run_state,event,finish_source
    from torchgwas.planning_session import IncrementalPlanningBudget
    timing=Clock()
    monkeypatch.setattr('torchgwas.planning_session.time.perf_counter',lambda:timing.wall)
    monkeypatch.setattr('torchgwas.planning_session.time.thread_time',lambda:timing.cpu)
    run=run_state(budget=IncrementalPlanningBudget(max_cpu_seconds=.08,max_steps=5))
    probe=controller(tmp_path)
    def measure(state):
        def one():
            timing.advance(.01,.015);return sample(clock)
        probe.advance(one)
        return dict(chunk_size=state['current_chunk_size'],baseline_seconds=0.,candidate_seconds=0.)
    costs=dict(remaining_seconds=10.,expected_cpu_seconds=.001,expected_wall_seconds=.002)
    first=run.planning_step(measure,**costs)
    assert not first['evaluated'] and probe.snapshot()['state']=='pending'
    outcomes=[]
    for _ in range(3):
        run.output_written(event());outcomes.append(run.planning_step(measure,**costs))
    assert not outcomes[-1]['usable_for_decision'] and not outcomes[-1]['applied']
    assert outcomes[-1]['total_tuning_cost_seconds']==pytest.approx(.06)
    assert run.snapshot()['planning']['cpu_seconds']==pytest.approx(.09)
    assert len(probe.snapshot()['samples'])==6 and not records(tmp_path)
    assert run.for_partition('a')(0,64,16)==4
    finish_source(run)
    assert run.finish(successful=True)['finished']['successful']
