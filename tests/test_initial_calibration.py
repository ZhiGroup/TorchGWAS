"""Cache hits receive live checks; only changed GPUs spend the refresh budget."""
from dataclasses import asdict, replace
from pathlib import Path
import time
from unittest.mock import patch

import pytest

from torchgwas.adaptive_chunks import ChunkObservation, ChunkDeviceTiming
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.initial_calibration import (CACHE_NAME, InitialCalibrationController,
                                           compare_component_windows)


DEPS=dict(protocol='synthetic-controller-test',source='v1',shape=[129,10,5],input='fixture')


def observation(device, start, *, scale=1., width=10):
    t=time.perf_counter()
    return ChunkObservation(start,start+width,10,device,t,t+.01*scale,t+.012*scale,
        t+.02*scale,t+.03*scale,.01*scale,1,width*40,
        cuda=ChunkDeviceTiming(.002*scale,.001*scale,.005*scale,.0003*scale),
        consumer_cpu_seconds=.003*scale,consumer_probe_wall_seconds=.0002*scale)


def rows(device, *, scale=1.):
    offset=0 if device=='cuda:0' else 1000
    return [asdict(observation(device,offset+i*10,scale=scale)) for i in range(8)]


def seed(cache, device, *, scale=1., age=10.):
    return cache.store('stage_observations',CACHE_NAME+':'+device,
        dict(device=device,observations=rows(device,scale=scale),
             consumer_probe_protocol='thread_cpu_schedstat_probe_wall_v2'),dependencies=DEPS,
        provenance=dict(source='synthetic control'),max_age_seconds=300.,
        observed_unix_seconds=time.time()-age)


def controller(cache, devices=('cuda:0',), **kwargs):
    defaults=dict(max_chunks_per_device=8,validation_chunks_per_device=2,warmup_chunks=0,
                  stride=1,max_window_seconds=100.,max_age_seconds=300.)
    defaults.update(kwargs)
    return InitialCalibrationController(devices,cache=cache,dependencies=DEPS,
        provenance=dict(job='unit-test'),**defaults)


def feed(control, device, start, scale=1.):
    if control.reserve_read(start,start+10,device):
        control(observation(device,start,scale=scale))
        return True
    return False


def test_matched_work_check_detects_changed_spans_without_inferring_capacity():
    old=rows('cuda:0')[:2]
    assert compare_component_windows(old,rows('cuda:0')[:2])['status']=='consistent'
    changed=compare_component_windows(old,rows('cuda:0',scale=4.)[:2])
    assert changed['status']=='drift'
    assert all(row['drift'] for row in changed['metrics'].values())
    assert changed['matched_samples']==2
    changed_shape=[asdict(observation('cuda:0',start,width=7)) for start in (0,10)]
    assert compare_component_windows(old,changed_shape)['status']=='incomparable'
    changed_payload=rows('cuda:0')[:2];changed_payload[0]['result_bytes']+=4
    assert compare_component_windows(old,changed_payload)['changed_output_ranges']
    unknown=rows('cuda:0')[:2];unknown[0]['cuda']=None
    assert compare_component_windows(old,unknown)['status']=='incomparable'


def test_tiny_intervals_do_not_trigger_arbitrary_ratios():
    old=rows('cuda:0')[:2];fresh=rows('cuda:0')[:2]
    for a,b in zip(old,fresh):
        a['cuda']['h2d']=0.;b['cuda']['h2d']=1e-7
    report=compare_component_windows(old,fresh)
    assert report['status']=='consistent' and report['metrics']['h2d']['ratio'] is None


def test_reader_cpu_drift_refreshes_a_matched_cached_window():
    old=rows('cuda:0')[:2];fresh=rows('cuda:0')[:2]
    for row in old: row['read_cpu_seconds']=.004
    for row in fresh: row['read_cpu_seconds']=.020
    report=compare_component_windows(old,fresh)
    assert report['status']=='drift'
    assert report['metrics']['read_thread_cpu']['drift']
    assert not report['metrics']['read_decode']['drift']
    for row in old: row['read_cpu_seconds']=None
    assert compare_component_windows(old,fresh)['incomparable_metric']=='read_thread_cpu'
    fresh[0]['read_cpu_seconds']=-1.
    with pytest.raises(ValueError,match='reader thread CPU'):
        compare_component_windows(old,fresh)


def test_downstream_cpu_drift_refreshes_without_calling_it_writer_service():
    old=rows('cuda:0')[:2];fresh=rows('cuda:0')[:2]
    for row in fresh:row['consumer_cpu_seconds']=.03
    report=compare_component_windows(old,fresh)
    assert report['status']=='drift'
    assert report['metrics']['consumer_thread_cpu']['drift']
    assert not report['metrics']['read_decode']['drift']
    del old[0]['consumer_cpu_seconds']
    assert compare_component_windows(old,fresh)['incomparable_metric']=='consumer_thread_cpu'
    fresh[0]['consumer_cpu_seconds']=-1.
    with pytest.raises(ValueError,match='consumer thread CPU'):
        compare_component_windows(old,fresh)
    fresh[0]['consumer_cpu_seconds']=.003
    fresh[0]['consumer_probe_wall_seconds']=-1.
    with pytest.raises(ValueError,match='consumer_probe_wall_seconds'):
        compare_component_windows(old,fresh)


def test_per_device_refresh_does_not_renew_unchanged_device_record(tmp_path):
    cache=CalibrationParameterCache(tmp_path/'cache')
    initial={d:seed(cache,d) for d in ('cuda:0','cuda:1')}
    original={d:Path(row['path']).read_bytes() for d,row in initial.items()}
    control=controller(cache,('cuda:0','cuda:1'))
    for i in range(8):
        feed(control,'cuda:0',10*i)
        feed(control,'cuda:1',1000+10*i,scale=4.)
    snapshot=control.snapshot()
    assert snapshot['reserved']=={'cuda:0':2,'cuda:1':8}
    assert snapshot['current_limits']=={'cuda:0':2,'cuda:1':8}
    assert snapshot['refresh_decisions']['cuda:0']['state']=='cached_consistent'
    assert snapshot['refresh_decisions']['cuda:1']['state']=='refreshed'
    report=control.finish(successful=True)
    assert not report['incomplete_devices'] and len(report['publications'])==2
    assert control.finish(successful=True)==report
    for device in ('cuda:0','cuda:1'):
        assert Path(initial[device]['path']).read_bytes()==original[device]
        current=cache.lookup('stage_observations',CACHE_NAME+':'+device,dependencies=DEPS)
        assert current['hit']
        assert (current['record_sha256']==initial[device]['record_sha256'])==(device=='cuda:0')
    assert cache.lookup('stage_observations',CACHE_NAME+':cuda:0:validation',dependencies=DEPS)['hit']


@pytest.mark.parametrize('reason',['missing','expired','invalid'])
def test_unusable_prior_collects_full_bounded_window(tmp_path, reason):
    cache=CalibrationParameterCache(tmp_path/'cache')
    if reason=='expired':seed(cache,'cuda:0',age=301.)
    if reason=='invalid':
        cache.store('stage_observations',CACHE_NAME+':cuda:0',dict(device='cuda:0',observations=[]),
            dependencies=DEPS,provenance=dict(job='malformed'),max_age_seconds=300.)
    control=controller(cache,max_chunks_per_device=3)
    for i in range(20):feed(control,'cuda:0',i*10)
    snapshot=control.snapshot()
    assert snapshot['reserved']=={'cuda:0':3}
    assert snapshot['refresh_decisions']['cuda:0']['state']=='collected'
    assert not snapshot['refresh_decisions']['cuda:0']['cache_hit']
    assert not control.finish(successful=True)['incomplete_devices']
    assert cache.lookup('stage_observations',CACHE_NAME+':cuda:0',dependencies=DEPS)['hit']


@pytest.mark.parametrize('version',['v1','v2'])
def test_old_probe_protocol_cannot_validate_new_component_window(tmp_path,version):
    cache=CalibrationParameterCache(tmp_path/'cache')
    old=cache.store('stage_observations','initial_chunk_components.'+version+':cuda:0',
        dict(device='cuda:0',observations=rows('cuda:0')),
        dependencies=DEPS,provenance=dict(job='old-protocol'),max_age_seconds=300.)
    old_bytes=Path(old['path']).read_bytes()
    control=controller(cache,max_chunks_per_device=2)
    decision=control.snapshot()['refresh_decisions']['cuda:0']
    assert not decision['cache_hit'] and decision['state']=='collecting'
    feed(control,'cuda:0',0);feed(control,'cuda:0',10)
    report=control.finish(successful=True)
    assert report['publications']['cuda:0']['record_sha256']!=old['record_sha256']
    assert Path(old['path']).read_bytes()==old_bytes


def test_missing_consumer_probe_binding_forces_refresh(tmp_path):
    cache=CalibrationParameterCache(tmp_path/'cache')
    cache.store('stage_observations',CACHE_NAME+':cuda:0',
        dict(device='cuda:0',observations=rows('cuda:0')),
        dependencies=DEPS,provenance=dict(job='missing-protocol'),max_age_seconds=300.)
    control=controller(cache,max_chunks_per_device=2)
    decision=control.snapshot()['refresh_decisions']['cuda:0']
    assert not decision['cache_hit'] and decision['cache_reason']=='invalid_component_window'


def test_short_or_failed_job_does_not_publish_a_fresh_baseline(tmp_path):
    cache=CalibrationParameterCache(tmp_path/'cache')
    prior=seed(cache,'cuda:0')
    short=controller(cache)
    feed(short,'cuda:0',0)
    report=short.finish(successful=True)
    assert report['incomplete_devices']==['cuda:0'] and not report['publications']
    assert cache.lookup('stage_observations',CACHE_NAME+':cuda:0',dependencies=DEPS)['record_sha256']==prior['record_sha256']
    failed=controller(cache)
    feed(failed,'cuda:0',0);feed(failed,'cuda:0',10)
    assert not failed.finish(successful=False)['publications']


def test_drift_cannot_expand_beyond_window_deadline(tmp_path):
    cache=CalibrationParameterCache(tmp_path/'cache');seed(cache,'cuda:0')
    control=controller(cache,max_window_seconds=1.)
    with patch('torchgwas.adaptive_chunks.time.perf_counter',return_value=10.):
        assert feed(control,'cuda:0',0,scale=4.)
        assert feed(control,'cuda:0',10,scale=4.)
    with patch('torchgwas.adaptive_chunks.time.perf_counter',return_value=11.):
        assert not feed(control,'cuda:0',20,scale=4.)
        assert control.snapshot()['new_measurements_stopped']
    # Synthetic counter values have no epoch meaning: no incomplete publication.
    report=control.finish(successful=True)
    assert report['incomplete_devices']==['cuda:0'] and not report['publications']


def test_cache_io_failure_does_not_invalidate_completed_results(tmp_path):
    cache=CalibrationParameterCache(tmp_path/'cache')
    control=controller(cache,max_chunks_per_device=2)
    feed(control,'cuda:0',0);feed(control,'cuda:0',10)
    with patch.object(cache,'store',side_effect=OSError('cache is read only')):
        report=control.finish(successful=True)
    assert report['successful']
    assert report['publications']['cuda:0']==dict(saved=False,error='cache is read only')


def test_decision_cpu_budget_stops_new_reservations_without_canceling_data(tmp_path):
    cache=CalibrationParameterCache(tmp_path/'cache');seed(cache,'cuda:0')
    control=controller(cache,max_decision_cpu_seconds=.01)
    assert control.reserve_read(0,10,'cuda:0')
    assert control.reserve_read(10,20,'cuda:0')
    with patch('torchgwas.initial_calibration.time.thread_time',side_effect=[0.,.02]):
        control(observation('cuda:0',0))
    assert not control.reserve_read(20,30,'cuda:0')
    # Already issued observations/results must still complete after the cutoff.
    control(observation('cuda:0',10))
    snapshot=control.snapshot()
    assert snapshot['decision_budget_exhausted'] and snapshot['decision_cpu_seconds']>=.02
    assert not snapshot['pending'] and len(snapshot['observations'])==2
    assert control.finish(successful=True)['decision_budget_exhausted']


@pytest.mark.parametrize('provenance',[None,{},'job'])
def test_invalid_provenance_fails_before_scan(tmp_path,provenance):
    with pytest.raises(ValueError,match='provenance'):
        InitialCalibrationController(['cuda:0'],cache=CalibrationParameterCache(tmp_path),
            dependencies=DEPS,provenance=provenance)
