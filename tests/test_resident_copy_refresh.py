"""Public refresh never renews old evidence or uses it before checking."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import time

import pytest

from torchgwas.calibration_cache import CalibrationParameterCache,_digest,read_calibration_record
from torchgwas.detailed_calibration import bind_detailed_profile,sha256_file,validate_detailed_profile
from torchgwas.initial_chunk_autotune import PublicInitialChunkTuning
from torchgwas.price_binding import validate_price_bindings
from torchgwas.resident_copy_refresh import ResidentCopyRefresh,ELEMENTS
from torchgwas.sumstats import DenseWriteProgress
from test_initial_chunk_autotune import config as initial_config


@pytest.fixture
def seed(tmp_path,monkeypatch):
    now=[100.];monkeypatch.setattr('time.time',lambda:now[0])
    execution=dict(devices={'cuda:0':dict(uuid='fixture')});source={'test':'fixed'}
    deps=dict(source_sha256=source,execution_context=execution,measurement_protocol={'operation':'old-copy-protocol'})
    record=CalibrationParameterCache(tmp_path/'old').store('cpu_capacity','old-copy',dict(rate=1e-10),
        dependencies=deps,provenance={'scope':'Synthetic test'},max_age_seconds=20.,observed_unix_seconds=100.)
    contexts=[dict(name='fixture',devices=['cuda:0'],profiles={'cuda:0':dict(
        return_beta=True,process_units=dict(numpy_copy_bytes=1e-10),
        owned_result_copy_scenario=dict(resident_cpu_seconds_per_byte=1e-10))})]
    binding=dict(artifact=str(Path(record['path']).resolve()),kind='cpu_capacity',name='old-copy',dependencies=deps,
        max_age_seconds=None,targets=[dict(context_path=[0,'profiles','cuda:0','process_units','numpy_copy_bytes'],value_path=['rate'])])
    profile=bind_detailed_profile(contexts,execution,sources=source,limitations=['Synthetic accounting test'],
        component_artifacts={record['path']:sha256_file(record['path'])},price_bindings=[binding])
    cfg=dict(cache_dir=str(tmp_path/'fresh'),profile_dir=str(tmp_path/'profiles'),binding_indexes=[0],
        max_age_seconds=20.,expected_cpu_seconds=.01,expected_wall_seconds=.01)
    return SimpleNamespace(now=now,profile=profile,config=cfg,source=source,execution=execution,record=record,tmp=tmp_path)


def sample(e,cpu=.01):
    observed=e.now[0];e.now[0]+=.01
    return dict(cpu_seconds=cpu,wall_seconds=.02,work_units=ELEMENTS*4,observed_unix_seconds=observed)


def collect(e,profile=None,cpu=.01):
    refresh=ResidentCopyRefresh(e.profile if profile is None else profile,e.config)
    refresh.sample=lambda:sample(e,cpu)
    for event in range(1,6):
        result=refresh.advance(written_events=event,issued_revision=event+2)
        if result is not None:return refresh
    pytest.fail('Refresh exceeded bounded sample count')


def test_only_named_expiry_can_be_deferred_and_never_validates_a_timing_price(seed):
    e=seed;before=deepcopy(e.profile);e.now[0]=150.
    with pytest.raises(ValueError,match='[Ee]xpired'):validate_price_bindings(e.profile)
    report=validate_price_bindings(e.profile,deferred_bindings=[0])
    assert report['status']=='pending_refresh' and report['verified_targets']==0
    assert report['bindings'][0]['fresh'] is False and report['bindings'][0]['age_seconds']==50.
    assert validate_detailed_profile(e.profile,e.execution,sources=e.source,deferred_bindings=[0])['context_matches']
    assert e.profile==before
    with pytest.raises(ValueError,match='[Ee]xpired'):
        read_calibration_record(e.record['path'],kind='cpu_capacity',name='old-copy',dependencies=e.profile['price_bindings'][0]['dependencies'])
    # Deferral is per binding, even when another target shares the same record.
    second=deepcopy(e.profile);other=deepcopy(second['price_bindings'][0])
    other['targets'][0]['context_path']=[0,'profiles','cuda:0','owned_result_copy_scenario','resident_cpu_seconds_per_byte']
    second['price_bindings'].append(other)
    with pytest.raises(ValueError,match='[Ee]xpired'):validate_price_bindings(second,deferred_bindings=[0])
    assert validate_price_bindings(second,deferred_bindings=[0,1])['status']=='pending_refresh'


@pytest.mark.parametrize('fault',['artifact','value','dependency','future','index','kind'])
def test_deferral_does_not_bypass_integrity_or_context(seed,fault):
    e=seed;p=deepcopy(e.profile);e.now[0]=150.;indexes=[0]
    if fault=='artifact':Path(e.record['path']).write_text('changed')
    if fault=='value':p['contexts'][0]['profiles']['cuda:0']['process_units']['numpy_copy_bytes']=2.
    if fault=='dependency':p['price_bindings'][0]['dependencies']['source_sha256']={}
    if fault=='future':e.now[0]=50.
    if fault=='index':indexes=[1]
    if fault=='kind':p['price_bindings'][0]['kind']='stage_observations'
    with pytest.raises(ValueError):validate_price_bindings(p,deferred_bindings=indexes)


def test_expired_original_refreshes_in_four_batches_without_upfront_io(seed,monkeypatch):
    e=seed;e.now[0]=150.;before=Path(e.record['path']).read_bytes();profile_before=deepcopy(e.profile)
    probe=ResidentCopyRefresh(e.profile,e.config)
    assert not (e.tmp/'fresh').exists() and not (e.tmp/'profiles').exists()
    probe.sample=lambda:sample(e)
    for i in range(4):
        result=probe.advance(written_events=i+1,issued_revision=i+4)
        assert (result is not None)==(i==3)
        assert len(list((e.tmp/'fresh').glob('*/*.json')))==int(i==3)
    report=probe.snapshot();bound=report['published']['price_evidence']['bindings'][0]
    assert [batch['samples'] for batch in report['batches']]==[2,4,6,7]
    assert bound['observed_unix_seconds']==150. and bound['created_unix_seconds']>150.
    assert report['buffers_released'] and not report['pending']
    assert Path(e.record['path']).read_bytes()==before and e.profile==profile_before
    assert str(Path(e.record['path']).resolve()) not in probe.profile['component_artifacts']
    assert _digest(probe.profile)!=_digest(e.profile)
    assert validate_price_bindings(probe.profile)['status']=='declared_targets_verified'


def test_matching_later_check_reuses_exact_measurement_and_profile(seed):
    e=seed;e.now[0]=150.;first=collect(e);old=first.snapshot()
    original=deepcopy(first.profile);files={p:p.read_bytes() for root in ('fresh','profiles') for p in (e.tmp/root).rglob('*.json')}
    e.now[0]=155.;later=collect(e,first.profile,cpu=.011);report=later.snapshot()
    assert len(report['controller']['samples'])==2 and report['controller']['result']['status']=='reused_original'
    assert later.profile==original
    assert report['published']['profile_sha256']==old['published']['profile_sha256']
    assert report['controller']['result']['record_sha256']==old['controller']['result']['record_sha256']
    assert report['controller']['result']['record']['observed_unix_seconds']==150.
    assert files=={p:p.read_bytes() for root in ('fresh','profiles') for p in (e.tmp/root).rglob('*.json')}


@pytest.mark.parametrize('reason',['drift','expiry'])
def test_drift_or_expiry_creates_new_evidence_and_keeps_old_bytes(seed,reason):
    e=seed;e.now[0]=150.;first=collect(e);old=first.snapshot()['controller']['result']
    old_bytes=Path(old['path']).read_bytes();e.now[0]=155. if reason=='drift' else 175.
    later=collect(e,first.profile,cpu=.04 if reason=='drift' else .01)
    report=later.snapshot()['controller'];new=report['result']
    assert new['record_sha256']!=old['record_sha256'] and Path(old['path']).read_bytes()==old_bytes
    assert len(report['samples'])==7 and new['record']['observed_unix_seconds']>old['record']['observed_unix_seconds']
    assert (report['check'] is not None)==(reason=='drift')


def test_partial_or_failed_probe_never_publishes_a_profile(seed):
    e=seed;e.now[0]=150.;probe=ResidentCopyRefresh(e.profile,e.config);probe.sample=lambda:sample(e)
    assert probe.advance(written_events=1,issued_revision=3) is None
    def failed():raise OSError('probe failed')
    probe.sample=failed
    with pytest.raises(OSError,match='probe failed'):probe.advance(written_events=2,issued_revision=4)
    probe.close()
    assert probe.snapshot()['closed'] and probe.snapshot()['buffers_released']
    assert not list((e.tmp/'profiles').glob('*.json')) and not list((e.tmp/'fresh').glob('*/*.json'))


def test_changed_context_after_last_sample_prevents_record_and_profile_publication(seed):
    e=seed;e.now[0]=150.;probe=ResidentCopyRefresh(e.profile,e.config);probe.sample=lambda:sample(e)
    for i in range(3):probe.advance(written_events=i+1,issued_revision=i+3)
    def changed():raise ValueError('context changed during probe')
    with pytest.raises(ValueError,match='context changed'):
        probe.advance(written_events=4,issued_revision=6,validate=changed)
    assert probe.snapshot()['controller']['state']=='failed'
    assert probe.snapshot()['buffers_released']
    assert not list((e.tmp/'profiles').glob('*.json')) and not list((e.tmp/'fresh').glob('*/*.json'))


@pytest.mark.parametrize('field',['path','target','duplicate','age','cost'])
def test_invalid_refresh_contract_is_rejected(seed,field):
    e=seed;profile=deepcopy(e.profile);config=deepcopy(e.config)
    if field=='path':config['cache_dir']=''
    if field=='target':profile['price_bindings'][0]['targets'][0]['context_path']=[0,'profiles','cuda:0','process_units','gpu_bandwidth']
    if field=='duplicate':config['binding_indexes']=[0,0]
    if field=='age':config['max_age_seconds']=float('inf')
    if field=='cost':config['expected_wall_seconds']=0.
    with pytest.raises(ValueError):ResidentCopyRefresh(profile,config)


@pytest.mark.parametrize('stop_early',[False,True])
def test_productive_controller_charges_probe_then_forecast_and_preserves_source(seed,monkeypatch,stop_early):
    e=seed;e.now[0]=150.;refresh=ResidentCopyRefresh(e.profile,e.config);refresh.sample=lambda:sample(e)
    cfg=initial_config('trait',3);cfg['budget']['max_steps']=1 if stop_early else 5
    cfg['resident_copy_refresh']=e.config
    start=dict(partitions=[dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])],
        chunk_sizes=[4,8],initial_size=4,context='fixture',
        candidate=dict(tiles=[dict(device='cuda:0',profile=deepcopy(e.profile['contexts'][0]['profiles']['cuda:0']))]))
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),profile=deepcopy(e.profile),refresh=refresh)
    tuner=PublicInitialChunkTuning(owner,start,{},{});seen=[]
    monkeypatch.setattr(tuner,'_check',lambda:validate_price_bindings(owner.profile,deferred_bindings=refresh.deferred_bindings))
    def propose(snapshot,size):
        assert not refresh.pending
        assert len(refresh.snapshot()['batches'])==4
        assert start['candidate']['tiles'][0]['profile']['process_units']['numpy_copy_bytes']==owner.profile['contexts'][0]['profiles']['cuda:0']['process_units']['numpy_copy_bytes']
        seen.append(snapshot['written_events'])
        return dict(chunk_size=size,baseline_seconds=90.,candidate_seconds=0.)
    monkeypatch.setattr(tuner,'_propose',propose)
    bound=tuner.for_partition('cuda:0',[0,257],[0,3]);cursor=0
    for _ in range(5):
        count=bound(cursor,257,8);cursor+=count
        tuner.output_written(DenseWriteProgress(cursor-count,cursor,(0,3),12*count,cursor,time.perf_counter(),'test'))
    state=tuner.run.snapshot()
    assert state['current_chunk_size']==(4 if stop_early else 8)
    assert seen==([] if stop_early else [5])
    assert len(state['planning']['steps'])==(1 if stop_early else 5)
    assert all(not row['applied'] for row in state['decisions'][:(1 if stop_early else 4)])
    if not stop_early:assert state['decisions'][-1]['total_tuning_cost_seconds']>=state['planning']['wall_seconds']
    while cursor<257:cursor+=bound(cursor,257,8)
    result=tuner.finish(successful=True)
    assert result['resident_copy_refresh']['buffers_released'] and result['resident_copy_refresh']['closed']
    assert result['resident_copy_refresh']['pending']==stop_early
