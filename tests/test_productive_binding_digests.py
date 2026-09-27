"""Productive digest reuse must keep source, price and live checks effective."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import torch

from torchgwas import binding_digests,detailed_calibration as calibration
from torchgwas.analytical_plan_cache import input_identity
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.initial_chunk_autotune import PublicInitialChunkTuning,validate_initial_chunk_config
from test_initial_chunk_autotune import config


@pytest.fixture
def environment(tmp_path,monkeypatch):
    source=tmp_path/'source.py';library=tmp_path/'library.so';data=tmp_path/'input.pgen'
    for path,body in ((source,'source'),(library,'binary'),(data,'input')):path.write_text(body)
    now=[100.];monkeypatch.setattr('time.time',lambda:now[0])
    monkeypatch.setattr(binding_digests,'_stable',lambda row:True)
    monkeypatch.setattr(binding_digests,'_host_identity',lambda:'test-boot')
    live=dict(threads=2,gpu=1000,host=1000)
    def sources(*,digest_cache=None):
        values=([calibration.sha256_file(source)] if digest_cache is None else digest_cache.digests([source]))
        return {'source.py':values[0]}
    def context(*a,digest_cache=None,**kw):
        values=([calibration.sha256_file(library)] if digest_cache is None else digest_cache.digests([library]))
        return dict(devices={'cuda:0':dict(uuid='fixture')},library=values[0],threads=live['threads'])
    monkeypatch.setattr(calibration,'source_identity',sources)
    monkeypatch.setattr(calibration,'execution_context',context)
    monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(live['gpu'],2000))
    monkeypatch.setattr(torch.cuda,'memory_reserved',lambda d:0)
    monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:live['host'])
    deps=dict(source_sha256=sources(),execution_context=context(),measurement_protocol={'operation':'test'})
    record=CalibrationParameterCache(tmp_path/'prices').store('cpu_capacity','test',dict(rate=3.),
        dependencies=deps,provenance={'scope':'Synthetic test'},observed_unix_seconds=90.,max_age_seconds=30.)
    profile=calibration.bind_detailed_profile(
        [dict(name='fixture',devices=['cuda:0'],profiles={'cuda:0':dict(rate=3.)})],context(),sources=sources(),
        component_artifacts={record['path']:calibration.sha256_file(record['path'])},limitations=['Synthetic test'],
        price_bindings=[dict(artifact=record['path'],kind='cpu_capacity',name='test',dependencies=deps,max_age_seconds=None,
            targets=[dict(context_path=[0,'profiles','cuda:0','rate'],value_path=['rate'])])])
    def make(enabled=True):
        cfg=config();cfg.update(binding_digest_cache_dir=str(tmp_path/'digests'),reuse_binding_digests=enabled)
        owner=SimpleNamespace(config=dict(initial_chunks=cfg),profile=deepcopy(profile),reduction=None,
            devices=['cuda:0'],input_path=str(data),output_path=str(tmp_path))
        start=dict(partitions=[dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])],
            chunk_sizes=[4,8],initial_size=4,input_file_identity=input_identity(data),
            memory=dict(device_bytes={'cuda:0':100},host_bytes=100))
        return PublicInitialChunkTuning(owner,start,{},dict.fromkeys(owner.devices,0))
    return SimpleNamespace(make=make,now=now,live=live,source=source,library=library,data=data,
        directory=tmp_path/'digests',record=record,profile=profile)


def complete(tuner):
    bound=tuner.for_partition('cuda:0',[0,257],[0,3]);cursor=0
    while cursor<257:cursor+=bound(cursor,257,8)
    return tuner.finish(successful=True)


def test_digest_reuse_is_lazy_cross_job_and_preserves_measurement_age(environment,monkeypatch):
    e=environment;tuner=e.make();original=Path(e.record['path']).read_bytes()
    assert tuner._binding_digests is None and not e.directory.exists()
    tuner._check();tuner._check()
    assert not e.directory.exists()
    first=complete(tuner);files={p:p.read_bytes() for p in e.directory.rglob('*.json')}
    assert first['binding_digests']['state']['hashed_files']==2
    assert first['binding_digests']['state']['memory_hits']==2
    assert len(first['binding_digests']['publication']['stored'])==2
    assert first['binding_digests']['state']['closed'] and first['validation']['checks']==2
    later=e.make();e.now[0]=110.;real=calibration.sha256_file
    def hash_uncached_only(path):
        assert Path(path) not in (e.source,e.library)
        return real(path)
    monkeypatch.setattr(calibration,'sha256_file',hash_uncached_only)
    later._check();second=complete(later)
    assert second['binding_digests']['state']['disk_hits']==2
    assert second['binding_digests']['state']['hashed_files']==0
    assert second['binding_digests']['publication']['stored']==[]
    assert files=={p:p.read_bytes() for p in e.directory.rglob('*.json')}
    assert Path(e.record['path']).read_bytes()==original and later.owner.profile==e.profile
    assert second['binding_digests']['publication_wall_seconds']>=0


@pytest.mark.parametrize('fault',['source','library','input','threads','gpu','host','expiry','artifact'])
def test_cached_digests_never_bypass_live_price_or_source_validation(environment,fault):
    e=environment;tuner=e.make();tuner._check()
    if fault=='source':e.source.write_text('change')
    if fault=='library':e.library.write_text('change')
    if fault=='input':e.data.write_text('other')
    if fault=='threads':e.live['threads']=4
    if fault=='gpu':e.live['gpu']=1
    if fault=='host':e.live['host']=1
    if fault=='expiry':e.now[0]=121.
    if fault=='artifact':Path(e.record['path']).write_text('changed')
    with pytest.raises(ValueError):tuner._check()
    assert tuner.run.snapshot()['current_chunk_size']==4
    report=tuner.finish(successful=False)
    assert report['binding_digests']['publication']['stored']==[]
    assert report['binding_digests']['state']['closed'] and not e.directory.exists()


def test_digest_publication_rechecks_files_after_the_last_productive_check(environment):
    e=environment;tuner=e.make();tuner._check();e.library.write_text('change')
    report=complete(tuner)
    assert report['finished']['successful']
    assert report['binding_digests']['publication']['status']=='files_changed'
    assert not e.directory.exists() and report['binding_digests']['state']['closed']


def test_optional_cache_failure_falls_back_to_hashing_or_preserves_completed_output(environment,monkeypatch):
    e=environment;tuner=e.make()
    def fail(*a,**kw):raise OSError('unavailable')
    monkeypatch.setattr(binding_digests,'_host_identity',fail)
    tuner._check();report=complete(tuner)
    assert report['finished']['successful'] and not report['binding_digests']['enabled']
    assert 'unavailable' in report['binding_digests']['initialization_error']
    monkeypatch.setattr(binding_digests,'_host_identity',lambda:'test-boot')
    later=e.make();later._check();monkeypatch.setattr(later._binding_digests,'publish',fail)
    report=complete(later)
    assert report['finished']['successful'] and report['binding_digests']['state']['closed']
    assert report['binding_digests']['publication']['status']=='cache_error'


def test_explicit_disable_keeps_full_byte_validation(environment,monkeypatch):
    e=environment;tuner=e.make(enabled=False);hashed=[];real=calibration.sha256_file
    def track(path):hashed.append(Path(path));return real(path)
    monkeypatch.setattr(calibration,'sha256_file',track)
    tuner._check();tuner._check();report=complete(tuner)
    assert hashed.count(e.source)==2 and hashed.count(e.library)==2
    assert not report['binding_digests']['enabled'] and not e.directory.exists()


@pytest.mark.parametrize('field,value',[('reuse_binding_digests',1),('binding_digest_cache_dir',''),('binding_digest_cache_dir',None)])
def test_invalid_digest_controls_fail_before_work(field,value):
    cfg=config();cfg[field]=value
    with pytest.raises(ValueError):
        validate_initial_chunk_config(cfg,dict(chunks=[4,8],trait_blocks=[3]),[dict(name='fixture')],None,
            dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')}))
