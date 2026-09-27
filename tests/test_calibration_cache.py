"""A previous job can contribute evidence without freezing volatile capacity."""
import json
from pathlib import Path
from unittest.mock import patch
import pytest
from torchgwas.calibration_cache import CalibrationParameterCache


@pytest.fixture
def cache(tmp_path):
    return CalibrationParameterCache(tmp_path/'calibration')


DEPENDENCIES = dict(source_sha256='source-v1', cuda_library_sha256='library-v1',
                    gpu_uuid='physical-0', shape=[129, 17, 5, 2], protocol='v1')
PROVENANCE = dict(job='test-job', artifact='measured-component.json')


def test_structural_records_reuse_without_age_expiry_but_require_exact_dependencies(cache):
    with patch('torchgwas.calibration_cache.time.time', return_value=10.):
        first = cache.store('kernel_geometry', 'statistics', {'block': [256, 1, 1]},
                            dependencies=DEPENDENCIES, provenance=PROVENANCE)
    with patch('torchgwas.calibration_cache.time.time', return_value=1e10):
        hit = cache.lookup('kernel_geometry', 'statistics', dependencies=DEPENDENCIES)
        assert hit['hit'] and hit['record_sha256'] == first['record_sha256']
        for key, value in [('source_sha256','source-v2'), ('cuda_library_sha256','library-v2'),
                           ('gpu_uuid','physical-1'), ('shape',[129, 7, 5, 2]), ('protocol','v2')]:
            result = cache.lookup('kernel_geometry', 'statistics', dependencies=dict(DEPENDENCIES, **{key:value}))
            assert not result['hit']


def test_empirical_lifetime_cannot_be_extended_by_a_later_job(cache):
    with patch('torchgwas.calibration_cache.time.time', return_value=10.):
        cache.store('storage_capacity', 'read', 1e9, dependencies=DEPENDENCIES,
                    provenance=PROVENANCE, max_age_seconds=30.)
    with patch('torchgwas.calibration_cache.time.time', return_value=35.):
        assert cache.lookup('storage_capacity','read',dependencies=DEPENDENCIES)['hit']
        assert cache.lookup('storage_capacity','read',dependencies=DEPENDENCIES,max_age_seconds=20.)['reason']=='expired'
    with patch('torchgwas.calibration_cache.time.time', return_value=50.):
        assert cache.lookup('storage_capacity','read',dependencies=DEPENDENCIES,max_age_seconds=100.)['reason']=='expired'
    with patch('torchgwas.calibration_cache.time.time', return_value=5.):
        assert cache.lookup('storage_capacity','read',dependencies=DEPENDENCIES)['reason']=='future_timestamp'


@pytest.mark.parametrize('kind', ['available_memory','contention'])
def test_live_observations_are_never_reusable(cache, kind):
    record = cache.store(kind, 'cuda:0', 123, dependencies=DEPENDENCIES, provenance=PROVENANCE)
    assert not record['reusable']
    assert cache.lookup(kind,'cuda:0',dependencies=DEPENDENCIES)==dict(hit=False,reason='requires_live_observation')


def test_refresh_appends_and_corruption_is_not_accepted(cache):
    with patch('torchgwas.calibration_cache.time.time', return_value=10.):
        first=cache.store('cpu_capacity','decode',1e8,dependencies=DEPENDENCIES,
                          provenance=PROVENANCE,max_age_seconds=100.)
    old=open(first['path']).read()
    with patch('torchgwas.calibration_cache.time.time', return_value=20.):
        second=cache.store('cpu_capacity','decode',2e8,dependencies=DEPENDENCIES,
                           provenance=PROVENANCE,max_age_seconds=100.)
    with patch('torchgwas.calibration_cache.time.time', return_value=30.):
        assert cache.lookup('cpu_capacity','decode',dependencies=DEPENDENCIES)['record']['value']==2e8
        assert open(first['path']).read()==old and first['path']!=second['path']
        with open(second['path'],'w') as stream: stream.write('{"tampered":true}')
        result=cache.lookup('cpu_capacity','decode',dependencies=DEPENDENCIES)
        assert result['hit'] and result['record']['value']==1e8 and result['invalid_records']==1


@pytest.mark.parametrize('kind,age', [('gpu_capacity',None), ('gpu_capacity',0),
    ('gpu_capacity',True), ('gpu_capacity',float('nan')), ('workspace',10.)])
def test_invalid_freshness_policy_is_rejected(cache, kind, age):
    with pytest.raises(ValueError):
        cache.store(kind,'value',1,dependencies=DEPENDENCIES,provenance=PROVENANCE,max_age_seconds=age)


def test_same_record_publication_is_idempotent(cache):
    with patch('torchgwas.calibration_cache.time.time',return_value=10.):
        first=cache.store('source_work','ops',{'flops':20},dependencies=DEPENDENCIES,provenance=PROVENANCE)
        second=cache.store('source_work','ops',{'flops':20},dependencies=DEPENDENCIES,provenance=PROVENANCE)
    assert first==second
    assert len(list(cache.directory.glob('*/*.json')))==1


def test_delayed_publication_does_not_make_old_measurements_fresh(cache):
    with patch('torchgwas.calibration_cache.time.time',return_value=1000.):
        cache.store('gpu_capacity','old',1e9,dependencies=DEPENDENCIES,
            provenance=PROVENANCE,max_age_seconds=30.,observed_unix_seconds=900.)
        assert cache.lookup('gpu_capacity','old',dependencies=DEPENDENCIES)['reason']=='expired'
        cache.store('gpu_capacity','recent',1e9,dependencies=DEPENDENCIES,
            provenance=PROVENANCE,max_age_seconds=30.,observed_unix_seconds=990.)
        assert cache.lookup('gpu_capacity','recent',dependencies=DEPENDENCIES)['age_seconds']==10.
    with patch('torchgwas.calibration_cache.time.time',return_value=1020.):
        assert cache.lookup('gpu_capacity','recent',dependencies=DEPENDENCIES)['reason']=='expired'


@pytest.mark.parametrize('kind,observed',[('workspace',1.),('gpu_capacity',float('nan')),
    ('gpu_capacity',True),('gpu_capacity',2000.)])
def test_invalid_observation_time_is_rejected(cache, kind, observed):
    with patch('torchgwas.calibration_cache.time.time',return_value=1000.):
        with pytest.raises(ValueError):
            cache.store(kind,'invalid',1,dependencies=DEPENDENCIES,provenance=PROVENANCE,
                max_age_seconds=30. if kind=='gpu_capacity' else None,observed_unix_seconds=observed)

@pytest.mark.parametrize('recent_has_observation_time', [True, False])
def test_out_of_order_jobs_reuse_newest_measurement_not_latest_publication(cache, recent_has_observation_time):
    # The shorter job observes later and finishes first. Legacy records use
    # publication time as their observation time, so mixed caches also work.
    with patch('torchgwas.calibration_cache.time.time', return_value=100.):
        extra = dict(observed_unix_seconds=90.) if recent_has_observation_time else {}
        recent = cache.store('cpu_capacity', 'decode', 2e8, dependencies=DEPENDENCIES,
            provenance=dict(PROVENANCE, job='short-job'), max_age_seconds=100., **extra)
    with patch('torchgwas.calibration_cache.time.time', return_value=110.):
        old = cache.store('cpu_capacity', 'decode', 1e8, dependencies=DEPENDENCIES,
            provenance=dict(PROVENANCE, job='long-job'), max_age_seconds=100.,
            observed_unix_seconds=50.)
        result = cache.lookup('cpu_capacity', 'decode', dependencies=DEPENDENCIES)
        assert result['hit'] and result['record_sha256'] == recent['record_sha256']
        assert result['age_seconds'] == (20. if recent_has_observation_time else 10.)
        assert old['path'] != recent['path']
        assert len(list(cache.directory.glob('*/*.json'))) == 2


@pytest.mark.parametrize('kind', ['cpu_capacity', 'gpu_capacity', 'transfer_capacity',
    'storage_capacity', 'stage_observations'])
@pytest.mark.parametrize('requested_age', [None, 20.])
def test_lookup_cannot_return_a_record_expiring_during_cache_reads(cache, monkeypatch, kind, requested_age):
    now = [100.]
    monkeypatch.setattr('torchgwas.calibration_cache.time.time', lambda: now[0])
    saved = cache.store(kind, 'service', 3., dependencies=DEPENDENCIES,
        provenance=PROVENANCE, max_age_seconds=30., observed_unix_seconds=90.)
    original = Path(saved['path']).read_bytes()
    read = Path.read_text
    def delayed(path, *args, **kwargs):
        contents = read(path, *args, **kwargs)
        now[0] = 90. + (30. if requested_age is None else requested_age)
        return contents
    monkeypatch.setattr(Path, 'read_text', delayed)
    result = cache.lookup(kind, 'service', dependencies=DEPENDENCIES, max_age_seconds=requested_age)
    assert not result['hit'] and result['reason'] == 'expired'
    assert result['expired_records'] == 1
    assert Path(saved['path']).read_bytes() == original
    assert len(list(cache.directory.glob('*/*.json'))) == 1


def test_lookup_reports_age_after_cache_reads_without_renewing_observation(cache, monkeypatch):
    now = [100.]
    monkeypatch.setattr('torchgwas.calibration_cache.time.time', lambda: now[0])
    saved = cache.store('cpu_capacity', 'service', 3., dependencies=DEPENDENCIES,
        provenance=PROVENANCE, max_age_seconds=30., observed_unix_seconds=90.)
    original = Path(saved['path']).read_bytes()
    read = Path.read_text
    def delayed(path, *args, **kwargs):
        contents = read(path, *args, **kwargs)
        now[0] = 115.
        return contents
    monkeypatch.setattr(Path, 'read_text', delayed)
    hit = cache.lookup('cpu_capacity', 'service', dependencies=DEPENDENCIES)
    assert hit['hit'] and hit['record_sha256'] == saved['record_sha256']
    assert hit['age_seconds'] == 25.
    assert hit['record']['observed_unix_seconds'] == 90.
    assert hit['record']['created_unix_seconds'] == 100.
    assert Path(saved['path']).read_bytes() == original


@pytest.mark.parametrize('reverse_order', [False, True])
def test_record_expiring_during_lookup_does_not_hide_older_fresh_measurement(cache, monkeypatch, reverse_order):
    now = [80.]
    monkeypatch.setattr('torchgwas.calibration_cache.time.time', lambda: now[0])
    old = cache.store('gpu_capacity', 'service', 2., dependencies=DEPENDENCIES,
        provenance=PROVENANCE, max_age_seconds=100., observed_unix_seconds=70.)
    now[0] = 95.
    recent = cache.store('gpu_capacity', 'service', 3., dependencies=DEPENDENCIES,
        provenance=PROVENANCE, max_age_seconds=20., observed_unix_seconds=90.)
    paths = [Path(old['path']), Path(recent['path'])]
    original = {p:p.read_bytes() for p in paths}
    now[0] = 100.
    read = Path.read_text
    def delayed(path, *args, **kwargs):
        contents = read(path, *args, **kwargs)
        now[0] = 115.
        return contents
    monkeypatch.setattr(Path, 'read_text', delayed)
    monkeypatch.setattr(Path, 'glob', lambda path, pattern: iter(paths[::-1] if reverse_order else paths))
    hit = cache.lookup('gpu_capacity', 'service', dependencies=DEPENDENCIES)
    assert hit['hit'] and hit['record_sha256'] == old['record_sha256']
    assert hit['age_seconds'] == 45.
    assert {p:p.read_bytes() for p in paths} == original


@pytest.mark.parametrize('kind', ['cpu_capacity', 'source_work'])
def test_lookup_rejects_a_backwards_clock_step_during_read(cache, monkeypatch, kind):
    now = [100.]
    monkeypatch.setattr('torchgwas.calibration_cache.time.time', lambda: now[0])
    extra = dict(max_age_seconds=30., observed_unix_seconds=90.) if kind == 'cpu_capacity' else {}
    saved = cache.store(kind, 'service', 3., dependencies=DEPENDENCIES, provenance=PROVENANCE, **extra)
    original = Path(saved['path']).read_bytes()
    read = Path.read_text
    def delayed(path, *args, **kwargs):
        contents = read(path, *args, **kwargs)
        now[0] = 105.
        return contents
    now[0] = 110.
    monkeypatch.setattr(Path, 'read_text', delayed)
    result = cache.lookup(kind, 'service', dependencies=DEPENDENCIES)
    assert not result['hit'] and result['reason'] == 'clock_moved_backwards'
    assert Path(saved['path']).read_bytes() == original


def test_structural_cache_reads_do_not_create_an_empirical_expiry(cache, monkeypatch):
    now = [100.]
    monkeypatch.setattr('torchgwas.calibration_cache.time.time', lambda: now[0])
    saved = cache.store('source_work', 'ops', 3., dependencies=DEPENDENCIES, provenance=PROVENANCE)
    read = Path.read_text
    def delayed(path, *args, **kwargs):
        contents = read(path, *args, **kwargs)
        now[0] = 1e10
        return contents
    monkeypatch.setattr(Path, 'read_text', delayed)
    hit = cache.lookup('source_work', 'ops', dependencies=DEPENDENCIES)
    assert hit['hit'] and hit['record_sha256'] == saved['record_sha256']
    assert hit['age_seconds'] == 1e10-100.
