"""Price identity and freshness survive profile publication and reuse."""
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pytest

from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import (PRICE_SCHEMA,bind_detailed_profile,
    read_detailed_profile,sha256_file,validate_detailed_profile,write_detailed_profile)
from torchgwas.price_binding import validate_price_bindings


@pytest.fixture
def evidence(tmp_path,monkeypatch):
    now=[100.];monkeypatch.setattr('time.time',lambda:now[0])
    source={'model.py':'actual-source'}
    execution=dict(devices={'cuda:0':dict(uuid='card-zero')},affinity=[1,2])
    contexts=[dict(name='one',devices=['cuda:0'],profiles={'cuda:0':dict(rate=3.,other=8.)})]
    deps=dict(source_sha256=source,execution_context=execution,measurement_protocol=dict(operation='independent-copy'))
    cache=CalibrationParameterCache(tmp_path/'cache')
    record=cache.store('cpu_capacity','copy',dict(cpu_seconds_per_byte=3.),dependencies=deps,
        provenance=dict(job='original'),max_age_seconds=30.,observed_unix_seconds=90.)
    binding=dict(artifact=record['path'],kind='cpu_capacity',name='copy',dependencies=deps,max_age_seconds=None,
        targets=[dict(context_path=[0,'profiles','cuda:0','rate'],value_path=['cpu_seconds_per_byte'])])
    def bind(bindings=None):
        return bind_detailed_profile(contexts,execution,sources=source,limitations=['Explicit test prices'],
            component_artifacts={record['path']:sha256_file(record['path'])},
            price_bindings=[binding] if bindings is None else bindings)
    return dict(now=now,source=source,execution=execution,contexts=contexts,record=record,
        binding=binding,bind=bind,cache=cache,deps=deps)


def validate(profile,e):
    return validate_detailed_profile(profile,e['execution'],sources=e['source'])


def test_profile_links_exact_price_and_preserves_age_through_rebinding(evidence,tmp_path):
    e=evidence;p=e['bind']();original=Path(e['record']['path']).read_bytes()
    assert p['schema']==PRICE_SCHEMA
    e['now'][0]=110.
    new=e['bind']();report=validate(new,e)['price_evidence'];row=report['bindings'][0]
    assert report['status']=='declared_targets_verified' and report['verified_targets']==1
    assert row['observed_unix_seconds']==90. and row['age_seconds']==20.
    assert row['expires_unix_seconds']==120.
    assert row['record_sha256']==e['record']['record_sha256']
    path=tmp_path/'profile.json';write_detailed_profile(new,path)
    assert read_detailed_profile(path)==new and Path(e['record']['path']).read_bytes()==original
    e['now'][0]=120.
    with pytest.raises(ValueError,match='[Ee]xpired'):validate(p,e)
    with pytest.raises(ValueError,match='[Ee]xpired'):e['bind']()


def test_newer_record_does_not_replace_bound_price(evidence):
    e=evidence;p=e['bind']();e['now'][0]=110.
    e['cache'].store('cpu_capacity','copy',dict(cpu_seconds_per_byte=9.),dependencies=e['deps'],
        provenance=dict(job='new'),max_age_seconds=30.,observed_unix_seconds=110.)
    row=validate(p,e)['price_evidence']['bindings'][0]
    assert row['record_sha256']==e['record']['record_sha256'] and row['age_seconds']==20.
    e['now'][0]=121.
    assert e['cache'].lookup('cpu_capacity','copy',dependencies=e['deps'])['hit']
    with pytest.raises(ValueError,match='[Ee]xpired'):validate(p,e)


def test_caller_can_only_shorten_producer_age(evidence):
    e=evidence;p=e['bind']();p['price_bindings'][0]['max_age_seconds']=1000.
    assert validate(p,e)['price_evidence']['bindings'][0]['max_age_seconds']==30.
    with pytest.raises(ValueError,match='expired'):validate_price_bindings(p,max_age_seconds=10.)
    p['price_bindings'][0]['max_age_seconds']=10.
    with pytest.raises(ValueError,match='expired'):validate(p,e)


@pytest.mark.parametrize('change,match',[
    (lambda p:p['contexts'][0]['profiles']['cuda:0'].update(rate=4.),'Calculator price'),
    (lambda p:p['price_bindings'][0]['dependencies'].update(source_sha256={}),'source'),
    (lambda p:p['price_bindings'][0]['dependencies'].update(execution_context={}), 'execution'),
    (lambda p:p['price_bindings'][0]['dependencies'].update(measurement_protocol={}), 'protocol'),
    (lambda p:p['price_bindings'][0].update(kind='stage_observations'),'independent'),
    (lambda p:p['price_bindings'][0].update(kind='available_memory'),'independent'),
    (lambda p:p['price_bindings'][0].update(kind='source_work'),'independent'),
    (lambda p:p['price_bindings'][0]['targets'][0].update(value_path=['missing']),'path'),
    (lambda p:p['price_bindings'][0]['targets'][0].update(context_path=[True,'profiles']),'path'),
    (lambda p:p['price_bindings'][0]['targets'][0].update(context_path=[-1,'profiles']),'path'),
    (lambda p:p['price_bindings'][0]['targets'][0].update(context_path=[0,'profiles','missing']),'path'),
    (lambda p:p['price_bindings'][0].update(targets=[]),'targets'),
    (lambda p:p['price_bindings'][0].update(max_age_seconds=True),'age'),
    (lambda p:p.update(price_bindings=[]),'bindings'),
    (lambda p:p['price_bindings'].append(deepcopy(p['price_bindings'][0])),'Overlapping'),
])
def test_changed_evidence_and_unbound_values_are_refused(evidence,change,match):
    p=evidence['bind']();change(p)
    with pytest.raises(ValueError,match=match):validate(p,evidence)


def test_artifact_identity_check_cannot_be_disabled_for_prices(evidence):
    e=evidence;p=e['bind']();Path(e['record']['path']).write_text('changed')
    with pytest.raises(ValueError,match='artifact changed'):
        validate_detailed_profile(p,e['execution'],sources=e['source'],verify_artifacts=False)


def test_record_must_already_belong_to_bound_profile(evidence,tmp_path):
    e=evidence;p=e['bind']();other=tmp_path/'other.json';other.write_bytes(Path(e['record']['path']).read_bytes())
    p['price_bindings'][0]['artifact']=str(other)
    with pytest.raises(ValueError,match='absent'):validate(p,e)


def test_final_age_check_catches_expiry_while_evidence_is_read(evidence):
    import torchgwas.price_binding as module
    e=evidence;p=e['bind']();read=module.read_calibration_record
    def slow(*args,**kwargs):
        value=read(*args,**kwargs);e['now'][0]=120.;return value
    with patch.object(module,'read_calibration_record',side_effect=slow):
        with pytest.raises(ValueError,match='expired during'):validate(p,e)


def test_profile_without_declared_prices_remains_explicitly_unchecked(evidence):
    e=evidence;p=e['bind']();p.pop('price_bindings');p['schema']='torchgwas.detailed_calibration.v1'
    row=validate(p,e)['price_evidence']
    assert row['status']=='undeclared' and row['verified_targets']==0
    p['schema']=PRICE_SCHEMA
    with pytest.raises(ValueError,match='fields'):validate(p,e)


def test_bindings_and_reports_do_not_alias_caller_state(evidence):
    e=evidence;p=e['bind']();report=validate(p,e)
    e['binding']['targets'].clear();report['price_evidence']['bindings'][0]['targets'].clear()
    assert validate(p,e)['price_evidence']['verified_targets']==1


def test_target_subtrees_may_not_overlap(evidence):
    e=evidence;p=e['bind']()
    p['price_bindings'][0]['targets'].append(dict(context_path=[0,'profiles','cuda:0'],value_path=['cpu_seconds_per_byte']))
    with pytest.raises(ValueError,match='Overlapping'):validate(p,e)


def test_reusing_setup_profile_preserves_price_record_and_original_age(evidence,tmp_path,monkeypatch):
    from test_covariate_setup_calibration import bank
    import torchgwas.setup_calibration as module
    e=evidence;e['contexts'][0]['profiles']['cuda:0']['setup_primitives']=bank(8)
    profile=e['bind']();e['now'][0]=110.
    monkeypatch.setattr(module,'execution_context',lambda *a,**kw:e['execution'])
    monkeypatch.setattr(module,'source_identity',lambda:e['source'])
    monkeypatch.setattr(module,'validate_detailed_profile',lambda p,c:validate_detailed_profile(p,c,sources=e['source']))
    monkeypatch.setattr(module,'bind_detailed_profile',lambda *a,**kw:bind_detailed_profile(*a,**kw,sources=e['source']))
    monkeypatch.setattr(module.subprocess,'run',lambda *a,**kw:pytest.fail('References already available'))
    module.complete_setup_primitives(profile,[8],input_path='unused',output_path='unused',collection_dir=tmp_path/'derived')
    derived=read_detailed_profile(tmp_path/'derived'/'profile.json')
    assert derived['price_bindings']==profile['price_bindings']
    row=validate(derived,e)['price_evidence']['bindings'][0]
    assert row['observed_unix_seconds']==90. and row['age_seconds']==20.
    assert row['record_sha256']==e['record']['record_sha256']
