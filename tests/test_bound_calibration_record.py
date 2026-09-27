"""A bound empirical artifact cannot silently change identity or renew age."""
import json
from pathlib import Path
from unittest.mock import patch
import pytest
from torchgwas.calibration_cache import CalibrationParameterCache,read_calibration_record

ARGS=dict(kind='cpu_capacity',name='significant_host_components',dependencies={'source':'one'})


def saved(tmp_path, observed=90., lifetime=20.):
    cache=CalibrationParameterCache(tmp_path)
    with patch('torchgwas.calibration_cache.time.time',return_value=100.):
        value=cache.store(ARGS['kind'],ARGS['name'],{'price':3.},dependencies=ARGS['dependencies'],
            provenance={'job':'independent'},observed_unix_seconds=observed,max_age_seconds=lifetime)
    return cache,Path(value['path'])


def test_exact_bound_record_stays_immutable_when_newer_record_is_published(tmp_path):
    cache,path=saved(tmp_path)
    original=path.read_bytes()
    with patch('torchgwas.calibration_cache.time.time',return_value=101.):
        cache.store(ARGS['kind'],ARGS['name'],{'price':9.},dependencies=ARGS['dependencies'],
            provenance={'job':'next'},observed_unix_seconds=101.,max_age_seconds=20.)
        record=read_calibration_record(path,**ARGS)
        assert record['record']['value']['price']==3. and record['age_seconds']==11.
        assert cache.lookup(**ARGS)['record']['value']['price']==9.
    assert path.read_bytes()==original


@pytest.mark.parametrize('now,override',[(110.,None),(111.,1000.),(100.,5.)])
def test_observation_expiry_cannot_be_extended_by_rebinding_or_override(tmp_path,now,override):
    _,path=saved(tmp_path)
    with patch('torchgwas.calibration_cache.time.time',return_value=now):
        with pytest.raises(ValueError,match='Expired'):
            read_calibration_record(path,**ARGS,max_age_seconds=override)


@pytest.mark.parametrize('fault',['kind','name','dependency','hash','legacy','future','duplicate'])
def test_incompatible_or_unauditable_bound_artifacts_fail(tmp_path,fault):
    _,path=saved(tmp_path,observed=None if fault=='legacy' else 90.)
    args=dict(ARGS)
    if fault=='kind':args['kind']='gpu_capacity'
    if fault=='name':args['name']='different'
    if fault=='dependency':args['dependencies']={'source':'changed'}
    if fault=='hash':
        value=json.loads(path.read_text());value['value']['price']=999.
        path.write_text(json.dumps(value))
    if fault=='duplicate':path.write_text(path.read_text().replace('"value":','"value":{},"value":'))
    with patch('torchgwas.calibration_cache.time.time',return_value=99. if fault=='future' else 101.):
        with pytest.raises(ValueError):read_calibration_record(path,**args)
