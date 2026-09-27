"""The execution bridge rechecks general component ages, including plan hits."""
import pytest

from test_jagwas_autotune import planning_fixture,controller
from test_detailed_autotune import source
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import bind_detailed_profile,sha256_file,validate_detailed_profile


def fixture(tmp_path,monkeypatch):
    import torchgwas.detailed_autotune as module
    now=[100.];monkeypatch.setattr('time.time',lambda:now[0])
    profile,config,path,calls,select,values=planning_fixture(tmp_path,monkeypatch)
    context=module.execution_context(None)
    profile['contexts'][0]['name']='test';profile['contexts'][0]['profiles']['cuda:1']['copy_cpu']=3.
    deps=dict(source_sha256=profile['source_sha256'],execution_context=context,
        measurement_protocol=dict(operation='independent-copy'))
    record=CalibrationParameterCache(tmp_path/'copy').store('cpu_capacity','copy',dict(rate=3.),
        dependencies=deps,provenance=dict(test='synthetic clock'),observed_unix_seconds=90.,max_age_seconds=15.)
    artifacts=dict(profile['component_artifacts'],**{record['path']:sha256_file(record['path'])})
    bound=bind_detailed_profile(profile['contexts'],context,sources=profile['source_sha256'],
        component_artifacts=artifacts,limitations=['Synthetic numerical controls'],price_bindings=[
            dict(artifact=record['path'],kind='cpu_capacity',name='copy',dependencies=deps,max_age_seconds=None,
                targets=[dict(context_path=[0,'profiles','cuda:1','copy_cpu'],value_path=['rate'])])])
    profile.clear();profile.update(bound)
    monkeypatch.setattr(module,'validate_detailed_profile',lambda p,c:
        validate_detailed_profile(p,c,sources=profile['source_sha256']))
    return now,profile,config,path,calls,select,values


@pytest.mark.parametrize('during',['qc','planning'])
def test_expired_nonreduction_price_prevents_selection(tmp_path,monkeypatch,during):
    import torchgwas.detailed_autotune as module
    now,p,c,path,calls,select,values=fixture(tmp_path,monkeypatch)
    now[0]=101.;tuner=controller(p,c,path)
    if during=='qc':now[0]=106.
    else:
        original=module.bounded_jagwas_plan
        def slow(*args,**kwargs):
            result=original(*args,**kwargs);now[0]=106.;return result
        monkeypatch.setattr(module,'bounded_jagwas_plan',slow)
    y,cov,qc,output=values
    with pytest.raises(ValueError,match='[Ee]xpired'):
        tuner.select(source(path),y,cov,qc,output=output)
    assert len(calls)==(during=='planning') and not (tmp_path/'output').exists()


def test_plan_hit_checks_original_record_and_reports_current_age(tmp_path,monkeypatch):
    now,p,c,path,calls,select,values=fixture(tmp_path,monkeypatch)
    now[0]=101.;first,audit,_=select();now[0]=104.;second,hit,_=select()
    assert first==second and len(calls)==1 and hit['analytical_cache']['status']=='hit'
    before=audit['price_evidence']['bindings'][0];after=hit['price_evidence']['bindings'][0]
    assert before['record_sha256']==after['record_sha256']
    assert before['age_seconds']==11. and after['age_seconds']==14.
    now[0]=106.
    with pytest.raises(ValueError,match='[Ee]xpired'):select()
    assert len(calls)==1
