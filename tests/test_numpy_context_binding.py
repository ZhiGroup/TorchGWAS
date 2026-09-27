"""Same-version NumPy builds and CPU features must not share old prices."""
from copy import deepcopy
import json
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas import binding_digests,detailed_calibration as calibration
from torchgwas.binding_digests import BindingDigestCache
from test_detailed_calibration import fixture as profile_fixture


@pytest.fixture
def core(tmp_path,monkeypatch):
    path=tmp_path/'numpy_core.so';path.write_bytes(b'first-build')
    value=SimpleNamespace(__file__=str(path),__cpu_features__={'SSE2':True,'AVX2':True},
        __cpu_baseline__=['SSE2'],__cpu_dispatch__=['AVX2'])
    monkeypatch.setattr(np._core,'_multiarray_umath',value)
    monkeypatch.setattr(binding_digests,'_host_identity',lambda:'fixture-boot')
    monkeypatch.setattr(binding_digests,'_stable',lambda row:True)
    return value,path


def test_cross_job_numpy_digest_reuse_keeps_runtime_metadata_fresh(core,tmp_path):
    module,path=core;directory=tmp_path/'digests'
    first=BindingDigestCache(directory)
    original=calibration._numpy_core_context(digest_cache=first)
    assert original==calibration._numpy_core_context()
    first.publish(successful=True);first.close()
    artifacts={p:p.read_bytes() for p in directory.rglob('*.json')}
    later=BindingDigestCache(directory)
    with patch.object(calibration,'sha256_file',side_effect=AssertionError('unchanged binary reread')):
        assert calibration._numpy_core_context(digest_cache=later)==original
        module.__cpu_features__['AVX2']=False
        changed=calibration._numpy_core_context(digest_cache=later)
    assert original['cpu_features']['AVX2'] is True and changed['cpu_features']['AVX2'] is False
    assert later.snapshot()['disk_hits']==later.snapshot()['memory_hits']==1
    assert later.publish(successful=True)['stored']==[]
    assert artifacts=={p:p.read_bytes() for p in directory.rglob('*.json')}


def test_same_version_binary_replacement_invalidates_reused_price(core,tmp_path):
    module,path=core;cache=BindingDigestCache(tmp_path/'digests')
    profile,execution,sources,_=profile_fixture(tmp_path)
    execution.update(numpy_version=np.__version__,numpy_core=calibration._numpy_core_context(digest_cache=cache))
    profile['execution_context']=deepcopy(execution)
    saved=json.dumps(profile,sort_keys=True)
    replacement=path.with_suffix('.new');replacement.write_bytes(b'other-build');replacement.replace(path)
    current=dict(execution,numpy_core=calibration._numpy_core_context(digest_cache=cache))
    with pytest.raises(ValueError,match='numpy_core.library_sha256'):
        calibration.validate_detailed_profile(profile,current,sources=sources)
    assert cache.publish(successful=True)['status']=='files_changed'
    assert json.dumps(profile,sort_keys=True)==saved


@pytest.mark.parametrize('field,value,changed',[
    ('__cpu_features__',{'SSE2':True,'AVX2':False},'cpu_features.AVX2'),
    ('__cpu_baseline__',['SSE2','SSE3'],'cpu_baseline'),
    ('__cpu_dispatch__',[],'cpu_dispatch')])
def test_cpu_dispatch_changes_invalidate_same_binary_profile(core,tmp_path,field,value,changed):
    module,_=core;profile,execution,sources,_=profile_fixture(tmp_path)
    execution['numpy_core']=calibration._numpy_core_context()
    profile['execution_context']=deepcopy(execution)
    setattr(module,field,value)
    execution['numpy_core']=calibration._numpy_core_context()
    with pytest.raises(ValueError,match='numpy_core.'+changed):
        calibration.validate_detailed_profile(profile,execution,sources=sources)


def test_legacy_profile_without_binary_context_is_not_silently_upgraded(core,tmp_path):
    profile,execution,sources,_=profile_fixture(tmp_path)
    execution['numpy_core']=calibration._numpy_core_context()
    with pytest.raises(ValueError,match='context.numpy_core'):
        calibration.validate_detailed_profile(profile,execution,sources=sources)
    assert 'numpy_core' not in profile['execution_context']


@pytest.mark.parametrize('field,value',[
    ('__file__',None),('__file__',''),('__cpu_features__',None),('__cpu_features__',{}),
    ('__cpu_features__',{'AVX2':1}),('__cpu_features__',{'':True}),
    ('__cpu_baseline__',None),('__cpu_dispatch__','AVX2'),('__cpu_dispatch__',[1])])
def test_missing_or_unrecognized_runtime_metadata_fails_closed(core,field,value):
    module,_=core;setattr(module,field,value)
    with pytest.raises(ValueError,match='NumPy'):
        calibration._numpy_core_context()
