"""Cross-job operation reuse never reuses a priced result or refreshes evidence."""
from copy import deepcopy
import json
from pathlib import Path
from unittest.mock import patch

import pytest
import torch

from torchgwas import structural_tensor_cache as saved
from torchgwas import tensor_work
from torchgwas.reduction_tensor_work import jagwas_tensor_work
from torchgwas.structural_tensor_cache import StructuralTensorWorkCache


def statistics():return tensor_work.eager_statistics_work(32,4,3,8,True)


def traces():
    return [statistics(),jagwas_tensor_work(32,4,3,phase='prepare'),
            jagwas_tensor_work(32,4,3,phase='reduce')]


def files(directory):return {str(p):p.read_bytes() for p in directory.glob('*/*.json')}


def test_structural_traces_are_staged_then_reused_by_another_job_without_refresh(tmp_path):
    expected=traces();first=StructuralTensorWorkCache(tmp_path)
    with first.activate():
        actual=traces();assert actual==expected
        actual[0]['steps'].clear()
        assert statistics()==expected[0]
    assert not files(tmp_path)
    assert first.snapshot()['pending']==3 and first.snapshot()['memory_hits']==1
    published=first.publish(successful=True);assert len(published['stored'])==3
    before=files(tmp_path);first.close()
    assert all(json.loads(v)['kind']=='source_work' and json.loads(v)['max_age_seconds'] is None for v in before.values())
    later=StructuralTensorWorkCache(tmp_path)
    with patch('torchgwas.calibration_cache.time.time',return_value=1e11):
        with later.activate():assert traces()==expected
        assert later.snapshot()['disk_hits']==3 and later.snapshot()['misses']==0
        assert later.publish(successful=True)['stored']==[]
    assert files(tmp_path)==before
    assert all('observed_unix_seconds' not in json.loads(v) for v in before.values())


@pytest.mark.parametrize('shape',[(33,4,3,8,True),(32,2,3,8,True),(32,4,2,8,True),
    (32,4,3,7,True),(32,4,3,8,False)])
def test_changed_shape_or_validation_has_its_own_source_ledger(tmp_path,shape):
    first=StructuralTensorWorkCache(tmp_path)
    with first.activate():statistics()
    first.publish(successful=True)
    later=StructuralTensorWorkCache(tmp_path)
    with later.activate():actual=tensor_work.eager_statistics_work(*shape)
    assert later.snapshot()['misses']==1 and later.snapshot()['disk_hits']==0
    assert actual==tensor_work.eager_statistics_work(*shape)


def test_dtype_library_source_and_implementation_changes_invalidate_reuse(tmp_path):
    first=StructuralTensorWorkCache(tmp_path)
    with first.activate():statistics()
    first.publish(successful=True)
    old=torch.get_default_dtype()
    try:
        torch.set_default_dtype(torch.float64);later=StructuralTensorWorkCache(tmp_path)
        with later.activate():actual=statistics()
        assert later.snapshot()['misses']==1 and actual==statistics()
    finally:torch.set_default_dtype(old)
    runtime=saved._runtime();runtime['torch_git']='different-build'
    with patch.object(saved,'_runtime',return_value=runtime):
        later=StructuralTensorWorkCache(tmp_path)
        with later.activate():assert statistics()==statistics()
        assert later.snapshot()['misses']==1 and later.snapshot()['memory_hits']==1
    sources=dict(saved._sources(),**{'linear.py':'different-source'})
    with patch.object(saved,'_sources',return_value=sources):
        later=StructuralTensorWorkCache(tmp_path)
        with later.activate():statistics()
        assert later.snapshot()['misses']==1
    # A runtime replacement absent from the source hash is also incompatible.
    later=StructuralTensorWorkCache(tmp_path)
    with patch.object(tensor_work,'_trace_statistics_work',wraps=tensor_work._trace_statistics_work) as trace:
        with later.activate():statistics();statistics()
        assert trace.call_count==2 and later.snapshot()['bypasses']==2
    assert len(files(tmp_path))==1


def test_jagwas_phase_and_compute_dtype_are_distinct(tmp_path):
    first=StructuralTensorWorkCache(tmp_path)
    with first.activate():
        a=jagwas_tensor_work(32,4,3,phase='prepare')
        b=jagwas_tensor_work(32,4,3,phase='reduce')
        c=jagwas_tensor_work(32,4,3,phase='reduce',compute_dtype='float64')
    assert first.snapshot()['misses']==3 and a!=b and b!=c
    assert len(first.publish(successful=True)['stored'])==3


def test_implementation_identity_is_stable_across_trace_execution():
    from torchgwas.reduction_tensor_work import _trace_jagwas_tensor_work
    from torchgwas.reduce import JagwasReduction
    functions=(_trace_jagwas_tensor_work,JagwasReduction.prepare,JagwasReduction.reduce)
    before=saved._dependencies(saved._sources(),dict(samples=32,markers=4,traits=3),functions)
    for _ in range(3):
        traces()
        assert saved._dependencies(saved._sources(),dict(samples=32,markers=4,traits=3),functions)==before


def test_service_prices_are_recalculated_with_saved_source_work(tmp_path):
    from test_mechanistic_shapes import fixture,component
    from torchgwas.mechanistic_torch import torch_scan_work
    data,profile=fixture();first=StructuralTensorWorkCache(tmp_path)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component) as price:
        with first.activate():torch_scan_work(data,profile)
        assert price.call_count==2
        first.publish(successful=True);before=files(tmp_path)
        later=StructuralTensorWorkCache(tmp_path);changed=deepcopy(profile)
        changed['gpu_resources']['hbm_bytes_per_second']/=2
        with later.activate():torch_scan_work(data,changed)
        assert price.call_count==4 and later.snapshot()['disk_hits']==2
        later.publish(successful=True)
        assert files(tmp_path)==before
    assert all('gpu_resources' not in json.loads(raw)['dependencies'] for raw in before.values())


def test_bounds_eviction_and_close_drop_unpublished_work(tmp_path):
    bounded=StructuralTensorWorkCache(tmp_path,max_entries=1)
    with bounded.activate():
        statistics();jagwas_tensor_work(32,4,3,phase='reduce')
    assert bounded.snapshot()['entries']==1 and bounded.snapshot()['evictions']==1
    assert len(bounded.publish(successful=True)['stored'])==1
    tiny=StructuralTensorWorkCache(tmp_path/'tiny',max_bytes=1024)
    with tiny.activate():statistics()
    assert tiny.snapshot()['entries']==0 and tiny.snapshot()['bypasses']==1
    assert tiny.publish(successful=True)['stored']==[]
    pending=StructuralTensorWorkCache(tmp_path/'unpublished')
    with pending.activate():statistics()
    pending.close()
    assert pending.snapshot()['estimated_retained_bytes']==0 and not files(tmp_path/'unpublished')
    with pytest.raises(ValueError),pending.activate():pass
    with pytest.raises(ValueError):pending.publish(successful=True)


def test_failed_work_source_changes_and_cache_failures_cannot_publish_good_looking_records(tmp_path):
    cache=StructuralTensorWorkCache(tmp_path)
    with patch.object(cache.cache,'lookup',side_effect=OSError('unavailable')):
        with cache.activate():assert statistics()==statistics()
    assert cache.snapshot()['cache_errors']==1
    assert cache.publish(successful=False)['status']=='unsuccessful' and not files(tmp_path)
    with patch.object(saved,'_sources',return_value={}):
        assert cache.publish(successful=True)['status']=='source_changed' and not files(tmp_path)
    with patch.object(cache.cache,'store',side_effect=OSError('unavailable')):
        assert cache.publish(successful=True)['status']=='partial' and not files(tmp_path)
    assert cache.snapshot()['pending']==1 and cache.snapshot()['cache_errors']==2
    assert len(cache.publish(successful=True)['stored'])==1


def test_corrupt_saved_ledger_is_recomputed_without_modifying_the_old_record(tmp_path):
    first=StructuralTensorWorkCache(tmp_path)
    with first.activate():expected=statistics()
    published=first.publish(successful=True);path=Path(published['stored'][0]['path'])
    path.write_text('{"invalid":true}');damaged=path.read_bytes()
    later=StructuralTensorWorkCache(tmp_path)
    with later.activate():assert statistics()==expected
    assert later.snapshot()['misses']==1 and later.snapshot()['disk_hits']==0
    assert len(later.publish(successful=True)['stored'])==1
    assert path.read_bytes()==damaged and len(files(tmp_path))==2


def test_nested_scopes_and_failed_trace_do_not_leak_or_cache_failure(tmp_path):
    outer=StructuralTensorWorkCache(tmp_path/'outer');inner=StructuralTensorWorkCache(tmp_path/'inner')
    with outer.activate():
        statistics()
        with pytest.raises(ValueError),inner.activate():
            jagwas_tensor_work(32,4,3,phase='invalid')
        assert inner.snapshot()['entries']==0
        statistics()
    assert outer.snapshot()['memory_hits']==1 and saved._ACTIVE.get() is None


@pytest.mark.parametrize('options',[dict(max_entries=True),dict(max_entries=0),dict(max_bytes=0)])
def test_invalid_cache_limits_are_rejected(tmp_path,options):
    with pytest.raises(ValueError):StructuralTensorWorkCache(tmp_path,**options)
