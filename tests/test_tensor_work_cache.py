import copy
from unittest.mock import patch

import pytest
import torch

import torchgwas.tensor_work as work


def test_scoped_trace_reuse_preserves_full_ledger_and_mutation_independence():
    expected=work.eager_statistics_work(32,4,3,8,True)
    with patch.object(work,'_trace_statistics_work',wraps=work._trace_statistics_work) as trace:
        with work.tensor_work_cache():
            first=work.eager_statistics_work(32,4,3,8,True)
            second=work.eager_statistics_work(32,4,3,8,True)
            assert first==second==expected and trace.call_count==1
            first['steps'][0]['inputs'][0]['shape'][0]=-1
            first['source_sha256'].clear()
            assert work.eager_statistics_work(32,4,3,8,True)==second==expected
        assert work.eager_statistics_work(32,4,3,8,True)==expected
        assert trace.call_count==2


def test_trace_cache_is_bounded_and_distinguishes_geometry_validation_and_dtype():
    def fake(n,b,k,c,validate):
        return dict(n=n,b=b,k=k,c=c,validate=validate,dtype=str(torch.get_default_dtype()))
    with patch.object(work,'_trace_statistics_work',side_effect=fake) as trace:
        with work.tensor_work_cache(max_entries=2):
            for n,b,k,c,v in [(32,4,3,8,False),(32,4,3,8,True),(32,4,4,8,True)]:
                work.eager_statistics_work(n,b,k,c,v)
            work.eager_statistics_work(32,4,3,8,False)
            assert trace.call_count==4  # the first shape was evicted
            old=torch.get_default_dtype()
            try:
                torch.set_default_dtype(torch.float64)
                assert work.eager_statistics_work(32,4,3,8,False)['dtype']=='torch.float64'
                assert trace.call_count==5
            finally:torch.set_default_dtype(old)


def test_nested_and_failed_plans_cannot_leak_trace_cache():
    with patch.object(work,'_trace_statistics_work',return_value={'steps':[]}) as trace:
        with work.tensor_work_cache():
            work.eager_statistics_work(32,4)
            with pytest.raises(RuntimeError),work.tensor_work_cache():
                work.eager_statistics_work(32,4)
                raise RuntimeError('failed nested plan')
            work.eager_statistics_work(32,4)
            assert trace.call_count==2
        work.eager_statistics_work(32,4)
        assert trace.call_count==3
    assert work._STATISTICS_WORK_CACHE.get() is None


def test_decorated_plan_has_a_fresh_cache_on_every_call_and_on_failure():
    @work.reuse_tensor_work
    def plan(fail=False):
        value=work.eager_statistics_work(32,4)
        work.eager_statistics_work(32,4)
        if fail:raise RuntimeError('invalid candidate')
        return value
    with patch.object(work,'_trace_statistics_work',return_value={'steps':[]}) as trace:
        plan();plan()
        with pytest.raises(RuntimeError):plan(True)
        assert trace.call_count==3
    assert work._STATISTICS_WORK_CACHE.get() is None


def test_invalid_cache_bounds():
    for maximum in [True,0,-1,1.5]:
        with pytest.raises(ValueError),work.tensor_work_cache(maximum):pass