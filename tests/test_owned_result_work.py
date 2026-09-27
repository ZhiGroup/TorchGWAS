import pytest
from torchgwas.owned_result_work import owned_result_work,owned_result_copy_service


def test_owned_result_layout_and_explicit_copy_states():
    work=owned_result_work(1024,512)
    assert work['allocation_calls']==4
    assert work['copy_bytes']==(8*512+5)*1024
    scenario=dict(resident_cpu_seconds_per_byte=1.,fresh_cpu_seconds_per_byte=6.,fresh_fraction=0.)
    resident=owned_result_copy_service(work,scenario)
    scenario['fresh_fraction']=1.
    fresh=owned_result_copy_service(work,scenario)
    assert fresh['additional_copy_cpu_seconds']==6*resident['additional_copy_cpu_seconds']
    assert fresh['additional_copy_bytes']==work['copy_bytes']-416
    assert not fresh['prediction_complete']
    assert owned_result_copy_service(owned_result_work(1,1),scenario)['additional_copy_cpu_seconds']==0


@pytest.mark.parametrize('fraction',[-1,1.1,float('nan'),True])
def test_invalid_memory_state_is_not_silently_clamped(fraction):
    with pytest.raises(ValueError):
        owned_result_copy_service(owned_result_work(1,1),dict(resident_cpu_seconds_per_byte=1.,fresh_cpu_seconds_per_byte=6.,fresh_fraction=fraction))


def test_numpy_advice_threshold_is_per_allocation_and_inclusive():
    from torchgwas.owned_result_work import owned_result_allocation_policy
    policy=dict(numpy_version='2.2.6',numpy_madvise_hugepage=True,linux_thp_enabled='madvise')
    below=owned_result_allocation_policy(owned_result_work(1024,512),policy)
    assert below['advice_eligible_array_bytes']=={}
    assert below['advice_compaction_exposure'] is False
    at=owned_result_allocation_policy(owned_result_work(2048,512),policy)
    assert at['advice_eligible_array_bytes']==dict(beta=4*1024**2,t=4*1024**2)
    assert at['advice_compaction_exposure'] is True
    assert at['compaction_latency_bound_seconds'] is None


def test_unknown_allocator_version_or_policy_is_not_assumed_safe():
    from torchgwas.owned_result_work import owned_result_allocation_policy
    work=owned_result_work(4096,512)
    assert owned_result_allocation_policy(work)['advice_compaction_exposure'] is None
    unknown=owned_result_allocation_policy(work,dict(numpy_version='future',numpy_madvise_hugepage=True,linux_thp_enabled='madvise'))
    assert unknown['verified_advice_threshold_bytes'] is None
    assert unknown['advice_compaction_exposure'] is None
    disabled=owned_result_allocation_policy(work,dict(numpy_madvise_hugepage=False))
    assert disabled['advice_compaction_exposure'] is False
    assert not disabled['prediction_complete']
    with pytest.raises(ValueError):
        owned_result_allocation_policy(work,dict(numpy_madvise_hugepage='0'))

def test_numeric_copy_gil_threshold_and_unknown_build():
    from torchgwas.owned_result_work import owned_result_copy_threading
    policy=dict(numpy_version='2.2.6',numpy_allow_threads=True)
    at=owned_result_copy_threading(owned_result_work(500,1),policy)
    above=owned_result_copy_threading(owned_result_work(501,1),policy)
    assert set(at['gil_released_by_array'].values())=={False}
    assert set(above['gil_released_by_array'].values())=={True}
    mixed=owned_result_copy_threading(owned_result_work(1,501),policy)
    assert mixed['gil_released_by_array']==dict(beta=True,t=True,status=False,df=False)
    assert set(owned_result_copy_threading(owned_result_work(501,1))['gil_released_by_array'].values())=={None}
    with pytest.raises(ValueError):
        owned_result_copy_threading(owned_result_work(501,1),dict(numpy_allow_threads='yes'))