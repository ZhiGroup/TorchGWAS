from torchgwas.tensor_memory import eager_scan_memory

def test_more_outstanding_chunks_add_staging_and_retained_outputs():
    a=eager_scan_memory(32,8,3,2,2)
    b=eager_scan_memory(32,8,3,2,3)
    assert b['distinct_temporary_bytes']==a['distinct_temporary_bytes']
    assert b['tensor_storage_budget']-a['tensor_storage_budget']==32*8+a['retained_output_bytes']

def test_hidden_promotion_and_views_are_accounted_as_storage():
    a=eager_scan_memory(32,8,3,2,2)
    assert a['distinct_temporary_bytes']>=8*32*8
    assert a['distinct_temporary_storages']>0
    assert a['tensor_storage_budget']==sum(a[k] for k in ['resident_bytes','staging_bytes','distinct_temporary_bytes','retained_output_bytes','previous_genotype_bytes'])

def test_reference_lifetimes_reduce_the_allocation_sum():
    a=eager_scan_memory(128,32,17,8,3)
    assert a['lifetime_trace']['temporary_live_bytes']>8*128*32
    assert a['lifetime_candidate_bytes']<a['tensor_storage_budget']
    assert a['lifetime_trace']['peak_operation']


def test_library_composition_requires_explicit_profile():
    import pytest
    from torchgwas.tensor_memory import eager_memory_plan
    with pytest.raises(ValueError,match='Missing'):
        eager_memory_plan(128,32,17,8,3,{})
    profile=dict(torch_version='2.5.1+cu124',sm_count=132,max_threads_per_sm=2048,
        compute_capability=[9,0],cublas_workspace_config=None,cublas_handle_stream_pairs=1)
    one=eager_memory_plan(128,32,17,8,3,profile)
    profile['cublas_handle_stream_pairs']=2
    two=eager_memory_plan(128,32,17,8,3,profile)
    assert two['device_bytes']-one['device_bytes']==32*1024**2
    assert not two['prediction_complete']
    assert two['lifetime_candidate_bytes']<=two['device_bytes']
    profile['torch_version']='2.6.0'
    with pytest.raises(ValueError,match='2.5.1'):
        eager_memory_plan(128,32,17,8,3,profile)
