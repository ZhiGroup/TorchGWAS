import pytest
from torchgwas.reduction_memory import column_sum_workspace

def test_output_vectorized_global_buffer():
    r=column_sum_workspace(4096,2048,sm_count=132,max_threads_per_sm=2048)
    assert r['block']==[32,4] and r['grid']==[16,64]
    assert r['workspace_bytes']==67108864
    assert r['semaphore_bytes']==64

def test_small_input_needs_no_global_reduction():
    r=column_sum_workspace(32,17,sm_count=132,max_threads_per_sm=2048)
    assert r['workspace_bytes']==0 and r['semaphore_bytes']==0
    assert r['output_vector_width']==1

def test_split_tensor_is_not_silently_extrapolated():
    with pytest.raises(ValueError,match='splitting'):
        column_sum_workspace(2**20,2**12,sm_count=132,max_threads_per_sm=2048)

def test_setup_includes_reduction_workspace_only_with_explicit_device():
    from torchgwas.setup_work import setup_memory
    old=setup_memory(4096,2048,8)
    fixed=setup_memory(4096,2048,8,sm_count=132,max_threads_per_sm=2048)
    assert old['design_reduction_workspace_bytes'] is None
    assert fixed['design_reduction_workspace_bytes']==67108864+64
    assert fixed['design_device_live_bytes_upper']==old['design_device_live_bytes_upper']+67108864+64
    with pytest.raises(ValueError,match='Both'):
        setup_memory(4096,2048,8,sm_count=132)


def test_reduction_allocation_is_not_reduction_traffic():
    r=column_sum_workspace(4096,2048,sm_count=132,max_threads_per_sm=2048)
    assert r['partial_sum_bytes']==4*16*64*32*4
    assert r['partial_sum_read_write_bytes']==1048576
    assert r['workspace_bytes']==128*r['partial_sum_bytes']


def test_residual_std_workspace_is_larger_than_sum():
    from torchgwas.reduction_memory import column_std_workspace
    from torchgwas.setup_work import setup_memory
    from torchgwas.tensor_memory import eager_memory_plan
    # Actual allocator-history audit found this 512 MiB request at std().
    std=column_std_workspace(4096,4096,sm_count=132,max_threads_per_sm=2048)
    assert std['workspace_bytes']==512*1024**2
    assert std['accumulator_bytes']==16
    assert std['partial_sum_bytes']*128==std['workspace_bytes']
    fixed=setup_memory(4096,4096,8,sm_count=132,max_threads_per_sm=2048)
    old=setup_memory(4096,4096,8)
    assert fixed['residual_device_live_bytes_upper']==old['residual_device_live_bytes_upper']+512*1024**2+128
    profile=dict(torch_version='2.5.1',sm_count=132,max_threads_per_sm=2048,
        compute_capability=[9,0],cublas_workspace_config=None,cublas_handle_stream_pairs=2)
    plan=eager_memory_plan(4096,128,4096,8,4,profile)
    assert plan['device_bytes']>=fixed['residual_device_live_bytes_upper']+64*1024**2
    assert plan['device_bytes']>=637682176  # held-out allocated peak, not a model coefficient