"""Shared-ring bounds must cover mixed lifetimes and shifted LD restarts."""
import copy
from unittest.mock import patch
import numpy as np
import pytest
from torchgwas.adaptive_chunks import AlignedChunkSizeControl,aligned_chunk_shapes
from torchgwas.adaptive_candidate import adaptive_tensor_memory,adaptive_candidate_memory
from torchgwas.decoder_work import native_read_layout,native_reader_workspace
from torchgwas.pgen_work_census import census
from torchgwas.tensor_memory import eager_memory_plan
from test_jagwas_actual_candidate import actual_candidate
from test_pgen_native_reader import write_pgen_mixed


def memory_profile():
    return dict(torch_version='2.5.1+cu124',sm_count=108,max_threads_per_sm=2048,
        compute_capability=[8,0],cublas_workspace_config=None,cublas_handle_stream_pairs=2)


@pytest.mark.parametrize('sizes,length',[
    ([4,8,16],137),([4,12,20],39),([128,256,512],1025),([128,256,512],1),
    ([4,8,16],16),([4,8,16],32)])
def test_every_reachable_transition_keeps_known_shapes_and_exact_extent(sizes,length):
    # Enumerate reachable states, not just a handpicked transition sequence.
    control=AlignedChunkSizeControl(sizes,initial=sizes[0]);pending=[0];seen=set()
    shapes=set(aligned_chunk_shapes(sizes,length))
    while pending:
        start=pending.pop()
        if start in seen:continue
        seen.add(start)
        if start==length:continue
        assert start<length
        for size in sizes:
            control.set_size(size)
            count=control(start+7,length+7,max(sizes))
            assert count in shapes and count<=length-start
            assert start%sizes[0]==0
            pending.append(start+count)
    assert length in seen
    # An oversize final request splits into known sizes instead of shape 9.
    if sizes==[4,8,16] and length==137:
        control.set_size(16);assert control(128,137,16)==8


def test_nondivisible_grid_refused_without_changing_private_general_control():
    with pytest.raises(ValueError,match='multiples'):AlignedChunkSizeControl([3,17],initial=3)


@pytest.mark.parametrize('reduction',[None,'jagwas'])
def test_retained_capacity_is_not_the_current_small_chunk(reduction):
    kwargs=dict(device_profile=memory_profile(),reduction=reduction)
    bound=adaptive_tensor_memory(129,7,2,4,chunk_sizes=[4,8,16],markers=137,**kwargs)
    small=eager_memory_plan(129,4,7,2,4,**kwargs)
    assert bound['tensor_allocations']['staging_bytes']==4*129*16
    assert bound['device_bytes']>=small['device_bytes']
    assert bound['pinned']['requested_bytes']>small['scan']['staging_bytes']
    assert bound['work_shapes']==[1,4,8,16] and not bound['prediction_complete']


def test_mixed_previous_output_and_current_workspace_need_componentwise_bound():
    # A smaller work shape can use a different library/source branch. A maximum
    # of the two separate run peaks misses its workspace plus a previous result.
    def plan(n,b,k,c,depth,profile,**kwargs):
        scan=dict(resident_bytes=10,distinct_temporary_bytes=100 if b==4 else 20,
            retained_output_bytes=10 if b==4 else 90,previous_genotype_bytes=4*b)
        return dict(scan=scan,setup_bytes=0,cublas={'total_bytes':0},unresolved_memory_terms=[])
    with patch('torchgwas.adaptive_candidate.eager_memory_plan',side_effect=plan):
        bound=adaptive_tensor_memory(32,2,0,2,chunk_sizes=[4,8],markers=32,device_profile={})
    separate=[10+t+r+4*b+2*32*b for b,t,r in [(4,100,10),(8,20,90)]]
    assert bound['device_bytes']>max(separate)


def joint_fixture(tmp_path):
    path=tmp_path/'shifted.pgen';n,m=2049,1025
    calls=np.tile((np.arange(n)%3).astype(np.uint8),(m,1))
    for row in range(m):calls[row,row%n]=(calls[row,row%n]+1)%3
    forms=[0 if i%512==0 else 2 for i in range(m)]
    write_pgen_mixed(path,calls,forms)
    candidate=actual_candidate(path,512,1)
    return path,candidate,census(path,128,include_chunks=True)


def test_shifted_ld_prefix_and_grow_only_reader_memory_cover_all_possible_starts(tmp_path):
    path,candidate,fine=joint_fixture(tmp_path);before=copy.deepcopy((candidate,fine))
    bound=adaptive_candidate_memory(candidate,chunk_sizes=[128,256,512],source_census=fine,
        reduction='jagwas',device_memory_profiles={'cuda:0':memory_profile()})
    row=bound['tiles'][0]
    old=native_reader_workspace(candidate['tiles'][0]['data']['encoded']['chunks'])
    assert row['reader_payload_upper_bytes']>old[0]
    assert bound['decoder_extra_bytes_by_device']['cuda:0']>0
    assert not bound['missing_geometry'] and row['work_shapes']==[1,128,256,512]
    for start in range(0,1025,128):
        for size in [128,256,512]:
            end=min(start+size,1025)
            encoded=census(path,size,variant_range=(start,end),include_chunks=True)
            layout=native_read_layout(encoded['chunks'][0])
            assert layout['read_bytes']<=row['reader_payload_upper_bytes']
            assert layout['extra_workspace_bytes']<=row['reader_replay_scratch_upper_bytes']
    assert (candidate,fine)==before
    assert bound['host_bytes']>bound['base_fixed_capacity_memory']['host_bytes']


@pytest.mark.parametrize('fault',['capacity','source','gap','budget','shapes','devices'])
def test_adaptive_preflight_rejects_unbound_or_unbounded_input_before_tensor_work(tmp_path,fault):
    _,candidate,fine=joint_fixture(tmp_path)
    options=dict(chunk_sizes=[128,256,512],source_census=fine,reduction='jagwas',
        device_memory_profiles={'cuda:0':memory_profile()})
    if fault=='capacity':options['chunk_sizes']=[128,256]
    if fault=='source':fine['path']='different.pgen'
    if fault=='gap':fine['chunks'][1]['variant_range'][0]+=1
    if fault=='budget':options['max_census_chunks']=1
    if fault=='shapes':options['max_kernel_shapes']=1
    if fault=='devices':options['device_memory_profiles']={}
    with patch('torchgwas.adaptive_candidate.adaptive_tensor_memory',side_effect=AssertionError('tensor work started')):
        with pytest.raises(ValueError):adaptive_candidate_memory(candidate,**options)


def test_missing_smaller_geometry_is_reported_without_claiming_admission(tmp_path):
    _,candidate,fine=joint_fixture(tmp_path)
    candidate['tiles'][0]['profile']['joint_kernel_geometry']=[
        r for r in candidate['tiles'][0]['profile']['joint_kernel_geometry'] if r['B']!=256]
    bound=adaptive_candidate_memory(candidate,chunk_sizes=[128,256,512],source_census=fine,
        reduction='jagwas',device_memory_profiles={'cuda:0':memory_profile()})
    assert bound['missing_geometry']==[dict(device='cuda:0',shape=[2049,256,512,2],reduction='jagwas')]


@pytest.mark.parametrize('reduction',[None,'significant'])
def test_wide_output_admission_keeps_original_writer_and_selection_reserves(tmp_path,reduction):
    from test_trait_tiling_model import candidate as trait_candidate
    from test_pgen_native_reader import write_pgen
    path=tmp_path/'wide.pgen';write_pgen(path,np.arange(10*32,dtype=np.uint8).reshape(10,32)%4)
    candidate=trait_candidate(path,width=2,count=2,block_bytes=None)
    before=copy.deepcopy(candidate);profiles={d:memory_profile() for d in candidate['devices']}
    bound=adaptive_candidate_memory(candidate,chunk_sizes=[2,4],source_census=census(path,2,include_chunks=True),
        reduction=reduction,device_memory_profiles=profiles,host_reserve_bytes=1234,device_reserve_bytes=5678)
    base=bound['base_fixed_capacity_memory']
    assert bound['host_bytes']>=base['host_bytes'] and candidate==before
    assert all(bound['device_bytes'][d]>=base['device_bytes'][d] for d in profiles)
    if reduction=='significant':
        assert base['queued_result_bytes']>0 and base['selection_bytes_by_device']
    assert all(row['tensor_memory']['capacity']==4 for row in bound['tiles'])
