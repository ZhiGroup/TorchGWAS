"""Workspace contracts and synchronous native storage, independent of durations."""
import copy
import pytest
from torchgwas.nonzero_memory import device_nonzero_workspace,device_selection_memory


def census_fixture(cells):
    context=dict(torch_version='2.5.1+cu124',cuda_runtime='12.4',device_uuid='test-device',
        compute_capability=[8,0],sm_count=108,library_sha256='test-binary')
    rows=[]
    for retained in sorted({0,1,cells}):
        private=[dict(action='alloc',addr=1,size=4),dict(action='alloc',addr=2,size=4095),dict(action='alloc',addr=3,size=1023)]
        coordinates=[dict(action='alloc',addr=4,size=16*retained)] if retained else []
        frees=[dict(action='free_requested',addr=e['addr'],size=e['size']) for e in [private[1],private[2],private[0]]]
        trace=private[:2]+coordinates+private[2:]+frees
        rows.append(dict(cells=cells,retained=retained,shape=[1,cells],dtype='torch.bool',
            output_ptr=4 if retained else 0,output_bytes=16*retained,output_allocation_bytes=((16*retained+511)//512)*512,
            extra_allocated_peak_bytes=5632+((16*retained+511)//512)*512,
            trace=trace,duration_fields_recorded=False))
    return dict(context=context,rows=rows,duration_fields_recorded=False)


def test_private_workspace_includes_overlap_and_excludes_coordinates():
    census=census_fixture(4093)
    work=device_nonzero_workspace(census,4093,census['context'])
    assert work['private_requested_bytes']==5122
    assert work['private_rounded_bytes']==5632


@pytest.mark.parametrize('fault',['binary','extent','time','peak','survivor','release'])
def test_unknown_or_inconsistent_allocation_evidence_is_refused(fault):
    census=census_fixture(4093);context=copy.deepcopy(census['context']);cells=4093
    if fault=='binary':context['library_sha256']='other'
    elif fault=='extent':cells=4092
    elif fault=='time':census['rows'][0]['trace'][0]['time_us']=123
    elif fault=='peak':census['rows'][0]['extra_allocated_peak_bytes']-=512
    elif fault=='survivor':census['rows'].pop()
    elif fault=='release':census['rows'][0]['trace'][-1]['addr']=999
    with pytest.raises(ValueError):device_nonzero_workspace(census,cells,context)


def test_one_cell_case_has_two_distinct_occupancy_controls():
    census=census_fixture(1)
    assert device_nonzero_workspace(census,1,census['context'])['private_rounded_bytes']==5632


def test_wide_selector_bound_depends_on_blocks_not_whole_trait_axis():
    census=census_fixture(17)
    one=device_selection_memory(40,1,17,census=census,context=census['context'],max_cells=17)
    many=device_selection_memory(40,100,17*7,census=census,context=census['context'],max_cells=17)
    assert many['selection_gpu_bytes']==one['selection_gpu_bytes']
    assert many['maximum_selection_cells']==17
    assert many['selected_payload_bytes_per_block']==340


def test_tail_extent_requires_its_own_installed_workspace_request():
    census=census_fixture(17)
    with pytest.raises(ValueError,match='controls'):
        device_selection_memory(40,1,18,census=census,context=census['context'],max_cells=17)
    census['rows']+=census_fixture(1)['rows']
    work=device_selection_memory(40,1,18,census=census,context=census['context'],max_cells=17)
    assert {row['cells'] for row in work['block_shapes']}=={1,17}


def test_memory_admission_uses_production_selection_geometry():
    census=census_fixture(8)
    census['rows']+=census_fixture(4)['rows']
    work=device_selection_memory(40,4,9,census=census,context=census['context'],max_cells=8)
    assert {(row['rows'],row['traits'],row['cells']) for row in work['block_shapes']}=={(4,2,8),(4,1,4)}
    assert work['maximum_selection_cells']==8
    assert 'selection_geometry.py' in work['source_sha256']


def test_old_geometry_workspace_census_cannot_admit_new_tail_shape():
    census=census_fixture(8)
    census['rows']+=census_fixture(1)['rows']
    with pytest.raises(ValueError,match='controls'):
        device_selection_memory(40,4,9,census=census,context=census['context'],max_cells=8)


def test_device_significance_stages_only_input_pins_and_retains_one_previous_dense_result():
    from torchgwas.pinned_work import pinned_scan_work
    from torchgwas.tensor_memory import eager_scan_memory
    pins=pinned_scan_work(128,32,17,3,reduction='device_significant')
    assert [row['name'] for row in pins['allocations']]==['genotype']
    assert pins['requested_bytes']==128*32*3
    two=eager_scan_memory(128,32,17,8,2,reduction='device_significant')
    four=eager_scan_memory(128,32,17,8,4,reduction='device_significant')
    assert two['retained_output_bytes']==0
    assert two['previous_dense_result_bytes']>8*32*17
    assert four['tensor_storage_budget']-two['tensor_storage_budget']==2*128*32
    assert two['critical_device_bytes']>=4*129
    assert not two['selection_memory_included']


def test_device_control_has_no_dense_result_ring_submission_or_future_resolution():
    from torchgwas.native_control_work import native_control_work
    work=native_control_work(reduction='device_significant')
    assert work['result_submit']=={'event_record':1,'executor_submit':1}
    assert work['finish']==work['resolve']=={}
    assert work['release']=={'event_synchronize_ready':1,'queue_put':1}


def device_candidate(input_path,**kwargs):
    from test_trait_tiling_model import candidate
    value=candidate(input_path,block_bytes=None,**kwargs)
    for tile in value['tiles']:
        tile['profile'].update(reduction='device_significant',result_ownership='owned')
        tile['data']['phenotype_complete']=True
    return value


def memory_profiles(candidate):
    cells=set()
    for tile in candidate['tiles']:
        b=tile['profile']['chunk_markers'];m=tile['data']['markers'];k=tile['data']['traits_analyzed']
        cells.add(min(b,m)*k)
        if m%b:cells.add((m%b)*k)
    census=census_fixture(next(iter(cells)))
    census['rows']=[row for cell in sorted(cells) for row in census_fixture(cell)['rows']]
    profile=dict(torch_version='2.5.1+cu124',sm_count=108,max_threads_per_sm=2048,
        compute_capability=[8,0],cublas_workspace_config=None,cublas_handle_stream_pairs=1,
        nonzero_workspace_census=census,nonzero_workspace_context=census['context'])
    return {device:copy.deepcopy(profile) for device in candidate['devices']}


def test_device_candidate_counts_global_queue_once_and_uses_only_input_pins(tmp_path):
    from test_trait_tiling_model import write_pgen
    import numpy as np
    from torchgwas.significant_device_model import significant_device_memory
    path=tmp_path/'input.pgen';write_pgen(path,np.arange(10*32,dtype=np.uint8).reshape(10,32)%4)
    candidate=device_candidate(path)
    profiles=memory_profiles(candidate)
    before=copy.deepcopy(candidate)
    first=significant_device_memory(candidate,device_memory_profiles=profiles)
    assert candidate==before
    candidate['output']['queue_depth']=3
    second=significant_device_memory(candidate,device_memory_profiles=profiles)
    assert second['host_bytes']-first['host_bytes']==2*28*8
    assert second['device_bytes']==first['device_bytes']
    assert first['maximum_selection_cells']==8
    assert first['pinned_cache_bytes_by_device']=={'cuda:0':256,'cuda:1':256}
    assert first['count_pinned_allocator_bytes']==8
    single=device_candidate(path,count=1)
    memory=significant_device_memory(single,device_memory_profiles=memory_profiles(single))
    assert memory['queued_result_bytes']==0


def test_device_memory_refuses_missing_full_or_tail_workspace_profile(tmp_path):
    from test_trait_tiling_model import write_pgen
    import numpy as np
    from torchgwas.significant_device_model import significant_device_memory
    path=tmp_path/'input.pgen';write_pgen(path,np.arange(10*32,dtype=np.uint8).reshape(10,32)%4)
    candidate=device_candidate(path,width=2,traits=2,count=1)
    profiles=memory_profiles(candidate)
    profiles['cuda:0']['nonzero_workspace_census']['rows']=[r for r in profiles['cuda:0']['nonzero_workspace_census']['rows'] if r['cells']!=4]
    with pytest.raises(ValueError,match='controls'):
        significant_device_memory(candidate,device_memory_profiles=profiles)


def test_cached_coordinate_block_is_not_misclassified_as_private_cub_workspace():
    census=census_fixture(4093)
    census['rows'][-1]['output_allocation_bytes']+=12288
    census['rows'][-1]['extra_allocated_peak_bytes']+=12288
    result=device_nonzero_workspace(census,4093,census['context'])
    assert result['private_rounded_bytes']==5632
    assert result['observed_coordinate_block_excess_bytes']==[0,0,12288]
