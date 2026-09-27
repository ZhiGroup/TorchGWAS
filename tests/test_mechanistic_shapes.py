from unittest.mock import patch
import pytest
from torchgwas.mechanistic_torch import torch_scan_work,torch_scan_runtime


def fixture(k=3,c=2):
    n,m=32,10
    data=dict(samples=n,markers=m,covariates=c,traits_analyzed=k,matching_sample_order=True,
        encoded=dict(samples=n,markers=m,record_form_counts={'0':m},ld_records_at_chunk_starts=0,record_payload_bytes=80))
    profile=dict(chunk_markers=8,depth=2,decode_workers=2,cpu_fraction=1.,cpu_available_cores=2.,
        read_bytes_per_second=1e8,shared_dram_bytes_per_second=1e9,h2d_bytes_per_second=1e8,d2h_bytes_per_second=1e8,
        process_units=dict(pgen_index_records=1e-8,finish_fixed_calls=1e-6,numpy_copy_bytes=1e-10),
        decode_units={'copy_packed_byte':1e-9,'expand4_int8':1e-9},executor_cpu_seconds=1e-6,stream_submission_cpu_seconds=1e-6,
        gpu_resources=dict(hbm_bytes_per_second=1e12,l2_bytes_per_second=2e12,fp32_flops_per_second=1e13,
            kernel_launch_seconds=1e-6,host_dispatch_cpu_seconds=2e-6,available_l2_bytes=1<<20,sm_count=100),host_primitives={},
        kernel_geometry=[dict(N=n,B=b,K=k,C=c,kernels=[]) for b in (8,2)])
    return data,profile


def component(*args,**kwargs):
    return dict(host_dispatch_cpu_seconds=1e-5,estimated_span_seconds=2e-5,kernel_service_seconds=2e-5,
        operations=[dict(host_submit_finish=1e-5,kernel_service_seconds=1e-5)],gemm={},unpriced_terms=[])


def test_exact_trait_covariate_shape_and_tail_are_counted():
    data,profile=fixture()
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(data,profile)
        assert len(work['blocks'])==2
        assert [b['h2d_bytes'] for b in work['blocks']]==[256,64]
        assert [b['d2h_bytes'] for b in work['blocks']]==[8*29,2*29]
        assert work['host_workspace']['packed_decoder_bytes_upper']==2*8*8
        assert work['host_workspace']['packed_decoder_allocation_calls_upper']==2
        result=torch_scan_runtime(data,profile)
        assert result['estimated_scan_seconds']>0
        assert not result['prediction_complete']


def test_geometry_does_not_reuse_k1_for_other_trait_widths():
    data,profile=fixture()
    for row in profile['kernel_geometry']:row['K']=1
    with pytest.raises(ValueError,match='geometry'):
        torch_scan_work(data,profile)


@pytest.mark.parametrize('setting',[{'borrow_results':True},{'result_ownership':'borrowed'}])
def test_owned_calculator_rejects_unmodeled_borrowed_lifetimes(setting):
    data,profile=fixture();profile.update(setting)
    with pytest.raises(ValueError,match='acknowledgement'):
        torch_scan_work(data,profile)


def test_borrowed_scan_has_no_owned_copy_free_or_duplicate_finish_control():
    from torchgwas.native_control_work import native_control_work
    data,profile=fixture(k=512)
    profile.update(result_ownership='borrowed',result_finish_service=dict(cpu_seconds=.01,serial_cpu_seconds=.008,
        baseline_copy_bytes=0,replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False),
        control_primitives={name:1e-6 for counts in native_control_work().values() for name in counts},
        owned_result_allocator={'unused':'deliberately not valid for owned results'})
    profile['process_units']['finish_fixed_calls']=99.
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(data,profile)
    assert work['result_ownership']=='borrowed'
    for block in work['blocks']:
        assert block['finish_seconds']==pytest.approx(.010001)
        assert sum(p['seconds'] for p in block['finish_operations'])==pytest.approx(.010001)
        assert block.get('discard_seconds',0.)==0.
    for result in work['owned_result_work'].values():
        assert result['work']['copy_bytes']==result['work']['allocation_calls']==0
        assert result['work']['borrowed_array_bytes']['beta']>0
        assert 'allocator' not in result


def test_range_validation_requires_matching_compiled_geometry():
    data,profile=fixture();profile['validate_range']=True
    with pytest.raises(ValueError,match='geometry'):
        torch_scan_work(data,profile)


def test_census_mismatch_and_integer_ld_keys_fail_closed():
    data,profile=fixture();data['encoded']['markers']=11
    with pytest.raises(ValueError,match='census dimensions'):
        torch_scan_work(data,profile)
    data,profile=fixture();data['encoded']['record_form_counts']={2:10}
    with pytest.raises(ValueError,match='LD replay'):
        torch_scan_work(data,profile)

def test_scan_retains_unpriced_geometry_terms():
    data,profile=fixture()
    def unresolved(*args,**kwargs):
        result=component();result['unpriced_terms']=['unknown synchronization']
        return result
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=unresolved):
        result=torch_scan_runtime(data,profile)
    assert 'unknown synchronization' in result['unpriced_terms']
    assert not result['prediction_complete']

def test_full_runtime_exposes_unpriced_scan_costs():
    from torchgwas.mechanistic_torch import torch_runtime
    data,profile=fixture(k=1,c=8)
    data.update(traits_in_file=1,tables={'pvar':{'bytes':200,'field_characters':100},'psam':{'bytes':100}})
    data['encoded']['file_bytes']=100
    profile.update(write_bytes_per_second=1e8,npy_itemsize=8,initial_index_parses=2,
        timing_boundary='environment-ready',pin_cpu_seconds_per_page=1e-6,
        pin_driver_seconds_per_page=1e-6,tiny_setup_cpu_seconds=1e-5,
        tiny_setup_non_cpu_seconds=1e-5,first_use_seconds=.01,design_first_use_seconds=.1,fsync_seconds=.001)
    profile['gpu_resources']['gpu_fraction']=1.
    for name in ['pvar_rows','psam_rows','phenotype_qc_cells','covariate_qc_cells',
                 'covariate_basis_work','bytearray_zero_bytes','variant_id_rows']:
        profile['process_units'][name]=1e-8
    def unresolved(*args,**kwargs):
        result=component();result['unpriced_terms']=['unknown synchronization']
        return result
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=unresolved):
        result=torch_runtime(data,profile)
    assert result['estimated_seconds']>0
    assert result['stage_seconds']['design_first_use']==.1
    assert not result['prediction_complete']
    assert 'unknown synchronization' in result['unpriced_terms']



def test_wide_runtime_requires_and_uses_blocked_setup_primitives():
    from torchgwas.mechanistic_torch import torch_runtime
    data,profile=fixture(k=3,c=8)
    with pytest.raises(ValueError,match='blocked setup primitives'):
        torch_runtime(data,profile)
    data.update(traits_in_file=3,tables={'pvar':{'bytes':200,'field_characters':100},'psam':{'bytes':100}})
    data['encoded']['file_bytes']=100
    profile.update(write_bytes_per_second=1e8,npy_itemsize=8,initial_index_parses=2,
        timing_boundary='environment-ready',pin_cpu_seconds_per_page=1e-6,
        pin_driver_seconds_per_page=1e-6,tiny_setup_cpu_seconds=1e-5,
        tiny_setup_non_cpu_seconds=1e-5,first_use_seconds=.01,design_first_use_seconds=.1,fsync_seconds=.001)
    profile['gpu_resources']['gpu_fraction']=1.
    profile['setup_primitives']={name:dict(reference_shape=[32,1,8],cpu_seconds=1e-5,non_cpu_seconds=1e-6)
        for name in ['residual_common','residual_block','design_common','design_block']}
    for name in ['pvar_rows','psam_rows','phenotype_qc_cells','covariate_qc_cells',
                 'covariate_basis_work','bytearray_zero_bytes','variant_id_rows']:
        profile['process_units'][name]=1e-8
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        result=torch_runtime(data,profile)
    assert result['source_work']['pinned_bytes']==2*(32*8+(8*3+5)*8)
    assert result['source_work']['setup']['traits']==3
    assert result['source_work']['setup']['h2d_bytes']==4*32*(2*3+2*8)
    assert result['schedule']['writer_payload_bytes']==12*10*3+4*10
    assert result['stage_seconds']['gpu_design']==result['source_work']['setup_service']['seconds']
    assert not result['prediction_complete']


def test_shard_scan_accounts_for_full_file_reader_index():
    data,profile=fixture()
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        baseline=torch_scan_work(data,profile)
        data['encoded']['file_markers']=110
        shard=torch_scan_work(data,profile)
    # Two reader workers each parse all110 index entries rather than10.
    assert shard['cpu_work_seconds']-baseline['cpu_work_seconds']==pytest.approx(2*100*profile['process_units']['pgen_index_records'])
    assert [b['decode_seconds'] for b in shard['blocks']]==[b['decode_seconds'] for b in baseline['blocks']]
    assert sum(b.get('reader_init_seconds',0.) for b in shard['blocks'])==pytest.approx(2*110*profile['process_units']['pgen_index_records'])


def test_owned_result_scenario_replaces_bulk_copy_price_without_double_counting():
    data,profile=fixture(k=64)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        original=torch_scan_work(data,profile)
        profile['owned_result_copy_scenario']=dict(resident_cpu_seconds_per_byte=2e-10,
            fresh_cpu_seconds_per_byte=8e-10,fresh_fraction=1.)
        fresh=torch_scan_work(data,profile)
    extra_bytes=sum(max(0,(8*64+5)*b-416) for b in [8,2])
    assert fresh['cpu_work_seconds']-original['cpu_work_seconds']==pytest.approx(extra_bytes*(8e-10-1e-10))
    assert fresh['owned_result_work'][8]['work']['allocation_calls']==4
    assert any('page-reuse' in term for term in fresh['unpriced_terms'])


def test_allocation_exposure_propagates_per_chunk_without_fitted_latency():
    from torchgwas.mechanistic_torch import torch_multigpu_scan_runtime
    data,profile=fixture(k=512)
    data['markers']=2049
    data['encoded'].update(markers=2049,record_form_counts={'0':2049},record_payload_bytes=8*2049)
    profile['chunk_markers']=2048
    profile['kernel_geometry']=[dict(N=32,B=b,K=512,C=2,kernels=[]) for b in (2048,1)]
    policy=dict(numpy_version='2.2.6',numpy_madvise_hugepage=True,linux_thp_enabled='madvise')
    profile['owned_result_allocation_policy']=policy
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        exposed=torch_scan_runtime(data,profile)
        multi=torch_multigpu_scan_runtime([dict(device='cuda:0',data=data,profile=profile)],dict(cpu=2.,dram=1e9,input=1e8))
        policy['numpy_madvise_hugepage']=False
        disabled=torch_scan_runtime(data,profile)
    allocation=exposed['source_work']['owned_result_work']
    assert allocation[2048]['allocation_policy']['advice_compaction_exposure'] is True
    assert allocation[1]['allocation_policy']['advice_compaction_exposure'] is False
    assert any('hugepage' in term for term in exposed['unpriced_terms'])
    assert any('hugepage' in term for term in multi['unpriced_terms'])
    assert multi['shards'][0]['owned_result_work'][2048]['allocation_policy']['advice_compaction_exposure'] is True
    assert not any('hugepage' in term for term in disabled['unpriced_terms'])
    assert exposed['estimated_scan_seconds']==disabled['estimated_scan_seconds']
    assert not disabled['prediction_complete']

def test_input_read_cpu_is_charged_once_and_propagates_to_scan_graph():
    data,profile=fixture()
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        baseline=torch_scan_work(data,profile)
        profile['input_read_cpu_prices']=dict(cpu_seconds_per_byte=1e-6,cpu_seconds_per_call=2e-6)
        with_read=torch_scan_work(data,profile)
    expected=80e-6+2*2e-6
    assert with_read['cpu_work_seconds']-baseline['cpu_work_seconds']==pytest.approx(expected)
    assert sum(block['decode_read_seconds']*block['read_resources']['cpu'] for block in with_read['blocks'])==pytest.approx(expected)
    assert sum(block['decode_seconds']*block['decode_resources']['cpu'] for block in with_read['blocks'])==pytest.approx(sum(block['decode_seconds']*block['decode_resources']['cpu'] for block in baseline['blocks']))
    assert 'input-read syscall/copy CPU service' in baseline['unpriced_terms']
    assert 'input-read syscall/copy CPU service' not in with_read['unpriced_terms']


def test_exact_chunk_work_changes_local_service_without_changing_totals(tmp_path):
    from test_pgen_native_reader import write_pgen_mixed
    from torchgwas.pgen_work_census import census
    from torchgwas.decoder_work import decoder_work
    import numpy as np
    data, profile = fixture()
    path = tmp_path/'unequal.pgen'
    values = np.zeros((10, 32), dtype=np.uint8)
    values[:8] = np.arange(32, dtype=np.uint8)%4
    write_pgen_mixed(path, values, [0]*8+[4]*2)
    data['encoded'] = census(path, 8, include_chunks=True)
    units = decoder_work(data['encoded'], 'torch_native_int8')['source_units']
    profile['decode_units'] = {name:1e-6 for name in units}
    profile['decode_units']['copy_packed_byte'] = 4e-6
    with patch('torchgwas.mechanistic_torch.tensor_stage_service', side_effect=component):
        exact = torch_scan_work(data, profile)
        runtime = torch_scan_runtime(data, profile)
        del data['encoded']['chunks']
        uniform = torch_scan_work(data, profile)
    assert exact['encoded_work_distribution'] == 'exact_encoded_chunks'
    assert runtime['source_work']['encoded_work_distribution'] == 'exact_encoded_chunks'
    assert uniform['encoded_work_distribution'] == 'uniform_aggregate_scenario'
    assert exact['cpu_work_seconds'] == pytest.approx(uniform['cpu_work_seconds'])
    assert sum(b['decode_read_seconds'] for b in exact['blocks']) == pytest.approx(sum(b['decode_read_seconds'] for b in uniform['blocks']))
    assert exact['blocks'][0]['decode_seconds'] > uniform['blocks'][0]['decode_seconds']
    assert exact['blocks'][1]['decode_seconds'] < uniform['blocks'][1]['decode_seconds']
    assert exact['blocks'][0]['decode_read_seconds'] > uniform['blocks'][0]['decode_read_seconds']
    assert [b['h2d_bytes'] for b in exact['blocks']] == [b['h2d_bytes'] for b in uniform['blocks']]
