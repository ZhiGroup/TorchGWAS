"""The calculator counts the actual JAGWAS operations, including FP64 work."""
from unittest.mock import patch
import pytest
from torchgwas.reduction_tensor_work import jagwas_tensor_work
from torchgwas.reduce import JagwasReduction


@pytest.mark.parametrize('traits',[1,7,128])
def test_joint_factor_and_projection_types_shapes_and_aliases(traits):
    n,b=257,13
    preparation=jagwas_tensor_work(n,b,traits,phase='prepare')
    projection=jagwas_tensor_work(n,b,traits,phase='reduce')
    assert preparation['initial_storages'][0]['bytes']==4*n*traits
    assert preparation['result_storages'][0]['bytes']==8*traits*traits
    assert preparation['result_storages'][0]['dtype']=='torch.float64'
    assert sum(s['bytes'] for s in projection['result_storages'])==17*b
    assert projection['result_storages'][1]['dtype']=='torch.float32'
    # The factor's Gram is FP64 (one sample block of the FP32 panel at n=257).
    for ledger,flops,dtype in [(preparation,2*n*traits*traits,'torch.float64'),
                               (projection,2*b*traits*traits,'torch.float64')]:
        products=[row for row in ledger['steps'] if 'matmul_flops' in row]
        assert len(products)==1
        assert products[0]['matmul_flops']==flops and products[0]['matmul_dtype']==dtype
        assert ledger['distinct_temporary_bytes']>0
        assert ledger['logical_bytes']==sum(s['logical_bytes'] for s in ledger['steps'])
        assert all(s['logical_bytes']==0 for s in ledger['steps'] if s['alias_only'])
        assert ledger['unpriced_terms']
    # The default eigen factor: eigh of R, then QR of the scaled kept directions.
    assert any('linalg_eigh' in row['op'] for row in preparation['steps'])
    assert any('linalg_qr' in row['op'] for row in preparation['steps'])
    rounding=jagwas_tensor_work(n,b,traits,phase='prepare',method='rounding')
    assert any('cholesky' in row['op'] for row in rounding['steps'])
    assert any('solve_triangular' in row['op'] for row in rounding['steps'])
    # Status and df are passed through rather than allocated by the reducer.
    for result,initial in zip(projection['result_storages'][-2:],projection['initial_storages'][2:4]):
        assert result['storage']==initial['storage']


def test_a_real_source_operation_change_changes_the_work_ledger():
    from torchgwas.jagwas_projection import JagwasReduction
    traced=JagwasReduction
    expected=jagwas_tensor_work(257,13,128,phase='reduce')
    original=traced.reduce
    def extra(self,beta,t,status,df,width):
        retained=t.clone()
        result=original(self,beta,t,status,df,width)
        assert retained.shape==t.shape
        return result
    with patch.object(traced,'reduce',extra):
        changed=jagwas_tensor_work(257,13,128,phase='reduce')
    assert changed['distinct_temporary_bytes']>expected['distinct_temporary_bytes']
    assert changed['logical_bytes']==expected['logical_bytes']+2*4*13*128
    assert len(changed['steps'])==len(expected['steps'])+1


def test_voxel_scale_requests_allocate_only_meta_storage():
    ledger=jagwas_tensor_work(22250,512,2075298,phase='prepare')
    assert ledger['result_storages'][0]['bytes']==8*2075298**2
    assert ledger['distinct_temporary_bytes']>28*2075298**2
    # Beyond the factorization, the projection's row blocks are views only.
    assert sum(not step['alias_only'] for step in ledger['steps'])<30


@pytest.mark.parametrize('change',[dict(samples=0),dict(markers=True),dict(traits=-1),dict(phase='scan')])
def test_invalid_requests_are_refused(change):
    args=dict(samples=32,markers=4,traits=3,phase='prepare');args.update(change)
    with pytest.raises(ValueError):jagwas_tensor_work(**args)


@pytest.mark.parametrize('phase',['prepare','reduce'])
def test_fp64_inputs_change_real_work_without_changing_factor_precision(phase):
    ledger=jagwas_tensor_work(257,13,7,phase=phase,compute_dtype='float64')
    assert ledger['compute_dtype']=='float64'
    assert ledger['initial_storages'][0]['dtype']=='torch.float64'
    assert ledger['initial_storages'][0]['bytes']==8*(257 if phase=='prepare' else 13)*7
    products=[row for row in ledger['steps'] if 'matmul_flops' in row]
    assert len(products)==1 and products[0]['matmul_dtype']=='torch.float64'
    assert ledger['result_storages'][0]['dtype']=='torch.float64'
    if phase=='reduce':
        assert ledger['result_storages'][1]['bytes']==8*13
        assert ledger['result_storages'][-1]['dtype']=='torch.float64'  # df passes through in the compute dtype


def test_unknown_compute_precision_is_refused():
    with pytest.raises(ValueError,match='compute dtype'):
        jagwas_tensor_work(32,4,3,phase='prepare',compute_dtype='float16')


@pytest.mark.parametrize('compute_dtype,element_bytes',[('float32',4),('float64',8)])
def test_factor_capacity_counts_simultaneously_live_matrices(compute_dtype,element_bytes):
    from torchgwas.reduction_tensor_work import jagwas_factor_memory_floor
    n,k=257,128
    work=jagwas_factor_memory_floor(n,k,compute_dtype=compute_dtype)
    # Eigen (default): eigenvectors, values and their scaled copy beside the
    # FP64 correlation exceed the Gram stage (correlation and one FP64 block)
    # at this shape.
    assert work['method']=='eigen' and work['explicit_live_bytes']==element_bytes*n*k+24*k*k+8*k
    assert work['arrays']['gram_block']==(8*n*k if compute_dtype=='float32' else 0)
    wide=jagwas_factor_memory_floor(8192,64,compute_dtype=compute_dtype)
    assert wide['explicit_live_bytes']==element_bytes*8192*64+max(24*64*64+8*64,8*64*64+wide['arrays']['gram_block'])
    # Rounding cutoff: the solve stage (correlation, factor, identity, inverse).
    rounding=jagwas_factor_memory_floor(n,k,compute_dtype=compute_dtype,method='rounding')
    assert rounding['explicit_live_bytes']==element_bytes*n*k+32*k*k
    if compute_dtype=='float32':assert wide['arrays']['gram_block']==8*4096*64>24*64*64
    assert work['persistent_factor_bytes']==8*k*k
    trace=jagwas_tensor_work(n,3,k,phase='prepare',compute_dtype=compute_dtype)
    assert work['explicit_live_bytes']<=element_bytes*n*k+trace['distinct_temporary_bytes']


def test_factor_capacity_probes_each_device_and_counts_reusable_cache():
    from torchgwas.reduction_tensor_work import require_jagwas_factor_capacity,jagwas_factor_memory_floor
    work=jagwas_factor_memory_floor(257,128);needed=work['explicit_live_bytes']
    with patch('torch.cuda.mem_get_info',side_effect=[(needed,needed*2),(needed-1,needed*2)]) as probe, \
         patch('torch.cuda.memory_reserved',return_value=4096), \
         patch('torch.cuda.memory_allocated',return_value=4096):
        with pytest.raises(ValueError,match='cuda:3'):
            require_jagwas_factor_capacity(257,128,['cpu','cuda:1','cuda:3'])
        assert [str(call.args[0]) for call in probe.call_args_list]==['cuda:1','cuda:3']
    with patch('torch.cuda.mem_get_info',return_value=(needed-1,needed*2)), \
         patch('torch.cuda.memory_reserved',return_value=4096), \
         patch('torch.cuda.memory_allocated',return_value=4095):
        assert require_jagwas_factor_capacity(257,128,['cuda:1'])==work


@pytest.mark.parametrize('retained',[0,1,13,4096])
def test_joint_writer_census_matches_real_npz_parts(tmp_path,retained):
    import json
    import numpy as np
    from torchgwas.reduced_output_work import jagwas_writer_work
    from torchgwas.sumstats_indexed import write_indexed_sumstats
    b=4096;t=np.full((b,1),np.nan,np.float32);t[:retained,0]=np.arange(retained)+.5
    beta=np.full_like(t,np.nan);index=np.zeros((b,1),np.int32)
    chunks=iter([(0,b,beta,t,None,index)])
    count,summary=write_indexed_sumstats(tmp_path,[str(i) for i in range(b)],['a','b'],129,chunks,
        kind='jagwas',df=125,chi2_df=2,fsync=False)
    work=jagwas_writer_work(b,retained,fsync=False)
    assert count==retained and work['part_fsync_calls']==0
    manifest=json.loads((tmp_path/'manifest.json').read_text())
    assert len(manifest['parts'])==int(retained>0)
    assert work['part']['file_bytes']==sum(p.stat().st_size for p in tmp_path.glob('part_*.npz'))
    assert work['indexed_array_payload_bytes']==16*retained
    if retained:
        with np.load(tmp_path/manifest['parts'][0]['file']) as part:
            assert [(name,part[name].dtype.str) for name in part.files]==[('variant_index','<i8'),('chi2','<f8')]
            np.testing.assert_array_equal(part['variant_index'],np.arange(retained))
            np.testing.assert_array_equal(part['chi2'],t[:retained,0].astype(np.float64))


@pytest.mark.parametrize('retained',[0,17])
def test_joint_writer_counts_transfer_and_queue_arrays_separately(retained):
    from torchgwas.reduced_output_work import jagwas_writer_work
    work=jagwas_writer_work(17,retained)
    assert work['native_result_d2h_bytes']==17*17
    assert work['owned_chunk_array_bytes']==12*17
    assert work['part_fsync_calls']==int(retained>0)
    assert work['writer_array_bytes_upper']>=25*17+24*retained


@pytest.mark.parametrize('rows',[-1,True,134217728])
def test_joint_writer_refuses_invalid_or_unbounded_zip_members(rows):
    from torchgwas.reduced_output_work import jagwas_indexed_part_work
    with pytest.raises(ValueError):jagwas_indexed_part_work(rows)


@pytest.mark.parametrize('chunk,depth',[(1,1),(511,3),(8192,4)])
def test_joint_pinned_requests_match_actual_reducer_buffers(chunk,depth):
    from torchgwas.pinned_work import pinned_scan_work
    work=pinned_scan_work(257,chunk,128,depth,reduction='jagwas')
    buffers=JagwasReduction().host_buffers(chunk,1,pin_memory=False)
    requests={r['name']:r for r in work['allocations']}
    assert work['allocation_count']==6*depth
    for name,buffer in zip(['beta','tstat','trait_index','flags','df'],buffers):
        request=requests[name]
        assert request['requested_bytes']==buffer.numel()*buffer.element_size()
        assert request['count']==depth
        assert request['allocator_bytes']>=request['requested_bytes']
    assert work['requested_bytes']==depth*chunk*(257+17)
    # Full-K compute does not imply a full-K result ring.
    assert pinned_scan_work(257,chunk,128,depth)['requested_bytes']>work['requested_bytes']


def test_joint_scan_memory_keeps_dense_work_but_retains_narrow_outputs():
    from torchgwas.tensor_memory import eager_scan_memory
    n,b,k,c=257,511,128,2
    dense=eager_scan_memory(n,b,k,c,3)
    joint=eager_scan_memory(n,b,k,c,3,reduction='jagwas')
    deeper=eager_scan_memory(n,b,k,c,4,reduction='jagwas')
    narrow=sum(((size+511)//512)*512 for size in [4*b,4*b,4*b,b,4*b])
    assert joint['retained_output_bytes']==2*narrow
    assert joint['retained_output_bytes']<dense['retained_output_bytes']
    assert joint['distinct_temporary_bytes']>dense['distinct_temporary_bytes']
    assert joint['resident_bytes']==dense['resident_bytes']+8*k*k
    assert joint['previous_genotype_bytes']==dense['previous_genotype_bytes']
    assert joint['previous_reduced_result_bytes']==narrow
    assert deeper['tensor_storage_budget']-joint['tensor_storage_budget']==n*b+narrow
    assert 'jagwas_projection.py' in joint['source_sha256']


def test_joint_memory_plan_includes_preparation_and_persistent_factor(monkeypatch):
    # The census is xpotrf's: the rounding cutoff's Cholesky factor.
    monkeypatch.setenv('TORCHGWAS_JAGWAS_RCOND', '0')
    from torchgwas.tensor_memory import eager_memory_plan
    profile=dict(torch_version='2.5.1',sm_count=108,max_threads_per_sm=2048,
        compute_capability=[8,0],cublas_workspace_config=None,cublas_handle_stream_pairs=2)
    n,b,k,c=4097,512,2048,27
    plan=eager_memory_plan(n,b,k,c,3,profile,reduction='jagwas')
    assert plan['prediction_complete'] is False
    assert plan['setup_bytes']>=4*n*k+plan['factor_setup']['distinct_temporary_bytes']+plan['cublas']['total_bytes']
    assert plan['setup_bytes']>=plan['setup']['design_device_live_bytes_upper']+8*k*k+plan['cublas']['total_bytes']
    assert plan['device_bytes']>=plan['scan']['tensor_storage_budget']+plan['cublas']['total_bytes']
    assert any('Cholesky' in term for term in plan['unresolved_memory_terms'])
    with pytest.raises(ValueError,match='whole phenotype'):
        eager_memory_plan(n,b,k,c,3,profile,reduction='jagwas',preprocessing_traits=k+1)


@pytest.mark.parametrize('function',['pinned','device'])
def test_unknown_reduction_memory_layout_cannot_fall_back_to_dense(function):
    from torchgwas.pinned_work import pinned_scan_work
    from torchgwas.tensor_memory import eager_scan_memory
    with pytest.raises(ValueError,match='Unsupported'):
        if function=='pinned':pinned_scan_work(257,13,7,3,reduction='topk')
        else:eager_scan_memory(257,13,7,2,3,reduction='topk')
