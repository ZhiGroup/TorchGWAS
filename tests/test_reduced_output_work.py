import numpy as np
import pytest
from torchgwas.reduced_output_work import significant_output_work,jagwas_output_work
from torchgwas.sumstats_indexed import write_indexed_sumstats,open_indexed_sumstats


@pytest.mark.parametrize('backend',['host','device'])
@pytest.mark.parametrize('store_beta',[True,False])
def test_source_payload_matches_real_indexed_parts(tmp_path,backend,store_beta):
    args=(40,7,11,4)
    shape=significant_output_work(*args,backend=backend,max_selection_cells=5)
    counts=[0 if i%3==0 else b['cells'] for i,b in enumerate(shape['blocks'])]
    work=significant_output_work(*args,backend=backend,max_selection_cells=5,retained_per_block=counts,store_beta=store_beta)
    chunks=[]
    for b,count in zip(work['blocks'],counts):
        vi,ti=np.indices((b['variant_range'][1]-b['variant_range'][0],b['trait_range'][1]-b['trait_range'][0]))
        chunks.append((*b['variant_range'],vi.ravel()[:count]+b['variant_range'][0],ti.ravel()[:count]+b['trait_range'][0],
            np.ones(count,np.float32),np.full(count,3.,np.float32),np.full(count,38.,np.float32)))
    rows,_=write_indexed_sumstats(tmp_path,list(map(str,range(7))),list(map(str,range(11))),40,chunks,
        kind='significant',df=38,store_beta=store_beta)
    manifest,parts=open_indexed_sumstats(tmp_path);parts=list(parts)
    assert rows==work['retained_pairs'] and len(parts)==work['nonempty_parts']
    assert sum(a.nbytes for p in parts for a in p.values())==work['indexed_array_payload_bytes']
    assert sum(b['cells'] for b in work['blocks'])==7*11
    if backend=='device':assert work['maximum_selection_block_cells']<=5
    assert work['selection_count_d2h_bytes'] == (4*work['selection_blocks']
                                                 if backend=='device' else 0)


def test_unknown_sparse_counts_remain_bounds_and_distribution_changes_fsync():
    unknown=significant_output_work(40,8,11,4,backend='device',max_selection_cells=11)
    assert unknown['retained_pairs'] is None and unknown['result_payload_d2h_bytes'] is None
    assert unknown['retained_pair_bounds']==[0,88]
    a=significant_output_work(40,8,11,4,backend='device',max_selection_cells=11,retained_per_block=[8,0,0,0,0,0,0,0])
    b=significant_output_work(40,8,11,4,backend='device',max_selection_cells=11,retained_per_block=[1]*8)
    assert a['indexed_array_payload_bytes']==b['indexed_array_payload_bytes']
    assert a['part_fsync_calls']==1 and b['part_fsync_calls']==8
    with pytest.raises(ValueError):significant_output_work(40,8,11,4,backend='device',retained_per_block=[89])


def test_jagwas_counts_match_narrow_native_ring_and_joint_rank_limit():
    w=jagwas_output_work(40,100,11,7,covariate_rank=2,retained_variants=91)
    assert w['native_result_payload_d2h_bytes']==100*(4+4+4+1+4)
    assert w['indexed_array_payload_bytes']==91*(8+8)
    assert w['persistent_inverse_cholesky_bytes']==8*11**2
    assert not w['trait_separable'] and w['full_rank_possible']
    assert not jagwas_output_work(40,100,38,7,covariate_rank=2)['full_rank_possible']

def test_voxel_scale_aggregate_work_does_not_expand_millions_of_events():
    args=(22250,8931083,2075298,512)
    with pytest.raises(ValueError,match='max_blocks'):
        significant_output_work(*args,backend='device',max_blocks=10)
    work=significant_output_work(*args,backend='device',include_blocks=False)
    # One selection block per chunk: 512 x 2,075,298 cells is below CUDA
    # nonzero's INT_MAX limit.
    assert work['blocks'] is None and work['selection_blocks']==-(-8931083//512)
    assert work['selection_cells']==8931083*2075298
    assert work['maximum_selection_block_cells']==512*2075298
    bounded=significant_output_work(*args,backend='device',include_blocks=False,max_selection_cells=1<<20)
    assert bounded['selection_blocks']>1000000 and bounded['maximum_selection_block_cells']<=1<<20

def test_inclusive_threshold_one_skips_inverse_cdf_service():
    for backend in ['host','device']:
        assert significant_output_work(40,7,11,4,backend=backend,significance_threshold=1.)['critical_value_evaluations']==0
    with pytest.raises(ValueError):significant_output_work(40,7,11,4,backend='host',significance_threshold=0.)


def test_significant_layout_caps_active_devices_and_global_readers():
    from torchgwas.reduced_output_work import significant_execution_layout
    layout = significant_execution_layout(17, 6, ['cuda:1', 'cuda:2'], 3, queue_depth=1)
    assert layout['tiles'] == 3 and layout['readers_per_device'] == [2, 1]
    assert layout['queue_depth'] == 1 and layout['producer_pending_slots'] == 2
    one = significant_execution_layout(3, 6, ['cuda:2', 'cuda:1'], 1)
    assert one['devices'] == ['cuda:2'] and one['queue_depth'] == 0
    with pytest.raises(ValueError, match='at least one'):
        significant_execution_layout(17, 6, ['cuda:1', 'cuda:2'], 1)
    with pytest.raises(ValueError, match='Unique'):
        significant_execution_layout(17, 6, ['cuda:1', 'cuda:1'], 3)
    with pytest.raises(ValueError, match='queue_depth'):
        significant_execution_layout(17, 6, ['cuda:1', 'cuda:2'], 3, queue_depth=0)

def test_api_shared_critical_table_removes_repeated_inverse_cdf_work():
    for backend in ['host', 'device']:
        work = significant_output_work(40000, 1000000, 8192, 512, backend=backend,
            critical_table_reused=True, include_blocks=False)
        assert work['critical_value_evaluations'] == 0
        assert work['critical_lookup_entries'] == 1000000
