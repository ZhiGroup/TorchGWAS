"""Untimed installed launches and independently supplied resource arithmetic."""
import copy
import json
from pathlib import Path
import pytest
from torchgwas.device_significance_work import device_significant_tensor_work
from torchgwas.device_selection_gpu_work import selection_gpu_census,selection_gpu_work,selection_gpu_service
from torchgwas.tensor_service import DeviceService


def fixture(index):
    census=json.loads((Path(__file__).parent/'fixtures/device_selection_geometry_cuda0.json').read_text())
    row=census['rows'][index]
    work=device_significant_tensor_work(row['N'],row['B'],row['K'],[b['retained'] for b in row['blocks']],max_cells=row['max_cells'])
    kernels=selection_gpu_census(work,census,census['context'])
    return census,work,kernels


def resources():
    return DeviceService(hbm_bytes_per_second=1e12,l2_bytes_per_second=2e12,fp32_flops_per_second=1e13,
        kernel_launch_seconds=1e-6,host_dispatch_cpu_seconds=0.,available_l2_bytes=40<<20,sm_count=108)


@pytest.mark.parametrize('index',range(12))
def test_all_real_a100_selector_launches_are_assigned_once(index):
    census,work,kernels=fixture(index)
    ledger=selection_gpu_work(work,kernels,compute_capability=census['context']['compute_capability'])
    used=[i for indices in ledger['operations'].values() for i in indices]
    used += [i for group in ledger['nonzero_groups'] for indices in group.values() for i in indices]
    assert sorted(used)==list(range(len(kernels)))
    for block,phases in zip(work['blocks'],ledger['nonzero_groups']):
        assert len(phases['count'])==(2 if block['cells']>4096 else 1)
        assert len(phases['flagged_select'])==2
        assert ('coordinate_scatter' in phases)==bool(block['retained'])
    service=selection_gpu_service(ledger,resources(),traffic_mode='logical_hbm',int64_divmod_per_second=1e10)
    assert sum(row['seconds'] for row in service['operation_services'].values())+sum(row['seconds'] for group in service['nonzero_services'] for row in group.values())==pytest.approx(service['gpu_service_seconds'])
    assert service['gpu_service_seconds']>=len(kernels)*1e-6


def test_empty_results_keep_selection_and_both_count_passes():
    _,work,kernels=fixture(8)
    ledger=selection_gpu_work(work,kernels,compute_capability=[8,0])
    assert ledger['kernel_count']==28
    assert [len(g['count']) for g in ledger['nonzero_groups']]==[2,1]
    assert not any(r.get('int64_divmod_pairs') for r in ledger['kernels'])
    first=ledger['kernels'][ledger['nonzero_groups'][0]['flagged_select'][0]]
    assert first['write_bytes']==8*(455+32)+4
    selection_gpu_service(ledger,resources(),traffic_mode='logical_l2')


def test_nonempty_coordinates_require_a_separate_integer_capacity():
    _,work,kernels=fixture(10)
    ledger=selection_gpu_work(work,kernels,compute_capability=[8,0])
    with pytest.raises(ValueError,match='int64 div/mod'):
        selection_gpu_service(ledger,resources(),traffic_mode='logical_hbm')
    first=selection_gpu_service(ledger,resources(),traffic_mode='logical_hbm',int64_divmod_per_second=1e9)
    second=selection_gpu_service(ledger,resources(),traffic_mode='logical_l2',int64_divmod_per_second=1e9)
    assert first['gpu_service_seconds']>=second['gpu_service_seconds']


@pytest.mark.parametrize('field,value',[('cuda_runtime','12.9'),('sm_count',80),('device_uuid','other'),('library_sha256','wrong')])
def test_census_cannot_transfer_to_another_installed_context(field,value):
    census,work,_=fixture(8);context=copy.deepcopy(census['context']);context[field]=value
    with pytest.raises(ValueError,match='context'):selection_gpu_census(work,census,context)


def test_census_rejects_changed_selection_geometry_source():
    census,work,_=fixture(0)
    census['source_sha256']['selection_geometry.py']='stale'
    with pytest.raises(ValueError,match='source changed'):
        selection_gpu_census(work,census,census['context'])


@pytest.mark.parametrize('change,match',[
    ('extra_kernel','Unassigned'),('missing_kernel','match source'),('duration','Duration-free'),
    ('count_block','count-reduction'),('select_grid','architecture policy'),('unknown_arch','architecture')])
def test_missing_or_incompatible_kernel_work_is_refused(change,match):
    _,work,kernels=fixture(8);arch=[8,0]
    if change=='extra_kernel':kernels.append(copy.deepcopy(kernels[-1]))
    elif change=='missing_kernel':kernels.pop()
    elif change=='duration':kernels[0]['duration']=1.
    elif change=='count_block':kernels[17]['geometry']['block']=[128,1,1]
    elif change=='select_grid':kernels[20]['geometry']['grid']=[1,1,1]
    elif change=='unknown_arch':arch=[7,0]
    with pytest.raises(ValueError,match=match):selection_gpu_work(work,kernels,compute_capability=arch)


def test_lookback_scenario_changes_source_traffic_without_changing_launches():
    _,work,kernels=fixture(8)
    a=selection_gpu_work(work,kernels,compute_capability=[8,0],lookback_windows=1)
    b=selection_gpu_work(work,kernels,compute_capability=[8,0],lookback_windows=3)
    assert a['kernel_count']==b['kernel_count']
    assert b['logical_bytes']-a['logical_bytes']==2*32*8*((455-1)+(2-1))


def whole_chunk(index):
    census=json.loads((Path(__file__).parent/'fixtures/device_selection_whole_chunk_a100.json').read_text())
    row=census['rows'][index]
    work=device_significant_tensor_work(row['N'],row['B'],row['K'],[b['retained'] for b in row['blocks']],max_cells=row['max_cells'])
    return census,row,work,selection_gpu_census(work,census,census['context'])


@pytest.mark.parametrize('index',range(7))
def test_whole_chunk_blocks_above_1m_cells_are_assigned_once(index):
    # One nonzero per chunk (DEVICE_SELECTION_MAX_CELLS): 2.1M to 33.5M cells.
    census,row,work,kernels=whole_chunk(index)
    assert len(work['blocks'])==1 and work['blocks'][0]['cells']==row['B']*row['K']>1<<20
    ledger=selection_gpu_work(work,kernels,compute_capability=census['context']['compute_capability'])
    used=[i for indices in ledger['operations'].values() for i in indices]
    used+=[i for group in ledger['nonzero_groups'] for indices in group.values() for i in indices]
    assert sorted(used)==list(range(len(kernels)))
    sweep=ledger['kernels'][ledger['nonzero_groups'][0]['flagged_select'][1]]
    assert sweep['cub_tiles']==-(-row['B']*row['K']//(384*6))


def test_count_grid_saturates_at_the_cub_occupancy_bound():
    census,row,work,kernels=whole_chunk(6)
    assert (row['B'],row['K'],row['mode'])==(4096,8192,'dense')
    ledger=selection_gpu_work(work,kernels,compute_capability=[8,0])
    count=ledger['kernels'][ledger['nonzero_groups'][0]['count'][0]]
    # 108 SMs x 8 resident blocks x subscription 5, below ceil(cells / 4096) = 8192.
    assert count['kernel']['geometry']['grid'][0]==4320
    assert count['read_bytes']==row['B']*row['K'] and count['write_bytes']==4*4320


def test_blocks_beyond_the_cuda_nonzero_limit_are_refused():
    census,row,work,kernels=whole_chunk(0)
    work=copy.deepcopy(work);work['blocks'][0]['cells']=1<<31
    with pytest.raises(ValueError,match='nonzero limit'):
        selection_gpu_work(work,kernels,compute_capability=[8,0])
