"""Source payload, dense memory admission and shared-graph candidate controls."""
import copy
from unittest.mock import patch
import numpy as np
import pytest
from test_trait_tiling_model import candidate, input_path
from test_mechanistic_shapes import component
from torchgwas.significant_host_work import (indexed_part_work, host_significant_selection_work,
    host_selection_service, significant_host_memory)
from torchgwas.significant_host_model import significant_host_runtime
from torchgwas.significant_host_plan import detailed_significant_host_plan
from torchgwas.numpy_nonzero_work import nonzero_protocol


def bank():
    return dict(host_selector='bounded_flat_v2',predicate_max_cells=1<<20,nonzero_protocol=nonzero_protocol(),prices={name:dict(call_cpu_seconds=1e-6,unit_cpu_seconds=1e-9,dram_bytes_per_unit=8)
        for name in ['critical_one','critical_lookup','critical_round','mask_allocate','predicate_block','flatnonzero_empty','flatnonzero_sparse','flatnonzero_dense','coordinate_divmod','matrix_gather_flat','df_gather_row','inplace_index_add','index_cast','index_add']},
        archive={str(b):dict(call_cpu_seconds=1e-6,byte_cpu_seconds=1e-9) for b in [False,True]},
        queue_cpu_seconds=dict(put=1e-6,get=1e-6))


def runtime(c, occupancy='dense', **kwargs):
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return significant_host_runtime(c,bank(),occupancy=occupancy,host_serial_fraction=.5,**kwargs)


@pytest.mark.parametrize('rows',[1,37,4096])
@pytest.mark.parametrize('store_beta',[True,False])
def test_indexed_part_header_and_payload_bytes_match_numpy(tmp_path,rows,store_beta):
    work = indexed_part_work(rows,store_beta=store_beta)
    values = {array['field']:np.arange(rows,dtype=array['dtype']) for array in work['arrays']}
    path = tmp_path/'part.npz'
    with path.open('wb') as handle:
        np.savez(handle,**values)
    assert work['file_bytes'] == path.stat().st_size
    assert work['array_payload_bytes'] == sum(value.nbytes for value in values.values())
    assert work['fsync_calls'] == 1
    assert indexed_part_work(0)['file_bytes'] == 0


def test_large_archive_contract_is_refused_without_allocating_payload():
    with pytest.raises(ValueError,match='ZIP'):
        indexed_part_work(1 << 29)


def test_dense_memory_admission_counts_shared_queue_without_multiplying_per_device(input_path):
    c = candidate(input_path,width=2,count=2,block_bytes=None)
    first = significant_host_memory(c,host_reserve_bytes=1000,device_reserve_bytes=2000)
    c['output']['queue_depth'] = 3
    second = significant_host_memory(c,host_reserve_bytes=1000,device_reserve_bytes=2000)
    assert second['host_bytes']-first['host_bytes'] == 2*28*4*2
    assert first['device_bytes'] == second['device_bytes']
    assert first['shared_critical_table_bytes'] == 8*(32+1)
    single = significant_host_memory(candidate(input_path,width=2,count=1,block_bytes=None))
    assert single['queued_result_bytes'] == 0 and single['layout']['queue_depth'] == 0
    assert first['occupancy'].startswith('all variant-trait pairs')


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
def test_reduced_bridge_keeps_setup_reader_lifetime_and_durable_writer(input_path,policy):
    c = candidate(input_path,block_bytes=None)
    original = copy.deepcopy(c)
    graph = runtime(c,return_graph=True,host_serial_policy=policy)
    result = graph.solve()
    assert c == original
    assert result['start']['tile:0:submit_decode:0'] >= result['end']['tile:0:prepare:complete']
    assert result['start']['tile:2:prepare:covariate_basis'] >= result['end']['tile:0:significant:producer_complete']
    assert result['seconds'] == result['end']['significant:finalize:done']
    for resource,capacity in graph.capacities.items():
        work = sum(graph.nodes[name][0]*demands.get(resource,0.) for name,demands in graph.demands.items())
        assert result['seconds']+1e-10 >= work/capacity
    report = runtime(c)
    expected = sum(indexed_part_work(b*k)['file_bytes'] for k in [2,2,1] for b in [4,4,2])
    assert report['indexed_part_bytes'] == expected and report['parts'] == 9
    assert any('page residency' in term for term in report['unpriced_terms'])
    empty = runtime(c,occupancy='empty')
    assert empty['indexed_part_bytes'] == 0 and empty['parts'] == 0
    assert empty['estimated_tile_seconds'] < report['estimated_tile_seconds']


def test_shared_read_write_and_pcie_resources_are_charged_once(input_path):
    c = candidate(input_path,block_bytes=None)
    c['shared_storage_bytes_per_second'] = 1.
    c['shared_links'] = [dict(devices=c['devices'],h2d_bytes_per_second=1.,d2h_bytes_per_second=1.)]
    graph = runtime(c,return_graph=True)
    result = graph.solve()
    for resource in ['storage','link:0:h2d','link:0:d2h']:
        work = sum(graph.nodes[name][0]*demands.get(resource,0.) for name,demands in graph.demands.items())
        assert work>0 and result['seconds']+1e-8 >= work


def test_explicit_counts_preserve_bursts_and_threshold_one_uses_its_own_primitive(input_path):
    c = candidate(input_path,block_bytes=None)
    counts = [[1,0,0],[0,1,0],[0,0,1]]
    result = runtime(c,occupancy=counts)
    assert result['parts'] == 3 and result['indexed_part_bytes'] == 3*indexed_part_work(1)['file_bytes']
    with pytest.raises(ValueError,match='One exact retained'):
        runtime(c,occupancy=[[0],[0],[0]])
    prices = bank()['prices']
    prices.pop('critical_lookup')
    work = host_significant_selection_work(4,3,0,threshold_one=True)
    assert len(host_selection_service(work,prices,cpu_fraction=1.,dram_bytes_per_second=1e9,host_serial_fraction=0.)) == 11


def options():
    return dict(significance_threshold=.05,occupancy_scenarios=dict(null='empty',dense='dense'),
        host_scenarios=dict(fluid=dict(host_serial_fraction=0.)),cpu_workers=2,
        host_memory_bytes=1<<30,device_memory_bytes={'cuda:0':1<<30,'cuda:1':1<<30})


def test_bounded_host_selection_keeps_statistic_and_checks_memory_before_graph(input_path):
    candidates = [candidate(input_path,width=2,count=2,block_bytes=None),candidate(input_path,width=5,count=1,block_bytes=None)]
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        result = detailed_significant_host_plan(candidates,bank(),**options())
    assert result['selected']['worst_supplied_scenario_seconds'] == min(r['worst_supplied_scenario_seconds'] for r in result['candidates'])
    selected = result['selected']
    assert selected['api_kwargs']['reduce'] == 'significant' and selected['api_kwargs']['significance_threshold'] == .05
    assert selected['required_environment']['TORCHGWAS_SIGNIFICANCE_BACKEND'] == 'host'
    assert not result['selection_validated']
    args = options(); args['host_memory_bytes'] = 1
    with patch('torchgwas.significant_host_plan.significant_host_runtime',side_effect=AssertionError('must reject first')):
        with pytest.raises(ValueError,match='host_memory'):
            detailed_significant_host_plan(candidates,bank(),**args)

def test_oversized_archive_proposal_reaches_memory_rejection(input_path):
    # No phenotype or selected arrays are allocated. Dense output would exceed
    # the priced ZIP offset contract, but this candidate already exceeds RAM.
    huge = candidate(input_path,width=1<<25,count=1,traits=1<<25,block_bytes=None)
    memory = significant_host_memory(huge)
    assert memory['archive_buffer_bytes'] == 32 << 20
    assert memory['host_bytes'] > 1 << 30
    with patch('torchgwas.significant_host_plan.significant_host_runtime',side_effect=AssertionError('Reject before pricing')):
        with pytest.raises(ValueError,match='host_memory'):
            detailed_significant_host_plan([huge],bank(),**options())


def test_old_host_prices_cannot_silently_price_new_selector(input_path):
    c = candidate(input_path, block_bytes=None)
    prices = bank(); prices.pop('host_selector')
    with pytest.raises(ValueError, match='host selector'):
        significant_host_runtime(c, prices, occupancy='empty', host_serial_fraction=0.)
    prices = bank(); prices['predicate_max_cells'] //= 2
    with pytest.raises(ValueError, match='predicate limit'):
        significant_host_runtime(c, prices, occupancy='empty', host_serial_fraction=0.)


@pytest.mark.parametrize('shape', [(512,2048),(1024,8193),(3,(1<<20)+5)])
def test_host_ledger_tracks_predicate_block_tail_and_owned_coordinates(shape):
    from torchgwas.host_significance import predicate_block_shape
    work = host_significant_selection_work(*shape, retained=7)
    height, width, calls = predicate_block_shape(*shape)
    assert work['predicate_block_shape'] == [height,width]
    assert work['predicate_calls'] == calls
    assert work['nonzero_output_indices'] == work['coordinate_divmod_elements'] == 7
    assert work['index_cast_elements'] == work['inplace_index_add_elements'] == 7
    prices = bank()['prices']
    prices['predicate_block']['call_cpu_seconds'] = .25
    prices['predicate_block']['unit_cpu_seconds'] = 0.
    steps = host_selection_service(work, prices, cpu_fraction=1., dram_bytes_per_second=1e30, host_serial_fraction=0.)
    assert steps[3]['seconds'] == calls*.25
    assert sum(s['seconds']*s['resources']['dram'] for s in steps) > 25*shape[0]*shape[1]


@pytest.mark.parametrize('return_beta',[False,True])
@pytest.mark.parametrize('step',[1,7,None])
def test_selected_allocation_extents_match_production_owned_storage(return_beta,step):
    from torchgwas.host_significance import select_host_pairs
    values=np.zeros((13,19),np.float32)
    if step is not None:values.ravel()[::step]=3.
    beta=np.full_like(values,2.) if return_beta else None
    result=select_host_pairs(beta,values,np.full((13,1),97.,np.float32),np.full((13,1),2.,np.float64))
    work=host_significant_selection_work(13,19,len(result[0]),return_beta=return_beta)
    actual={name:array for name,array in zip(['variant_index','trait_index','beta','t_stat','df'],result) if array is not None}
    ledger={a['name']:a for a in work['allocations'] if a['lifetime']=='selected_result'}
    assert set(ledger)==set(actual)
    for name,array in actual.items():
        assert ledger[name]['bytes']==array.nbytes
        assert not np.shares_memory(array,values)
    assert sum(a.nbytes for a in actual.values())==work['selected_array_bytes']
    assert work['allocation_data_bytes']==values.size+sum(a.nbytes for a in actual.values())+8*len(result[0])


def test_trait_rebase_prices_owned_copy_and_inplace_add_without_allocating_add():
    prices=bank()['prices'];prices.pop('index_add')
    prices['inplace_index_add']=dict(call_cpu_seconds=.002,unit_cpu_seconds=.003,dram_bytes_per_unit=16)
    work=host_significant_selection_work(13,19,37)
    steps=host_selection_service(work,prices,cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)
    assert steps[8]['seconds']==steps[10]['seconds']==.002+37*.003
    assert steps[10]['seconds']*steps[10]['resources']['dram']==pytest.approx(16*37)
    assert [(a['name'],a['bytes']) for a in work['allocations'] if a['lifetime']=='rebased_result']==[('rebased_trait_index',8*37)]
    assert not any(a['lifetime']=='outer_rebase_temporary' for a in work['allocations'])
