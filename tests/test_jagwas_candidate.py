"""Real encoded extents, captured projection geometry, synthetic service controls."""
import copy
from unittest.mock import patch
import numpy as np
import pytest
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.jagwas_candidate import jagwas_candidate_shape,jagwas_candidate_memory,jagwas_candidate_runtime
from torchgwas.linear import multigpu_variant_ranges
from torchgwas.pgen_work_census import census
from torchgwas.reduced_output_work import jagwas_indexed_part_work
from test_jagwas_scan_work import joint_fixture,scan_component
from test_jagwas_writer_service import bank,archive
from test_pgen_native_reader import write_pgen


@pytest.fixture
def input_path(tmp_path):
    data,profile=joint_fixture();n=data['samples'];m=4*profile['chunk_markers']
    path=tmp_path/'joint.pgen'
    write_pgen(path,(np.arange(m*n,dtype=np.uint32).reshape(m,n)%3).astype(np.uint8))
    return path


def candidate(path,count=2):
    data,profile=joint_fixture();b=profile['chunk_markers'];m=census(path,b)['markers'];tiles=[]
    caps=dict(cpu=3.,dram=1e9,input=1e8,output=1e8)
    for index,span in enumerate(multigpu_variant_ranges(m,b,count)):
        d,p=copy.deepcopy(data),copy.deepcopy(profile)
        d.update(markers=span[1]-span[0],phenotype_complete=True)
        d['encoded']=census(path,b,variant_range=span,include_chunks=True)
        p.update(result_ownership='owned',validate_range=True,event_wait_cpu_fraction=1.,
            decode_workers=3//count+(index<3%count),cpu_available_cores=caps['cpu'],
            write_bytes_per_second=caps['output'],fsync_seconds=1e-6,
            writeback_service=dict(pagecache_seconds_per_byte=1e-9,storage_seconds_per_byte=1e-8))
        for row in p['kernel_geometry']:row['validate_range']=True
        tiles.append(dict(device='cuda:'+str(index),data=d,profile=p,trait_range=[0,d['traits_analyzed']],variant_range=list(span)))
    return dict(tiles=tiles,trait_block=data['traits_analyzed'],devices=[t['device'] for t in tiles],
        partition_axis='variant',shared_capacities=caps,
        output=dict(block_bytes=None,queue_depth=2,store_beta=False,fsync=True))


def preparation(c):
    shape=jagwas_candidate_shape(c);shared=ExecutionGraph()
    shared.add('residual',.01,resources=dict(cpu=1.,host_serial=.5))
    device_graphs={};cleanup={}
    for device in shape['devices']:
        graph=ExecutionGraph();graph.add('factor',.02)
        graph.add('design',.001,after=['factor'],resources={device+':h2d':100.})
        device_graphs[device]=graph;cleanup[device]=[dict(seconds=.0001,resources=dict(cpu=1.,host_serial=.5))]
    columns=c['tiles'][0]['data'].get('covariate_columns',shape['covariates'])
    return dict(dimensions=[shape['samples'],shape['traits'],shape['covariates'],columns,shape['chunk_size'],shape['depth']],shared_graph=shared,
        device_graphs=device_graphs,cleanup=cleanup,finalize=[dict(seconds=.001)],unpriced_terms=['synthetic preparation control'])


def runtime(c,occupancy='dense',**kwargs):
    prepare=kwargs.pop('preparation',preparation(c))
    prices=dict(prices=bank(),archive=archive(),queue_cpu_seconds=dict(put=1e-6,get=1e-6))
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=scan_component):
        return jagwas_candidate_runtime(c,prices,preparation=prepare,occupancy=occupancy,host_serial_fraction=.5,**kwargs)


def test_shared_memory_and_narrow_queues_do_not_scale_with_trait_width(input_path):
    c=candidate(input_path);old=copy.deepcopy(c);memory=jagwas_candidate_memory(c)
    assert c==old
    assert memory['queued_result_bytes']==2*12*128
    assert memory['consumer_result_bytes']==2*12*128
    assert memory['shared_host_arrays']['phenotype_cast_and_result']==8*2049*512
    c['output']['queue_depth']=5
    more=jagwas_candidate_memory(c)
    assert more['host_bytes']-memory['host_bytes']==3*12*128
    assert more['device_bytes']==memory['device_bytes']
    c['trait_block']=7
    for tile in c['tiles']:
        tile['data']['traits_analyzed']=7;tile['trait_range']=[0,7]
    narrow=jagwas_candidate_memory(c)
    assert narrow['queued_result_bytes']==more['queued_result_bytes']
    assert narrow['writer_array_bytes']==more['writer_array_bytes']
    assert narrow['device_bytes']['cuda:0']<more['device_bytes']['cuda:0']
    assert all(row['persistent_factor_bytes']>=8*7*7 for row in narrow['devices'].values())
    one=jagwas_candidate_memory(candidate(input_path,1))
    assert one['queued_result_bytes']==0
    assert one['shared_host_bytes']==memory['shared_host_bytes']
    assert not one['prediction_complete']


def test_pinned_owned_depth_and_reserves_are_counted(input_path):
    c=candidate(input_path);old=jagwas_candidate_memory(c)
    for tile in c['tiles']:tile['profile']['depth']=3
    new=jagwas_candidate_memory(c,host_reserve_bytes=1000,device_reserve_bytes=2000)
    expected=1000+sum(row['pinned']['allocator_bytes']//2+29*128 for row in old['devices'].values())
    assert new['host_bytes']-old['host_bytes']==expected
    assert all(new['device_bytes'][device]>=old['device_bytes'][device]+2000 for device in c['devices'])
    assert all(row['pinned']['allocation_count']==18 for row in new['devices'].values())


@pytest.mark.parametrize('fault',['trait','borrowed','borrow_flag','incomplete','dtype','beta','coalescing','rank','readers','range'])
def test_invalid_executor_contract_refused_before_scan(fault,input_path):
    c=candidate(input_path);tile=c['tiles'][0]
    if fault=='trait':c['partition_axis']='trait'
    if fault=='borrowed':tile['profile']['result_ownership']='borrowed'
    if fault=='borrow_flag':tile['profile']['borrow_results']=True
    if fault=='incomplete':tile['data']['phenotype_complete']=False
    if fault=='dtype':tile['profile']['compute_dtype']='float64'
    if fault=='beta':c['output']['store_beta']=True
    if fault=='coalescing':c['output']['block_bytes']=1
    if fault=='rank':
        c['trait_block']=tile['data']['samples']
        for t in c['tiles']:t['data']['traits_analyzed']=c['trait_block'];t['trait_range']=[0,c['trait_block']]
    if fault=='readers':
        tile['profile']['decode_workers']=1;c['tiles'][1]['profile']['decode_workers']=2
    if fault=='range':tile['variant_range'][0]+=1
    with pytest.raises(ValueError):jagwas_candidate_memory(c)


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
def test_joint_graph_conserves_shared_resources_preparation_and_output(input_path,policy):
    c=candidate(input_path);old=copy.deepcopy(c)
    graph=runtime(c,return_graph=True,host_serial_policy=policy);result=graph.solve()
    assert c==old
    for index in range(2):
        assert result['start'][f'tile:{index}:prepare:factor']>=result['end']['shared_prepare:complete']
        assert result['start'][f'tile:{index}:submit_decode:0']>=result['end'][f'tile:{index}:prepare:complete']
    assert result['seconds']==result['end']['jagwas:finalize:done']
    for resource,capacity in graph.capacities.items():
        total=sum(graph.nodes[name][0]*demands.get(resource,0.) for name,demands in graph.demands.items())
        assert result['seconds']+1e-10>=total/capacity
    report=runtime(c)
    from torchgwas.executor_timing import SETUP_SCAN_WRITE_METRIC, SETUP_SCAN_WRITE_BOUNDARY
    assert report['observed_metric'] == SETUP_SCAN_WRITE_METRIC
    assert report['timing_boundary'] == SETUP_SCAN_WRITE_BOUNDARY
    assert report['resource_balance']['scheduled_seconds'] == report['estimated_seconds']
    assert 0 < report['resource_balance']['resource_lower_bound_seconds'] <= report['estimated_seconds']
    assert report['parts']==4 and report['retained_variants']==512
    assert report['indexed_part_bytes']==4*jagwas_indexed_part_work(128)['file_bytes']
    assert sum(row['d2h_bytes'] for row in report['shards'])==17*512
    assert report['shared_preprocessing_passes']==report['genotype_passes']==1
    assert not report['automatic_selection_ready']
    empty=runtime(c,occupancy='empty')
    assert empty['parts']==empty['indexed_part_bytes']==0
    assert empty['estimated_seconds']<report['estimated_seconds']


def test_single_device_and_burst_counts_use_the_same_consumer(input_path):
    report=runtime(candidate(input_path,1),occupancy=[[0,1,0,2]])
    assert report['parts']==2 and report['retained_variants']==3
    assert report['queue_depth']==0
    assert report['indexed_part_bytes']==sum(jagwas_indexed_part_work(i)['file_bytes'] for i in [1,2])
    with pytest.raises(ValueError,match='retained count'):runtime(candidate(input_path),occupancy=[[0],[0]])
    with pytest.raises(ValueError,match='Retained'):runtime(candidate(input_path),occupancy=[[129,0],[0,0]])


def test_physical_storage_and_shared_links_include_preparation_once(input_path):
    c=candidate(input_path);c['shared_storage_bytes_per_second']=1e6
    c['shared_links']=[dict(devices=c['devices'],h2d_bytes_per_second=1e7,d2h_bytes_per_second=1e7)]
    graph=runtime(c,return_graph=True);result=graph.solve()
    for key in ['storage','link:0:h2d','link:0:d2h']:
        amount=sum(graph.nodes[name][0]*d.get(key,0.) for name,d in graph.demands.items())
        assert amount>0 and result['seconds']+1e-10>=amount/graph.capacities[key]
    assert graph.demands['tile:0:prepare:design']['link:0:h2d']==100.


@pytest.mark.parametrize('fault',['missing','dimensions','device','empty','finalize','unbounded','limit'])
def test_missing_preparation_or_graph_expansion_refused(fault,input_path):
    c=candidate(input_path);prep=preparation(c);kwargs={}
    if fault=='missing':prep={}
    if fault=='dimensions':prep['dimensions'][1]+=1
    if fault=='device':prep['device_graphs'].pop('cuda:0')
    if fault=='empty':prep['shared_graph']=ExecutionGraph()
    if fault=='finalize':prep['finalize']=[]
    if fault=='unbounded':prep['shared_graph'].demands['residual']['unnamed_gpu']=1.
    if fault=='limit':kwargs['max_source_chunks']=3
    with patch('torchgwas.jagwas_candidate.jagwas_candidate_shape',wraps=jagwas_candidate_shape):
        with pytest.raises(ValueError):runtime(c,preparation=prep,**kwargs)


def plan(candidates,**kwargs):
    from torchgwas.jagwas_candidate import detailed_jagwas_plan
    defaults=dict(preparations=[preparation(c) for c in candidates],
        occupancy_scenarios={'empty':'empty','full':'dense'},host_scenarios={'half':dict(host_serial_fraction=.5)},
        cpu_workers=3,host_memory_bytes=1<<40,device_memory_bytes={'cuda:0':1<<40,'cuda:1':1<<40})
    defaults.update(kwargs)
    prices=dict(prices=bank(),archive=archive(),queue_cpu_seconds=dict(put=1e-6,get=1e-6))
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=scan_component):
        return detailed_jagwas_plan(candidates,prices,**defaults)


def test_bounded_optimizer_preserves_memory_under_empty_output_and_api_mapping(input_path):
    choices=[candidate(input_path,1),candidate(input_path,2)]
    result=plan(choices)
    scores=[max(row['estimate']['estimated_seconds'] for row in item['scenarios']) for item in result['candidates']]
    assert scores==sorted(scores)
    chosen=result['selected'];kwargs=chosen['api_kwargs']
    assert kwargs['variant_devices']==chosen['devices']
    assert kwargs['reduce']=='jagwas' and 'trait_block' not in kwargs
    assert not result['automatic_selection_ready'] and not result['selection_validated']
    before=[jagwas_candidate_memory(c)['host_bytes'] for c in choices]
    empty=plan(choices,occupancy_scenarios={'none':'empty'})
    for row in empty['candidates']:assert row['memory']['host_bytes']==before[row['candidate_index']]
    limit=min(before)-1
    with patch('torchgwas.jagwas_candidate.jagwas_candidate_runtime',side_effect=AssertionError('scored infeasible geometry')):
        with pytest.raises(ValueError,match='No feasible'):plan(choices,host_memory_bytes=limit)


def test_bounded_optimizer_refuses_budget_overflow_and_missing_capacity(input_path):
    choices=[candidate(input_path,1),candidate(input_path,2)]
    for options,message in [({'max_candidates':1},'bounded'),({'max_scenario_evaluations':3},'max_scenario'),
        ({'device_memory_bytes':{'cuda:0':1<<40}},'Missing'),({'cpu_workers':2},'No feasible')]:
        with pytest.raises(ValueError,match=message):plan(choices,**options)


def test_named_scenario_pairs_cannot_collide(input_path):
    result=plan([candidate(input_path,1)],occupancy_scenarios={'a:b':'empty','a':'dense'},
        host_scenarios={'c':dict(host_serial_fraction=0.),'b:c':dict(host_serial_fraction=1.)})
    assert len(result['selected']['scenarios'])==4


def test_owned_release_service_moves_to_the_single_consumer(input_path):
    from torchgwas.mechanistic_torch import torch_scan_work
    c=candidate(input_path)
    def with_release(data,profile):
        work=torch_scan_work(data,profile)
        for block in work['blocks']:block['discard_seconds']=.007
        return work
    for occupancy,stage in [('dense','write:3'),('empty','select:5')]:
        with patch('torchgwas.mechanistic_torch.torch_scan_work',side_effect=with_release):
            graph=runtime(c,occupancy=occupancy,return_graph=True)
        names=[name for name in graph.nodes if name.startswith('tile:') and name.endswith(stage)]
        assert len(names)==4
        assert sum(graph.nodes[name][0]*graph.demands[name]['cpu'] for name in names)==pytest.approx(.028)
        assert all(graph.demands[name]['host_serial']==1. for name in names)
        assert all(':jagwas:' in name for name in names)


def test_candidate_memory_includes_queried_host_workspace_per_factor(input_path, monkeypatch):
    # The census is xpotrf's: the rounding cutoff's Cholesky factor.
    monkeypatch.setenv('TORCHGWAS_JAGWAS_RCOND', '0')
    from test_cusolver_memory import evidence
    c=candidate(input_path);before=jagwas_candidate_memory(c)
    census,profile=evidence();census['rows'][0].update(traits=512,device_workspace_bytes=1152,host_workspace_bytes=32768)
    profile['jagwas_factor_workspace_census']=census
    after=jagwas_candidate_memory(c,device_memory_profiles={device:profile for device in c['devices']})
    assert after['host_bytes']-before['host_bytes']==2*32768
    assert all(row['host_arrays']['factor_host_workspace']==32768 for row in after['devices'].values())


def test_exact_tail_census_is_kept_and_missing_tail_geometry_is_refused(tmp_path):
    path=tmp_path/'tail.pgen'
    write_pgen(path,np.zeros((513,2049),np.uint8))
    c=candidate(path)
    assert [tile['variant_range'] for tile in c['tiles']]==[[0,384],[384,513]]
    memory=jagwas_candidate_memory(c)
    assert memory['shape']['markers']==513
    assert memory['queued_result_bytes']==2*12*128
    with pytest.raises(ValueError,match='geometry'):runtime(c)


@pytest.mark.parametrize('dimension',[3,4,5])
def test_preparation_cannot_reuse_another_covariate_chunk_or_ring_shape(input_path,dimension):
    c=candidate(input_path);prep=preparation(c);prep['dimensions'][dimension]+=1
    with pytest.raises(ValueError,match='Preparation dimensions'):runtime(c,preparation=prep)


def test_optimization_cannot_change_input_covariate_column_count(input_path):
    first=candidate(input_path,1);second=candidate(input_path,2)
    for tile in second['tiles']:tile['data']['covariate_columns']=3
    with pytest.raises(ValueError,match='same statistical workload'):plan([first,second])
