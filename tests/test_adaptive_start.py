"""One admitted starting layout must not hide an upfront performance search."""
import copy
from unittest.mock import patch
import pytest

from torchgwas.adaptive_start import prepare_adaptive_start
from torchgwas.pgen_work_census import census
from torchgwas.productive_run import ProductiveTuningRun
from test_adaptive_candidate import memory_profile,joint_fixture
from test_trait_candidate_space import spec as dense_spec


def options(context,*,reduction=None,axis='trait',width=2,sizes=(2,4)):
    return dict(chunk_sizes=list(sizes),initial_size=sizes[0],partition_axis=axis,
        trait_block=width,reduction=reduction,cpu_workers=4,host_memory_bytes=1<<40,
        device_memory_bytes={d:1<<40 for d in context['devices']},
        device_memory_profiles={d:memory_profile() for d in context['devices']},
        host_reserve_bytes=1<<20,device_reserve_bytes=1<<20)


@pytest.mark.parametrize('reduction',[None,'significant'])
def test_single_tiled_start_builds_no_runtime_candidate_and_keeps_all_work_shapes(tmp_path,reduction):
    spec=dense_spec(tmp_path);context=spec['contexts'][1];output=dict(spec['output'],block_bytes=None)
    before=copy.deepcopy((spec,context,output))
    retained=[];certified=[]
    with patch('torchgwas.trait_candidate_space.census',wraps=census) as collect,\
         patch('torchgwas.execution_graph.ExecutionGraph.solve',side_effect=AssertionError('upfront graph solve')),\
         patch('torchgwas.mechanistic_torch.torch_scan_work',side_effect=AssertionError('upfront scan timing')):
        result=prepare_adaptive_start(spec['workload'],context,output=output,
            **options(context,reduction=reduction),significance_threshold=.01 if reduction else None,
            _header_receiver=retained.append,_index_receiver=certified.append)
    assert len(retained)==1 and retained[0][0]==result['input_file_identity']
    assert len(certified)==1 and certified[0][0]==retained[0][0]
    assert certified[0][1] is retained[0][1]
    assert certified[0][2].dtype.name=='uint32'
    assert not certified[0][2].flags.writeable
    assert not retained[0][1].record_offsets.flags.writeable
    assert collect.call_count==0
    assert result['structural_work']['header_passes']==1
    assert result['structural_work']['census_passes']==0
    assert result['source_census'] is None
    assert result['structural_work']['layouts_built']==1
    assert result['structural_work']['runtime_candidates_evaluated']==0
    assert result['capacity']==4 and result['initial_size']==2
    assert result['api_kwargs']['chunk_size']==4
    assert [p['trait_range'] for p in result['partitions']]==[[0,2],[2,4],[4,5]]
    assert all(row['variant_range']==[0,10] for row in result['partitions'])
    assert not result['memory']['missing_geometry']
    for tile in result['candidate']['tiles']:
        assert {row['B'] for row in tile['profile']['kernel_geometry']}=={2,4}
    if reduction:
        assert result['api_kwargs']['reduce']=='significant'
        assert result['api_kwargs']['significance_threshold']==.01
        assert result['memory']['base_fixed_capacity_memory']['selection_bytes_by_device']
    assert (spec,context,output)==before
    # The execution driver can begin work at the smaller size without changing
    # the allocation cap; none of these productive reservations runs a planner.
    run=ProductiveTuningRun(result['partitions'],chunk_sizes=result['chunk_sizes'],initial=result['initial_size'])
    assert run.for_partition('0')(0,10,4)==2
    assert run.snapshot()['planning']['steps']==[]


def test_admission_reuses_source_header_without_reopening_index(tmp_path):
    from torchgwas.analytical_plan_cache import input_identity
    from torchgwas.pgen_reader import read_header
    spec=dense_spec(tmp_path);context=spec['contexts'][1]
    path=spec['workload']['genotype']
    prepared=(input_identity(path),read_header(path))
    with patch('torchgwas.trait_candidate_space.read_header',
               side_effect=AssertionError('index reread')):
        reused=prepare_adaptive_start(spec['workload'],context,output=spec['output'],
            **options(context),_prepared_header=prepared)
    fresh=prepare_adaptive_start(spec['workload'],context,output=spec['output'],
        **options(context))
    assert reused['memory']==fresh['memory']
    assert reused['partitions']==fresh['partitions']
    assert reused['source_layout']==fresh['source_layout']


@pytest.mark.parametrize('reduction',[None,'significant'])
def test_compact_start_matches_full_memory_admission(tmp_path,reduction):
    spec=dense_spec(tmp_path);context=spec['contexts'][1]
    output=dict(spec['output'],block_bytes=None) if reduction else spec['output']
    kw=options(context,reduction=reduction)
    if reduction:kw['significance_threshold']=.01
    full=prepare_adaptive_start(spec['workload'],context,output=output,**kw)
    compact=prepare_adaptive_start(spec['workload'],context,output=output,**kw,
        _compact_memory=True)
    assert compact['memory']==full['memory']
    assert compact['partitions']==full['partitions']
    assert compact['api_kwargs']==full['api_kwargs']
    assert compact['source_layout']['logical_chunks']==len(full['source_layout']['chunks'])


def test_compact_jagwas_start_matches_full_memory_admission(tmp_path):
    path,candidate,_=joint_fixture(tmp_path)
    context=dict(name='single',devices=candidate['devices'],
        profiles={t['device']:t['profile'] for t in candidate['tiles']},
        shared_capacities=candidate['shared_capacities'])
    workload=dict(genotype=str(path),samples=2049,markers=1025,traits=512,covariates=2,
        matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
    kw=options(context,reduction='jagwas',axis='variant',width=None,sizes=(128,256,512))
    full=prepare_adaptive_start(workload,context,output=candidate['output'],**kw)
    compact=prepare_adaptive_start(workload,context,output=candidate['output'],**kw,
        _compact_memory=True)
    assert compact['memory']==full['memory']
    assert compact['partitions']==full['partitions']
    assert compact['api_kwargs']==full['api_kwargs']


def test_jagwas_start_keeps_joint_panel_and_admits_shifted_ld_reads(tmp_path):
    path,candidate,source=joint_fixture(tmp_path)
    context=dict(name='single',devices=candidate['devices'],
        profiles={t['device']:t['profile'] for t in candidate['tiles']},
        shared_capacities=candidate['shared_capacities'])
    workload=dict(genotype=str(path),samples=2049,markers=1025,traits=512,covariates=2,
        matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
    with patch('torchgwas.trait_candidate_space.census',wraps=census) as collect:
        result=prepare_adaptive_start(workload,context,output=candidate['output'],
            **options(context,reduction='jagwas',axis='variant',width=None,sizes=(128,256,512)))
    assert collect.call_count==0
    from torchgwas.decoder_work import native_read_layout
    for exact,metadata in zip(source['chunks'],result['source_layout']['chunks']):
        a,b=native_read_layout(exact),native_read_layout(metadata)
        assert a['read_bytes']==b['read_bytes']
        assert a['extra_workspace_bytes']==b['extra_workspace_bytes']
    assert result['partitions']==[dict(id='0',device='cuda:0',trait_range=[0,512],variant_range=[0,1025])]
    assert 'trait_block' not in result['api_kwargs']
    assert result['api_kwargs']['variant_devices']==['cuda:0']
    assert result['memory']['decoder_extra_bytes_by_device']['cuda:0']>0
    assert result['memory']['required_geometry']
    assert result['memory']['host_bytes']>result['memory']['base_fixed_capacity_memory']['host_bytes']


@pytest.mark.parametrize('fault',['initial','axis','jagwas_tiles','significant_variants','threshold','zero_reserve',
    'workers','budget_devices','grid','census_budget','tile_budget'])
def test_invalid_start_fails_before_source_census(tmp_path,fault):
    spec=dense_spec(tmp_path);context=spec['contexts'][1];args=options(context)
    if fault=='initial':args['initial_size']=3
    if fault=='axis':args['partition_axis']='sample'
    if fault=='jagwas_tiles':args['reduction']='jagwas'
    if fault=='significant_variants':args.update(reduction='significant',partition_axis='variant',trait_block=None)
    if fault=='threshold':args['significance_threshold']=.05
    if fault=='zero_reserve':args['host_reserve_bytes']=0
    if fault=='workers':args['cpu_workers']=1
    if fault=='budget_devices':args['device_memory_bytes'].pop('cuda:1')
    if fault=='grid':args['chunk_sizes']=[2,3]
    if fault=='census_budget':args['max_census_chunks']=3
    if fault=='tile_budget':args['max_tiles']=1
    with patch('torchgwas.trait_candidate_space.census',side_effect=AssertionError('premature source census')):
        with pytest.raises(ValueError):prepare_adaptive_start(spec['workload'],context,output=spec['output'],**args)


@pytest.mark.parametrize('fault',['host_capacity','device_capacity','input_change'])
def test_unadmitted_layout_is_never_returned(tmp_path,fault):
    spec=dense_spec(tmp_path);context=spec['contexts'][1];args=options(context)
    if fault=='host_capacity':args['host_memory_bytes']=1
    if fault=='device_capacity':args['device_memory_bytes']['cuda:0']=1
    from torchgwas.adaptive_start import input_identity
    identity=input_identity(spec['workload']['genotype'])
    def changed(path):
        changed.calls+=1
        return dict(identity,mtime_ns=identity['mtime_ns']+changed.calls-1) if fault=='input_change' else identity
    changed.calls=0
    with patch('torchgwas.adaptive_start.input_identity',side_effect=changed):
        with pytest.raises(ValueError):prepare_adaptive_start(spec['workload'],context,output=spec['output'],**args)


@pytest.mark.parametrize('fault',['missing','ambiguous'])
def test_missing_timing_geometry_does_not_prevent_memory_admitted_work(tmp_path,fault):
    spec=dense_spec(tmp_path);context=spec['contexts'][1]
    for profile in context['profiles'].values():
        profile['kernel_geometry']=[] if fault=='missing' else profile['kernel_geometry']*2
    result=prepare_adaptive_start(spec['workload'],context,output=spec['output'],**options(context))
    assert result['timing_geometry_missing'] and not result['runtime_work_ready']
    assert result['memory']['device_bytes'] and result['source_census'] is None
    from torchgwas.mechanistic_torch import torch_scan_work
    first=result['candidate']['tiles'][0]
    with pytest.raises(ValueError,match='Memory-only'):
        torch_scan_work(first['data'],first['profile'])


def test_finer_admission_grid_is_bounded_without_adding_a_performance_candidate(tmp_path):
    from torchgwas.trait_candidate_space import prepare_trait_candidates
    spec=dense_spec(tmp_path);context=spec['contexts'][0]
    result=prepare_trait_candidates(spec['workload'],[context],chunks=[4],trait_blocks=[2],
        output=spec['output'],census_chunk_size=2)
    assert len(result['candidates'])==1 and result['census_passes']==1
    assert result['source_census']['chunk_markers']==2 and len(result['source_census']['chunks'])==5
    assert result['census_chunks']==8
    with patch('torchgwas.trait_candidate_space.census',side_effect=AssertionError('budget must reject first')):
        with pytest.raises(ValueError):prepare_trait_candidates(spec['workload'],[context],chunks=[4],trait_blocks=[2],
            output=spec['output'],census_chunk_size=2,max_census_chunks=7)
