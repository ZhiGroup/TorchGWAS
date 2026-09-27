"""Bounded JAGWAS proposals with real encoded extents and captured kernels.

Prices are synthetic resource controls; these tests do not certify ranking.
"""
import copy
import json
from unittest.mock import patch

import pytest

from test_jagwas_actual_candidate import input_path, actual_candidate, writer_prices
from test_jagwas_preparation import profile as prepared_profile
from torchgwas.jagwas_candidate import detailed_jagwas_plan
from torchgwas.jagwas_candidate_space import prepare_jagwas_candidates, bounded_jagwas_plan
from torchgwas.jagwas_preparation import build_jagwas_preparation
from torchgwas.pgen_work_census import census
from torchgwas.trait_candidate_space import prepare_trait_candidates


def specification(path):
    contexts=[];services={}
    hosts={'parallel':dict(host_serial_fraction=0.), 'serial':dict(host_serial_fraction=1.)}
    for count in (1,2):
        example=actual_candidate(path,128,count)
        name='gpus'+str(count)
        contexts.append(dict(name=name,devices=example['devices'],
            profiles={tile['device']:prepared_profile(tile['profile']) for tile in example['tiles']},
            shared_capacities=example['shared_capacities']))
        services[name]={key:dict(library_arithmetic={d:'scalar' for d in example['devices']},
            shared_cpu_steps=[dict(seconds=.001,resources=dict(cpu=1.,host_serial=host['host_serial_fraction']))],
            finalize=[dict(seconds=.001,resources=dict(cpu=1.,host_serial=host['host_serial_fraction']))])
            for key,host in hosts.items()}
    return dict(workload=dict(genotype=str(path),samples=2049,markers=1025,traits=512,covariates=2,
            matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True),
        contexts=contexts,bounds=dict(chunks=[128,256,512]),output=example['output'],prices=writer_prices(),
        preparation_services=services,joint=dict(occupancy_scenarios={'none':'empty','all':'dense'},
            host_scenarios=hosts,cpu_workers=3,host_memory_bytes=1<<40,
            device_memory_bytes={'cuda:0':1<<40,'cuda:1':1<<40},max_scenario_evaluations=24))


def proposed(spec):
    return prepare_jagwas_candidates(spec['workload'],spec['contexts'],output=spec['output'],**spec['bounds'])


def test_full_panel_space_reuses_one_census_and_preserves_exact_shards_and_tails(input_path):
    spec=specification(input_path);before=copy.deepcopy(spec)
    with patch('torchgwas.trait_candidate_space.census',wraps=census) as collect:
        space=proposed(spec)
    assert collect.call_count==1
    assert space['raw_candidates']==len(space['candidates'])==6
    assert space['reduction']=='jagwas' and not space['missing_geometry']
    assert spec==before
    for candidate in space['candidates']:
        assert candidate['partition_axis']=='variant' and candidate['trait_block']==512
        assert sum(t['data']['markers'] for t in candidate['tiles'])==1025
        for tile in candidate['tiles']:
            assert tile['trait_range']==[0,512] and tile['data']['phenotype_complete'] is True
            assert tile['profile']['result_ownership']=='owned'
            assert tile['profile']['reduction']=='jagwas'
            b=tile['profile']['chunk_markers'];m=tile['data']['markers']
            widths={min(b,m)} | ({m%b} if m%b else set())
            assert {r['B'] for r in tile['profile']['joint_kernel_geometry']}==widths
            assert {r['B'] for r in tile['profile']['kernel_geometry']}==widths
    assert any(row.get('reduction')=='jagwas' and row['shape'][1]==1 for row in space['required_geometry'])


@pytest.mark.parametrize('fault',['trait_partition','partial_width','rank','borrowed','wrong_reduction','beta','coalescing','budget'])
def test_invalid_joint_scope_and_budgets_fail_before_reading_pgen(input_path,fault):
    spec=specification(input_path)
    options=dict(chunks=[128],trait_blocks=[512],partition_axes=['variant'],reduction='jagwas')
    if fault=='trait_partition':options['partition_axes']=['trait']
    if fault=='partial_width':options['trait_blocks']=[256]
    if fault=='rank':spec['workload']['samples']=512
    if fault=='borrowed':spec['contexts'][0]['profiles']['cuda:0']['result_ownership']='borrowed'
    if fault=='wrong_reduction':spec['contexts'][0]['profiles']['cuda:0']['reduction']='device_significant'
    if fault=='beta':spec['output']['store_beta']=True
    if fault=='coalescing':spec['output']['block_bytes']=1024
    if fault=='budget':options['max_candidates']=1
    with patch('torchgwas.trait_candidate_space.read_header') as header, patch('torchgwas.trait_candidate_space.census') as collect:
        with pytest.raises(ValueError):
            prepare_trait_candidates(spec['workload'],spec['contexts'],output=spec['output'],**options)
        header.assert_not_called();collect.assert_not_called()


def test_generated_plan_matches_explicit_source_preparation_for_every_scenario(input_path):
    spec=specification(input_path);before=copy.deepcopy(spec);space=proposed(spec)
    assignments={row['candidate_index']:row['context'] for row in space['assignments']}
    calls=[]
    def build(index,candidate,name,host):
        calls.append((index,name))
        return build_jagwas_preparation(candidate,host_serial_fraction=host['host_serial_fraction'],
            **spec['preparation_services'][assignments[index]][name])['preparation']
    expected=detailed_jagwas_plan(space['candidates'],spec['prices'],preparation_factory=build,**spec['joint'])
    actual=bounded_jagwas_plan(**spec)
    assert len(calls)==12 and len(set(calls))==12  # Reused across both occupancies.
    assert actual['selected']==expected['selected']
    assert actual['candidates']==expected['candidates'] and spec==before
    assert actual['preparation_policy']=='per_host_scenario'
    assert actual['candidates_feasible']==6 and not actual['automatic_selection_ready']
    assert actual['selected']['api_kwargs']['reduce']=='jagwas'
    assert 'trait_block' not in actual['selected']['api_kwargs']
    for candidate in actual['candidates']:
        assert len(candidate['scenarios'])==4
        assert {r['estimate']['retained_variants'] for r in candidate['scenarios']}=={0,1025}
        assert all(r['estimate']['shared_preprocessing_passes']==1 for r in candidate['scenarios'])


def test_preparation_factory_runs_after_memory_and_evaluation_admission(input_path):
    spec=specification(input_path);space=proposed(spec)
    options=copy.deepcopy(spec['joint'])
    with patch('torchgwas.jagwas_candidate_space.build_jagwas_preparation',wraps=build_jagwas_preparation) as build:
        too_small=copy.deepcopy(spec);too_small['joint']['max_scenario_evaluations']=23
        with pytest.raises(ValueError,match='max_scenario_evaluations'):bounded_jagwas_plan(**too_small)
        build.assert_not_called()
        too_many_chunks=copy.deepcopy(spec);too_many_chunks['joint']['max_source_chunks']=1
        with pytest.raises(ValueError,match='source-chunk expansion'):bounded_jagwas_plan(**too_many_chunks)
        build.assert_not_called()
        impossible=copy.deepcopy(spec);impossible['joint']['host_memory_bytes']=1
        with pytest.raises(ValueError,match='No feasible'):bounded_jagwas_plan(**impossible)
        build.assert_not_called()
    # GPU 1 cannot hold a full panel. Its missing geometry must not require a
    # capture, and no factor/setup service is built for that rejected context.
    spec['joint']['device_memory_bytes']['cuda:1']=1
    spec['contexts'][1]['profiles']['cuda:1']['joint_kernel_geometry']=[]
    with patch('torchgwas.jagwas_candidate_space.build_jagwas_preparation',wraps=build_jagwas_preparation) as build:
        result=bounded_jagwas_plan(**spec)
        assert build.call_count==6
        assert all(call.args[0]['devices']==['cuda:0'] for call in build.call_args_list)
    assert result['candidates_feasible']==3 and len(result['rejected'])==3
    assert all(row['reason']=='device_memory' for row in result['rejected'])


def test_feasible_missing_projection_geometry_is_not_silently_skipped(input_path):
    spec=specification(input_path)
    p=spec['contexts'][1]['profiles']['cuda:1']
    p['joint_kernel_geometry']=[row for row in p['joint_kernel_geometry'] if row['B']!=1]
    space=proposed(spec)
    assert space['missing_geometry']
    assert all(row.get('reduction')=='jagwas' and row['shape']==[2049,1,512] for row in space['missing_geometry'])
    with pytest.raises(ValueError,match='joint geometry'):bounded_jagwas_plan(**spec)


@pytest.mark.parametrize('fault',['trait_bound','missing_context','missing_host','missing_final','missing_gpu'])
def test_serializable_services_cannot_omit_preparation_before_census(input_path,fault):
    spec=specification(input_path)
    if fault=='trait_bound':spec['bounds']['trait_blocks']=[256]
    if fault=='missing_context':spec['preparation_services'].pop('gpus2')
    if fault=='missing_host':spec['preparation_services']['gpus1'].pop('parallel')
    if fault=='missing_final':spec['preparation_services']['gpus1']['parallel']['finalize']=[]
    if fault=='missing_gpu':spec['preparation_services']['gpus2']['parallel']['library_arithmetic'].pop('cuda:1')
    with patch('torchgwas.trait_candidate_space.read_header') as header:
        with pytest.raises(ValueError):bounded_jagwas_plan(**spec)
        header.assert_not_called()


def test_cli_plans_from_json_without_manual_graphs(input_path,tmp_path,capsys):
    from torchgwas.autotune import main
    spec=specification(input_path);spec['model']='detailed_jagwas_space'
    path=tmp_path/'space.json';path.write_text(json.dumps(spec))
    with patch('sys.argv',['autotune',str(path)]):main()
    result=json.loads(capsys.readouterr().out)
    assert result['candidates_evaluated']==6
    assert result['search_space']['census_passes']==1
    assert result['preparation_policy']=='per_host_scenario'
    assert result['selected']['api_kwargs']['variant_devices'] in [['cuda:0'],['cuda:0','cuda:1']]
