import copy
import json
from unittest.mock import patch

import numpy as np
import pytest

from test_pgen_native_reader import write_pgen
from test_trait_tiling_model import candidate, component, options
from torchgwas.pgen_work_census import census
from torchgwas.trait_candidate_space import prepare_trait_candidates, bounded_trait_plan


def spec(tmp_path):
    path=tmp_path/'input.pgen'
    write_pgen(path,np.arange(320,dtype=np.uint8).reshape(10,32)%4)
    contexts=[]
    for count in (1,2):
        example=candidate(path,count=count)
        profiles={tile['device']:copy.deepcopy(tile['profile']) for tile in example['tiles']}
        for profile in profiles.values():
            profile['kernel_geometry']=[dict(N=32,B=b,K=k,C=8,validate_range=True,kernels=[])
                                        for b in (2,4) for k in (1,2,5)]
        contexts.append(dict(name=f'gpus{count}',devices=example['devices'],profiles=profiles,
                             shared_capacities=example['shared_capacities']))
    return dict(workload=dict(genotype=str(path),samples=32,markers=10,traits=5,covariates=8,
                             matching_sample_order=True,complete_phenotypes=True),contexts=contexts,
                bounds=dict(chunks=[2,4],trait_blocks=[2,5]),joint=options(),output=example['output'])


def prepare(value,**kwargs):
    return prepare_trait_candidates(value['workload'],value['contexts'],output=value['output'],**(value['bounds']|kwargs))


def test_bounded_cartesian_space_preserves_context_and_all_tails(tmp_path):
    value=spec(tmp_path);before=copy.deepcopy(value)
    with patch('torchgwas.trait_candidate_space.census',wraps=census) as collect:
        space=prepare(value)
    assert collect.call_count==1
    assert space['raw_candidates']==8 and len(space['candidates'])==6
    assert len(space['excluded'])==2 and all(r['reason']=='idle_devices' for r in space['excluded'])
    assert space['census_chunks']==8 and space['census_views']==2 and not space['missing_geometry']
    assert value==before
    small=next(c for c in space['candidates'] if c['trait_block']==2 and len(c['devices'])==2)
    assert [tile['trait_range'] for tile in small['tiles']]==[[0,2],[2,4],[4,5]]
    assert [tile['device'] for tile in small['tiles']]==['cuda:0','cuda:1','cuda:0']
    assert [tile['profile']['decode_workers'] for tile in small['tiles']]==[1,1,1]
    single=next(c for c in space['candidates'] if len(c['devices'])==1)
    assert single['tiles'][0]['profile']['decode_workers']==2
    for c in space['candidates']:
        encoded=c['tiles'][0]['data']['encoded']
        assert encoded['variant_range']==[0,10]
        assert all(tile['data']['encoded'] is encoded for tile in c['tiles'])
    assert {tuple(row['shape']) for row in space['required_geometry']} >= {(32,2,1,8),(32,4,2,8),(32,2,5,8)}


def test_optional_shared_transfer_prices_do_not_change_old_layouts(tmp_path):
    value=spec(tmp_path)
    original=prepare(value)
    for context in value['contexts']:
        context['shared_transfer_capacities']=dict(h2d=1e10,d2h=8e9)
    with_transfer=prepare(value)
    assert with_transfer['candidates']==original['candidates']
    value['contexts'][0]['shared_transfer_capacities']['h2d']=0.
    with pytest.raises(ValueError,match='shared H2D/D2H'):
        prepare(value)


def test_explicit_source_layout_is_preserved_only_for_whole_panel(tmp_path):
    value=spec(tmp_path)
    for flag in [True,False]:
        value['workload']['phenotype_c_contiguous']=flag
        space=prepare(value)
        for c in space['candidates']:
            assert all(t['data']['phenotype_c_contiguous']==(flag and c['trait_block']==5) for t in c['tiles'])
    value['workload']['phenotype_c_contiguous']=1
    with pytest.raises(ValueError,match='boolean'):prepare(value)


def test_profile_projection_keeps_only_required_geometry_and_independent_copies(tmp_path):
    value=spec(tmp_path);before=copy.deepcopy(value)
    projected=prepare(value)
    for c in projected['candidates']:
        for tile in c['tiles']:
            p=tile['profile'];k=tile['data']['traits_analyzed'];b=p['chunk_markers']
            assert {tuple(row[key] for key in ['N','B','K','C']) for row in p['kernel_geometry']}=={(32,size,k,8) for size in {min(b,10),10%b}-{0}}
    repeated=next(c for c in projected['candidates'] if c['trait_block']==2 and len(c['devices'])==1)
    repeated['tiles'][0]['profile']['kernel_geometry'][0]['kernels'].append({'test':'mutation'})
    assert not repeated['tiles'][1]['profile']['kernel_geometry'][0]['kernels']
    assert value==before


def test_projected_and_complete_geometry_banks_give_identical_plans(tmp_path):
    value=spec(tmp_path)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        projected=bounded_trait_plan(**value)
        with patch('torchgwas.trait_candidate_space._tile_profile',side_effect=lambda template,*args:copy.deepcopy(template)):
            complete=bounded_trait_plan(**value)
    assert projected==complete


@pytest.mark.parametrize('limit',[
    dict(max_candidates=7),dict(max_candidate_tiles=15),dict(max_census_chunks=7)])
def test_search_budget_failure_precedes_any_file_read(tmp_path,limit):
    value=spec(tmp_path)
    with patch('torchgwas.trait_candidate_space.read_header') as header, patch('torchgwas.trait_candidate_space.census') as collect:
        with pytest.raises(ValueError,match='before census'):prepare(value,**limit)
        header.assert_not_called();collect.assert_not_called()


@pytest.mark.parametrize('bad',['dimension','complete','covariates','context_device','workers','shared','duplicate','wide'])
def test_malformed_workload_and_context_are_refused(tmp_path,bad):
    value=spec(tmp_path)
    if bad=='dimension':value['workload']['markers']=11
    elif bad=='complete':value['workload']['complete_phenotypes']=False
    elif bad=='covariates':value['workload']['covariates']=-1
    elif bad=='context_device':value['contexts'][0]['devices'].append('cuda:1')
    elif bad=='workers':value['contexts'][1]['profiles']['cuda:1']['decode_workers']=2
    elif bad=='shared':value['contexts'][0]['shared_capacities']['input']*=1.01
    elif bad=='duplicate':value['bounds']['chunks']=[4,4]
    else:value['bounds']['trait_blocks']=[6]
    with pytest.raises(ValueError):prepare(value)


def test_missing_or_ambiguous_geometry_is_reported_and_never_silently_skipped(tmp_path):
    value=spec(tmp_path)
    profile=value['contexts'][0]['profiles']['cuda:0']
    profile['kernel_geometry']=[row for row in profile['kernel_geometry'] if row['K']!=1]
    profile['kernel_geometry'].append(copy.deepcopy(profile['kernel_geometry'][0]))
    space=prepare(value)
    assert {row['reason'] for row in space['missing_geometry']}=={'missing','ambiguous'}
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        with pytest.raises(ValueError,match='lacks model coverage'):bounded_trait_plan(**value)


def test_infeasible_full_panel_does_not_require_an_oom_geometry_capture(tmp_path):
    from torchgwas.trait_tiling_model import trait_tiled_memory
    value=spec(tmp_path)
    value['workload']['traits']=33
    value['contexts']=value['contexts'][:1]
    value['bounds']=dict(chunks=[4],trait_blocks=[2,33])
    space=prepare(value)
    assert space['missing_geometry'] and all(row['shape'][2]==33 for row in space['missing_geometry'])
    small,large=space['candidates']
    low=trait_tiled_memory(small)['device_bytes']['cuda:0']
    high=trait_tiled_memory(large)['device_bytes']['cuda:0']
    assert low<high
    value['joint']['device_memory_bytes']['cuda:0']=(low+high)//2
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        plan=bounded_trait_plan(**value)
    assert plan['selected']['trait_block']==2
    assert plan['candidates_feasible']==1
    assert plan['rejected'][0]['reason']=='device_memory'


def test_rejects_file_changes_between_census_passes(tmp_path):
    import os
    from pathlib import Path
    value=spec(tmp_path)
    def changing(path,*args,**kwargs):
        result=census(path,*args,**kwargs)
        stat=Path(path).stat()
        os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns+1))
        return result
    with patch('torchgwas.trait_candidate_space.census',side_effect=changing):
        with pytest.raises(ValueError,match='changed during'):prepare(value)


def test_generated_space_uses_the_same_detailed_graph(tmp_path):
    from torchgwas.trait_tiling_plan import detailed_trait_plan
    value=spec(tmp_path)
    space=prepare(value)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        expected=detailed_trait_plan(space['candidates'],**value['joint'])
        actual=bounded_trait_plan(**value)
    assert actual['selected']==expected['selected']
    assert actual['selected']['api_kwargs']['trait_block'] in [2,5]
    assert actual['search_space']['census_passes']==1
    assert actual['search_space']['census_views']==2
    assert actual['selection_validated'] is False


def test_cli_accepts_bounded_trait_space(tmp_path,capsys):
    from torchgwas.autotune import main
    value=spec(tmp_path);value['model']='detailed_traits_space'
    path=tmp_path/'space.json';path.write_text(json.dumps(value))
    with patch('sys.argv',['autotune',str(path)]),patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        main()
    result=json.loads(capsys.readouterr().out)
    assert result['candidates_evaluated']==6
    assert result['search_space']['raw_candidates']==8
