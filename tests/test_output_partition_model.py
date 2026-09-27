"""Both output partition axes use the same component prices and execution graph."""
import copy
from unittest.mock import patch

import pytest

from test_trait_candidate_space import spec
from test_trait_tiling_model import component
from torchgwas.trait_candidate_space import prepare_trait_candidates,bounded_trait_plan
from torchgwas.trait_tiling_model import trait_tiled_shape,trait_tiled_memory,torch_trait_tiled_runtime


def candidates(value):
    return prepare_trait_candidates(value['workload'],value['contexts'],output=value['output'],
        partition_axes=['trait','variant'],**value['bounds'])


def test_joint_axes_preserve_trait_proposals_and_exact_variant_censuses(tmp_path):
    value=spec(tmp_path);value['workload']['phenotype_c_contiguous']=True
    old=prepare_trait_candidates(value['workload'],value['contexts'],output=value['output'],**value['bounds'])
    before=copy.deepcopy(value);mixed=candidates(value)
    assert value==before and mixed['candidates'][:6]==old['candidates']
    assert mixed['raw_candidates']==12 and len(mixed['candidates'])==10
    assert mixed['census_passes']==1 and mixed['census_views']==6 and mixed['census_chunks']==16
    assert not mixed['missing_geometry']
    for candidate in mixed['candidates'][6:]:
        shape=trait_tiled_shape(candidate)
        assert shape['partition_axis']=='variant' and shape['traits']==5 and shape['markers']==10
        cursor=0
        for tile in candidate['tiles']:
            lo,hi=tile['variant_range'];data=tile['data'];encoded=data['encoded']
            assert lo==cursor and data['markers']==hi-lo and encoded['variant_range']==[lo,hi]
            assert encoded['file_markers']==10 and data['traits_analyzed']==5
            assert data['phenotype_c_contiguous'] is True
            cursor=hi
        assert cursor==10


def test_single_device_full_panel_has_equal_priced_work_on_either_axis(tmp_path):
    value=spec(tmp_path);value['workload']['phenotype_c_contiguous']=True
    mixed=candidates(value)['candidates']
    trait=next(c for c in mixed if c.get('partition_axis') is None and c['trait_block']==5 and len(c['devices'])==1)
    variant=next(c for c in mixed if c.get('partition_axis')=='variant' and len(c['devices'])==1)
    assert trait_tiled_memory(trait)==trait_tiled_memory(variant)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        for fraction in (0.,1.):
            for policy in ('held-first','held-last'):
                t=torch_trait_tiled_runtime(trait,host_serial_fraction=fraction,host_serial_policy=policy)
                v=torch_trait_tiled_runtime(variant,host_serial_fraction=fraction,host_serial_policy=policy)
                assert t['estimated_tile_seconds']==v['estimated_tile_seconds']
                assert t['binary_payload_bytes']==v['binary_payload_bytes']==440
                assert v['genotype_passes']==1


def test_variant_graph_conserves_shared_work_and_does_not_reread_payload(tmp_path):
    value=spec(tmp_path);mixed=candidates(value)['candidates']
    candidate=next(c for c in mixed if c.get('partition_axis')=='variant' and len(c['devices'])==2)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        graph=torch_trait_tiled_runtime(candidate,host_serial_fraction=.5,host_serial_policy='held-last',return_graph=True)
        result=graph.solve();report=torch_trait_tiled_runtime(candidate,host_serial_fraction=.5)
    assert report['binary_payload_bytes']==440 and report['df_payload_bytes']==40
    assert report['genotype_passes']==1 and sum(r['blocks'] for r in report['tiles'])==5
    for resource,capacity in graph.capacities.items():
        work=sum(graph.nodes[name][0]*demand.get(resource,0.) for name,demand in graph.demands.items())
        assert result['seconds']+1e-10>=work/capacity
    assert result['start']['tile1:prepare:covariate_basis']<result['end']['tile0:complete']


def test_joint_planner_emits_executable_variant_settings_and_hard_bounds(tmp_path):
    value=spec(tmp_path);value['bounds']['partition_axes']=['trait','variant']
    value['joint']['shortlist_size']=20
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        plan=bounded_trait_plan(value['workload'],value['contexts'],bounds=value['bounds'],joint=value['joint'],output=value['output'])
    rows=[r for r in plan['shortlist'] if r.get('partition_axis')=='variant']
    assert len(rows)==4
    for row in rows:
        assert row['genotype_passes']==1 and row['trait_block'] is None
        assert row['api_kwargs']['variant_devices']==row['devices']
        assert not {'trait_block','trait_devices'}&set(row['api_kwargs'])
    for limit,value_limit in [('max_candidates',11),('max_candidate_tiles',21),('max_census_chunks',15)]:
        with patch('torchgwas.trait_candidate_space.census',side_effect=AssertionError('must refuse before census')):
            with pytest.raises(ValueError,match=limit):
                prepare_trait_candidates(value['workload'],value['contexts'],output=value['output'],
                    **dict(value['bounds'],**{limit:value_limit}))


@pytest.mark.parametrize('change',[
    lambda c:c['tiles'][1].update(variant_range=[5,10]),
    lambda c:c['tiles'][0].update(trait_range=[1,6]),
    lambda c:c['tiles'][0]['data']['encoded'].update(variant_range=[0,10]),
    lambda c:c['tiles'][0]['data'].update(markers=3),
    lambda c:c['tiles'][0]['profile'].update(decode_workers=3),
])
def test_variant_candidates_cannot_drift_from_executor_partition(tmp_path,change):
    value=spec(tmp_path);c=next(c for c in candidates(value)['candidates'] if c.get('partition_axis')=='variant' and len(c['devices'])==2)
    change(c)
    with pytest.raises(ValueError):trait_tiled_shape(c)
