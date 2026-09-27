"""Live-range binding and ragged multi-device remaining-work scenarios."""
from copy import deepcopy
import time
import pytest

from torchgwas.productive_forecast import productive_window_proposal
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.sumstats_indexed import IndexedChunkWrite


OPTIONS=dict(boundary_adjustments=dict(baseline=[0.,0.],candidate=[0.,0.]),
    relative_model_error=0.,max_slope_change=.1,max_extrapolation=8.,
    assumptions=dict(source_work='Declared uniform source-work scenario.',
        output_occupancy='Dense cells stay fixed.',resource_capacity='Independent supplied prices stay valid.',
        partition_balance='Use the balanced core/envelope for exact unequal remaining extents.'))


def run_fixture():
    run=ProductiveTuningRun([dict(id=str(i),device='cuda:'+str(i),variant_range=[100,164],trait_range=[i,i+1]) for i in (0,1)],
        chunk_sizes=[2,4],initial=2,budget=IncrementalPlanningBudget(max_steps=3,max_cpu_seconds=10.,max_window_seconds=60.))
    run.for_partition('0')(100,164,4);run.for_partition('0')(102,164,4)
    run.for_partition('1')(100,164,4)
    now=time.perf_counter();run.output_written(IndexedChunkWrite(100,102,'significant',0,0,None,now,now,False))
    return run


def series(snapshot):
    windows=[dict(device=p['device'],trait_range=p['trait_range'],variant_start=p['cursor'],
        issued_chunks=p['issued_chunks'],chunk_markers=2,profile_sha256='all-profile',
        fixed_data_sha256='fixed-data',chunk_invariant_profile_sha256='fixed-profile')
        for p in snapshot['partitions'] if p['cursor']<p['variant_range'][1]]
    config=dict(partition_axis='trait',windows=windows);alternative=deepcopy(config)
    for w in alternative['windows']:w['chunk_markers']=4;w['profile_sha256']='different-chunk-geometry'
    contract=dict(schema='torchgwas.prepared_comparison.v1',model_identity={'source':'synthetic-independent-controls'},
        total_traits=2,reduction=None,output=dict(block_bytes=64),baseline=config,candidate=alternative)
    from torchgwas.window_model import _coverage
    reports=[]
    for width in (4,8,12):
        values={};pairs=width*len(windows)
        for name,rate in [('baseline',.1),('candidate',.05)]:
            rows=[dict(device=w['device'],trait_range=w['trait_range'],variant_range=[w['variant_start'],w['variant_start']+width],
                issued_chunks=w['issued_chunks'],writer_streams=dict(t=dict(block_bytes=64,payload_bytes=4*width,
                    writeback_interval_bytes=64<<20,writeback_submit_calls=0))) for w in windows]
            values[name]=dict(estimated_window_seconds=1+rate*pairs,unpriced_terms=[],windows=rows)
        coverage=_coverage([dict(trait_range=w['trait_range'],data=dict(encoded=dict(variant_range=w['variant_range']))) for w in rows])
        reports.append(dict(comparison_contract=deepcopy(contract),coverage=coverage,**values,
            calculation_cpu_seconds=.001,calculation_wall_seconds=.001))
    return reports


def propose(state,rows=None,**options):
    return productive_window_proposal(state,series(state) if rows is None else rows,chunk_sizes=[2,4],**dict(OPTIONS,**options))


def test_exact_nonzero_unissued_extents_include_prefetch_and_unequal_gpu_progress():
    state=run_fixture().snapshot();before=deepcopy(state);proposal,audit=propose(state)
    assert state==before and audit['remaining_pairs']==122
    assert [r['variant_range'] for r in audit['remaining_partitions']]==[[104,164],[102,164]]
    assert audit['partition_extent_scenarios']['baseline']==dict(lower_pairs=120,upper_pairs=124,lower_ratio=5.,upper_ratio=124/24)
    assert proposal==dict(chunk_size=4,baseline_seconds=pytest.approx(13.),candidate_seconds=pytest.approx(7.2))
    assert not audit['selection_validated']


def test_forecast_step_applies_only_to_future_reservations_and_retains_audit():
    run=run_fixture();before=run.snapshot();result=run.forecast_step(series,forecast_options=OPTIONS,
        remaining_seconds=100.,expected_cpu_seconds=.001,expected_wall_seconds=.001,expected_gain_seconds=5.)
    assert result['applied'] and result['forecast_audit']['issued_revision']==3
    assert run.snapshot()['partitions']==before['partitions']
    assert run.snapshot()['decisions'][-1]['forecast_audit']==result['forecast_audit']
    assert run.for_partition('0')(104,164,4)==4
    assert run.for_partition('1')(102,164,4)==4
    run.finish(successful=False)


def test_no_comparison_is_built_before_productive_output():
    run=ProductiveTuningRun([dict(id='0',device='cuda:0',variant_range=[0,100],trait_range=[0,2])],chunk_sizes=[2,4],initial=2)
    result=run.forecast_step(lambda snapshot:pytest.fail('upfront model construction'),forecast_options=OPTIONS,
        remaining_seconds=100.,expected_cpu_seconds=.001,expected_wall_seconds=.001)
    assert not result['evaluated'] and result['reason']=='no_useful_output_yet'
    run.finish(successful=False)


def test_jagwas_refuses_partial_phenotype_panel():
    state=run_fixture().snapshot();state['partitions']=state['partitions'][:1];state['issued_revision']=2
    rows=series(state)
    for report in rows:
        contract=report['comparison_contract'];contract['reduction']='jagwas'
        for name in ('baseline','candidate'):contract[name]['partition_axis']='variant'
    with pytest.raises(ValueError,match='complete phenotype'):propose(state,rows)


def test_changed_frontier_refuses_old_comparisons_and_scientific_scan_continues():
    run=run_fixture();old=series(run.snapshot());run.for_partition('0')(104,164,4)
    result=run.forecast_step(lambda snapshot:old,forecast_options=OPTIONS,remaining_seconds=100.,
        expected_cpu_seconds=.001,expected_wall_seconds=.001)
    assert result['error']=='ValueError' and not result['applied']
    assert run.snapshot()['stop_reason']=='planning_error'
    assert run.for_partition('0')(106,164,4)==2
    run.finish(successful=False)


@pytest.mark.parametrize('fault',['missing_partition','cursor','ordinal','current_size','unadmitted','mixed_size','device',
    'data','prices','axis','past_end','revision','overlap','gap','unadmitted_prefix','budget','future_horizon'])
def test_bad_model_or_prefix_cannot_authorize_an_action(fault):
    state=run_fixture().snapshot();rows=series(state);options={}
    if fault=='missing_partition':
        for r in rows:r['comparison_contract']['candidate']['windows'].pop()
    elif fault in ('cursor','ordinal','current_size','unadmitted','mixed_size','device','data','prices'):
        name='baseline' if fault=='current_size' else 'candidate'
        field,value={'cursor':('variant_start',100),'ordinal':('issued_chunks',0),'current_size':('chunk_markers',4),
            'unadmitted':('chunk_markers',3),'mixed_size':('chunk_markers',2),'device':('device','cuda:2'),
            'data':('fixed_data_sha256','changed'),'prices':('chunk_invariant_profile_sha256','changed')}[fault]
        for r in rows:r['comparison_contract'][name]['windows'][0][field]=value
    elif fault=='axis':
        for r in rows:r['comparison_contract']['candidate']['partition_axis']='variant'
    elif fault=='past_end':state['partitions'][0]['cursor']=165
    elif fault=='revision':state['issued_revision']+=1
    elif fault=='overlap':state['partitions'][1]['trait_range']=[0,1]
    elif fault=='gap':state['partitions'][0]['ranges'][1][0]+=1
    elif fault=='unadmitted_prefix':state['partitions'][0].update(ranges=[[100,103]],cursor=103,issued_chunks=1);state['issued_revision']=2
    elif fault=='budget':options['max_partitions']=1
    else:
        for p in state['partitions']:p['variant_range'][1]=p['cursor']+10
    with pytest.raises(ValueError):propose(state,rows,**options)


def test_extrapolation_limit_applies_to_longest_partition_not_average_extent():
    state=run_fixture().snapshot()
    # Mean ratio is 122/24, below this limit; the envelope is 124/24, above it.
    with pytest.raises(ValueError,match='extrapolation limit'):propose(state,max_extrapolation=123/24)


def test_stable_but_unprofitable_and_unstable_forecasts_keep_the_current_size():
    run=run_fixture()
    costly=run.forecast_step(series,forecast_options=OPTIONS,remaining_seconds=100.,
        expected_cpu_seconds=.001,expected_wall_seconds=.001,publication_seconds=10.)
    assert costly['evaluated'] and not costly['applied']
    def unstable(state):
        rows=series(state);rows[-1]['candidate']['estimated_window_seconds']+=100.;return rows
    result=run.forecast_step(unstable,forecast_options=OPTIONS,remaining_seconds=100.,
        expected_cpu_seconds=.001,expected_wall_seconds=.001)
    assert result['evaluated'] and not result['applied']
    assert result['forecast_audit']['forecast_status']=='unstable_marginal_cost'
    assert run.snapshot()['current_chunk_size']==2
    run.finish(successful=False)


def test_fully_issued_partition_is_excluded_without_claiming_its_output_is_complete():
    run=run_fixture();control=run.for_partition('0')
    for start in range(104,164,2):control(start,164,4)
    state=run.snapshot();proposal,audit=propose(state)
    assert audit['remaining_pairs']==62 and [r['id'] for r in audit['remaining_partitions']]==['1']
    assert state['written_events']==1 and state['written_rows']==0
    assert proposal['chunk_size']==4
    run.finish(successful=False)


@pytest.mark.parametrize('changed',['decode_workers','fsync_seconds'])
def test_source_generated_profile_binding_detects_settings_and_price_changes(tmp_path,changed):
    import numpy as np
    from unittest.mock import patch
    from test_pgen_native_reader import write_pgen
    from test_trait_tiling_model import candidate
    from test_mechanistic_shapes import component
    from torchgwas.pgen_work_bounds import PgenHeaderWork
    from torchgwas.window_model import compare_prepared_windows
    path=tmp_path/'input.pgen';write_pgen(path,(np.arange(32*32).reshape(32,32)%3).astype(np.uint8))
    choice=candidate(path,traits=2,width=2,count=1,block_bytes=32);header=PgenHeaderWork(path);layouts=[]
    for size in (2,4):
        w=deepcopy(choice['tiles'][0]);w['issued_chunks']=2;w['profile']['chunk_markers']=size
        w['data'].update(markers=8,encoded=header.window(4,12,size))
        if size==4:w['profile'][changed]=1 if changed=='decode_workers' else .001
        layouts.append(dict(windows=[w],partition_axis='trait'))
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        comparison=compare_prepared_windows(*layouts,total_traits=2,reduction=None,output=choice['output'],
            shared_capacities=choice['shared_capacities'],endpoint='upper',host_serial_fraction=.5,
            model_identity=dict(source='unit source graph',prices='synthetic controls'))
    state=dict(prefix_complete=True,finished=None,first_written=1.,issued_revision=2,current_chunk_size=2,
        partitions=[dict(id='0',device='cuda:0',variant_range=[0,32],trait_range=[0,2],ranges=[[0,2],[2,4]],cursor=4,issued_chunks=2)])
    with pytest.raises(ValueError,match='fixed data, prices or execution settings'):
        propose(state,[comparison]*3)
