"""Remaining-work scenarios expose their assumptions and all tuning costs."""
from copy import deepcopy
import json
from unittest.mock import patch
import numpy as np
import pytest

from torchgwas.window_forecast import forecast_remaining_windows,remaining_forecast_payback
from torchgwas.window_model import compare_prepared_windows
from torchgwas.pgen_work_bounds import PgenHeaderWork
from test_pgen_native_reader import write_pgen
from test_trait_tiling_model import candidate
from test_mechanistic_shapes import component


ASSUMPTIONS=dict(source_work='Uniform source primitive work beyond the sampled prefix.',
    output_occupancy='Dense output has fixed pair coverage.',resource_capacity='Bound independent capacities remain available.',
    partition_balance='Remaining work has the modeled per-partition proportions.')


def series():
    configuration=dict(partition_axis='trait',windows=[dict(device='cuda:0',trait_range=[0,2],variant_start=0,
        issued_chunks=4,chunk_markers=2,profile_sha256='toy-independent-profile')])
    contract=dict(schema='torchgwas.prepared_comparison.v1',model_identity={'source':'toy-graph-v1'},
        reduction=None,output=dict(block_bytes=64),baseline=configuration,candidate=deepcopy(configuration))
    result=[]
    for pairs in (10,20,30):
        values={}
        for name,seconds in [('baseline',2+.1*pairs),('candidate',3+.05*pairs)]:
            values[name]=dict(estimated_window_seconds=seconds,unpriced_terms=[],windows=[dict(device='cuda:0',
                trait_range=[0,2],variant_range=[0,pairs//2],issued_chunks=4,writer_streams=dict(t=dict(block_bytes=64,
                    payload_bytes=4*pairs,writeback_interval_bytes=64<<20,writeback_submit_calls=0)))])
        result.append(dict(comparison_contract=deepcopy(contract),coverage=[dict(variant_range=[0,pairs//2],trait_ranges=[[0,2]])],
            **values,calculation_cpu_seconds=.01,calculation_wall_seconds=.02))
    return result


def forecast(rows=None,**kw):
    options=dict(remaining_pairs=100,boundary_adjustments=dict(baseline=[0.,0.],candidate=[0.,0.]),
        relative_model_error=0.,max_slope_change=.1,max_extrapolation=100.,assumptions=ASSUMPTIONS)
    options.update(kw)
    return forecast_remaining_windows(series() if rows is None else rows,**options)


def test_affine_pipeline_extrapolation_preserves_one_fill_and_drain():
    rows=series();original=deepcopy(rows);result=forecast(rows)
    assert rows==original and result['status']=='stable_scenario'
    assert result['forecasts']['baseline']['lower_seconds']==pytest.approx(12.)
    assert result['forecasts']['candidate']['upper_seconds']==pytest.approx(8.)
    assert result['scenario_gain_floor_seconds']==pytest.approx(4.)
    assert result['calculation_cpu_seconds']==pytest.approx(.03)
    assert result['calculation_wall_seconds']==pytest.approx(.06)
    assert not result['selection_validated']
    # Repeating an isolated window would incorrectly pay fill/drain each time.
    assert result['forecasts']['baseline']['upper_seconds']<rows[-1]['baseline']['estimated_window_seconds']*100/30


def test_boundary_adjustment_and_declared_model_error_remain_visible():
    result=forecast(relative_model_error=.1,boundary_adjustments=dict(baseline=[-2.,-1.],candidate=[1.,2.]))
    assert result['forecasts']['baseline']['lower_seconds']==pytest.approx(8.8)
    assert result['forecasts']['baseline']['upper_seconds']==pytest.approx(12.2)
    assert result['forecasts']['candidate']['lower_seconds']==pytest.approx(8.2)
    assert result['forecasts']['candidate']['upper_seconds']==pytest.approx(10.8)
    assert result['scenario_gain_floor_seconds']==pytest.approx(-2.)


def test_publication_and_reserve_can_make_a_window_win_unprofitable():
    prediction=forecast()
    yes=remaining_forecast_payback(prediction,planning_seconds=.5,switching_seconds=.25,publication_seconds=.5,reserve_seconds=.25)
    assert yes['worthwhile'] and yes['net_scenario_gain_seconds']==pytest.approx(2.5)
    no=remaining_forecast_payback(prediction,planning_seconds=.5,switching_seconds=.25,publication_seconds=3.,reserve_seconds=.25)
    assert not no['worthwhile'] and no['reason']=='gain_does_not_repay_tuning'
    assert not yes['selection_validated']


def test_unstable_marginal_work_cannot_pass_the_payback_prerequisite():
    rows=series()
    for row,seconds in zip(rows,[10.,11.,20.]):row['baseline']['estimated_window_seconds']=seconds
    result=forecast(rows)
    assert result['status']=='unstable_marginal_cost' and not result['forecasts']['baseline']['stable']
    gate=remaining_forecast_payback(result,planning_seconds=0.,switching_seconds=0.,publication_seconds=0.,reserve_seconds=0.)
    assert not gate['worthwhile'] and gate['reason']=='unstable_marginal_cost'


def test_large_remaining_extent_does_not_expand_source_or_pair_arrays():
    result=forecast(remaining_pairs=10**15,max_extrapolation=10**15)
    assert result['remaining_pairs']==10**15 and result['horizon_pairs']==[10,20,30]
    assert result['forecasts']['baseline']['upper_seconds']==pytest.approx(2+.1*10**15)
    assert len(result['forecasts']['baseline']['marginal_seconds_per_pair'])==2
    assert result['status']=='unmodeled_writer_regime'
    assert len(result['writer_regime_warnings'])==2
    gate=remaining_forecast_payback(result,planning_seconds=0.,switching_seconds=0.,publication_seconds=0.,reserve_seconds=0.)
    assert not gate['worthwhile'] and gate['reason']=='unmodeled_writer_regime'


def test_a_sampled_periodic_writeback_cycle_does_not_trigger_unseen_regime_warning():
    rows=series()
    for report in rows:
        for name in ('baseline','candidate'):
            stream=report[name]['windows'][0]['writer_streams']['t'];stream['writeback_interval_bytes']=48
            stream['writeback_submit_calls']=stream['payload_bytes']//48
    result=forecast(rows,remaining_pairs=100)
    assert result['status']=='stable_scenario' and result['writer_regime_warnings']==[]


@pytest.mark.parametrize('mutation',[
    lambda r:r[1]['comparison_contract'].update(model_identity={'source':'changed'}),
    lambda r:r[1]['comparison_contract']['output'].update(block_bytes=32),
    lambda r:r[1]['comparison_contract']['candidate']['windows'][0].update(profile_sha256='changed-rate'),
    lambda r:r[1]['comparison_contract']['candidate']['windows'][0].update(issued_chunks=0),
    lambda r:r[1]['candidate']['windows'][0].update(issued_chunks=0),
    lambda r:r[1]['candidate']['windows'][0].update(device='cuda:1'),
    lambda r:r[1]['coverage'][0]['variant_range'].__setitem__(1,11),
    lambda r:r[1]['baseline'].update(estimated_window_seconds=100.),
    lambda r:r.reverse()])
def test_mismatched_bindings_coverage_and_nonpositive_slopes_refuse_forecasts(mutation):
    rows=series();mutation(rows)
    with pytest.raises(ValueError):forecast(rows)


def test_unknown_identity_and_unresolved_dense_writer_configuration_are_rejected():
    for change in ('identity','writer'):
        rows=series()
        for row in rows:
            if change=='identity':row['comparison_contract']['model_identity']=None
            else:row['comparison_contract']['output']['block_bytes']=None
        with pytest.raises(ValueError,match='identity' if change=='identity' else 'writer block_bytes'):forecast(rows)


def test_changed_writer_configuration_is_not_hidden_by_a_stable_slope():
    rows=series();rows[1]['baseline']['windows'][0]['writer_streams']['t']['block_bytes']=32
    with pytest.raises(ValueError,match='Writer configuration changed'):forecast(rows)


def test_unequal_partition_growth_is_rejected_even_with_matching_pair_totals():
    rows=series()
    for index,row in enumerate(rows):
        # Two disjoint one-trait windows. Total coverage is still 10,20,30
        # pairs, but the middle horizon moves work to the first partition.
        counts=[5,5] if index==0 else [12,8] if index==1 else [15,15]
        row['coverage']=[dict(variant_range=[0,min(counts)],trait_ranges=[[0,2]])]
        if counts[0]!=counts[1]:row['coverage'].append(dict(variant_range=[min(counts),max(counts)],trait_ranges=[[0,1]]))
        for name in ('baseline','candidate'):
            bound=row['comparison_contract'][name]['windows'][0]
            row['comparison_contract'][name]['windows']=[dict(bound,device='cuda:'+str(i),trait_range=[i,i+1]) for i in (0,1)]
            window=row[name]['windows'][0]
            row[name]['windows']=[dict(deepcopy(window),device='cuda:'+str(i),trait_range=[i,i+1],variant_range=[0,count])
                for i,count in enumerate(counts)]
    with pytest.raises(ValueError,match='preserve partition balance'):forecast(rows)


@pytest.mark.parametrize('options',[dict(remaining_pairs=29),dict(max_extrapolation=2.),dict(relative_model_error=1.),
    dict(relative_model_error=True),dict(max_slope_change=-1.),dict(assumptions={}),dict(boundary_adjustments={}),
    dict(boundary_adjustments=dict(baseline=[1.,0.],candidate=[0.,0.])),dict(remaining_pairs=10**400,max_extrapolation=1e300)])
def test_invalid_or_excessive_forecast_requests_fail(options):
    with pytest.raises(ValueError):forecast(**options)


def test_actual_header_window_models_carry_consistent_horizon_contracts(tmp_path):
    path=tmp_path/'input.pgen';write_pgen(path,(np.arange(16*32).reshape(16,32)%3).astype(np.uint8))
    choice=candidate(path,traits=2,width=2,count=1,block_bytes=32);header=PgenHeaderWork(path);rows=[]
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        for count in (4,8,12):
            layouts=[]
            for size in (2,4):
                window=deepcopy(choice['tiles'][0]);window['issued_chunks']=4
                window['data'].update(markers=count,encoded=header.window(0,count,size))
                window['profile']['chunk_markers']=size
                layouts.append(dict(windows=[window],partition_axis='trait'))
            rows.append(compare_prepared_windows(*layouts,total_traits=2,reduction=None,output=choice['output'],
                shared_capacities=choice['shared_capacities'],endpoint='upper',host_serial_fraction=.5,
                model_identity=dict(source='test-fixture',independent_prices='synthetic controls',runtime=('python','torch'))))
    assert rows[0]['comparison_contract']==rows[1]['comparison_contract']==rows[2]['comparison_contract']
    rows[1]=json.loads(json.dumps(rows[1]))
    assert rows[0]['comparison_contract']==rows[1]['comparison_contract']
    result=forecast(rows,remaining_pairs=32,max_slope_change=.9)
    assert result['horizon_pairs']==[8,16,24] and result['extrapolation_ratio']==pytest.approx(4/3)
    assert all(r['lower_seconds']<=r['upper_seconds'] for r in result['forecasts'].values())
