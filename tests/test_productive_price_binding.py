"""Expired or substituted component evidence cannot authorize a live switch."""
from copy import deepcopy
import time
import pytest

from test_price_binding import evidence
from test_productive_forecast import OPTIONS,series
from torchgwas.calibration_cache import _digest
from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.sumstats_indexed import IndexedChunkWrite


def run():
    value=ProductiveTuningRun([dict(id='0',device='cuda:0',variant_range=[100,164],trait_range=[0,1])],
        chunk_sizes=[2,4],initial=2,budget=IncrementalPlanningBudget(max_steps=2,max_cpu_seconds=10.,max_window_seconds=60.))
    value.for_partition('0')(100,164,4);value.for_partition('0')(102,164,4)
    now=time.perf_counter();value.output_written(IndexedChunkWrite(100,102,'significant',0,0,None,now,now,False))
    return value


def bound_series(state,profile):
    rows=series(state)
    fixed=profile['contexts'][0]['profiles']['cuda:0']
    for row in rows:
        contract=row['comparison_contract']
        contract['model_identity']=dict(source_sha256=profile['source_sha256'],price_profile_sha256=_digest(profile))
        for name in ('baseline','candidate'):
            for window in contract[name]['windows']:window['chunk_invariant_profile_sha256']=_digest(fixed)
    return rows


def step(value,builder,profile):
    return value.forecast_step(builder,price_profile=profile,forecast_options=OPTIONS,
        remaining_seconds=100.,expected_cpu_seconds=.001,expected_wall_seconds=.001,expected_gain_seconds=2.)


def test_valid_price_record_is_checked_and_audited_before_switch(evidence):
    profile=evidence['bind']();value=run();before=deepcopy(profile)
    result=step(value,lambda state:bound_series(state,profile),profile)
    assert result['applied'] and profile==before
    row=result['forecast_audit']['price_evidence']['bindings'][0]
    assert row['observed_unix_seconds']==90. and row['age_seconds']==10.
    assert row['record_sha256']==evidence['record']['record_sha256']
    assert value.for_partition('0')(104,164,4)==4
    value.finish(successful=False)


@pytest.mark.parametrize('fault',['already_expired','expires_during_model','different_profile','different_prices','different_source','undeclared','variable_field'])
def test_invalid_optional_price_calculation_keeps_valid_scan(evidence,fault):
    profile=evidence['bind']();value=run();calls=[]
    if fault=='already_expired':evidence['now'][0]=120.
    if fault=='undeclared':profile.pop('price_bindings')
    if fault=='variable_field':
        profile['contexts'][0]['profiles']['cuda:0']['chunk_markers']=3.
        profile['price_bindings'][0]['targets'][0]['context_path'][-1]='chunk_markers'
    def build(state):
        calls.append(True);rows=bound_series(state,profile)
        if fault=='expires_during_model':evidence['now'][0]=120.
        if fault=='different_profile':
            for row in rows:row['comparison_contract']['model_identity']['price_profile_sha256']='other'
        if fault=='different_source':
            for row in rows:row['comparison_contract']['model_identity']['source_sha256']={}
        if fault=='different_prices':
            for row in rows:
                for name in ('baseline','candidate'):
                    row['comparison_contract'][name]['windows'][0]['chunk_invariant_profile_sha256']='other'
        return rows
    result=step(value,build,profile)
    assert result['error']=='ValueError' and not result['applied']
    assert bool(calls)==(fault not in ('already_expired','undeclared'))
    assert value.snapshot()['planning']['cpu_seconds']>0
    assert value.for_partition('0')(104,164,4)==2
    value.finish(successful=False)


def test_resident_copy_price_changes_actual_jagwas_finish_service():
    from unittest.mock import patch
    from test_jagwas_scan_work import joint_fixture,scan_component
    from torchgwas.mechanistic_torch import torch_scan_work
    data,profile=joint_fixture();changed=deepcopy(profile)
    changed['owned_result_copy_scenario']['resident_cpu_seconds_per_byte']*=3
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=scan_component):
        before=torch_scan_work(data,profile);after=torch_scan_work(data,changed)
    for block,altered in zip(before['blocks'],after['blocks']):
        old=before['owned_result_work'][block['markers']]['service']['additional_copy_cpu_seconds']
        new=after['owned_result_work'][block['markers']]['service']['additional_copy_cpu_seconds']
        assert old>0 and new==pytest.approx(3*old)
        assert altered['finish_seconds']-block['finish_seconds']==pytest.approx((new-old)/profile['cpu_fraction'])
