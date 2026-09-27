import copy
import pytest
from torchgwas.gil_service import gil_service_prices


def probe():
    controls={name:dict(cpu_seconds=.02,detached_cpu_seconds=.019 if count else 0.,
                        detached_intervals=count,balance_errors=0) for name,count in [('native',1),('nested',3),('python',0)]}
    meter_rows=[dict(repeat=r,recording=on,intervals=20000,detached_intervals=20000 if on else 0,
                     cpu_seconds=.04 if on else .02,detached_cpu_seconds=.01 if on else 0.,balance_errors=0)
                for r in range(3) for on in [False,True]]
    rows=[dict(primitive='gemm',repeat=r,recording=on,sample_count=100,cpu_seconds_per_call=10e-6,
               raw_cpu_seconds_per_call=11e-6,raw_cpu_sum_seconds=.0011,empty_cpu_seconds_per_call=1e-6,
               detached_cpu_seconds_per_call=6e-6 if on else 0.,detached_intervals=200 if on else 0,balance_errors=0)
          for r in range(3) for on in [False,True]]
    return dict(context_verified=True,torch_version='test',affinity=[1,2],results=[dict(mode='single1',workers=[dict(
        device=1,context_verified=True,controls=controls,meter_rows=meter_rows,rows=rows)])])


def load(value):
    return gil_service_prices(value,mode='single1',device=1,torch_version='test',cpu_affinity=[1,2])


def test_meter_correction_conserves_total_and_detached_cpu():
    result=load(probe())
    assert result['cpu_primitives']['gemm']==pytest.approx(8e-6)
    assert result['measurement_controls_passed']
    assert not result['total_cpu_transfer_qualified']
    assert not result['serial_transfer_qualified']
    assert not result['prediction_complete']
    assert result['serial_primitives']['gemm']==pytest.approx(3e-6)
    assert result['meter_cpu_seconds_per_interval']==pytest.approx(1e-6)
    assert result['meter_detached_cpu_seconds_per_interval']==pytest.approx(.5e-6)


@pytest.mark.parametrize('change',['context','version','affinity','control','nested','meter_count',
    'duplicate_meter','missing_pair','raw_sum','mean','unbalanced','bad_correction','duplicate_api'])
def test_corrupt_or_unverified_evidence_rejected(change):
    value=probe();worker=value['results'][0]['workers'][0]
    if change=='context':value['context_verified']=False
    if change=='version':value['torch_version']='different'
    if change=='affinity':value['affinity']=[1]
    if change=='control':worker['controls']['native']['detached_cpu_seconds']=0.
    if change=='nested':worker['controls']['nested']['detached_intervals']=1
    if change=='meter_count':worker['meter_rows'][1]['detached_intervals']=19999
    if change=='duplicate_meter':worker['meter_rows'].append(worker['meter_rows'][0])
    if change=='missing_pair':worker['rows'].pop()
    if change=='raw_sum':worker['rows'][0]['raw_cpu_sum_seconds']=0.
    if change=='mean':worker['rows'][0]['cpu_seconds_per_call']=1e-6
    if change=='unbalanced':worker['rows'][1]['balance_errors']=1
    if change=='bad_correction':worker['rows'][1]['detached_cpu_seconds_per_call']=10e-6
    if change=='duplicate_api':worker['rows'].append(worker['rows'][0])
    with pytest.raises(ValueError):load(value)
