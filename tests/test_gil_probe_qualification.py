import copy
import pytest
from test_gil_service import probe
from torchgwas.gil_probe_qualification import compare_gil_probe_bracket


def evidence():
    audited=probe()
    audited.update(host='same-host',devices={'1':dict(name='GPU',capability=[8,0])},
        primitive_bank='jagwas',fixed_shape=[32,32],source_sha256={
            'benchmarks/direct_jagwas_host_primitives.py':'bank',
            'benchmarks/direct_calculator_gil_audit.c':'C',
            '/build/libtorchgwas_gil_audit.so':'lib'})
    before=copy.deepcopy(audited);before['ld_audit']=None
    worker=before['results'][0]['workers'][0];worker['complete']=True
    for row in worker['rows']:
        row['phase']=int(row.pop('recording'))
        row.update(cpu_seconds_per_call=8e-6,raw_cpu_seconds_per_call=9e-6,raw_cpu_sum_seconds=.0009)
    return before,audited,copy.deepcopy(before)


def compare(values):
    return compare_gil_probe_bracket(*values,maximum_relative_difference=.1)


def scale(probe, scale):
    for row in probe['results'][0]['workers'][0]['rows']:
        row['cpu_seconds_per_call']*=scale
        row['raw_cpu_seconds_per_call']=row['cpu_seconds_per_call']+row['empty_cpu_seconds_per_call']
        row['raw_cpu_sum_seconds']=row['raw_cpu_seconds_per_call']*row['sample_count']


def test_matched_total_cpu_does_not_certify_serial_fraction_or_prediction():
    result=compare(evidence())
    assert result['total_cpu_compatible']
    assert not result['serial_transfer_qualified']
    assert not result['prediction_complete']
    row=result['contexts']['single1:1']['primitives']['gemm']
    assert row['unhooked_cpu_seconds']==pytest.approx(8e-6)
    assert row['audit_relative_delta_interval']==pytest.approx([0,0])
    assert row['disabled_relative_delta_interval']==pytest.approx([.25,.25])


def test_a_broad_drift_envelope_is_not_acceptance():
    values=evidence();scale(values[0],.5);scale(values[2],1.5)
    row=compare(values)['contexts']['single1:1']['primitives']['gemm']
    assert row['audit_relative_delta_interval'][0]<0<row['audit_relative_delta_interval'][1]
    assert not row['bracket_stable']
    assert not row['total_cpu_compatible']


def test_stable_but_perturbed_audit_is_rejected():
    values=evidence();scale(values[0],2);scale(values[2],2)
    row=compare(values)['contexts']['single1:1']['primitives']['gemm']
    assert row['bracket_stable']
    assert not row['total_cpu_compatible']


@pytest.mark.parametrize('change',['host','affinity','devices','bank','source','driver','hooked','incomplete',
    'duplicate','raw','negative','primitive','repeat','sample_count','audited_controls'])
def test_mismatched_or_corrupt_controls_rejected(change):
    values=evidence();before,audited,after=values;worker=before['results'][0]['workers'][0]
    if change=='host':before['host']='different'
    if change=='affinity':before['affinity']=[1]
    if change=='devices':before['devices']['1']['name']='different'
    if change=='bank':before['primitive_bank']='scan'
    if change=='source':before['source_sha256']['benchmarks/direct_jagwas_host_primitives.py']='different'
    if change=='driver':before['source_sha256']['extra_driver']='new'
    if change=='hooked':before['ld_audit']='lib'
    if change=='incomplete':worker['complete']=False
    if change=='duplicate':worker['rows'].append(worker['rows'][0])
    if change=='raw':worker['rows'][0]['raw_cpu_sum_seconds']=0
    if change=='negative':worker['rows'][0]['cpu_seconds_per_call']=-1
    if change=='primitive':worker['rows'][0]['primitive']='unknown'
    if change=='repeat':worker['rows'].pop()
    if change=='sample_count':
        worker['rows'][0]['sample_count']=200
        worker['rows'][0]['raw_cpu_sum_seconds']*=2
    if change=='audited_controls':audited['context_verified']=False
    with pytest.raises(ValueError):compare(values)


@pytest.mark.parametrize('tolerance',[True,-1,float('inf'),float('nan')])
def test_explicit_finite_tolerance_required(tolerance):
    with pytest.raises(ValueError):compare_gil_probe_bracket(*evidence(),maximum_relative_difference=tolerance)


def test_zero_after_empty_control_is_unresolved_not_an_error_or_free_api():
    values=evidence();scale(values[0],0);scale(values[2],0)
    row=compare(values)['contexts']['single1:1']['primitives']['gemm']
    assert row['resolution_limited']
    assert row['audit_relative_delta_interval'] is None
    assert row['bracket_relative_drift'] is None
    assert not row['total_cpu_compatible']
