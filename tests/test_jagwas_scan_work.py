"""Joint projection integration uses real untimed kernels and synthetic prices."""
import copy
from unittest.mock import patch
import pytest
from torchgwas.execution_graph import torch_scan_schedule
from torchgwas.mechanistic_torch import torch_scan_work,torch_scan_runtime,torch_runtime
from torchgwas.native_control_work import native_control_work
from torchgwas.reduction_tensor_work import append_jagwas_projection
from test_mechanistic_shapes import fixture
from test_jagwas_tensor_service import CAPTURE,prices,resources


def joint_fixture():
    from dataclasses import asdict
    row=next(r for r in CAPTURE['rows'] if r['K']==512 and r['compute_dtype']=='float32')
    n,b,k=row['N'],row['B'],row['K'];m=2*b
    data,profile=fixture(k=k,c=2)
    data.update(samples=n,markers=m)
    data['encoded'].update(samples=n,markers=m,record_form_counts={'0':m},record_payload_bytes=((n+3)//4)*m)
    profile['decode_units']['expand_int8_tail_sample']=1e-9
    gpu=asdict(resources());gpu.pop('host_cpu_fraction')
    profile.update(reduction='jagwas',chunk_markers=b,gpu_resources=gpu,
        kernel_geometry=[dict(N=n,B=b,K=k,C=2,kernels=[])],
        joint_kernel_geometry=[row],joint_host_primitives=prices(row),
        result_finish_service=dict(cpu_seconds=1e-5,serial_cpu_seconds=8e-6,
            baseline_copy_bytes=544,reduction='jagwas',result_arrays=5,baseline_rows=32,
            return_df=False,replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False),
        control_primitives={key:1e-6 for counts in native_control_work(reduction='jagwas').values() for key in counts},
        owned_result_copy_scenario=dict(resident_cpu_seconds_per_byte=1e-10,fresh_cpu_seconds_per_byte=8e-10,fresh_fraction=0.))
    return data,profile


def scan_component(*args,**kwargs):
    return dict(host_dispatch_cpu_seconds=1e-5,estimated_span_seconds=2e-5,kernel_service_seconds=1e-5,
        host_calls=[dict(call_id=0,primitive='scan_api',cpu_seconds=1e-5,submit_finish=1e-5)],
        operations=[dict(op='scan',phase='statistics',host_submit_finish=1e-5,kernel_service_seconds=1e-5)],
        gemm={},unpriced_terms=[],kernel_count=1,logical_bytes=19,modeled_hbm_bytes=23)


def test_native_joint_payload_and_projection_precede_result_transfer():
    data,profile=joint_fixture();snapshot=copy.deepcopy(profile)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=scan_component):
        work=torch_scan_work(data,profile)
        runtime=torch_scan_runtime(data,profile)
    assert profile==snapshot
    assert sum(block['d2h_bytes'] for block in work['blocks'])==17*data['markers']
    assert sum(block['h2d_bytes'] for block in work['blocks'])==data['samples']*data['markers']
    for component in work['components'].values():
        joint=component['reduction_component']
        assert component['kernel_count']==1+joint['kernel_count']
        assert component['operations'][0]['phase']=='statistics'
        assert all(op['phase']=='joint_projection' for op in component['operations'][1:])
        assert component['host_dispatch_cpu_seconds']==pytest.approx(1e-5+joint['host_dispatch_cpu_seconds'])
        assert component['logical_bytes']==19+joint['logical_bytes']
        assert component['modeled_hbm_bytes']==23+joint['modeled_hbm_bytes']
        assert joint['gemm']['useful_flops']==2*profile['chunk_markers']*data['traits_analyzed']**2
    for result in work['owned_result_work'].values():
        assert result['work']['allocation_calls']==5
        assert result['service']['additional_copy_bytes']==17*profile['chunk_markers']-544
    slower=copy.deepcopy(work['blocks'])
    for block in slower:
        block['operations'][-1]['kernel_service_seconds']+=1.
    schedule=lambda blocks:torch_scan_schedule(blocks,depth=2,decode_workers=2,
        shared_capacities=dict(cpu=profile['cpu_available_cores'],dram=profile['shared_dram_bytes_per_second'],input=profile['read_bytes_per_second']))
    assert schedule(slower)['seconds']>schedule(work['blocks'])['seconds']+1.
    assert runtime['prediction_complete'] is False
    assert 'per-device factor preparation' in runtime['scope']
    assert any('cache state' in term for term in runtime['unpriced_terms'])


def test_joint_host_serial_work_is_conserved_without_dense_finish_reuse():
    data,profile=joint_fixture()
    profile['host_serial_primitives']={'scan_api':3e-6,**{name:value/2 for name,value in profile['joint_host_primitives'].items()}}
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=scan_component):
        work=torch_scan_work(data,profile)
    component=next(iter(work['components'].values()))
    expected=3e-6+component['reduction_component']['host_dispatch_cpu_seconds']/2
    assert component['host_dispatch_serial_cpu_seconds']==pytest.approx(expected)
    assert work['blocks'][0]['host_submit_serial_cpu_seconds']==pytest.approx(expected)
    assert work['blocks'][0]['result_submit_seconds']==sum(
        count*profile['control_primitives'][name] for name,count in native_control_work(reduction='jagwas')['result_submit'].items())


@pytest.mark.parametrize('fault,message',[
    ('dense_finish','five-array'),('wrong_bytes','coverage'),('no_finish','five-array'),
    ('borrowed','owned results'),('fp64_stats','float32 statistics'),('rank','rank exceeds'),
    ('missing_joint','joint geometry'),('missing_price','joint host primitives'),('duplicate_joint','Duplicate joint')])
def test_joint_scan_refuses_unmodeled_paths(fault,message):
    data,profile=joint_fixture()
    if fault=='dense_finish':profile['result_finish_service'].pop('reduction')
    if fault=='wrong_bytes':profile['result_finish_service']['baseline_copy_bytes']=416
    if fault=='no_finish':profile.pop('result_finish_service')
    if fault=='borrowed':profile['result_ownership']='borrowed'
    if fault=='fp64_stats':profile['compute_dtype']='float64'
    if fault=='rank':data['traits_analyzed']=data['samples']
    if fault=='missing_joint':profile['joint_kernel_geometry']=[]
    if fault=='missing_price':profile.pop('joint_host_primitives')
    if fault=='duplicate_joint':profile['joint_kernel_geometry']*=2
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=scan_component),pytest.raises(ValueError,match=message):
        torch_scan_work(data,profile)


def test_dense_full_runtime_cannot_silently_write_joint_results_as_dense():
    data,profile=joint_fixture()
    with pytest.raises(ValueError,match='indexed writer'):torch_runtime(data,profile)
    profile.pop('reduction')
    with pytest.raises(ValueError,match='reduction mismatch'):torch_scan_work(data,profile)


def test_eager_composition_overlaps_host_submission_with_previous_gpu_work():
    scan=scan_component()
    scan.update(host_dispatch_cpu_seconds=2.,kernel_service_seconds=10.)
    scan['host_calls'][0].update(cpu_seconds=2.,submit_finish=2.)
    scan['operations'][0].update(host_submit_finish=2.,kernel_service_seconds=10.)
    projection=dict(reduction='jagwas',host_dispatch_cpu_seconds=3.,kernel_service_seconds=7.,
        host_calls=[dict(call_id=0,cpu_seconds=3.,submit_finish=3.,primitive='joint')],
        operations=[dict(host_submit_finish=3.,kernel_service_seconds=7.)],kernel_count=1,
        modeled_hbm_bytes=31,logical_bytes=29,source_sha256={},unpriced_terms=[])
    value=append_jagwas_projection(scan,projection,1.)
    assert value['estimated_span_seconds']==19.
    assert value['operations'][-1]['host_submit_finish']==5.
    assert value['host_dispatch_cpu_seconds']==5.
    assert [call['call_id'] for call in value['host_calls']]==[0,1]
    assert scan['operations'][0].get('gpu_finish') is None
