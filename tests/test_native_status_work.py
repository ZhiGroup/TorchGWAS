"""Independent status service and its native scan/selector ownership boundary."""
import copy
from unittest.mock import patch
import pytest
from test_mechanistic_shapes import fixture,component
from test_device_significance_service import selector_fixture
from torchgwas.native_control_work import native_control_work
from torchgwas.native_status_work import native_status_service
from torchgwas.mechanistic_torch import torch_scan_work,torch_scan_runtime
from torchgwas.indexed_schedule import significant_trait_schedule


def status_prices():
    return dict(dtype='uint8',copy=dict(before_cpu_seconds=2e-6,after_cpu_seconds=3e-6),
        transfer=dict(latency_seconds=4e-6,bytes_per_second=2e8),host_serial_fraction=.5,wait_cpu_fraction=.25,
        qc={name:dict(call_cpu_seconds=1e-6,row_cpu_seconds=2e-9) for name in ['malformed_empty','count_status']})


def device_fixture():
    data,profile=fixture(k=7)
    profile.update(reduction='device_significant',device_status_service=status_prices(),event_wait_cpu_fraction=0.,
        control_primitives={key:1e-6 for phase in native_control_work(reduction='device_significant').values() for key in phase})
    return data,profile


def test_status_cpu_and_traffic_conservation():
    model=native_status_service(32,status_prices(),cpu_fraction=.5,dram_bytes_per_second=1e9)
    assert model['status_bytes']==32
    assert model['d2h_seconds']==pytest.approx(4e-6+32/2e8)
    assert model['cpu_seconds']==pytest.approx(8e-6+3*32*2e-9)
    assert model['qc_work']['logical_dram_bytes']==9*32
    assert sum(p['seconds'] for p in model['finish_operations'])==model['finish_seconds']
    assert model['result_wait_resources']=={'cpu':.125}


@pytest.mark.parametrize('change,match',[
    ({'dtype':'float32'},'uint8'),({'copy':{'cpu_seconds':1.}},'partitioned'),
    ({'transfer':{}},'pageable'),({'wait_cpu_fraction':-1.},'bounded'),
    ({'qc':{}},'malformed')])
def test_missing_or_incompatible_independent_prices_are_refused(change,match):
    prices=status_prices();prices.update(change)
    with pytest.raises(ValueError,match=match):
        native_status_service(32,prices,cpu_fraction=1.,dram_bytes_per_second=1e9)


def test_native_device_path_does_not_price_dense_transfer_or_result_copy():
    data,profile=device_fixture();original=copy.deepcopy(profile)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(data,profile)
        # These dense-only parameters are not consulted by this branch.
        profile['process_units'].update(finish_fixed_calls=999.,numpy_copy_bytes=999.)
        same=torch_scan_work(data,profile)
    assert work['blocks']==same['blocks']
    assert work['cpu_work_seconds']==same['cpu_work_seconds']
    assert sum(b['d2h_bytes'] for b in work['blocks'])==10
    for rows,block in zip([8,2],work['blocks']):
        assert block['resolve_seconds']==0.
        assert block['result_submit_seconds']==2e-6
        assert block['d2h_seconds']==pytest.approx(4e-6+rows/2e8)
        assert work['owned_result_work'][rows]['work']['allocation_calls']==0
    profile['process_units']=original['process_units']
    assert profile==original


@pytest.mark.parametrize('change,match',[
    ({'result_ownership':'borrowed'},'owned'),({'compute_dtype':'float64'},'float32'),
    ({'result_finish_service':{}},'finish'),({'control_primitives':None},'control'),
    ({'device_status_service':None},'status'),({'owned_result_copy_scenario':{}},'dense result copier')])
def test_native_device_mode_refuses_dense_contracts(change,match):
    data,profile=device_fixture();profile.update(change)
    with pytest.raises(ValueError,match=match):torch_scan_work(data,profile)


def test_scan_only_api_cannot_report_partial_device_selection_as_complete_scan():
    data,profile=device_fixture()
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        with pytest.raises(ValueError,match='synchronous indexed selection'):
            torch_scan_runtime(data,profile)


def test_native_status_composes_with_selector_and_releases_input_independently():
    data,profile=device_fixture()
    # Large status dispatch exposes its ordering independently of all tiny
    # source/native services. Input release must not await this blocking call.
    profile['device_status_service']['copy']['before_cpu_seconds']=5.
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        work=torch_scan_work(data,profile)
    selections=[selector_fixture([0]*rows)[2] for rows in [8,2]]
    for model in selections:model['graph'].capacities['cpu']=4.
    tile=dict(device='cuda:0',backend='device',blocks=work['blocks'],depth=2,decode_workers=2,
        outputs=[[dict(cells=7,retained=0,selection=[],writer=[]) for _ in range(rows)] for rows in [8,2]],
        selection_graphs=selections,cleanup=[])
    graph=significant_trait_schedule([tile],queue_depth=0,queue_service=None,finalize=[],
        shared_capacities={'cpu':4.,'host_serial':1.,'dram':1e9,'input':1e8,'d2h':4.},return_graph=True)
    solved=graph.solve();assert solved==graph._solve_shared_python()
    assert solved['end']['tile:0:release:0']<solved['end']['tile:0:status_submit:0']
    assert solved['start']['tile:0:d2h:0']>=solved['end']['tile:0:status_submit:0']
    assert solved['start']['tile:0:significant:0:device:begin']>=solved['end']['tile:0:finish:0']
    assert solved['start']['tile:0:host_start:1']>=solved['end']['tile:0:consume:0']
    assert solved['start']['tile:0:kernel:1:0']>=solved['end']['tile:0:significant:0:device:'+selections[0]['gpu_tail']]
