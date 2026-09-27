import copy
import pytest
from torchgwas.host_serial_work import attach_host_serial_work
from torchgwas.execution_graph import torch_scan_schedule, torch_multigpu_schedule


def component():
    return dict(host_dispatch_cpu_seconds=6., host_calls=[
        dict(primitive='view',cpu_seconds=1.,submit_finish=2.),
        dict(primitive='matmul',cpu_seconds=4.,submit_finish=10.),
        dict(primitive='view',cpu_seconds=1.,submit_finish=12.)], operations=[
        dict(host_submit_finish=10.,kernel_service_seconds=0.),
        dict(host_submit_finish=10.,kernel_service_seconds=0.)])


def test_views_and_multi_kernel_apis_conserve_serial_cpu_once():
    original=component()
    result=attach_host_serial_work(original,{'view':.5,'matmul':1.},.5)
    assert result['host_dispatch_serial_cpu_seconds']==2.
    assert [op['host_serial_cpu_finish'] for op in result['operations']]==[1.5,1.5]
    assert 'host_serial_cpu_finish' not in original['operations'][0]


@pytest.mark.parametrize('prices',[{'view':1.}, {'view':True,'matmul':1.},
    {'view':-.1,'matmul':1.},{'view':2.,'matmul':1.},{'view':float('nan'),'matmul':1.}])
def test_bad_primitive_service_rejected(prices):
    with pytest.raises(ValueError):attach_host_serial_work(component(),prices,.5)


def test_inconsistent_source_ledgers_rejected():
    for mutation in ['offset','total','kernel']:
        value=component()
        if mutation=='offset':value['host_calls'][0]['submit_finish']=3.
        if mutation=='total':value['host_dispatch_cpu_seconds']=10.
        if mutation=='kernel':value['operations'][0]['host_submit_finish']=3.
        with pytest.raises(ValueError):attach_host_serial_work(value,{'view':.5,'matmul':1.},.5)


def block(serial=2.):
    return dict(decode_seconds=0.,h2d_seconds=0.,d2h_seconds=0.,finish_seconds=0.,consumer_seconds=0.,
        host_resources={'cpu':1.},host_submit_seconds=10.,host_submit_serial_cpu_seconds=serial,
        operations=[dict(host_submit_finish=10.,host_serial_cpu_finish=serial,kernel_service_seconds=0.)])


def test_per_api_service_overrides_unknown_host_fraction_in_shared_schedule():
    def run(serial):
        shards=[dict(device=str(d),blocks=[block(serial)],depth=2,decode_workers=1) for d in [1,2]]
        return torch_multigpu_schedule(shards,{'cpu':2.},host_serial_fraction=1.)
    assert run(2.)['seconds']==pytest.approx(10.)
    assert run(10.)['seconds']==pytest.approx(20.)


def test_tail_keeps_unknown_result_control_separate_and_conserves_work():
    b=block(1.5)
    b.update(host_submit_seconds=6.,result_submit_seconds=4.,host_resources={'cpu':1.,'host_serial':1.})
    b['operations']=[dict(host_submit_finish=4.,host_serial_cpu_finish=1.,kernel_service_seconds=0.)]
    g=torch_scan_schedule([b],depth=2,decode_workers=1,shared_capacities={'cpu':2.},return_graph=True)
    work=sum(g.nodes[name][0]*demands.get('host_serial',0.) for name,demands in g.demands.items())
    assert work==pytest.approx(1.5+4.)


@pytest.mark.parametrize('total,offset',[(11.,1.), (2.,3.), (2.,None), (2.,-1.), (2.,True)])
def test_inconsistent_graph_serial_endpoints_rejected(total,offset):
    b=block(total);b['operations'][0]['host_serial_cpu_finish']=offset
    with pytest.raises(ValueError):torch_scan_schedule([b],depth=2,decode_workers=1,shared_capacities={'cpu':2.})
