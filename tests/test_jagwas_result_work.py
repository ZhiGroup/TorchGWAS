"""Reconcile the result ledger against real buffers and the exact finish body."""
import ast
import inspect
import numpy as np
import pytest
import torch
from torchgwas.native_scan import dosage_cuda_iterator
from torchgwas.reduce import JagwasReduction
from torchgwas.owned_result_work import owned_result_work,owned_result_copy_service
from torchgwas.native_control_work import native_control_work,native_control_service
from torchgwas.result_service import result_finish_prices
from torchgwas.allocator_service import owned_allocator_service
from test_result_service import fixture


@pytest.mark.parametrize('markers',[1,32,501,4096])
def test_source_host_buffers_determine_copy_extents_and_control_calls(markers):
    buffers=JagwasReduction().host_buffers(markers,1,pin_memory=False)
    work=owned_result_work(markers,2048,reduction='jagwas')
    assert list(work['array_bytes'].values())==[t.numel()*t.element_size() for t in buffers]
    assert list(work['array_elements'].values())==[t.numel() for t in buffers]
    assert work['allocation_calls']==len(buffers)
    counts=native_control_work(reduction='jagwas')
    assert {counts['result_submit'][key] for key in ['tensor_slice','copy_d2h','tensor_record_stream']}=={len(buffers)}
    assert counts['finish']['tensor_slice']==counts['finish']['tensor_numpy']==len(buffers)
    prices={key:1. for phase in counts.values() for key in phase}
    dense=native_control_service(prices)['cpu_seconds']
    joint=native_control_service(prices,reduction='jagwas')['cpu_seconds']
    assert joint['result_submit']-dense['result_submit']==3.
    assert joint['finish']-dense['finish']==2.
    with pytest.raises(ValueError,match='matching finish'):owned_result_copy_service(work,{})


@pytest.mark.parametrize('borrowed',[False,True])
def test_exact_finish_returns_index_releases_df_and_invalidates_rows(borrowed):
    source=ast.parse(inspect.getsource(dosage_cuda_iterator))
    body=next(node for node in source.body[0].body if isinstance(node,ast.FunctionDef) and node.name=='finish')
    buffers=JagwasReduction().host_buffers(4,1,pin_memory=False)
    for value in buffers:value.fill_(1)
    buffers[3].copy_(torch.tensor([0,1,2,0],dtype=torch.uint8))
    buffers[4].copy_(torch.tensor([30.,29.,28.,27.]))
    class Ready:
        def synchronize(self):pass
    namespace=dict(np=np,result_done=[Ready()],result_buffers=[buffers],borrow_results=borrowed,
        reduction=JagwasReduction(),compute_p_values=False,profiling=False,return_df=False)
    exec(compile(ast.Module(body=[body],type_ignores=[]),'exact_native_finish','exec'),namespace)
    emitted,missing,invariant,times=namespace['finish'](0,45,49)
    assert emitted[:2]==(45,49) and emitted[4] is None
    assert emitted[2].flags.owndata is (not borrowed)
    assert len(emitted)==6 and emitted[3].shape==emitted[5].shape==(4,1)
    assert np.isnan(emitted[2][1:3]).all() and np.isnan(emitted[3][1:3]).all()
    assert (missing,invariant,times)==(1,1,None)
    assert sum(a.nbytes for a in [emitted[2],emitted[3],emitted[5]])==12*4
    assert owned_result_work(4,999,reduction='jagwas')['copy_bytes']==17*4


def test_finish_evidence_cannot_cross_reduction_or_returned_df_contexts():
    probe,kwargs=fixture()
    with pytest.raises(ValueError):result_finish_prices(probe,reduction='jagwas',**kwargs)
    probe.update(reduction='jagwas',return_df=False,
        result_layout=owned_result_work(32,1,reduction='jagwas')['array_bytes'],
        timing_contract='Exact source finish(0,0,32), K=1, status clear, JAGWAS reduction, no p-values.')
    with pytest.raises(ValueError):result_finish_prices(probe,**kwargs)
    prices=result_finish_prices(probe,reduction='jagwas',**kwargs)['prices']
    assert prices['owned']['baseline_copy_bytes']==544
    assert prices['owned']['result_arrays']==5 and prices['owned']['return_df'] is False
    assert prices['borrowed']['baseline_copy_bytes']==0
    probe['return_df']=True
    with pytest.raises(ValueError):result_finish_prices(probe,reduction='jagwas',**kwargs)


def test_joint_index_is_retained_but_status_and_df_are_released_by_worker():
    work=owned_result_work(1024,512,reduction='jagwas')
    geometry=dict(array_bytes={str(size):dict(route='mmap',mapped_bytes=4096) for size in work['array_bytes'].values()})
    prices=dict(page_bytes=4096,allocate_cpu_seconds_per_mmap=1.,release_cpu_seconds_per_mapped_page=2.)
    service=owned_allocator_service(work,prices,geometry)
    assert service['worker_release_cpu_seconds']==4.
    assert service['consumer_release_cpu_seconds']==6.
    assert service['arrays']['status']['release_owner']=='finish_worker'
    assert service['arrays']['trait_index']['release_owner']=='consumer'
    assert service['arrays']['df']['release_owner']=='finish_worker'
