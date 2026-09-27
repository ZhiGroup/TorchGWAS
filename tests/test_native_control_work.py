import pytest
from torchgwas.native_control_work import native_control_work,native_control_service
from torchgwas.execution_graph import torch_scan_schedule


def test_control_reuse_has_only_previous_slot_waits():
    cold=native_control_work(False);hot=native_control_work(True)
    assert cold['result_submit']==hot['result_submit']
    assert cold['result_submit']['copy_d2h']==4
    assert cold['result_submit']['tensor_record_stream']==4
    assert cold['result_submit']['executor_submit']==2
    assert cold['transfer_submit']['future_result']==0
    assert hot['transfer_submit']['future_result']==1
    assert hot['transfer_submit']['stream_wait_event']==2
    prices={key:1. for counts in hot.values() for key in counts}
    a=native_control_service(prices,False)['cpu_seconds']
    b=native_control_service(prices,True)['cpu_seconds']
    assert sum(b.values())-sum(a.values())==2.
    with pytest.raises(ValueError,match='Missing native control'):
        native_control_service({})


def test_release_callback_gates_pinned_slot_reuse():
    block=dict(decode_seconds=1.,h2d_seconds=1.,operations=[],host_submit_seconds=0.,
               d2h_seconds=0.,finish_seconds=0.,consumer_seconds=0.,release_seconds=5.)
    result=torch_scan_schedule([block]*3,depth=2,decode_workers=2)
    # One release thread serializes its callbacks; slot 0 cannot refill before
    # its first H2D completes plus the callback's independently priced service.
    assert result['start']['submit_decode:2']==result['end']['release:0']
    assert result['start']['release:1']>=result['end']['release:0']
    assert result['end']['release:0']-result['end']['h2d:0']==5.


def test_publish_and_resolve_cpu_are_on_the_consumer_path():
    block=dict(decode_seconds=1.,h2d_seconds=1.,operations=[],host_submit_seconds=0.,
               d2h_seconds=0.,finish_seconds=0.,consumer_seconds=0.,
               publish_seconds=2.,resolve_seconds=3.)
    result=torch_scan_schedule([block],depth=2,decode_workers=1)
    assert result['seconds']==7.
    assert result['end']['publish:0']==3.
    assert result['end']['consume:0']==7.


def test_release_cannot_start_before_the_main_thread_submits_it():
    block=dict(decode_seconds=1.,h2d_seconds=1.,
               operations=[dict(host_submit_finish=10.,kernel_service_seconds=1.)],
               host_submit_seconds=10.,d2h_seconds=0.,finish_seconds=0.,consumer_seconds=0.)
    result=torch_scan_schedule([block]*3,depth=2,decode_workers=2)
    assert result['end']['h2d:0']<result['end']['host_done:0']
    assert result['start']['release:0']==result['end']['host_done:0']
    assert result['start']['submit_decode:2']==result['end']['release:0']


def test_early_release_refills_input_during_statistics_submission():
    block=dict(decode_seconds=1.,h2d_seconds=1.,
               operations=[dict(host_submit_finish=10.,kernel_service_seconds=1.)],
               host_submit_seconds=10.,d2h_seconds=0.,finish_seconds=0.,
               consumer_seconds=0.,release_after_transfer=True)
    result=torch_scan_schedule([block]*3,depth=2,decode_workers=2)
    assert result['start']['release:0']==result['end']['h2d:0']
    assert result['start']['submit_decode:2']==result['end']['release:0']
    assert result['end']['decode:2']<result['end']['host_done:0']
    # An unfinished DMA remains the boundary even with an early callback.
    block.update(h2d_seconds=20.,host_submit_seconds=0.,operations=[],
                 event_wait_resources={'cpu':1.})
    result=torch_scan_schedule([block]*3,depth=2,decode_workers=2,
                               shared_capacities={'cpu':4.})
    assert result['end']['begin_release_wait:0']<result['end']['h2d:0']
    assert result['start']['release:0']==result['end']['h2d:0']


def test_reader_initialization_precedes_payload_read_and_is_paid_once():
    block=dict(decode_seconds=2.,decode_read_seconds=3.,h2d_seconds=0.,
               operations=[],host_submit_seconds=0.,d2h_seconds=0.,
               finish_seconds=0.,consumer_seconds=0.)
    result=torch_scan_schedule([dict(block,reader_init_seconds=5.),block],
                               depth=2,decode_workers=1)
    assert result['start']['decode_read:0']==5.
    assert result['start']['decode:0']==8.
    assert result['end']['decode:1']==15.
    assert 'reader_init:1' not in result['start']


def test_pending_window_already_gates_previous_worker_job_for_unequal_chunks():
    blocks=[dict(decode_seconds=value,decode_read_seconds=.2,h2d_seconds=.1,operations=[],
                 host_submit_seconds=0.,d2h_seconds=0.,finish_seconds=0.,consumer_seconds=0.)
            for value in [1.,.01,2.,.1,3.,.001,1.5]]
    graph=torch_scan_schedule(blocks,depth=4,decode_workers=2,return_graph=True)
    original=graph.solve()
    for i in range(2,len(blocks)):
        name=f'decode_read:{i}'
        seconds,deps=graph.nodes[name]
        assert f'decode:{i-2}' in deps
        graph.nodes[name]=(seconds,tuple(dep for dep in deps if dep!=f'decode:{i-2}'))
    reduced=graph.solve()
    assert reduced['start']==pytest.approx(original['start'])
    assert reduced['end']==pytest.approx(original['end'])