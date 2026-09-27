import copy
import pytest
from torchgwas.result_service import result_finish_prices, acknowledgement_prices


def observation(cpu, detached=0., intervals=0):
    return dict(cpu_seconds=cpu, detached_cpu_seconds=detached,
                detached_intervals=intervals, balance_errors=0)


def fixture():
    worker = dict(worker=0, context_verified=True,
        controls=dict(native=observation(1., .95, 1), python=observation(1.)), meter_rows=[], rows=[])
    for repeat in range(3):
        for record in [False, True]:
            worker['meter_rows'].append(dict(observation(40. if record else 20., 10. if record else 0., 20000 if record else 0),
                repeat=repeat, recording=record, intervals=20000))
            for mode in ['owned', 'borrowed']:
                cpu = (.1 if record else .09)+(.01 if mode=='owned' else 0.)
                worker['rows'].append(dict(repeat=repeat, ownership=mode, recording=record, worker=0,
                    samples=[observation(cpu, .012 if record else 0., 8 if record else 0) for _ in range(3)],
                    empty_samples=[observation(.001) for _ in range(3)],
                    dummy_event_samples=[observation(.002) for _ in range(3)]))
    probe = dict(context_verified=True, numpy_version='np', torch_version='torch', python_version='py',
        affinity=[1, 2], numpy_madvise_hugepage=False, source_sha256={'src/torchgwas/native_scan.py':'abc'},
        timing_contract='Exact source finish(0,0,32), K=1, status clear, no reduction or p-values.',
        results=[dict(worker_count=1, workers=[worker])])
    kw = dict(workers=1, numpy_version='np', torch_version='torch', python_version='py', cpu_affinity=[1,2], source_sha256='abc')
    return probe, kw


def test_finish_replaces_tiny_control_and_keeps_expensive_calls():
    probe, kw = fixture()
    value = result_finish_prices(probe, **kw)['prices']['borrowed']
    assert value['cpu_seconds'] == pytest.approx(.09)
    assert value['serial_cpu_seconds'] == pytest.approx(.082)
    assert value['baseline_copy_bytes'] == 0
    for row in probe['results'][0]['workers'][0]['rows']:
        if row['recording'] and row['ownership']=='borrowed':
            row['samples'][-1]['cpu_seconds'] += .3
    assert result_finish_prices(probe, **kw)['prices']['borrowed']['cpu_seconds'] == pytest.approx(.19)


@pytest.mark.parametrize('fault', ['source','runtime','duplicate','missing','meter','balance','dummy','disabled','negative'])
def test_bad_finish_evidence_is_rejected(fault):
    probe, kw = fixture(); worker = probe['results'][0]['workers'][0]
    if fault=='source': kw['source_sha256']='different'
    if fault=='runtime': kw['numpy_version']='different'
    if fault=='duplicate': worker['rows'].append(copy.deepcopy(worker['rows'][0]))
    if fault=='missing': worker['rows'].pop()
    if fault=='meter': worker['meter_rows'][0]['detached_intervals']=1
    if fault=='balance': worker['rows'][0]['samples'][0]['balance_errors']=1
    if fault=='dummy': worker['rows'][2]['dummy_event_samples'][0].update(detached_intervals=1,detached_cpu_seconds=.001)
    if fault=='disabled': worker['rows'][0]['samples'][0]['detached_intervals']=1
    if fault=='negative':
        for row in worker['rows']:
            row['dummy_event_samples']=[observation(1.) for _ in range(3)]
    with pytest.raises(ValueError): result_finish_prices(probe, **kw)


def ack_fixture():
    ready=dict(python_version='py',affinity=[1,2],results=[dict(pairs=1,workers=[dict(worker=0,rows=[
        dict(repeat=i,calls=300,loop_cpu_seconds=.001,publish_cpu_seconds=.002,receive_cpu_seconds=.003)
        for i in range(3)])])])
    blocked=dict(python_version='py',affinity=[1,2],acknowledgement=[dict(pairs=1,repeat=i,workers=[dict(worker=0,rows=[
        dict(blocked_verified=True,create_cpu_seconds=.0005,publisher_cpu_seconds=.002,
             receiver_cpu_seconds=.003,publication_to_return_seconds=.010) for _ in range(3)])]) for i in range(3)])
    return blocked,ready,dict(pairs=1,python_version='py',cpu_affinity=[1,2])


def test_acknowledgement_uses_ready_costs_and_only_extra_blocked_latency():
    blocked,ready,kw=ack_fixture();value=acknowledgement_prices(blocked,ready,**kw)
    assert value['create_cpu_seconds']==pytest.approx(.0005)
    assert value['publish_cpu_seconds']==pytest.approx(.001)
    assert value['receive_cpu_seconds']==pytest.approx(.002)
    assert value['wakeup_seconds']==pytest.approx(.007)
    blocked['acknowledgement'][0]['workers'][0]['rows'][0]['blocked_verified']=False
    with pytest.raises(ValueError,match='verified'):acknowledgement_prices(blocked,ready,**kw)


def test_acknowledgement_rejects_duplicate_and_mismatched_context():
    blocked,ready,kw=ack_fixture();ready['affinity']=[3]
    with pytest.raises(ValueError,match='context'):acknowledgement_prices(blocked,ready,**kw)
    blocked,ready,kw=ack_fixture();blocked['acknowledgement'].append(copy.deepcopy(blocked['acknowledgement'][0]))
    with pytest.raises(ValueError,match='coverage'):acknowledgement_prices(blocked,ready,**kw)
