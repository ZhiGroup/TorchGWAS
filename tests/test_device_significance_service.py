"""Every current device selector API has an independent typed probe."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from torchgwas.device_significance_service import device_selection_host_primitives
from torchgwas.device_significance_work import device_significant_tensor_work


@pytest.mark.parametrize('shape,counts', [((13,7,17),[0]*7),((13,7,17),[1]*7),((1,37,11),[0,1,2,4])])
def test_all_source_calls_resolve_to_available_fixed_primitives(shape, counts):
    build_bank = pytest.importorskip('direct_device_significance_primitives').build_bank
    b,k,limit = shape
    work = device_significant_tensor_work(40,b,k,counts,max_cells=limit)
    resolved = device_selection_host_primitives(work)
    assert len(resolved['calls']) == len(work['host_calls'])
    bank = build_bank('cpu')
    assert set(resolved['counts']) <= set(bank)
    copies = sum(resolved['counts'].get(name,0) for name in ['copy_cpu_numpy_fp32','copy_cpu_numpy_int64','copy_cpu_numpy_int32'])
    assert copies == work['selected_copy_calls']
    assert sum(resolved['counts'].get(name,0) for name in ['nonzero_empty','nonzero_nonempty']) == len(counts)
    for name in resolved['counts']:
        value = bank[name]()
        assert value is not None


def test_unknown_source_call_is_not_assigned_a_generic_price():
    work = device_significant_tensor_work(40,1,7,[1])
    call = next(c for c in work['host_calls'] if c['name']=='isfinite')
    work['steps'][call['step_indices'][0]]['op'] = 'unknown'
    with pytest.raises(ValueError,match='Unpriced'):
        device_selection_host_primitives(work)


def test_fixed_cpu_bank_matches_source_aten_dispatch_for_noncopy_calls():
    import torch
    from torch.utils._python_dispatch import TorchDispatchMode
    build_bank = pytest.importorskip('direct_device_significance_primitives').build_bank

    class Record(TorchDispatchMode):
        def __init__(self): self.ops = []
        def __torch_dispatch__(self, function, types, args=(), kwargs=None):
            self.ops.append(str(function))
            return function(*args, **(kwargs or {}))

    work = device_significant_tensor_work(40,13,7,[1]*7,max_cells=17)
    resolved = device_selection_host_primitives(work)
    bank = build_bank('cpu')
    for call in resolved['calls']:
        if call['primitive'].startswith('copy_'): continue
        observer = Record()
        with observer: bank[call['primitive']]()
        expected = [work['steps'][i]['op'] for i in call['step_indices']]
        assert observer.ops == expected, (call['primitive'], observer.ops, expected)


def selector_fixture(counts, *, cells=7, count_seconds=2., select_seconds=10., scatter_seconds=2.):
    from torchgwas.device_significance_service import device_selection_graph
    work = device_significant_tensor_work(40,len(counts),cells,counts,max_cells=cells)
    primitives = device_selection_host_primitives(work)['counts']
    prices = {name:({'before_cpu_seconds':0.,'after_cpu_seconds':0.}
        if name.startswith(('copy_cpu_numpy_','nonzero_')) else {'cpu_seconds':0.}) for name in primitives}
    services = {str(i):dict(op=s['op'],seconds=0.) for i,s in enumerate(work['steps'])
        if not(s['alias_only'] or s['allocation_only'] or s['host_copy'] or s['op']=='aten.nonzero.default')}
    phases = [dict(count=dict(seconds=count_seconds),flagged_select=dict(seconds=select_seconds),
        **({'coordinate_scatter':dict(seconds=scatter_seconds)} if count else {})) for count in counts]
    transfers = {name:dict(latency_seconds=0.,bytes_per_second=4.,resources=['d2h']) for name in ['count','payload']}
    kwargs = dict(host_prices=prices,operation_services=services,nonzero_services=phases,
        transfer_prices=transfers,yield_cpu_seconds=0.,capacities={'cpu':1.,'d2h':4.})
    return work,kwargs,device_selection_graph(work,**kwargs)


def test_empty_yield_precedes_flagged_select_drain_and_next_count_stays_on_stream():
    _,_,model = selector_fixture([0,0])
    solved = model['graph'].solve()
    assert solved['end'][model['yield_nodes'][0]] == 3.
    assert solved['end'][model['block_gpu_tails'][0]] == 13.
    assert solved['end'][model['yield_nodes'][1]] == 16.
    assert solved['seconds'] == 26.
    assert model['nonzero_count_d2h_bytes'] == 8
    assert model['selected_payload_d2h_bytes'] == 0


def test_nonempty_yield_waits_for_coordinates_and_one_packed_payload_copy():
    _,_,model = selector_fixture([1])
    solved = model['graph'].solve()
    # count 2 + count copy 1 + select 10 + scatter 2 + 20 packed bytes at 4 B/s.
    assert solved['end'][model['yield_nodes'][0]] == 20.
    assert model['selected_payload_d2h_bytes'] == 20
    assert model['nonzero_count_d2h_bytes'] == 4


def test_writer_pause_can_overlap_empty_selection_gpu_tail():
    _,_,model = selector_fixture([0,0])
    graph = model['graph']
    writer = graph.add('external_writer',1.,after=[model['yield_nodes'][0]])
    resume = model['resume_nodes'][0]
    seconds,deps = graph.nodes[resume]; graph.nodes[resume] = (seconds,(*deps,writer))
    solved = graph.solve()
    assert solved['end'][writer] == 4.
    assert solved['seconds'] == 26.


def test_whole_blocking_api_timing_cannot_silently_double_count_gpu_work():
    from torchgwas.device_significance_service import device_selection_graph
    work,kwargs,_ = selector_fixture([0])
    kwargs['host_prices']['nonzero_empty'] = {'cpu_seconds':0.0001}
    with pytest.raises(ValueError,match='partitioned around the barrier'):
        device_selection_graph(work,**kwargs)


def test_unknown_gpu_operation_service_is_refused():
    from torchgwas.device_significance_service import device_selection_graph
    work,kwargs,_ = selector_fixture([0])
    next(iter(kwargs['operation_services'].values()))['op'] = 'unknown'
    with pytest.raises(ValueError,match='source-operation mismatch'):
        device_selection_graph(work,**kwargs)


def test_indexed_pipeline_preserves_pending_selector_stream_work_between_chunks():
    from test_significant_schedule import tile, run
    candidate = tile(chunks=2,retained=0)
    candidate['selection_graphs'] = [selector_fixture([0])[2] for _ in range(2)]
    for group in candidate['outputs']:
        group[0].update(cells=7,selection=[])
    graph = run([candidate],return_graph=True)
    solved = graph.solve()
    assert solved == graph._solve_shared_python()
    # First chunk yields at 8 while its GPU tail continues until 18. The
    # next H2D can start at 8, but its statistics kernel shares that stream.
    assert solved['end']['tile:0:consume:0'] == 8.
    assert solved['start']['tile:0:h2d:1'] == 8.
    assert solved['start']['tile:0:kernel:1:0'] == 18.
    assert solved['seconds'] == 36.


def test_single_writer_resume_port_delays_next_selection_and_preserves_owned_output():
    from test_significant_schedule import tile, run
    candidate = tile(chunks=1,retained=1,writer=7.)
    candidate['selection_graphs'] = [selector_fixture([1,1])[2]]
    candidate['outputs'] = [[dict(cells=7,retained=1,selection=[],writer=[dict(seconds=7.)]) for _ in range(2)]]
    solved = run([candidate])
    prefix = 'tile:0:significant:0:device:'
    first_yield,second_yield = candidate['selection_graphs'][0]['yield_nodes']
    # Each nonempty block copies one 20-byte packed payload (5 s at 4 B/s).
    assert solved['end'][prefix+first_yield] == 25.
    assert solved['end'][prefix+'block:0:resume'] == 32.
    assert solved['end'][prefix+second_yield] == 52.
    assert solved['seconds'] == 61.


def test_device_graph_rejects_wrong_occupancy_or_duplicate_selection_service():
    from test_significant_schedule import tile, run
    candidate = tile(chunks=1,retained=0)
    candidate['selection_graphs'] = [selector_fixture([1])[2]]
    with pytest.raises(ValueError,match='occupancy'):
        run([candidate])
    candidate['selection_graphs'] = [selector_fixture([0])[2]]
    with pytest.raises(ValueError,match='twice'):
        run([candidate])


def test_explicit_wait_cpu_conserves_resource_work_without_becoming_dispatch():
    from torchgwas.device_significance_service import device_selection_graph
    from torchgwas.resource_balance import resource_balance
    work,kwargs,_ = selector_fixture([1])
    kwargs['wait_cpu_fraction'] = 1.
    model = device_selection_graph(work,**kwargs)
    graph = model['graph']; solved = graph.solve()
    assert solved == graph._solve_shared_python()
    assert solved['seconds'] == 20.
    assert model['host_cpu_seconds'] == 0.
    balance = resource_balance(graph,solved)
    assert balance['resources']['cpu']['work'] == 20.


def test_two_device_source_graphs_share_one_writer_and_resume_after_queue_admission():
    from test_significant_schedule import tile,run
    candidates=[]
    for device in ['cuda:1','cuda:2']:
        candidate=tile(device=device,chunks=2,retained=1,writer=20.)
        candidate['selection_graphs']=[selector_fixture([1,1])[2] for _ in range(2)]
        candidate['outputs']=[[dict(cells=7,retained=1,selection=[],writer=[dict(seconds=20.)]) for _ in range(2)] for _ in range(2)]
        candidates.append(candidate)
    graph=run(candidates,return_graph=True)
    solved=graph.solve()
    assert solved==graph._solve_shared_python()
    writes=sorted((solved['start'][n],solved['end'][n]) for n in graph.nodes if n.endswith(':write:0'))
    assert len(writes)==8
    assert all(a[1]<=b[0] for a,b in zip(writes,writes[1:]))
    for device in range(2):
        for chunk in range(2):
            for part in range(2):
                prefix=f'tile:{device}:significant:{chunk}:'
                assert solved['end'][prefix+f'device:block:{part}:resume']==solved['end'][prefix+f'{part}:put:done']
