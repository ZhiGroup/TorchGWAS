"""Hand-calculated dependency/occupancy controls, not GWAS timing fits."""
import copy
import pytest
from torchgwas.selection_geometry import DEVICE_SELECTION_MAX_CELLS
from torchgwas.significant_schedule import significant_trait_schedule


def step(seconds):
    return [dict(seconds=seconds)]


def tile(device='cuda:1', backend='device', chunks=2, kernel=2., selection=3., writer=5., retained=1):
    block = dict(decode_seconds=1., h2d_seconds=1., operations=[dict(host_submit_finish=0., kernel_service_seconds=kernel)],
        host_submit_seconds=0., d2h_seconds=1., finish_seconds=0., consumer_seconds=0.)
    output = dict(cells=8, retained=retained, selection=step(selection), writer=step(writer) if retained else [])
    return dict(device=device, backend=backend, blocks=[copy.deepcopy(block) for _ in range(chunks)],
        outputs=[[copy.deepcopy(output)] for _ in range(chunks)], depth=2, decode_workers=2, cleanup=[])


def run(tiles, **kwargs):
    multiple = len({t['device'] for t in tiles}) > 1
    opts = dict(queue_depth=1 if multiple else 0, shared_capacities={},
        queue_service=dict(put=step(0.), get=step(0.)) if multiple else None, finalize=step(2.))
    opts.update(kwargs)
    return significant_trait_schedule(tiles, **opts)


def test_synchronous_device_results_hold_next_submission_until_writer_returns():
    result = run([tile()])
    # Chunk 0: decode 1, H2D 1, compute 2, status D2H 1,
    # select 3, durable writer 5 = 13. Chunk 1 was decoded ahead.
    assert result['seconds'] == 27.  # 13 + (1+2+1+3+5) + metadata 2
    assert result['start']['tile:0:host_start:1'] == result['end']['tile:0:consume:0'] == 13.
    assert result['end']['tile:0:decode:1'] < 13.
    assert result['selected_pairs'] == 2 and result['parts'] == 2


def test_host_dense_ring_keeps_the_existing_overlap():
    result = run([tile(backend='host')])
    assert result['start']['tile:0:host_start:1'] < result['end']['tile:0:consume:0']
    assert result['seconds'] < run([tile()])['seconds']


def test_shared_writer_uses_arrival_order_not_device_or_tile_order():
    tiles = [tile(chunks=1, kernel=10.), tile(device='cuda:2', chunks=1, kernel=1.)]
    result = run(tiles)
    start, end = result['start'], result['end']
    fast = 'tile:1:significant:0:0:write'
    slow = 'tile:0:significant:0:0:write'
    assert end[fast + ':done'] <= start[slow + ':0']
    assert start[fast + ':0'] < start[slow + ':0']
    assert result['seconds'] == end['significant:finalize:done']


def test_one_slot_queue_limits_admission_and_writer_never_overlaps():
    result = run([tile(chunks=4, writer=20.), tile(device='cuda:2', chunks=4, writer=20.)])
    start, end = result['start'], result['end']
    writes = sorted((start[name], end[name]) for name in start if name.endswith(':write:0'))
    assert len(writes) == 8
    assert all(before[1] <= after[0] for before, after in zip(writes, writes[1:]))
    events = []
    for name in start:
        if name.endswith(':acquire'):
            events.append((start[name], 1))
        if name.endswith(':get_start'):
            events.append((start[name], -1))
    # Resolve simultaneous release/acquire together, so zero-service transfers
    # are not assigned an artificial order by the test.
    count = 0
    for when in sorted(set(time for time, _ in events)):
        count += sum(delta for time, delta in events if time == when)
        assert 0 <= count <= 1
    assert count == 0


def test_next_tile_can_prepare_while_previous_parts_are_still_writing():
    tiles = [tile(chunks=1, writer=30.), tile(device='cuda:2', chunks=1, kernel=100.),
             tile(chunks=1, writer=30.)]
    result = run(tiles)
    assert result['start']['tile:2:submit_decode:0'] < result['end']['tile:0:significant:0:0:write:done']
    assert result['start']['tile:2:submit_decode:0'] >= result['end']['tile:0:significant:producer_complete']


def test_empty_selection_does_not_create_part_but_still_selects_and_queues():
    result = run([tile(retained=0), tile(device='cuda:2', retained=0)])
    assert result['parts'] == 0 and result['selected_pairs'] == 0 and result['selections'] == 4
    assert len([name for name in result['start'] if name.endswith(':get_start')]) == 6


def test_refuse_unknown_occupancy_missing_stages_or_unbounded_block():
    for field, value, message in [('retained', None, 'retained'), ('selection', [], 'Selection'),
                                   ('writer', [], 'writer'), ('cells', DEVICE_SELECTION_MAX_CELLS + 1, 'bound')]:
        candidate = tile()
        candidate['outputs'][0][0][field] = value
        with pytest.raises(ValueError, match=message):
            run([candidate])
    with pytest.raises(ValueError, match='queue put/get'):
        run([tile(), tile(device='cuda:2')], queue_service=None)


def test_existing_native_and_reference_solver_agree_on_reduced_queue_graph():
    graph = run([tile(chunks=3), tile(device='cuda:2', chunks=3)], return_graph=True)
    assert graph.solve() == graph._solve_shared_python()

def test_finite_schedule_rejects_expansion_beyond_explicit_bounds():
    with pytest.raises(ValueError, match='max_source_chunks'):
        run([tile(chunks=2)], max_source_chunks=1)
    with pytest.raises(ValueError, match='max_selection_blocks'):
        run([tile(chunks=2)], max_selection_blocks=1)