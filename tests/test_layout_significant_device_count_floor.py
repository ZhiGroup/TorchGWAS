"""Device count barriers follow selection blocks and per-GPU serial ownership."""
from copy import deepcopy

import pytest

from test_layout_frontier import productive_snapshot
from test_source_layout_floor import part, source
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.layout_significant_device_count_floor import (
    native_layout_significant_device_count_floor)
from torchgwas.source_layout_floor import native_layout_source_floor


def fixture(tmp_path, devices):
    floor = source(tmp_path)(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', devices[0], (0, 2), floor),
        part('high', devices[1], (2, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(
        layout, significant_backend='device', device_selection_max_cells=2,
        retained_ranges={'low': [0, 0], 'high': [0, 0]},
        shared_d2h_bytes_per_second=1e9,
        per_device_d2h_bytes_per_second={device: 1e9 for device in set(devices)},
        output_bytes_per_second=1e9)
    return layout, output


def price():
    # The source graph prices every count transfer as latency + 4/rate.
    return dict(latency_seconds=2., bytes_per_second=4., resources=['d2h'])


def test_empty_device_blocks_still_serialize_count_barriers_per_gpu(tmp_path):
    layout, output = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    report = native_layout_significant_device_count_floor(
        layout, output, {'cuda:0': price(), 'cuda:1': price()})
    assert [row['selection_blocks'] for row in report['partitions']] == [8, 12]
    assert [row['count_d2h_bytes'] for row in report['partitions']] == [32, 48]
    assert report['per_device_serial_seconds'] == {'cuda:0': 24., 'cuda:1': 36.}
    assert report['count_barrier_floor_seconds'] == 36.
    assert output['total_d2h_payload_bytes'] == [8 * 2 + 4 * (8 + 12)] * 2
    serial_layout, serial_output = fixture(tmp_path, ('cuda:0', 'cuda:0'))
    serial = native_layout_significant_device_count_floor(
        serial_layout, serial_output, {'cuda:0': price()})
    assert serial['count_barrier_floor_seconds'] == 60.


def test_compact_count_latency_is_a_floor_of_the_source_selector_graph(tmp_path):
    from test_device_significance_service import selector_fixture
    from torchgwas.device_significance_service import device_selection_graph

    floor = source(tmp_path)(0, 2, 2)
    layout = native_layout_source_floor([
        part('tile', 'cuda:0', (0, 3), floor)],
        total_traits=3, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(
        layout, significant_backend='device', device_selection_max_cells=3,
        retained_ranges={'tile': [0, 0]},
        shared_d2h_bytes_per_second=1e9,
        per_device_d2h_bytes_per_second={'cuda:0': 1e9},
        output_bytes_per_second=1e9)
    compact = native_layout_significant_device_count_floor(
        layout, output, {'cuda:0': price()})
    work, kwargs, _ = selector_fixture([0, 0], cells=3,
                                        count_seconds=0., select_seconds=0.)
    kwargs['transfer_prices']['count'] = price()
    exact = device_selection_graph(work, **kwargs)['graph'].solve()
    assert compact['count_barrier_floor_seconds'] == 6.
    assert compact['count_barrier_floor_seconds'] <= exact['seconds']


def test_count_floor_binds_to_unissued_partial_envelope(tmp_path):
    layout, output = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    count = native_layout_significant_device_count_floor(
        layout, output, {'cuda:0': price(), 'cuda:1': price()})
    frontier = unissued_frontier(productive_snapshot(reduction='significant'),
        source_identity=layout['input_identity'], reduction='significant',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1e9,
        per_device_h2d_bytes_per_second={'cuda:0': 1e9, 'cuda:1': 1e9},
        peak_fp32_flops_per_second={'cuda:0': 1e9, 'cuda:1': 1e9})
    report = native_layout_partial_envelope(
        frontier, layout, compute, output, significant_device_count=count)
    assert report['stage_floor_seconds']['significant_device_count_barrier'] == [36., 36.]
    assert all(value >= 36. for value in report['partial_floor_seconds'])
    changed = deepcopy(output)
    changed['partitions'][0]['device_selection_blocks'] += 1
    with pytest.raises(ValueError, match='Matching device count-transfer floor'):
        native_layout_partial_envelope(
            frontier, layout, compute, changed, significant_device_count=count)


@pytest.mark.parametrize('damage', ['price', 'missing', 'mode', 'geometry'])
def test_count_floor_rejects_unbound_price_or_geometry(tmp_path, damage):
    layout, output = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    prices = {'cuda:0': price(), 'cuda:1': price()}
    if damage == 'price':
        prices['cuda:0']['latency_seconds'] = -1.
    elif damage == 'missing':
        prices.pop('cuda:1')
    elif damage == 'mode':
        output['significant_backend'] = 'host'
    else:
        output['partitions'][1]['device_selection_blocks'] -= 1
    with pytest.raises(ValueError):
        native_layout_significant_device_count_floor(layout, output, prices)


def test_count_floor_keeps_an_immutable_price_snapshot(tmp_path):
    layout, output = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    prices = {'cuda:0': price(), 'cuda:1': price()}
    report = native_layout_significant_device_count_floor(layout, output, prices)
    prices['cuda:0']['latency_seconds'] = 100.
    prices['cuda:1']['resources'].append('changed')
    assert report['count_transfer_prices_by_device']['cuda:0']['latency_seconds'] == 2.
    assert report['count_transfer_prices_by_device']['cuda:1']['resources'] == ['d2h']
