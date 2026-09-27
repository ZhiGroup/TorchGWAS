"""Compact selector launches match source block geometry and join count waits."""
from copy import deepcopy

import pytest

from test_layout_significant_device_count_floor import fixture, price as count_price
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.layout_significant_device_count_floor import (
    native_layout_significant_device_count_floor)
from torchgwas.layout_significant_device_launch_floor import (
    _chunk_launch_counts, native_layout_significant_device_launch_floor)
from torchgwas.reduced_output_work import significant_output_work


def profile():
    return dict(torch_version='2.5.1+cu124', cuda_runtime='12.4',
                compute_capability=[8, 0], kernel_launch_seconds=2.,
                gpu_fraction=1.)


@pytest.mark.parametrize('markers,traits,limit', [
    (13, 7, 17), (512, 600000, 1 << 20), (4100, 2, 1 << 20),
    (5000, 5000, 8192), (1, 5000, 4097)])
def test_compact_launch_classes_equal_expanded_source_blocks(markers, traits, limit):
    blocks, two_count, largest = _chunk_launch_counts(markers, traits, limit)
    source = significant_output_work(129, markers, traits, markers,
        backend='device', max_selection_cells=limit, max_blocks=100000)
    assert blocks == len(source['blocks'])
    assert two_count == sum(row['cells'] > 4096 for row in source['blocks'])
    assert largest == max(row['cells'] for row in source['blocks'])


def test_empty_device_launches_join_serial_count_barriers(tmp_path):
    source, output = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    profiles = {'cuda:0': profile(), 'cuda:1': profile()}
    launches = native_layout_significant_device_launch_floor(source, output, profiles)
    profiles['cuda:0']['compute_capability'][0] = 9
    assert launches['launch_profiles_by_device']['cuda:0'][
        'compute_capability'] == [8, 0]
    count = native_layout_significant_device_count_floor(
        source, output, {'cuda:0': count_price(), 'cuda:1': count_price()})
    assert [row['selection_blocks'] for row in launches['partitions']] == [8, 12]
    assert all(row['coordinate_scatter_launches'] == [0, 0]
               for row in launches['partitions'])
    assert launches['selector_launch_floor_seconds'] == [72., 72.]
    from test_layout_frontier import productive_snapshot
    from torchgwas.layout_frontier import unissued_frontier
    frontier = unissued_frontier(productive_snapshot(reduction='significant'),
        source_identity=source['input_identity'], reduction='significant',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(source, covariate_rank=3,
        shared_h2d_bytes_per_second=1e9,
        per_device_h2d_bytes_per_second={'cuda:0': 1e9, 'cuda:1': 1e9},
        peak_fp32_flops_per_second={'cuda:0': 1e9, 'cuda:1': 1e9})
    bound = native_layout_partial_envelope(frontier, source, compute, output,
        significant_device_count=count, significant_device_launch=launches)
    stages = bound['stage_floor_seconds']
    assert stages['significant_device_count_barrier'] == [36., 36.]
    assert stages['significant_device_launch_service'] == [72., 72.]
    assert stages['significant_device_launch_count_serial'] == [108., 108.]
    assert bound['partial_floor_seconds'][0] >= 108.
    other = tmp_path / 'same_gpu'
    other.mkdir()
    serial_source, serial_output = fixture(other, ('cuda:0', 'cuda:0'))
    serial_launch = native_layout_significant_device_launch_floor(
        serial_source, serial_output, {'cuda:0': profile()})
    serial_count = native_layout_significant_device_count_floor(
        serial_source, serial_output, {'cuda:0': count_price()})
    assert serial_launch['selector_launch_floor_seconds'] == [120., 120.]
    assert serial_count['count_barrier_floor_seconds'] == 60.
    changed = deepcopy(output)
    changed['partitions'][0]['device_selection_blocks'] += 1
    with pytest.raises(ValueError, match='Matching device selector launch'):
        native_layout_partial_envelope(frontier, source, compute, changed,
            significant_device_launch=launches)


def test_retained_intervals_bound_scatter_launches(tmp_path):
    source, _ = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    output = native_layout_output_floor(source,
        significant_backend='device', device_selection_max_cells=2,
        retained_ranges={'low': [1, 4], 'high': [0, 5]},
        shared_d2h_bytes_per_second=1e9,
        per_device_d2h_bytes_per_second={'cuda:0': 1e9, 'cuda:1': 1e9},
        output_bytes_per_second=1e9)
    launches = native_layout_significant_device_launch_floor(source, output,
        {'cuda:0': profile(), 'cuda:1': profile()})
    low, high = launches['partitions']
    assert low['coordinate_scatter_launches'] == [1, 4]
    assert high['coordinate_scatter_launches'] == [0, 5]
    assert launches['selector_launch_floor_seconds'] == [72., 82.]


@pytest.mark.parametrize('damage', ['version', 'build', 'capability', 'fraction', 'missing'])
def test_unbound_or_unsupported_launch_price_rejected(tmp_path, damage):
    source, output = fixture(tmp_path, ('cuda:0', 'cuda:1'))
    prices = {'cuda:0': profile(), 'cuda:1': profile()}
    if damage == 'version': prices['cuda:0']['torch_version'] = '2.6.0'
    elif damage == 'build': prices['cuda:0']['torch_version'] = '2.5.1+cu121'
    elif damage == 'capability': prices['cuda:0']['compute_capability'] = [7, 5]
    elif damage == 'fraction': prices['cuda:0']['gpu_fraction'] = 0.
    else: prices.pop('cuda:1')
    with pytest.raises(ValueError):
        native_layout_significant_device_launch_floor(source, output, prices)