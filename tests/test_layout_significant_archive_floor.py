"""Compact host significant archive encloses per-chunk finite services."""
from copy import deepcopy

import pytest

from test_layout_frontier import productive_snapshot
from test_source_layout_floor import part, source
from torchgwas.selection_geometry import DEVICE_SELECTION_MAX_CELLS
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.layout_significant_archive_floor import native_layout_significant_archive_floor
from torchgwas.significant_host_model import _writer_service
from torchgwas.significant_host_work import indexed_part_work
from torchgwas.source_layout_floor import native_layout_source_floor


def prices():
    return {str(beta): dict(call_cpu_seconds=.001,
                            byte_cpu_seconds=.00001)
            for beta in (False, True)}


def profile():
    return dict(cpu_fraction=.5, shared_dram_bytes_per_second=1000.,
                fsync_seconds=.003,
                writeback_service=dict(pagecache_seconds_per_byte=.00001,
                                       storage_seconds_per_byte=.0001),
                process_units=dict(numpy_copy_bytes=.00001))


def single(tmp_path, retained, *, beta=True, backend='host', fsync=True,
           selection_max_cells=DEVICE_SELECTION_MAX_CELLS):
    floor = source(tmp_path)(4, 9, 2)
    layout = native_layout_source_floor([
        part('tile', 'cuda:0', (0, 3), floor)],
        total_traits=3, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(layout, store_beta=beta,
        significant_backend=backend,
        device_selection_max_cells=selection_max_cells,
        retained_ranges={'tile': retained}, significant_writer_fsync=fsync,
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400.},
        output_bytes_per_second=200.)
    return layout, output


@pytest.mark.parametrize('beta,retained,counts', [
    (True, [0, 0], [0, 0, 0]),
    (False, [0, 0], [0, 0, 0]),
    (True, [15, 15], [6, 6, 3]),
    (False, [15, 15], [6, 6, 3]),
    (True, [1, 8], [2, 0, 1]),
    (False, [1, 8], [0, 5, 1]),
])
def test_compact_archive_bounds_per_chunk_writer(tmp_path, beta, retained, counts):
    layout, output = single(tmp_path, retained, beta=beta)
    settings = profile()
    report = native_layout_significant_archive_floor(layout, output,
        prices(), {'cuda:0': settings})
    actual_cpu = actual_dram = actual_file = actual_serial = actual_parts = 0.
    for count in counts:
        part_work = indexed_part_work(count, store_beta=beta)
        steps = _writer_service(part_work, prices(), settings, serial=0.)
        actual_cpu += sum(row['seconds'] * row.get('resources', {}).get('cpu', 0.)
                          for row in steps)
        actual_dram += sum(row['seconds'] * row.get('resources', {}).get('dram', 0.)
                           for row in steps)
        actual_file += part_work['file_bytes']
        actual_serial += sum(row['seconds'] for row in steps)
        actual_parts += bool(count)
    for key, value in [('total_cpu_seconds', actual_cpu),
                       ('total_logical_dram_bytes', actual_dram),
                       ('total_file_bytes', actual_file),
                       ('total_nonempty_parts', actual_parts),
                       ('single_consumer_floor_seconds', actual_serial)]:
        low, high = report[key]
        assert low <= value + 1e-12, key
        assert value <= high + 1e-12, key
    assert report['store_beta'] is beta


@pytest.mark.parametrize('backend', ['host', 'device'])
def test_archive_enters_matching_significant_partial_envelope(tmp_path, backend):
    floor = source(tmp_path)(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(layout,
        significant_backend=backend, significant_writer_fsync=True,
        device_selection_max_cells=3 if backend == 'device' else DEVICE_SELECTION_MAX_CELLS,
        retained_ranges={'low': [0, 16], 'high': [1, 24]},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    archive = native_layout_significant_archive_floor(layout, output,
        prices(), {'cuda:0': profile(), 'cuda:1': profile()})
    frontier = unissued_frontier(productive_snapshot(reduction='significant'),
        source_identity=layout['input_identity'], reduction='significant',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.})
    envelope = native_layout_partial_envelope(frontier, layout, compute,
        output, significant_archive=archive)
    stages = envelope['stage_floor_seconds']
    assert stages['significant_archive_service'] == archive['archive_service_floor_seconds']
    assert stages['significant_single_consumer'] == (
        archive['single_consumer_floor_seconds'])
    assert stages['combined_shared_cpu'] == [
        (layout['resource_work']['cpu_seconds'][i] +
         archive['total_cpu_seconds'][i]) /
        layout['shared_capacities']['cpu'] for i in (0, 1)]
    assert stages['combined_shared_dram'] == [
        (layout['resource_work']['dram_bytes'] + compute['total_h2d_bytes'] +
         output['total_d2h_payload_bytes'][i] +
         archive['total_logical_dram_bytes'][i]) /
        layout['shared_capacities']['dram'] for i in (0, 1)]
    damaged = deepcopy(output)
    damaged['partitions'][0]['retained_rows'] = [1, 16]
    with pytest.raises(ValueError, match='Matching significant archive'):
        native_layout_partial_envelope(frontier, layout, compute, damaged,
                                       significant_archive=archive)


@pytest.mark.parametrize('beta,retained,counts', [
    (True, [0, 0], [0, 0, 0, 0, 0]),
    (False, [15, 15], [3, 3, 3, 3, 3]),
    (True, [1, 8], [0, 3, 0, 1, 0]),
])
def test_device_archive_bounds_per_selection_block(tmp_path, beta, retained, counts):
    layout, output = single(tmp_path, retained, beta=beta, backend='device',
                            selection_max_cells=4)
    report = native_layout_significant_archive_floor(
        layout, output, prices(), {'cuda:0': profile()})
    assert report['partitions'][0]['possible_parts'] == 5
    assert report['partitions'][0]['maximum_part_pairs'] == 3
    assert output['partitions'][0]['d2h_payload_bytes'] == [
        5 + 4 * 5 + 28 * value for value in retained]
    actual = {key: 0. for key in ('total_cpu_seconds', 'total_logical_dram_bytes',
                                  'total_file_bytes', 'total_nonempty_parts',
                                  'single_consumer_floor_seconds')}
    for count in counts:
        work = indexed_part_work(count, store_beta=beta)
        steps = _writer_service(work, prices(), profile(), serial=0.)
        actual['total_cpu_seconds'] += sum(
            row['seconds'] * row.get('resources', {}).get('cpu', 0.) for row in steps)
        actual['total_logical_dram_bytes'] += sum(
            row['seconds'] * row.get('resources', {}).get('dram', 0.) for row in steps)
        actual['total_file_bytes'] += work['file_bytes']
        actual['total_nonempty_parts'] += bool(count)
        actual['single_consumer_floor_seconds'] += sum(row['seconds'] for row in steps)
    for key, value in actual.items():
        low, high = report[key]
        assert low <= value + 1e-12, key
        assert value <= high + 1e-12, key
    changed = deepcopy(output)
    changed['partitions'][0]['device_selection_blocks'] += 1
    with pytest.raises(ValueError, match='Device selection blocks'):
        native_layout_significant_archive_floor(layout, changed, prices(),
                                                 {'cuda:0': profile()})


def test_nondurable_writer_needs_separate_service(tmp_path):
    layout, output = single(tmp_path, [0, 15], fsync=False)
    with pytest.raises(ValueError, match='Matching durable significant'):
        native_layout_significant_archive_floor(layout, output, prices(),
                                                 {'cuda:0': profile()})


def test_t_only_significant_output_still_transports_scan_beta(tmp_path):
    layout, output = single(tmp_path, [1, 8], beta=False)
    assert output['partitions'][0]['d2h_payload_bytes'] == [
        8 * 5 * 3 + 5 * 5] * 2
    assert output['partitions'][0]['output_array_payload_bytes'] == [24, 192]
    archive = native_layout_significant_archive_floor(layout, output,
        prices(), {'cuda:0': profile()})
    assert archive['store_beta'] is False
