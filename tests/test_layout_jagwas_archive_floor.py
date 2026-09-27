"""Compact indexed-writer intervals enclose the public per-part service."""
from copy import deepcopy

import pytest

from test_layout_frontier import productive_snapshot
from test_layout_jagwas_selection_floor import prices as selector_prices
from test_source_layout_floor import part, source
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.layout_jagwas_archive_floor import native_layout_jagwas_archive_floor
from torchgwas.layout_jagwas_selection_floor import native_layout_jagwas_selection_floor
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.reduced_output_work import (jagwas_archive_service,
                                           jagwas_indexed_part_work,
                                           jagwas_writer_work)
from torchgwas.source_layout_floor import native_layout_source_floor


def archive_price():
    return dict(field_schema=[['variant_index', '<i8'], ['chi2', '<f8']],
                call_cpu_seconds=.001, byte_cpu_seconds=.00001)


def profile():
    return dict(cpu_fraction=.5, shared_dram_bytes_per_second=1000.,
                fsync_seconds=.003,
                writeback_service=dict(pagecache_seconds_per_byte=.00001,
                                       storage_seconds_per_byte=.0001))


def single(tmp_path, retained, *, start=4, stop=9, chunk=2, fsync=True):
    floor = source(tmp_path)(start, stop, chunk)
    layout = native_layout_source_floor([
        part('shard', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    output = native_layout_output_floor(layout,
        retained_ranges={'shard': retained}, jagwas_writer_fsync=fsync,
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400.},
        output_bytes_per_second=200.)
    return layout, output


@pytest.mark.parametrize('retained,counts', [
    ([0, 0], [0, 0, 0]),
    ([5, 5], [2, 2, 1]),
    ([1, 4], [1, 0, 1]),
    ([1, 4], [0, 2, 0]),
])
def test_archive_bounds_existing_part_service(tmp_path, retained, counts):
    layout, output = single(tmp_path, retained)
    settings = profile()
    report = native_layout_jagwas_archive_floor(layout, output,
        archive_price(), {'cuda:0': settings})
    actual_cpu = actual_dram = actual_file = actual_serial = actual_parts = 0.
    for markers, keep in zip([2, 2, 1], counts):
        work = jagwas_writer_work(markers, keep)
        steps = jagwas_archive_service(work, archive_price(), settings,
                                       host_serial_fraction=0.)
        actual_cpu += sum(row['seconds'] * row.get('resources', {}).get('cpu', 0.)
                          for row in steps)
        actual_dram += sum(row['seconds'] * row.get('resources', {}).get('dram', 0.)
                           for row in steps)
        actual_file += work['part']['file_bytes']
        actual_serial += sum(row['seconds'] for row in steps)
        actual_parts += bool(keep)
    for key, value in [('total_cpu_seconds', actual_cpu),
                       ('total_logical_dram_bytes', actual_dram),
                       ('total_file_bytes', actual_file),
                       ('total_nonempty_parts', actual_parts),
                       ('single_consumer_floor_seconds', actual_serial)]:
        low, high = report[key]
        assert low <= value + 1e-12, key
        assert value <= high + 1e-12, key
    if retained in ([0, 0], [5, 5]):
        assert report['total_file_bytes'] == [actual_file, actual_file]


def test_archive_and_selector_share_consumer_cpu_dram_and_output(tmp_path):
    make = source(tmp_path)
    left, right = make(4, 8, 2), make(8, 12, 2)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), left),
        part('right', 'cuda:1', (0, 5), right)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    output = native_layout_output_floor(layout,
        retained_ranges={'left': [0, 4], 'right': [1, 4]},
        jagwas_writer_fsync=True,
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    selector = native_layout_jagwas_selection_floor(layout, output,
        selector_prices(), {'cuda:0': .5, 'cuda:1': .5})
    archive = native_layout_jagwas_archive_floor(layout, output,
        archive_price(), {'cuda:0': profile(), 'cuda:1': profile()})
    frontier = unissued_frontier(productive_snapshot(reduction='jagwas'),
        source_identity=layout['input_identity'], reduction='jagwas',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.},
        peak_fp64_flops_per_second={'cuda:0': 1000., 'cuda:1': 1000.})
    envelope = native_layout_partial_envelope(frontier, layout, compute,
        output, jagwas_selection=selector, jagwas_archive=archive)
    stages = envelope['stage_floor_seconds']
    assert stages['jagwas_archive_service'] == archive['archive_service_floor_seconds']
    assert stages['jagwas_single_consumer'] == [
        selector['single_consumer_cpu_floor_seconds'][i] +
        archive['single_consumer_floor_seconds'][i] for i in (0, 1)]
    assert stages['combined_shared_cpu'] == [
        (layout['resource_work']['cpu_seconds'][i] +
         selector['total_cpu_seconds'][i] + archive['total_cpu_seconds'][i]) /
        layout['shared_capacities']['cpu'] for i in (0, 1)]
    assert stages['combined_shared_dram'] == [
        (layout['resource_work']['dram_bytes'] + compute['total_h2d_bytes'] +
         output['total_d2h_payload_bytes'][i] +
         selector['total_logical_dram_bytes'][i] +
         archive['total_logical_dram_bytes'][i]) /
        layout['shared_capacities']['dram'] for i in (0, 1)]
    assert all(envelope['partial_floor_seconds'][i] >=
               stages['jagwas_single_consumer'][i] for i in (0, 1))
    changed = deepcopy(output)
    changed['partitions'][0]['retained_rows'] = [1, 4]
    with pytest.raises(ValueError, match='Matching durable JAGWAS archive'):
        native_layout_partial_envelope(frontier, layout, compute, changed,
                                       jagwas_archive=archive)


def test_archive_requires_durable_boundary_and_schema(tmp_path):
    layout, output = single(tmp_path, [0, 5], fsync=False)
    with pytest.raises(ValueError, match='Matching durable'):
        native_layout_jagwas_archive_floor(layout, output, archive_price(),
                                            {'cuda:0': profile()})
    layout, output = single(tmp_path, [0, 5])
    wrong = archive_price()
    wrong['field_schema'].append(['df', '<f4'])
    with pytest.raises(ValueError, match='two-array'):
        native_layout_jagwas_archive_floor(layout, output, wrong,
                                            {'cuda:0': profile()})


def test_two_array_header_overhead_is_monotone_at_shape_digit_boundaries():
    sizes = [1, 9, 10, 99, 100, 999, 1000, 9999, 10000, 65536]
    overhead = [jagwas_indexed_part_work(rows)['file_bytes'] - 16 * rows
                for rows in sizes]
    assert overhead == sorted(overhead)


def test_sparse_occupancy_caps_largest_possible_part(tmp_path):
    layout, output = single(tmp_path, [0, 1], chunk=64)
    report = native_layout_jagwas_archive_floor(layout, output,
        archive_price(), {'cuda:0': profile()})
    assert report['total_file_bytes'][1] == jagwas_indexed_part_work(1)['file_bytes']
    assert report['total_nonempty_parts'] == [0., 1.]
