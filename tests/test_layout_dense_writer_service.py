"""Compact writer service conserves the finite graph's priced CPU work."""
from copy import deepcopy

import pytest

from test_layout_partial_envelope import case
from test_source_layout_floor import part, source
from torchgwas.binary_output_work import binary_output_work
from torchgwas.binary_schedule import BinaryWriterSchedule
from torchgwas.execution_graph import ExecutionGraph
from torchgwas.layout_dense_writer_service import native_layout_dense_writer_service_floor
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.source_layout_floor import native_layout_source_floor


OPTIONS = dict(block_bytes=32, queue_depth=2, borrow_chunks=False,
               fsync=True, writeback_bytes=64, sync_file_range=True,
               store_variant_df=True)


def profile():
    return dict(cpu_fraction=.5,
                writer_copy_service=dict(cpu_seconds_per_byte=.0001,
                                         cpu_seconds_per_call=.001),
                process_units=dict(bytearray_zero_bytes=.00001),
                executor_cpu_seconds=.002,
                writeback_service=dict(pagecache_seconds_per_byte=.0001,
                                       storage_seconds_per_byte=.0002,
                                       submit_seconds=.003, wait_seconds=.004,
                                       fadvise_seconds=.005,
                                       fadvise_eviction_seconds_per_byte=.00001),
                fsync_seconds=.01)


def priced_output(layout):
    devices = {row['device'] for row in layout['partitions']}
    return native_layout_output_floor(layout,
        shared_d2h_bytes_per_second=1000.,
        per_device_d2h_bytes_per_second={device: 1000. for device in devices},
        output_bytes_per_second=1000., dense_writer_options=OPTIONS)


def test_compact_cpu_work_matches_expanded_writer_graph(tmp_path):
    floor = source(tmp_path)(0, 12, 3)
    layout = native_layout_source_floor([
        part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction=None, partition_axis='variant')
    output = priced_output(layout)
    settings = profile()
    report = native_layout_dense_writer_service_floor(
        layout, output, {'cuda:0': settings})
    expanded = binary_output_work(12, 5, 3, store_beta=True, **OPTIONS)
    q = settings['cpu_fraction']
    writeback = {key: value if key == 'storage_seconds_per_byte' else value / q
                 for key, value in settings['writeback_service'].items()}
    writer = BinaryWriterSchedule(expanded,
        copy_seconds_per_byte=settings['writer_copy_service']['cpu_seconds_per_byte'] / q,
        copy_seconds_per_call=settings['writer_copy_service']['cpu_seconds_per_call'] / q,
        zero_seconds_per_byte=settings['process_units']['bytearray_zero_bytes'] / q,
        write_seconds_per_byte=1.,
        fsync_seconds_per_array=settings['fsync_seconds'],
        append_seconds=settings['executor_cpu_seconds'] / q,
        handoff_seconds=settings['executor_cpu_seconds'] / q,
        cpu_fraction=q, write_capacity=5000.,
        writeback_bytes=OPTIONS['writeback_bytes'], writeback_service=writeback)
    graph = ExecutionGraph()
    graph.capacities = {'cpu': layout['shared_capacities']['cpu'], 'output': 5000.}
    last = graph.add('start')
    for index in range(4):
        last = writer.append(graph, index, [last])
    writer.close(graph, [last])
    cpu = sum(seconds * graph.demands.get(name, {}).get('cpu', 0.)
              for name, (seconds, _) in graph.nodes.items())
    assert report['total_writer_cpu_seconds'] == pytest.approx(cpu)
    assert report['writer_service_floor_seconds'] <= graph.solve()['seconds']
    assert report['partitions'][0]['write_calls_minimum'] == writer.write_calls


def test_writer_service_enters_matching_unissued_envelope(tmp_path):
    frontier, layout, compute, _ = case(tmp_path)
    output = priced_output(layout)
    settings = {device: profile() for device in ('cuda:0', 'cuda:1')}
    writer = native_layout_dense_writer_service_floor(layout, output, settings)
    report = native_layout_partial_envelope(frontier, layout, compute, output,
                                            writer_service=writer)
    assert report['stage_floor_seconds']['dense_writer_service'] == [
        writer['writer_service_floor_seconds']] * 2
    assert report['partial_floor_seconds'][0] >= writer['writer_service_floor_seconds']
    damaged = deepcopy(output)
    damaged['dense_writer_work'][0]['work']['staging_copy_calls'] += 1
    with pytest.raises(ValueError):
        native_layout_partial_envelope(frontier, layout, compute, damaged,
                                       writer_service=writer)
    changed_block = deepcopy(output)
    changed_block['dense_writer_work'][0]['work']['block_bytes'] += 1
    with pytest.raises(ValueError, match='Matching dense writer'):
        native_layout_partial_envelope(frontier, layout, compute, changed_block,
                                       writer_service=writer)


def test_writer_service_requires_explicit_prices_and_matching_output(tmp_path):
    floor = source(tmp_path)(0, 12, 3)
    layout = native_layout_source_floor([
        part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction=None, partition_axis='variant')
    output = priced_output(layout)
    settings = profile()
    settings.pop('writer_copy_service')
    with pytest.raises(ValueError):
        native_layout_dense_writer_service_floor(layout, output, {'cuda:0': settings})
    settings['process_units']['numpy_copy_bytes'] = .0001
    fallback = native_layout_dense_writer_service_floor(
        layout, output, {'cuda:0': settings})
    assert fallback['partitions'][0]['copy_price_source'] == (
        'legacy_numpy_copy_bytes_without_call_price')
    assert fallback['unpriced_terms']
    settings = profile()
    settings['writeback_service']['storage_seconds_per_byte'] = 0.
    with pytest.raises(ValueError):
        native_layout_dense_writer_service_floor(layout, output, {'cuda:0': settings})
    bare = native_layout_output_floor(layout,
        shared_d2h_bytes_per_second=1000.,
        per_device_d2h_bytes_per_second={'cuda:0': 1000.},
        output_bytes_per_second=1000.)
    with pytest.raises(ValueError):
        native_layout_dense_writer_service_floor(layout, bare, {'cuda:0': profile()})
