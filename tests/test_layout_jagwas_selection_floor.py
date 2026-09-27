"""The compact JAGWAS selector bounds the existing per-chunk service."""
from copy import deepcopy

import pytest

from test_layout_frontier import productive_snapshot
from test_source_layout_floor import part, source
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.layout_jagwas_selection_floor import native_layout_jagwas_selection_floor
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.reduced_output_work import (jagwas_host_selection_service,
                                           jagwas_writer_work)
from torchgwas.source_layout_floor import native_layout_source_floor


def prices():
    names = ('fp32_to_fp64_view', 'finite_fp64', 'flatnonzero_empty',
             'flatnonzero_nonempty', 'index_add', 'fp64_gather')
    return {name: dict(call_cpu_seconds=(index + 1) * .001,
                       unit_cpu_seconds=(index + 1) * .00001)
            for index, name in enumerate(names)}


def case(tmp_path, ranges):
    make = source(tmp_path)
    left, right = make(4, 8, 2), make(8, 12, 2)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), left),
        part('right', 'cuda:1', (0, 5), right)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    output = native_layout_output_floor(layout, retained_ranges=ranges,
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    return layout, output


def expanded_service_cpu(markers_per_chunk, retained_per_chunk, bank, fraction):
    cpu = dram = 0.
    for markers, retained in zip(markers_per_chunk, retained_per_chunk):
        work = jagwas_writer_work(markers, retained)
        steps = jagwas_host_selection_service(work, bank,
            cpu_fraction=fraction, dram_bytes_per_second=1e8,
            host_serial_fraction=0.)
        cpu += sum(step['seconds'] * step['resources']['cpu'] for step in steps)
        dram += sum(step['seconds'] * step['resources']['dram'] for step in steps)
    return cpu, dram


@pytest.mark.parametrize('ranges,retained_chunks', [
    ({'left': [0, 0], 'right': [0, 0]}, ([0, 0], [0, 0])),
    ({'left': [4, 4], 'right': [4, 4]}, ([2, 2], [2, 2])),
    ({'left': [1, 3], 'right': [0, 2]}, ([1, 1], [0, 1])),
])
def test_compact_selector_bounds_expanded_service(tmp_path, ranges, retained_chunks):
    layout, output = case(tmp_path, ranges)
    bank = prices()
    fractions = {'cuda:0': .5, 'cuda:1': .25}
    report = native_layout_jagwas_selection_floor(layout, output, bank,
                                                    fractions)
    actual_cpu = actual_dram = actual_serial = 0.
    for partition, counts in zip(layout['partitions'], retained_chunks):
        cpu, dram = expanded_service_cpu([2, 2], counts, bank,
                                          fractions[partition['device']])
        actual_cpu += cpu
        actual_dram += dram
        actual_serial += cpu / fractions[partition['device']]
    assert report['total_cpu_seconds'][0] <= actual_cpu + 1e-12
    assert actual_cpu <= report['total_cpu_seconds'][1] + 1e-12
    assert report['total_logical_dram_bytes'][0] <= actual_dram
    assert actual_dram <= report['total_logical_dram_bytes'][1]
    assert report['single_consumer_cpu_floor_seconds'][0] <= actual_serial + 1e-12
    assert actual_serial <= report['single_consumer_cpu_floor_seconds'][1] + 1e-12
    if ranges['left'] == [0, 0] or ranges['left'] == [4, 4]:
        assert report['total_cpu_seconds'][0] == pytest.approx(actual_cpu)
        assert report['total_cpu_seconds'][1] == pytest.approx(actual_cpu)


def test_jagwas_selector_binds_unissued_partial_envelope(tmp_path):
    layout, output = case(tmp_path, {'left': [0, 4], 'right': [0, 4]})
    frontier = unissued_frontier(productive_snapshot(reduction='jagwas'),
        source_identity=layout['input_identity'], reduction='jagwas',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.},
        peak_fp64_flops_per_second={'cuda:0': 1000., 'cuda:1': 1000.})
    selector = native_layout_jagwas_selection_floor(layout, output, prices(),
        {'cuda:0': .5, 'cuda:1': .5})
    envelope = native_layout_partial_envelope(frontier, layout, compute, output,
                                               jagwas_selection=selector)
    assert envelope['stage_floor_seconds']['jagwas_selection_service'] == (
        selector['selector_service_floor_seconds'])
    assert envelope['partial_floor_seconds'][0] >= (
        selector['selector_service_floor_seconds'][0])
    assert envelope['stage_floor_seconds']['combined_shared_dram'] == [
        (layout['resource_work']['dram_bytes'] + compute['total_h2d_bytes'] +
         output['total_d2h_payload_bytes'][i] +
         selector['total_logical_dram_bytes'][i]) /
        layout['shared_capacities']['dram'] for i in (0, 1)]
    changed = deepcopy(output)
    changed['partitions'][0]['retained_rows'] = [1, 4]
    with pytest.raises(ValueError, match='Matching JAGWAS selector'):
        native_layout_partial_envelope(frontier, layout, compute, changed,
                                       jagwas_selection=selector)


def test_short_final_chunk_keeps_empty_and_nonempty_bounds(tmp_path):
    floor = source(tmp_path)(4, 9, 2)
    layout = native_layout_source_floor([
        part('tail', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    output = native_layout_output_floor(layout,
        retained_ranges={'tail': [1, 3]},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400.},
        output_bytes_per_second=200.)
    bank = prices()
    report = native_layout_jagwas_selection_floor(layout, output, bank,
                                                    {'cuda:0': .5})
    cpu, dram = expanded_service_cpu([2, 2, 1], [1, 0, 1], bank, .5)
    assert report['partitions'][0]['chunks'] == 3
    assert report['total_cpu_seconds'][0] <= cpu <= report['total_cpu_seconds'][1]
    assert report['total_logical_dram_bytes'][0] <= dram <= (
        report['total_logical_dram_bytes'][1])


def test_selector_requires_complete_independent_prices_and_full_panel(tmp_path):
    layout, output = case(tmp_path, {'left': [0, 4], 'right': [0, 4]})
    incomplete = prices()
    incomplete.pop('flatnonzero_empty')
    with pytest.raises(ValueError, match='Complete independent'):
        native_layout_jagwas_selection_floor(layout, output, incomplete,
                                              {'cuda:0': .5, 'cuda:1': .5})
    with pytest.raises(ValueError, match='One measured'):
        native_layout_jagwas_selection_floor(layout, output, prices(),
                                              {'cuda:0': .5})
    damaged = deepcopy(layout)
    damaged['partitions'][0]['trait_range'] = [0, 4]
    with pytest.raises(ValueError, match='geometry'):
        native_layout_jagwas_selection_floor(damaged, output, prices(),
                                              {'cuda:0': .5, 'cuda:1': .5})
