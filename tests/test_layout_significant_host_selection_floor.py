"""Compact host-selector intervals enclose the public per-chunk primitives."""
from copy import deepcopy

import pytest

from test_layout_frontier import productive_snapshot
from test_source_layout_floor import part, source
from test_layout_significant_archive_floor import (prices as archive_prices,
                                                   profile as archive_profile)
from test_significant_host_model import bank as price_bank
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_significant_archive_floor import native_layout_significant_archive_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.layout_significant_host_selection_floor import (
    native_layout_significant_host_selection_floor)
from torchgwas.significant_host_work import (host_significant_selection_work,
                                             host_selection_service)
from torchgwas.source_layout_floor import native_layout_source_floor


def profile():
    return dict(cpu_fraction=.5, shared_dram_bytes_per_second=1e8)


def single(tmp_path, retained, *, threshold_one=False, store_beta=True):
    floor = source(tmp_path)(4, 9, 2)
    layout = native_layout_source_floor([
        part('tile', 'cuda:0', (0, 3), floor)],
        total_traits=3, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(layout, store_beta=store_beta,
        significant_backend='host', significant_threshold_one=threshold_one,
        retained_ranges={'tile': retained},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400.},
        output_bytes_per_second=200.)
    return layout, output


@pytest.mark.parametrize('retained,counts', [
    ([0, 0], [0, 0, 0]),
    ([15, 15], [6, 6, 3]),
    ([1, 8], [1, 0, 2]),
    ([1, 8], [0, 5, 1]),
])
@pytest.mark.parametrize('threshold_one', [False, True])
def test_compact_selector_bounds_existing_host_steps(tmp_path, retained,
                                                      counts, threshold_one):
    layout, output = single(tmp_path, retained,
                            threshold_one=threshold_one)
    evidence = price_bank()
    report = native_layout_significant_host_selection_floor(layout, output,
        evidence, {'cuda:0': profile()})
    actual_cpu = actual_dram = 0.
    for markers, keep in zip([2, 2, 1], counts):
        work = host_significant_selection_work(markers, 3, keep,
            threshold_one=threshold_one, return_beta=True)
        steps = host_selection_service(work, evidence['prices'],
            cpu_fraction=.5, dram_bytes_per_second=1e8,
            host_serial_fraction=0.)
        actual_cpu += sum(row['seconds'] * row['resources']['cpu']
                          for row in steps)
        actual_dram += sum(row['seconds'] * row['resources']['dram']
                           for row in steps)
    assert report['total_cpu_seconds'][0] <= actual_cpu + 1e-12
    assert actual_cpu <= report['total_cpu_seconds'][1] + 1e-12
    assert report['total_logical_dram_bytes'][0] <= actual_dram + 1e-12
    assert actual_dram <= report['total_logical_dram_bytes'][1] + 1e-12
    if retained in ([0, 0], [15, 15]):
        assert report['total_cpu_seconds'][0] == pytest.approx(actual_cpu)
        assert report['total_cpu_seconds'][1] == pytest.approx(actual_cpu)


def test_t_only_writer_keeps_full_public_host_selection(tmp_path):
    layout, output = single(tmp_path, [1, 8], store_beta=False)
    report = native_layout_significant_host_selection_floor(layout, output,
        price_bank(), {'cuda:0': profile()})
    assert report['return_beta'] is True
    assert output['partitions'][0]['d2h_payload_bytes'] == [145, 145]


def test_host_selector_binds_partial_envelope_and_shared_resources(tmp_path):
    floor = source(tmp_path)(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(layout,
        significant_backend='host', significant_threshold_one=False,
        significant_writer_fsync=True,
        retained_ranges={'low': [0, 16], 'high': [1, 24]},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    selector = native_layout_significant_host_selection_floor(layout, output,
        price_bank(), {'cuda:0': profile(), 'cuda:1': profile()})
    archive = native_layout_significant_archive_floor(layout, output,
        archive_prices(), {'cuda:0': archive_profile(),
                           'cuda:1': archive_profile()})
    frontier = unissued_frontier(productive_snapshot(reduction='significant'),
        source_identity=layout['input_identity'], reduction='significant',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.})
    envelope = native_layout_partial_envelope(frontier, layout, compute,
        output, significant_host_selection=selector,
        significant_archive=archive)
    stages = envelope['stage_floor_seconds']
    assert stages['significant_host_selector_service'] == (
        selector['selector_service_floor_seconds'])
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
    changed = deepcopy(output)
    changed['significant_threshold_one'] = True
    with pytest.raises(ValueError, match='Matching host significant selector'):
        native_layout_partial_envelope(frontier, layout, compute, changed,
                                       significant_host_selection=selector)


def test_selector_requires_bound_threshold_and_host_backend(tmp_path):
    layout, output = single(tmp_path, [0, 15])
    output['significant_threshold_one'] = None
    with pytest.raises(ValueError, match='threshold-bound'):
        native_layout_significant_host_selection_floor(layout, output,
            price_bank(), {'cuda:0': profile()})
    output['significant_threshold_one'] = False
    output['significant_backend'] = 'device'
    with pytest.raises(ValueError, match='threshold-bound'):
        native_layout_significant_host_selection_floor(layout, output,
            price_bank(), {'cuda:0': profile()})


def test_sparse_and_dense_nonzero_prices_bound_mixed_short_tail(tmp_path):
    floor = source(tmp_path)(4, 9, 2)
    layout = native_layout_source_floor([
        part('tile', 'cuda:0', (0, 8), floor)],
        total_traits=8, reduction='significant', partition_axis='trait')
    output = native_layout_output_floor(layout,
        significant_backend='host', significant_threshold_one=False,
        retained_ranges={'tile': [3, 3]},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400.},
        output_bytes_per_second=200.)
    evidence = price_bank()
    evidence['prices']['flatnonzero_sparse']['call_cpu_seconds'] = 5e-6
    evidence['prices']['flatnonzero_dense']['call_cpu_seconds'] = 2e-6
    report = native_layout_significant_host_selection_floor(layout, output,
        evidence, {'cuda:0': profile()})
    actual_cpu = actual_dram = 0.
    for markers, keep in zip([2, 2, 1], [1, 0, 2]):
        work = host_significant_selection_work(markers, 8, keep,
            threshold_one=False, return_beta=True)
        steps = host_selection_service(work, evidence['prices'],
            cpu_fraction=.5, dram_bytes_per_second=1e8,
            host_serial_fraction=0.)
        actual_cpu += sum(row['seconds'] * row['resources']['cpu']
                          for row in steps)
        actual_dram += sum(row['seconds'] * row['resources']['dram']
                           for row in steps)
    assert report['total_cpu_seconds'][0] <= actual_cpu <= report['total_cpu_seconds'][1]
    assert (report['total_logical_dram_bytes'][0] <= actual_dram <=
            report['total_logical_dram_bytes'][1])
    regimes = report['partitions'][0]['possible_shapes'][0]['possible_nonzero_regimes']
    assert 'flatnonzero_sparse' in regimes and 'flatnonzero_dense' in regimes
