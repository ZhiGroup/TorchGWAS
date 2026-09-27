"""Only matching, complete unissued work can share one partial resource floor."""
from copy import deepcopy

import pytest

from test_layout_frontier import productive_snapshot
from test_source_layout_floor import part, source
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.layout_output_floor import native_layout_output_floor
from torchgwas.layout_partial_envelope import native_layout_partial_envelope
from torchgwas.source_layout_floor import native_layout_source_floor


def case(tmp_path):
    floor = source(tmp_path)(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    frontier = unissued_frontier(productive_snapshot(),
        source_identity=floor['input_identity'], reduction=None,
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.})
    output = native_layout_output_floor(layout,
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    return frontier, layout, compute, output


def test_complete_unissued_layout_combines_only_necessary_stage_floors(tmp_path):
    frontier, layout, compute, output = case(tmp_path)
    report = native_layout_partial_envelope(frontier, layout, compute, output)
    expected = [max(layout['source_stage_floor_seconds'][i],
        compute['compute_transfer_floor_seconds'], output['payload_floor_seconds'][i],
        (layout['resource_work']['dram_bytes'] + compute['total_h2d_bytes'] +
         output['total_d2h_payload_bytes'][i]) / layout['shared_capacities']['dram'])
        for i in (0, 1)]
    assert report['partial_floor_seconds'] == expected
    assert report['stage_floor_seconds']['combined_shared_dram'][0] > (
        layout['resource_work']['dram_bytes'] / layout['shared_capacities']['dram'])
    assert report['required_cells'] == 8 * 5
    assert report['coverage']['exact_coverage']
    assert not report['prediction_complete'] and not report['selection_validated']


def test_shared_dram_counts_decode_and_both_dma_directions_once(tmp_path):
    floor = source(tmp_path, capacities=dict(cpu=1., dram=100., input=5e6))(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    frontier = unissued_frontier(productive_snapshot(),
        source_identity=floor['input_identity'], reduction=None,
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1e9,
        per_device_h2d_bytes_per_second={'cuda:0': 1e9, 'cuda:1': 1e9},
        peak_fp32_flops_per_second={'cuda:0': 1e9, 'cuda:1': 1e9})
    output = native_layout_output_floor(layout,
        shared_d2h_bytes_per_second=1e9,
        per_device_d2h_bytes_per_second={'cuda:0': 1e9, 'cuda:1': 1e9},
        output_bytes_per_second=1e9)
    report = native_layout_partial_envelope(frontier, layout, compute, output)
    assert report['partial_floor_seconds'][0] == (
        layout['resource_work']['dram_bytes'] + compute['total_h2d_bytes'] +
        output['total_d2h_payload_bytes'][0]) / 100.
    assert report['partial_floor_seconds'][0] > layout['source_stage_floor_seconds'][0]


def test_empty_first_significant_part_still_binds_full_future_scenarios(tmp_path):
    floor = source(tmp_path)(4, 12, 2)
    layout = native_layout_source_floor([
        part('low', 'cuda:0', (0, 2), floor),
        part('high', 'cuda:1', (2, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    snapshot = productive_snapshot(reduction='significant')
    assert snapshot['first_written'] is not None and snapshot['first_material_part'] is None
    frontier = unissued_frontier(snapshot, source_identity=floor['input_identity'],
        reduction='significant', total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.})
    output = native_layout_output_floor(layout, significant_backend='host',
        retained_ranges={'low': [0, 16], 'high': [0, 24]},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    report = native_layout_partial_envelope(frontier, layout, compute, output)
    assert report['coverage']['exact_coverage']
    assert report['stage_floor_seconds']['d2h_and_output'][1] > 0
    assert not report['selection_validated']


def test_full_panel_jagwas_variant_shards_bind_reduced_output(tmp_path):
    floor = source(tmp_path)
    left, right = floor(4, 8, 2), floor(8, 12, 2)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), left),
        part('right', 'cuda:1', (0, 5), right)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    frontier = unissued_frontier(productive_snapshot(reduction='jagwas'),
        source_identity=left['input_identity'], reduction='jagwas',
        total_traits=5, job_variant_range=[0, 12])
    compute = native_layout_compute_floor(layout, covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
        peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.},
        peak_fp64_flops_per_second={'cuda:0': 1000., 'cuda:1': 1000.})
    output = native_layout_output_floor(layout,
        retained_ranges={'left': [0, 4], 'right': [0, 4]},
        shared_d2h_bytes_per_second=500.,
        per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
        output_bytes_per_second=200.)
    report = native_layout_partial_envelope(frontier, layout, compute, output)
    assert report['required_cells'] == 8 * 5
    assert report['stage_floor_seconds']['d2h_and_output'][0] > 0


@pytest.mark.parametrize('damage', ['source', 'compute_id', 'compute_samples',
                                    'output_rows', 'output_mode', 'missing_pair',
                                    'links'])
def test_partial_envelope_rejects_mixed_or_incomplete_work(tmp_path, damage):
    frontier, layout, compute, output = case(tmp_path)
    layout, compute, output = map(deepcopy, (layout, compute, output))
    if damage == 'source':
        compute['input_identity']['bytes'] += 1
    elif damage == 'compute_id':
        compute['partitions'][0]['id'] = 'other'
    elif damage == 'compute_samples':
        compute['samples'] += 1
    elif damage == 'output_rows':
        output['partitions'].pop()
    elif damage == 'output_mode':
        output['reduction'] = 'jagwas'
    elif damage == 'links':
        output['shared_d2h_link_loads'][0]['declarations'].append(
            dict(devices=['cuda:0'], h2d_bytes_per_second=1., d2h_bytes_per_second=1.))
    else:
        layout['partitions'][1]['trait_range'] = [3, 5]
    with pytest.raises(ValueError):
        native_layout_partial_envelope(frontier, layout, compute, output)
