"""Native int8 transfer and matrix-product floors follow fixed GPU partitions."""
from copy import deepcopy

import pytest

from test_source_layout_floor import part, source
from torchgwas.layout_compute_floor import native_layout_compute_floor
from torchgwas.source_layout_floor import native_layout_source_floor


def priced(layout, **updates):
    devices = {row['device'] for row in layout['partitions']}
    settings = dict(covariate_rank=3,
        shared_h2d_bytes_per_second=1000.,
        per_device_h2d_bytes_per_second={d: 800. for d in devices},
        peak_fp32_flops_per_second={d: 100000. for d in devices})
    if layout['reduction'] == 'jagwas':
        settings['peak_fp64_flops_per_second'] = {d: 1000. for d in devices}
    settings.update(updates)
    return native_layout_compute_floor(layout, **settings)


def test_trait_tiling_rereads_int8_source_and_counts_intercept(tmp_path):
    floor = source(tmp_path)(0, 12)
    partitions = [part('a', 'cuda:0', (0, 3), floor),
                  part('b', 'cuda:0', (3, 5), floor)]
    serial = native_layout_source_floor(partitions,
        total_traits=5, reduction=None, partition_axis='trait')
    assert serial['samples'] == 129
    result = priced(serial)
    assert result['total_h2d_bytes'] == 2 * 12 * 129
    assert result['per_device_work']['cuda:0']['fp32_gemm_flops'] == (
        2 * 12 * 129 * ((3 + 3 + 1) + (2 + 3 + 1)))
    assert result['h2d_payload_floor_seconds'] == 2 * 12 * 129 / 800
    assert result['matrix_product_floor_seconds'] == (
        result['per_device_work']['cuda:0']['fp32_gemm_flops'] / 100000)
    assert not result['selection_validated']
    partitions[1]['device'] = 'cuda:1'
    parallel = native_layout_source_floor(partitions,
        total_traits=5, reduction=None, partition_axis='trait')
    other = priced(parallel)
    assert other['total_h2d_bytes'] == result['total_h2d_bytes']
    assert other['matrix_product_floor_seconds'] < result['matrix_product_floor_seconds']
    assert other['h2d_payload_floor_seconds'] == 2 * 12 * 129 / 1000


def test_jagwas_shards_duplicate_full_panel_setup_width_and_project(tmp_path):
    floor = source(tmp_path)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), floor(0, 6)),
        part('right', 'cuda:1', (0, 5), floor(6, 12))],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    report = priced(layout)
    assert report['total_h2d_bytes'] == 12 * 129
    assert report['per_device_work']['cuda:0']['fp32_gemm_flops'] == 2 * 6 * 129 * 9
    assert report['per_device_work']['cuda:1']['fp64_projection_flops'] == 2 * 6 * 25
    assert report['matrix_product_floor_seconds'] == 2 * 6 * 129 * 9 / 100000 + 2 * 6 * 25 / 1000


def test_shared_h2d_link_changes_multi_gpu_necessary_floor(tmp_path):
    floor = source(tmp_path)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), floor(0, 6)),
        part('right', 'cuda:1', (0, 5), floor(6, 12))],
        total_traits=5, reduction=None, partition_axis='variant')
    together = priced(layout, shared_links=[
        dict(devices=['cuda:0', 'cuda:1'], h2d_bytes_per_second=300.,
             d2h_bytes_per_second=300.)])
    apart = priced(layout, shared_links=[
        dict(devices=['cuda:0'], h2d_bytes_per_second=300., d2h_bytes_per_second=300.),
        dict(devices=['cuda:1'], h2d_bytes_per_second=300., d2h_bytes_per_second=300.)])
    assert together['shared_h2d_link_loads']['links'][0]['bytes'] == 12 * 129
    assert together['h2d_payload_floor_seconds'] == 12 * 129 / 300
    assert apart['h2d_payload_floor_seconds'] == 6 * 129 / 300
    assert together['total_h2d_bytes'] == apart['total_h2d_bytes']


@pytest.mark.parametrize('damage', ['packed', 'rank', 'missing_gpu', 'fp64',
                                    'sample', 'panel', 'stale'])
def test_compute_floor_rejects_unbound_or_changed_execution(tmp_path, damage):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    settings = {}
    if damage == 'packed':
        settings['transport'] = 'pgen_2bit'
    elif damage == 'rank':
        settings['covariate_rank'] = 127
    elif damage == 'missing_gpu':
        settings['per_device_h2d_bytes_per_second'] = {'cuda:1': 800.}
    elif damage == 'fp64':
        settings['peak_fp64_flops_per_second'] = None
    elif damage == 'sample':
        layout = deepcopy(layout)
        layout['samples'] = 0
    elif damage == 'panel':
        layout = deepcopy(layout)
        layout['partitions'][0]['trait_range'] = [0, 4]
    else:
        from pathlib import Path
        path = Path(layout['input_identity']['path'])
        path.write_bytes(path.read_bytes() + b'changed')
    with pytest.raises(ValueError):
        priced(layout, **settings)
