"""Full/tail GPU service must match the actual unissued tile ownership."""
from copy import deepcopy
from unittest.mock import patch

import pytest

from test_productive_source_floor import setup
from torchgwas.layout_gpu_shape_service import native_layout_gpu_shape_service
from torchgwas.productive_source_floor import (
    productive_partial_floor, productive_source_floor)


def _profile(samples, rank, widths, sizes, reduction):
    profile = dict(reduction=reduction, cpu_fraction=.5,
                   gpu_resources=dict(hbm_bytes_per_second=1e12,
                       l2_bytes_per_second=2e12,
                       fp32_flops_per_second=1e13,
                       kernel_launch_seconds=1e-6,
                       host_dispatch_cpu_seconds=2e-6,
                       available_l2_bytes=1 << 20, sm_count=100),
                   host_primitives={'synthetic': 1e-6},
                   kernel_geometry=[dict(N=samples, B=size, K=width, C=rank,
                                         kernels=[])
                       for width in widths for size in sizes])
    if reduction == 'jagwas':
        profile['joint_host_primitives'] = {'synthetic': 1e-6}
        profile['joint_kernel_geometry'] = [
            dict(N=samples, B=size, K=width, compute_dtype='float32',
                 kernels=[])
            for width in widths for size in sizes]
    return profile


def _component(samples, size, width, rank, profile, gpu, statistics, joint):
    assert statistics['N'] == samples and statistics['B'] == size
    assert statistics['K'] == width and statistics['C'] == rank
    assert (joint is not None) == (profile['reduction'] == 'jagwas')
    return dict(kernel_service_seconds=.1 * size * width,
                host_dispatch_cpu_seconds=.01 * size,
                host_dispatch_serial_cpu_seconds=.002 * size,
                source_sha256='test-source', unpriced_terms=[])


def test_trait_tiles_and_tail_use_distinct_shapes_without_chunk_graph(tmp_path):
    header, source_profile, frontier = setup(tmp_path)
    partitions = [
        dict(id='low', device='cuda:0', variant_range=[4, 12],
             trait_range=[0, 2]),
        dict(id='high', device='cuda:1', variant_range=[4, 12],
             trait_range=[2, 5])]
    source = productive_source_floor(frontier, header, partitions,
        chunk_markers=3, profiles={row['id']: source_profile for row in partitions},
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis='trait')['source']
    profiles = {device: _profile(129, 3, [width], [3, 2], None)
                for device, width in [('cuda:0', 2), ('cuda:1', 3)]}
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component) as priced:
        report = native_layout_gpu_shape_service(
            source, covariate_rank=3, profiles=profiles)
    assert priced.call_count == 4
    assert report['distinct_shapes'] == 4
    assert report['per_device_service']['cuda:0']['chunks'] == 3
    assert report['per_device_service']['cuda:1']['chunks'] == 3
    assert report['per_device_service']['cuda:0']['kernel_service_seconds'] == pytest.approx(1.6)
    assert report['per_device_service']['cuda:1']['kernel_service_seconds'] == pytest.approx(2.4)
    assert report['total_host_work']['host_dispatch_cpu_seconds'] == pytest.approx(.16)
    assert report['calculation_wall_seconds'] >= 0
    assert not report['prediction_complete']

    with pytest.raises(ValueError, match='budget'):
        native_layout_gpu_shape_service(source, covariate_rank=3,
                                        profiles=profiles, max_shapes=3)
    stale = deepcopy(profiles)
    stale['cuda:1']['kernel_geometry'] = stale['cuda:1']['kernel_geometry'][:1]
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component), pytest.raises(
                   ValueError, match='statistics geometry'):
        native_layout_gpu_shape_service(source, covariate_rank=3,
                                        profiles=stale)


def test_jagwas_shape_service_keeps_full_panel_on_each_gpu(tmp_path):
    header, source_profile, frontier = setup(tmp_path, reduction='jagwas')
    partitions = [
        dict(id='left', device='cuda:0', variant_range=[4, 8],
             trait_range=[0, 5]),
        dict(id='right', device='cuda:1', variant_range=[8, 12],
             trait_range=[0, 5])]
    source = productive_source_floor(frontier, header, partitions,
        chunk_markers=3, profiles={row['id']: source_profile for row in partitions},
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis='variant')['source']
    profiles = {device: _profile(129, 2, [5], [3, 1], 'jagwas')
                for device in ('cuda:0', 'cuda:1')}
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component) as priced:
        report = native_layout_gpu_shape_service(
            source, covariate_rank=2, profiles=profiles)
    assert priced.call_count == 4
    assert [row['traits'] for row in report['partitions']] == [5, 5]
    assert all(row['chunks'] == 2 for row in report['partitions'])
    changed = deepcopy(source)
    changed['partitions'][0]['trait_range'] = [0, 4]
    with pytest.raises(ValueError, match='JAGWAS panel'):
        native_layout_gpu_shape_service(changed, covariate_rank=2,
                                        profiles=profiles)


@pytest.mark.parametrize('scan_mode', [None, 'device_significant'])
def test_significant_output_accepts_its_two_real_scan_modes(tmp_path,
                                                               scan_mode):
    header, source_profile, frontier = setup(tmp_path,
                                              reduction='significant')
    partition = dict(id='all', device='cuda:0', variant_range=[4, 12],
                     trait_range=[0, 5])
    source_report = productive_source_floor(frontier, header, [partition],
        chunk_markers=4, profiles={'all': source_profile},
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis='trait')
    source = source_report['source']
    profile = _profile(129, 3, [5], [4], scan_mode)
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component):
        report = native_layout_gpu_shape_service(source,
            covariate_rank=3, profiles={'cuda:0': profile})
    assert report['per_device_service']['cuda:0']['chunks'] == 2
    assert report['reduction'] == 'significant'
    with pytest.raises(ValueError, match='selector backend'):
        productive_partial_floor(frontier, source_report,
            compute_options=dict(covariate_rank=3,
                shared_h2d_bytes_per_second=1000.,
                per_device_h2d_bytes_per_second={'cuda:0': 800.},
                peak_fp32_flops_per_second={'cuda:0': 100000.}),
            output_options=dict(store_beta=True,
                significant_backend=('device' if scan_mode is None else 'host'),
                shared_d2h_bytes_per_second=500.,
                per_device_d2h_bytes_per_second={'cuda:0': 400.},
                output_bytes_per_second=200.),
            occupancy_scenario='empty',
            shape_profiles={'cuda:0': profile})
