"""Payload floors must distinguish dense, selected and JAGWAS execution."""
from copy import deepcopy

import pytest

from test_source_layout_floor import part, source
from torchgwas.source_layout_floor import native_layout_source_floor
from torchgwas.layout_output_floor import native_layout_output_floor


def priced(layout, **kwargs):
    devices = {row['device'] for row in layout['partitions']}
    settings = dict(
        shared_d2h_bytes_per_second=200.,
        per_device_d2h_bytes_per_second={device: 100. for device in devices},
        output_bytes_per_second=80.)
    settings.update(kwargs)
    return native_layout_output_floor(layout, **settings)


def test_dense_tiles_count_each_result_transfer_but_one_statistic_matrix(tmp_path):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([
        part('a', 'cuda:0', (0, 3), floor),
        part('b', 'cuda:1', (3, 5), floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    report = priced(layout)
    assert report['total_d2h_payload_bytes'] == [8 * 12 * 5 + 5 * 12 * 2] * 2
    assert report['total_output_array_payload_bytes'] == [8 * 12 * 5] * 2
    assert report['per_device_d2h_payload_bytes']['cuda:0'] == [8 * 12 * 3 + 5 * 12] * 2
    assert report['payload_floor_seconds'] == [480 / 80] * 2
    assert not report['selection_validated']


def test_dense_writer_counts_df_and_fixed_chunk_writeback_without_expansion(tmp_path):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([
        part('a', 'cuda:0', (0, 3), floor),
        part('b', 'cuda:1', (3, 5), floor)],
        total_traits=5, reduction=None, partition_axis='trait')
    writer = dict(block_bytes=32, queue_depth=2, borrow_chunks=False,
                  fsync=True, writeback_bytes=64, sync_file_range=True,
                  store_variant_df=True)
    report = priced(layout, dense_writer_options=writer)
    assert report['total_output_array_payload_bytes'] == [12 * 12 * 5 + 4 * 12 * 2] * 2
    assert report['payload_floor_seconds'] == [(12 * 12 * 5 + 4 * 12 * 2) / 80] * 2
    assert len(report['dense_writer_work']) == 2
    assert sum(row['work']['payload_bytes'] for row in report['dense_writer_work']) == 816
    assert all(row['work']['fsync_calls'] == 4 for row in report['dense_writer_work'])


def test_selected_output_rejects_dense_writer_settings(tmp_path):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    with pytest.raises(ValueError):
        priced(layout, significant_backend='host', dense_writer_options={})


def test_shared_d2h_link_changes_multi_gpu_payload_floor(tmp_path):
    floor = source(tmp_path)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), floor(0, 6)),
        part('right', 'cuda:1', (0, 5), floor(6, 12))],
        total_traits=5, reduction=None, partition_axis='variant')
    together = priced(layout, shared_links=[
        dict(devices=['cuda:0', 'cuda:1'], h2d_bytes_per_second=20.,
             d2h_bytes_per_second=20.)])
    apart = priced(layout, shared_links=[
        dict(devices=['cuda:0'], h2d_bytes_per_second=20., d2h_bytes_per_second=20.),
        dict(devices=['cuda:1'], h2d_bytes_per_second=20., d2h_bytes_per_second=20.)])
    assert together['shared_d2h_link_loads'][0]['links'][0]['bytes'] == 8 * 12 * 5 + 5 * 12
    assert together['payload_floor_seconds'] == [27.] * 2
    assert apart['payload_floor_seconds'] == [13.5] * 2


def test_significant_host_and_device_selector_have_distinct_transport(tmp_path):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([
        part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    retained = {'a': [2, 7]}
    host = priced(layout, significant_backend='host', retained_ranges=retained)
    device = priced(layout, significant_backend='device', retained_ranges=retained)
    assert host['total_d2h_payload_bytes'] == [8 * 60 + 5 * 12] * 2
    assert device['total_d2h_payload_bytes'] == [12 + 4 * 4 + 28 * 2,
                                                   12 + 4 * 4 + 28 * 7]
    assert host['total_output_array_payload_bytes'] == [56, 196]
    assert device['total_output_array_payload_bytes'] == host['total_output_array_payload_bytes']
    assert device['payload_floor_seconds'][0] < host['payload_floor_seconds'][0]
    unknown = priced(layout, significant_backend='device')
    assert unknown['total_output_array_payload_bytes'] == [0, 28 * 60]
    assert unknown['total_d2h_payload_bytes'] == [12 + 4 * 4,
                                                   12 + 4 * 4 + 28 * 60]


def test_device_selection_count_transfer_uses_exact_chunk_and_tile_geometry(tmp_path):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([
        part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    report = priced(layout, significant_backend='device',
                    device_selection_max_cells=7, retained_ranges={'a': [0, 0]})
    assert report['device_selection_max_cells'] == 7
    assert report['partitions'][0]['device_selection_blocks'] == 12
    assert report['partitions'][0]['maximum_selection_block_cells'] == 5
    assert report['total_d2h_payload_bytes'] == [12 + 4 * 12] * 2
    with pytest.raises(ValueError, match='device selection cell limit'):
        priced(layout, significant_backend='device', device_selection_max_cells=0)
    with pytest.raises(ValueError, match='device selection cell limit'):
        priced(layout, significant_backend='host', device_selection_max_cells=7)


def test_jagwas_variant_shards_keep_full_panel_and_count_selected_rows(tmp_path):
    floor = source(tmp_path)
    layout = native_layout_source_floor([
        part('left', 'cuda:0', (0, 5), floor(0, 6)),
        part('right', 'cuda:1', (0, 5), floor(6, 12))],
        total_traits=5, reduction='jagwas', partition_axis='variant')
    report = priced(layout, retained_ranges={'left': [1, 4], 'right': [0, 3]})
    assert report['total_d2h_payload_bytes'] == [17 * 12] * 2
    assert report['total_output_array_payload_bytes'] == [16, 112]
    assert report['per_device_d2h_payload_bytes'] == {
        'cuda:0': [17 * 6] * 2, 'cuda:1': [17 * 6] * 2}
    assert report['payload_floor_seconds'] == [102 / 100, 112 / 80]


@pytest.mark.parametrize('damage', ['backend', 'retained', 'capacity', 'stale',
                                    'chunk', 'beta', 'mode'])
def test_output_floor_rejects_unbound_or_changed_layout(tmp_path, damage):
    floor = source(tmp_path)(0, 12)
    layout = native_layout_source_floor([part('a', 'cuda:0', (0, 5), floor)],
        total_traits=5, reduction='significant', partition_axis='trait')
    kwargs = dict(significant_backend='device', retained_ranges={'a': [0, 60]})
    if damage == 'backend':
        kwargs['significant_backend'] = None
    elif damage == 'retained':
        kwargs['retained_ranges'] = {'a': [0, 61]}
    elif damage == 'capacity':
        kwargs['per_device_d2h_bytes_per_second'] = {'cuda:1': 100.}
    elif damage == 'stale':
        from pathlib import Path
        path = Path(layout['input_identity']['path'])
        path.write_bytes(path.read_bytes() + b'changed')
    elif damage == 'chunk':
        layout = deepcopy(layout)
        layout['partitions'][0]['chunks'] += 1
    elif damage == 'beta':
        kwargs['store_beta'] = 'yes'
    else:
        layout = deepcopy(layout)
        layout['reduction'] = 'jagwas'
    with pytest.raises(ValueError):
        priced(layout, **kwargs)
