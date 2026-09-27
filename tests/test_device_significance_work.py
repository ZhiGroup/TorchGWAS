"""Declared nonzero extents retain source control flow and gather payloads."""
import pytest
from torchgwas.device_significance_work import device_significant_tensor_work
from torchgwas.reduced_output_work import significant_output_work


@pytest.mark.parametrize('shape', [(13, 7, 17), (5, 37, 11), (257, 4093, 1 << 20)])
@pytest.mark.parametrize('occupancy', ['empty', 'dense', 'one_each'])
def test_source_blocks_count_actual_copies_and_preserve_tail(shape, occupancy):
    b, k, limit = shape
    blocks = significant_output_work(40, b, k, b, backend='device', max_selection_cells=limit)['blocks']
    counts = [0 if occupancy == 'empty' else row['cells'] if occupancy == 'dense' else 1 for row in blocks]
    report = device_significant_tensor_work(40, b, k, counts, max_cells=limit)
    assert report['nonzero_calls'] == report['selection_blocks'] == len(blocks)
    # One packed int32 copy (row, trait, beta, t, df) per nonempty block.
    assert report['selected_copy_calls'] == sum(v > 0 for v in counts)
    assert report['selected_payload_d2h_bytes'] == 20 * sum(counts)
    assert [v['retained'] for v in report['blocks']] == counts
    assert report['retained_pairs'] == sum(counts)
    assert all(v['cells'] <= limit for v in report['blocks'])
    gathers = [v for v in report['steps'] if v['op'] == 'aten.index.Tensor']
    assert all(v['read_bytes'] == sum(t['bytes'] for t in v['outputs']) + sum(t['bytes'] for t in v['inputs'][1:]) for v in gathers)
    assert not report['prediction_complete']


@pytest.mark.parametrize('counts', [None, [], [True], [-1], [8]])
def test_missing_invalid_or_excess_survivors_are_refused(counts):
    with pytest.raises(ValueError):
        device_significant_tensor_work(40, 1, 7, counts)


def test_voxel_shape_requires_bounded_block_expansion():
    with pytest.raises(ValueError, match='max_blocks'):
        device_significant_tensor_work(22250, 1024, 2075298, [], max_cells=1 << 20, max_blocks=8)


def test_default_block_is_the_whole_chunk_below_the_nonzero_limit():
    from torchgwas.selection_geometry import DEVICE_SELECTION_MAX_CELLS
    report = device_significant_tensor_work(40, 1024, 8192, [3])
    assert report['selection_blocks'] == report['nonzero_calls'] == report['selected_copy_calls'] == 1
    assert report['blocks'][0]['cells'] == 1024 * 8192 <= DEVICE_SELECTION_MAX_CELLS
    # Predicate temporaries stay within the host predicate bound (1M cells):
    # eight 128-row strips, each one abs, one compare and one finite test.
    assert sum(v['op'] == 'aten.lt.Scalar' for v in report['steps']) == 8
    assert report['selected_payload_d2h_bytes'] == 60


def test_row_setup_is_reused_across_trait_windows():
    report = device_significant_tensor_work(40, 1, 37, [0] * 4, max_cells=11)
    # df cast/clamp/critical gather runs once, while each trait window selects.
    assert sum(v['op'] == 'aten.clamp.default' for v in report['steps']) == 1
    assert sum(v['op'] == 'aten.nonzero.default' for v in report['steps']) == 4
    assert len(report['blocks']) == 4


def test_nonzero_keeps_the_second_mask_scan_for_empty_results():
    from torchgwas.device_significance_work import cuda_nonzero_work
    empty = cuda_nonzero_work(100, 0, torch_version='2.5.1+cu124')
    assert [p['phase'] for p in empty['phases']] == ['count', 'count_to_host', 'flagged_select']
    assert empty['gpu_logical_bytes'] == 212
    assert empty['count_d2h_bytes'] == 4 and empty['output_bytes'] == 0
    dense = cuda_nonzero_work(100, 100, torch_version='2.5.1+cu124')
    assert dense['gpu_logical_bytes'] == 3412
    assert dense['output_bytes'] == 1600 and dense['flat_indices_alias_output']
    assert dense['phases'][-1]['grid'] == [1, 1, 1]


def test_nonzero_refuses_an_unverified_backend_or_overflow():
    from torchgwas.device_significance_work import cuda_nonzero_work
    with pytest.raises(ValueError, match='2.5.1'):
        cuda_nonzero_work(100, 1, torch_version='2.6.0')
    with pytest.raises(ValueError, match='INT_MAX'):
        cuda_nonzero_work((1 << 31)-1, 0, torch_version='2.5.1')
