"""The selector and calculator share the minimum-block bounded geometry."""
from torchgwas.selection_geometry import device_selection_shape
from torchgwas.reduced_output_work import significant_output_work


def test_shape_minimizes_blocking_counts_against_all_small_feasible_heights():
    for rows in range(1, 18):
        for traits in range(1, 25):
            for cap in range(1, 33):
                width, height, count = device_selection_shape(rows, traits, cap)
                assert 1 <= width <= traits and 1 <= height <= rows
                assert width * height <= cap
                brute = min(
                    ((rows + h - 1) // h) *
                    ((traits + min(traits, cap // h) - 1) // min(traits, cap // h))
                    for h in range(1, min(rows, cap) + 1))
                assert count == brute


def test_large_tile_reduces_count_barriers_without_expanding_blocks():
    rows, traits, cap = 512, 600000, 1 << 20
    width, height, count = device_selection_shape(rows, traits, cap)
    old_width = min(traits, cap)
    old_height = max(1, cap // old_width)
    old_count = ((rows + old_height - 1) // old_height) * (
        (traits + old_width - 1) // old_width)
    assert count < old_count
    assert (width, height, count) == (2048, 512, 293)
    work = significant_output_work(22250, 512 * 100 + 7, traits, rows,
        backend='device', max_selection_cells=cap, include_blocks=False)
    tail = device_selection_shape(7, traits, cap)
    assert work['selection_blocks'] == 100 * count + tail[2]
    assert work['selection_count_d2h_bytes'] == 4 * work['selection_blocks']
    assert work['maximum_selection_block_cells'] <= cap
