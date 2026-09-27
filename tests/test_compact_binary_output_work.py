"""Compact dense-writer counts must match the source-expanded ledger."""
import itertools

import pytest

from torchgwas.binary_output_work import binary_output_work, compact_binary_output_work


@pytest.mark.parametrize('markers,traits,chunk,block,borrow,df', [
    (markers, traits, chunk, block, borrow, df)
    for markers, traits, chunk, block, borrow, df in itertools.product(
        (0, 1, 2, 7, 15, 33), (1, 3), (1, 4, 8),
        (8, 32, 64), (False, True), (False, True))
    if markers <= 15 or chunk == 4
])
def test_compact_matches_expanded_fixed_chunks(markers, traits, chunk, block, borrow, df):
    kwargs = dict(markers=markers, traits=traits, chunk_markers=chunk,
                  block_bytes=block, queue_depth=2, borrow_chunks=borrow,
                  store_beta=True, fsync=True, writeback_bytes=32,
                  sync_file_range=True, store_variant_df=df)
    compact = compact_binary_output_work(**kwargs)
    expanded = binary_output_work(**kwargs)
    for short, long in (
        ('allocated_staging_bytes', 'allocated_staging_bytes'),
        ('zero_initialization_bytes', 'zero_initialization_bytes'),
        ('payload_bytes', 'binary_payload_bytes'),
        ('staging_copy_bytes', 'staging_copy_bytes'),
        ('staging_copy_calls', 'staging_copy_calls'),
        ('write_calls_minimum', 'binary_write_calls_minimum'),
        ('payload_queued_only_at_close_bytes', 'payload_queued_only_at_close_bytes'),
        ('fsync_calls', 'fsync_calls'),
    ):
        assert compact[short] == expanded[long], (kwargs, short)
    for name, traits_for_stream in [('beta', traits), ('t', traits)] + ([('df', 1)] if df else []):
        size = block if name != 'df' else min(block, 1 << 20)
        one = binary_output_work(markers, traits_for_stream, chunk, size, 2,
                                 borrow, False, True, 32, True)
        row = compact['streams'][name]
        assert row['payload_bytes'] == one['binary_payload_bytes']
        assert row['staging_copy_bytes'] == one['staging_copy_bytes']
        assert row['staging_copy_calls'] == one['staging_copy_calls']
        assert row['write_calls_minimum'] == one['binary_write_calls_minimum']
        assert row['writeback_submit_calls'] == one['writeback_per_array']['submit_calls']
        assert row['writeback_wait_calls'] == one['writeback_per_array']['wait_calls']
        assert row['writeback_unsubmitted_tail_bytes'] == one['writeback_per_array']['unsubmitted_tail_bytes']


def test_large_dense_output_is_compact():
    result = compact_binary_output_work(8_086_101, 128, 1024, block_bytes=16 << 20,
                                        borrow_chunks=False, store_variant_df=True)
    assert result['payload_bytes'] == 8_086_101 * (128 * 8 + 4)
    assert len(result['streams']) == 3
    assert result['write_calls_minimum'] > 0


@pytest.mark.parametrize('writeback,sync', [(0, True), (32, False), (32, True)])
def test_auto_blocks_and_optional_writeback(writeback, sync):
    settings = dict(markers=19, traits=3, chunk_markers=7,
                    block_bytes=None, borrow_chunks=False,
                    store_beta=False, store_variant_df=True,
                    writeback_bytes=writeback, sync_file_range=sync)
    compact = compact_binary_output_work(**settings)
    expanded = binary_output_work(**settings)
    assert compact['block_bytes'] == expanded['block_bytes']
    assert compact['payload_bytes'] == expanded['binary_payload_bytes']
    assert compact['staging_copy_calls'] == expanded['staging_copy_calls']
    assert compact['write_calls_minimum'] == expanded['binary_write_calls_minimum']
    assert compact['writeback_submit_calls'] == (
        sum(row['writeback_submit_calls'] for row in compact['streams'].values()))
    if not sync or not writeback:
        assert compact['writeback_submit_calls'] == 0


def test_invalid_compact_inputs():
    with pytest.raises(ValueError):
        compact_binary_output_work(10, writeback_bytes=-1)
    with pytest.raises(ValueError):
        compact_binary_output_work(10, store_beta=1)
