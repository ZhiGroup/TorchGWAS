"""Price live accepted writer bytes without calling them a job checkpoint."""
from copy import deepcopy
import time

import numpy as np
import pytest

from torchgwas.productive_dense_writer_queue_service import (price_dense_writer_queue_observation,
    price_bracketed_dense_writer_queues)
from torchgwas.productive_dense_writer_queue_bracket import bracket_dense_writer_queues
from torchgwas.sumstats import BinarySumstatsWriter


def profile():
    return dict(cpu_fraction=.5,
                writer_copy_service=dict(cpu_seconds_per_byte=.01,
                                         cpu_seconds_per_call=.001),
                process_units=dict(bytearray_zero_bytes=.01),
                executor_cpu_seconds=.001,
                writeback_service=dict(pagecache_seconds_per_byte=.1,
                                       storage_seconds_per_byte=.2,
                                       submit_seconds=.01, wait_seconds=.01,
                                       fadvise_seconds=.01),
                fsync_seconds=.1)


def test_staged_live_writer_is_priced_from_existing_independent_service(tmp_path):
    writer = BinarySumstatsWriter(tmp_path, 2, ['a', 'b'], 129, 123,
        block_bytes=32, writeback_bytes=0, store_variant_df=True,
        on_write_progress=lambda event: None)
    try:
        values = np.ones((1, 2), np.float32)
        writer.write_chunk(0, 1, values, values, np.ones((1, 1), np.float32))
        observation = writer.queue_snapshot()
        priced = price_dense_writer_queue_observation(observation, profile())
        assert observation['atomic_writer_streams']
        assert set(priced['streams']) == {'beta', 't_stat', 'df'}
        assert priced['pending_write_bytes_interval'] == [20, 20]
        assert priced['pagecache_cpu_seconds_at_fixed_price'] == pytest.approx([2., 2.])
        assert priced['storage_seconds_at_fixed_price'] == pytest.approx([4., 4.])
        assert priced['streams']['df']['pending_write_bytes_interval'] == [4, 4]
        assert priced['streams']['df']['stream_serial_pagecache_seconds_at_fixed_fraction'] == pytest.approx([.8, .8])
        assert not priced['prediction_complete'] and not priced['selection_validated']
        writer.write_chunk(1, 2, values, values, np.ones((1, 1), np.float32))
        writer.close()
    finally:
        writer.abort()


def test_active_block_counts_only_in_upper_byte_endpoint(tmp_path):
    writer = BinarySumstatsWriter(tmp_path, 1, ['a'], 129, 123,
        block_bytes=4, writeback_bytes=0, on_write_progress=lambda event: None)
    try:
        row = dict(accepted_bytes=4, written_bytes=0, staging_bytes=0,
                   queued_bytes=0, active_bytes=4, error=False,
                   pending_write_bytes_interval=[0, 4])
        observed = writer.queue_snapshot()
        observed['streams'] = {'t_stat': row}
        priced = price_dense_writer_queue_observation(observed, profile())
        assert priced['pending_write_bytes_interval'] == [0, 4]
        assert priced['pagecache_cpu_seconds_at_fixed_price'] == pytest.approx([0, .4])
        for edit in (lambda q: q.update(atomic_writer_streams=False),
                     lambda q: q['streams']['t_stat'].update(pending_write_bytes_interval=[4, 4]),
                     lambda q: q['streams']['t_stat'].update(active_bytes=3)):
            broken = deepcopy(observed); edit(broken)
            with pytest.raises(ValueError):
                price_dense_writer_queue_observation(broken, profile())
    finally:
        writer.abort()


def test_two_writers_have_one_common_anchor_and_distinct_device_prices(tmp_path):
    writers = [BinarySumstatsWriter(tmp_path/name, 2, ['a'], 129, 123,
        block_bytes=32, writeback_bytes=0, on_write_progress=lambda event: None)
        for name in ('one', 'two')]
    try:
        for writer in writers:
            row = np.ones((1, 1), np.float32)
            writer.write_chunk(0, 1, row, row)
        keys = [str(writer.directory.resolve()) for writer in writers]
        def pass_once():
            return {key:dict(device='cuda:'+str(index),
                             observation=writer.queue_snapshot())
                    for index,(key,writer) in enumerate(zip(keys,writers))}
        first = pass_once()
        anchor = time.perf_counter()
        second = pass_once()
        bracket = bracket_dense_writer_queues(first, second, anchor)
        assert bracket['logical_pending_bytes_interval'] == [16, 16]
        assert bracket['os_write_bytes_upper'] == 16
        other = deepcopy(profile())
        other['writeback_service']['pagecache_seconds_per_byte'] = .2
        other['writeback_service']['storage_seconds_per_byte'] = .4
        priced = price_bracketed_dense_writer_queues(
            bracket, {'cuda:0':profile(), 'cuda:1':other})
        assert priced['pagecache_cpu_seconds_upper_at_fixed_price'] == pytest.approx(2.4)
        assert priced['storage_seconds_upper_at_fixed_price'] == pytest.approx(4.8)
        stale = deepcopy(second)
        stale[keys[0]]['observation']['streams']['beta']['accepted_bytes'] = 0
        with pytest.raises(ValueError):
            bracket_dense_writer_queues(first, stale, anchor)
        stale = deepcopy(second)
        stale[keys[1]]['observation']['capture_started_seconds'] = anchor - 1.
        with pytest.raises(ValueError):
            bracket_dense_writer_queues(first, stale, anchor)
    finally:
        for writer in writers:
            writer.abort()
