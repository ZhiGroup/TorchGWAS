"""Written prefixes, coalescing, concurrent arrays and GPU/tile ownership."""
from dataclasses import FrozenInstanceError
from pathlib import Path
import threading
import time
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.sumstats import BinarySumstatsWriter, DenseWriteProgress, open_binary_sumstats, open_binary_df
from torchgwas.sumstats_tiled import write_trait_tiled_sumstats
from torchgwas.sumstats_sharded import write_variant_sharded_sumstats
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.productive_boundary import ProductiveBoundaryProgress


@pytest.mark.parametrize('block', [10, 128, None])
@pytest.mark.parametrize('borrow', [False, True])
@pytest.mark.parametrize('beta', [False, True])
@pytest.mark.parametrize('variant_df', [False, True])
def test_only_written_common_prefixes_and_exact_data(tmp_path, block, borrow, beta, variant_df):
    n, k = 17, 3; values = np.arange(n*k, dtype=np.float32).reshape(n, k); events = []
    dfs = np.arange(n, dtype=np.float32).reshape(n, 1)+100
    def observe(event):
        assert isinstance(event, DenseWriteProgress)
        assert event.start == (events[-1].end if events else 0)
        assert event.start < event.end <= n and event.trait_range == (0, k)
        assert not (tmp_path/'manifest.json').exists()
        for name, expected in [('tstat.f32', values+1)] + ([('beta.f32', values)] if beta else []):
            assert (tmp_path/name).stat().st_size >= event.end*k*4
            actual = np.fromfile(tmp_path/name, dtype='<f4', count=event.end*k).reshape(event.end, k)
            np.testing.assert_array_equal(actual, expected[:event.end])
        assert event.variant_df_complete_to is None if not variant_df else 0 <= event.variant_df_complete_to <= n
        events.append(event)
    writer = BinarySumstatsWriter(tmp_path, n, list('abc'), 129, 123, block_bytes=block,
        queue_depth=1, writeback_bytes=0, borrow_chunks=borrow, store_beta=beta,
        store_variant_df=variant_df, on_write_progress=observe)
    for start in range(0, n, 4):
        stop = min(start+4, n)
        writer.write_chunk(start, stop, values[start:stop], (values+1)[start:stop],
                           dfs[start:stop] if variant_df else None)
    summary = writer.close()
    assert events[-1].end == n
    assert sum(event.rows for event in events) == n*k
    assert sum(event.statistic_bytes for event in events) == 4*n*k*(1+int(beta))
    assert all(a.completed <= b.completed for a, b in zip(events, events[1:]))
    if variant_df:
        np.testing.assert_array_equal(open_binary_df(tmp_path), dfs)
    stored_beta, stored_t, _ = open_binary_sumstats(tmp_path)
    np.testing.assert_array_equal(stored_t, values+1)
    if beta: np.testing.assert_array_equal(stored_beta, values)
    assert summary['payload_bytes'] == 4*n*k*(1+int(beta)) + (4*n if variant_df else 0)
    assert summary['progress_callback_seconds'] > 0
    for event in events:
        queue = event.writer_queue
        assert queue['valid'] and queue['kind'] == 'torchgwas.dense_writer_queue_observation.v1'
        assert queue['atomic_writer_streams']
        assert set(queue['streams']) == ({'t_stat'} | ({'beta'} if beta else set()) |
                                         ({'df'} if variant_df else set()))
        for state in queue['streams'].values():
            assert state['accepted_bytes'] - state['written_bytes'] == (
                state['staging_bytes'] + state['queued_bytes'] + state['active_bytes'])
            assert state['pending_write_bytes_interval'] == [
                state['staging_bytes'] + state['queued_bytes'],
                state['accepted_bytes'] - state['written_bytes']]
            assert state['queued_blocks'] >= 0
    final = writer.queue_snapshot()
    assert final['valid']
    assert all(row['accepted_bytes'] == row['written_bytes'] for row in final['streams'].values())
    for stream in (writer._beta, writer._tstat, writer._df):
        if stream is not None:
            assert not stream._thread.is_alive() and stream._on_written is None
    with pytest.raises(FrozenInstanceError): events[0].end = 99


def test_staging_and_one_completed_array_are_not_full_progress(tmp_path):
    from torchgwas import sumstats
    events = []; beta_written = threading.Event(); release_t = threading.Event(); progress = threading.Event()
    writer = BinarySumstatsWriter(tmp_path, 4, ['a', 'b'], 129, 123, block_bytes=32,
        writeback_bytes=0, on_write_progress=lambda event: (events.append(event), progress.set()))
    original = sumstats.os.write; syncs = []
    def write(fd, data):
        if fd == writer._tstat._fd:
            assert release_t.wait(5), 'test failed to release t writer'
        result = original(fd, data)
        if fd == writer._beta._fd: beta_written.set()
        return result
    try:
        with patch('torchgwas.sumstats.os.write', side_effect=write), \
             patch('torchgwas.sumstats.os.fsync', side_effect=lambda fd: syncs.append(fd)):
            values = np.ones((4, 2), np.float32)
            writer.write_chunk(0, 1, values[:1], values[:1])
            assert events == [] and writer._beta._offset == writer._tstat._offset == 0
            writer.write_chunk(1, 4, values[1:], values[1:])
            assert beta_written.wait(5) and events == []
            queue = writer.queue_snapshot()
            assert queue['valid'] and queue['streams']['t_stat']['active_bytes'] == 32
            release_t.set()
            assert progress.wait(5) and events[-1].end == 4 and syncs == []
            writer.close()
            assert syncs  # fsync policy remains at close, not at each notification.
    finally:
        release_t.set(); writer.abort()


def test_native_dense_progress_binds_queue_and_unique_producer(tmp_path):
    partitions=[dict(id='only',device='cuda:0',variant_range=[0,8],trait_range=[0,2])]
    run=ProductiveTuningRun(partitions,chunk_sizes=[2,4],initial=4)
    boundary=ProductiveBoundaryProgress(partitions,None)
    control=run.for_partition('only')
    def observe(event):
        run.output_written(event)
        boundary.observe(event)
    writer=BinarySumstatsWriter(tmp_path,8,['a','b'],129,123,block_bytes=32,
                                writeback_bytes=0,on_write_progress=observe)
    values=np.ones((4,2),np.float32)
    for first in (0,4):
        assert control(first,8,4)==4
        writer.write_chunk(first,first+4,values,values)
    writer.close()
    observed=boundary.bind(run.snapshot())
    assert observed['valid'] and observed['partitions'][0]['matrix_written_to']==8
    queue=observed['writer_queue_observation']
    assert queue['written_event']==run.snapshot()['written_events']
    assert queue['observation']['valid']
    assert all(row['accepted_bytes']==row['written_bytes']
               for row in queue['observation']['streams'].values())
    run.finish(successful=False)


def test_short_writes_do_not_report_incomplete_rows(tmp_path):
    from torchgwas import sumstats
    original = sumstats.os.write; writes = []; events = []
    def short(fd, data):
        result = original(fd, data[:3]); writes.append(result); return result
    writer = BinarySumstatsWriter(tmp_path, 3, ['a', 'b'], 129, 123, block_bytes=11,
                                  on_write_progress=events.append, writeback_bytes=0)
    with patch('torchgwas.sumstats.os.write', side_effect=short):
        values = np.ones((3, 2), np.float32)
        writer.write_chunk(0, 3, values, values); writer.close()
    assert len(writes) > 8 and events[-1].end == 3
    assert sum(e.rows for e in events) == 6
    assert all(e.writer_queue['valid'] for e in events)


def test_failed_callback_prevents_manifest_and_releases_writer_cycles(tmp_path):
    events = []
    def fail(event):
        events.append(event); raise RuntimeError('intentional callback failure')
    writer = BinarySumstatsWriter(tmp_path, 8, ['a'], 129, 123, block_bytes=8,
                                  on_write_progress=fail, writeback_bytes=0)
    try:
        with pytest.raises(RuntimeError, match='sumstats write'):
            writer.write_chunk(0, 8, np.ones((8, 1)), np.ones((8, 1))); writer.close()
    finally:
        writer.abort()
    assert len(events) == 1 and not (tmp_path/'manifest.json').exists()
    assert all(not s._thread.is_alive() and s._on_written is None for s in (writer._beta, writer._tstat))


def test_callback_disabled_does_not_install_stream_hooks(tmp_path):
    writer = BinarySumstatsWriter(tmp_path, 2, ['a'], 129, 123, block_bytes=4)
    with patch.object(writer, '_stream_written', side_effect=AssertionError('disabled progress')):
        writer.write_chunk(0, 2, np.ones((2, 1)), np.ones((2, 1))); summary=writer.close()
    assert writer._progress_lock is None
    assert writer._beta._queue_state_lock is None
    assert writer._tstat._queue_state_lock is None
    assert summary['progress_callback_seconds'] == 0


@pytest.mark.parametrize('axis', ['trait', 'variant'])
def test_partition_writers_label_real_global_ranges(tmp_path, axis):
    n, k = 11, 5; values = np.arange(n*k, dtype=np.float32).reshape(n, k); events = []
    def chunks(start, stop, first, width):
        for lo in range(start, stop, 3):
            hi = min(lo+3, stop); block = values[lo:hi, first:first+width]
            yield lo-start, hi-start, block, block+1, None, np.full((hi-lo, 1), 123, np.float32)
    opened = []
    kwargs = dict(n_variants=n, trait_names=list('abcde'), n_samples=129, df=123,
                  devices=['cuda:0', 'cuda:1'], reader_workers=2, block_bytes=24,
                  on_write_progress=events.append,
                  on_writer_open=lambda writer, device: opened.append((str(writer.directory), device)))
    if axis == 'trait':
        result = write_trait_tiled_sumstats(tmp_path, trait_block=2,
            scan_factory=lambda first, width, device, workers: chunks(0, n, first, width), **kwargs)
        partitions = [(tuple(row['trait_range']), (0, n), row['device']) for row in result['tiles']]
    else:
        result = write_variant_sharded_sumstats(tmp_path, chunk_size=4,
            scan_factory=lambda start, stop, device, workers: chunks(start, stop, 0, k), **kwargs)
        partitions = [((0, k), tuple(row['variant_range']), row['device']) for row in result['shards']]
    for traits, span, device in partitions:
        matching = sorted([event for event in events if event.trait_range == traits and event.device == device
                           and span[0] <= event.start < event.end <= span[1]], key=lambda e: e.start)
        cursor = span[0]
        for event in matching:
            assert event.start == cursor; cursor = event.end
        assert cursor == span[1]
    assert all(event.writer_queue['valid'] for event in events)
    assert len(opened) == len(partitions)
    assert {device for _, device in opened} == {device for _, _, device in partitions}
    assert sum(e.rows for e in events) == n*k
    b, t, _ = open_binary_sumstats(tmp_path)
    np.testing.assert_array_equal(b, values); np.testing.assert_array_equal(t, values+1)


def test_dense_progress_starts_tuning_without_claiming_chunk_or_fsync(tmp_path):
    run = ProductiveTuningRun([dict(id='0', device='cuda:0', variant_range=[0, 20], trait_range=[0, 3])],
                             chunk_sizes=[2, 4], initial=2)
    control = run.for_partition('0'); assert control(0, 20, 4) == 2
    before = run.planning_step(lambda state: pytest.fail('before output'), remaining_seconds=10.,
                              expected_cpu_seconds=.001, expected_wall_seconds=.002)
    assert before['reason'] == 'no_useful_output_yet'
    now = time.perf_counter()
    event = DenseWriteProgress(0, 1, (0, 3), 24, 0, now, str(tmp_path))
    run.output_written(event)
    step = run.planning_step(lambda state: dict(chunk_size=4, baseline_seconds=10., candidate_seconds=5.),
                            remaining_seconds=10., expected_cpu_seconds=.001, expected_wall_seconds=.002)
    assert step['applied'] and control(2, 20, 4) == 4
    snapshot = run.snapshot()
    assert snapshot['first_written'] == snapshot['first_material_output'] == now
    assert snapshot['first_fsynced_part'] is snapshot['first_material_part'] is None
    assert snapshot['written_chunks'] == 0 and snapshot['written_events'] == 1
    assert snapshot['dense_statistic_bytes'] == 24 and snapshot['part_bytes'] == 0
    run.finish(successful=False)
