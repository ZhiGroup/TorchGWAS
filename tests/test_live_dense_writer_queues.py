"""Live dense queue capture detects issue/output races without blocking writers."""
import threading
from types import SimpleNamespace

import numpy as np

from torchgwas.initial_chunk_autotune import PublicInitialChunkTuning
from torchgwas.productive_dense_writer_queue_service import price_bracketed_dense_writer_queues
from torchgwas.sumstats import BinarySumstatsWriter


def test_live_capture_has_current_staging_and_rejects_advanced_issue(tmp_path):
    cfg = dict(chunk_size=4, window_markers=[4, 8, 12],
               budget=dict(max_steps=1, max_cpu_seconds=1., max_window_seconds=10.))
    owner = SimpleNamespace(config=dict(initial_chunks=cfg), reduction=None,
                            refresh=None)
    partition = dict(id='only', device='cuda:0', variant_range=[0, 12],
                     trait_range=[0, 2])
    start = dict(partitions=[partition], chunk_sizes=[4], initial_size=4,
                 memory=dict(retained_index_bases_bytes=0))
    tuner = PublicInitialChunkTuning(owner, start, {}, {})
    tuner._planning_done = True
    control = tuner.for_partition('cuda:0', [0, 12], [0, 2])
    first_output = threading.Event()

    def observe(event):
        tuner.output_written(event)
        first_output.set()

    writer = BinarySumstatsWriter(tmp_path, 12, ['a', 'b'], 129, 123,
                                  block_bytes=32, writeback_bytes=0,
                                  on_write_progress=observe)
    tuner.register_writer(writer, 'cuda:0')
    try:
        values = np.ones((4, 2), np.float32)
        assert control(0, 12, 4) == 4
        writer.write_chunk(0, 4, values, values)
        assert first_output.wait(5.)
        assert control(4, 12, 4) == 4
        writer.write_chunk(4, 5, values[:1], values[:1])
        live = tuner.capture_live_writer_queues()
        assert live['observation_valid'] and live['stable_revision']
        assert live['registered_writers'] == 1
        assert live['bracket']['logical_pending_bytes_interval'] == [16, 16]
        only = next(iter(live['writers'].values()))
        assert only['device'] == 'cuda:0'
        assert only['first']['observation']['streams']['beta']['staging_bytes'] == 8
        assert only['second']['observation']['streams']['t_stat']['staging_bytes'] == 8
        original = writer.queue_snapshot
        calls = 0

        def accept_during_capture():
            nonlocal calls
            calls += 1
            if calls == 2:
                writer.write_chunk(5, 6, values[:1], values[:1])
            return original()

        writer.queue_snapshot = accept_during_capture
        changing = tuner.capture_live_writer_queues()
        assert changing['observation_valid'] and changing['stable_revision']
        assert changing['bracket']['logical_pending_bytes_interval'] == [16, 32]
        assert changing['bracket']['os_write_bytes_upper'] == 32
        profile = dict(cpu_fraction=.5,
            writer_copy_service=dict(cpu_seconds_per_byte=0.,cpu_seconds_per_call=0.),
            process_units=dict(bytearray_zero_bytes=0.),executor_cpu_seconds=0.,
            writeback_service=dict(pagecache_seconds_per_byte=.1,
                storage_seconds_per_byte=.2,submit_seconds=0.,wait_seconds=0.,
                fadvise_seconds=0.),fsync_seconds=0.)
        priced = price_bracketed_dense_writer_queues(changing['bracket'],
                                                     {'cuda:0': profile})
        assert priced['os_write_bytes_upper'] == 32
        assert priced['pagecache_cpu_seconds_upper_at_fixed_price'] == 3.2
        assert priced['storage_seconds_upper_at_fixed_price'] == 6.4
        advanced = False

        def advance_issue_during_capture():
            nonlocal advanced
            if not advanced:
                assert control(8, 12, 4) == 4
                advanced = True
            return original()

        writer.queue_snapshot = advance_issue_during_capture
        stale = tuner.capture_live_writer_queues()
        assert not stale['observation_valid'] and not stale['stable_revision']
        assert stale['issue_output_token']['issued_revision'] + 1 == tuner.run.revision_token()['issued_revision']
    finally:
        writer.abort()
        tuner.finish(successful=False)
