"""Only issued chunks without output completion belong to the held backlog."""
import queue
import time

import numpy as np
import pytest

from test_pgen_work_bounds import write_records
from torchgwas.layout_frontier import unissued_frontier
from torchgwas.pgen_reader import pack_genovec
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.productive_boundary import ProductiveBoundaryProgress
from torchgwas.productive_issued_work import productive_issued_work
from torchgwas.productive_indexed_result_queue import snapshot_indexed_result_queue
from torchgwas.productive_indexed_queue_join import bind_jagwas_result_queue
from torchgwas.productive_issued_queue_refinement import refine_issued_jagwas_with_queue
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.sumstats import DenseWriteProgress
from torchgwas.sumstats_indexed import IndexedChunkWrite, IndexedOutputPartition


def header_for(tmp_path):
    samples = 129
    path = tmp_path / 'issued.pgen'
    records = [pack_genovec(np.full(samples, i % 3, dtype=np.uint8),
                            samples).tobytes() for i in range(12)]
    write_records(path, samples, [0] * 12, records)
    return PgenHeaderWork(path)


def frontier_for(snapshot, header, reduction):
    return unissued_frontier(snapshot, source_identity=header.input_identity,
                             reduction=reduction, total_traits=5,
                             job_variant_range=[0, 12])


def test_dense_partial_prefix_replays_full_issued_chunks(tmp_path, monkeypatch):
    header = header_for(tmp_path)
    partitions = [dict(id='original', device='cuda:0',
                       variant_range=[0, 12], trait_range=[0, 5])]
    run = ProductiveTuningRun(partitions, chunk_sizes=[2, 4], initial=4)
    progress = ProductiveBoundaryProgress(partitions, None)
    reserve = run.for_partition('original')
    assert reserve(0, 12, 4) == 4
    assert reserve(4, 12, 4) == 4
    event = DenseWriteProgress(0, 2, (0, 5), 80, 2,
                               time.perf_counter(), 'sumstats', 'cuda:0')
    run.output_written(event)
    progress.observe(event)
    snapshot = run.snapshot()
    frontier = frontier_for(snapshot, header, None)
    boundary = progress.bind(snapshot)
    assert boundary['valid']
    report = productive_issued_work(boundary, frontier, header,
                                    covariate_rank=3)
    assert [row['variant_range'] for row in report['chunks']] == [[0, 4], [4, 8]]
    assert report['pending_source_records'] == 8
    assert report['total_work']['read_bytes'] == (
        header.bounds(0, 4)['read_bytes'] + header.bounds(4, 8)['read_bytes'])
    assert report['total_work']['h2d_bytes'] == 8 * 129
    assert report['total_work']['fp32_gemm_flops'] == 2 * 8 * 129 * 9
    with monkeypatch.context() as patch:
        patch.setattr(header, 'bounds',
                      lambda *args, **kwargs: pytest.fail('Budget checked after source read'))
        with pytest.raises(ValueError, match='pending-chunk|bounded budget'):
            productive_issued_work(boundary, frontier, header,
                                   covariate_rank=3, max_pending_chunks=1)
    run.finish(successful=False)


def test_jagwas_completed_empty_part_excludes_only_its_source_chunk(tmp_path):
    header = header_for(tmp_path)
    partitions = [
        dict(id='left', device='cuda:0', variant_range=[0, 6], trait_range=[0, 5]),
        dict(id='right', device='cuda:1', variant_range=[6, 12], trait_range=[0, 5])]
    run = ProductiveTuningRun(partitions, chunk_sizes=[3, 6], initial=3)
    progress = ProductiveBoundaryProgress(partitions, 'jagwas')
    assert run.for_partition('left')(0, 6, run.capacity) == 3
    assert run.for_partition('right')(6, 12, run.capacity) == 3
    owner = IndexedOutputPartition('cuda:0', (0, 6), (0, 5))
    now = time.perf_counter()
    event = IndexedChunkWrite(0, 3, 'jagwas', 0, 0, None, now, now,
                              False, owner, (0, 3))
    run.output_written(event)
    progress.observe(event)
    snapshot = run.snapshot()
    frontier = frontier_for(snapshot, header, 'jagwas')
    boundary = progress.bind(snapshot)
    report = productive_issued_work(boundary, frontier, header,
                                    covariate_rank=2)
    assert report['pending_chunks'] == 1
    assert report['chunks'][0]['variant_range'] == [6, 9]
    assert report['chunks'][0]['device'] == 'cuda:1'
    assert report['total_work']['fp64_projection_flops'] == 2 * 3 * 5 * 5
    stale = dict(boundary, written_events=boundary['written_events'] + 1)
    with pytest.raises(ValueError, match='bound source and output'):
        productive_issued_work(stale, frontier, header, covariate_rank=2)
    run.finish(successful=False)


def test_significant_trait_tiles_keep_only_unfinished_producer(tmp_path):
    header = header_for(tmp_path)
    partitions = [
        dict(id='low', device='cuda:0', variant_range=[0, 12], trait_range=[0, 2]),
        dict(id='high', device='cuda:1', variant_range=[0, 12], trait_range=[2, 5])]
    run = ProductiveTuningRun(partitions, chunk_sizes=[2, 4], initial=4)
    progress = ProductiveBoundaryProgress(partitions, 'significant')
    assert run.for_partition('low')(0, 12, run.capacity) == 4
    assert run.for_partition('high')(0, 12, run.capacity) == 4
    owner = IndexedOutputPartition('cuda:0', (0, 12), (0, 2))
    now = time.perf_counter()
    event = IndexedChunkWrite(0, 4, 'significant', 0, 0, None, now, now,
                              False, owner, (0, 4))
    run.output_written(event)
    progress.observe(event)
    snapshot = run.snapshot()
    frontier = frontier_for(snapshot, header, 'significant')
    report = productive_issued_work(progress.bind(snapshot), frontier, header,
                                    covariate_rank=2)
    assert report['pending_chunks'] == 1
    assert report['chunks'][0]['id'] == 'high'
    assert report['chunks'][0]['trait_range'] == [2, 5]
    assert report['total_work']['h2d_bytes'] == 4 * 129
    assert report['total_work']['fp32_gemm_flops'] == 2 * 4 * 129 * 6
    run.finish(successful=False)


def test_jagwas_queued_result_removes_completed_producer_work(tmp_path):
    header = header_for(tmp_path)
    partitions = [
        dict(id='left', device='cuda:0', variant_range=[0, 6], trait_range=[0, 5]),
        dict(id='right', device='cuda:1', variant_range=[6, 12], trait_range=[0, 5])]
    run = ProductiveTuningRun(partitions, chunk_sizes=[3, 6], initial=3)
    progress = ProductiveBoundaryProgress(partitions, 'jagwas')
    assert run.for_partition('left')(0, 6, run.capacity) == 3
    assert run.for_partition('right')(6, 12, run.capacity) == 3
    now = time.perf_counter()
    event = IndexedChunkWrite(0, 3, 'jagwas', 0, 0, None, now, now, False,
        IndexedOutputPartition('cuda:0', (0, 6), (0, 5)), (0, 3))
    run.output_written(event)
    progress.observe(event)
    held = run.snapshot()
    frontier = frontier_for(held, header, 'jagwas')
    boundary = progress.bind(held)
    issued = productive_issued_work(boundary, frontier, header, covariate_rank=2)
    source = queue.Queue(maxsize=2)
    finished = object()
    source.put((6, 9, None, np.arange(15, dtype=np.float64), None, None))
    observed = snapshot_indexed_result_queue(source, [(0, 6), (6, 12)],
        ['cuda:0', 'cuda:1'], finished)
    joined = bind_jagwas_result_queue(boundary, observed)
    refined = refine_issued_jagwas_with_queue(issued, joined)
    assert refined['queued_chunks_with_completed_producer'] == 1
    assert refined['upstream_or_unresolved_chunks'] == 0
    assert refined['queued_producer_work_completed'] == issued['total_work']
    assert all(value == 0 for value in
        refined['upstream_or_unresolved_full_chunk_work_upper'].values())
    assert refined['per_device_upstream_or_unresolved_upper'] == {}
    assert not refined['prediction_complete'] and not refined['selection_validated']
    stale = dict(joined, issued_revision=joined['issued_revision'] + 1)
    with pytest.raises(ValueError, match='Same-revision'):
        refine_issued_jagwas_with_queue(issued, stale)
    wrong = dict(joined, queued_source_chunks=[dict(joined['queued_source_chunks'][0],
        variant_range=[9, 12])])
    with pytest.raises(ValueError, match='Queued result differs'):
        refine_issued_jagwas_with_queue(issued, wrong)
    run.finish(successful=False)


def test_jagwas_active_writer_and_queue_exclude_both_producers(tmp_path):
    header = header_for(tmp_path)
    partitions = [
        dict(id='left', device='cuda:0', variant_range=[0, 6], trait_range=[0, 5]),
        dict(id='right', device='cuda:1', variant_range=[6, 12], trait_range=[0, 5])]
    run = ProductiveTuningRun(partitions, chunk_sizes=[3, 6], initial=3)
    progress = ProductiveBoundaryProgress(partitions, 'jagwas')
    left = run.for_partition('left')
    assert left(0, 6, run.capacity) == 3
    assert left(3, 6, run.capacity) == 3
    assert run.for_partition('right')(6, 12, run.capacity) == 3
    now = time.perf_counter()
    event = IndexedChunkWrite(0, 3, 'jagwas', 0, 0, None, now, now, False,
        IndexedOutputPartition('cuda:0', (0, 6), (0, 5)), (0, 3))
    run.output_written(event)
    progress.observe(event)
    held = run.snapshot()
    frontier = frontier_for(held, header, 'jagwas')
    boundary = progress.bind(held)
    issued = productive_issued_work(boundary, frontier, header, covariate_rank=2)
    source = queue.Queue(maxsize=2)
    source.put((6, 9, None, np.arange(15, dtype=np.float64), None, None))
    observed = snapshot_indexed_result_queue(source, [(0, 6), (6, 12)],
        ['cuda:0', 'cuda:1'], object())
    joined = bind_jagwas_result_queue(boundary, observed, active_writer=dict(
        device='cuda:0', variant_range=[3, 6], trait_range=[0, 5]))
    refined = refine_issued_jagwas_with_queue(issued, joined)
    assert refined['queued_chunks_with_completed_producer'] == 1
    assert refined['active_writer_chunks_with_completed_producer'] == 1
    assert refined['upstream_or_unresolved_chunks'] == 0
    assert refined['observed_producer_work_completed'] == issued['total_work']
    with pytest.raises(ValueError, match='not one unfinished issued'):
        bind_jagwas_result_queue(boundary, observed, active_writer=dict(
            device='cuda:1', variant_range=[6, 9], trait_range=[0, 5]))
    run.finish(successful=False)
