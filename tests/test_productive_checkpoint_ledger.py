"""A candidate may move unissued work, never its already-issued chunks."""
from copy import deepcopy
import queue
import time
from unittest.mock import patch

import numpy as np
import pytest

from test_productive_issued_work import header_for, frontier_for
from test_layout_gpu_shape_service import _component, _profile
from torchgwas.layout_gpu_shape_service import productive_issued_gpu_shape_service
from torchgwas.productive_boundary import ProductiveBoundaryProgress
from torchgwas.productive_checkpoint_ledger import productive_checkpoint_ledger
from torchgwas.productive_issued_work import productive_issued_work
from torchgwas.productive_indexed_result_queue import snapshot_indexed_result_queue
from torchgwas.productive_indexed_queue_join import bind_jagwas_result_queue
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.productive_source_floor import (
    productive_partial_floor, productive_source_floor)
from torchgwas.sumstats import DenseWriteProgress
from torchgwas.sumstats_indexed import IndexedChunkWrite, IndexedOutputPartition


def test_checkpoint_keeps_issued_gpu_when_future_moves_to_another_gpu(tmp_path):
    header = header_for(tmp_path)
    original = [dict(id='original', device='cuda:0',
                     variant_range=[0, 12], trait_range=[0, 5])]
    run = ProductiveTuningRun(original, chunk_sizes=[2, 4], initial=4)
    progress = ProductiveBoundaryProgress(original, None)
    reserve = run.for_partition('original')
    assert reserve(0, 12, run.capacity) == 4
    assert reserve(4, 12, run.capacity) == 4
    event = DenseWriteProgress(0, 2, (0, 5), 80, 2,
                               time.perf_counter(), 'sumstats', 'cuda:0')
    run.output_written(event)
    progress.observe(event)
    snapshot = run.snapshot()
    frontier = frontier_for(snapshot, header, None)
    boundary = progress.bind(snapshot)
    issued = productive_issued_work(boundary, frontier, header,
                                    covariate_rank=3)
    schedule = header.schedule_bounds(8, 12, 2)
    profile = dict(
        decode_units={name: 1e-8 for name in schedule['source_units']},
        cpu_fraction=.5, depth=2, decode_workers=2,
        cpu_available_cores=2.,
        shared_dram_bytes_per_second=1e8,
        read_bytes_per_second=1e7,
        input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-8,
                                   cpu_seconds_per_call=1e-5))
    candidate = [dict(id='moved', device='cuda:1',
                      variant_range=[8, 12], trait_range=[0, 5])]
    source = productive_source_floor(
        frontier, header, candidate, chunk_markers=2,
        profiles={'moved': profile},
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis='trait')
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component):
        partial = productive_partial_floor(
            frontier, source,
            compute_options=dict(covariate_rank=3,
                shared_h2d_bytes_per_second=1000.,
                per_device_h2d_bytes_per_second={'cuda:1': 800.},
                peak_fp32_flops_per_second={'cuda:1': 100000.}),
            output_options=dict(store_beta=True,
                shared_d2h_bytes_per_second=500.,
                per_device_d2h_bytes_per_second={'cuda:1': 400.},
                output_bytes_per_second=200.),
            output_boundary=boundary,
            shape_profiles={'cuda:1': _profile(129, 3, [5], [2], None)})
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component):
        issued_gpu = productive_issued_gpu_shape_service(
            issued, frontier,
            profiles={'cuda:0': _profile(129, 3, [5], [4], None)})
    ledger = productive_checkpoint_ledger(
        frontier, partial, issued, issued_gpu_shape_service=issued_gpu)
    assert ledger['issued_pending_chunks'] == 2
    assert ledger['unissued_required']['h2d_bytes'] == 4 * 129
    assert ledger['issued_full_chunk_upper']['h2d_bytes'] == 8 * 129
    assert ledger['nominal_subset_work_upper']['h2d_bytes'] == 12 * 129
    assert ledger['nominal_subset_work_upper']['output_array_payload_bytes_upper'] == 400
    assert ledger['per_device_nominal_subset_upper']['cuda:0']['h2d_bytes'] == 8 * 129
    assert ledger['per_device_nominal_subset_upper']['cuda:1']['h2d_bytes'] == 4 * 129
    assert ledger['conditional_unissued_gpu_shape_service']['distinct_shapes'] == 1
    assert ledger['conditional_unissued_gpu_shape_service']['per_device']['cuda:1']['chunks'] == 2
    assert ledger['conditional_issued_gpu_shape_service']['per_device']['cuda:0']['chunks'] == 2
    assert ledger['conditional_issued_gpu_shape_service']['per_device']['cuda:0']['kernel_service_seconds'] == pytest.approx(4.)
    assert 'producer stage state' in ledger['missing_completion_terms'][1]
    assert ledger['missing_completion_terms'] and not ledger['prediction_complete']
    assert ledger['missing_mode_service_floors'] == ['dense_writer_service']
    stale = deepcopy(issued)
    stale['issued_revision'] += 1
    with pytest.raises(ValueError, match='revision'):
        productive_checkpoint_ledger(frontier, partial, stale)
    changed = deepcopy(partial)
    changed['gpu_shape_service']['partitions'][0]['markers'] += 1
    with pytest.raises(ValueError, match='GPU shape partition'):
        productive_checkpoint_ledger(frontier, changed, issued)
    changed_gpu = deepcopy(issued_gpu)
    changed_gpu['written_events'] += 1
    with pytest.raises(ValueError, match='Issued GPU shape'):
        productive_checkpoint_ledger(
            frontier, partial, issued, issued_gpu_shape_service=changed_gpu)
    run.finish(successful=False)


def test_jagwas_checkpoint_uses_one_full_panel_occupancy_scenario(tmp_path):
    header = header_for(tmp_path)
    original = [
        dict(id='left', device='cuda:0', variant_range=[0, 6], trait_range=[0, 5]),
        dict(id='right', device='cuda:1', variant_range=[6, 12], trait_range=[0, 5])]
    run = ProductiveTuningRun(original, chunk_sizes=[3, 6], initial=3)
    progress = ProductiveBoundaryProgress(original, 'jagwas')
    assert run.for_partition('left')(0, 6, run.capacity) == 3
    assert run.for_partition('right')(6, 12, run.capacity) == 3
    now = time.perf_counter()
    event = IndexedChunkWrite(0, 3, 'jagwas', 0, 0, None, now, now,
        False, IndexedOutputPartition('cuda:0', (0, 6), (0, 5)), (0, 3))
    run.output_written(event)
    progress.observe(event)
    snapshot = run.snapshot()
    frontier = frontier_for(snapshot, header, 'jagwas')
    boundary = progress.bind(snapshot)
    issued = productive_issued_work(boundary, frontier, header,
                                    covariate_rank=2)
    future = [
        dict(id='future_left', device='cuda:0',
             variant_range=[3, 6], trait_range=[0, 5]),
        dict(id='future_right', device='cuda:1',
             variant_range=[9, 12], trait_range=[0, 5])]
    schedules = [header.schedule_bounds(*row['variant_range'], 3)
                 for row in future]
    profiles = {row['id']: dict(
        decode_units={name: 1e-8 for name in schedule['source_units']},
        cpu_fraction=.5, depth=2, decode_workers=2,
        cpu_available_cores=2., shared_dram_bytes_per_second=1e8,
        read_bytes_per_second=1e7,
        input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-8,
                                   cpu_seconds_per_call=1e-5))
        for row, schedule in zip(future, schedules)}
    source = productive_source_floor(
        frontier, header, future, chunk_markers=3, profiles=profiles,
        shared_capacities=dict(cpu=1., dram=5e7, input=5e6),
        partition_axis='variant')
    partial = productive_partial_floor(
        frontier, source,
        compute_options=dict(covariate_rank=2,
            shared_h2d_bytes_per_second=1000.,
            per_device_h2d_bytes_per_second={'cuda:0': 800., 'cuda:1': 800.},
            peak_fp32_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.},
            peak_fp64_flops_per_second={'cuda:0': 100000., 'cuda:1': 100000.}),
        output_options=dict(store_beta=True,
            shared_d2h_bytes_per_second=500.,
            per_device_d2h_bytes_per_second={'cuda:0': 400., 'cuda:1': 400.},
            output_bytes_per_second=200.),
        occupancy_scenario='dense', output_boundary=boundary)
    ledger = productive_checkpoint_ledger(frontier, partial, issued)
    assert ledger['output_scenario'] == 'dense'
    assert ledger['issued_pending_chunks'] == 1
    assert ledger['nominal_subset_work_upper']['fp64_projection_flops'] == 450
    assert ledger['nominal_subset_work_upper']['output_array_payload_bytes_upper'] == 144
    assert ledger['missing_mode_service_floors'] == [
        'jagwas_selection_service', 'jagwas_archive_service']
    result_queue = queue.Queue(maxsize=2)
    finished = object()
    result_queue.put((6, 9, None, np.arange(15, dtype=np.float64), None, None))
    queue_observation = snapshot_indexed_result_queue(
        result_queue, [(0, 6), (6, 12)], ['cuda:0', 'cuda:1'], finished)
    joined_queue = bind_jagwas_result_queue(boundary, queue_observation)
    queue_ledger = productive_checkpoint_ledger(frontier, partial, issued,
        bound_jagwas_result_queue=joined_queue)
    assert queue_ledger['issued_remaining_producer_chunks'] == 0
    assert queue_ledger['issued_queue_refinement']['queued_chunks_with_completed_producer'] == 1
    assert queue_ledger['nominal_subset_work_upper']['indexed_read_bytes'] == (
        ledger['unissued_required']['indexed_read_bytes'])
    assert queue_ledger['nominal_subset_work_upper']['h2d_bytes'] == (
        ledger['unissued_required']['h2d_bytes'])
    assert queue_ledger['nominal_subset_work_upper']['fp64_projection_flops'] == 300
    assert queue_ledger['issued_output_array_payload_upper_bytes'] == (
        ledger['issued_output_array_payload_upper_bytes'])
    result_queue.get_nowait()
    empty_observation = snapshot_indexed_result_queue(
        result_queue, [(0, 6), (6, 12)], ['cuda:0', 'cuda:1'], finished)
    active_join = bind_jagwas_result_queue(boundary, empty_observation,
        active_writer=dict(device='cuda:1', variant_range=[6, 9],
                           trait_range=[0, 5]))
    active_ledger = productive_checkpoint_ledger(frontier, partial, issued,
        bound_jagwas_result_queue=active_join)
    assert active_ledger['issued_queue_refinement'][
        'active_writer_chunks_with_completed_producer'] == 1
    assert active_ledger['nominal_subset_work_upper'] == (
        queue_ledger['nominal_subset_work_upper'])
    with patch('torchgwas.mechanistic_torch._shape_component',
               side_effect=_component):
        issued_gpu = productive_issued_gpu_shape_service(
            issued, frontier,
            profiles={'cuda:1': _profile(129, 2, [5], [3], 'jagwas')})
    assert issued_gpu['pending_chunks'] == 1
    assert issued_gpu['per_device_service']['cuda:1']['kernel_service_seconds'] == pytest.approx(1.5)
    with pytest.raises(ValueError, match='without full issued GPU shape service'):
        productive_checkpoint_ledger(frontier, partial, issued,
            issued_gpu_shape_service=issued_gpu,
            bound_jagwas_result_queue=joined_queue)
    changed = deepcopy(partial)
    changed['issued_output_backlog']['occupancy_scenario'] = 'empty'
    with pytest.raises(ValueError, match='scenarios differ'):
        productive_checkpoint_ledger(frontier, changed, issued)
    run.finish(successful=False)
