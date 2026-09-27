"""Join fixed issued work and candidate unissued work at one JIT checkpoint.

The joined counts are deliberately narrower than a completion prediction.
Issued chunks may already have finished source/GPU stages, while unissued
chunks must still run. Keeping those scopes separate prevents an in-flight
upper workload from being mistaken for a necessary resource floor.
"""

from copy import deepcopy
import math

from .layout_frontier import KIND as FRONTIER_KIND


def _count(value, name):
    if (type(value) is float and math.isfinite(value) and
            0 <= value <= 2**53 and value.is_integer()):
        return int(value)
    if type(value) is not int or value < 0:
        raise ValueError('Nonnegative integer ' + name + ' required')
    return value


def productive_checkpoint_ledger(frontier, partial, issued, *,
                                 issued_gpu_shape_service=None, bound_jagwas_result_queue=None):
    """Compose disjoint work ledgers without converting any to elapsed time."""
    if (not isinstance(frontier, dict) or
            frontier.get('kind') != FRONTIER_KIND or
            not isinstance(partial, dict) or
            partial.get('kind') != 'torchgwas.productive_partial_floor.v1' or
            not isinstance(issued, dict) or
            issued.get('kind') != 'torchgwas.productive_issued_work.v1' or
            not isinstance(partial.get('envelope'), dict) or
            partial['envelope'].get('coverage', {}).get('exact_coverage') is not True or
            not isinstance(partial.get('issued_output_backlog'), dict)):
        raise ValueError('Exact priced unissued and bound issued work required')
    source, compute, output = (partial.get(name) for name in
                               ('source', 'compute', 'output'))
    backlog = partial['issued_output_backlog']
    if any(not isinstance(row, dict) for row in (source, compute, output)):
        raise ValueError('Typed source, GPU and output work required')
    for row in (source, compute, output, issued):
        if (row.get('input_identity') != frontier.get('input_identity') or
                row.get('reduction') != frontier.get('reduction')):
            raise ValueError('Checkpoint source/output mode differs')
    for row in (partial['envelope'], issued, backlog):
        if (row.get('issued_revision') != frontier.get('issued_revision') or
                row.get('written_events') != frontier.get('written_events')):
            raise ValueError('Checkpoint issue/output revision differs')
    if (backlog.get('kind') != 'torchgwas.productive_output_backlog.v1' or
            backlog.get('reduction') != frontier['reduction'] or
            backlog.get('store_beta') != output.get('store_beta') or
            compute.get('samples') != issued.get('samples') or
            compute.get('covariate_rank') != issued.get('covariate_rank') or
            partial['envelope'].get('required_cells') !=
            frontier.get('required_cells')):
        raise ValueError('Checkpoint workload contract differs')
    survivors = partial.get('survivors')
    expected_scenario = None if survivors is None else survivors.get('scenario')
    if backlog.get('occupancy_scenario') != expected_scenario:
        raise ValueError('Issued and unissued output scenarios differ')
    if (not isinstance(issued.get('original_partitions'), dict) or
            not isinstance(issued.get('chunks'), list) or
            issued.get('pending_chunks') != len(issued['chunks'])):
        raise ValueError('Bound original issued producer ledger required')
    for rectangle in frontier['rectangles']:
        original = issued['original_partitions'].get(rectangle['id'])
        if (not isinstance(original, dict) or
                original.get('device') != rectangle['device'] or
                original.get('trait_range') != rectangle['trait_range'] or
                [original.get('issued_to'),
                 original.get('variant_range', [None, None])[1]] !=
                rectangle['variant_range']):
            raise ValueError('Issued producer differs from unissued rectangle')
    fields = ('read_bytes', 'h2d_bytes', 'fp32_gemm_flops',
              'fp64_projection_flops')
    issued_totals = issued.get('total_work')
    if not isinstance(issued_totals, dict):
        raise ValueError('Issued full-chunk workload required')
    for name in fields:
        amount = _count(issued_totals.get(name), name)
        if amount != sum(_count(chunk.get(name), name)
                         for chunk in issued['chunks']):
            raise ValueError('Issued per-chunk work differs from total')
    queue_refinement = None
    counted_issued = issued_totals
    if bound_jagwas_result_queue is not None:
        if frontier['reduction'] != 'jagwas' or issued_gpu_shape_service is not None:
            raise ValueError('Queue refinement requires JAGWAS without full issued GPU shape service')
        from .productive_issued_queue_refinement import refine_issued_jagwas_with_queue
        queue_refinement = refine_issued_jagwas_with_queue(
            issued, bound_jagwas_result_queue)
        counted_issued = queue_refinement['upstream_or_unresolved_full_chunk_work_upper']
    future_read = _count(source['resource_work'].get('input_bytes'),
                         'unissued indexed read bytes')
    future_h2d = _count(compute.get('total_h2d_bytes'),
                        'unissued H2D bytes')
    future_output = output.get('total_output_array_payload_bytes')
    if (not isinstance(future_output, list) or len(future_output) != 2 or
            _count(future_output[0], 'unissued output bytes') >
            _count(future_output[1], 'unissued output bytes')):
        raise ValueError('Bound unissued output array payload required')
    pending_output = _count(backlog.get('array_payload_bytes_upper'),
                            'issued output array payload')
    if partial.get('array_payload_work_upper_bytes') != future_output[1] + pending_output:
        raise ValueError('Issued/unissued output upper work differs')
    future_device = compute.get('per_device_work')
    pending_device = (issued.get('per_device_work') if queue_refinement is None else
        queue_refinement['per_device_upstream_or_unresolved_upper'])
    if not isinstance(future_device, dict) or not isinstance(pending_device, dict):
        raise ValueError('Per-device GPU work required')
    devices = set(future_device) | set(pending_device)
    per_device = {}
    for device in sorted(devices):
        future = future_device.get(device, {})
        pending = pending_device.get(device, {})
        per_device[device] = {}
        for name in ('h2d_bytes', 'fp32_gemm_flops',
                     'fp64_projection_flops'):
            per_device[device][name] = (
                _count(future.get(name, 0), device + ' unissued ' + name) +
                _count(pending.get(name, 0), device + ' issued ' + name))
    if (sum(row['h2d_bytes'] for row in per_device.values()) !=
            future_h2d + counted_issued['h2d_bytes']):
        raise ValueError('Shared and per-device H2D work differ')
    total = dict(
        indexed_read_bytes=future_read + counted_issued['read_bytes'],
        h2d_bytes=future_h2d + counted_issued['h2d_bytes'],
        fp32_gemm_flops=sum(row['fp32_gemm_flops'] for row in per_device.values()),
        fp64_projection_flops=sum(row['fp64_projection_flops']
                                  for row in per_device.values()),
        output_array_payload_bytes_upper=future_output[1] + pending_output)
    stage_floors = partial['envelope'].get('stage_floor_seconds')
    if not isinstance(stage_floors, dict):
        raise ValueError('Explicit incomplete stage-floor audit required')
    mode = frontier['reduction']
    if mode is None:
        required_floors = ('dense_writer_service',)
    elif mode == 'jagwas':
        required_floors = ('jagwas_selection_service',
                           'jagwas_archive_service')
    elif output.get('significant_backend') == 'host':
        required_floors = ('significant_host_selector_service',
                           'significant_archive_service')
    else:
        required_floors = ('significant_device_count_barrier',
                           'significant_archive_service')
    missing_floors = [name for name in required_floors
                      if stage_floors.get(name) is None]
    gpu_shape = partial.get('gpu_shape_service')
    if gpu_shape is not None:
        if (not isinstance(gpu_shape, dict) or
                gpu_shape.get('kind') !=
                'torchgwas.pgen_layout_gpu_shape_service.v1' or
                gpu_shape.get('input_identity') != source['input_identity'] or
                gpu_shape.get('reduction') != source['reduction'] or
                gpu_shape.get('samples') != compute['samples'] or
                gpu_shape.get('covariate_rank') != compute['covariate_rank'] or
                gpu_shape.get('chunk_markers') != source['chunk_markers'] or
                not isinstance(gpu_shape.get('partitions'), list) or
                len(gpu_shape['partitions']) != len(source['partitions'])):
            raise ValueError('GPU shape service differs from checkpoint')
        for bound, priced in zip(source['partitions'],
                                 gpu_shape['partitions']):
            if (priced.get('id') != bound['id'] or
                    priced.get('device') != bound['device'] or
                    priced.get('markers') !=
                    bound['variant_range'][1] - bound['variant_range'][0] or
                    priced.get('traits') !=
                    bound['trait_range'][1] - bound['trait_range'][0] or
                    priced.get('chunks') != bound['chunks']):
                raise ValueError('GPU shape partition differs from checkpoint')
        if (set(gpu_shape.get('per_device_service', {})) !=
                {row['device'] for row in source['partitions']}):
            raise ValueError('GPU shape device set differs from checkpoint')
    if issued_gpu_shape_service is not None:
        row = issued_gpu_shape_service
        if (not isinstance(row, dict) or
                row.get('kind') !=
                'torchgwas.productive_issued_gpu_shape_service.v1' or
                any(row.get(key) != issued.get(key) for key in
                    ('input_identity', 'issued_revision', 'written_events',
                     'reduction', 'samples', 'covariate_rank',
                     'pending_chunks', 'pending_source_records')) or
                set(row.get('per_device_service', {})) !=
                set(issued['per_device_work'])):
            raise ValueError('Issued GPU shape service differs from checkpoint')
    return dict(kind='torchgwas.productive_checkpoint_ledger.v1',
                input_identity=deepcopy(frontier['input_identity']),
                issued_revision=frontier['issued_revision'],
                written_events=frontier['written_events'],
                reduction=frontier['reduction'],
                output_scenario=deepcopy(expected_scenario),
                required_unissued_cells=frontier['required_cells'],
                issued_pending_chunks=issued['pending_chunks'],
                issued_remaining_producer_chunks=(issued['pending_chunks'] if
                    queue_refinement is None else
                    queue_refinement['upstream_or_unresolved_chunks']),
                unissued_required=dict(indexed_read_bytes=future_read,
                    h2d_bytes=future_h2d,
                    output_array_payload_bytes=list(future_output)),
                issued_full_chunk_upper=dict(read_bytes=issued_totals['read_bytes'],
                    h2d_bytes=issued_totals['h2d_bytes'],
                    fp32_gemm_flops=issued_totals['fp32_gemm_flops'],
                    fp64_projection_flops=issued_totals['fp64_projection_flops']),
                issued_output_array_payload_upper_bytes=pending_output,
                issued_queue_refinement=deepcopy(queue_refinement),
                nominal_subset_work_upper=total,
                per_device_nominal_subset_upper=per_device,
                available_stage_floor_seconds=deepcopy(stage_floors),
                conditional_unissued_gpu_shape_service=(
                    None if gpu_shape is None else dict(
                        distinct_shapes=gpu_shape['distinct_shapes'],
                        per_device=deepcopy(gpu_shape['per_device_service']),
                        total_host_work=deepcopy(gpu_shape['total_host_work']),
                        profile_sha256=deepcopy(gpu_shape['profile_sha256']))),
                conditional_issued_gpu_shape_service=(
                    None if issued_gpu_shape_service is None else dict(
                        distinct_shapes=issued_gpu_shape_service['distinct_shapes'],
                        per_device=deepcopy(
                            issued_gpu_shape_service['per_device_service']),
                        total_host_work=deepcopy(
                            issued_gpu_shape_service['total_host_work']),
                        profile_sha256=deepcopy(
                            issued_gpu_shape_service['profile_sha256']))),
                missing_mode_service_floors=missing_floors,
                missing_completion_terms=[
                    'issued decoder and read service still in flight' if queue_refinement is None else
                    'upstream or active-writer issued source stage is uncertain',
                    ('complete shape-specific GPU statistics/reduction service'
                     if gpu_shape is None else
                     'issued GPU shape service and producer stage state'
                     if issued_gpu_shape_service is None else
                     'producer stage state and unpriced GPU operations'),
                    'issued D2H and selection service',
                    'live result queues and writer staging' if queue_refinement is None else
                    'queued JAGWAS results and active writer still need selection and durable output',
                    'writeback, fsync and manifest/directory publication',
                    'loaded shared-capacity and scheduler service'],
                prediction_complete=False, selection_validated=False,
                scope='Exact unissued coverage plus conservative issued producer workload, optionally refined by a bound JAGWAS result queue, and conditional output-array workload at one source/output revision. Totals are upper workloads for a subset of necessary operations, not a lower resource floor or elapsed completion ceiling. Common issued work remains on original devices while candidate unissued work may move. No JIT switch authorization.')
