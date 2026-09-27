"""Bound compact whole-source work to a written productive issue frontier.

This is the first post-output pass for a finite JIT continuation. It can
compare chunk sizes and alternative trait/variant assignments, but it only
reports necessary read/decode work and exact unissued coverage. It neither
prices the rest of the pipeline nor authorizes a layout change.
"""
import time

from .calibration_cache import _digest
from .layout_frontier import KIND, _rectangles_match, bind_layout_to_frontier
from .pgen_work_bounds import PgenHeaderWork, native_schedule_source_floor
from .source_layout_floor import native_layout_source_floor


def productive_source_floor(frontier, header, partitions, *, chunk_markers,
                            profiles, shared_capacities, partition_axis,
                            max_partitions=16, max_unique_records=20_000_000,
                            max_chunks_per_partition=100_000,
                            staged_source=None):
    """Price exact unissued PGEN source schedules with bounded header work.

    Each candidate partition has an id, device, variant_range and trait_range.
    Profiles are independently measured component prices keyed by partition
    id. Identical source schedules are inspected once even when multiple trait
    tiles reread the same records; their physical work is still charged once
    per tile by the layout floor.
    """
    for name, value in (('max_partitions', max_partitions),
                        ('max_unique_records', max_unique_records),
                        ('max_chunks_per_partition', max_chunks_per_partition),
                        ('chunk_markers', chunk_markers)):
        if type(value) is not int or value < 1:
            raise ValueError('Positive bounded ' + name + ' required')
    if (not isinstance(header, PgenHeaderWork) or
            not isinstance(frontier, dict) or
            frontier.get('kind') != KIND or
            frontier.get('input_identity') != header.input_identity):
        raise ValueError('Source header differs from the productive frontier')
    staged_worker_cost = None
    if staged_source is not None:
        from .incremental_pgen_schedule import IncrementalPgenSchedule
        from .productive_source_stage import ProductiveSourceStage
        if isinstance(staged_source, ProductiveSourceStage):
            staged_source, staged_worker_cost = staged_source.completed_ledger()
            if staged_worker_cost['input_identity'] != header.input_identity:
                raise ValueError('Productive source stage differs from this PGEN')
        if (not isinstance(staged_source, IncrementalPgenSchedule) or
                staged_source.input_identity != header.input_identity or
                not staged_source.snapshot()['complete']):
            raise ValueError('Complete staged source for this PGEN required')
        staged_cost = staged_source.snapshot()
    else:
        staged_cost = None
    if (not isinstance(partitions, list) or
            not 0 < len(partitions) <= max_partitions or
            not isinstance(profiles, dict) or
            not isinstance(shared_capacities, dict) or
            set(shared_capacities) != {'cpu', 'dram', 'input'}):
        raise ValueError('Bounded source layout and capacities required')
    ids = [row.get('id') for row in partitions if isinstance(row, dict)]
    if (len(ids) != len(partitions) or
            any(not isinstance(key, str) or not key for key in ids) or
            len(set(ids)) != len(ids) or
            set(ids) != set(profiles)):
        raise ValueError('Exactly one priced profile per source partition required')
    if (partition_axis not in ('trait', 'variant') or
            (frontier['reduction'] == 'jagwas' and partition_axis != 'variant') or
            (frontier['reduction'] == 'significant' and partition_axis != 'trait')):
        raise ValueError('Source partition axis differs from output mode')

    # Reject an oversized plan before parsing any whole-source schedule.
    unique = set()
    for row in partitions:
        if set(row) != {'id', 'device', 'variant_range', 'trait_range'}:
            raise ValueError('Explicit candidate source partition required')
        span = row['variant_range']
        if (not isinstance(span, (list, tuple)) or len(span) != 2 or
                any(type(value) is not int for value in span) or
                not 0 <= span[0] < span[1] <= header._header.variant_ct):
            raise ValueError('Candidate source range differs from PGEN header')
        trait = row['trait_range']
        if (not isinstance(trait, (list, tuple)) or len(trait) != 2 or
                any(type(value) is not int for value in trait) or
                not 0 <= trait[0] < trait[1] <= frontier['total_traits']):
            raise ValueError('Candidate phenotype range differs from frontier')
        if (frontier['reduction'] == 'jagwas' and
                list(trait) != [0, frontier['total_traits']]):
            raise ValueError('JAGWAS requires the complete phenotype panel')
        if (span[1] - span[0] + chunk_markers - 1) // chunk_markers > max_chunks_per_partition:
            raise ValueError('Whole-source chunk budget exceeded')
        unique.add(tuple(span))
    if sum(hi - lo for lo, hi in unique) > max_unique_records:
        raise ValueError('Whole-source record budget exceeded')
    # Coverage is a cheap rectangle sweep. Reject a stale/overlapping tile or
    # GPU proposal before a multi-million-record schedule walk.
    _rectangles_match(frontier['rectangles'], partitions, 1_000_000)

    started = time.perf_counter()
    cpu_started = time.thread_time()
    schedules = {}
    priced_floors = {}
    priced = []
    for row in partitions:
        span = tuple(row['variant_range'])
        if span not in schedules:
            schedules[span] = (header.schedule_bounds(
                *span, chunk_markers,
                max_records=max_unique_records,
                max_chunks=max_chunks_per_partition)
                if staged_source is None else
                staged_source.rebase(span[0], chunk_markers, stop=span[1]))
        # Trait tiles reread the same physical source, so their work is
        # charged once per tile below. The immutable source price and exact
        # schedule need only be calculated once for identical profiles.
        price_key = (span, _digest(profiles[row['id']]))
        if price_key not in priced_floors:
            priced_floors[price_key] = native_schedule_source_floor(
                schedules[span], profiles[row['id']], shared_capacities)
        priced.append(dict(id=row['id'], device=row['device'],
                           trait_range=list(row['trait_range']),
                           floor=priced_floors[price_key]))
    layout = native_layout_source_floor(
        priced, total_traits=frontier['total_traits'],
        reduction=frontier['reduction'], partition_axis=partition_axis,
        max_partitions=max_partitions)
    coverage = bind_layout_to_frontier(frontier, layout)
    return dict(source=layout, coverage=coverage,
                source_schedule_method=('direct_header' if staged_source is None
                                        else 'staged_primary_rebase'),
                staged_segments_used=(None if staged_source is None else
                                      staged_source.segments),
                prior_staged_calculation=(None if staged_cost is None else
                    dict(cpu_seconds=(staged_cost['calculation_cpu_seconds']
                                      if staged_worker_cost is None else
                                      staged_worker_cost['cpu_seconds']),
                         wall_seconds=(staged_cost['calculation_wall_seconds']
                                       if staged_worker_cost is None else
                                       staged_worker_cost['wall_seconds']),
                         scope=('Incremental schedule calculation only; a public worker must additionally charge header and observation work.'
                                if staged_worker_cost is None else
                                staged_worker_cost['scope']))),
                productive_stage_worker_cost=staged_worker_cost,
                unique_schedules=len(schedules),
                unique_priced_source_floors=len(priced_floors),
                unique_source_records=sum(hi - lo for lo, hi in unique),
                calculation_wall_seconds=time.perf_counter() - started,
                calculation_cpu_seconds=time.thread_time() - cpu_started,
                prediction_complete=False, selection_validated=False,
                scope='Exact unissued coverage and conditional whole-source read/decode floor after written output. Optional staged primary work is reported separately from this rebase call so all actual planning cost can be charged once per job. No GPU, transfer, output, queue, drain or live transition price; no switch authorization.')


def productive_partial_floor(frontier, source_report, *, compute_options,
                              output_options, occupancy_scenario=None,
                              stage_service_reports=None, stage_service_factory=None,
                              output_boundary=None, shape_profiles=None):
    """Compose necessary whole-job resource loads for one bound JIT scenario.

    This feeds the existing analytical components rather than extrapolating
    short windows. All capacities and optional service reports must be
    independently supplied and current. The returned high endpoint bounds
    only the partial *floor*, never candidate completion.
    """
    from .layout_compute_floor import native_layout_compute_floor
    from .layout_output_floor import native_layout_output_floor
    from .layout_partial_envelope import native_layout_partial_envelope
    from .layout_gpu_shape_service import native_layout_gpu_shape_service
    from .productive_occupancy import whole_layout_survivors
    from .productive_output_backlog import productive_output_backlog

    if (not isinstance(source_report, dict) or
            not isinstance(source_report.get('source'), dict) or
            not isinstance(source_report.get('coverage'), dict) or
            source_report['coverage'].get('exact_coverage') is not True or
            not isinstance(compute_options, dict) or
            not isinstance(output_options, dict) or
            'retained_ranges' in output_options):
        raise ValueError('Bound source and explicit compute/output inputs required')
    if stage_service_factory is not None and (stage_service_reports is not None or
                                               not callable(stage_service_factory)):
        raise ValueError('One bound stage-service source required')
    service = {} if stage_service_reports is None else stage_service_reports
    allowed = {'writer_service', 'jagwas_selection', 'jagwas_archive',
               'significant_archive', 'significant_host_selection',
               'significant_device_count', 'significant_device_launch'}
    if not isinstance(service, dict) or set(service) - allowed:
        raise ValueError('Known bound optional stage services required')
    source = source_report['source']
    coverage = source_report['coverage']
    if (not isinstance(frontier, dict) or
            any(coverage.get(name) != frontier.get(name)
                for name in ('input_identity', 'issued_revision',
                             'written_events', 'required_cells')) or
            source.get('input_identity') != frontier.get('input_identity')):
        raise ValueError('Source floor differs from the held productive frontier')
    reduction = source.get('reduction')
    if (reduction is None) != (occupancy_scenario is None):
        raise ValueError('Reduced output requires one explicit occupancy scenario')
    started = time.perf_counter()
    cpu_started = time.thread_time()
    survivors = (None if reduction is None else
                 whole_layout_survivors(source, occupancy_scenario))
    compute = native_layout_compute_floor(source, **compute_options)
    output = native_layout_output_floor(
        source, retained_ranges=(None if survivors is None else
                                  survivors['retained_ranges']),
        **output_options)
    if stage_service_factory is not None:
        service = stage_service_factory(source, output)
        if not isinstance(service, dict) or set(service) - allowed:
            raise ValueError('Known bound generated stage services required')
    if reduction == 'significant' and shape_profiles is not None:
        expected_scan_mode = (None if output['significant_backend'] == 'host'
                              else 'device_significant')
        if (not isinstance(shape_profiles, dict) or
                any(not isinstance(profile, dict) or
                    profile.get('reduction') != expected_scan_mode
                    for profile in shape_profiles.values())):
            raise ValueError('Significant selector backend differs from GPU shape scan mode')
    gpu_shape = (None if shape_profiles is None else
                 native_layout_gpu_shape_service(
                     source, covariate_rank=compute['covariate_rank'],
                     profiles=shape_profiles))
    backlog = None
    if output_boundary is not None:
        if (not isinstance(output_boundary, dict) or
                output_boundary.get('issued_revision') != frontier['issued_revision'] or
                output_boundary.get('written_events') != frontier['written_events'] or
                output_boundary.get('reduction') != reduction):
            raise ValueError('Issued output boundary differs from the held productive frontier')
        dense_writer = output_options.get('dense_writer_options')
        backlog = productive_output_backlog(
            output_boundary, total_traits=frontier['total_traits'],
            store_beta=output_options.get('store_beta', True),
            store_variant_df=(False if dense_writer is None else
                              dense_writer['store_variant_df']),
            occupancy_scenario=occupancy_scenario)
    envelope = native_layout_partial_envelope(
        frontier, source, compute, output, **service)
    return dict(kind='torchgwas.productive_partial_floor.v1',
                source=source, compute=compute, output=output,
                gpu_shape_service=gpu_shape,
                included_stage_services=sorted(service),
                survivors=survivors, envelope=envelope,
                issued_output_backlog=backlog,
                array_payload_work_upper_bytes=(None if backlog is None else
                    output['total_output_array_payload_bytes'][1] +
                    backlog['array_payload_bytes_upper']),
                calculation_wall_seconds=time.perf_counter() - started,
                calculation_cpu_seconds=time.thread_time() - cpu_started,
                prediction_complete=False, selection_validated=False,
                scope='Exact unissued coverage and conditional necessary whole-source, H2D/GEMM and output loads, with optional independently priced exact GPU full/tail shape service. An optional synchronized issued-output backlog adds only an upper array-payload workload, not an elapsed-time ceiling. Unpriced stage work, queued work, live capacity and final drain prevent completion bounds or a JIT switch.')
