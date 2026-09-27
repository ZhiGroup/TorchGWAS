"""Bounded post-output layout screen using one completed PGEN source stage.

This joins exact unissued coverage, independently priced component work and
mode-specific output service. Its floors cannot authorize a JIT switch: queued
work, live availability, final publication and a completion ceiling remain.
"""
import math
import time

from .productive_source_floor import productive_partial_floor, productive_source_floor
from .productive_source_stage import ProductiveSourceStage


def _mode_services(source, output, options):
    """Build only mode-matching, source-bound service reports."""
    if not isinstance(options, dict):
        raise ValueError('Explicit output-mode service prices required')
    reduction = source['reduction']
    if reduction is None:
        from .layout_dense_writer_service import native_layout_dense_writer_service_floor
        if set(options) != {'writer_profiles'}:
            raise ValueError('Dense writer service profiles required')
        return dict(writer_service=native_layout_dense_writer_service_floor(
            source, output, options['writer_profiles']))
    if reduction == 'jagwas':
        from .layout_jagwas_selection_floor import native_layout_jagwas_selection_floor
        from .layout_jagwas_archive_floor import native_layout_jagwas_archive_floor
        if set(options) != {'selection_prices', 'cpu_fraction_by_device',
                            'archive_price', 'archive_profiles'}:
            raise ValueError('Complete JAGWAS selection/archive prices required')
        return dict(
            jagwas_selection=native_layout_jagwas_selection_floor(
                source, output, options['selection_prices'],
                options['cpu_fraction_by_device']),
            jagwas_archive=native_layout_jagwas_archive_floor(
                source, output, options['archive_price'],
                options['archive_profiles']))
    if output['significant_backend'] == 'host':
        from .layout_significant_host_selection_floor import (
            native_layout_significant_host_selection_floor)
        from .layout_significant_archive_floor import native_layout_significant_archive_floor
        if set(options) != {'selection_prices', 'selection_profiles',
                            'archive_prices', 'archive_profiles'}:
            raise ValueError('Complete host significant selection/archive prices required')
        return dict(
            significant_host_selection=native_layout_significant_host_selection_floor(
                source, output, options['selection_prices'],
                options['selection_profiles']),
            significant_archive=native_layout_significant_archive_floor(
                source, output, options['archive_prices'],
                options['archive_profiles']))
    from .layout_significant_device_count_floor import native_layout_significant_device_count_floor
    from .layout_significant_device_launch_floor import native_layout_significant_device_launch_floor
    from .layout_significant_archive_floor import native_layout_significant_archive_floor
    if set(options) != {'count_transfer_prices', 'launch_profiles',
                        'archive_prices', 'archive_profiles'}:
        raise ValueError('Complete device significant count/archive prices required')
    return dict(
        significant_device_count=native_layout_significant_device_count_floor(
            source, output, options['count_transfer_prices']),
        significant_device_launch=native_layout_significant_device_launch_floor(
            source, output, options['launch_profiles']),
        significant_archive=native_layout_significant_archive_floor(
            source, output, options['archive_prices'],
            options['archive_profiles']))


def productive_staged_partial_screen(frontier, stage, candidates, *,
                                      shared_source_capacities,
                                      occupancy_scenario=None,
                                      output_boundary=None,
                                      max_candidates=4,
                                      max_partitions=16,
                                      max_unique_records=20_000_000,
                                      max_chunks_per_partition=100_000,
                                      max_cpu_seconds=10.,
                                      max_wall_seconds=30.):
    """Screen explicit admitted layouts after useful output, with no switch.

    Each candidate declares id, chunk_markers, partitions, partition_axis,
    source_profiles, compute_options, output_options, mode_service_options,
    and optional shape_profiles. The same source/output scenario and shared
    source capacities apply to all candidates. The caller must still prove
    current price bindings, memory admission and a live transition before a
    future comparison. A cooperative CPU/wall budget stops after a candidate;
    it cannot preempt one in-progress metadata rebase.
    """
    for name, value in (('max_candidates', max_candidates),
                        ('max_partitions', max_partitions),
                        ('max_unique_records', max_unique_records),
                        ('max_chunks_per_partition', max_chunks_per_partition)):
        if type(value) is not int or value < 1:
            raise ValueError('Positive bounded ' + name + ' required')
    for name, value in (('max_cpu_seconds', max_cpu_seconds),
                        ('max_wall_seconds', max_wall_seconds)):
        if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or value <= 0):
            raise ValueError('Positive finite screen ' + name + ' required')
    if (not isinstance(frontier, dict) or
            not isinstance(stage, ProductiveSourceStage) or
            not isinstance(candidates, list) or
            not 0 < len(candidates) <= max_candidates or
            not isinstance(shared_source_capacities, dict) or
            set(shared_source_capacities) != {'cpu', 'dram', 'input'}):
        raise ValueError('Completed productive stage and bounded candidates required')
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or
           not math.isfinite(value) or value <= 0
           for value in shared_source_capacities.values()):
        raise ValueError('Positive finite shared source capacities required')
    ledger, prior = stage.completed_ledger()
    header = ledger.header
    if (frontier.get('input_identity') != prior['input_identity'] or
            frontier.get('input_identity') != header.input_identity):
        raise ValueError('Staged source differs from the held issue frontier')
    ids = [row.get('id') for row in candidates if isinstance(row, dict)]
    if (len(ids) != len(candidates) or
            any(not isinstance(key, str) or not key for key in ids) or
            len(set(ids)) != len(ids)):
        raise ValueError('Unique named screen candidates required')
    required = {'id', 'chunk_markers', 'partitions', 'partition_axis',
                'source_profiles', 'compute_options', 'output_options',
                'mode_service_options'}
    # Validate schema before touching any multi-million-record source schedule.
    if any(set(row) not in (required, required | {'shape_profiles'})
           for row in candidates):
        raise ValueError('Explicit priced screen candidate required')
    started = time.perf_counter()
    cpu_started = time.thread_time()
    reports = []
    stop_reason = 'complete'
    for row in candidates:
        if (reports and (time.thread_time() - cpu_started >= max_cpu_seconds or
                         time.perf_counter() - started >= max_wall_seconds)):
            stop_reason = 'screen_budget'
            break
        source_report = productive_source_floor(
            frontier, header, row['partitions'],
            chunk_markers=row['chunk_markers'],
            profiles=row['source_profiles'],
            shared_capacities=shared_source_capacities,
            partition_axis=row['partition_axis'],
            max_partitions=max_partitions,
            max_unique_records=max_unique_records,
            max_chunks_per_partition=max_chunks_per_partition,
            staged_source=stage)
        common = dict(compute_options=row['compute_options'],
                      output_options=row['output_options'],
                      occupancy_scenario=occupancy_scenario,
                      output_boundary=output_boundary,
                      shape_profiles=row.get('shape_profiles'))
        partial = productive_partial_floor(frontier, source_report,
            stage_service_factory=lambda source, output: _mode_services(
                source, output, row['mode_service_options']), **common)
        if partial['envelope']['coverage']['exact_coverage'] is not True:
            raise ValueError('Candidate changes the held association frontier')
        reports.append(dict(id=row['id'], chunk_markers=row['chunk_markers'],
            source_schedule_method=source_report['source_schedule_method'],
            unique_source_schedules=source_report['unique_schedules'],
            unique_priced_source_floors=source_report['unique_priced_source_floors'],
            source_floor_cpu_seconds=source_report['calculation_cpu_seconds'],
            source_floor_wall_seconds=source_report['calculation_wall_seconds'],
            partial=partial))
    return dict(kind='torchgwas.productive_staged_partial_screen.v1',
        input_identity=prior['input_identity'],
        issued_revision=frontier['issued_revision'],
        written_events=frontier['written_events'],
        reduction=frontier['reduction'],
        occupancy_scenario=occupancy_scenario,
        prior_stage_cost_once=prior,
        candidates=reports,
        requested_candidates=len(candidates),
        evaluated_candidates=len(reports),
        stop_reason=stop_reason,
        screen_cpu_seconds=time.thread_time() - cpu_started,
        screen_wall_seconds=time.perf_counter() - started,
        prediction_complete=False, selection_validated=False,
        scope='Exact post-output unissued coverage and conditional source, compute, output, and mode-specific service floors from one completed staged ledger. Prior stage work is charged once. Budget is cooperative after each candidate; partial floors are not elapsed completion bounds, and no layout switch is authorized.')
