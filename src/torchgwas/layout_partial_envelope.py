"""Bind exact unissued coverage to compact, incomplete resource floors."""
import math

from .calibration_cache import _digest
from .layout_frontier import bind_layout_to_frontier


def _bounds(value, name):
    if (not isinstance(value, list) or len(value) != 2 or
            any(isinstance(x, bool) or not isinstance(x, (int, float)) or
                not math.isfinite(x) or x < 0 for x in value) or value[0] > value[1]):
        raise ValueError('Finite nonnegative ' + name + ' floor interval required')
    return value


def native_layout_partial_envelope(frontier, source, compute, output,
                                   *, writer_service=None,
                                   jagwas_selection=None,
                                   jagwas_archive=None,
                                   significant_archive=None,
                                   significant_host_selection=None,
                                   significant_device_count=None,
                                   significant_device_launch=None,
                                   max_coverage_visits=1000000):
    """Compose only known necessary stage loads for the same unissued pairs.

    The lower endpoint can later be compared with a qualified, feasible
    baseline ceiling to reject an impossible candidate. The upper endpoint
    bounds the possible *partial floor*, not the completion time. Neither
    endpoint authorizes a chunk, tile or GPU switch.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            not isinstance(compute, dict) or
            compute.get('kind') != 'torchgwas.pgen_layout_compute_floor.v1' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1'):
        raise ValueError('Typed source, compute and output floors required')
    coverage = bind_layout_to_frontier(frontier, source,
                                        max_coverage_visits=max_coverage_visits)
    for name, floor in [('compute', compute), ('output', output)]:
        for field in ('input_identity', 'reduction'):
            if floor.get(field) != source.get(field):
                raise ValueError(name + ' floor differs from source ' + field)
    if compute.get('samples') != source.get('samples'):
        raise ValueError('Compute sample count differs from bound source')
    h2d_links = compute.get('shared_h2d_link_loads')
    d2h_links = output.get('shared_d2h_link_loads')
    if (not isinstance(h2d_links, dict) or not isinstance(d2h_links, list) or
            len(d2h_links) != 2 or any(not isinstance(row, dict) or
                row.get('declarations') != h2d_links.get('declarations')
                for row in d2h_links)):
        raise ValueError('Compute/output shared-link declarations differ')
    rows = source['partitions']
    if (not isinstance(compute.get('partitions'), list) or
            not isinstance(output.get('partitions'), list) or
            len(rows) != len(compute['partitions']) or
            len(rows) != len(output['partitions'])):
        raise ValueError('One compute/output work row per source partition required')
    for src, gpu, writer in zip(rows, compute['partitions'], output['partitions']):
        markers = src['variant_range'][1] - src['variant_range'][0]
        traits = src['trait_range'][1] - src['trait_range'][0]
        for name, other in [('compute', gpu), ('output', writer)]:
            if (other.get('id') != src['id'] or other.get('device') != src['device'] or
                    other.get('markers') != markers):
                raise ValueError(name + ' partition differs from source layout')
        if gpu.get('traits') != traits or writer.get('cells') != markers * traits:
            raise ValueError('Compute/output phenotype extent differs from source layout')
    read = _bounds(source.get('source_stage_floor_seconds'), 'source')
    result = _bounds(output.get('payload_floor_seconds'), 'output')
    matrix = compute.get('compute_transfer_floor_seconds')
    if (isinstance(matrix, bool) or not isinstance(matrix, (int, float)) or
            not math.isfinite(matrix) or matrix < 0):
        raise ValueError('Finite nonnegative compute/transfer floor required')
    try:
        source_dram = source['resource_work']['dram_bytes']
        h2d = compute['total_h2d_bytes']
        d2h = output['total_d2h_payload_bytes']
        dram_capacity = source['shared_capacities']['dram']
    except (KeyError, TypeError):
        raise ValueError('Complete bound shared host-memory work required') from None
    if (isinstance(source_dram, bool) or not isinstance(source_dram, (int, float)) or
            not math.isfinite(source_dram) or source_dram < 0 or
            type(h2d) is not int or h2d < 0 or
            not isinstance(d2h, list) or len(d2h) != 2 or
            any(type(value) is not int or value < 0 for value in d2h) or
            isinstance(dram_capacity, bool) or
            not isinstance(dram_capacity, (int, float)) or
            not math.isfinite(dram_capacity) or dram_capacity <= 0):
        raise ValueError('Finite nonnegative shared host-memory work and capacity required')
    # Decoder writes the host genotype staging; H2D subsequently reads it,
    # and D2H writes a distinct host result buffer. These transfers share
    # the same declared host-memory service even when GPUs use separate links.
    writer_floor = 0.
    writer_cpu = 0.
    if writer_service is not None:
        if (source['reduction'] is not None or not isinstance(writer_service, dict) or
                writer_service.get('kind') != 'torchgwas.pgen_layout_dense_writer_service_floor.v1' or
                writer_service.get('input_identity') != source['input_identity'] or
                writer_service.get('reduction') is not None or
                writer_service.get('shared_cpu_capacity') != source['shared_capacities']['cpu'] or
                not isinstance(writer_service.get('partitions'), list) or
                len(writer_service['partitions']) != len(rows) or
                not isinstance(output.get('dense_writer_work'), list) or
                len(output['dense_writer_work']) != len(rows) or
                writer_service.get('output_work_sha256') !=
                _digest(output['dense_writer_work'])):
            raise ValueError('Matching dense writer service and source capacity required')
        for src, priced, declared in zip(rows, writer_service['partitions'],
                                         output['dense_writer_work']):
            work = declared['work']
            if (priced.get('id') != src['id'] or priced.get('device') != src['device'] or
                    priced.get('markers') != work.get('markers') or
                    priced.get('traits') != work.get('traits') or
                    priced.get('chunks') != src['chunks'] or
                    any(priced.get(name) != work.get(name) for name in
                        ('payload_bytes', 'staging_copy_calls',
                         'write_calls_minimum', 'fsync_calls'))):
                raise ValueError('Dense writer service differs from bound output work')
        writer_floor = writer_service.get('writer_service_floor_seconds')
        if (isinstance(writer_floor, bool) or not isinstance(writer_floor, (int, float)) or
                not math.isfinite(writer_floor) or writer_floor < 0):
            raise ValueError('Finite nonnegative dense writer floor required')
        writer_cpu = writer_service.get('total_writer_cpu_seconds')
        if (isinstance(writer_cpu, bool) or not isinstance(writer_cpu, (int, float)) or
                not math.isfinite(writer_cpu) or writer_cpu < 0):
            raise ValueError('Finite nonnegative dense writer CPU work required')
    selector_floor = [0., 0.]
    selector_dram = [0., 0.]
    selector_cpu = [0., 0.]
    if jagwas_selection is not None:
        if (source['reduction'] != 'jagwas' or
                not isinstance(jagwas_selection, dict) or
                jagwas_selection.get('kind') !=
                'torchgwas.pgen_layout_jagwas_selection_floor.v1' or
                jagwas_selection.get('input_identity') != source['input_identity'] or
                jagwas_selection.get('reduction') != 'jagwas' or
                jagwas_selection.get('shared_cpu_capacity') !=
                source['shared_capacities']['cpu'] or
                jagwas_selection.get('shared_dram_capacity') != dram_capacity or
                jagwas_selection.get('output_partitions_sha256') !=
                _digest(output['partitions']) or
                not isinstance(jagwas_selection.get('partitions'), list) or
                len(jagwas_selection['partitions']) != len(rows)):
            raise ValueError('Matching JAGWAS selector service required')
        for src, priced, declared in zip(rows, jagwas_selection['partitions'],
                                         output['partitions']):
            if (priced.get('id') != src['id'] or
                    priced.get('device') != src['device'] or
                    priced.get('markers') != declared['markers'] or
                    priced.get('chunks') != src['chunks'] or
                    priced.get('retained_rows') != declared['retained_rows']):
                raise ValueError('JAGWAS selector service differs from bound output')
        selector_floor = _bounds(
            jagwas_selection.get('selector_service_floor_seconds'),
            'JAGWAS selector')
        selector_dram = _bounds(
            jagwas_selection.get('total_logical_dram_bytes'),
            'JAGWAS selector DRAM work')
        selector_cpu = _bounds(jagwas_selection.get('total_cpu_seconds'),
                               'JAGWAS selector CPU work')
        if any(not math.isclose(selector_dram[i], math.fsum(
                row['logical_dram_bytes'][i]
                for row in jagwas_selection['partitions']), rel_tol=1e-12,
                abs_tol=1e-9) for i in (0, 1)):
            raise ValueError('JAGWAS selector DRAM total differs from partitions')
    archive_floor = [0., 0.]
    archive_cpu = [0., 0.]
    archive_dram = [0., 0.]
    if jagwas_archive is not None:
        if (source['reduction'] != 'jagwas' or
                output.get('jagwas_writer_fsync') is not True or
                not isinstance(jagwas_archive, dict) or
                jagwas_archive.get('kind') !=
                'torchgwas.pgen_layout_jagwas_archive_floor.v1' or
                jagwas_archive.get('input_identity') != source['input_identity'] or
                jagwas_archive.get('reduction') != 'jagwas' or
                jagwas_archive.get('jagwas_writer_fsync') is not True or
                jagwas_archive.get('output_partitions_sha256') !=
                _digest(output['partitions']) or
                jagwas_archive.get('shared_cpu_capacity') !=
                source['shared_capacities']['cpu'] or
                jagwas_archive.get('shared_dram_capacity') != dram_capacity or
                jagwas_archive.get('shared_output_capacity') !=
                output['capacities']['output_bytes_per_second'] or
                not isinstance(jagwas_archive.get('partitions'), list) or
                len(jagwas_archive['partitions']) != len(rows)):
            raise ValueError('Matching durable JAGWAS archive service required')
        for src, priced, declared in zip(rows, jagwas_archive['partitions'],
                                         output['partitions']):
            if (priced.get('id') != src['id'] or
                    priced.get('device') != src['device'] or
                    priced.get('markers') != declared['markers'] or
                    priced.get('chunks') != src['chunks'] or
                    priced.get('retained_rows') != declared['retained_rows']):
                raise ValueError('JAGWAS archive differs from bound output')
        archive_floor = _bounds(
            jagwas_archive.get('archive_service_floor_seconds'),
            'JAGWAS archive service')
        archive_cpu = _bounds(jagwas_archive.get('total_cpu_seconds'),
                              'JAGWAS archive CPU work')
        archive_dram = _bounds(jagwas_archive.get('total_logical_dram_bytes'),
                               'JAGWAS archive DRAM work')
        for key, total in (('cpu_seconds', archive_cpu),
                           ('logical_dram_bytes', archive_dram)):
            if any(not math.isclose(total[i], math.fsum(
                    row[key][i] for row in jagwas_archive['partitions']),
                    rel_tol=1e-12, abs_tol=1e-9) for i in (0, 1)):
                raise ValueError('JAGWAS archive work differs from partitions')
    significant_floor = [0., 0.]
    significant_cpu = [0., 0.]
    significant_dram = [0., 0.]
    significant_consumer = None
    if significant_archive is not None:
        if (source['reduction'] != 'significant' or
                output.get('significant_backend') not in ('host', 'device') or
                output.get('significant_writer_fsync') is not True or
                not isinstance(significant_archive, dict) or
                significant_archive.get('kind') !=
                'torchgwas.pgen_layout_significant_archive_floor.v1' or
                significant_archive.get('input_identity') != source['input_identity'] or
                significant_archive.get('reduction') != 'significant' or
                significant_archive.get('significant_backend') != output['significant_backend'] or
                significant_archive.get('device_selection_max_cells') !=
                output.get('device_selection_max_cells') or
                significant_archive.get('store_beta') != output['store_beta'] or
                significant_archive.get('significant_writer_fsync') is not True or
                significant_archive.get('output_partitions_sha256') !=
                _digest(output['partitions']) or
                significant_archive.get('shared_cpu_capacity') !=
                source['shared_capacities']['cpu'] or
                significant_archive.get('shared_dram_capacity') != dram_capacity or
                significant_archive.get('shared_output_capacity') !=
                output['capacities']['output_bytes_per_second'] or
                not isinstance(significant_archive.get('partitions'), list) or
                len(significant_archive['partitions']) != len(rows)):
            raise ValueError('Matching significant archive service required')
        for src, priced, declared in zip(rows, significant_archive['partitions'],
                                         output['partitions']):
            if (priced.get('id') != src['id'] or
                    priced.get('device') != src['device'] or
                    priced.get('markers') != declared['markers'] or
                    priced.get('traits') !=
                    src['trait_range'][1] - src['trait_range'][0] or
                    priced.get('chunks') != src['chunks'] or
                    priced.get('retained_rows') != declared['retained_rows']):
                raise ValueError('Significant archive differs from bound output')
        significant_floor = _bounds(
            significant_archive.get('archive_service_floor_seconds'),
            'significant archive service')
        significant_cpu = _bounds(significant_archive.get('total_cpu_seconds'),
                                   'significant archive CPU work')
        significant_dram = _bounds(
            significant_archive.get('total_logical_dram_bytes'),
            'significant archive DRAM work')
        significant_consumer = _bounds(
            significant_archive.get('single_consumer_floor_seconds'),
            'significant archive consumer')
        for key, total in (('cpu_seconds', significant_cpu),
                           ('logical_dram_bytes', significant_dram)):
            if any(not math.isclose(total[i], math.fsum(
                    row[key][i] for row in significant_archive['partitions']),
                    rel_tol=1e-12, abs_tol=1e-9) for i in (0, 1)):
                raise ValueError('Significant archive work differs from partitions')
    host_selector_floor = [0., 0.]
    host_selector_cpu = [0., 0.]
    host_selector_dram = [0., 0.]
    if significant_host_selection is not None:
        if (source['reduction'] != 'significant' or
                output.get('significant_backend') != 'host' or
                type(output.get('significant_threshold_one')) is not bool or
                not isinstance(significant_host_selection, dict) or
                significant_host_selection.get('kind') !=
                'torchgwas.pgen_layout_significant_host_selection_floor.v1' or
                significant_host_selection.get('input_identity') !=
                source['input_identity'] or
                significant_host_selection.get('reduction') != 'significant' or
                significant_host_selection.get('significant_backend') != 'host' or
                significant_host_selection.get('return_beta') is not True or
                significant_host_selection.get('significant_threshold_one') !=
                output['significant_threshold_one'] or
                significant_host_selection.get('output_partitions_sha256') !=
                _digest(output['partitions']) or
                significant_host_selection.get('shared_cpu_capacity') !=
                source['shared_capacities']['cpu'] or
                significant_host_selection.get('shared_dram_capacity') !=
                dram_capacity or
                not isinstance(significant_host_selection.get('partitions'), list) or
                len(significant_host_selection['partitions']) != len(rows)):
            raise ValueError('Matching host significant selector service required')
        for src, priced, declared in zip(rows,
                                         significant_host_selection['partitions'],
                                         output['partitions']):
            if (priced.get('id') != src['id'] or
                    priced.get('device') != src['device'] or
                    priced.get('markers') != declared['markers'] or
                    priced.get('traits') !=
                    src['trait_range'][1] - src['trait_range'][0] or
                    priced.get('chunks') != src['chunks'] or
                    priced.get('retained_rows') != declared['retained_rows']):
                raise ValueError('Host significant selector differs from output')
        host_selector_floor = _bounds(
            significant_host_selection.get('selector_service_floor_seconds'),
            'host significant selector service')
        host_selector_cpu = _bounds(
            significant_host_selection.get('total_cpu_seconds'),
            'host significant selector CPU work')
        host_selector_dram = _bounds(
            significant_host_selection.get('total_logical_dram_bytes'),
            'host significant selector DRAM work')
        for key, total in (('cpu_seconds', host_selector_cpu),
                           ('logical_dram_bytes', host_selector_dram)):
            if any(not math.isclose(total[i], math.fsum(
                    row[key][i]
                    for row in significant_host_selection['partitions']),
                    rel_tol=1e-12, abs_tol=1e-9) for i in (0, 1)):
                raise ValueError('Host significant selector work differs from partitions')
    device_count_floor = 0.
    device_count_totals = None
    if significant_device_count is not None:
        if (source['reduction'] != 'significant' or
                output.get('significant_backend') != 'device' or
                not isinstance(significant_device_count, dict) or
                significant_device_count.get('kind') !=
                'torchgwas.pgen_layout_significant_device_count_floor.v1' or
                significant_device_count.get('input_identity') != source['input_identity'] or
                significant_device_count.get('reduction') != 'significant' or
                significant_device_count.get('significant_backend') != 'device' or
                significant_device_count.get('device_selection_max_cells') !=
                output.get('device_selection_max_cells') or
                significant_device_count.get('output_partitions_sha256') !=
                _digest(output['partitions']) or
                not isinstance(significant_device_count.get('partitions'), list) or
                len(significant_device_count['partitions']) != len(rows)):
            raise ValueError('Matching device count-transfer floor required')
        device_totals = {row['device']: 0. for row in rows}
        for src, priced, declared in zip(rows, significant_device_count['partitions'],
                                         output['partitions']):
            if (priced.get('id') != src['id'] or
                    priced.get('device') != src['device'] or
                    priced.get('markers') != declared['markers'] or
                    priced.get('traits') !=
                    src['trait_range'][1] - src['trait_range'][0] or
                    priced.get('chunks') != src['chunks'] or
                    priced.get('selection_blocks') !=
                    declared.get('device_selection_blocks') or
                    priced.get('count_d2h_bytes') !=
                    4 * declared.get('device_selection_blocks', -1)):
                raise ValueError('Device count transfer differs from bound output')
            duration = priced.get('serial_count_transfer_seconds')
            if (isinstance(duration, bool) or not isinstance(duration, (int, float)) or
                    not math.isfinite(duration) or duration < 0):
                raise ValueError('Finite device count-transfer service required')
            device_totals[src['device']] += duration
        declared_totals = significant_device_count.get('per_device_serial_seconds')
        device_count_floor = significant_device_count.get('count_barrier_floor_seconds')
        if (not isinstance(declared_totals, dict) or
                set(declared_totals) != set(device_totals) or
                any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                    not math.isfinite(value) or value < 0
                    for value in declared_totals.values()) or
                any(not math.isclose(device_totals[device], declared_totals[device],
                                     rel_tol=1e-12, abs_tol=1e-9)
                    for device in device_totals) or
                isinstance(device_count_floor, bool) or
                not isinstance(device_count_floor, (int, float)) or
                not math.isfinite(device_count_floor) or
                not math.isclose(device_count_floor, max(device_totals.values()),
                                 rel_tol=1e-12, abs_tol=1e-9)):
            raise ValueError('Device count-transfer floor differs from partitions')
        device_count_totals = device_totals
    device_launch_floor = None
    device_selector_serial = None
    if significant_device_launch is not None:
        launch = significant_device_launch
        if (source['reduction'] != 'significant' or
                output.get('significant_backend') != 'device' or
                not isinstance(launch, dict) or
                launch.get('kind') !=
                'torchgwas.pgen_layout_significant_device_launch_floor.v1' or
                launch.get('input_identity') != source['input_identity'] or
                launch.get('reduction') != 'significant' or
                launch.get('significant_backend') != 'device' or
                launch.get('device_selection_max_cells') !=
                output.get('device_selection_max_cells') or
                launch.get('output_partitions_sha256') !=
                _digest(output['partitions']) or
                not isinstance(launch.get('partitions'), list) or
                len(launch['partitions']) != len(rows)):
            raise ValueError('Matching device selector launch floor required')
        launch_totals = {row['device']: [0., 0.] for row in rows}
        for src, priced, declared in zip(rows, launch['partitions'],
                                         output['partitions']):
            if (priced.get('id') != src['id'] or
                    priced.get('device') != src['device'] or
                    priced.get('markers') != declared['markers'] or
                    priced.get('traits') !=
                    src['trait_range'][1] - src['trait_range'][0] or
                    priced.get('chunks') != src['chunks'] or
                    priced.get('selection_blocks') !=
                    declared.get('device_selection_blocks') or
                    priced.get('maximum_selection_block_cells') !=
                    declared.get('maximum_selection_block_cells') or
                    priced.get('retained_rows') != declared.get('retained_rows')):
                raise ValueError('Device selector launches differ from bound output')
            seconds = _bounds(priced.get('serial_launch_seconds'),
                              'device selector launch service')
            for index in (0, 1):
                launch_totals[src['device']][index] += seconds[index]
        declared_totals = launch.get('per_device_serial_seconds')
        device_launch_floor = _bounds(launch.get('selector_launch_floor_seconds'),
                                      'device selector launch')
        if (not isinstance(declared_totals, dict) or
                set(declared_totals) != set(launch_totals)):
            raise ValueError('Device selector launch totals differ from partitions')
        for value in declared_totals.values():
            _bounds(value, 'per-device selector launch')
        if (
                any(not math.isclose(launch_totals[device][index],
                                     declared_totals[device][index],
                                     rel_tol=1e-12, abs_tol=1e-9)
                    for device in launch_totals for index in (0, 1)) or
                any(not math.isclose(device_launch_floor[index],
                                     max(value[index]
                                         for value in launch_totals.values()),
                                     rel_tol=1e-12, abs_tol=1e-9)
                    for index in (0, 1))):
            raise ValueError('Device selector launch floor differs from partitions')
        device_selector_serial = [max(
            launch_totals[device][index] +
            (0. if device_count_totals is None else
             device_count_totals[device])
            for device in launch_totals) for index in (0, 1)]
        if not all(math.isfinite(value) for value in device_selector_serial):
            raise ValueError('Device selector serial service overflow')
    shared_dram = [(source_dram + h2d + d2h[i] + selector_dram[i] +
                    archive_dram[i] + significant_dram[i] +
                    host_selector_dram[i]) /
                   dram_capacity for i in (0, 1)]
    if not all(math.isfinite(value) for value in shared_dram):
        raise ValueError('Shared host-memory floor overflow')
    source_cpu = _bounds(source['resource_work']['cpu_seconds'], 'source CPU work')
    shared_cpu = [(source_cpu[i] + writer_cpu + selector_cpu[i] +
                   archive_cpu[i] + significant_cpu[i] +
                   host_selector_cpu[i]) /
                  source['shared_capacities']['cpu'] for i in (0, 1)]
    if not all(math.isfinite(value) for value in shared_cpu):
        raise ValueError('Shared CPU floor overflow')
    joint_consumer = None
    if jagwas_selection is not None and jagwas_archive is not None:
        selection_serial = _bounds(
            jagwas_selection.get('single_consumer_cpu_floor_seconds'),
            'JAGWAS selector consumer')
        archive_serial = _bounds(
            jagwas_archive.get('single_consumer_floor_seconds'),
            'JAGWAS archive consumer')
        joint_consumer = [selection_serial[i] + archive_serial[i]
                          for i in (0, 1)]
    envelope = [max(read[i], matrix, result[i], shared_dram[i], writer_floor,
                    selector_floor[i], archive_floor[i], significant_floor[i],
                    host_selector_floor[i], device_count_floor,
                    0. if device_selector_serial is None else device_selector_serial[i],
                    shared_cpu[i],
                    0. if joint_consumer is None else joint_consumer[i])
                for i in (0, 1)]
    return dict(kind='torchgwas.pgen_layout_partial_envelope.v1',
                input_identity=dict(source['input_identity']),
                reduction=source['reduction'],
                issued_revision=coverage['issued_revision'],
                written_events=coverage['written_events'],
                required_cells=coverage['required_cells'],
                coverage=coverage,
                stage_floor_seconds=dict(source=list(read),
                                         h2d_and_matrix=[matrix, matrix],
                                         d2h_and_output=list(result),
                                         combined_shared_dram=shared_dram,
                                         combined_shared_cpu=shared_cpu,
                                         dense_writer_service=([writer_floor] * 2
                                                               if writer_service is not None else None),
                                         jagwas_selection_service=(list(selector_floor)
                                                                   if jagwas_selection is not None else None),
                                         jagwas_archive_service=(list(archive_floor)
                                                                 if jagwas_archive is not None else None),
                                         jagwas_single_consumer=joint_consumer,
                                         significant_archive_service=(list(significant_floor)
                                                                      if significant_archive is not None else None),
                                         significant_single_consumer=significant_consumer,
                                         significant_host_selector_service=(list(host_selector_floor)
                                                                            if significant_host_selection is not None else None),
                                         significant_device_count_barrier=([device_count_floor] * 2
                                                                           if significant_device_count is not None else None),
                                         significant_device_launch_service=device_launch_floor,
                                         significant_device_launch_count_serial=device_selector_serial),
                partial_floor_seconds=envelope,
                prediction_complete=False, selection_validated=False,
                scope='Exact unissued coverage plus conditional necessary read/decode, unpacked PGEN H2D/GEMM and native result/output payload floors. Shared CPU and host-memory constraints sum distinct source, DMA, selector and writer work. Optional bound dense-writer, full-panel JAGWAS selector/archive, significant selector/archive and device count-transfer/mandatory CUDA launch terms add independently priced service. Remaining device kernel bodies and synchronization service, queueing, filesystem throttling, in-flight work and final drain can raise completion time. The high endpoint is not a completion ceiling; no switch is authorized.')
