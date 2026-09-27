"""Compact indexed-archive service for significant variant-trait pairs."""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .reduced_output_work import significant_output_work
from .significant_host_work import indexed_part_work


def _number(name, value, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value < 0 or (positive and value == 0)):
        raise ValueError('Finite independent significant archive price required: ' + name)
    return value


def _profile(row):
    if not isinstance(row, dict):
        raise ValueError('Explicit significant archive profile required')
    try:
        q = _number('cpu_fraction', row['cpu_fraction'], positive=True)
        dram = _number('shared_dram_bytes_per_second',
                       row['shared_dram_bytes_per_second'], positive=True)
        fsync = _number('fsync_seconds', row['fsync_seconds'])
        page = _number('pagecache_seconds_per_byte',
                       row['writeback_service']['pagecache_seconds_per_byte'])
        storage = _number('storage_seconds_per_byte',
                          row['writeback_service']['storage_seconds_per_byte'],
                          positive=True)
        memcpy = _number('numpy_copy_bytes',
                         row['process_units']['numpy_copy_bytes'])
    except (KeyError, TypeError):
        raise ValueError('Complete significant archive service profile required') from None
    if q > 1:
        raise ValueError('Significant archive CPU fraction exceeds one')
    return q, dram, fsync, page, storage, memcpy


def native_layout_significant_archive_floor(source, output, archive_prices,
                                            profiles_by_device):
    """Bound durable NPZ parts from retained-pair intervals.

    Host selection writes at most one part per source chunk; device selection
    writes one per nonempty selection block. Unknown retained distribution
    widens framing, fsync and serial-writer service without expanding events.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            source.get('reduction') != 'significant' or
            source.get('partition_axis') != 'trait' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1' or
            output.get('reduction') != 'significant' or
            output.get('significant_backend') not in ('host', 'device') or
            output.get('significant_writer_fsync') is not True or
            output.get('input_identity') != source.get('input_identity')):
        raise ValueError('Matching durable significant source/output required')
    identity = source['input_identity']
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed before significant archive pricing')
    if not isinstance(archive_prices, dict) or set(archive_prices) != {'True', 'False'}:
        raise ValueError('Complete independent significant archive price bank required')
    beta = output.get('store_beta')
    if type(beta) is not bool:
        raise ValueError('Explicit significant beta-output setting required')
    price = archive_prices[str(beta)]
    if not isinstance(price, dict) or set(price) != {'call_cpu_seconds',
                                                     'byte_cpu_seconds'}:
        raise ValueError('Complete significant archive primitive required')
    fixed = _number('archive call CPU', price['call_cpu_seconds'])
    bulk = _number('archive byte CPU', price['byte_cpu_seconds'])
    partitions = source.get('partitions')
    output_rows = output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(output_rows, list) or len(output_rows) != len(partitions)):
        raise ValueError('One output row per significant trait tile required')
    devices = {row['device'] for row in partitions}
    if not isinstance(profiles_by_device, dict) or set(profiles_by_device) != devices:
        raise ValueError('One significant archive profile per active GPU required')
    profiles = {device: _profile(profile)
                for device, profile in profiles_by_device.items()}
    cpu_cap = _number('shared CPU capacity',
                      source['shared_capacities']['cpu'], positive=True)
    dram_cap = _number('shared DRAM capacity',
                       source['shared_capacities']['dram'], positive=True)
    output_cap = _number('shared output capacity',
                         output['capacities']['output_bytes_per_second'],
                         positive=True)
    rows = []
    totals = {key: [[], []] for key in ('cpu_seconds', 'logical_dram_bytes',
                                        'file_bytes', 'nonempty_parts',
                                        'serial_consumer_floor_seconds')}
    payload_per_pair = 28 if beta else 24
    backend = output['significant_backend']
    if backend == 'device' and (type(output.get('device_selection_max_cells')) is not int or
                                output['device_selection_max_cells'] < 1):
        raise ValueError('Bound device selection cell limit required')
    for partition, result in zip(partitions, output_rows):
        markers = partition['variant_range'][1] - partition['variant_range'][0]
        traits = partition['trait_range'][1] - partition['trait_range'][0]
        cells = markers * traits
        block = source['chunk_markers']
        chunks = (markers + block - 1) // block
        retained = result.get('retained_rows')
        if (result.get('id') != partition['id'] or
                result.get('device') != partition['device'] or
                result.get('markers') != markers or
                result.get('cells') != cells or
                partition['chunks'] != chunks or
                not isinstance(retained, list) or len(retained) != 2 or
                any(type(value) is not int for value in retained) or
                not 0 <= retained[0] <= retained[1] <= cells):
            raise ValueError('Significant archive geometry or occupancy differs from output')
        if result.get('output_array_payload_bytes') != [
                payload_per_pair * value for value in retained]:
            raise ValueError('Significant array payload differs from output')
        if backend == 'device':
            selection = significant_output_work(
                source['samples'], markers, traits, block, backend='device',
                max_selection_cells=output['device_selection_max_cells'],
                include_blocks=False)
            maximum_part_pairs = selection['maximum_selection_block_cells']
            possible_parts = selection['selection_blocks']
            if (result.get('device_selection_blocks') != possible_parts or
                    result.get('maximum_selection_block_cells') != maximum_part_pairs):
                raise ValueError('Device selection blocks differ from bound output')
        else:
            maximum_part_pairs = min(block, markers) * traits
            possible_parts = chunks
        max_rows = max(1, min(maximum_part_pairs, retained[1]))
        smallest = indexed_part_work(1, store_beta=beta)
        largest = indexed_part_work(max_rows, store_beta=beta)
        local = sum(row['local_header_bytes'] for row in smallest['arrays'])
        overhead = [smallest['file_bytes'] - payload_per_pair,
                    largest['file_bytes'] - payload_per_pair * max_rows]
        parts = [(retained[0] + maximum_part_pairs - 1) // maximum_part_pairs,
                 min(retained[1], possible_parts)]
        file_bytes = [payload_per_pair * retained[i] + parts[i] * overhead[i]
                      for i in (0, 1)]
        submitted = [file_bytes[i] + parts[i] * local for i in (0, 1)]
        q, host_dram, fsync, page, storage, memcpy = profiles[partition['device']]
        cpu = [parts[i] * fixed + bulk * payload_per_pair * retained[i] +
               page * submitted[i] + 4 * retained[i] * memcpy
               for i in (0, 1)]
        dram = [5 * payload_per_pair * retained[i] + 8 * retained[i] +
                2 * (submitted[i] - payload_per_pair * retained[i])
                for i in (0, 1)]
        serial = [max(cpu[0] / q, dram[0] / host_dram) +
                  file_bytes[0] * storage + parts[0] * fsync,
                  cpu[1] / q + dram[1] / host_dram +
                  file_bytes[1] * storage + parts[1] * fsync]
        if not all(math.isfinite(value) for value in (*cpu, *serial)):
            raise ValueError('Significant archive service overflow')
        data = dict(id=partition['id'], device=partition['device'],
                    markers=markers, traits=traits, chunks=chunks,
                    possible_parts=possible_parts,
                    maximum_part_pairs=maximum_part_pairs,
                    retained_rows=list(retained), nonempty_parts=parts,
                    file_bytes=file_bytes,
                    local_header_rewrite_bytes_per_part=local,
                    cpu_seconds=cpu, logical_dram_bytes=dram,
                    serial_consumer_floor_seconds=serial)
        rows.append(data)
        for key in totals:
            for i in (0, 1):
                totals[key][i].append(data[key][i])
    sums = {key: [math.fsum(values) for values in pair]
            for key, pair in totals.items()}
    floor = [max(sums['cpu_seconds'][i] / cpu_cap,
                 sums['logical_dram_bytes'][i] / dram_cap,
                 sums['file_bytes'][i] / output_cap,
                 sums['serial_consumer_floor_seconds'][i])
             for i in (0, 1)]
    if not all(math.isfinite(value) for value in floor):
        raise ValueError('Significant archive floor overflow')
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during significant archive pricing')
    return dict(kind='torchgwas.pgen_layout_significant_archive_floor.v1',
                input_identity=dict(identity), reduction='significant',
                significant_backend=output['significant_backend'],
                device_selection_max_cells=output['device_selection_max_cells'],
                store_beta=beta, significant_writer_fsync=True,
                output_partitions_sha256=_digest(output_rows), partitions=rows,
                total_cpu_seconds=sums['cpu_seconds'],
                total_logical_dram_bytes=sums['logical_dram_bytes'],
                total_file_bytes=sums['file_bytes'],
                total_nonempty_parts=sums['nonempty_parts'],
                single_consumer_floor_seconds=sums['serial_consumer_floor_seconds'],
                shared_cpu_capacity=cpu_cap, shared_dram_capacity=dram_cap,
                shared_output_capacity=output_cap,
                archive_service_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Conditional significant-pair NPZ archive CPU, logical DRAM, durable file bytes and one-consumer service for explicit retained-count intervals. Host output has one possible part per source chunk; device output has one per selection block. Selector service is separate. Upper bounds only this partial writer floor; queueing, throttling, final metadata and in-flight work can raise completion time. No JIT switch is authorized.')
