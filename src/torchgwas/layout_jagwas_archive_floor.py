"""Compact conditional indexed-archive service for JAGWAS variant shards."""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .reduced_output_work import jagwas_indexed_part_work


_SCHEMA = [['variant_index', '<i8'], ['chi2', '<f8']]


def _price(name, value, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value < 0 or (positive and value == 0)):
        raise ValueError('Finite independent JAGWAS archive price required: ' + name)
    return value


def _profile(row):
    if not isinstance(row, dict):
        raise ValueError('Explicit JAGWAS archive profile required')
    try:
        q = _price('cpu_fraction', row['cpu_fraction'], positive=True)
        dram = _price('shared_dram_bytes_per_second',
                      row['shared_dram_bytes_per_second'], positive=True)
        fsync = _price('fsync_seconds', row['fsync_seconds'])
        writeback = row['writeback_service']
        page = _price('pagecache_seconds_per_byte',
                      writeback['pagecache_seconds_per_byte'])
        storage = _price('storage_seconds_per_byte',
                         writeback['storage_seconds_per_byte'], positive=True)
    except (KeyError, TypeError):
        raise ValueError('Complete JAGWAS archive service profile required') from None
    if q > 1:
        raise ValueError('JAGWAS archive CPU fraction exceeds one')
    return q, dram, fsync, page, storage


def native_layout_jagwas_archive_floor(source, output, archive_price,
                                       profiles_by_device):
    """Bound existing two-array NPZ service from shard-level retained ranges.

    The number of nonempty chunks and their NPY header lengths can vary with
    the unobserved distribution of retained variants. Bounds use the smallest
    and largest possible nonempty part overhead, not a selectivity estimate.
    All indexed parts use the public one-consumer path and durable fsync.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            source.get('reduction') != 'jagwas' or
            source.get('partition_axis') != 'variant' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1' or
            output.get('reduction') != 'jagwas' or
            output.get('jagwas_writer_fsync') is not True or
            output.get('input_identity') != source.get('input_identity')):
        raise ValueError('Matching durable full-panel JAGWAS source/output required')
    identity = source['input_identity']
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed before JAGWAS archive pricing')
    if (not isinstance(archive_price, dict) or
            set(archive_price) != {'field_schema', 'call_cpu_seconds',
                                   'byte_cpu_seconds'} or
            archive_price['field_schema'] != _SCHEMA):
        raise ValueError('Independent two-array JAGWAS archive price required')
    fixed = _price('archive call CPU', archive_price['call_cpu_seconds'])
    bulk = _price('archive byte CPU', archive_price['byte_cpu_seconds'])
    partitions = source.get('partitions')
    output_rows = output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(output_rows, list) or len(output_rows) != len(partitions)):
        raise ValueError('One output row per JAGWAS variant shard required')
    devices = {row['device'] for row in partitions}
    if not isinstance(profiles_by_device, dict) or set(profiles_by_device) != devices:
        raise ValueError('One JAGWAS archive service profile per active GPU required')
    profiles = {device: _profile(profile)
                for device, profile in profiles_by_device.items()}
    cpu_cap = _price('shared CPU capacity',
                     source['shared_capacities']['cpu'], positive=True)
    dram_cap = _price('shared DRAM capacity',
                      source['shared_capacities']['dram'], positive=True)
    output_cap = _price('shared output capacity',
                        output['capacities']['output_bytes_per_second'], positive=True)
    rows = []
    totals = {key: [[], []] for key in ('cpu_seconds', 'logical_dram_bytes',
                                        'file_bytes', 'nonempty_parts',
                                        'serial_consumer_floor_seconds')}
    for partition, result in zip(partitions, output_rows):
        markers = partition['variant_range'][1] - partition['variant_range'][0]
        block = source['chunk_markers']
        chunks = (markers + block - 1) // block
        retained = result.get('retained_rows')
        if (result.get('id') != partition['id'] or
                result.get('device') != partition['device'] or
                result.get('markers') != markers or
                result.get('cells') != markers * source['total_traits']):
            raise ValueError('JAGWAS archive geometry differs from bound output')
        if (partition['trait_range'] != [0, source['total_traits']] or
                partition['chunks'] != chunks or
                not isinstance(retained, list) or len(retained) != 2 or
                any(type(value) is not int for value in retained) or
                not 0 <= retained[0] <= retained[1] <= markers):
            raise ValueError('JAGWAS archive occupancy differs from bound output')
        if result.get('output_array_payload_bytes') != [16 * value for value in retained]:
            raise ValueError('JAGWAS array payload differs from bound output')
        # A shard-level occupancy cap also bounds any one nonempty part.
        max_rows = max(1, min(markers, block, retained[1]))
        smallest = jagwas_indexed_part_work(1)
        largest = jagwas_indexed_part_work(max_rows)
        local = sum(row['local_header_bytes'] for row in smallest['arrays'])
        overhead = [smallest['file_bytes'] - 16,
                    largest['file_bytes'] - 16 * max_rows]
        # NPY header length is nondecreasing with the printed shape length;
        # ZIP framing and each local header are constant for the two fields.
        parts = [(retained[0] + block - 1) // block,
                 min(retained[1], chunks)]
        file_bytes = [16 * retained[i] + parts[i] * overhead[i]
                      for i in (0, 1)]
        submitted = [file_bytes[i] + parts[i] * local for i in (0, 1)]
        q, host_dram, fsync, page, storage = profiles[partition['device']]
        cpu = [parts[i] * fixed + 16 * retained[i] * bulk +
               submitted[i] * page for i in (0, 1)]
        dram = [5 * 16 * retained[i] +
                2 * (submitted[i] - 16 * retained[i]) for i in (0, 1)]
        # Sum(max(CPU, DRAM)) is bounded below by max of the sums, and above
        # by their sum. Storage and fsync are serial after each NPZ CPU step.
        serial = [max(cpu[0] / q, dram[0] / host_dram) +
                  file_bytes[0] * storage + parts[0] * fsync,
                  cpu[1] / q + dram[1] / host_dram +
                  file_bytes[1] * storage + parts[1] * fsync]
        if not all(math.isfinite(value) for value in (*cpu, *serial)):
            raise ValueError('JAGWAS archive service overflow')
        data = dict(id=partition['id'], device=partition['device'],
                    markers=markers, chunks=chunks, retained_rows=list(retained),
                    nonempty_parts=parts, file_bytes=file_bytes,
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
        raise ValueError('JAGWAS archive floor overflow')
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during JAGWAS archive pricing')
    return dict(kind='torchgwas.pgen_layout_jagwas_archive_floor.v1',
                input_identity=dict(identity), reduction='jagwas',
                output_partitions_sha256=_digest(output_rows),
                jagwas_writer_fsync=True, partitions=rows,
                total_cpu_seconds=sums['cpu_seconds'],
                total_logical_dram_bytes=sums['logical_dram_bytes'],
                total_file_bytes=sums['file_bytes'],
                total_nonempty_parts=sums['nonempty_parts'],
                single_consumer_floor_seconds=sums['serial_consumer_floor_seconds'],
                shared_cpu_capacity=cpu_cap, shared_dram_capacity=dram_cap,
                shared_output_capacity=output_cap,
                archive_service_floor_seconds=floor,
                prediction_complete=False, selection_validated=False,
                scope='Conditional two-array JAGWAS NPZ archive CPU, logical DRAM, durable file bytes and one-consumer service from independent primitive prices. Shard-level retained intervals bound possible nonempty part and NPY-header counts; upper bounds this partial service only. Queues, filesystem throttling, final metadata and in-flight work remain outside; no JIT switch is authorized.')
