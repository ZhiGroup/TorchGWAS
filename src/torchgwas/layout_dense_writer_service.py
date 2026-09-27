"""Compact mandatory CPU and serial-close work of the native dense writer."""
import math

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .output_write_work import writer_copy_cost


def _price(name, value, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value < 0 or (positive and value == 0)):
        raise ValueError('Finite nonnegative independent writer price required: ' + name)
    return value


def _service(profile):
    if not isinstance(profile, dict):
        raise ValueError('Explicit dense writer profile required')
    try:
        q = _price('cpu_fraction', profile['cpu_fraction'], positive=True)
        if q > 1:
            raise ValueError('Writer CPU fraction exceeds one')
        copy_source = ('writer_copy_service' if profile.get('writer_copy_service') is not None
                       else 'legacy_numpy_copy_bytes_without_call_price')
        copy = writer_copy_cost(profile)
        zero = profile['process_units']['bytearray_zero_bytes']
        executor = profile['executor_cpu_seconds']
        writeback = profile['writeback_service']
        required = {'pagecache_seconds_per_byte', 'storage_seconds_per_byte',
                    'submit_seconds', 'wait_seconds', 'fadvise_seconds'}
        if (not isinstance(writeback, dict) or
                set(writeback) not in (required, required | {'fadvise_eviction_seconds_per_byte'})):
            raise ValueError('Complete existing writeback service required')
        fsync = profile['fsync_seconds']
    except (KeyError, TypeError):
        raise ValueError('Complete dense writer service profile required') from None
    for name, value in copy.items():
        _price('writer_copy_service.' + name, value)
    _price('bytearray_zero_bytes', zero)
    _price('executor_cpu_seconds', executor)
    for name, value in writeback.items():
        _price('writeback_service.' + name, value,
               positive=name == 'storage_seconds_per_byte')
    _price('fsync_seconds', fsync)
    return q, copy, zero, executor, writeback, fsync, copy_source


def native_layout_dense_writer_service_floor(source, output, profiles_by_device):
    """Price fixed dense writer counts without expanding output events.

    Prices are the existing finite graph's independent primitive inputs. The
    CPU work sum and per-device serial writer chains are necessary conditional
    floors; storage transfer is already in the output payload report. They
    omit queue stalls, dirty-page coupling and final metadata publication.
    """
    if (not isinstance(source, dict) or
            source.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            not isinstance(output, dict) or
            output.get('kind') != 'torchgwas.pgen_layout_output_floor.v1' or
            source.get('reduction') is not None or output.get('reduction') is not None or
            source.get('input_identity') != output.get('input_identity')):
        raise ValueError('Matching bound dense source and output reports required')
    identity = source['input_identity']
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed before dense writer pricing')
    partitions = source.get('partitions')
    writers = output.get('dense_writer_work')
    output_rows = output.get('partitions')
    if (not isinstance(partitions, list) or not partitions or
            not isinstance(writers, list) or len(writers) != len(partitions) or
            not isinstance(output_rows, list) or len(output_rows) != len(partitions)):
        raise ValueError('One explicit compact dense writer per source partition required')
    devices = {row['device'] for row in partitions}
    if not isinstance(profiles_by_device, dict) or set(profiles_by_device) != devices:
        raise ValueError('One dense writer service profile per active device required')
    service = {device: _service(profile) for device, profile in profiles_by_device.items()}
    capacity = _price('shared CPU capacity', source['shared_capacities']['cpu'], positive=True)
    rows = []
    cpu_parts = []
    chains = {device: [] for device in devices}
    for partition, declared, output_row in zip(partitions, writers, output_rows):
        work = declared.get('work') if isinstance(declared, dict) else None
        if (not isinstance(work, dict) or
                work.get('kind') != 'torchgwas.compact_binary_output_work.v1' or
                declared.get('id') != partition['id'] or
                declared.get('device') != partition['device'] or
                output_row.get('id') != partition['id'] or
                output_row.get('device') != partition['device']):
            raise ValueError('Dense writer differs from bound source partition')
        markers = partition['variant_range'][1] - partition['variant_range'][0]
        traits = partition['trait_range'][1] - partition['trait_range'][0]
        chunks = partition['chunks']
        if (work.get('markers') != markers or work.get('traits') != traits or
                work.get('chunk_markers') != source['chunk_markers'] or
                chunks != (markers + work['chunk_markers'] - 1) // work['chunk_markers'] or
                output_row.get('output_array_payload_bytes') != [work.get('payload_bytes')] * 2 or
                work.get('fsync_calls') != work.get('arrays')):
            raise ValueError('Dense writer geometry or closing fsync differs from output')
        q, copy, zero, executor, writeback, fsync, copy_source = service[partition['device']]
        setup_and_consumer = (
            work['zero_initialization_bytes'] * zero +
            work['staging_copy_bytes'] * copy['cpu_seconds_per_byte'] +
            work['staging_copy_calls'] * copy['cpu_seconds_per_call'] +
            (chunks + work['write_calls_minimum']) * executor)
        stream_rows = {}
        for name, stream in work['streams'].items():
            writer_cpu = (
                stream['payload_bytes'] * writeback['pagecache_seconds_per_byte'] +
                stream['writeback_submit_calls'] * writeback['submit_seconds'] +
                stream['writeback_wait_calls'] * writeback['wait_seconds'] +
                stream['writeback_fadvise_calls'] * writeback['fadvise_seconds'] +
                stream['writeback_waited_bytes'] *
                writeback.get('fadvise_eviction_seconds_per_byte', 0.))
            stream_rows[name] = dict(cpu_seconds=writer_cpu,
                                     serial_seconds=writer_cpu / q + fsync)
        cpu = setup_and_consumer + sum(row['cpu_seconds'] for row in stream_rows.values())
        # Each stream writes serially. Their writes may overlap, but close
        # calls are sequential and the next tile on a GPU starts after drain.
        serial = max(max(row['serial_seconds'] for row in stream_rows.values()),
                     len(stream_rows) * fsync)
        if not math.isfinite(cpu) or not math.isfinite(serial):
            raise ValueError('Dense writer service overflow')
        cpu_parts.append(cpu)
        chains[partition['device']].append(serial)
        rows.append(dict(id=partition['id'], device=partition['device'],
                         markers=markers, traits=traits, chunks=chunks,
                         payload_bytes=work['payload_bytes'],
                         staging_copy_calls=work['staging_copy_calls'],
                         write_calls_minimum=work['write_calls_minimum'],
                         fsync_calls=work['fsync_calls'],
                         copy_price_source=copy_source,
                         cpu_seconds=cpu, consumer_cpu_seconds=setup_and_consumer,
                         stream_service=stream_rows,
                         serial_writer_floor_seconds=serial))
    cpu_total = math.fsum(cpu_parts)
    cpu_floor = cpu_total / capacity
    device_chains = {device: math.fsum(values) for device, values in chains.items()}
    floor = max(cpu_floor, *device_chains.values())
    if input_identity(identity['path']) != identity:
        raise ValueError('PGEN input changed during dense writer pricing')
    return dict(kind='torchgwas.pgen_layout_dense_writer_service_floor.v1',
                input_identity=dict(identity), reduction=None,
                output_work_sha256=_digest(writers),
                partitions=rows, total_writer_cpu_seconds=cpu_total,
                shared_cpu_capacity=capacity,
                shared_cpu_floor_seconds=cpu_floor,
                per_device_serial_writer_floor_seconds=device_chains,
                writer_service_floor_seconds=floor,
                unpriced_terms=(['Writer staging-copy dispatch uses the existing zero-call-cost NumPy fallback; a source-matched writer-specific call price is absent.']
                                if any(row['copy_price_source'] != 'writer_copy_service'
                                       for row in rows) else []),
                prediction_complete=False, selection_validated=False,
                scope='Conditional necessary source-matched dense writer CPU-equivalent work and per-device serial write/close chains, using the finite graph primitive price fields and compact exact counts. Storage payload is in the separate output floor. Queueing, dirty-page coupling, short writes, live writer state and final metadata can raise completion time. Not a complete completion bound or switch authorization.')
