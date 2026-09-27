"""Price output already accepted by one observed native dense writer."""
import math

from .calibration_cache import _digest
from .layout_dense_writer_service import _service


def _bytes(value, name):
    if type(value) is not int or value < 0:
        raise ValueError('Nonnegative integer writer ' + name + ' required')
    return value


def price_dense_writer_queue_observation(observation, profile):
    """Price the accepted bytes awaiting writes at one writer-local instant.

    Interval endpoints bound bytes, not elapsed seconds: coefficients are
    independent calibrated service inputs, rather than performance guarantees.
    The event may be historical by the time its boundary is inspected.
    """
    if (not isinstance(observation, dict) or
            observation.get('kind') != 'torchgwas.dense_writer_queue_observation.v1' or
            observation.get('valid') is not True or
            observation.get('atomic_writer_streams') is not True or
            not isinstance(observation.get('streams'), dict) or
            not observation['streams'] or
            set(observation['streams']) - {'beta', 't_stat', 'neg_log10_p', 'df'} or
            't_stat' not in observation['streams']):
        raise ValueError('One valid atomic native dense writer observation required')
    started = observation.get('capture_started_seconds')
    finished = observation.get('capture_finished_seconds')
    if (any(type(value) not in (int, float) or not math.isfinite(value)
            for value in (started, finished)) or started > finished):
        raise ValueError('Finite ordered writer capture times required')
    q, _, _, _, writeback, _, _ = _service(profile)
    pagecache = writeback['pagecache_seconds_per_byte']
    storage = writeback['storage_seconds_per_byte']
    streams = {}
    for name, state in observation['streams'].items():
        if not isinstance(state, dict) or state.get('error') is not False:
            raise ValueError('Valid native writer stream counters required')
        accepted, written, staging, queued, active = (
            _bytes(state.get(field), field) for field in
            ('accepted_bytes', 'written_bytes', 'staging_bytes',
             'queued_bytes', 'active_bytes'))
        if (written > accepted or accepted-written != staging+queued+active or
                state.get('pending_write_bytes_interval') !=
                [staging+queued, accepted-written]):
            raise ValueError('Writer accepted/pending byte invariant failed')
        low, high = staging+queued, accepted-written
        streams[name] = dict(
            accepted_bytes=accepted, completely_written_bytes=written,
            pending_write_bytes_interval=[low, high],
            pagecache_cpu_seconds_at_fixed_price=[low*pagecache, high*pagecache],
            storage_seconds_at_fixed_price=[low*storage, high*storage],
            stream_serial_pagecache_seconds_at_fixed_fraction=
                [low*pagecache/q, high*pagecache/q])
    low = sum(row['pending_write_bytes_interval'][0] for row in streams.values())
    high = sum(row['pending_write_bytes_interval'][1] for row in streams.values())
    cpu = [math.fsum(row['pagecache_cpu_seconds_at_fixed_price'][i]
                     for row in streams.values()) for i in (0, 1)]
    disk = [math.fsum(row['storage_seconds_at_fixed_price'][i]
                      for row in streams.values()) for i in (0, 1)]
    if not all(math.isfinite(value) for value in (*cpu, *disk)):
        raise ValueError('Dense writer queue service overflow')
    return dict(kind='torchgwas.productive_dense_writer_queue_service.v1',
                capture_started_seconds=started,
                capture_finished_seconds=finished,
                profile_sha256=_digest(profile),
                streams=streams,
                pending_write_bytes_interval=[low, high],
                pagecache_cpu_seconds_at_fixed_price=cpu,
                storage_seconds_at_fixed_price=disk,
                prediction_complete=False, selection_validated=False,
                scope='One historical writer-local atomic observation. Lower bytes are staged or queued and have not entered os.write; upper bytes also include an active block that may already be partly written. Independent fixed page-cache and storage prices convert this accepted-byte interval to nominal work only. It omits unaccepted output, producer/GPU work, dirty data already written, writeback calls, fsync, other writers and capacity contention. It is neither a synchronized checkpoint nor a completion-time bound.')


def price_bracketed_dense_writer_queues(bracket, profiles_by_device):
    """Price the common-instant upper accepted write workload across writers."""
    if (not isinstance(bracket, dict) or
            bracket.get('kind') != 'torchgwas.bracketed_dense_writer_queues.v1' or
            not isinstance(bracket.get('writers'), dict) or not bracket['writers'] or
            not isinstance(profiles_by_device, dict) or
            type(bracket.get('anchor_seconds')) not in (int,float) or
            not math.isfinite(bracket['anchor_seconds']) or
            any(not isinstance(row,dict) or
                not isinstance(row.get('device'),str) or not row['device']
                for row in bracket['writers'].values())):
        raise ValueError('Common-instant writer bracket and profiles required')
    devices = {row['device'] for row in bracket['writers'].values()}
    if (set(profiles_by_device) != devices or
            not isinstance(bracket.get('logical_pending_bytes_interval'), list)):
        raise ValueError('One independent profile per active writer device required')
    services = {device: _service(profiles_by_device[device])
                for device in devices}
    writers = {}
    total_bytes = total_cpu = total_storage = 0
    for directory, row in bracket['writers'].items():
        if (not isinstance(directory, str) or not directory or
                not isinstance(row, dict) or not isinstance(row.get('streams'), dict) or
                not row['streams']):
            raise ValueError('Bound writer stream workload required')
        q, _, _, _, writeback, _, _ = services[row['device']]
        streams = {}
        for name, pair in row['streams'].items():
            if (name not in ('beta', 't_stat', 'neg_log10_p', 'df') or
                    not isinstance(pair, list) or len(pair) != 2 or
                    any(type(value) is not int or value < 0 for value in pair) or
                    pair[0] > pair[1]):
                raise ValueError('Bound logical writer byte interval required')
            upper = pair[1]
            streams[name] = dict(logical_pending_bytes_interval=list(pair),
                os_write_bytes_upper=upper,
                pagecache_cpu_seconds_upper_at_fixed_price=
                    upper * writeback['pagecache_seconds_per_byte'],
                storage_seconds_upper_at_fixed_price=
                    upper * writeback['storage_seconds_per_byte'],
                stream_serial_pagecache_seconds_upper_at_fixed_fraction=
                    upper * writeback['pagecache_seconds_per_byte'] / q)
        bytes_upper = sum(value['os_write_bytes_upper'] for value in streams.values())
        if row.get('logical_pending_bytes_interval') != [
                sum(value['logical_pending_bytes_interval'][i]
                    for value in streams.values()) for i in (0, 1)]:
            raise ValueError('Writer bracket aggregate differs from streams')
        cpu = math.fsum(value['pagecache_cpu_seconds_upper_at_fixed_price']
                        for value in streams.values())
        storage = math.fsum(value['storage_seconds_upper_at_fixed_price']
                            for value in streams.values())
        total_bytes += bytes_upper
        total_cpu += cpu
        total_storage += storage
        writers[directory] = dict(device=row['device'],
            profile_sha256=_digest(profiles_by_device[row['device']]),
            streams=streams,os_write_bytes_upper=bytes_upper,
            pagecache_cpu_seconds_upper_at_fixed_price=cpu,
            storage_seconds_upper_at_fixed_price=storage)
    if (total_bytes != bracket.get('os_write_bytes_upper') or
            not all(math.isfinite(value) for value in (total_cpu, total_storage))):
        raise ValueError('Writer bracket or fixed service price overflow')
    return dict(kind='torchgwas.priced_bracketed_dense_writer_queues.v1',
        anchor_seconds=bracket['anchor_seconds'],writers=writers,
        os_write_bytes_upper=total_bytes,
        pagecache_cpu_seconds_upper_at_fixed_price=total_cpu,
        storage_seconds_upper_at_fixed_price=total_storage,
        prediction_complete=False,selection_validated=False,
        scope='Common-instant upper workload for bytes already accepted by registered dense writers, priced with independent fixed page-cache and storage coefficients. Active blocks may already be partly written. This excludes unaccepted output, dirty data previously written, writeback calls, final fsync, producer/GPU work and contention. It is not an elapsed completion ceiling or JIT switch authorization.')
