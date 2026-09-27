"""Independent output storage service, distinct from buffered-write CPU work."""
import math
import statistics

from .first_principles import positive


def direct_write_capacity(probe, *, workers, cpu_affinity, filesystem):
    """Read a verified fixed-size generic O_DIRECT service measurement.

    Write completion excludes the separately recorded final fsync. Both the
    transfer accounting and full readback must pass; buffered observations are
    diagnostics and cannot supply a storage price. This does not characterize
    mixed read/write service, new-extent allocation or dirty-page throttling.
    """
    if isinstance(workers, bool) or not isinstance(workers, int) or workers not in (1, 4):
        raise ValueError('Unmeasured direct-write worker context')
    if probe.get('affinity') != list(cpu_affinity):
        raise ValueError('Direct-write CPU affinity mismatch')
    inputs = probe['input']
    if inputs.get('filesystem') != filesystem:
        raise ValueError('Direct-write filesystem mismatch')
    block, count, amount = inputs.get('block_bytes'), inputs.get('blocks'), inputs.get('bytes')
    if (block, count, amount) != (16 << 20, 64, 1 << 30) or inputs.get('preallocated') is not True:
        raise ValueError('Unverified direct-write extent or allocation policy')
    rows = [row for row in probe['rows'] if row.get('mode') == 'direct' and row.get('workers') == workers]
    if len(rows) != 5 or {row['repeat'] for row in rows} != set(range(5)):
        raise ValueError('Five complete unique direct-write repeats required')
    rates = []
    for row in rows:
        if (row.get('bytes') != amount or row.get('resident_at_launch') != [0]*workers
                or row.get('resident_after') != [0]*workers or row.get('verified_payload') is not True
                or row.get('write_io_delta', {}).get('write_bytes') != amount
                or row.get('durable_io_delta', {}).get('write_bytes') != amount):
            raise ValueError('Direct-write bypass, transferred bytes or payload verification failed')
        seconds = positive('direct write elapsed', row['seconds'])
        positive('write fsync elapsed', row['fsync_seconds'], True)
        calls = row['calls']
        coverage = [(call['worker'], call['part']) for call in calls]
        if len(calls) != count or set(coverage) != {
                (worker, part) for worker in range(workers) for part in range(count//workers)}:
            raise ValueError('Direct-write request coverage mismatch')
        for call in calls:
            if call['bytes'] != block or call['offset'] != call['part']*block:
                raise ValueError('Direct-write request geometry mismatch')
            start = positive('write start', call['start'], True)
            end = positive('write end', call['end'])
            positive('write CPU service', call['cpu_seconds'], True)
            if end <= start or not math.isclose(call['wall_seconds'], end-start, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError('Direct-write interval mismatch')
        for worker in range(workers):
            sequential = sorted((call for call in calls if call['worker'] == worker), key=lambda c: c['part'])
            if any(left['end'] > right['start'] for left, right in zip(sequential, sequential[1:])):
                raise ValueError('Sequential worker requests overlap')
        if not math.isclose(seconds, max(call['end'] for call in calls), rel_tol=1e-12):
            raise ValueError('Direct-write timing boundary mismatch')
        if not math.isclose(row['sum_cpu_seconds'], sum(call['cpu_seconds'] for call in calls), rel_tol=1e-12):
            raise ValueError('Direct-write CPU accounting mismatch')
        rates.append(amount/seconds)
    return dict(bytes_per_second=statistics.median(rates), workers=workers,
        repeat_bytes_per_second=rates, filesystem=filesystem,
        service_kind='independent_direct_write_aggregate',
        unpriced_terms=['Transfer of direct service to buffered kernel writeback, extent allocation, fragmentation, '
                        'mixed reads/writes, device/controller cache and dirty-page throttling'],
        scope='Fixed 16MiB writes,1GiB total, preallocated private files; verified O_DIRECT page bypass, '
              'process write-byte accounting and full payload readback. Wall service at the requested concurrency; '
              'final fsync and buffered CPU copying are separate. No association times or fitted multiplier.')


def with_direct_write_capacity(profile, probe, *, cpu_affinity, filesystem):
    """Return a profile with per-stream service and four-stream shared capacity.

    CPU page-cache and clean-syscall prices are preserved verbatim. Independent
    one-worker and aggregate four-worker prices have different roles; the
    shared output resource prevents streams each receiving aggregate capacity.
    """
    import copy
    single = direct_write_capacity(probe, workers=1, cpu_affinity=cpu_affinity, filesystem=filesystem)
    aggregate = direct_write_capacity(probe, workers=4, cpu_affinity=cpu_affinity, filesystem=filesystem)
    result = copy.deepcopy(profile)
    result['writeback_service']['storage_seconds_per_byte'] = 1/single['bytes_per_second']
    result['write_bytes_per_second'] = aggregate['bytes_per_second']
    result['output_storage_provenance'] = dict(single_stream=single, shared=aggregate)
    return result


def staging_copy_prices(probe, *, cpu_affinity, python_version, numpy_version, mechanism='numpy_copyto'):
    """Price view dispatch plus one bulk copy from fixed generic controls.

    The empty-copy CPU control prices each call; the independent 64MiB control
    supplies bulk CPU service after subtracting that one call. The 4KiB and
    two-worker observations remain diagnostics, not shape interpolation.
    """
    if (probe.get('affinity') != list(cpu_affinity) or probe.get('python_version') != python_version
            or probe.get('numpy_version') != numpy_version):
        raise ValueError('Staging-copy runtime context mismatch')
    if mechanism not in ('numpy_copyto', 'memoryview_slice'):
        raise ValueError('Unmodeled staging-copy mechanism')
    service = {}
    for size in (0, 64 << 20):
        rows = [r for r in probe['rows'] if r.get('mechanism') == mechanism
                and r.get('workers') == 1 and r.get('payload_bytes') == size]
        if len(rows) != 5 or {r['repeat'] for r in rows} != set(range(5)):
            raise ValueError('Incomplete generic staging-copy controls')
        expected_loops = 10000 if not size else 16
        values = []
        for row in rows:
            if row.get('loops') != expected_loops or len(row['threads']) != 1 or row['threads'][0]['worker'] != 0:
                raise ValueError('Staging-copy control geometry differs')
            cpu = positive('staging CPU', row['cpu_seconds'])
            positive('staging elapsed', row['seconds'])
            if not math.isclose(cpu, row['threads'][0]['cpu_seconds'], rel_tol=1e-12):
                raise ValueError('Staging-copy CPU accounting mismatch')
            values.append(cpu/expected_loops)
        service[size] = statistics.median(values)
    bulk = positive('staging bulk CPU', service[64 << 20]-service[0])/(64 << 20)
    return dict(cpu_seconds_per_call=service[0], cpu_seconds_per_byte=bulk, mechanism=mechanism,
        scope='Independent one-worker empty and64MiB CPU controls for '+mechanism+'. '
              'Fixed bulk service is a DRAM-sized scenario; cache residency, view keyword dispatch and '
              'loaded multi-caller ordering remain approximate. No scan times, shape interpolation or fitted residual.')


def writer_copy_cost(profile):
    """Validate an explicit copy service or retain the archived byte-only input."""
    price = profile.get('writer_copy_service')
    if price is None:
        price = dict(cpu_seconds_per_byte=profile['process_units']['numpy_copy_bytes'], cpu_seconds_per_call=0.)
    return {key: positive(key, price[key], True) for key in ['cpu_seconds_per_byte', 'cpu_seconds_per_call']}
