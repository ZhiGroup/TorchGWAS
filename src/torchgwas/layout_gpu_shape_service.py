"""Compact, source-bound GPU shape service for a fixed unissued PGEN layout.

Only distinct (device, chunk rows, trait width) tensor shapes are priced. A
multi-million-record schedule therefore does not create one Python GPU graph
per chunk. These are conditional component services, not completion bounds.
"""

import math
import time

from .analytical_plan_cache import input_identity
from .calibration_cache import _digest
from .planning_session import planning_work_scope
from .tensor_service import DeviceService


def _positive(value, name):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value <= 0):
        raise ValueError('Positive finite ' + name + ' required')
    return value


def _geometry(rows, key, name):
    if not isinstance(rows, list):
        raise ValueError('Explicit ' + name + ' geometry required')
    matching = [row for row in rows if isinstance(row, dict) and key(row)]
    if len(matching) != 1:
        raise ValueError('Exactly one matching ' + name + ' geometry required')
    return matching[0]


def _price_shape_groups(groups, profiles, samples, covariate_rank, mode):
    from .mechanistic_torch import _shape_component

    started = time.perf_counter()
    cpu_started = time.thread_time()
    devices = {device for device, _, _ in groups}
    if not isinstance(profiles, dict) or set(profiles) != devices:
        raise ValueError('One shape profile per active device required')
    priced = []
    per_device = {device: dict(chunks=0, kernel_service_seconds=0.,
                               host_dispatch_cpu_seconds=0.,
                               host_dispatch_serial_cpu_seconds=0.)
                  for device in devices}
    with planning_work_scope():
        for (device, size, width), count in sorted(groups.items()):
            profile = profiles[device]
            scan_modes = ({None, 'device_significant'} if mode == 'significant'
                          else {mode})
            if (not isinstance(profile, dict) or
                    profile.get('reduction') not in scan_modes or
                    not isinstance(profile.get('gpu_resources'), dict) or
                    not isinstance(profile.get('host_primitives'), dict)):
                raise ValueError('Matching independent tensor profile required')
            q = _positive(profile.get('cpu_fraction'), 'host CPU fraction')
            if q > 1:
                raise ValueError('Host CPU fraction exceeds one')
            statistics = _geometry(profile.get('kernel_geometry'),
                lambda row: (row.get('N'), row.get('B'), row.get('K', 1),
                             row.get('C', 8), row.get('validate_range', False)) ==
                    (samples, size, width, covariate_rank,
                     profile.get('validate_range', False)), 'statistics')
            joint = None
            if mode == 'jagwas':
                if not isinstance(profile.get('joint_host_primitives'), dict):
                    raise ValueError('Independent JAGWAS host primitives required')
                joint = _geometry(profile.get('joint_kernel_geometry'),
                    lambda row: (row.get('N'), row.get('B'), row.get('K'),
                                 row.get('compute_dtype')) ==
                        (samples, size, width, 'float32'), 'JAGWAS projection')
            gpu = DeviceService(**dict(profile['gpu_resources'],
                                       host_cpu_fraction=q))
            component = _shape_component(samples, size, width,
                                         covariate_rank, profile, gpu,
                                         statistics, joint)
            if component.get('status') == 'zero_available_capacity':
                raise ValueError('Zero available GPU tensor capacity')
            values = [component.get('kernel_service_seconds'),
                      component.get('host_dispatch_cpu_seconds'),
                      component.get('host_dispatch_serial_cpu_seconds', 0.)]
            if any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                   not math.isfinite(value) or value < 0 for value in values):
                raise ValueError('Finite nonnegative tensor shape service required')
            device_work = per_device[device]
            device_work['chunks'] += count
            for key, value in zip(('kernel_service_seconds',
                                   'host_dispatch_cpu_seconds',
                                   'host_dispatch_serial_cpu_seconds'), values):
                device_work[key] += count * value
            priced.append(dict(device=device, markers=size, traits=width,
                               chunks=count, kernel_service_seconds=values[0],
                               host_dispatch_cpu_seconds=values[1],
                               host_dispatch_serial_cpu_seconds=values[2],
                               component_source_sha256=component.get('source_sha256'),
                               unpriced_terms=component.get('unpriced_terms', [])))
    totals = dict(host_dispatch_cpu_seconds=math.fsum(
                      row['host_dispatch_cpu_seconds'] for row in per_device.values()),
                  host_dispatch_serial_cpu_seconds=math.fsum(
                      row['host_dispatch_serial_cpu_seconds'] for row in per_device.values()))
    return dict(distinct_shapes=len(groups), shapes=priced,
                per_device_service=per_device, total_host_work=totals,
                profile_sha256={device: _digest(profiles[device])
                                for device in sorted(devices)},
                calculation_wall_seconds=time.perf_counter() - started,
                calculation_cpu_seconds=time.thread_time() - cpu_started)


def native_layout_gpu_shape_service(layout, *, covariate_rank, profiles,
                                    max_shapes=64):
    """Sum exact full/tail tensor shapes using the scan calculator's prices.

    Profiles contain independently measured device resources, host primitive
    prices and duration-free compiled geometry. Pricing starts after useful
    output in the productive planner; it does no cold, full-job graph build.
    """
    if (not isinstance(layout, dict) or
            layout.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            type(covariate_rank) is not int or
            type(max_shapes) is not int or max_shapes < 1):
        raise ValueError('Bound native source layout and shape budget required')
    samples = layout.get('samples')
    if (type(samples) is not int or samples < 3 or
            not 0 <= covariate_rank < samples - 2):
        raise ValueError('Bound sample count and residual covariate rank required')
    source = layout.get('input_identity')
    if (not isinstance(source, dict) or 'path' not in source or
            input_identity(source['path']) != source):
        raise ValueError('PGEN input changed before GPU shape pricing')
    rows = layout.get('partitions')
    chunk = layout.get('chunk_markers')
    mode = layout.get('reduction')
    if (not isinstance(rows, list) or not rows or
            type(chunk) is not int or chunk < 1 or
            mode not in (None, 'significant', 'jagwas')):
        raise ValueError('Fixed nonempty native source layout required')
    if any(not isinstance(row, dict) for row in rows):
        raise ValueError('Complete source partitions required')
    devices = {row.get('device') for row in rows}
    if (len(devices) == 0 or any(not isinstance(device, str) for device in devices) or
            not isinstance(profiles, dict) or set(profiles) != devices):
        raise ValueError('One shape profile per active device required')
    groups = {}
    partitions = []
    for row in rows:
        if (not isinstance(row, dict) or
                set(row) != {'id', 'device', 'variant_range', 'trait_range', 'chunks'}):
            raise ValueError('Complete source partition required')
        variant, trait = row['variant_range'], row['trait_range']
        if (not isinstance(variant, list) or len(variant) != 2 or
                not isinstance(trait, list) or len(trait) != 2 or
                any(type(value) is not int for value in variant + trait) or
                not 0 <= variant[0] < variant[1] or
                not 0 <= trait[0] < trait[1]):
            raise ValueError('Nonempty source shape ranges required')
        markers = variant[1] - variant[0]
        width = trait[1] - trait[0]
        full, tail = divmod(markers, chunk)
        counts = {chunk: full} if full else {}
        if tail:
            counts[tail] = counts.get(tail, 0) + 1
        if (row['chunks'] != full + bool(tail) or
                (mode == 'jagwas' and trait != [0, layout['total_traits']])):
            raise ValueError('Source chunk or JAGWAS panel geometry changed')
        for size, count in counts.items():
            key = (row['device'], size, width)
            groups[key] = groups.get(key, 0) + count
        partitions.append(dict(id=row['id'], device=row['device'],
                               markers=markers, traits=width,
                               chunks=full + bool(tail),
                               chunk_shape_counts=dict(counts)))
    if len(groups) > max_shapes:
        raise ValueError('Distinct GPU shape budget exceeded')
    service = _price_shape_groups(groups, profiles, samples,
                                  covariate_rank, mode)
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed during GPU shape pricing')
    return dict(kind='torchgwas.pgen_layout_gpu_shape_service.v1',
                input_identity=dict(source), reduction=mode,
                samples=samples, covariate_rank=covariate_rank,
                chunk_markers=chunk, partitions=partitions,
                **service,
                prediction_complete=False, selection_validated=False,
                scope='Conditional GPU kernel and host dispatch service from independently priced exact full/tail tensor shapes. Per-device kernels are serialized and shared host work is counted once. This is not a necessary capacity lower bound or an elapsed completion ceiling; transfer, output, queues, live load and final drain are separate.')


def productive_issued_gpu_shape_service(issued, frontier, *, profiles,
                                        max_shapes=64):
    """Conservatively price full issued chunks on their original devices."""
    from .layout_frontier import KIND as FRONTIER_KIND

    if (not isinstance(issued, dict) or
            issued.get('kind') != 'torchgwas.productive_issued_work.v1' or
            not isinstance(frontier, dict) or frontier.get('kind') != FRONTIER_KIND or
            any(issued.get(key) != frontier.get(key) for key in
                ('input_identity', 'issued_revision', 'written_events',
                 'reduction')) or
            type(max_shapes) is not int or max_shapes < 1 or
            not isinstance(issued.get('chunks'), list) or
            issued.get('pending_chunks') != len(issued['chunks'])):
        raise ValueError('Matching bound issued chunks and frontier required')
    source = issued['input_identity']
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed before issued GPU shape pricing')
    samples = issued['samples']
    rank = issued['covariate_rank']
    mode = issued['reduction']
    groups = {}
    records = 0
    for row in issued['chunks']:
        variant, trait = row['variant_range'], row['trait_range']
        if (not isinstance(variant, list) or len(variant) != 2 or
                not isinstance(trait, list) or len(trait) != 2 or
                any(type(value) is not int for value in variant + trait) or
                not 0 <= variant[0] < variant[1] or
                not 0 <= trait[0] < trait[1] <= frontier['total_traits'] or
                (mode == 'jagwas' and trait != [0, frontier['total_traits']])):
            raise ValueError('Bound issued full-chunk geometry required')
        size = variant[1] - variant[0]
        width = trait[1] - trait[0]
        records += size
        key = (row['device'], size, width)
        groups[key] = groups.get(key, 0) + 1
    if (records != issued['pending_source_records'] or
            len(groups) > max_shapes):
        raise ValueError('Issued GPU shape budget or source count differs')
    service = _price_shape_groups(groups, profiles, samples, rank, mode)
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed during issued GPU shape pricing')
    return dict(kind='torchgwas.productive_issued_gpu_shape_service.v1',
                input_identity=dict(source), reduction=mode,
                issued_revision=issued['issued_revision'],
                written_events=issued['written_events'],
                samples=samples, covariate_rank=rank,
                pending_chunks=issued['pending_chunks'],
                pending_source_records=records,
                **service,
                prediction_complete=False, selection_validated=False,
                scope='Conditional full-chunk GPU service for issued work without output completion, on original devices. Some or all service may already be in flight or complete. This is nominal replay work, not remaining time or an elapsed completion ceiling.')
