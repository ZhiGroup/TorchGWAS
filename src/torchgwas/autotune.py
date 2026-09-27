"""Bounded joint geometry/device search using the shared pipeline accounting.

Detailed scan selection uses the calculator's finite execution graph through
``detailed_joint_plan``. The older ``joint_plan`` remains a coarse advisory
resource report for other paths. Neither fits association timings.
"""
from __future__ import annotations

from dataclasses import replace
import itertools
import math

from .pipeline_model import Workload, InputProfile, Hardware, PipelinePlan, estimate
from .linear import multigpu_variant_ranges
from .pinned_work import pinned_scan_work
from .setup_work import setup_memory
from .tensor_memory import eager_scan_memory,eager_memory_plan
from .mechanistic_plan import detailed_joint_plan
from .trait_tiling_plan import detailed_trait_plan
from .trait_candidate_space import prepare_trait_candidates, bounded_trait_plan


def _positive_int(name, value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return value


def joint_plan(workload: Workload, input_profile: InputProfile, hardware: Hardware,
               *, devices=('cuda:0',), device_hardware=None,
               chunks=(256, 1024, 2048, 4096), trait_blocks=None, device_memory_profiles=None,
               workers=(1, 2, 4), depths=(2, 4), device_sets=None,
               modes=('variant', 'trait'), reduction=None,
               host_reserve_bytes=0, device_reserve_bytes=0,
               result_queue_depth=4, shared_links=(),
               max_candidates=10000, shortlist_size=5, tie_fraction=0.0):
    """Minimize a conditional resource score over a finite feasible set.

    workers is PER DEVICE; the sum cannot exceed hardware.cpu_workers.
    hardware storage/host-memory/CPU fields describe SHARED resources.
    device_hardware can override per-device GPU/link/memory capacities.
    shared_links contains {devices, h2d_bytes_per_second, d2h_bytes_per_second}.
    Multiple groups can overlap (e.g. a PCIe switch and its upstream root).
    Reserves cover live phenotype arrays, allocator and decoder workspace,
    retained reductions and writer state absent from the ring accounting.

    Trait sharding maps to run_linear_gwas(reduce='significant'), which rereads
    all variants per block. Full-output variant sharding maps to
    linear_scan_multigpu(ordered=False). No hybrid executor is assumed.
    """
    devices = tuple(devices)
    if not devices or len(set(devices)) != len(devices):
        raise ValueError('devices must be nonempty and unique')
    if any(not isinstance(d, str) or not d.startswith('cuda:') or not d[5:].isdigit() for d in devices):
        raise ValueError('explicit CUDA device names required')
    for name, value in [('max_candidates', max_candidates), ('shortlist_size', shortlist_size),
                        ('result_queue_depth', result_queue_depth)]:
        _positive_int(name, value)
    for name, value in [('host_reserve_bytes', host_reserve_bytes),
                        ('device_reserve_bytes', device_reserve_bytes), ('tie_fraction', tie_fraction)]:
        if not math.isfinite(value) or value < 0:
            raise ValueError(f'{name} must be finite and nonnegative')
    if reduction not in (None, 'significant'):
        raise ValueError('joint planning supports full or significant output; JAGWAS is not trait-separable')
    if not set(modes) <= {'variant', 'trait'} or not modes:
        raise ValueError('modes must contain variant and/or trait')
    if input_profile.cpu_decode_core_seconds_total is not None or input_profile.genotype_record_read_bytes_total is not None:
        raise ValueError('exact decode-work census must be recomputed per shard/geometry; aggregate overrides unsupported')
    per_device = {d: (device_hardware or {}).get(d, hardware) for d in devices}
    if any(h.device_count != 1 for h in per_device.values()):
        raise ValueError('supply one Hardware record per device')
    links = []
    for group in shared_links:
        members = tuple(group['devices'])
        if not members or len(set(members)) != len(members) or not set(members) <= set(devices):
            raise ValueError('shared link members must be unique known devices')
        for direction in ('h2d', 'd2h'):
            rate = group[direction + '_bytes_per_second']
            if not math.isfinite(rate) or rate <= 0:
                raise ValueError('shared link rates must be finite and positive')
        links.append(group)
    if device_sets is None:
        device_sets = [group for count in range(1, len(devices) + 1)
                       for group in itertools.combinations(devices, count)]
    device_sets = [tuple(group) for group in device_sets]
    for group in device_sets:
        if not group or len(set(group)) != len(group) or not set(group) <= set(devices):
            raise ValueError('device_sets must contain nonempty subsets of devices')
    axes = []
    for name, axis in [('chunks', chunks), ('workers', workers), ('depths', depths)]:
        axis = sorted(set(axis))
        if not axis:
            raise ValueError(name + ' cannot be empty')
        for value in axis:
            _positive_int(name, value)
        axes.append(axis)
    chunks, workers, depths = axes
    if min(depths) < 2:
        raise ValueError('depths must be at least two')
    if trait_blocks is None:
        trait_blocks = sorted({workload.traits} | {
            min(workload.traits, 64 * 2**i) for i in range(max(1, workload.traits.bit_length() - 5))})
    trait_blocks = sorted(set(trait_blocks))
    for width in trait_blocks:
        _positive_int('trait_blocks', width)
        if width > workload.traits:
            raise ValueError('trait block cannot exceed phenotype width')
    geometries = []
    if 'variant' in modes and reduction is None:
        geometries.append(('variant', workload.traits))
    if 'trait' in modes and reduction == 'significant':
        geometries.extend(('trait', width) for width in trait_blocks)
    if not geometries:
        raise ValueError('no executable mode for requested reduction')
    count = len(geometries) * len(device_sets) * len(chunks) * len(workers) * len(depths)
    if count > max_candidates:
        raise ValueError(f'{count} candidates exceeds max_candidates={max_candidates}; narrow explicit bounds')

    ranked, rejected, unknown = [], {}, set()
    for (mode, width), group, chunk, worker, depth in itertools.product(
            geometries, device_sets, chunks, workers, depths):
        reason = None
        if worker * len(group) > hardware.cpu_workers:
            reason = 'shared_cpu_budget'
        if mode == 'variant':
            ranges = multigpu_variant_ranges(workload.variants, chunk, len(group))
            tasks = [(i, end-start, width, (start, end)) for i, (start, end) in enumerate(ranges)]
        else:
            tasks = [(i % len(group), workload.variants, min(width, workload.traits-offset),
                      (offset, min(offset+width, workload.traits)))
                     for i, offset in enumerate(range(0, workload.traits, width))]
        if len({task[0] for task in tasks}) != len(group):
            reason = 'idle_devices'
        if reason:
            rejected[reason] = rejected.get(reason, 0) + 1
            continue
        shared = dict(storage=0., cpu_decode=0., host_memory=0.)
        device_cost = {d: dict(gpu=0., h2d=0., d2h=0., cpu_decode=0.) for d in group}
        transfer = {d: dict(h2d=0., d2h=0.) for d in group}
        host_peak = dict.fromkeys(group, 0)
        gpu_peak = dict.fromkeys(group, 0)
        overhead = dict.fromkeys(group, 0.)
        assignments = []
        try:
            for index, markers, traits, span in tasks:
                d = group[index]
                local_h = replace(per_device[d], cpu_workers=hardware.cpu_workers,
                    disk_bytes_per_second=hardware.disk_bytes_per_second,
                    output_bytes_per_second=hardware.output_bytes_per_second,
                    host_bytes_per_second=hardware.host_bytes_per_second,
                    host_memory_bytes=hardware.host_memory_bytes,
                    device_memory_bytes=max(0, per_device[d].device_memory_bytes-device_reserve_bytes))
                local_w = replace(workload, variants=markers, traits=traits)
                local_p = replace(input_profile, stored_bytes=math.ceil(input_profile.stored_bytes * markers / workload.variants))
                cost = estimate(local_w, local_p, local_h,
                    PipelinePlan(chunk, chunk, chunk, worker, depth))
                for resource in shared:
                    if resource == 'cpu_decode':
                        # Sum core-seconds, then divide by allocated global cores.
                        shared[resource] += cost['resource_seconds'][resource] * min(worker, depth)
                    else:
                        shared[resource] += cost['resource_seconds'][resource]
                for resource in device_cost[d]:
                    device_cost[d][resource] += cost['resource_seconds'][resource]
                # The design is resident for a task and uploaded once per task.
                design_bytes = 4 * workload.samples * (traits + workload.covariates + 1)
                design_seconds = design_bytes / local_h.h2d_bytes_per_second
                device_cost[d]['h2d'] += design_seconds
                overhead[d] += design_seconds
                shared['host_memory'] += design_bytes / hardware.host_bytes_per_second
                transfer[d]['h2d'] += cost['input_transfer_bytes'] + design_bytes
                transfer[d]['d2h'] += cost['result_d2h_bytes']
                # Independent source scans run consecutively on each GPU.
                overhead[d] += cost['fill_drain_seconds']
                # Queue plus producer-held item, with full-width result upper estimate.
                queue_bytes = (max(result_queue_depth, 4)+1) * chunk * (8*traits+5)
                host_buffers = cost['host_buffer_bytes']
                if mode == 'variant' and input_profile.direct_native_fill and not input_profile.decode_on_gpu:
                    pins = pinned_scan_work(workload.samples, chunk, traits, depth,
                        transfer_bytes_per_variant=math.ceil(input_profile.transfer_bytes_per_variant))
                    host_buffers += pins['allocator_bytes'] - pins['requested_bytes']
                    unknown.add('pinned allocator cache history and driver backing beyond rounded requests')
                host_peak[d] = max(host_peak[d], host_buffers + queue_bytes)
                setup = setup_memory(workload.samples, traits, workload.covariates)
                full_setup = setup_memory(workload.samples, workload.traits, workload.covariates)
                unknown.update(setup['unresolved_memory_terms'])
                gpu_peak[d] = max(gpu_peak[d], cost['device_buffer_bytes'],
                    setup['design_device_live_bytes_upper'],
                    full_setup['residual_device_live_bytes_upper'])
                if input_profile.direct_native_fill and not input_profile.decode_on_gpu and input_profile.transfer_bytes_per_variant == workload.samples:
                    eager = eager_scan_memory(workload.samples,chunk,traits,workload.covariates,depth)
                    gpu_peak[d] = max(gpu_peak[d],eager['tensor_storage_budget'])
                    if device_memory_profiles and d in device_memory_profiles:
                        memory = eager_memory_plan(workload.samples,chunk,traits,workload.covariates,depth,
                            device_memory_profiles[d],preprocessing_traits=workload.traits)
                        gpu_peak[d] = max(gpu_peak[d],memory['device_bytes'])
                        unknown.update(memory['unresolved_memory_terms'])
                unknown.add('device allocator reservations and library workspace beyond tensor storage accounting')
                if gpu_peak[d]+device_reserve_bytes > per_device[d].device_memory_bytes:
                    reason = 'device_memory'
                assignments.append(dict(device=d, axis=mode, span=list(span), traits=traits, variants=markers))
                unknown.update(cost['unresolved_timing_terms'])
        except ValueError as error:
            # Shape-specific service/census must never silently extrapolate.
            reason = 'unsupported_component: ' + str(error)
        if reason:
            rejected[reason] = rejected.get(reason, 0) + 1
            continue
        host_bytes = sum(host_peak.values()) + host_reserve_bytes
        if host_bytes > hardware.host_memory_bytes:
            rejected['host_memory'] = rejected.get('host_memory', 0) + 1
            continue
        shared['cpu_decode'] /= len(group) * min(worker, depth)
        resource = dict(shared)
        for d, costs in device_cost.items():
            for key, value in costs.items():
                resource[d + ':' + key] = value
        for index, link in enumerate(links):
            for direction in ('h2d', 'd2h'):
                resource[f'link{index}:{direction}'] = sum(
                    transfer[d][direction] for d in group if d in link['devices']) / link[direction+'_bytes_per_second']
        bound = max(resource.values())
        # A coarse finite-batch correction, NOT an execution-graph prediction.
        score = bound + max(overhead.values())
        common = dict(chunk_size=chunk, reader_workers=worker, prefetch_chunks=depth)
        if mode == 'variant':
            entrypoint = 'torchgwas.linear.linear_scan_multigpu'
            kwargs = dict(common, devices=list(group), reader_workers=worker*len(group),
                          ordered=False, result_queue_depth=result_queue_depth, compute_p_values=False)
        else:
            entrypoint = 'torchgwas.api.run_linear_gwas'
            kwargs = dict(common, device=group[0], trait_devices=list(group),
                          trait_block=width, reduce='significant')
        ranked.append(dict(mode=mode, devices=list(group), chunk_size=chunk, trait_block=width,
            workers_per_device=worker, prefetch_chunks=depth, resource_score_seconds=score,
            resource_floor_seconds=bound, resource_seconds=resource,
            bottleneck=max(resource, key=resource.get), host_bytes=host_bytes,
            device_bytes={d: gpu_peak[d]+device_reserve_bytes for d in group},
            genotype_passes=1 if mode == 'variant' else math.ceil(workload.traits/width),
            assignments=assignments, entrypoint=entrypoint, scan_kwargs=kwargs,
            source_controls=dict(read_variants=chunk, decode_variants=chunk)))
    if not ranked:
        raise ValueError('no feasible supported candidate: ' + str(rejected))
    ranked.sort(key=lambda row: row['resource_score_seconds'])
    best_score = ranked[0]['resource_score_seconds']
    tied = [row for row in ranked if row['resource_score_seconds'] <= best_score*(1+tie_fraction)]
    chosen = min(tied, key=lambda row: (len(row['devices']), row['host_bytes']+sum(row['device_bytes'].values()),
                                      row['workers_per_device'], row['resource_score_seconds']))
    shortlist = [chosen] + [row for row in ranked if row is not chosen][:shortlist_size-1]
    return dict(status='advisory_resource_plan', workload=dict(variants=workload.variants, samples=workload.samples, traits=workload.traits, covariates=workload.covariates), selected=chosen, shortlist=shortlist,
        candidates_evaluated=count, candidates_feasible=len(ranked), rejected=rejected,
        predicted_runtime_seconds=None, runtime_prediction_validated=False,
        unresolved_timing_terms=sorted(unknown | {
            'chunk/trait shape efficiency beyond supplied component coverage',
            'preprocessing, per-block setup, exact significance/reduction and writer service',
            'finite output queues, host scheduling and device placement contention',
            'source controls must be checked against the actual adapter before execution'}),
        assumptions=[
            'Score reuses pipeline_model resource accounting; it is not the newer detailed K1 runtime model.',
            'GPU logical-traffic/roofline estimates are incomplete; no optimal measured runtime claim.',
            'All repeated genotype passes charged to the supplied input bandwidth; cache reuse is not assumed.',
            'Variable-record bytes and decode demand distributed uniformly across variant shards.',
            'Per-device link capacities require explicit shared_links for common upstream links.',
            'Memory is an estimate: reserves must cover decoder, phenotype, writer and retained reduction state.',
            'Trait significance D2H is conservatively priced as full-width results; output density remains unknown.',
            'Near-tie tolerance is a selection policy, not measured uncertainty.',
            'Completion-order variant output requires a consumer that uses absolute coordinates.'])


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profile', help='JSON workload/input/hardware plus joint search bounds')
    parser.add_argument('--output')
    args = parser.parse_args()
    with open(args.profile, encoding='utf-8-sig') as handle:
        data = json.load(handle)
    options = data.get('joint', {})
    if data.get('model') == 'detailed_traits_space':
        if set(data) != {'model', 'workload', 'contexts', 'bounds', 'joint', 'output'}:
            raise ValueError('Detailed trait space needs workload, contexts, bounds, joint and output')
        result = bounded_trait_plan(data['workload'], data['contexts'], bounds=data['bounds'],
                                    joint=options, output=data['output'])
    elif data.get('model') == 'detailed_significant_host_space':
        from .significant_candidate_space import bounded_significant_host_plan
        if set(data)-{'significance_threshold'} != {'model','workload','contexts','bounds','joint','output','prices'}:
            raise ValueError('Detailed significant space needs workload, contexts, bounds, joint, output, prices and optional significance_threshold')
        result = bounded_significant_host_plan(**{key:value for key,value in data.items() if key!='model'})
    elif data.get('model') == 'detailed_jagwas_space':
        from .jagwas_candidate_space import bounded_jagwas_plan
        if set(data) != {'model','workload','contexts','bounds','joint','output','prices','preparation_services'}:
            raise ValueError('Detailed JAGWAS space needs workload, contexts, bounds, joint, output, prices and preparation_services')
        result = bounded_jagwas_plan(**{key:value for key,value in data.items() if key!='model'})
    elif data.get('model') in ('detailed_scan','detailed_traits'):
        if set(data) - {'model', 'candidates', 'joint'}:
            raise ValueError('Detailed scan input accepts model, candidates and joint only')
        planner = detailed_joint_plan if data['model']=='detailed_scan' else detailed_trait_plan
        result = planner(data['candidates'], **options)
    else:
        if data.get('model') not in (None, 'coarse_resource'):
            raise ValueError('Unknown planning model')
        if 'device_hardware' in options:
            options['device_hardware'] = {key: Hardware(**value) for key, value in options['device_hardware'].items()}
        result = joint_plan(Workload(**data['workload']), InputProfile(**data['input']),
                            Hardware(**data['hardware']), **options)
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()

