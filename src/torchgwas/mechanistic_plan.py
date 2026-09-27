"""Bounded selection using the same finite execution graph as the calculator.

Inputs contain source censuses, untimed kernel geometry and independent prices.
They never contain measured scan durations. This first execution contract is a
complete-phenotype native-int8 PGEN scan with borrowed, immediately consumed
results. Setup and durable output are outside the timing objective.
"""
from __future__ import annotations

import math

from .gpu_identity import canonical_cuda_device as _device
from .linear import multigpu_variant_ranges
from .mechanistic_torch import torch_multigpu_scan_runtime
from .pinned_work import pinned_scan_work
from .setup_work import setup_memory
from .tensor_memory import eager_scan_memory, eager_memory_plan


def _integer(name, value, minimum=1):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def _shape(candidate):
    """Check the candidate maps exactly to the public variant-shard executor."""
    allowed = {'shards', 'shared_capacities', 'shared_links', 'queue_service', 'result_queue_depth'}
    if set(candidate) - allowed:
        raise ValueError('Unknown candidate fields: ' + str(sorted(set(candidate) - allowed)))
    shards = candidate['shards']
    if not isinstance(shards, list) or not shards:
        raise ValueError('Nonempty shard list required')
    devices = [_device(s['device']) for s in shards]
    if len(set(devices)) != len(devices):
        raise ValueError('Duplicate device')
    first = shards[0]
    n, k, c = (first['data'][key] for key in ('samples', 'traits_analyzed', 'covariates'))
    for name, value, minimum in [('samples', n, 32), ('traits', k, 1), ('covariates', c, 0)]:
        _integer(name, value, minimum)
    b, depth = (first['profile'][key] for key in ('chunk_markers', 'depth'))
    _integer('chunk_markers', b)
    _integer('depth', depth, 2)
    m = sum(_integer('markers', s['data']['markers']) for s in shards)
    ranges = multigpu_variant_ranges(m, b, len(shards))
    if len(ranges) != len(shards):
        raise ValueError('Idle device shard')
    capacities = candidate['shared_capacities']
    if set(capacities) != {'cpu', 'dram', 'input'}:
        raise ValueError('Declare aggregate cpu, dram and input capacities')
    for key, value in capacities.items():
        if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise ValueError('Invalid shared capacity: ' + key)
    workers = []
    event_wait = first['profile'].get('event_wait_cpu_fraction')
    if isinstance(event_wait, bool) or event_wait not in (0., 1.):
        raise ValueError('Executable selection requires explicit spin or blocking completion events')
    for shard, span in zip(shards, ranges):
        data, profile = shard['data'], shard['profile']
        if (data['samples'], data['traits_analyzed'], data['covariates']) != (n, k, c):
            raise ValueError('Shard phenotype dimensions differ')
        if data.get('matching_sample_order') is not True:
            raise ValueError('Matching full sample order required')
        if (profile['chunk_markers'], profile['depth']) != (b, depth):
            raise ValueError('Executor requires a common chunk and prefetch depth')
        if profile.get('result_ownership') != 'borrowed':
            raise ValueError('Detailed selection currently requires borrowed results')
        if profile.get('validate_range') is not True:
            raise ValueError('Public native scan requires genotype range validation')
        if profile.get('event_wait_cpu_fraction') != event_wait:
            raise ValueError('Completion-event policy must agree across shards')
        encoded = data['encoded']
        if encoded.get('variant_range') != list(span) or encoded.get('file_markers') != m:
            raise ValueError('Census ranges must match full-file executor partition')
        if encoded.get('path') != first['data']['encoded'].get('path') or not encoded.get('path'):
            raise ValueError('All shards must census the same input path')
        if not encoded.get('chunks') or encoded.get('chunk_markers') != b:
            raise ValueError('Exact per-chunk encoded census required')
        for resource, field in [('cpu', 'cpu_available_cores'), ('dram', 'shared_dram_bytes_per_second'),
                                ('input', 'read_bytes_per_second')]:
            if profile[field] != capacities[resource]:
                raise ValueError('Shared capacity differs from shard profile: ' + resource)
        workers.append(_integer('decode_workers', profile['decode_workers']))
    total_workers = sum(workers)
    per_device, remainder = divmod(total_workers, len(shards))
    if workers != [per_device + (i < remainder) for i in range(len(shards))]:
        raise ValueError('Reader allocation differs from public executor')
    queue_depth = _integer('result_queue_depth', candidate.get('result_queue_depth', 4))
    if len(shards) > 1 and candidate.get('queue_service') is None:
        raise ValueError('Independent multiGPU queue service required')
    return dict(samples=n, markers=m, traits=k, covariates=c, chunk_size=b,
                depth=depth, workers=workers, devices=devices, reader_workers=total_workers,
                result_queue_depth=queue_depth, input_path=first['data']['encoded']['path'],
                blocking_events=event_wait == 0.)


def _memory(candidate, shape, host_reserve, device_reserve, memory_profiles):
    from .decoder_work import native_reader_workspace
    n, k, c, b, depth = (shape[key] for key in ('samples', 'traits', 'covariates', 'chunk_size', 'depth'))
    pins = pinned_scan_work(n, b, k, depth)
    setup = setup_memory(n, k, c)
    eager = eager_scan_memory(n, b, k, c, depth)
    device_bytes, host_shards = {}, []
    unknown = set(setup['unresolved_memory_terms'])
    for shard in candidate['shards']:
        device, encoded = shard['device'], shard['data']['encoded']
        active = min(shard['profile']['decode_workers'], depth, math.ceil(shard['data']['markers']/b))
        # Both decoder arrays are grow-only. Charge the largest encoded chunk
        # to every reader, plus its packed workspace. Index/metadata Python
        # representations and construction peaks remain in explicit reserves.
        packed = active * min(b, shard['data']['markers']) * ((n+3)//4)
        payload,replay_workspace = native_reader_workspace(encoded['chunks'])
        records = active * payload
        host_shards.append(dict(device=device, pinned_bytes=pins['allocator_bytes'],
                                packed_decoder_bytes=packed, record_buffer_bytes=records))
        if replay_workspace:host_shards[-1]['ld_replay_workspace_bytes']=active*replay_workspace
        budget = max(eager['tensor_storage_budget'], setup['device_live_bytes_upper'])
        if device in memory_profiles:
            memory = eager_memory_plan(n, b, k, c, depth, memory_profiles[device])
            budget = max(budget, memory['device_bytes'])
            unknown.update(memory['unresolved_memory_terms'])
        else:
            unknown.add('CUDA library workspaces require a device memory profile or explicit reserve')
        device_bytes[device] = budget + device_reserve
    # Shared residual phenotype and covariate basis persist while workers scan.
    shared_host = 4*n*(k+c)
    host_bytes = shared_host + host_reserve + sum(
        row['pinned_bytes']+row['packed_decoder_bytes']+row['record_buffer_bytes']+row.get('ld_replay_workspace_bytes',0) for row in host_shards)
    unknown.update(['input arrays, index/metadata objects and their construction peaks require host reserve',
                    'allocator caches, driver allocations and unlisted library workspace require reserves'])
    return dict(host_bytes=host_bytes, device_bytes=device_bytes, host_shards=host_shards,
                shared_host_bytes=shared_host, host_reserve_bytes=host_reserve,
                device_reserve_bytes=device_reserve, unresolved_memory_terms=sorted(unknown))


def detailed_joint_plan(candidates, *, host_scenarios, cpu_workers, host_memory_bytes,
                        device_memory_bytes, host_reserve_bytes=0, device_reserve_bytes=0,
                        device_memory_profiles=None, max_candidates=1000,
                        max_scenario_evaluations=4000, shortlist_size=5,
                        max_slowdown_fraction=0.):
    """Minimize max_s T_graph(candidate, s) over explicitly bounded candidates.

    All scenarios must be supplied; their extrema are not confidence bounds.
    Memory constraints use source storage accounting plus caller reserves.
    max_slowdown_fraction=0 selects predicted throughput. A positive, explicit
    allowance instead prefers fewer devices within that relative slowdown.
    No scan is run, no profile is modified and no coarse score is consulted.
    """
    for name, value in [('cpu_workers', cpu_workers), ('host_memory_bytes', host_memory_bytes),
                        ('max_candidates', max_candidates), ('max_scenario_evaluations', max_scenario_evaluations),
                        ('shortlist_size', shortlist_size)]:
        _integer(name, value)
    for name, value in [('host_reserve_bytes', host_reserve_bytes), ('device_reserve_bytes', device_reserve_bytes)]:
        _integer(name, value, 0)
    if isinstance(max_slowdown_fraction, bool) or not math.isfinite(max_slowdown_fraction) or max_slowdown_fraction < 0:
        raise ValueError('max_slowdown_fraction must be finite and nonnegative')
    if not isinstance(candidates, (list, tuple)) or not candidates or len(candidates) > max_candidates:
        raise ValueError('Nonempty candidate list within max_candidates required')
    if not isinstance(host_scenarios, dict) or not host_scenarios:
        raise ValueError('Explicit named host scenarios required')
    if len(candidates)*len(host_scenarios) > max_scenario_evaluations:
        raise ValueError('Search exceeds max_scenario_evaluations')
    for name, scenario in host_scenarios.items():
        if not isinstance(name, str) or not name or set(scenario) != {'host_serial_fraction', 'host_serial_policy'}:
            raise ValueError('Each named scenario needs host_serial_fraction and host_serial_policy')
        fraction = scenario['host_serial_fraction']
        if isinstance(fraction, bool) or not math.isfinite(fraction) or not 0 <= fraction <= 1:
            raise ValueError('Invalid host_serial_fraction')
        if scenario['host_serial_policy'] not in ('fluid', 'held-first', 'held-last'):
            raise ValueError('Invalid host_serial_policy')
    for device, value in device_memory_bytes.items():
        _device(device)
        _integer('device_memory_bytes', value)
    shapes = [_shape(candidate) for candidate in candidates]
    identity_keys = ('samples', 'markers', 'traits', 'covariates', 'input_path')
    identity = {key: shapes[0][key] for key in identity_keys}
    if any(any(shape[key] != identity[key] for key in identity_keys) for shape in shapes):
        raise ValueError('Candidates must represent the same complete workload')
    ranked, rejected, unknown = [], [], set()
    for index, (candidate, shape) in enumerate(zip(candidates, shapes)):
        if not set(shape['devices']) <= set(device_memory_bytes):
            raise ValueError('Missing device memory budget')
        if shape['reader_workers'] > cpu_workers:
            rejected.append(dict(candidate_index=index, reason='shared_reader_budget'))
            continue
        memory = _memory(candidate, shape, host_reserve_bytes, device_reserve_bytes, device_memory_profiles or {})
        if memory['host_bytes'] > host_memory_bytes:
            rejected.append(dict(candidate_index=index, reason='host_memory', memory=memory))
            continue
        if any(memory['device_bytes'][d] > device_memory_bytes[d] for d in shape['devices']):
            rejected.append(dict(candidate_index=index, reason='device_memory', memory=memory))
            continue
        scenarios = {}
        try:
            for name, scenario in host_scenarios.items():
                result = torch_multigpu_scan_runtime(candidate['shards'], candidate['shared_capacities'],
                    candidate.get('shared_links', ()), ordered=False,
                    queue_service=candidate.get('queue_service'), result_queue_depth=shape['result_queue_depth'],
                    **scenario)
                seconds = result.get('estimated_scan_seconds')
                if isinstance(seconds, bool) or seconds is None or not math.isfinite(seconds) or seconds <= 0:
                    raise ValueError('No finite positive execution-graph estimate')
                scenarios[name] = result
        except ValueError as error:
            rejected.append(dict(candidate_index=index, reason='unsupported_model_context', detail=str(error)))
            continue
        for result in scenarios.values():
            unknown.update(result['unpriced_terms'])
        kwargs = dict(devices=shape['devices'], chunk_size=shape['chunk_size'],
                      reader_workers=shape['reader_workers'], prefetch_chunks=shape['depth'],
                      result_queue_depth=shape['result_queue_depth'], ordered=False,
                      compute_dtype='float32', compute_p_values=False, borrow_results=True)
        # Native readers prefer source.decode_workers. Setting only the scan
        # argument would leave single-device candidates on the source default.
        source_kwargs = dict(genotype_path=shape['input_path'], genotype_format='pgen', pgen_mode='hardcall',
                             reader_workers=shape['reader_workers'], pgen_decode_workers=shape['reader_workers'],
                             prefetch_chunks=shape['depth'])
        ranked.append(dict(candidate_index=index, devices=shape['devices'], chunk_size=shape['chunk_size'],
            memory=memory, scenarios=scenarios,
            worst_supplied_scenario_seconds=max(r['estimated_scan_seconds'] for r in scenarios.values()),
            entrypoint='torchgwas.linear.linear_scan_multigpu', scan_kwargs=kwargs,
            source_entrypoint='torchgwas.io.load_genotype', source_kwargs=source_kwargs,
            required_environment=dict(TORCHGWAS_NATIVE_STATS='0', TORCHGWAS_PGEN_PACKED='0',
                TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_SCAN_PROFILE='0',
                TORCHGWAS_BLOCKING_EVENTS='1' if shape['blocking_events'] else '0'),
            required_torch_settings=dict(allow_tf32=False),
            execution_requirements=['Reuse the CPU affinity, thread and library/device context of the independent prices.',
                'Verify input/source fingerprints and complete phenotypes with the stated covariate rank before execution.',
                'Consume borrowed views before requesting another chunk; copying, retaining or writing output changes the objective.']))
    if not ranked:
        raise ValueError('No feasible supported detailed candidate: ' + str(rejected))
    ranked.sort(key=lambda row: (row['worst_supplied_scenario_seconds'], row['candidate_index']))
    best = ranked[0]['worst_supplied_scenario_seconds']
    eligible = [row for row in ranked if row['worst_supplied_scenario_seconds'] <= best*(1+max_slowdown_fraction)]
    selected = min(eligible, key=lambda row: (len(row['devices']), row['worst_supplied_scenario_seconds'], row['candidate_index']))
    return dict(status='conditional_detailed_scan_plan', workload=identity, selected=selected,
        shortlist=[selected]+[row for row in ranked if row is not selected][:shortlist_size-1],
        candidates_evaluated=len(candidates), candidates_feasible=len(ranked), rejected=rejected,
        objective='minimize worst supplied execution-graph scenario', max_slowdown_fraction=max_slowdown_fraction,
        predicted_selection_penalty_fraction=selected['worst_supplied_scenario_seconds']/best-1,
        runtime_prediction_validated=False, prediction_complete=False,
        unresolved_timing_terms=sorted(unknown),
        scope='Complete-phenotype, matching-order native-int8 PGEN variant scans to immediate borrowed-result consumption. '
              'Setup, phenotype tiling, p-values, retained results and durable output are excluded. '
              'Memory checks are source storage estimates plus explicit reserves, not allocator guarantees.')
