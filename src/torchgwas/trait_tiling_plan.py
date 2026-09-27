"""Explicitly bounded development selection for durable trait-tiled output."""
import math

from .mechanistic_plan import _device, _integer
from .trait_tiling_model import (trait_tiled_shape, trait_tiled_memory, torch_trait_tiled_runtime,
                               reuse_trait_work, trait_tiled_graph_chunks)


@reuse_trait_work
def detailed_trait_plan(candidates, *, host_scenarios, cpu_workers, host_memory_bytes,
                        device_memory_bytes, host_reserve_bytes=0, device_reserve_bytes=0,
                        device_memory_profiles=None, max_candidates=1000,
                        max_scenario_evaluations=4000, max_tiles=10000, max_chunk_evaluations=1000000,
                        shortlist_size=5, max_slowdown_fraction=0.):
    """Minimize the worst supplied tile graph, with hard search/resource budgets.

    Limited measured selection checks exist; general selection remains unvalidated.
    It never collects scan timings, fits coefficients or uses the coarse model.
    """
    for name,value in [('cpu_workers',cpu_workers),('host_memory_bytes',host_memory_bytes),
                       ('max_candidates',max_candidates),('max_scenario_evaluations',max_scenario_evaluations),
                       ('max_tiles',max_tiles),('max_chunk_evaluations',max_chunk_evaluations),('shortlist_size',shortlist_size)]:
        _integer(name,value)
    for name,value in [('host_reserve_bytes',host_reserve_bytes),('device_reserve_bytes',device_reserve_bytes)]:
        _integer(name,value,0)
    if isinstance(max_slowdown_fraction,bool) or not math.isfinite(max_slowdown_fraction) or max_slowdown_fraction<0:
        raise ValueError('Invalid max_slowdown_fraction')
    if not isinstance(candidates,(list,tuple)) or not candidates or len(candidates)>max_candidates:
        raise ValueError('Nonempty candidate list within max_candidates required')
    if not isinstance(host_scenarios,dict) or not host_scenarios:
        raise ValueError('Explicit named host scenarios required')
    for name,scenario in host_scenarios.items():
        if not isinstance(name,str) or not name or set(scenario)-{'pageable_arena_fresh_fraction'}!={'host_serial_fraction','host_serial_policy'}:
            raise ValueError('Each named scenario needs host_serial_fraction and host_serial_policy')
        value = scenario['host_serial_fraction']
        if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError('Invalid host_serial_fraction')
        if scenario['host_serial_policy'] not in ('fluid','held-first','held-last'):
            raise ValueError('Invalid host_serial_policy')
        if 'pageable_arena_fresh_fraction' in scenario:
            value=scenario['pageable_arena_fresh_fraction']
            if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=1:
                raise ValueError('Invalid pageable_arena_fresh_fraction')
    for device,value in device_memory_bytes.items():
        _device(device);_integer('device_memory_bytes',value)
    shapes = [trait_tiled_shape(candidate) for candidate in candidates]
    columns=[c['tiles'][0]['data'].get('covariate_columns',shape['covariates']) for c,shape in zip(candidates,shapes)]
    if any(value!=columns[0] for value in columns):
        raise ValueError('Candidates have different covariate column counts')
    keys = ['samples','markers','traits','covariates','input_path']
    identity = {key:shapes[0][key] for key in keys}
    if columns[0]!=identity['covariates']:identity['covariate_columns']=columns[0]
    output_identity = {key:candidates[0]['output'][key] for key in ['store_beta','fsync']}
    if any(any(shape[key]!=identity[key] for key in keys) for shape in shapes):
        raise ValueError('Candidates must represent the same complete workload')
    if any(any(candidate['output'][key]!=value for key,value in output_identity.items()) for candidate in candidates):
        raise ValueError('Candidates must preserve the same durable output fields')
    rows,rejected,unpriced,feasible = [],[],set(),[]
    for index,(candidate,shape) in enumerate(zip(candidates,shapes)):
        if not set(shape['devices'])<=set(device_memory_bytes):
            raise ValueError('Missing device memory budget')
        if shape['reader_workers']>cpu_workers:
            rejected.append(dict(candidate_index=index,reason='shared_reader_budget'));continue
        memory = trait_tiled_memory(candidate,host_reserve_bytes=host_reserve_bytes,
            device_reserve_bytes=device_reserve_bytes,device_memory_profiles=device_memory_profiles)
        if memory['host_bytes']>host_memory_bytes:
            rejected.append(dict(candidate_index=index,reason='host_memory',memory=memory));continue
        if any(memory['device_bytes'][d]>device_memory_bytes[d] for d in shape['devices']):
            rejected.append(dict(candidate_index=index,reason='device_memory',memory=memory));continue
        feasible.append((index,candidate,shape,memory))
    # These limits bound graph work, which rejected proposals never perform.
    # Check the entire feasible set before executing any scenario: truncating
    # it would silently change the optimization problem or make order matter.
    scenarios_count=len(host_scenarios)
    if len(feasible)*scenarios_count>max_scenario_evaluations:
        raise ValueError('Search exceeds max_scenario_evaluations for memory/reader-feasible candidates')
    if sum(len(c['tiles']) for _,c,_,_ in feasible)*scenarios_count>max_tiles:
        raise ValueError('Search exceeds max_tiles across feasible scenario evaluations')
    if sum(trait_tiled_graph_chunks(c) for _,c,_,_ in feasible)*scenarios_count>max_chunk_evaluations:
        raise ValueError('Search exceeds max_chunk_evaluations for memory/reader-feasible candidates')
    for index,candidate,shape,memory in feasible:
        scenarios = {}
        try:
            for name,scenario in host_scenarios.items():
                result = torch_trait_tiled_runtime(candidate,**scenario)
                seconds = result['estimated_tile_seconds']
                if isinstance(seconds,bool) or not math.isfinite(seconds) or seconds<=0:
                    raise ValueError('No positive finite tile graph estimate')
                scenarios[name] = result
        except (KeyError,ValueError) as error:
            rejected.append(dict(candidate_index=index,reason='unsupported_model_context',detail=str(error)));continue
        for result in scenarios.values():unpriced.update(result['unpriced_terms'])
        output = candidate['output']
        kwargs = dict(genotype=shape['input_path'],genotype_format='pgen',pgen_mode='hardcall',
            device=shape['devices'][0],trait_devices=shape['devices'],trait_block=shape['trait_block'],
            chunk_size=shape['chunk_size'],reader_workers=shape['reader_workers'],prefetch_chunks=shape['depth'],
            compute_dtype='float32',sumstats_format='binary',sumstats_fsync=True,
            sumstats_fields='beta+t' if output['store_beta'] else 't',sumstats_block_bytes=output['block_bytes'],
            sumstats_queue_depth=output['queue_depth'])
        variant_partition=shape.get('partition_axis')=='variant'
        if variant_partition:
            kwargs.pop('trait_devices');kwargs.pop('trait_block')
            kwargs['variant_devices']=shape['devices']
        rows.append(dict(candidate_index=index,devices=shape['devices'],trait_block=None if variant_partition else shape['trait_block'],
            chunk_size=shape['chunk_size'],genotype_passes=1 if variant_partition else len(candidate['tiles']),memory=memory,scenarios=scenarios,
            worst_supplied_scenario_seconds=max(r['estimated_tile_seconds'] for r in scenarios.values()),
            entrypoint='torchgwas.api.run_linear_gwas',api_kwargs=kwargs,
            required_environment=dict(TORCHGWAS_NATIVE_STATS='0',TORCHGWAS_PGEN_PACKED='0',
                TORCHGWAS_PGEN_BACKEND='native',TORCHGWAS_SCAN_PROFILE='0',
                TORCHGWAS_BLOCKING_EVENTS='1' if shape['blocking_events'] else '0'),
            required_torch_settings=dict(allow_tf32=False),
            execution_requirements=['Supply complete phenotype and covariate inputs plus output_dir.',
                'Verify source/input fingerprints, full sample order, covariate rank and independent price context.',
                'Reuse the CPU affinity, thread counts and library/device context of the supplied primitives.']))
        if variant_partition:rows[-1]['partition_axis']='variant'
    rejected.sort(key=lambda row:row['candidate_index'])
    if not rows:raise ValueError('No feasible supported detailed trait candidate: '+str(rejected))
    rows.sort(key=lambda row:(row['worst_supplied_scenario_seconds'],row['candidate_index']))
    best = rows[0]['worst_supplied_scenario_seconds']
    eligible = [row for row in rows if row['worst_supplied_scenario_seconds']<=best*(1+max_slowdown_fraction)]
    selected = min(eligible,key=lambda row:(len(row['devices']),row['worst_supplied_scenario_seconds'],row['candidate_index']))
    result=dict(status='development_detailed_trait_plan',workload=identity,selected=selected,
        shortlist=[selected]+[row for row in rows if row is not selected][:shortlist_size-1],
        admission_candidates=[{key:value for key,value in row.items() if key!='scenarios'}
            for row in rows],
        candidates_evaluated=len(candidates),candidates_feasible=len(rows),rejected=rejected,
        objective='minimize worst supplied tile execution-graph scenario',
        predicted_selection_penalty_fraction=selected['worst_supplied_scenario_seconds']/best-1,
        max_slowdown_fraction=max_slowdown_fraction,prediction_complete=False,
        runtime_prediction_validated=False,selection_validated=False,unresolved_timing_terms=sorted(unpriced),
        scope='Explicit full-output trait tiling, chunk and device selection over prepared independent profiles. '
              'Memory includes per-device pinned size-class retention and active writers. '
              'Tile-worker timing excludes input QC and final metadata publication; unpriced terms can change rankings.')
    if any(shape.get('partition_axis')=='variant' for shape in shapes):
        result.update(status='development_detailed_output_plan',
            objective='minimize worst supplied durable-output execution-graph scenario',
            scope='Joint trait/variant partition, chunk and device selection with one shared setup/scan/writer graph. Replicated phenotype setup, disjoint variant ranges or repeated full-input trait passes, aggregate resources, per-device pinned retention and durable writers are accounted. Input QC and final metadata remain outside the priced boundary.')
    return result
