"""Finite significant-host candidate selection over the shared analytical graph."""
import copy
import math
from .mechanistic_plan import _integer, _device
from .trait_tiling_model import trait_tiled_shape, reuse_trait_work
from .significant_host_work import significant_host_memory
from .significant_host_model import significant_host_runtime


@reuse_trait_work
def detailed_significant_host_plan(candidates, prices, *, significance_threshold,
        occupancy_scenarios, host_scenarios, cpu_workers, host_memory_bytes,
        device_memory_bytes, device_memory_profiles=None, host_reserve_bytes=0,
        device_reserve_bytes=0, max_candidates=64, max_scenario_evaluations=256,
        max_source_chunks=10000, max_selection_blocks=100000):
    """Minimize worst declared scenario; admit memory at dense occupancy.

    This development bridge does not authorize guarded API execution. Profile
    provenance, empirical selection coverage and device filtering are separate.
    """
    for name,value in [('cpu_workers',cpu_workers),('host_memory_bytes',host_memory_bytes),
        ('max_candidates',max_candidates),('max_scenario_evaluations',max_scenario_evaluations),
        ('max_source_chunks',max_source_chunks),('max_selection_blocks',max_selection_blocks)]:
        _integer(name,value)
    for name,value in [('host_reserve_bytes',host_reserve_bytes),('device_reserve_bytes',device_reserve_bytes)]:
        _integer(name,value,0)
    if not isinstance(device_memory_bytes,dict) or not device_memory_bytes:
        raise ValueError('Explicit nonempty per-device memory capacities required')
    for device,capacity in device_memory_bytes.items():
        _device(device);_integer('device capacity',capacity)
    if not isinstance(candidates,(list,tuple)) or not candidates or len(candidates)>max_candidates:
        raise ValueError('Nonempty bounded significant candidate list required')
    if significance_threshold is not None and (isinstance(significance_threshold,bool) or not isinstance(significance_threshold,(int,float)) or not math.isfinite(significance_threshold) or not 0<significance_threshold<=1):
        raise ValueError('Significant threshold must be in (0,1] or None for global Bonferroni')
    if not isinstance(occupancy_scenarios,dict) or not occupancy_scenarios or any(
            not isinstance(name,str) or not name or value not in ('empty','dense')
            for name,value in occupancy_scenarios.items()):
        raise ValueError('Named empty/dense occupancy scenarios required; counts may not be inferred from alpha')
    if not isinstance(host_scenarios,dict) or not host_scenarios:
        raise ValueError('Explicit host resource scenarios required')
    for name,host in host_scenarios.items():
        if (not isinstance(name,str) or not name or not isinstance(host,dict)
                or set(host)-{'host_serial_policy'} != {'host_serial_fraction'}):
            raise ValueError('Named host scenario requires serialization fraction and optional policy')
        value=host['host_serial_fraction']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError('Finite bounded host serialization fraction required')
        if host.get('host_serial_policy','fluid') not in ('fluid','held-first','held-last'):
            raise ValueError('Unknown host serialization policy')
    shapes = [trait_tiled_shape(candidate) for candidate in candidates]
    keys = ['samples','markers','traits','covariates','input_path']
    identity = {key:shapes[0][key] for key in keys}
    if any(any(shape[key]!=value for key,value in identity.items()) for shape in shapes):
        raise ValueError('Candidates must represent the same statistical workload')
    store_beta = candidates[0]['output']['store_beta']
    if any(candidate['output']['store_beta'] != store_beta for candidate in candidates):
        raise ValueError('Candidate output fields must be identical')
    rows, rejected, feasible = [], [], []
    for index,(candidate,shape) in enumerate(zip(candidates,shapes)):
        if shape['reader_workers']>cpu_workers:
            rejected.append(dict(candidate_index=index,reason='shared_reader_budget'));continue
        if not set(shape['devices'])<=set(device_memory_bytes):
            raise ValueError('Missing device memory capacity')
        memory = significant_host_memory(candidate,host_reserve_bytes=host_reserve_bytes,
            device_reserve_bytes=device_reserve_bytes,device_memory_profiles=device_memory_profiles)
        if memory['host_bytes']>host_memory_bytes:
            rejected.append(dict(candidate_index=index,reason='host_memory',memory=memory));continue
        if any(memory['device_bytes'][d]>device_memory_bytes[d] for d in shape['devices']):
            rejected.append(dict(candidate_index=index,reason='device_memory',memory=memory));continue
        feasible.append((index,candidate,shape,memory))
    if len(feasible)*len(occupancy_scenarios)*len(host_scenarios)>max_scenario_evaluations:
        raise ValueError('Significant search exceeds max_scenario_evaluations')
    unresolved = set()
    for index,candidate,shape,memory in feasible:
        scenarios = {}
        for occupancy_name,occupancy in occupancy_scenarios.items():
            for host_name,host in host_scenarios.items():
                result = significant_host_runtime(candidate,prices,occupancy=occupancy,significance_threshold=significance_threshold,**host,
                    max_source_chunks=max_source_chunks,max_selection_blocks=max_selection_blocks)
                seconds=result.get('estimated_tile_seconds')
                if isinstance(seconds,bool) or not isinstance(seconds,(int,float)) or not math.isfinite(seconds) or seconds<=0:
                    raise ValueError('No finite positive significant candidate score')
                scenarios[(occupancy_name,host_name)] = result
                unresolved.update(result['unpriced_terms'])
        score = max(value['estimated_tile_seconds'] for value in scenarios.values())
        if not math.isfinite(score) or score<=0:
            raise ValueError('No finite positive significant candidate score')
        rows.append(dict(candidate_index=index,devices=shape['devices'],chunk_size=shape['chunk_size'],
            trait_block=shape['trait_block'],memory=memory,
            scenarios=[dict(occupancy=label[0],host=label[1],estimate=value) for label,value in scenarios.items()],
            worst_supplied_scenario_seconds=score,api_kwargs=dict(genotype=shape['input_path'],
                genotype_format='pgen',pgen_mode='hardcall',device=shape['devices'][0],
                trait_devices=shape['devices'],trait_block=shape['trait_block'],chunk_size=shape['chunk_size'],
                reader_workers=shape['reader_workers'],prefetch_chunks=shape['depth'],compute_dtype='float32',
                reduce='significant',significance_threshold=significance_threshold,sumstats_format='binary',
                sumstats_fsync=True,sumstats_queue_depth=candidate['output']['queue_depth'],sumstats_fields='beta+t' if store_beta else 't'),
            required_environment=dict(TORCHGWAS_SIGNIFICANCE_BACKEND='host',TORCHGWAS_PGEN_BACKEND='native',
                TORCHGWAS_PGEN_PACKED='0',TORCHGWAS_NATIVE_STATS='0',TORCHGWAS_SCAN_PROFILE='0',
                TORCHGWAS_BLOCKING_EVENTS='1' if shape['blocking_events'] else '0')))
    if not rows:
        raise ValueError('No feasible significant host candidate: '+str(rejected))
    rows.sort(key=lambda row:(row['worst_supplied_scenario_seconds'],len(row['devices']),row['candidate_index']))
    return dict(status='development_detailed_significant_host_plan',selected=rows[0],candidates=rows,rejected=rejected,
        workload=identity,significance_threshold=significance_threshold,
        occupancy_scenarios=copy.deepcopy(occupancy_scenarios),host_scenarios=copy.deepcopy(host_scenarios),
        candidates_evaluated=len(candidates),candidates_feasible=len(rows),selection_validated=False,
        objective='minimize worst explicitly supplied survivor/resource scenario',unpriced_terms=sorted(unresolved),
        scope='One analytical setup/scan/host-selection/shared-writer graph with dense-occupancy array admission. Additional formats, device filtering and guarded API binding remain required; no GWAS timing is fitted.')