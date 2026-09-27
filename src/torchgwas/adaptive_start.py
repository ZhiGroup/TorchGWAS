"""Admit one productive starting layout without evaluating runtime candidates.

The caller supplies one dependency-validated hardware context and live resource
budgets. Only PGEN locators and mixed-size memory requests are built
here. Runtime ranking and association measurement cannot occur on this path.
"""
from copy import deepcopy
import time

from .adaptive_chunks import aligned_chunk_sizes,aligned_chunk_shapes
from .adaptive_candidate import adaptive_candidate_memory
from .analytical_plan_cache import input_identity
from .mechanistic_plan import _integer
from .trait_candidate_space import prepare_trait_candidates,_tile_profile
from .trait_tiling_model import trait_tiled_shape


def prepare_adaptive_start(workload, context, *, chunk_sizes, initial_size,
                           partition_axis, trait_block, reduction, output,
                           cpu_workers, host_memory_bytes, device_memory_bytes,
                           host_reserve_bytes, device_reserve_bytes,
                           device_memory_profiles, significance_threshold=None,
                           max_tiles=10000, max_census_chunks=1000000,
                           max_kernel_shapes=128, _header_receiver=None,
                           _prepared_header=None, _index_receiver=None,
                           _compact_memory=False, _compact_cache_dir=None,
                           _compact_cache_receiver=None, _defer_compact_bases=False):
    """Prepare exactly one fixed layout and admit all reachable chunk shapes.

    This does not choose the initial layout for the caller or read timing
    evidence. In particular it never calls the bounded performance optimizers.
    Chunk capacity is the largest choice, independent of the first work size.
    Host/device budgets must be capped by fresh available resources; the bridge
    must recheck them and the full input/context binding before starting scans.
    """
    started=time.perf_counter()
    sizes=aligned_chunk_sizes(chunk_sizes)
    if type(initial_size) is not int or initial_size not in sizes:
        raise ValueError('Initial size must be an admitted chunk choice')
    if reduction not in (None,'significant','jagwas'):
        raise ValueError('Adaptive start supports full, significant and jagwas output')
    if partition_axis not in ('trait','variant'):
        raise ValueError('Explicit trait or variant partition axis required')
    if reduction=='jagwas' and (partition_axis!='variant' or trait_block is not None):
        raise ValueError('JAGWAS requires variant partitioning and an unpartitioned phenotype panel')
    if reduction=='significant' and partition_axis!='trait':
        raise ValueError('Significant-pairs adaptive start requires trait partitioning')
    if partition_axis=='variant' and trait_block is not None:
        raise ValueError('Variant partitioning cannot also request a phenotype tile')
    width=workload['traits'] if partition_axis=='variant' else trait_block
    _integer('trait_block',width)
    if reduction=='significant':
        import math
        if significance_threshold is not None and (isinstance(significance_threshold,bool)
            or not isinstance(significance_threshold,(int,float)) or not math.isfinite(significance_threshold)
            or not 0<significance_threshold<=1):
            raise ValueError('Significance threshold must be in (0,1] or None')
    elif significance_threshold is not None:
        raise ValueError('A significance threshold applies only to significant pairs')
    if reduction is not None and output.get('block_bytes') is not None:
        raise ValueError('Reduced output requires indexed parts without dense coalescing')
    for name,value in [('cpu_workers',cpu_workers),('host_memory_bytes',host_memory_bytes),
        ('host_reserve_bytes',host_reserve_bytes),('device_reserve_bytes',device_reserve_bytes),
        ('max_tiles',max_tiles),('max_census_chunks',max_census_chunks),('max_kernel_shapes',max_kernel_shapes)]:
        _integer(name,value)
    devices=context['devices']
    if set(device_memory_bytes)!=set(devices) or set(device_memory_profiles)!=set(devices):
        raise ValueError('Exactly one budget and memory profile per starting device required')
    for value in device_memory_bytes.values():_integer('device_memory_bytes',value)
    readers=[context['profiles'][device]['decode_workers'] for device in devices]
    for value in readers:_integer('decode_workers',value)
    if sum(readers)>cpu_workers:
        raise ValueError('Starting readers exceed the shared CPU worker budget')
    identity=input_identity(workload['genotype'])
    # Only locators are needed to bound input/scratch allocations. A full
    # decoder census would read compressed payloads before useful output.
    space=prepare_trait_candidates(workload,[context],chunks=[sizes[-1]],trait_blocks=[width],
        partition_axes=[partition_axis],reduction='jagwas' if reduction=='jagwas' else None,
        output=output,max_candidates=1,max_candidate_tiles=max_tiles,
        max_census_chunks=max_census_chunks,census_chunk_size=sizes[0],memory_only=True,
        _header_receiver=_header_receiver,_prepared_header=_prepared_header,
        _index_receiver=_index_receiver,_compact_memory=_compact_memory,
        _compact_cache_dir=_compact_cache_dir,_compact_cache_receiver=_compact_cache_receiver,
        _defer_compact_bases=_defer_compact_bases)
    if len(space['candidates'])!=1:
        raise ValueError('Exactly one active starting layout required')
    candidate=space['candidates'][0];source=space['source_census']
    for tile in candidate['tiles']:
        data=tile['data'];n,m,k,c=[data[key] for key in ('samples','markers','traits_analyzed','covariates')]
        # Static preparation retains only the capacity/tail kernels. Restore
        # the finite set needed for every reachable adaptive work shape.
        profile=_tile_profile(context['profiles'][tile['device']],n,m,k,c,sizes[-1],
                              work_shapes=aligned_chunk_shapes(sizes,m))
        profile['return_beta']=output['store_beta'] if reduction!='jagwas' else True
        profile['chunk_markers']=sizes[-1];tile['profile']=profile
    shape=trait_tiled_shape(candidate,reduction='jagwas' if reduction=='jagwas' else None)
    if shape['reader_workers']>cpu_workers:
        raise ValueError('Starting readers exceed the shared CPU worker budget')
    memory=adaptive_candidate_memory(candidate,chunk_sizes=sizes,source_census=source,
        reduction=reduction,device_memory_profiles=device_memory_profiles,
        host_reserve_bytes=host_reserve_bytes,device_reserve_bytes=device_reserve_bytes,
        max_census_chunks=max_census_chunks,max_kernel_shapes=max_kernel_shapes)
    if memory['host_bytes']>host_memory_bytes:
        raise ValueError('Adaptive starting layout exceeds host memory budget')
    if any(value>device_memory_bytes[device] for device,value in memory['device_bytes'].items()):
        raise ValueError('Adaptive starting layout exceeds device memory budget')
    if input_identity(workload['genotype'])!=identity:
        raise ValueError('Input changed during adaptive startup admission')
    settings=dict(genotype=shape['input_path'],genotype_format='pgen',pgen_mode='hardcall',
        device=shape['devices'][0],chunk_size=sizes[-1],reader_workers=shape['reader_workers'],
        prefetch_chunks=shape['depth'],compute_dtype='float32',sumstats_format='binary',
        sumstats_fsync=True,sumstats_fields='beta+t' if output['store_beta'] else 't',
        sumstats_block_bytes=output['block_bytes'],sumstats_queue_depth=output['queue_depth'])
    if partition_axis=='variant':settings['variant_devices']=list(shape['devices'])
    else:settings.update(trait_devices=list(shape['devices']),trait_block=width)
    if reduction is not None:settings['reduce']=reduction
    if reduction=='significant':settings['significance_threshold']=significance_threshold
    partitions=[dict(id=str(i),device=tile['device'],trait_range=list(tile['trait_range']),
        variant_range=list(tile['data']['encoded']['variant_range'])) for i,tile in enumerate(candidate['tiles'])]
    environment=dict(TORCHGWAS_NATIVE_STATS='0',TORCHGWAS_PGEN_PACKED='0',
        TORCHGWAS_PGEN_BACKEND='native',TORCHGWAS_SCAN_PROFILE='0',
        TORCHGWAS_BLOCKING_EVENTS='1' if shape['blocking_events'] else '0')
    if reduction=='significant':environment['TORCHGWAS_SIGNIFICANCE_BACKEND']='host'
    return dict(candidate=candidate,source_layout=source,source_census=None,memory=memory,partitions=partitions,
        chunk_sizes=list(sizes),initial_size=initial_size,capacity=sizes[-1],api_kwargs=settings,
        required_environment=environment,required_torch_settings=dict(allow_tf32=False),
        input_file_identity=identity,context=context['name'],reduction=reduction,
        timing_geometry_missing=deepcopy(memory['missing_geometry']),runtime_work_ready=False,
        admission_seconds=time.perf_counter()-started,
        structural_work=dict(census_passes=space['census_passes'],header_passes=1,census_views=space['census_views'],
            census_chunks=space['census_chunks'],layouts_built=1,runtime_candidates_evaluated=0),
        resource_budgets=dict(cpu_workers=cpu_workers,host_memory_bytes=host_memory_bytes,
            device_memory_bytes=deepcopy(device_memory_bytes)),
        scope='One fixed starting layout, all allowed mixed chunk lifetimes, and exact source restart/workspace admission. No runtime ranking, graph solve, calibration run, or claim of optimal layout. The public execution bridge must still bind current inputs/source/prices and live resources, and attach a productive chunk controller to use initial_size.')
