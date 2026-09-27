"""Compose full-output trait tiles using the existing analytical scan graph.

No scan durations or fitted tile rates enter this model. This is a development
objective, from tile-worker entry through payload fsync. API input QC, file
page faults, manifests and final ID publication remain explicitly unpriced.
"""
from __future__ import annotations

import copy
import math
import json
from collections import OrderedDict
from contextvars import ContextVar
from functools import wraps

from .tensor_work import tensor_work_cache
from .planning_session import planning_work_scope
from .decoder_work import native_reader_workspace,census_chunk_ranges,scan_chunk_count,require_regular_memory

from .binary_output_work import binary_output_work, binary_output_memory
from .binary_schedule import BinaryWriterSchedule
from .execution_graph import ExecutionGraph, torch_scan_schedule
from .mechanistic_plan import _device, _integer
from .mechanistic_torch import torch_scan_work, handoff_summary
from .pinned_work import pinned_scan_work
from .output_write_work import writer_copy_cost
from .setup_work import setup_work, setup_service, setup_memory
from .tensor_memory import eager_scan_memory, eager_memory_plan


_SCAN_WORK_CACHE = ContextVar('torchgwas_trait_scan_work_cache', default=None)


def reuse_trait_work(function):
    @wraps(function)
    def run(*args,**kwargs):
        token=_SCAN_WORK_CACHE.set(OrderedDict())
        try:
            with tensor_work_cache(),planning_work_scope():return function(*args,**kwargs)
        finally:_SCAN_WORK_CACHE.reset(token)
    return run


def _cached_scan_work(data,profile):
    """Read-only reuse inside one plan; mutable graph blocks are copied below."""
    cache=_SCAN_WORK_CACHE.get()
    if cache is None:return torch_scan_work(data,profile)
    import torch
    try:
        key=json.dumps([data,profile,str(torch.get_default_dtype())],sort_keys=True,
                       separators=(',',':'),allow_nan=False)
    except (TypeError,ValueError):
        # Caching must not narrow the calculator's accepted Python input types.
        return torch_scan_work(data,profile)
    if key not in cache:
        result=torch_scan_work(data,profile)
        if len(cache)>=128:cache.popitem(last=False)
        cache[key]=result
    cache.move_to_end(key)
    return cache[key]


def trait_tiled_shape(candidate, *, reduction=None):
    if reduction not in (None,'jagwas','device_significant'):raise ValueError('Unsupported candidate reduction')
    allowed = {'tiles', 'trait_block', 'devices', 'shared_capacities', 'shared_links', 'shared_storage_bytes_per_second', 'output', 'partition_axis'}
    axis=candidate.get('partition_axis','trait')
    if axis not in ('trait','variant'):raise ValueError('Unknown output partition axis')
    if reduction=='jagwas' and axis!='variant':raise ValueError('JAGWAS requires variant partitioning')
    if reduction=='device_significant' and axis!='trait':raise ValueError('Device significance requires trait partitioning')
    if set(candidate) - allowed:
        raise ValueError('Unknown tiled candidate fields: ' + str(sorted(set(candidate)-allowed)))
    tiles = candidate['tiles']
    if not isinstance(tiles, list) or not tiles:
        raise ValueError('Nonempty tile list required')
    width = _integer('trait_block', candidate['trait_block'])
    devices = [_device(d) for d in candidate['devices']]
    if not devices or len(set(devices)) != len(devices) or len(devices) > len(tiles):
        raise ValueError('Unique active tile devices required')
    if axis=='variant' and len(devices)!=len(tiles):
        raise ValueError('Each variant shard requires one active device')
    first = tiles[0]
    n, m, c = [first['data'][key] for key in ['samples', 'markers', 'covariates']]
    if axis=='variant':m=sum(_integer('shard markers',tile['data']['markers']) for tile in tiles)
    for name, value, minimum in [('samples', n, 32), ('markers', m, 1), ('covariates', c, 0)]:
        _integer(name, value, minimum)
    columns=_integer('covariate_columns',first['data'].get('covariate_columns',c),c)
    if columns>=n-2:raise ValueError('Covariate dimensions require positive residual df')
    p = first['profile']
    chunk, depth = _integer('chunk_markers', p['chunk_markers']), _integer('depth', p['depth'], 2)
    if axis=='variant':
        from .linear import multigpu_variant_ranges
        spans=multigpu_variant_ranges(m,chunk,len(devices))
        if len(spans)!=len(tiles):raise ValueError('Idle devices in variant partition')
    event_wait = p.get('event_wait_cpu_fraction')
    if isinstance(event_wait, bool) or event_wait not in (0., 1.):
        raise ValueError('Explicit spin or blocking completion event policy required')
    output = candidate['output']
    if set(output) != {'block_bytes', 'queue_depth', 'store_beta', 'fsync'}:
        raise ValueError('Explicit block_bytes, queue_depth, store_beta and fsync required')
    if output['block_bytes'] is not None:
        _integer('block_bytes', output['block_bytes'])
    _integer('queue_depth', output['queue_depth'])
    if type(output['store_beta']) is not bool or output['fsync'] is not True:
        raise ValueError('Boolean store_beta and durable fsync output required')
    caps = candidate['shared_capacities']
    if set(caps) != {'cpu', 'dram', 'input', 'output'}:
        raise ValueError('Declare aggregate cpu, dram, input and output capacities')
    if any(isinstance(v, bool) or not math.isfinite(v) or v <= 0 for v in caps.values()):
        raise ValueError('Positive finite shared capacities required')
    workers = {}
    start = 0
    path = first['data']['encoded'].get('path')
    if not path:
        raise ValueError('Input census path required')
    for index, tile in enumerate(tiles):
        if set(tile) != ({'trait_range', 'device', 'data', 'profile'} | ({'variant_range'} if axis=='variant' else set())):
            raise ValueError('Tile must contain only trait_range, device, data and profile')
        data, profile = tile['data'], tile['profile']
        if 'phenotype_c_contiguous' in data and type(data['phenotype_c_contiguous']) is not bool:
            raise ValueError('Explicit boolean phenotype_c_contiguous required')
        if data.get('covariate_columns',c)!=columns:raise ValueError('Covariate column counts differ between tiles')
        k = _integer('tile traits', data['traits_analyzed'])
        if axis=='variant':
            if k!=width or tile['trait_range']!=[0,width] or tile['variant_range']!=list(spans[index]):
                raise ValueError('Variant shards must match chunk-balanced ranges and retain every trait')
        else:
            if k > width or (index < len(tiles)-1 and k != width) or tile['trait_range'] != [start, start+k]:
                raise ValueError('Tiles must cover consecutive full-width traits and one optional tail')
            start += k
        device = devices[index % len(devices)]
        if tile['device'] != device:
            raise ValueError('Tile devices must follow executor round-robin assignment')
        expected_markers=spans[index][1]-spans[index][0] if axis=='variant' else m
        if (data['samples'], data['markers'], data['covariates']) != (n, expected_markers, c) or data.get('matching_sample_order') is not True:
            raise ValueError('Matching complete samples and dimensions required')
        if (profile['chunk_markers'], profile['depth']) != (chunk, depth):
            raise ValueError('Common chunk and prefetch depth required')
        ownership='owned' if reduction in ('jagwas','device_significant') else 'borrowed'
        if profile.get('result_ownership') != ownership or profile.get('validate_range') is not True:
            raise ValueError('Tiled executor requires '+ownership+' results and range validation')
        if reduction=='jagwas':
            if profile.get('reduction')!='jagwas' or profile.get('compute_dtype','float32')!='float32':
                raise ValueError('JAGWAS candidate requires native FP32 statistics and joint reduction')
            if profile.get('borrow_results',False):raise ValueError('JAGWAS candidate cannot borrow results')
            if data.get('phenotype_complete') is not True:raise ValueError('Explicit complete phenotype contract required')
            if k>n-c-1:raise ValueError('JAGWAS trait count exceeds residual phenotype rank')
        if reduction=='device_significant':
            if profile.get('reduction')!=reduction or profile.get('compute_dtype','float32')!='float32':
                raise ValueError('Device significant candidate requires native FP32 selection')
            if profile.get('borrow_results',False):raise ValueError('Device significant results cannot be borrowed')
            if data.get('phenotype_complete') is not True:raise ValueError('Explicit complete phenotype contract required')
        if profile.get('event_wait_cpu_fraction') != event_wait:
            raise ValueError('Completion event policies differ')
        encoded = data['encoded']
        expected_span=list(spans[index]) if axis=='variant' else [0,m]
        if encoded.get('path') != path or encoded.get('variant_range') != expected_span or encoded.get('file_markers') != m:
            raise ValueError('Each tile must census the same full genotype input')
        if not encoded.get('chunks') or encoded.get('chunk_markers') != chunk:
            raise ValueError('Exact per-chunk input census required')
        for resource, field in [('cpu', 'cpu_available_cores'), ('dram', 'shared_dram_bytes_per_second'),
                                ('input', 'read_bytes_per_second'), ('output', 'write_bytes_per_second')]:
            # Shared availability can fall below the independently calibrated
            # service ceiling; leave the tile profile and its evidence intact.
            value = profile[field]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0 or caps[resource] > value:
                raise ValueError('Shared capacity exceeds tile profile: ' + resource)
        count = _integer('decode_workers', profile['decode_workers'])
        if device in workers and workers[device] != count:
            raise ValueError('Reader allocation changes between tiles on one device')
        workers[device] = count
    total = sum(workers.values())
    each, remainder = divmod(total, len(devices))
    if [workers[d] for d in devices] != [each+(i<remainder) for i in range(len(devices))]:
        raise ValueError('Reader allocation must match executor divmod budget')
    result=dict(samples=n, markers=m, covariates=c, traits=width if axis=='variant' else start, trait_block=width,
                devices=devices, reader_workers=total, chunk_size=chunk, depth=depth,
                input_path=path, blocking_events=event_wait == 0.)
    if axis=='variant':result['partition_axis']='variant'
    return result


def _output_work(tile, output):
    encoded=tile['data']['encoded']
    rows=([hi-lo for lo,hi in census_chunk_ranges(encoded,tile['profile']['chunk_markers'])]
          if 'chunk_ranges' in encoded else None)
    return binary_output_work(tile['data']['markers'], tile['data']['traits_analyzed'],
        tile['profile']['chunk_markers'], **output, borrow_chunks=False, store_variant_df=True,chunk_rows=rows)


def trait_tiled_memory(candidate, *, host_reserve_bytes=0, device_reserve_bytes=0,
                       device_memory_profiles=None):
    """Per-device peaks and retained pinned size classes, plus explicit reserves."""
    require_regular_memory(candidate)
    shape = trait_tiled_shape(candidate)
    for key, value in [('host_reserve_bytes', host_reserve_bytes), ('device_reserve_bytes', device_reserve_bytes)]:
        _integer(key, value, 0)
    host_peaks, gpu_peaks, pin_bins = {}, {}, {}
    # Source arrays and writer blocks depend on extents, not trait offsets.
    # Scope reuse to this call: input censuses and device profiles stay mutable
    # to callers, so identity-based reuse must not survive another request.
    tile_memory, census_max_payload = {}, {}
    unknown = set()
    for tile in candidate['tiles']:
        device, profile = tile['device'], tile['profile']
        n, m, k, c = [tile['data'][key] for key in ['samples', 'markers', 'traits_analyzed', 'covariates']]
        columns=tile['data'].get('covariate_columns',c)
        b, depth = profile['chunk_markers'], profile['depth']
        encoded=tile['data']['encoded']
        if id(encoded) not in census_max_payload:
            census_max_payload[id(encoded)]=native_reader_workspace(encoded['chunks'])
        payload,replay_workspace=census_max_payload[id(encoded)]
        key=(device,n,m,k,c,columns,b,depth,profile['decode_workers'],payload,replay_workspace,profile.get('return_beta',True))
        if key not in tile_memory:
            tile_memory[key]=_tile_memory(n,m,k,c,columns,b,depth,profile['decode_workers'],payload,
                candidate['output'],(device_memory_profiles or {}).get(device),replay_workspace,return_beta=profile.get('return_beta',True))
        host,gpu,needed,unresolved=tile_memory[key]
        cache = pin_bins.setdefault(device, {})
        for size, count in needed.items():cache[size] = max(cache.get(size, 0), count)
        host_peaks[device] = max(host_peaks.get(device, 0), host)
        gpu_peaks[device] = max(gpu_peaks.get(device, 0), gpu)
        unknown.update(unresolved)
    pinned = {device:sum(size*count for size,count in bins.items()) for device,bins in pin_bins.items()}
    unknown.update(['input mappings/page residency, metadata/index objects and QC peaks require host reserve',
                    'LAPACK workspace and NumPy reduction temporaries require host reserve',
                    'GPU allocator caches, fragmentation and driver allocations require device reserve',
                    'pinned size-class retention assumes reusable completed blocks; allocator fragmentation is not bounded'])
    return dict(host_bytes=sum(host_peaks.values())+sum(pinned.values())+host_reserve_bytes,
        device_bytes={d:b+device_reserve_bytes for d,b in gpu_peaks.items()},
        host_active_bytes_by_device=host_peaks, pinned_cache_bytes_by_device=pinned,
        host_reserve_bytes=host_reserve_bytes, device_reserve_bytes=device_reserve_bytes,
        unresolved_memory_terms=sorted(unknown))


def _tile_memory(n,m,k,c,columns,b,depth,decode_workers,payload,output,memory_profile,replay_workspace=0,*,reduction=None,return_beta=True):
    """One source extent; caller combines simultaneous device peaks and pins."""
    pins = pinned_scan_work(n, b, k, depth,reduction=reduction,return_beta=return_beta)
    needed = {}
    for row in pins['allocations']:
        needed[row['allocator_bytes']] = needed.get(row['allocator_bytes'], 0)+row['count']
    active = min(depth, decode_workers, (m+b-1)//b)
    decoder = active*(min(b, m)*((n+3)//4)+payload+replay_workspace)
    properties = ({} if memory_profile is None else
        {key:memory_profile[key] for key in ['sm_count','max_threads_per_sm']})
    setup = setup_memory(n, k, c, **properties)
    ledger = setup_work(n,k,c)
    # Sum explicitly sized arrays even where lifetimes can be disjoint.
    # FP32 tile cast and residual result; residual-block contiguous upload
    # and download; contiguous design upload; observed mask/counts. Basis
    # construction retains promoted/centered/scaled FP64 arrays and U,
    # plus FP32 Q. LAPACK workspace remains a reserve, not an invented rate.
    host_arrays = dict(tile_cast=4*n*k,residual_result=4*n*k,
        residual_host_blocks=8*n*ledger['residual_block_traits'],
        design_host_block=4*n*ledger['design_block_traits'],
        observed_counts=8*k,
        covariate_basis_arrays=36*n*columns+8*columns*columns+24*columns)
    writer = ({'allocated_staging_bytes': 0} if output is None else
        binary_output_memory(m,k,b,block_bytes=output['block_bytes'],
        queue_depth=output['queue_depth'],store_beta=output['store_beta'],store_variant_df=True))
    host = sum(host_arrays.values())+decoder+writer['allocated_staging_bytes']
    eager = eager_scan_memory(n, b, k, c, depth,reduction=reduction,return_beta=return_beta)
    gpu = max(setup['device_live_bytes_upper'], eager['tensor_storage_budget'])
    unknown=set(setup['unresolved_memory_terms'])
    if memory_profile is not None:
        detail = eager_memory_plan(n, b, k, c, depth, memory_profile,reduction=reduction,return_beta=return_beta)
        gpu = max(gpu, detail['device_bytes'])
        unknown.update(detail['unresolved_memory_terms'])
    else:
        unknown.add('CUDA library workspace requires an independent device profile or reserve')
    return host,gpu,needed,unknown


def _prepare_graph(work, profile, serial_fraction, pin_pages, cached_calls):
    """Ordered phase aggregates; intra-phase API/kernel overlap is unresolved."""
    estimate = setup_service(work, profile)
    q = profile['cpu_fraction']
    g = ExecutionGraph()
    g.capacities['host_serial'] = 1.
    cpu = {'cpu':q, 'host_serial':q*serial_fraction}
    last = None
    def add(name, seconds, resources=None):
        nonlocal last
        last = g.add(name, seconds, [] if last is None else [last], resources)
    n = work['samples']
    # The API constructs one rank-validated basis before the partition executor
    # begins, and every tile/shard reuses it. No worker-local SVD is performed.
    # Keep the boundary node for ordering; the shared API setup is unpriced here.
    add('covariate_basis', 0., cpu)
    for index, (phase, cost) in enumerate(zip(work['phases'], estimate['phases'])):
        # Each source phase is blocking. Spread its resource demand over its
        # aggregate duration; this is an explicit contention scenario.
        seconds = cost['seconds']
        resources = {'cpu':q*(cost['fixed_cpu_seconds']+cost['host_copy_seconds']+cost['host_page_seconds']+cost['pageable_d2h_cpu_seconds'])/seconds if seconds else 0.}
        # NumPy 2.x numeric array assignment releases the GIL above 500
        # elements. Bulk copy still consumes a CPU core and shared DRAM;
        # only its small-copy case belongs in the Python serialization pool.
        held_copy = cost['host_copy_seconds'] if n*phase['traits']<=500 else 0.
        resources['host_serial'] = q*((cost['fixed_cpu_seconds']+held_copy)*serial_fraction+cost['host_page_serial_seconds'])/seconds if seconds else 0.
        resources['dram'] = cost['host_dram_bytes']/seconds if seconds else 0.
        for direction in ['h2d','d2h']:
            resources['prep_'+direction] = phase[direction+'_bytes']/seconds if seconds else 0.
        add('phase:'+str(index), seconds, resources)
    fresh_cpu = profile['pin_cpu_seconds_per_page']*pin_pages/q
    fresh_non_cpu = profile['pin_driver_seconds_per_page']*pin_pages
    cached_cpu = profile['pin_cached_cpu_seconds_per_call']*cached_calls/q
    add('pin_cpu', fresh_cpu+cached_cpu, cpu)
    add('pin_driver', fresh_non_cpu)
    return g, estimate


def _tile_pin_state(tile, caches):
    """Advance the per-device completed-block cache, including size-class tails."""
    data,profile=tile['data'],tile['profile']
    pins=pinned_scan_work(data['samples'],profile['chunk_markers'],data['traits_analyzed'],profile['depth'],reduction=profile.get('reduction'),return_beta=profile.get('return_beta',True))
    needed={}
    for row in pins['allocations']:
        needed[row['allocator_bytes']]=needed.get(row['allocator_bytes'],0)+row['count']
    cache=caches.setdefault(tile['device'],{})
    fresh={size:max(0,count-cache.get(size,0)) for size,count in needed.items()}
    pages=sum(count*((size+4095)//4096) for size,count in fresh.items())
    cached_calls=sum(needed.values())-sum(fresh.values())
    for size,count in needed.items():cache[size]=max(cache.get(size,0),count)
    return pages,cached_calls


def _serial_tile_specs(candidate):
    """Call-scoped exact template keys; offsets are labels, all work is retained.

    Candidate inputs are read-only during evaluation. Memoize a shared census's
    serialization only within this traversal; never retain identities across
    calls. Non-JSON Python inputs remain supported without template reuse.
    """
    caches,encoded_keys={},{}
    for index,tile in enumerate(candidate['tiles']):
        pin_state=_tile_pin_state(tile,caches)
        data=tile['data'];encoded=data['encoded']
        try:
            if id(encoded) not in encoded_keys:
                encoded_keys[id(encoded)]=json.dumps(encoded,sort_keys=True,separators=(',',':'),allow_nan=False)
            key=(encoded_keys[id(encoded)],json.dumps(
                [tile['device'],{k:v for k,v in data.items() if k!='encoded'},tile['profile'],pin_state],
                sort_keys=True,separators=(',',':'),allow_nan=False))
        except (TypeError,ValueError):key=index
        yield tile,pin_state,key


def _fixed_transfer_capacities(candidate):
    # The composed graph has one fixed transfer capacity per device. Keep its
    # existing semantics for hand-written candidates with changing capacities.
    first={}
    for tile in candidate['tiles']:
        p=tile['profile'];device=tile['device']
        capacities=tuple(p.get(direction+'_bytes_per_second') for direction in ('h2d','d2h'))
        if device in first and first[device]!=capacities:return False
        first[device]=capacities
    return True


def _serial_tile_reuse(candidate):
    return len(candidate['devices'])==1 and _fixed_transfer_capacities(candidate)


def trait_tiled_graph_chunks(candidate):
    """Number of chunks actually expanded by a runtime call, before scenarios.

    Only a single active device permits exact sequential decomposition. All
    multi-device tiles still enter the joint shared-resource graph.
    """
    if not _serial_tile_reuse(candidate):
        return sum(scan_chunk_count(t['data'],t['profile'])
                   for t in candidate['tiles'])
    seen=set();chunks=0
    for tile,_,key in _serial_tile_specs(candidate):
        if key in seen:continue
        seen.add(key)
        chunks+=scan_chunk_count(tile['data'],tile['profile'])
    return chunks


def torch_trait_tiled_runtime(candidate, *, host_serial_fraction, host_serial_policy='fluid',
                              pageable_arena_fresh_fraction=None, return_graph=False):
    """Solve the same finite graph, reusing fully drained single-device tiles.

    A tile completion joins every scan, release and durable writer node. With
    one device, those cuts have no live resource, queue or token interaction.
    Solving each distinct tile from time zero and summing is therefore the
    same schedule up to floating-point accumulation. Pinned-cache state is
    carried explicitly. Multiple devices retain one shared event scheduler,
    admitting the next tile in each chain only after its predecessor drains.
    ``return_graph`` always returns the original fully expanded graph.
    """
    shape=trait_tiled_shape(candidate)
    options=dict(host_serial_fraction=host_serial_fraction,host_serial_policy=host_serial_policy,
                 pageable_arena_fresh_fraction=pageable_arena_fresh_fraction)
    if return_graph or not _serial_tile_reuse(candidate) or len(candidate['tiles'])==1:
        stream=not return_graph and len(shape['devices'])>1 and _fixed_transfer_capacities(candidate)
        return _trait_tiled_runtime(candidate,return_graph=return_graph,_stream_tiles=stream,**options)
    solved={};reports=[];durations=[];handoffs=[];unpriced=set()
    for tile,pin_state,key in _serial_tile_specs(candidate):
        if key not in solved:
            width=tile['data']['traits_analyzed']
            single=dict(candidate,trait_block=width,
                        tiles=[dict(tile,trait_range=[0,width])])
            solved[key]=_trait_tiled_runtime(single,_pin_states=[pin_state],**options)
        result=solved[key]
        report=copy.deepcopy(result['tiles'][0]);report['trait_range']=list(tile['trait_range'])
        reports.append(report);durations.append(result['estimated_tile_seconds'])
        handoffs.append(result['handoffs']);unpriced.update(result['unpriced_terms'])
    combined=dict(result,estimated_tile_seconds=math.fsum(durations),tiles=reports,
        genotype_passes=len(reports),binary_payload_bytes=sum(r['binary_payload_bytes'] for r in reports),
        df_payload_bytes=sum(r['df_payload_bytes'] for r in reports),unpriced_terms=sorted(unpriced))
    combined['handoffs']=dict(result['handoffs'],
        possible_waits=sum(h['possible_waits'] for h in handoffs),
        blocked_waits=sum(h['blocked_waits'] for h in handoffs),
        extra_elapsed_service_seconds=math.fsum(h['extra_elapsed_service_seconds'] for h in handoffs))
    return combined


def _trait_tiled_runtime(candidate, *, host_serial_fraction, host_serial_policy='fluid',
                         pageable_arena_fresh_fraction=None, return_graph=False, _pin_states=None,
                         _stream_tiles=False):
    """One shared graph, including repeated scans, setup, writer credits and fsync."""
    shape = trait_tiled_shape(candidate)
    if isinstance(host_serial_fraction, bool) or not math.isfinite(host_serial_fraction) or not 0<=host_serial_fraction<=1:
        raise ValueError('Invalid host_serial_fraction')
    if host_serial_policy not in ('fluid','held-first','held-last'):
        raise ValueError('Invalid host_serial_policy')
    if pageable_arena_fresh_fraction is not None and (isinstance(pageable_arena_fresh_fraction,bool)
        or not math.isfinite(pageable_arena_fresh_fraction) or not 0<=pageable_arena_fresh_fraction<=1):
        raise ValueError('Invalid pageable_arena_fresh_fraction')
    g = ExecutionGraph()
    g.capacities = dict(candidate['shared_capacities'], host_serial=1.)
    storage = candidate.get('shared_storage_bytes_per_second')
    if storage is not None:
        if isinstance(storage,bool) or not math.isfinite(storage) or storage<=0:
            raise ValueError('Positive finite shared storage capacity required')
        g.capacities['storage'] = storage
    links = candidate.get('shared_links', [])
    for index, link in enumerate(links):
        if not link['devices'] or not set(link['devices']) <= set(shape['devices']):
            raise ValueError('Unknown shared-link devices')
        for direction in ['h2d','d2h']:
            capacity = link[direction+'_bytes_per_second']
            if isinstance(capacity,bool) or not math.isfinite(capacity) or capacity<=0:
                raise ValueError('Positive finite shared-link capacity required')
            g.capacities[f'link:{index}:{direction}'] = capacity
    done, caches, reports, unpriced = {}, {}, [], set()
    streams={device:[] for device in shape['devices']};templates={}
    specs=iter(_serial_tile_specs(candidate)) if _stream_tiles else None
    for index, tile in enumerate(candidate['tiles']):
        device, data, profile = tile['device'], tile['data'], tile['profile']
        if _stream_tiles:
            _,pin_state,key=next(specs)
            if key in templates:
                local,report,terms=templates[key]
                streams[device].append((index,local));unpriced.update(terms)
                report=copy.deepcopy(report);report['trait_range']=list(tile['trait_range'])
                if shape.get('partition_axis')=='variant':report['variant_range']=list(tile['variant_range'])
                reports.append(report)
                continue
        if pageable_arena_fresh_fraction is not None:
            if profile.get('pageable_host_service') is None:
                raise ValueError('Arena scenario requires independent pageable host service')
            profile=dict(profile,pageable_host_service=dict(profile['pageable_host_service'],
                arena_fresh_fraction=pageable_arena_fresh_fraction))
        work = _cached_scan_work(data, profile)
        if work.get('status') == 'zero_available_capacity':
            raise ValueError('No available scan capacity')
        blocks = copy.deepcopy(work['blocks'])
        for block in blocks:
            block['host_resources']['host_serial'] = block['host_resources']['cpu']*host_serial_fraction
            for direction in ['h2d','d2h']:
                resources = block.setdefault(direction+'_resources', {})
                for link_index, link in enumerate(links):
                    if device in link['devices']:
                        seconds = block[direction+'_seconds']
                        resources[f'link:{link_index}:{direction}'] = block[direction+'_bytes']/seconds if seconds else 0.
        q = profile['cpu_fraction']
        units = profile['process_units']
        # Full durable output must price page-cache copies AND storage, even
        # below the periodic writeback threshold; no optimistic close shortcut.
        service = {key:value if key=='storage_seconds_per_byte' else value/q
                   for key,value in profile['writeback_service'].items()}
        writer_work = _output_work(tile, candidate['output'])
        copy_service = writer_copy_cost(profile)
        writer = BinaryWriterSchedule(writer_work,
            copy_seconds_per_byte=copy_service['cpu_seconds_per_byte']/q,
            copy_seconds_per_call=copy_service['cpu_seconds_per_call']/q,
            zero_seconds_per_byte=units['bytearray_zero_bytes']/q,
            write_seconds_per_byte=1/profile['write_bytes_per_second'], fsync_seconds_per_array=profile['fsync_seconds'],
            append_seconds=profile['executor_cpu_seconds']/q, handoff_seconds=profile['executor_cpu_seconds']/q,
            cpu_fraction=q, write_capacity=profile['write_bytes_per_second'], writeback_service=service)
        local = torch_scan_schedule(blocks, depth=work['depth'], decode_workers=work['workers'],
            consumer=writer, return_graph=True,
            consumer_open_before_scan=candidate['output']['block_bytes'] is not None)
        # Allocation/copy/writer Python serialization is currently a supplied
        # host scenario. Native NumPy copies may release the GIL.
        for name, resources in local.demands.items():
            if name.startswith('writer:') and resources.get('cpu'):
                resources['host_serial'] = resources['cpu']*host_serial_fraction
            if storage is not None:
                resources['storage'] = resources.get('input',0.)+resources.get('output',0.)
        n,k,c = data['samples'],data['traits_analyzed'],data['covariates']
        pages,cached_calls = (pin_state if _stream_tiles else
            _tile_pin_state(tile,caches) if _pin_states is None else _pin_states[index])
        prep, prep_estimate = _prepare_graph(setup_work(n,k,c,reuse_observed_counts=True,
            input_contiguous=data.get('phenotype_c_contiguous'),covariate_columns=data.get('covariate_columns')),profile,host_serial_fraction,pages,cached_calls)
        cleanup=prep_estimate['cleanup']
        # Generator exhaustion synchronizes its streams and joins workers before
        # releasing the CPU phenotype. The caller then closes the writer. Earlier
        # queued writes can still overlap this destruction.
        close_seconds,close_deps=local.nodes['writer:close:start']
        cleanup_seconds=cleanup['cpu_seconds']/q
        cleaned=local.add('prepare:host_cleanup',cleanup_seconds,
            list(close_deps)+[f'release:{len(blocks)-1}'],
            {'cpu':q,'host_serial':cleanup['serial_cpu_seconds']/cleanup_seconds if cleanup_seconds else 0.})
        local.nodes['writer:close:start']=(close_seconds,(cleaned,))
        for resources in prep.demands.values():
            for direction in ['h2d','d2h']:
                demand = resources.pop('prep_'+direction,0.)
                resources[device+':'+direction] = demand
                g.capacities[device+':'+direction] = profile[direction+'_bytes_per_second']
                for link_index,link in enumerate(links):
                    if device in link['devices']:resources[f'link:{link_index}:{direction}'] = demand
        original_roots = [name for name,(_,deps) in local.nodes.items() if not deps]
        early = candidate['output']['block_bytes'] is not None
        prepared = local.compose(prep, 'prepare:', ['writer:open'] if early else [])
        for name,(seconds,deps) in list(local.nodes.items()):
            if name.startswith('prepare:') or (early and name=='writer:open'):continue
            if (early and 'writer:open' in deps) or (not early and name in original_roots):
                local.nodes[name] = (seconds, tuple(d for d in deps if not early or d!='writer:open')+(prepared,))
        terms=work['unpriced_terms']+work['allocator_unpriced_terms']+prep_estimate['unpriced_terms']
        unpriced.update(terms)
        reports.append(dict(device=device,trait_range=tile['trait_range'],blocks=len(blocks),
            binary_payload_bytes=writer.payload,df_payload_bytes=writer_work['df_payload_bytes'],
            writer_copy_bytes=writer.copied,writer_copy_calls=writer.copy_calls,writer_write_calls=writer.write_calls,
            setup=prep_estimate,pin_fresh_pages=pages,pin_cached_calls=cached_calls))
        if shape.get('partition_axis')=='variant':reports[-1]['variant_range']=tile['variant_range']
        if _stream_tiles:
            if host_serial_policy!='fluid':local=local.with_serial_sections(host_serial_policy)
            templates[key]=(local,reports[-1],terms);streams[device].append((index,local))
        else:
            prefix = f'tile{index}:'
            done[device] = g.compose(local,prefix,[done[device]] if device in done else [])
    if _stream_tiles:
        result=g.solve_chains(streams.values(),shared_tokens={'exclusive:host_serial'},prefix='tile')
        handoffs=dict(handoff_summary({}),**result['conditional_summary'])
    else:
        if host_serial_policy!='fluid':g = g.with_serial_sections(host_serial_policy)
        if return_graph:return g
        result = g.solve();handoffs=handoff_summary(result)
    unpriced.update(['API input QC/loading, shared covariate basis and phenotype mmap faults/casts are outside the priced tile-worker boundary',
        'df return/validation, tile and root manifests, directory fsync, QC collection and variant IDs',
        'setup intra-phase CPU/GPU/transfer ordering is aggregated; host serialization is a supplied scenario',
        'fresh pinned fixed-call/driver service and actual allocator reuse across tile transitions',
        'kernel early writeback, dirty throttling and page-cache DRAM contention',
        'thread lifecycle, tile transition cleanup, cold library initialization and loaded driver contention'])
    if storage is None:unpriced.add('input/output physical storage sharing requires an explicit aggregate storage capacity')
    return dict(status='development_trait_tiled_candidate',estimated_tile_seconds=result['seconds'],
        prediction_complete=False,runtime_prediction_validated=False,tiles=reports,
        genotype_passes=1 if shape.get('partition_axis')=='variant' else len(reports),binary_payload_bytes=sum(r['binary_payload_bytes'] for r in reports),
        df_payload_bytes=sum(r['df_payload_bytes'] for r in reports),
        handoffs=handoffs,unpriced_terms=sorted(unpriced),
        scope=('Disjoint per-device variant ranges with replicated phenotype setup and durable shard writers. ' if shape.get('partition_axis')=='variant' else 'Ordered per-device full genotype passes with borrowed synchronous staging and durable payload close. ')+
              'Declared aggregate CPU/DRAM/input/output and PCIe sharing; no cross-device result queue. '
              'Phase setup aggregation and per-device pinned-cache reuse are explicit scenarios, not measured bounds.')
