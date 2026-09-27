"""JAGWAS candidate admission with shared preprocessing and narrow ownership.

This is the native-int8, complete FP32 phenotype contract. Tensor and array
requests plus caller reserves are not a bound for RSS or CUDA reservation.
"""
from .decoder_work import native_reader_workspace,scan_chunk_count,require_regular_memory
from .mechanistic_plan import _integer
from .pinned_work import pinned_scan_work
from .reduced_output_work import jagwas_writer_work
from .setup_work import setup_work
from .tensor_memory import eager_scan_memory,eager_memory_plan
from .trait_tiling_model import trait_tiled_shape,reuse_trait_work


def jagwas_candidate_shape(candidate):
    shape=trait_tiled_shape(candidate,reduction='jagwas')
    output=candidate['output']
    if output['block_bytes'] is not None or output['store_beta'] is not False:
        raise ValueError('JAGWAS requires two-array indexed parts without beta or dense block coalescing')
    shape['queue_depth']=output['queue_depth'] if len(shape['devices'])>1 else 0
    shape['reduction']='jagwas'
    shape['covariate_columns']=candidate['tiles'][0]['data'].get('covariate_columns',shape['covariates'])
    return shape


def jagwas_candidate_memory(candidate, *, host_reserve_bytes=0, device_reserve_bytes=0,
                             device_memory_profiles=None):
    """Count full-retention staging before considering any survivor scenario."""
    require_regular_memory(candidate)
    shape=jagwas_candidate_shape(candidate)
    _integer('host_reserve_bytes',host_reserve_bytes,0)
    _integer('device_reserve_bytes',device_reserve_bytes,0)
    n,k,c,b,depth=(shape[key] for key in ['samples','traits','covariates','chunk_size','depth'])
    columns=candidate['tiles'][0]['data'].get('covariate_columns',c)
    setup=setup_work(n,k,c,covariate_columns=columns)
    # Input cast and shared residual output; blocked residual copy/download;
    # covariate construction arrays and shared FP32 Q; observed-count vector.
    # Counting disjoint lifetimes together is intentional and conservative.
    shared=dict(phenotype_cast_and_result=8*n*k,
        residual_host_blocks=8*n*setup['residual_block_traits'],
        covariate_basis_arrays=36*n*columns+8*columns*columns+24*columns,
        shared_covariate_basis=4*n*c,observed_counts=8*k)
    devices={};unknown=set(setup['unresolved_memory_terms']);largest=0
    for tile in candidate['tiles']:
        device,data,profile=tile['device'],tile['data'],tile['profile']
        m=data['markers'];rows=min(b,m);largest=max(largest,rows)
        active=min(profile['decode_workers'],depth,(m+b-1)//b)
        payload,replay=native_reader_workspace(data['encoded']['chunks'])
        pins=pinned_scan_work(n,b,k,depth,reduction='jagwas')
        # Depth finish workers may each own five arrays (17 B/variant), a
        # malformed int64 index and masks; queued tuples retain only 12 B.
        arrays=dict(design_host_block=4*n*setup['design_block_traits'],
            phenotype_df_and_scale_vectors=20*k,
            decoder_workspace=active*(rows*((n+3)//4)+payload+replay),
            pinned_buffers=pins['allocator_bytes'],
            finish_and_pending_arrays=depth*(17+8+4)*rows,
            producer_pending_tuple=12*rows)
        memory=eager_scan_memory(n,b,k,c,depth,reduction='jagwas')
        profile_memory=(device_memory_profiles or {}).get(device)
        if profile_memory is None:
            # Factor preparation is not in eager_scan_memory. Require its
            # tensor sum even when an exact vendor workspace query is absent.
            from .reduction_tensor_work import jagwas_tensor_work
            factor=jagwas_tensor_work(n,b,k,phase='prepare')
            gpu=max(memory['tensor_storage_budget'],4*n*k+factor['distinct_temporary_bytes'],
                setup['design_device_live_bytes_upper']+memory['persistent_factor_bytes'],
                setup['residual_device_live_bytes_upper'])
            unknown.add('per-device CUDA library/reduction workspace requires a device profile or reserve')
        else:
            plan=eager_memory_plan(n,b,k,c,depth,profile_memory,reduction='jagwas')
            gpu=plan['device_bytes'];unknown.update(plan['unresolved_memory_terms'])
            arrays['factor_host_workspace']=plan.get('factor_host_workspace_bytes',0)
        devices[device]=dict(host_arrays=arrays,host_bytes=sum(arrays.values()),
            device_bytes=gpu+device_reserve_bytes,pinned=pins,persistent_factor_bytes=memory['persistent_factor_bytes'])
    queue=shape['queue_depth']*12*largest
    consumer=2*12*largest
    writer=jagwas_writer_work(largest,largest)['writer_array_bytes_upper']
    host=sum(shared.values())+sum(row['host_bytes'] for row in devices.values())+queue+consumer+writer+host_reserve_bytes
    unknown.update(['input mappings, IDs, genotype index/record metadata and NPZ-part metadata require host reserve',
        'QC construction peaks, LAPACK workspace, Python objects and allocator retention require host reserve',
        'CUDA driver/allocator reservation and pinned-cache fragmentation require explicit reserves'])
    return dict(host_bytes=host,device_bytes={d:row['device_bytes'] for d,row in devices.items()},
        shared_host_arrays=shared,shared_host_bytes=sum(shared.values()),devices=devices,
        queued_result_bytes=queue,consumer_result_bytes=consumer,writer_array_bytes=writer,
        host_reserve_bytes=host_reserve_bytes,device_reserve_bytes=device_reserve_bytes,
        occupancy='all variants retained; unaffected by predicted output occupancy',
        shape=shape,unresolved_memory_terms=sorted(unknown),prediction_complete=False,
        scope='Conservative explicit request accounting for one shared preprocessing pass, independent factors, owned narrow queues and one indexed writer. Reserves are mandatory for unlisted allocations; this is not a complete memory guarantee.')


def jagwas_candidate_runtime(candidate,prices,*,preparation,occupancy,host_serial_fraction,
                              host_serial_policy='fluid',return_graph=False,max_source_chunks=10000):
    """Compose source scan and writer work with explicit preparation services.

    Preparation supplies shared_graph, device_graphs, per-device cleanup,
    finalize, dimensions [N,K,rank,columns,chunk,depth], and unpriced_terms. There is deliberately no
    default empty setup or dense-output price reuse. Those graphs must come
    from independent component services; association durations are not inputs.
    """
    import copy
    import math
    from .execution_graph import ExecutionGraph
    from .indexed_schedule import jagwas_variant_schedule
    from .mechanistic_torch import torch_scan_work
    from .reduced_output_work import jagwas_host_selection_service,jagwas_archive_service
    shape=jagwas_candidate_shape(candidate)
    _integer('max_source_chunks',max_source_chunks)
    if isinstance(host_serial_fraction,bool) or not isinstance(host_serial_fraction,(int,float)) or not math.isfinite(host_serial_fraction) or not 0<=host_serial_fraction<=1:
        raise ValueError('Explicit finite host serialization fraction required')
    if not isinstance(preparation,dict) or set(preparation)!={'dimensions','shared_graph','device_graphs','cleanup','finalize','unpriced_terms'}:
        raise ValueError('Explicit complete preparation-service contract required')
    columns=candidate['tiles'][0]['data'].get('covariate_columns',shape['covariates'])
    if preparation['dimensions']!=[shape['samples'],shape['traits'],shape['covariates'],columns,shape['chunk_size'],shape['depth']]:
        raise ValueError('Preparation dimensions differ from candidate')
    devices=shape['devices']
    if set(preparation['device_graphs'])!=set(devices) or set(preparation['cleanup'])!=set(devices):
        raise ValueError('Independent preparation and cleanup for each active device required')
    graphs=[preparation['shared_graph'],*preparation['device_graphs'].values()]
    if any(not isinstance(g,ExecutionGraph) or not g.nodes for g in graphs):
        raise ValueError('Nonempty shared and per-device preparation graphs required')
    if not isinstance(preparation['finalize'],list) or not preparation['finalize']:
        raise ValueError('Explicit final metadata/commit service required')
    if not isinstance(preparation['unpriced_terms'],list) or any(not isinstance(v,str) or not v for v in preparation['unpriced_terms']):
        raise ValueError('Explicit preparation uncertainty list required')
    if isinstance(occupancy,str):
        if occupancy not in ('empty','dense'):raise ValueError('Unknown output occupancy scenario')
    elif not isinstance(occupancy,(list,tuple)) or len(occupancy)!=len(candidate['tiles']):
        raise ValueError('One explicit retained-count sequence per shard required')
    count=sum(scan_chunk_count(tile['data'],tile['profile']) for tile in candidate['tiles'])
    if count>max_source_chunks:raise ValueError('Candidate exceeds bounded source-chunk expansion')
    caps=dict(candidate['shared_capacities'],host_serial=1.)
    def positive(value,name,zero=False):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (not zero and not value):
            raise ValueError('Invalid independent '+name)
        return value
    storage=candidate.get('shared_storage_bytes_per_second')
    if storage is not None:caps['storage']=positive(storage,'shared storage capacity')
    links=candidate.get('shared_links',[])
    for index,link in enumerate(links):
        members=link['devices']
        if not members or len(set(members))!=len(members) or not set(members)<=set(devices):
            raise ValueError('Unique known shared-link devices required')
        for direction in ['h2d','d2h']:
            caps[f'link:{index}:{direction}']=positive(link[direction+'_bytes_per_second'],'link capacity')
    for tile in candidate['tiles']:
        for direction in ['h2d','d2h']:
            caps[tile['device']+':'+direction]=positive(tile['profile'][direction+'_bytes_per_second'],'device transfer capacity')
    for graph in graphs:
        available=set(caps)|set(graph.capacities)
        if any(set(demands)-available for demands in graph.demands.values()):
            raise ValueError('Preparation contains an unbounded/unnamed resource')
    shards=[];reports=[];unknown=set(preparation['unpriced_terms'])
    for index,tile in enumerate(candidate['tiles']):
        data,profile,device=tile['data'],tile['profile'],tile['device']
        work=torch_scan_work(data,profile)
        if work.get('status')=='zero_available_capacity':return work
        if work['result_ownership']!='owned':raise ValueError('Owned JAGWAS result contract required')
        blocks=copy.deepcopy(work['blocks']);m=data['markers'];b=profile['chunk_markers'];outputs=[]
        if not isinstance(occupancy,str) and len(occupancy[index])!=len(blocks):
            raise ValueError('One exact retained count per source chunk required')
        parts=payload=retained=0;q=profile['cpu_fraction']
        for chunk,block in enumerate(blocks):
            rows=block['markers']
            keep=(0 if occupancy=='empty' else rows if occupancy=='dense' else occupancy[index][chunk])
            _integer('retained variants',keep,0)
            output=jagwas_writer_work(rows,keep,fsync=candidate['output']['fsync'])
            selection=jagwas_host_selection_service(output,prices['prices'],cpu_fraction=q,
                dram_bytes_per_second=profile['shared_dram_bytes_per_second'],host_serial_fraction=host_serial_fraction)
            writer=jagwas_archive_service(output,prices['archive'],profile,host_serial_fraction=host_serial_fraction)
            if block.get('consumer_seconds',0.):raise ValueError('Unknown scan consumer service cannot be replaced')
            release=block.pop('discard_seconds',0.)
            if release:
                # Owned scan-to-discard frees move to the single consumer.
                # Exact prior-reference destruction around the next get is
                # separately disclosed as an ordering uncertainty.
                (writer if keep else selection).append(dict(seconds=release,resources=dict(cpu=q,host_serial=q)))
            outputs.append([dict(cells=rows,retained=keep,selection=selection,writer=writer)])
            parts+=bool(keep);payload+=output['part']['file_bytes'];retained+=keep
            block['host_resources']['host_serial']=block['host_resources']['cpu']*host_serial_fraction
            for direction in ['h2d','d2h']:
                rate=block[direction+'_bytes']/block[direction+'_seconds'] if block[direction+'_seconds'] else 0.
                resource=block.setdefault(direction+'_resources',{})
                resource[device+':'+direction]=rate
                for li,link in enumerate(links):
                    if device in link['devices']:resource[f'link:{li}:{direction}']=rate
        shards.append(dict(device=device,backend='host',blocks=blocks,outputs=outputs,depth=work['depth'],
            decode_workers=work['workers'],prepare=preparation['device_graphs'][device],cleanup=preparation['cleanup'][device]))
        reports.append(dict(device=device,variant_range=tile['variant_range'],source_chunks=len(blocks),
            indexed_part_bytes=payload,parts=parts,retained_variants=retained,
            d2h_bytes=sum(block['d2h_bytes'] for block in blocks)))
        unknown.update(work['unpriced_terms']+work['allocator_unpriced_terms'])
    queue=None
    if len(devices)>1:
        q=min(tile['profile']['cpu_fraction'] for tile in candidate['tiles'])
        queue={kind:[dict(seconds=positive(prices['queue_cpu_seconds'][kind],'queue CPU',True)/q,
            resources=dict(cpu=q,host_serial=q))] for kind in ['put','get']}
    graph=jagwas_variant_schedule(shards,shared_prepare=preparation['shared_graph'],
        queue_depth=shape['queue_depth'],shared_capacities=caps,queue_service=queue,
        finalize=preparation['finalize'],return_graph=True,max_source_chunks=max_source_chunks)
    # Shared physical resources also cover caller-supplied setup and final
    # graphs. Transfers keep device and every applicable upstream link demand.
    for demands in graph.demands.values():
        for li,link in enumerate(links):
            for direction in ['h2d','d2h']:
                key=f'link:{li}:{direction}'
                rate=sum(demands.get(device+':'+direction,0.) for device in link['devices'])
                if key in demands and demands[key]!=rate:
                    raise ValueError('Preparation link demand differs from its device transfers')
                if rate:demands[key]=rate
        if storage is not None:demands['storage']=demands.get('input',0.)+demands.get('output',0.)
    if host_serial_policy!='fluid':graph=graph.with_serial_sections(host_serial_policy)
    if return_graph:return graph
    from .resource_balance import resource_balance
    solved=graph.solve()
    unknown.update(['fixed primitive CPU/serial service transfer to loaded candidate context',
        'previous-result destruction placement around next queue get',
        'indexed-writer allocation/metadata and dirty-page overlap beyond supplied services',
        'blocked queue wakeup, timeout retries and worker lifecycle'])
    from .reduction_tensor_work import jagwas_factor_arithmetic_work
    factor=jagwas_factor_arithmetic_work(shape['samples'],shape['traits'])
    from .executor_timing import SETUP_SCAN_WRITE_METRIC, SETUP_SCAN_WRITE_BOUNDARY
    return dict(status='development_jagwas_candidate',estimated_seconds=solved['seconds'],
        observed_metric=SETUP_SCAN_WRITE_METRIC, timing_boundary=SETUP_SCAN_WRITE_BOUNDARY,
        resource_balance=resource_balance(graph, solved),
        factor_preparation=dict(instances=len(devices),per_device=factor,
            total_h2d_bytes=len(devices)*factor['h2d_bytes'],
            total_fp64_multiply_add_flops=len(devices)*factor['fp64_multiply_add_flops']),
        shards=reports,parts=sum(row['parts'] for row in reports),
        indexed_part_bytes=sum(row['indexed_part_bytes'] for row in reports),
        retained_variants=sum(row['retained_variants'] for row in reports),
        shared_preprocessing_passes=1,queue_depth=shape['queue_depth'],genotype_passes=1,
        unpriced_terms=sorted(unknown),prediction_complete=False,automatic_selection_ready=False,
        scope='One shared preparation, per-device factor/design graphs, native scan and projection, owned bounded result queue, one indexed consumer and explicit final commit. Conditional component services only; neither memory admission nor autotune qualification is implied.')


@reuse_trait_work
def detailed_jagwas_plan(candidates,prices,*,preparations=None,preparation_factory=None,occupancy_scenarios,host_scenarios,
                         cpu_workers,host_memory_bytes,device_memory_bytes,device_memory_profiles=None,
                         host_reserve_bytes=0,device_reserve_bytes=0,max_candidates=64,
                         max_scenario_evaluations=256,max_source_chunks=10000):
    """Minimize max_s T_graph over bounded, memory-admitted JAGWAS candidates.

    Scenarios are explicit modeling assumptions, not confidence bounds. This
    development plan does not bypass the public reduction-autotuning gate.
    preparation_factory(index, candidate, host_name, host) can supply preparation
    priced for that scenario; fixed preparations explicitly remain fixed while
    scan/writer host assumptions vary. The factory must not mutate candidates.
    """
    import math
    from .mechanistic_plan import _device
    for name,value in [('cpu_workers',cpu_workers),('host_memory_bytes',host_memory_bytes),
        ('max_candidates',max_candidates),('max_scenario_evaluations',max_scenario_evaluations),
        ('max_source_chunks',max_source_chunks)]:_integer(name,value)
    for name,value in [('host_reserve_bytes',host_reserve_bytes),('device_reserve_bytes',device_reserve_bytes)]:
        _integer(name,value,0)
    if not isinstance(candidates,(list,tuple)) or not candidates or len(candidates)>max_candidates:
        raise ValueError('Nonempty bounded JAGWAS candidate list required')
    if preparation_factory is not None:
        if preparations is not None or not callable(preparation_factory):
            raise ValueError('Supply either fixed preparations or a callable preparation_factory')
    elif not isinstance(preparations,(list,tuple)) or len(preparations)!=len(candidates):
        raise ValueError('One preparation-service contract per candidate required')
    if not isinstance(occupancy_scenarios,dict) or not occupancy_scenarios or any(
        not isinstance(name,str) or not name or value not in ('empty','dense') for name,value in occupancy_scenarios.items()):
        raise ValueError('Named explicit empty/dense survivor scenarios required')
    if not isinstance(host_scenarios,dict) or not host_scenarios:
        raise ValueError('Named explicit host scenarios required')
    for name,host in host_scenarios.items():
        if not isinstance(name,str) or not name or not isinstance(host,dict) or set(host)-{'host_serial_policy'}!={'host_serial_fraction'}:
            raise ValueError('Host scenario requires serialization fraction and optional policy')
        value=host['host_serial_fraction']
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError('Finite bounded host serialization fraction required')
        if host.get('host_serial_policy','fluid') not in ('fluid','held-first','held-last'):
            raise ValueError('Unknown host serialization policy')
    for device,capacity in device_memory_bytes.items():
        _device(device);_integer('device capacity',capacity)
    shapes=[jagwas_candidate_shape(candidate) for candidate in candidates]
    identity={key:shapes[0][key] for key in ['samples','markers','traits','covariates','covariate_columns','input_path']}
    if any(any(shape[key]!=value for key,value in identity.items()) for shape in shapes):
        raise ValueError('Candidates must represent the same statistical workload')
    feasible=[];rejected=[]
    for index,(candidate,shape) in enumerate(zip(candidates,shapes)):
        if shape['reader_workers']>cpu_workers:
            rejected.append(dict(candidate_index=index,reason='shared_reader_budget'));continue
        if not set(shape['devices'])<=set(device_memory_bytes):raise ValueError('Missing device memory capacity')
        memory=jagwas_candidate_memory(candidate,host_reserve_bytes=host_reserve_bytes,
            device_reserve_bytes=device_reserve_bytes,device_memory_profiles=device_memory_profiles)
        reason=('host_memory' if memory['host_bytes']>host_memory_bytes else
            'device_memory' if any(memory['device_bytes'][d]>device_memory_bytes[d] for d in shape['devices']) else None)
        if reason:rejected.append(dict(candidate_index=index,reason=reason,memory=memory))
        else:feasible.append((index,candidate,shape,memory))
    if len(feasible)*len(occupancy_scenarios)*len(host_scenarios)>max_scenario_evaluations:
        raise ValueError('JAGWAS search exceeds max_scenario_evaluations')
    if not feasible:raise ValueError('No feasible JAGWAS candidate: '+str(rejected))
    for _,candidate,shape,_ in feasible:
        count=sum(scan_chunk_count(tile['data'],tile['profile']) for tile in candidate['tiles'])
        if count>max_source_chunks:
            raise ValueError('Candidate exceeds bounded source-chunk expansion')
    rows=[];unpriced=set();memory_terms=set()
    for index,candidate,shape,memory in feasible:
        # Build once per host scenario, only after reader/memory admission and
        # the complete evaluation budget check. Reuse across survivor counts.
        prepared={name:(preparation_factory(index,candidate,name,dict(host))
                        if preparation_factory is not None else preparations[index])
                  for name,host in host_scenarios.items()}
        scenarios={}
        for occupancy_name,occupancy in occupancy_scenarios.items():
            for host_name,host in host_scenarios.items():
                # Tuple labels avoid collisions when caller names contain ':'.
                label=(occupancy_name,host_name)
                result=jagwas_candidate_runtime(candidate,prices,preparation=prepared[host_name],
                    occupancy=occupancy,max_source_chunks=max_source_chunks,**host)
                seconds=result.get('estimated_seconds')
                if isinstance(seconds,bool) or not isinstance(seconds,(int,float)) or not math.isfinite(seconds) or seconds<=0:
                    raise ValueError('No finite positive JAGWAS candidate score')
                scenarios[label]=result;unpriced.update(result['unpriced_terms'])
        score=max(result['estimated_seconds'] for result in scenarios.values())
        memory_terms.update(memory['unresolved_memory_terms'])
        rows.append(dict(candidate_index=index,devices=shape['devices'],chunk_size=shape['chunk_size'],
            memory=memory,scenarios=[dict(occupancy=label[0],host=label[1],estimate=result) for label,result in scenarios.items()],
            worst_supplied_scenario_seconds=score,api_kwargs=dict(genotype=shape['input_path'],
                genotype_format='pgen',pgen_mode='hardcall',device=shape['devices'][0],variant_devices=shape['devices'],
                chunk_size=shape['chunk_size'],reader_workers=shape['reader_workers'],prefetch_chunks=shape['depth'],
                compute_dtype='float32',reduce='jagwas',sumstats_format='binary',sumstats_fsync=True,
                sumstats_queue_depth=candidate['output']['queue_depth']),
            required_environment=dict(TORCHGWAS_PGEN_BACKEND='native',TORCHGWAS_PGEN_PACKED='0',
                TORCHGWAS_NATIVE_STATS='0',TORCHGWAS_SCAN_PROFILE='0',
                TORCHGWAS_BLOCKING_EVENTS='1' if shape['blocking_events'] else '0')))
    rows.sort(key=lambda row:(row['worst_supplied_scenario_seconds'],len(row['devices']),row['candidate_index']))
    return dict(status='development_detailed_jagwas_plan',selected=rows[0],candidates=rows,rejected=rejected,
        workload=identity,candidates_evaluated=len(candidates),candidates_feasible=len(rows),
        preparation_policy='per_host_scenario' if preparation_factory is not None else 'fixed_across_host_scenarios',
        occupancy_scenarios=dict(occupancy_scenarios),host_scenarios=dict(host_scenarios),
        unpriced_terms=sorted(unpriced),unresolved_memory_terms=sorted(memory_terms),
        prediction_complete=False,automatic_selection_ready=False,selection_validated=False,
        objective='minimize worst explicitly supplied survivor/resource scenario after full-retention memory admission',
        scope='Bounded optimization over the shared analytical execution graph. Exact full-file PGEN censuses and caller-supplied component preparation services; no trait tiling of the joint statistic, fitted GWAS durations or public autotune authorization.')
