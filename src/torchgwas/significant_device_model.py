"""Device-significant candidates in the shared source/resource calculator."""
import copy
import json
from .mechanistic_plan import _integer
from .trait_tiling_model import trait_tiled_shape,_tile_memory
from .decoder_work import native_reader_workspace
from .reduced_output_work import significant_execution_layout
from .nonzero_memory import device_selection_memory


def significant_device_memory(candidate, *, device_memory_profiles,
                              host_reserve_bytes=0,device_reserve_bytes=0):
    """Worst-case live/queued pairs, independent of predicted significance.

    Reuse native setup/decoder/eager memory with the synchronous selected-result
    contract. The global queue is counted once. CUB requests use exact installed
    allocation censuses for each distinct full/tail source selection extent.
    """
    shape=trait_tiled_shape(candidate,reduction='device_significant')
    _integer('host_reserve_bytes',host_reserve_bytes,0)
    _integer('device_reserve_bytes',device_reserve_bytes,0)
    if candidate['output']['block_bytes'] is not None:
        raise ValueError('Indexed selected parts do not use dense block coalescing')
    if set(device_memory_profiles)!=set(shape['devices']):
        raise ValueError('One explicit memory profile per active device required')
    layout=significant_execution_layout(shape['traits'],shape['trait_block'],shape['devices'],
        shape['reader_workers'],queue_depth=candidate['output']['queue_depth'])
    host_peaks={};gpu_peaks={};selection_peaks={};pins={};reports=[];unresolved=set();largest=0
    for tile in candidate['tiles']:
        data,profile,device=tile['data'],tile['profile'],tile['device']
        n,m,k,c=[data[key] for key in ['samples','markers','traits_analyzed','covariates']]
        b,depth=profile['chunk_markers'],profile['depth']
        memory_profile=device_memory_profiles[device]
        if not all(key in memory_profile for key in ['nonzero_workspace_census','nonzero_workspace_context']):
            raise ValueError('Exact installed nonzero workspace census and current context required')
        rows={min(b,m)}
        if m%b:rows.add(m%b)
        selection=[device_selection_memory(n,count,k,
            census=memory_profile['nonzero_workspace_census'],context=memory_profile['nonzero_workspace_context']) for count in sorted(rows)]
        cells=max(row['maximum_selection_cells'] for row in selection)
        workspace=max(row['selection_gpu_bytes'] for row in selection)
        largest=max(largest,cells)
        payload,replay=native_reader_workspace(data['encoded']['chunks'])
        host,gpu,needed,unknown=_tile_memory(n,m,k,c,data.get('covariate_columns',c),b,depth,
            profile['decode_workers'],payload,None,memory_profile,replay,reduction='device_significant')
        host_peaks[device]=max(host_peaks.get(device,0),host)
        # Adding to the conservative setup/scan maximum is deliberately safe
        # even though selection runs after design preparation has completed.
        gpu_peaks[device]=max(gpu_peaks.get(device,0),gpu+workspace)
        cache=pins.setdefault(device,{})
        for size,count in needed.items():cache[size]=max(cache.get(size,0),count)
        # Current/previous selected host payloads plus one owned int64 trait-
        # index copy rebased in place in api.emit. Status QC and threshold use
        # explicit bounded arrays; metadata/allocator effects remain reserves.
        local=2*(28+8)*cells+13*min(b,m)+25*n+16
        selection_peaks[device]=max(selection_peaks.get(device,0),local)
        unresolved.update(unknown)
        reports.append(dict(device=device,trait_range=tile['trait_range'],selection=selection,
            native_device_bytes=gpu,selection_gpu_bytes=workspace,producer_selected_host_bytes=local))
    pin_bytes={d:sum(size*count for size,count in sizes.items()) for d,sizes in pins.items()}
    count_pins=4*len(shape['devices'])
    queued=layout['queue_depth']*28*largest
    consumer=(2*28+4)*largest
    archive=2*min(8*largest,16<<20)
    critical=8*(shape['samples']+1)
    total=sum(host_peaks.values())+sum(pin_bytes.values())+sum(selection_peaks.values())
    total+=count_pins+queued+consumer+archive+critical+host_reserve_bytes
    unresolved.discard('Device selection intermediates, CUB scratch and owned host payloads require the selection ledger')
    unresolved.update(['input mmap residency, genotype/variant/part metadata and QC require host reserve',
        'Python objects, archive headers, allocator fragmentation and retained CPU pages require host reserve',
        'CUDA driver/allocator retention and pinned driver backing require explicit reserves'])
    return dict(host_bytes=total,device_bytes={d:value+device_reserve_bytes for d,value in gpu_peaks.items()},
        host_active_bytes_by_device=host_peaks,pinned_cache_bytes_by_device=pin_bytes,
        producer_selection_bytes_by_device=selection_peaks,count_pinned_allocator_bytes=count_pins,
        queued_result_bytes=queued,consumer_result_bytes=consumer,archive_buffer_bytes=archive,
        shared_critical_table_bytes=critical,maximum_selection_cells=largest,tiles=reports,
        host_reserve_bytes=host_reserve_bytes,device_reserve_bytes=device_reserve_bytes,
        unresolved_memory_terms=sorted(unresolved),layout=layout,
        occupancy='all variant-trait pairs retained in each bounded selection block',
        scope='Shared conservative source memory admission with exact installed nonzero requests. No dense pinned output ring; selected payloads and one global bounded queue are fully charged. Explicit reserves cover non-tensor/allocator effects; no duration or fitted memory multiplier is used.')


def _critical_prepare(prep,setup,profile,prices,serial):
    """Insert the cached-table rounding/upload before native design creation."""
    from .significant_host_model import _positive
    n=setup['samples'];q=profile['cpu_fraction']
    cpu=_positive('critical round call',prices['round_call_cpu_seconds'],True)
    cpu+=n*_positive('critical round row',prices['round_row_cpu_seconds'],True)
    # astype32, astype64, compare, nextafter, where and final concatenate.
    traffic=70*n+8;seconds=max(cpu/q,traffic/profile['shared_dram_bytes_per_second'])
    design=next(i for i,p in enumerate(setup['phases']) if p['phase']=='design_common')
    node=f'phase:{design}';duration,deps=prep.nodes[node]
    rounded=prep.add('critical:round',seconds,deps,dict(cpu=cpu/seconds if seconds else 0.,
        host_serial=serial*cpu/seconds if seconds else 0.,dram=traffic/seconds if seconds else 0.))
    before=_positive('critical transfer before',prices['copy']['before_cpu_seconds'],True)
    submitted=prep.add('critical:submit',before/q,[rounded],dict(cpu=q,host_serial=q*serial))
    transfer=prices['transfer'];size=4*(n+1)
    elapsed=_positive('critical H2D latency',transfer['latency_seconds'],True)+size/_positive('critical H2D capacity',transfer['bytes_per_second'])
    copied=prep.add('critical:h2d',elapsed,[submitted],dict(prep_h2d=size/elapsed))
    wait=prices['wait_cpu_fraction']
    if isinstance(wait,bool) or not isinstance(wait,(int,float)) or not 0<=wait<=1:
        raise ValueError('Explicit bounded critical wait CPU scenario required')
    if wait:prep.resource_waits['critical:wait']=dict(after=submitted,until=copied,resources=dict(cpu=q*wait))
    after=_positive('critical transfer after',prices['copy']['after_cpu_seconds'],True)
    ready=prep.add('critical:ready',after/q,[copied],dict(cpu=q,host_serial=q*serial))
    prep.nodes[node]=(duration,(ready,))
    return dict(h2d_bytes=size,round_logical_dram_bytes=traffic,host_cpu_seconds=cpu+before+after)


def _emit_rebase(model,counts,prices,profile,serial):
    """api.emit copies indices once, then rebases that owned copy in place."""
    from .significant_host_model import _positive
    graph=model['graph'];q=profile['cpu_fraction']
    for index,count in enumerate(counts):
        previous=model['yield_nodes'][index]
        for name in ['index_cast','inplace_index_add']:
            if name not in prices:raise ValueError('Missing independent emit primitive: '+name)
            price=prices[name]
            cpu=_positive(name+' call',price['call_cpu_seconds'],True)+count*_positive(name+' row',price['row_cpu_seconds'],True)
            traffic=16*count;seconds=max(cpu/q,traffic/profile['shared_dram_bytes_per_second'])
            previous=graph.add(f'emit:{index}:{name}',seconds,[previous],
                dict(cpu=cpu/seconds if seconds else 0.,host_serial=serial*cpu/seconds if seconds else 0.,dram=traffic/seconds if seconds else 0.))
        model['yield_nodes'][index]=previous
        resume=model['resume_nodes'][index];duration,deps=graph.nodes[resume]
        graph.nodes[resume]=(duration,(*deps,previous))


def significant_device_runtime(candidate, prices, *, selection_services, host_serial_fraction,
                               host_serial_policy='fluid',return_graph=False,
                               max_source_chunks=10000,max_selection_blocks=100000):
    """Compose exact source extents with independently supplied primitive services.

    selection_services[tile][chunk] contains retained_per_block and the inputs
    to device_selection_graph (without CPU/capacity overrides). Counts are
    declared scenarios, never inferred from p thresholds or observed GWAS time.
    Equal source/service records are reused within this call only.
    """
    from .mechanistic_torch import torch_scan_work
    from .trait_tiling_model import _prepare_graph,_tile_pin_state
    from .setup_work import setup_work
    from .device_significance_work import device_significant_tensor_work
    from .device_significance_service import device_selection_graph
    from .significant_host_model import _positive,_writer_service
    from .significant_host_work import indexed_part_work
    from .indexed_schedule import significant_trait_schedule
    from .resource_balance import resource_balance
    from .reduced_output_work import significant_output_work
    shape=trait_tiled_shape(candidate,reduction='device_significant')
    _integer('max_source_chunks',max_source_chunks);_integer('max_selection_blocks',max_selection_blocks)
    if candidate['output']['block_bytes'] is not None:raise ValueError('Indexed parts do not use dense coalescing')
    if isinstance(host_serial_fraction,bool) or not 0<=host_serial_fraction<=1:raise ValueError('Bounded host serialization required')
    if not isinstance(selection_services,(list,tuple)) or len(selection_services)!=len(candidate['tiles']):
        raise ValueError('One independent selection-service sequence per tile required')
    layout=significant_execution_layout(shape['traits'],shape['trait_block'],shape['devices'],shape['reader_workers'],
        queue_depth=candidate['output']['queue_depth'])
    chunks=sum((t['data']['markers']+t['profile']['chunk_markers']-1)//t['profile']['chunk_markers'] for t in candidate['tiles'])
    if chunks>max_source_chunks:raise ValueError('Device candidate exceeds bounded graph expansion')
    caps=dict(candidate['shared_capacities'],host_serial=1.);links=candidate.get('shared_links',[])
    storage=candidate.get('shared_storage_bytes_per_second')
    if storage is not None:caps['storage']=_positive('storage capacity',storage)
    for index,link in enumerate(links):
        if not set(link['devices'])<=set(shape['devices']):raise ValueError('Unknown link device')
        for direction in ['h2d','d2h']:caps[f'link:{index}:{direction}']=_positive('shared link',link[direction+'_bytes_per_second'])
    for tile in candidate['tiles']:
        for direction in ['h2d','d2h']:
            key=tile['device']+':'+direction;rate=_positive(key,tile['profile'][direction+'_bytes_per_second'])
            if key in caps and caps[key]!=rate:raise ValueError('Conflicting per-device transfer capacity')
            caps[key]=rate
    def names(device,direction):
        return [device+':'+direction]+[f'link:{i}:{direction}' for i,link in enumerate(links) if device in link['devices']]
    tiles=[];reports=[];pins={};cache={};unresolved=set();selection_count=0
    expected={'retained_per_block','host_prices','operation_services','nonzero_services','transfer_prices','yield_cpu_seconds','wait_cpu_fraction'}
    for index,tile in enumerate(candidate['tiles']):
        data,profile,device=tile['data'],tile['profile'],tile['device']
        n,m,k,c=(data[key] for key in ['samples','markers','traits_analyzed','covariates'])
        work=torch_scan_work(data,profile);blocks=copy.deepcopy(work['blocks'])
        if len(selection_services[index])!=len(blocks):raise ValueError('One selection service per source chunk required')
        models=[];outputs=[];parts=payload=selected=counts_bytes=0
        for chunk,block in enumerate(blocks):
            row=selection_services[index][chunk]
            if set(row)!=expected:raise ValueError('Explicit source selection occupancy and independent primitive services required')
            rows=min(profile['chunk_markers'],m-chunk*profile['chunk_markers'])
            count=significant_output_work(n,rows,k,rows,backend='device',include_blocks=False)['selection_blocks']
            selection_count+=count
            if selection_count>max_selection_blocks:raise ValueError('Device candidate exceeds bounded selection expansion')
            key=json.dumps([n,rows,k,row,device,profile['cpu_fraction']],sort_keys=True,allow_nan=False)
            if key not in cache:
                source=device_significant_tensor_work(n,rows,k,row['retained_per_block'])
                service=copy.deepcopy({name:value for name,value in row.items() if name!='retained_per_block'})
                for transfer in service['transfer_prices'].values():
                    if transfer['resources']!=['d2h']:raise ValueError('Selector transfers must declare their local d2h resource')
                    transfer['resources']=names(device,'d2h')
                model=device_selection_graph(source,**service,cpu_fraction=profile['cpu_fraction'],
                    host_serial_fraction=host_serial_fraction,capacities=caps)
                # Ordinary GPU steps cannot silently introduce an unshared
                # resource. Dedicated GPU stream ordering is already explicit.
                if any(set(r)-set(caps) for r in model['graph'].demands.values()):raise ValueError('Undeclared selector resource capacity')
                cache[key]=(source,model)
            source,template=cache[key];model=copy.deepcopy(template)
            _emit_rebase(model,row['retained_per_block'],prices['emit'],profile,host_serial_fraction)
            group=[]
            for selected_block in source['blocks']:
                count=selected_block['retained'];part=indexed_part_work(count,store_beta=candidate['output']['store_beta'])
                group.append(dict(cells=selected_block['cells'],retained=count,selection=[],
                    writer=_writer_service(part,prices['archive'],profile,host_serial_fraction)))
                parts+=bool(count);payload+=part['file_bytes'];selected+=count
            models.append(model);outputs.append(group);counts_bytes+=model['nonzero_count_d2h_bytes']
            block['host_resources']['host_serial']=profile['cpu_fraction']*host_serial_fraction
            for direction in ['h2d','d2h']:
                seconds=block[direction+'_seconds'];rate=block[direction+'_bytes']/seconds if seconds else 0.
                block[direction+'_resources']={name:rate for name in names(device,direction)}
        pages,cached=_tile_pin_state(tile,pins)
        setup=setup_work(n,k,c,reuse_observed_counts=True,input_contiguous=data.get('phenotype_c_contiguous'),covariate_columns=data.get('covariate_columns'))
        prep,estimate=_prepare_graph(setup,profile,host_serial_fraction,pages,cached)
        critical=_critical_prepare(prep,setup,profile,prices['critical'][device],host_serial_fraction)
        for resources in prep.demands.values():
            for direction in ['h2d','d2h']:
                rate=resources.pop('prep_'+direction,0.)
                resources.update({name:rate for name in names(device,direction)})
        clean=estimate['cleanup'];q=profile['cpu_fraction'];seconds=clean['cpu_seconds']/q
        cleanup=[dict(seconds=seconds,resources=dict(cpu=q,host_serial=clean['serial_cpu_seconds']/seconds if seconds else 0.))]
        tiles.append(dict(device=device,backend='device',blocks=blocks,outputs=outputs,selection_graphs=models,
            depth=work['depth'],decode_workers=work['workers'],cleanup=cleanup,prepare=prep))
        reports.append(dict(device=device,trait_range=tile['trait_range'],source_chunks=len(blocks),parts=parts,
            indexed_part_bytes=payload,status_d2h_bytes=m,selected_payload_d2h_bytes=28*selected,
            nonzero_count_d2h_bytes=counts_bytes,critical_prepare=critical,setup=estimate,pin_fresh_pages=pages,pin_cached_calls=cached))
        unresolved.update(work['unpriced_terms']+work['allocator_unpriced_terms']+estimate['unpriced_terms'])
    queue=None
    if len(shape['devices'])>1:
        q=min(t['profile']['cpu_fraction'] for t in candidate['tiles'])
        queue={name:[dict(seconds=_positive('queue '+name,prices['queue_cpu_seconds'][name],True)/q,
            resources=dict(cpu=q,host_serial=q))] for name in ['put','get']}
    graph=significant_trait_schedule(tiles,queue_depth=layout['queue_depth'],shared_capacities=caps,
        queue_service=queue,finalize=[],return_graph=True,max_source_chunks=max_source_chunks,max_selection_blocks=max_selection_blocks)
    if storage is not None:
        for resources in graph.demands.values():resources['storage']=resources.get('input',0.)+resources.get('output',0.)
    if host_serial_policy!='fluid':graph=graph.with_serial_sections(host_serial_policy)
    if return_graph:return graph
    solved=graph.solve()
    unresolved.discard('Device-selected tensor operations, count barriers and owned payloads require the indexed selection graph')
    unresolved.update(['Independent selector services require held-out transfer and candidate-ranking validation',
        'Shared FP64 critical table, API QC and final metadata/directory publication',
        'Selected output allocation/release beyond the declared primitive services',
        'NPZ extent allocation, dirty-page pressure and result-queue wakeup/timeout overhead'])
    return dict(status='development_significant_device_candidate',estimated_tile_seconds=solved['seconds'],
        resource_balance=resource_balance(graph,solved),tiles=reports,layout=layout,
        indexed_part_bytes=sum(r['indexed_part_bytes'] for r in reports),parts=sum(r['parts'] for r in reports),
        unpriced_terms=sorted(unresolved),prediction_complete=False,selection_validated=False,
        scope='Conditional source setup, native scan/status, device selection, trait-index rebasing, one bounded queue and durable indexed writer. Exact supplied occupancy; no measured GWAS duration or selectivity fit.')
