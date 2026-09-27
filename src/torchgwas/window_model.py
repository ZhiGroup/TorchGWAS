"""Finite prepared source windows with the existing output/resource models.

This is a bounded graph construction, not a live executor checkpoint or an
extrapolation to the rest of a job. Source/service validity belongs to callers.
"""
from copy import deepcopy
from pathlib import Path
import json
import math
import time

from .binary_output_work import binary_output_memory
from .binary_schedule import BinaryWriterSchedule
from .calibration_cache import _digest,_json
from .execution_graph import ExecutionGraph,torch_scan_schedule
from .gpu_identity import canonical_cuda_device
from .mechanistic_torch import torch_scan_header_work
from .output_write_work import writer_copy_cost
from .pgen_work_bounds import WINDOW_KIND
from .planning_session import planning_work_scope
from .resource_balance import resource_balance
from .trait_tiling_model import _output_work


def _integer(value,name,minimum=1):
    if type(value) is not int or value<minimum:raise ValueError('Invalid '+name)
    return value


def _number(value,name,zero=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (not zero and value==0):
        raise ValueError('Invalid '+name)
    return value


def prepared_source_window(template, header, *, start, stop, chunk_markers,
        issued_chunks, expected_input_identity=None, max_chunks=8, max_records=65536, max_signatures=256):
    """Own only the requested source window; never copy a discarded full index.

    The source template contributes fixed data and independent profile fields.
    Encoded ranges and marker count are rebuilt from the bound header. Returned
    data/profile/geometry are independently mutable, preserving caller isolation.
    Memory-only admission templates require the admission's input identity;
    they cannot silently bind to a replacement input inspected later.
    """
    from .pgen_work_bounds import PgenHeaderWork
    if not isinstance(template,dict) or not isinstance(header,PgenHeaderWork):
        raise ValueError('A source template and bound header are required')
    _integer(issued_chunks,'issued_chunks',0)
    _integer(start,'start',0);_integer(stop,'stop',1)
    data=template['data'];profile=template['profile']
    original=data['encoded']
    bound=original.get('input_identity')
    if expected_input_identity is None:expected_input_identity=bound
    if (expected_input_identity!=header.input_identity or
            (bound is not None and bound!=expected_input_identity) or
            Path(original['path']).resolve(strict=True)!=header.path):
        raise ValueError('Window header differs from the admitted source input')
    if data['samples']!=header._header.sample_ct or original['samples']!=data['samples']:
        raise ValueError('Window sample count differs from its source input')
    span=original['variant_range']
    if not isinstance(span,(tuple,list)) or len(span)!=2 or any(type(v) is not int for v in span):
        raise ValueError('Invalid fixed source partition')
    lo,hi=span
    if not 0<=lo<hi<=header._header.variant_ct:raise ValueError('Invalid fixed source partition')
    if not lo<=start<stop<=hi:raise ValueError('Window exceeds its fixed source partition')
    encoded=header.window(start,stop,chunk_markers,max_chunks=max_chunks,
        max_records=max_records,max_signatures=max_signatures)
    return dict(device=template['device'],trait_range=deepcopy(template['trait_range']),issued_chunks=issued_chunks,
        data=dict(deepcopy({key:value for key,value in data.items() if key not in ('encoded','markers')}),
            markers=stop-start,encoded=encoded),
        profile=dict(deepcopy({key:value for key,value in profile.items() if key!='chunk_markers'}),
            chunk_markers=chunk_markers))


def prepared_window_runtime(windows, *, total_traits, reduction, output, shared_capacities,
        endpoint, host_serial_fraction, prices=None, retained=None, significance_threshold=None, partition_axis=None,
        shared_links=(), shared_storage_bytes_per_second=None, host_serial_policy='fluid',
        max_windows=8, max_source_chunks=16, max_records=65536, max_output_blocks=256, max_writeback_actions=256,
        max_graph_nodes=100000, return_graph=False):
    """Price only these prepared windows, through payload/part drain and fsync.

    Every window has device, trait_range, data, profile and issued_chunks.
    Reduced outputs require one explicit survivor count per window/chunk.
    Dense windows use separate stores; indexed modes share one global writer.
    No phenotype/factor preparation, prior queued work, metadata publication,
    output-state continuation or future-source extrapolation is included.
    """
    for name,value in [('total_traits',total_traits),('max_windows',max_windows),
        ('max_source_chunks',max_source_chunks),('max_records',max_records),
        ('max_output_blocks',max_output_blocks),('max_writeback_actions',max_writeback_actions),('max_graph_nodes',max_graph_nodes)]:_integer(value,name)
    if reduction not in (None,'significant','jagwas'):raise ValueError('Unsupported window output mode')
    axis=('variant' if reduction=='jagwas' else 'trait') if partition_axis is None else partition_axis
    if axis not in ('trait','variant') or (reduction=='jagwas' and axis!='variant') or (reduction=='significant' and axis!='trait'):
        raise ValueError('Partition axis is not executable for this output mode')
    if endpoint not in ('lower','upper'):raise ValueError('Explicit decoder endpoint required')
    if type(return_graph) is not bool:raise ValueError('Boolean return_graph required')
    if host_serial_policy not in ('fluid','held-first','held-last'):raise ValueError('Unknown host serial policy')
    if _number(host_serial_fraction,'host_serial_fraction',True)>1:raise ValueError('Host serial fraction exceeds one')
    if not isinstance(windows,list) or not 0<len(windows)<=max_windows:raise ValueError('Bounded nonempty window list required')
    if set(output)!={'block_bytes','queue_depth','store_beta','fsync'} or type(output['store_beta']) is not bool or output['fsync'] is not True:
        raise ValueError('Explicit durable output contract required')
    _integer(output['queue_depth'],'queue_depth')
    if output['block_bytes'] is not None:_integer(output['block_bytes'],'block_bytes')
    if reduction and output['block_bytes'] is not None:raise ValueError('Indexed output cannot coalesce dense blocks')
    if reduction=='jagwas' and output['store_beta']:raise ValueError('JAGWAS does not store beta')
    if reduction is None:
        if retained is not None or prices is not None:raise ValueError('Dense output has no survivor or indexed price input')
    elif not isinstance(retained,list) or len(retained)!=len(windows) or not isinstance(prices,dict):
        raise ValueError('Explicit per-window survivors and independent indexed prices required')
    if significance_threshold is not None:
        if reduction!='significant' or not 0<_number(significance_threshold,'significance_threshold')<=1:
            raise ValueError('Significance threshold applies only to significant pairs')
    if reduction=='significant':
        from .numpy_nonzero_work import validate_host_price_protocol
        validate_host_price_protocol(prices)
    if set(shared_capacities)!={'cpu','dram','input','output'}:raise ValueError('Explicit shared CPU/DRAM/input/output capacities required')
    caps={key:_number(value,'shared '+key) for key,value in shared_capacities.items()};caps['host_serial']=1.
    if shared_storage_bytes_per_second is not None:caps['storage']=_number(shared_storage_bytes_per_second,'shared storage')
    devices=[];chunks=records=dense_blocks=writeback_actions=0;identity=None;dimensions=None;rectangles=[]
    for index,window in enumerate(windows):
        if set(window)!={'device','trait_range','data','profile','issued_chunks'}:raise ValueError('Explicit window fields required')
        device=canonical_cuda_device(window['device']);data=window['data'];profile=window['profile'];encoded=data['encoded']
        _integer(window['issued_chunks'],'issued_chunks',0)
        if device not in devices:devices.append(device)
        span=window['trait_range']
        if not isinstance(span,(list,tuple)) or len(span)!=2 or any(type(v) is not int for v in span) or not 0<=span[0]<span[1]<=total_traits:
            raise ValueError('Window phenotype range differs')
        if span[1]-span[0]!=data['traits_analyzed']:raise ValueError('Window phenotype width differs')
        if axis=='variant' and list(span)!=[0,total_traits]:raise ValueError('Variant shards retain the complete phenotype panel')
        if profile.get('compute_dtype','float32')!='float32' or profile.get('validate_range') is not True:
            raise ValueError('Public native FP32 range-validated scan required')
        ownership='owned' if reduction=='jagwas' else 'borrowed'
        if profile.get('result_ownership')!=ownership or profile.get('reduction')!=('jagwas' if reduction=='jagwas' else None):
            raise ValueError('Native result contract differs from output mode')
        if reduction=='jagwas' and (list(span)!=[0,total_traits] or data.get('phenotype_complete') is not True):
            raise ValueError('JAGWAS requires the complete phenotype panel')
        if encoded.get('kind')!=WINDOW_KIND or not isinstance(encoded.get('chunks'),list) or not encoded['chunks']:
            raise ValueError('Typed nonempty header window required')
        n=_integer(data['samples'],'samples');m=_integer(data['markers'],'markers')
        dims=(n,data['covariates'],data.get('covariate_columns',data['covariates']))
        if dimensions is not None and dims!=dimensions:raise ValueError('Window sample/covariate dimensions differ')
        dimensions=dims
        if identity is not None and encoded['input_identity']!=identity:raise ValueError('Windows must share an input binding')
        identity=encoded['input_identity']
        variant_span=encoded['variant_range']
        if not isinstance(variant_span,(list,tuple)) or len(variant_span)!=2 or any(type(v) is not int for v in variant_span) or variant_span[0]<0 or variant_span[1]-variant_span[0]!=m:
            raise ValueError('Window variant range differs')
        for previous_traits,previous_variants in rectangles:
            if axis=='trait' and max(span[0],previous_traits[0])<min(span[1],previous_traits[1]):
                raise ValueError('Trait partitions overlap')
            if max(span[0],previous_traits[0])<min(span[1],previous_traits[1]) and max(variant_span[0],previous_variants[0])<min(variant_span[1],previous_variants[1]):
                raise ValueError('Window association rectangles overlap')
        rectangles.append((span,variant_span))
        chunks+=len(encoded['chunks']);records+=m
        if chunks>max_source_chunks or records>max_records:raise ValueError('Source window expansion exceeds budget')
        # The profile fixes independently measured node service. A conditional
        # live-availability scenario may reduce aggregate capacity without
        # rewriting those prices or the profile digest used by JIT binding.
        for resource,field in [('cpu','cpu_available_cores'),('dram','shared_dram_bytes_per_second'),('input','read_bytes_per_second'),('output','write_bytes_per_second')]:
            if caps[resource]>_number(profile[field],field):
                raise ValueError('Shared '+resource+' exceeds window profile capacity')
        for direction in ('h2d','d2h'):
            key=device+':'+direction;value=_number(profile[direction+'_bytes_per_second'],key)
            if key in caps and caps[key]!=value:raise ValueError('Conflicting per-device transfer capacity')
            caps[key]=value
        if reduction is None:
            first=encoded['chunks'][0]['markers']
            memory=binary_output_memory(m,data['traits_analyzed'],first,output['block_bytes'],output['queue_depth'],output['store_beta'],True)
            for name,size in memory['block_bytes_by_stream'].items():
                amount=4*m*(1 if name=='df' else data['traits_analyzed'])
                dense_blocks+=(amount+size-1)//size
                submitted=amount//(64<<20)
                writeback_actions+=submitted+2*max(0,submitted-1)
            if dense_blocks>max_output_blocks:raise ValueError('Dense window output expansion exceeds budget')
            if writeback_actions>max_writeback_actions:raise ValueError('Dense window writeback expansion exceeds budget')
        else:
            counts=retained[index]
            if not isinstance(counts,list) or len(counts)!=len(encoded['chunks']):raise ValueError('One survivor count per source chunk required')
            for count,row in zip(counts,encoded['chunks']):
                _integer(count,'survivor count',0)
                if count>row['markers']*(data['traits_analyzed'] if reduction=='significant' else 1):raise ValueError('Survivors exceed source cells')
    if axis=='variant' and len(devices)!=len(windows):raise ValueError('Variant partitioning requires one window per active-device shard')
    for index,link in enumerate(shared_links):
        if set(link)!={'devices','h2d_bytes_per_second','d2h_bytes_per_second'} or not link['devices'] or len(set(link['devices']))!=len(link['devices']) or not set(link['devices'])<=set(devices):
            raise ValueError('Unique known devices required for a shared link')
        for direction in ('h2d','d2h'):caps[f'link:{index}:{direction}']=_number(link[direction+'_bytes_per_second'],'shared link')
    graph=ExecutionGraph();graph.capacities=dict(caps);done={};indexed=[];reports=[];unpriced=set()
    with planning_work_scope():
        for index,window in enumerate(windows):
            data,profile,device=window['data'],window['profile'],window['device'];q=profile['cpu_fraction']
            work=torch_scan_header_work(data,profile,endpoint=endpoint,issued_chunks=window['issued_chunks'],max_chunks=max_source_chunks,max_records=max_records)
            blocks=deepcopy(work['blocks']);unpriced.update(work['unpriced_terms']+work['allocator_unpriced_terms'])
            for block in blocks:
                block['host_resources']['host_serial']=block['host_resources']['cpu']*host_serial_fraction
                for direction in ('h2d','d2h'):
                    rate=block[direction+'_bytes']/block[direction+'_seconds'] if block[direction+'_seconds'] else 0.
                    resources=block.setdefault(direction+'_resources',{});resources[device+':'+direction]=rate
                    for li,link in enumerate(shared_links):
                        if device in link['devices']:resources[f'link:{li}:{direction}']=rate
            report=dict(device=device,trait_range=list(window['trait_range']),variant_range=list(data['encoded']['variant_range']),
                issued_chunks=window['issued_chunks'],decode_workers=work['workers'],depth=work['depth'],
                source_chunks=len(blocks),d2h_bytes=sum(b['d2h_bytes'] for b in blocks))
            if reduction is None:
                writer_work=_output_work(window,output);copy=writer_copy_cost(profile)
                writer=BinaryWriterSchedule(writer_work,copy_seconds_per_byte=copy['cpu_seconds_per_byte']/q,
                    copy_seconds_per_call=copy['cpu_seconds_per_call']/q,zero_seconds_per_byte=profile['process_units']['bytearray_zero_bytes']/q,
                    write_seconds_per_byte=1/profile['write_bytes_per_second'],fsync_seconds_per_array=profile['fsync_seconds'],
                    append_seconds=profile['executor_cpu_seconds']/q,handoff_seconds=profile['executor_cpu_seconds']/q,
                    cpu_fraction=q,write_capacity=profile['write_bytes_per_second'],
                    writeback_service={key:value if key=='storage_seconds_per_byte' else value/q for key,value in profile['writeback_service'].items()})
                local=torch_scan_schedule(blocks,depth=work['depth'],decode_workers=work['workers'],consumer=writer,
                    shared_capacities=caps,return_graph=True,consumer_open_before_scan=output['block_bytes'] is not None)
                for name,demands in local.demands.items():
                    if name.startswith('writer:') and demands.get('cpu'):demands['host_serial']=demands['cpu']*host_serial_fraction
                done[device]=graph.compose(local,f'window:{index}:',[done[device]] if device in done else [])
                report.update(payload_bytes=writer.payload,write_calls=writer.write_calls,fsync_calls=writer_work['fsync_calls'],
                    writer_streams={name:dict(block_bytes=stream['block_bytes'],
                        payload_bytes=stream['binary_payload_bytes']//stream['arrays'],
                        writeback_interval_bytes=stream['writeback_per_array']['interval_bytes'],
                        writeback_submit_calls=stream['writeback_per_array']['submit_calls'])
                        for name,stream in writer_work['stream_work'].items()})
            else:
                outputs=[];payload=parts=0
                for block,count in zip(blocks,retained[index]):
                    if reduction=='jagwas':
                        from .reduced_output_work import jagwas_writer_work,jagwas_host_selection_service,jagwas_archive_service
                        output_work=jagwas_writer_work(block['markers'],count,fsync=True)
                        selection=jagwas_host_selection_service(output_work,prices['prices'],cpu_fraction=q,dram_bytes_per_second=profile['shared_dram_bytes_per_second'],host_serial_fraction=host_serial_fraction)
                        writer=jagwas_archive_service(output_work,prices['archive'],profile,host_serial_fraction=host_serial_fraction)
                        part=output_work['part'];cells=block['markers']
                        release=block.pop('discard_seconds',0.)
                        if release:(writer if count else selection).append(dict(seconds=release,resources=dict(cpu=q,host_serial=q)))
                    else:
                        from .significant_host_work import host_significant_selection_work,host_selection_service,indexed_part_work
                        from .significant_host_model import _writer_service
                        cells=block['markers']*data['traits_analyzed']
                        selected=host_significant_selection_work(block['markers'],data['traits_analyzed'],count,threshold_one=significance_threshold==1.,return_beta=profile.get('return_beta',True))
                        unpriced.update(selected['allocation_unpriced_terms'])
                        selection=host_selection_service(selected,prices['prices'],cpu_fraction=q,dram_bytes_per_second=profile['shared_dram_bytes_per_second'],host_serial_fraction=host_serial_fraction)
                        part=indexed_part_work(count,store_beta=output['store_beta']);writer=_writer_service(part,prices['archive'],profile,host_serial_fraction)
                    outputs.append([dict(cells=cells,retained=count,selection=selection,writer=writer)])
                    payload+=part['file_bytes'];parts+=bool(count)
                indexed.append(dict(device=device,backend='host',blocks=blocks,outputs=outputs,depth=work['depth'],decode_workers=work['workers'],prepare=ExecutionGraph(),cleanup=[]))
                report.update(payload_bytes=payload,parts=parts,retained=sum(retained[index]),fsync_calls=parts)
            reports.append(report)
    if reduction is not None:
        from .indexed_schedule import jagwas_variant_schedule,significant_trait_schedule
        queue=None
        if len(devices)>1:
            q=min(w['profile']['cpu_fraction'] for w in windows)
            queue={kind:[dict(seconds=_number(prices['queue_cpu_seconds'][kind],'queue '+kind,True)/q,resources=dict(cpu=q,host_serial=q))] for kind in ('put','get')}
        options=dict(queue_depth=output['queue_depth'] if queue else 0,shared_capacities=caps,queue_service=queue,
            finalize=[],return_graph=True,max_source_chunks=max_source_chunks,max_selection_blocks=max_source_chunks)
        graph=(jagwas_variant_schedule(indexed,shared_prepare=ExecutionGraph(),**options) if reduction=='jagwas' else significant_trait_schedule(indexed,**options))
    for demands in graph.demands.values():
        if shared_storage_bytes_per_second is not None:demands['storage']=demands.get('input',0.)+demands.get('output',0.)
    if host_serial_policy!='fluid':graph=graph.with_serial_sections(host_serial_policy)
    if len(graph.nodes)>max_graph_nodes:raise ValueError('Window graph exceeds node budget')
    if return_graph:return graph
    solution=graph.solve()
    return dict(status='development_prepared_window',estimated_window_seconds=solution['seconds'],
        resource_balance=resource_balance(graph,solution),windows=reports,source_chunks=chunks,graph_nodes=len(graph.nodes),
        payload_bytes=sum(r['payload_bytes'] for r in reports),decoder_endpoint=endpoint,reduction=reduction,partition_axis=axis,
        input_identity=deepcopy(identity),unpriced_terms=sorted(unpriced),prediction_complete=False,selection_validated=False,
        scope='Isolated prepared source windows with empty initial queues/writer staging and final payload/part drain. No phenotype/factor setup, prior in-flight work, API metadata/directory publication, live checkpoint or forecast for unobserved variants/tiles. Decoder endpoint simulations are scenarios, not makespan bounds. Independent-price validity and live memory admission remain caller responsibilities.')


def _coverage(windows):
    """Canonical bounded rectangle union, without allocating any pair matrix."""
    for window in windows:
        for span in (window['trait_range'],window['data']['encoded']['variant_range']):
            if not isinstance(span,(list,tuple)) or len(span)!=2 or any(type(v) is not int for v in span) or not 0<=span[0]<span[1]:
                raise ValueError('Invalid comparison rectangle')
    endpoints=sorted({v for w in windows for v in w['data']['encoded']['variant_range']})
    bands=[]
    for lo,hi in zip(endpoints,endpoints[1:]):
        spans=sorted(tuple(w['trait_range']) for w in windows
            if w['data']['encoded']['variant_range'][0]<=lo and hi<=w['data']['encoded']['variant_range'][1])
        merged=[]
        for start,stop in spans:
            if merged and start<merged[-1][1]:raise ValueError('Association coverage overlaps')
            if merged and start==merged[-1][1]:merged[-1][1]=stop
            else:merged.append([start,stop])
        if not merged:continue
        if bands and bands[-1]['variant_range'][1]==lo and bands[-1]['trait_ranges']==merged:
            bands[-1]['variant_range'][1]=hi
        else:bands.append(dict(variant_range=[lo,hi],trait_ranges=merged))
    return bands


def _rebin_survivors(windows, evidence, coverage, reduction, total_traits, threshold, maximum):
    """Only split a count when it proves every cell or no cell survives."""
    if not isinstance(evidence,dict) or set(evidence)!={'input_identity','reduction','total_traits','significance_threshold','bins'}:
        raise ValueError('Explicit bound survivor evidence required')
    if (evidence['input_identity']!=windows[0]['data']['encoded']['input_identity'] or evidence['reduction']!=reduction
        or evidence['total_traits']!=total_traits or evidence['significance_threshold']!=threshold):
        raise ValueError('Survivor evidence binding differs')
    bins=evidence['bins']
    if not isinstance(bins,list) or not 0<len(bins)<=maximum:raise ValueError('Bounded survivor bins required')
    for row in bins:
        if set(row)!={'variant_range','trait_range','retained'}:raise ValueError('Explicit survivor rectangle and count required')
        for field in ('variant_range','trait_range'):
            span=row[field]
            if not isinstance(span,(list,tuple)) or len(span)!=2 or any(type(v) is not int for v in span) or not 0<=span[0]<span[1]:raise ValueError('Invalid survivor range')
        if row['trait_range'][1]>total_traits:raise ValueError('Survivor phenotype range exceeds the panel')
        if reduction=='jagwas' and list(row['trait_range'])!=[0,total_traits]:raise ValueError('JAGWAS survivor bins require the complete panel')
        cells=(row['variant_range'][1]-row['variant_range'][0])*(row['trait_range'][1]-row['trait_range'][0] if reduction=='significant' else 1)
        if _integer(row['retained'],'survivor count',0)>cells:raise ValueError('Survivors exceed evidence cells')
    if _coverage([dict(trait_range=r['trait_range'],data=dict(encoded=dict(variant_range=r['variant_range']))) for r in bins])!=coverage:
        raise ValueError('Survivor evidence must cover the exact comparison workload')
    result=[]
    for window in windows:
        traits=window['trait_range'];counts=[]
        for lo,hi in window['data']['encoded']['chunk_ranges']:
            count=0
            for row in bins:
                a,b=row['variant_range'];c,d=row['trait_range']
                variants=max(0,min(hi,b)-max(lo,a));width=max(0,min(traits[1],d)-max(traits[0],c))
                if not variants or not width:continue
                cells=(b-a)*(d-c if reduction=='significant' else 1)
                if lo<=a and b<=hi and traits[0]<=c and d<=traits[1]:count+=row['retained']
                elif row['retained']==0:pass
                elif row['retained']==cells:count+=variants*(width if reduction=='significant' else 1)
                else:raise ValueError('A partial survivor bin cannot be split; finer evidence is required')
            counts.append(count)
        result.append(counts)
    return result


def compare_prepared_windows(baseline, candidate, *, survivor_evidence=None, max_survivor_bins=128,
        model_identity=None, **common):
    """Compare two bounded layouts of the same associations/output contract.

    Each layout supplies windows and partition_axis. Reduction,
    threshold, writer settings and independent prices are common immutable
    inputs. Reduced layouts are rebinned from one common count ledger; partial
    bins cannot be split without finer evidence. Empty/full bins can be split.
    This is an isolated finite-window comparison, not a remaining-job proposal.
    """
    started=time.perf_counter();cpu=time.thread_time()
    if 'return_graph' in common or 'partition_axis' in common or 'retained' in common:
        raise ValueError('Comparison owns graph return, partition and survivor options')
    maximum=_integer(common.get('max_windows',8),'max_windows')
    chunk_limit=_integer(common.get('max_source_chunks',16),'max_source_chunks')
    _integer(max_survivor_bins,'max_survivor_bins')
    for layout in (baseline,candidate):
        if set(layout)!={'windows','partition_axis'} or not isinstance(layout['windows'],list) or not 0<len(layout['windows'])<=maximum:
            raise ValueError('Explicit bounded window layout required')
        if sum(len(w['data']['encoded']['chunk_ranges']) for w in layout['windows'])>chunk_limit:
            raise ValueError('Comparison source chunks exceed budget')
    coverage=_coverage(baseline['windows'])
    if coverage!=_coverage(candidate['windows']):raise ValueError('Layouts must cover identical variant-phenotype pairs')
    def binding(layout):
        w=layout['windows'][0];d=w['data']
        return (d['encoded']['input_identity'],d['samples'],d['covariates'],d.get('covariate_columns',d['covariates']))
    if binding(baseline)!=binding(candidate):raise ValueError('Comparison input/sample/covariate bindings differ')
    if model_identity is not None and (not isinstance(model_identity,dict) or not model_identity):
        raise ValueError('Explicit nonempty caller-validated model identity required')
    def configuration(layout):
        return dict(partition_axis=layout['partition_axis'],windows=[dict(device=w['device'],
            trait_range=list(w['trait_range']),variant_start=w['data']['encoded']['variant_range'][0],
            issued_chunks=w['issued_chunks'],chunk_markers=w['profile']['chunk_markers'],profile_sha256=_digest(w['profile']),
            chunk_invariant_profile_sha256=_digest({key:value for key,value in w['profile'].items()
                if key not in ('chunk_markers','kernel_geometry','joint_kernel_geometry')}),
            fixed_data_sha256=_digest({key:value for key,value in w['data'].items() if key not in ('markers','encoded')}))
            for w in layout['windows']])
    # A later horizon must price the same configuration and output contract.
    # Source/runtime validity is supplied by an already validated caller binding;
    # reading or hashing a package here would add hidden per-comparison I/O.
    contract=dict(schema='torchgwas.prepared_comparison.v1',model_identity=deepcopy(model_identity),
        input_binding=deepcopy(binding(baseline)),total_traits=common['total_traits'],reduction=common['reduction'],
        output=deepcopy(common['output']),significance_threshold=common.get('significance_threshold'),
        shared_capacities=deepcopy(common['shared_capacities']),shared_links=deepcopy(common.get('shared_links',())),
        shared_storage_bytes_per_second=common.get('shared_storage_bytes_per_second'),endpoint=common['endpoint'],
        host_serial_fraction=common['host_serial_fraction'],host_serial_policy=common.get('host_serial_policy','fluid'),
        indexed_prices_sha256=None if common.get('prices') is None else _digest(common['prices']),
        baseline=configuration(baseline),candidate=configuration(candidate))
    # Persisted and freshly constructed comparisons must have the same binding.
    # In particular, JSON represents both tuples and lists as arrays.
    contract=json.loads(_json(contract))
    counts=[None,None]
    if common.get('reduction') is not None:
        counts=[_rebin_survivors(layout['windows'],survivor_evidence,coverage,common['reduction'],common['total_traits'],
            common.get('significance_threshold'),max_survivor_bins) for layout in (baseline,candidate)]
    elif survivor_evidence is not None:raise ValueError('Dense comparisons have no survivor evidence')
    with planning_work_scope():
        before=prepared_window_runtime(baseline['windows'],partition_axis=baseline['partition_axis'],retained=counts[0],**common)
        after=prepared_window_runtime(candidate['windows'],partition_axis=candidate['partition_axis'],retained=counts[1],**common)
    gain=before['estimated_window_seconds']-after['estimated_window_seconds']
    return dict(baseline=before,candidate=after,coverage=coverage,comparison_contract=contract,
        estimated_window_gain_seconds=gain,
        preferred_window='candidate' if gain>0 else 'baseline',calculation_wall_seconds=time.perf_counter()-started,
        calculation_cpu_seconds=time.thread_time()-cpu,selection_validated=False,
        scope='Equivalent finite prepared-window model comparison. Calculation cost excludes caller header construction and parameter binding. No extrapolation, live-switch approval, actual speedup guarantee or claim of globally optimal layout.')
