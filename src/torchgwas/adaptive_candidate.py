"""Shared allocation admission for future chunk changes on fixed partitions.

This uses the existing source tensor, PGEN, pinned and output ledgers. It is
neither a timing fit nor a candidate-switching policy. Static chunk candidates
do not by themselves admit mixed sizes in one retained allocation.
"""
from __future__ import annotations

from collections import deque
import math

from .adaptive_chunks import aligned_chunk_sizes,aligned_chunk_shapes
from .decoder_work import native_read_layout,native_reader_workspace
from .mechanistic_plan import _integer
from .pinned_work import pinned_scan_work
from .tensor_memory import eager_scan_memory,eager_memory_plan
from .trait_tiling_model import trait_tiled_shape,trait_tiled_memory


def adaptive_tensor_memory(samples,traits,covariates,depth,*,chunk_sizes,markers,
                           device_profile,reduction=None,return_beta=True):
    """Component-wise maxima retain the fixed ring and mixed old/new outputs."""
    sizes=aligned_chunk_sizes(chunk_sizes);capacity=sizes[-1]
    shapes=aligned_chunk_shapes(sizes,markers)
    plans=[eager_memory_plan(samples,b,traits,covariates,depth,device_profile,reduction=reduction,return_beta=return_beta)
           for b in shapes]
    scans=[plan['scan'] for plan in plans]
    # Old conversion/results can belong to a different size than the currently
    # evaluated RHS. Maximize each storage class, not only whole-plan peaks.
    classes=['resident_bytes','distinct_temporary_bytes','retained_output_bytes',
             'previous_genotype_bytes']
    if reduction=='jagwas':classes.append('previous_reduced_result_bytes')
    if reduction=='device_significant':classes.append('previous_dense_result_bytes')
    if not return_beta:classes.append('previous_beta_bytes')
    allocations={name:max(scan[name] for scan in scans) for name in classes}
    allocations['staging_bytes']=depth*samples*capacity
    tensor_bytes=sum(allocations.values())
    library=max(plan['cublas']['total_bytes'] for plan in plans)
    setup=max(plan['setup_bytes'] for plan in plans)
    return dict(capacity=capacity,chunk_sizes=list(sizes),work_shapes=list(shapes),
        device_bytes=max(setup,tensor_bytes+library),tensor_allocations=allocations,
        tensor_storage_budget=tensor_bytes,library_bytes=library,setup_bytes=setup,
        pinned=pinned_scan_work(samples,capacity,traits,depth,reduction=reduction,return_beta=return_beta),
        unresolved_memory_terms=sorted({term for plan in plans for term in plan['unresolved_memory_terms']}),
        prediction_complete=False,
        scope='Fixed native ring plus component-wise maxima of source storage classes across admitted chunk/tail shapes. Explicit reserves are still required for unresolved allocations.')


def _fine_chunks(encoded,sizes,*,max_census_chunks):
    """Validate one complete small-grid census without materializing copies."""
    fine=sizes[0];parts=encoded.get('chunks')
    first,last=encoded['variant_range']
    if (first!=0 or last!=encoded['file_markers'] or encoded['markers']!=last
            or encoded['chunk_markers']!=fine):
        raise ValueError('A complete file census on the smallest chunk grid is required')
    if (not isinstance(parts,list) or len(parts)!=math.ceil(last/fine)
            or len(parts)>max_census_chunks):
        raise ValueError('Complete bounded fine-grid census required')
    context=['path','samples','file_markers','file_bytes','index_bytes']
    for index,part in enumerate(parts):
        lo=index*fine;hi=min(lo+fine,last)
        if (part.get('variant_range')!=[lo,hi] or part.get('markers')!=hi-lo
                or part.get('chunk_markers')!=fine
                or any(part.get(key)!=encoded.get(key) for key in context)):
            raise ValueError('Fine-grid source ranges or file context differ')
        if type(part['record_payload_bytes']) is not int or part['record_payload_bytes']<0:
            raise ValueError('Nonnegative exact source payload required')
    return parts


def _reader_envelope(parts,*,capacity,fine):
    # Each possible read starts on the fine grid. Its contiguous file payload
    # is at most the next capacity/fine records, plus the first part's LD
    # prefix. Interior fine-chunk restart prefixes must not be summed.
    window=deque();payload=0;read=0;scratch=0
    width=capacity//fine
    for part in reversed(parts):
        window.appendleft(part['record_payload_bytes']);payload+=window[0]
        if len(window)>width:payload-=window.pop()
        layout=native_read_layout(part)
        prefix=layout['read_bytes']-part['record_payload_bytes']
        read=max(read,payload+prefix)
        scratch=max(scratch,layout['extra_workspace_bytes'])
    return read,scratch


def future_chunk_candidate(candidate,*,source_census,chunk_sizes,next_size,issued_chunks,
                           reduction,issued_ranges=None,max_source_chunks=10000,max_census_chunks=1000000):
    """One bounded counterfactual, preserving every already-issued source chunk.

    issued_chunks is one count per tile/shard in candidate order, supplied by
    the executor or an analytical checkpoint adapter. It must include issued
    reads even when their results have not arrived. The full prefix is retained
    so the existing graph can resume with its queues and active services.
    This function neither reads input nor searches configurations. Admission of
    the shared allocation envelope and live empirical validity are separate.
    When issued_ranges is supplied it is the exact executor prefix per tile,
    which can differ from the old regular candidate. Such a new graph must not
    be resumed from an incompatible older analytical checkpoint.
    """
    from copy import deepcopy
    from .adaptive_chunks import AlignedChunkSizeControl
    from .decoder_work import census_chunk_ranges
    from .pgen_work_census import scheduled_census
    if source_census.get('kind') in ('torchgwas.pgen_memory_layout.v1',
            'torchgwas.pgen_compact_memory_envelope.v1'):
        raise ValueError('Memory-only layout cannot price a future candidate; exact decoder work evidence is required')
    if reduction not in (None,'significant','jagwas'):
        raise ValueError('Unsupported adaptive reduction')
    sizes=aligned_chunk_sizes(chunk_sizes)
    control=AlignedChunkSizeControl(sizes,initial=next_size)
    _integer('max_source_chunks',max_source_chunks);_integer('max_census_chunks',max_census_chunks)
    shape=trait_tiled_shape(candidate,reduction='jagwas' if reduction=='jagwas' else None)
    if shape['chunk_size']!=sizes[-1]:raise ValueError('Fixed capacity differs from admitted sizes')
    if (not isinstance(issued_chunks,(list,tuple)) or len(issued_chunks)!=len(candidate['tiles'])
        or any(type(value) is not int or value<0 for value in issued_chunks)):
        raise ValueError('One nonnegative issued-chunk count per tile required')
    if sum(issued_chunks)>max_source_chunks:
        raise ValueError('Issued prefix exceeds source-chunk budget')
    if issued_ranges is not None and (not isinstance(issued_ranges,(list,tuple))
            or len(issued_ranges)!=len(candidate['tiles'])):
        raise ValueError('One exact issued-range prefix per tile required')
    if (source_census.get('path')!=shape['input_path'] or source_census.get('samples')!=shape['samples']
        or source_census.get('markers')!=shape['markers']):
        raise ValueError('Source census differs from candidate')
    _fine_chunks(source_census,sizes,max_census_chunks=max_census_chunks)
    future=deepcopy({key:value for key,value in candidate.items() if key!='tiles'})
    future['tiles']=[];expanded=0
    for index,(tile,count) in enumerate(zip(candidate['tiles'],issued_chunks)):
        data,profile=tile['data'],tile['profile'];encoded=data['encoded']
        if any(encoded.get(key)!=source_census.get(key) for key in ('path','samples','file_markers','file_bytes','index_bytes')):
            raise ValueError('Candidate and fine-grid source identities differ')
        old=census_chunk_ranges(encoded,profile['chunk_markers'])
        if issued_ranges is None:
            if count>len(old):raise ValueError('Issued chunks exceed the source schedule')
            ranges=old[:count]
        else:
            prefix=issued_ranges[index]
            if not isinstance(prefix,(list,tuple)) or len(prefix)!=count:
                raise ValueError('Issued ranges differ from the supplied count')
            cursor=encoded['variant_range'][0];ranges=[]
            for span in prefix:
                if (not isinstance(span,(list,tuple)) or len(span)!=2 or any(type(x) is not int for x in span)
                        or span[0]!=cursor or not cursor<span[1]<=encoded['variant_range'][1]):
                    raise ValueError('Exact contiguous issued source ranges required')
                ranges.append(list(span));cursor=span[1]
        expanded+=count
        if expanded>max_source_chunks:raise ValueError('Future candidate exceeds source-chunk budget')
        lo,hi=encoded['variant_range'];cursor=ranges[-1][1] if ranges else lo
        shapes=set(aligned_chunk_shapes(sizes,data['markers']))
        if any(end-start not in shapes for start,end in ranges):
            raise ValueError('Issued shape is outside the admitted chunk set')
        while cursor<hi:
            if expanded>=max_source_chunks:raise ValueError('Future candidate exceeds source-chunk budget')
            end=cursor+control(cursor,hi,sizes[-1]);ranges.append([cursor,end]);cursor=end;expanded+=1
        updated=scheduled_census(source_census,sizes[-1],ranges,max_source_chunks=max_census_chunks)
        if issued_ranges is None and updated['chunks'][:count]!=encoded['chunks'][:count]:
            raise ValueError('Issued source counts or LD replay changed')
        copied=deepcopy({key:value for key,value in tile.items() if key!='data'})
        copied['data']=deepcopy({key:value for key,value in data.items() if key!='encoded'})
        copied['data']['encoded']=updated;future['tiles'].append(copied)
    return future


def adaptive_candidate_memory(candidate,*,chunk_sizes,source_census,reduction,
                              device_memory_profiles,host_reserve_bytes=0,
                              device_reserve_bytes=0,max_census_chunks=1000000,
                              max_kernel_shapes=128):
    """Admit an unchanged device/trait partition with one fixed largest ring.

    Use AlignedChunkSizeControl with the returned sizes/capacity. Devices,
    phenotype tiles, prefetch depth, reader budgets and variant intervals remain
    fixed. A later policy may change only the next not-yet-issued chunk size.
    The fine-grid census covers every possible native LD restart and source
    payload window, including starts absent from any regular coarse plan.
    """
    if reduction not in (None,'significant','jagwas'):
        raise ValueError('Adaptive candidate supports dense, host-significant and JAGWAS output')
    sizes=aligned_chunk_sizes(chunk_sizes);capacity=sizes[-1]
    for name,value in [('max_census_chunks',max_census_chunks),('max_kernel_shapes',max_kernel_shapes)]:
        _integer(name,value)
    shape=trait_tiled_shape(candidate,reduction='jagwas' if reduction=='jagwas' else None)
    if shape['chunk_size']!=capacity:
        raise ValueError('The fixed candidate must allocate the largest admitted chunk')
    if set(device_memory_profiles)!=set(shape['devices']):
        raise ValueError('Exactly one device memory profile per active adaptive device required')
    if (source_census.get('path')!=shape['input_path']
            or source_census.get('samples')!=shape['samples']
            or source_census.get('markers')!=shape['markers']):
        raise ValueError('Fine-grid source census differs from candidate')
    required={}
    for tile in candidate['tiles']:
        data=tile['data'];width=data['traits_analyzed']
        for b in aligned_chunk_shapes(sizes,data['markers']):
            key=(tile['device'],data['samples'],b,width,data['covariates'])
            required[key]=dict(device=tile['device'],shape=list(key[1:]),reduction=reduction)
    if len(required)>max_kernel_shapes:
        raise ValueError('Adaptive candidate exceeds max_kernel_shapes')
    from .pgen_memory_layout import COMPACT_KIND,compact_shifted_reader_envelope
    compact=source_census.get('kind')==COMPACT_KIND
    if compact:
        if (source_census.get('chunk_markers')!=sizes[0]
                or source_census.get('logical_chunks')!=math.ceil(shape['markers']/sizes[0])
                or source_census['logical_chunks']>max_census_chunks):
            raise ValueError('Complete bounded compact fine-grid source required')
        parts=None
    else:parts=_fine_chunks(source_census,sizes,max_census_chunks=max_census_chunks)
    options=dict(host_reserve_bytes=host_reserve_bytes,device_reserve_bytes=device_reserve_bytes,
                 device_memory_profiles=device_memory_profiles)
    if reduction=='jagwas':
        from .jagwas_candidate import jagwas_candidate_memory
        base=jagwas_candidate_memory(candidate,**options)
    elif reduction=='significant':
        from .significant_host_work import significant_host_memory
        base=significant_host_memory(candidate,**options)
    else:base=trait_tiled_memory(candidate,**options)
    rows=[];host_extra=dict.fromkeys(shape['devices'],0);gpu=dict(base['device_bytes'])
    unknown=set(base['unresolved_memory_terms']);missing=[];seen=set()
    tensor_cache={}
    for tile in candidate['tiles']:
        device,data,profile=tile['device'],tile['data'],tile['profile']
        n,m,k,c=[data[key] for key in ('samples','markers','traits_analyzed','covariates')]
        encoded=data['encoded']
        if any(encoded.get(key)!=source_census.get(key) for key in ('path','samples','file_markers','file_bytes','index_bytes')):
            raise ValueError('Candidate and fine-grid encoded file identities differ')
        lo,hi=tile.get('variant_range',[0,m])
        fine=sizes[0]
        if lo%fine or (hi%fine and hi!=source_census['markers']):
            raise ValueError('Fixed variant intervals must align to the smallest chunk grid')
        if compact:
            payload,scratch=compact_shifted_reader_envelope(source_census,(lo,hi),capacity)
        else:
            shard=parts[lo//fine:math.ceil(hi/fine)]
            payload,scratch=_reader_envelope(shard,capacity=capacity,fine=fine)
        old_payload,old_scratch=native_reader_workspace(encoded['chunks'])
        depth=profile['depth'];workers=profile['decode_workers'];packed=(n+3)//4
        old_active=min(depth,workers,math.ceil(m/capacity))
        active=min(depth,workers,math.ceil(m/fine))
        old_decoder=old_active*(min(capacity,m)*packed+old_payload+old_scratch)
        decoder=active*(min(capacity,m)*packed+payload+scratch)
        extra=max(0,decoder-old_decoder)
        host_extra[device]=max(host_extra[device],extra)
        tensor_key=(device,n,m,k,c,depth,profile.get('return_beta',True))
        if tensor_key not in tensor_cache:
            tensor_cache[tensor_key]=adaptive_tensor_memory(n,k,c,depth,chunk_sizes=sizes,markers=m,
                device_profile=device_memory_profiles[device],reduction='jagwas' if reduction=='jagwas' else None,
                return_beta=profile.get('return_beta',True))
        memory=tensor_cache[tensor_key]
        gpu[device]=max(gpu[device],memory['device_bytes']+device_reserve_bytes)
        unknown.update(memory['unresolved_memory_terms'])
        for b in memory['work_shapes']:
            key=(device,n,b,k,c)
            if key in seen:continue
            seen.add(key)
            statistics=[r for r in profile.get('kernel_geometry',[])
                if (r.get('N'),r.get('B'),r.get('K'),r.get('C'),r.get('validate_range'))==(n,b,k,c,True)]
            joint=[r for r in profile.get('joint_kernel_geometry',[])
                if (r.get('N'),r.get('B'),r.get('K'),r.get('compute_dtype'))==(n,b,k,'float32')]
            if len(statistics)!=1 or (reduction=='jagwas' and len(joint)!=1):
                missing.append(required[key])
        rows.append(dict(device=device,trait_range=list(tile['trait_range']),variant_range=[lo,hi],
            work_shapes=memory['work_shapes'],reader_active_upper=active,decoder_workspace_bytes=decoder,
            reader_payload_upper_bytes=payload,reader_replay_scratch_upper_bytes=scratch,
            previous_static_decoder_bytes=old_decoder,tensor_memory=memory))
    return dict(host_bytes=base['host_bytes']+sum(host_extra.values()),device_bytes=gpu,
        base_fixed_capacity_memory=base,decoder_extra_bytes_by_device=host_extra,tiles=rows,
        capacity=capacity,chunk_sizes=list(sizes),required_geometry=list(required.values()),
        missing_geometry=missing,unresolved_memory_terms=sorted(unknown),
        selection_validated=False,prediction_complete=False,
        scope='Shared allocation envelope for bounded future chunk changes on fixed partitions. Conservative source arrays and explicit reserves; no automatic ranking, transition policy or complete allocator bound.')
