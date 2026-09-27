"""Census encoded work without timing an association kernel.

Scope: variable-width, biallelic, unphased hardcall PGEN. LD-base obligations
are recorded at chunk starts, not assumed to mean the previous variant.
"""
import mmap
from collections import Counter
from pathlib import Path
import numpy as np
from .pgen_reader import read_header, PgenFormatError, _bytes_per_sample_id, _uleb128
from .first_principles import pgen_structure


def _counts():
    return dict(forms=Counter(), payload=0, entries=0, varints=0, groups=0,
                delta_integers=0, by_form={}, headers=Counter(), deltas=Counter(),
                tail_high=0, ld_starts=0, replay_bytes=0, native_replays=[])


def _accumulate_record(mm,h,i,id_bytes,targets):
    form=int(h.vrtypes[i]);offset=int(h.record_offsets[i]);length=int(h.record_lengths[i])
    for state in targets:
        state['forms'][form]+=1;state['payload']+=length
    if form==0:
        if length!=(h.sample_ct+3)//4:raise PgenFormatError('Unexpected plain record length')
        return
    prefix=1+(h.sample_ct+7)//8 if form==1 else 0
    if form==1 and h.sample_ct%8:
        # Native decoder's final conditional scalar set_category calls.
        if prefix>=length:raise PgenFormatError('Truncated one-bit record')
        tail_byte=mm[offset+prefix-1]&((1<<(h.sample_ct%8))-1)
        for state in targets:state['tail_high']+=int(tail_byte).bit_count()
    if prefix>=length:raise PgenFormatError('Missing difflist header')
    view=memoryview(mm)[offset:offset+length]
    tail=None
    try:
        entry_count,pos=_uleb128(view,prefix)
        if entry_count>h.sample_ct:raise PgenFormatError('Difflist exceeds cohort')
        count_header=pos-prefix;groups=(entry_count+63)//64
        payload_start=pos+groups*id_bytes+max(groups-1,0)+(entry_count+3)//4
        if payload_start>length:raise PgenFormatError('Truncated difflist')
        tail=np.frombuffer(view[payload_start:],dtype=np.uint8)
        expected_deltas=entry_count-groups
        endings=np.flatnonzero(tail<128)
        if len(endings)!=expected_deltas or (len(tail) and tail[-1]>=128):
            raise PgenFormatError('Difflist variable-integer count mismatch')
        lengths=np.diff(endings,prepend=-1)
        if len(lengths) and int(lengths.max())>5:
            raise PgenFormatError('Sample delta varint exceeds 32-bit encoding')
        # These validated lengths are in 1..5. A fixed-size count avoids a
        # per-record sort/unique and the duplicate predicate traversal.
        histogram={k:int(v) for k,v in enumerate(np.bincount(lengths,minlength=6)) if v}
        variable_bytes=count_header+length-payload_start
    finally:
        del tail
        view.release()
    for state in targets:
        state['headers'][count_header]+=1;state['deltas'].update(histogram)
        state['entries']+=entry_count;state['varints']+=variable_bytes
        state['groups']+=groups;state['delta_integers']+=expected_deltas
        row=state['by_form'].setdefault(form,dict(records=0,entries=0,varint_bytes=0,groups=0))
        row['records']+=1;row['entries']+=entry_count;row['varint_bytes']+=variable_bytes;row['groups']+=groups


def census(path, chunk_markers=2048, variant_range=None, *, include_chunks=False):
    if isinstance(chunk_markers, bool) or not isinstance(chunk_markers, int) or chunk_markers < 1:
        raise ValueError('chunk_markers must be positive integer')
    if type(include_chunks) is not bool:
        raise ValueError('include_chunks must be boolean')
    path=Path(path);h=read_header(path)
    if any(int(t) not in (0,1,2,3,4,6,7) for t in np.unique(h.vrtypes)):
        raise PgenFormatError('Census refuses reserved types, dosage, phase and multiallelic payloads')
    start,stop=(0,h.variant_ct) if variant_range is None else variant_range
    if any(isinstance(v,bool) or not isinstance(v,int) for v in (start,stop)) or not 0<=start<stop<=h.variant_ct:
        raise ValueError('variant_range must be a nonempty in-file half-open interval')
    id_bytes=_bytes_per_sample_id(h.sample_ct)
    total=_counts();parts=[];base=-1
    for prior in range(start-1,-1,-1):
        if int(h.vrtypes[prior]) not in (2,3):base=prior;break
    with path.open('rb') as file,mmap.mmap(file.fileno(),0,access=mmap.ACCESS_READ) as mm:
        for i in range(start,stop):
            if include_chunks and (i-start)%chunk_markers==0:
                current=_counts();parts.append((i,min(i+chunk_markers,stop),current))
            targets=(total,current) if include_chunks else (total,)
            form=int(h.vrtypes[i])
            if form in (2,3):
                if base<0:raise PgenFormatError('LD record has no earlier non-LD base')
                if (i-start)%chunk_markers==0:
                    replay=_counts()
                    _accumulate_record(mm,h,base,id_bytes,(replay,))
                    prefix_bytes=int(h.record_offsets[i])-int(h.record_offsets[base])
                    for state in targets:
                        state['ld_starts']+=1;state['replay_bytes']+=int(h.record_lengths[base])
                        state['native_replays'].append((i,base,prefix_bytes,replay))
            else:base=i
            _accumulate_record(mm,h,i,id_bytes,targets)
    file_bytes=path.stat().st_size
    index_bytes=file_bytes-int(h.record_lengths.astype(np.uint64).sum())

    def report(lo,hi,state):
        forms=dict(state['forms'])
        work=pgen_structure(h.sample_ct,forms,difflist_entries=state['entries'],varint_bytes=state['varints'])
        value=dict(path=str(path),samples=h.sample_ct,markers=hi-lo,file_markers=h.variant_ct,
            variant_range=[lo,hi],file_bytes=file_bytes,index_bytes=index_bytes,
            record_payload_bytes=state['payload'],record_form_counts=forms,difflist_by_form=state['by_form'],
            total_difflist_groups=state['groups'],sample_delta_integer_count=state['delta_integers'],source_work=work,
            native_onebit_tail_high_count=state['tail_high'],difflist_header_varint_lengths=dict(state['headers']),
            sample_delta_varint_lengths=dict(state['deltas']),total_varint_lengths=dict(state['headers']+state['deltas']),
            chunk_markers=chunk_markers,ld_records_at_chunk_starts=state['ld_starts'],
            additional_base_record_bytes_if_every_chunk_restarts=state['replay_bytes'],
            replay_interpretation='Potential base read obligations; actual reader-local cached-base reuse must be derived from scheduling',
            scope='Exact encoded counts for basic hardcall records; no elapsed-time coefficients, no genotype association')
        if state['native_replays']:
            value['native_ld_replays']=[dict(chunk_start=at,base_variant=base,
                read_prefix_bytes=amount,skipped_prefix_ld_records=at-base-1,
                base_record=report(base,base+1,counts))
                for at,base,amount,counts in state['native_replays']]
        return value

    result=report(start,stop,total)
    if include_chunks:
        result['chunks']=[report(lo,hi,state) for lo,hi,state in parts]
        result['chunk_census_scope']='One encoded-file pass, exact per-chunk counts and contiguous file-global ranges; no timing observations.'
    return result


def rechunk_census(encoded, chunk_markers, variant_range=None):
    """Derive a coarser aligned census from exact source-count chunks.

    Work counts are additive. LD restart work is not: keep only the restart
    ledger at each new chunk start, and retain that base's original source
    counts. No file reads, decode timings or uniform-work approximation enter.
    Returned records are independent copies; callers may mutate either view.
    """
    return _regroup_census(encoded,chunk_markers,variant_range,None)


def scheduled_census(encoded, chunk_markers, chunk_ranges, *, max_source_chunks=1000000):
    """Regroup a fine census for actual adaptive ranges without reopening input.

    The fixed ring capacity is chunk_markers. Ranges must be contiguous and
    aligned with the input census; only its final range may end off-grid.
    Both source and derived records remain independent, immutable-by-caller
    values. A subrange can describe remaining work, but does not imply that
    setup, in-flight work or reader initialization have already been paid.
    """
    if type(max_source_chunks) is not int or max_source_chunks<1:
        raise ValueError('Positive source-chunk budget required')
    if (not isinstance(chunk_ranges,list) or not chunk_ranges or len(chunk_ranges)>max_source_chunks
        or not isinstance(encoded.get('chunks'),list) or len(encoded['chunks'])>max_source_chunks):
        raise ValueError('Explicit schedule exceeds source-chunk budget or lacks source chunks')
    # Validate each row before indexing the first/last endpoints.
    for row in chunk_ranges:
        if not isinstance(row,(list,tuple)) or len(row)!=2 or any(type(v) is not int for v in row):
            raise ValueError('Explicit chunk ranges require integer endpoint pairs')
    return _regroup_census(encoded,chunk_markers,(chunk_ranges[0][0],chunk_ranges[-1][1]),chunk_ranges)


def _regroup_census(encoded,chunk_markers,variant_range,explicit_ranges):
    from copy import deepcopy
    from .decoder_work import census_chunk_ranges
    if 'chunk_ranges' in encoded:
        raise ValueError('Regrouping requires a regular fine source census')
    fine=encoded['chunk_markers'];parent_lo,parent_hi=encoded['variant_range']
    if encoded['markers']!=parent_hi-parent_lo:raise ValueError('Source census range does not match its marker count')
    if (type(chunk_markers) is not int or chunk_markers<1 or type(fine) is not int
        or fine<1 or chunk_markers%fine):
        raise ValueError('Rechunking requires a positive multiple of the source chunk size')
    lo,hi=(parent_lo,parent_hi) if variant_range is None else variant_range
    if (any(type(v) is not int for v in (lo,hi)) or not parent_lo<=lo<hi<=parent_hi
        or (lo-parent_lo)%fine or (hi!=parent_hi and (hi-parent_lo)%fine)):
        raise ValueError('Rechunking requires an aligned in-parent variant range')
    layout=dict(markers=hi-lo,variant_range=[lo,hi],chunk_markers=chunk_markers)
    if explicit_ranges is not None:layout['chunk_ranges']=explicit_ranges
    ranges=census_chunk_ranges(layout,chunk_markers)
    if any((start-parent_lo)%fine or (end!=parent_hi and (end-parent_lo)%fine) for start,end in ranges):
        raise ValueError('Scheduled chunks require aligned source-count boundaries')
    starts={start for start,_ in ranges}
    chunks=encoded.get('chunks')
    if not isinstance(chunks,list) or len(chunks)!=(parent_hi-parent_lo+fine-1)//fine:
        raise ValueError('Rechunking requires a complete source chunk census')
    context=('path','samples','file_markers','file_bytes','index_bytes','replay_interpretation','scope')
    scalars=('record_payload_bytes','total_difflist_groups','sample_delta_integer_count','native_onebit_tail_high_count')
    counters=('record_form_counts','difflist_header_varint_lengths','sample_delta_varint_lengths','total_varint_lengths')
    for i,part in enumerate(chunks):
        start=parent_lo+i*fine;end=min(start+fine,parent_hi)
        if (part.get('variant_range')!=[start,end] or part.get('markers')!=end-start
            or part.get('chunk_markers')!=fine or any(part.get(k)!=encoded.get(k) for k in context)):
            raise ValueError('Source chunk ranges or file context mismatch')
        if part.get('ld_records_at_chunk_starts',0) and not part.get('native_ld_replays'):
            raise ValueError('Exact native LD replay ledger required for rechunking')

    def merge(parts,start,end):
        value={key:deepcopy(encoded[key]) for key in context}
        value.update(markers=end-start,variant_range=[start,end],chunk_markers=chunk_markers)
        value.update({key:sum(part[key] for part in parts) for key in scalars})
        for key in counters:
            count=Counter()
            for part in parts:count.update(part[key])
            value[key]=dict(count)
        by_form={}
        for part in parts:
            for form,amounts in part['difflist_by_form'].items():
                row=by_form.setdefault(form,dict.fromkeys(('records','entries','varint_bytes','groups'),0))
                for key,amount in amounts.items():row[key]+=amount
        value['difflist_by_form']=by_form
        value['source_work']=pgen_structure(value['samples'],value['record_form_counts'],
            difflist_entries=sum(row['entries'] for row in by_form.values()),
            varint_bytes=sum(row['varint_bytes'] for row in by_form.values()))
        replays=[]
        for part in parts:
            # A constituent chunk's LD prefix is irrelevant unless the new
            # executor also starts a chunk at that exact position.
            if part['variant_range'][0] not in starts:continue
            for replay in part.get('native_ld_replays',[]):
                row=deepcopy(replay);row['base_record']['chunk_markers']=chunk_markers
                replays.append(row)
        value['ld_records_at_chunk_starts']=len(replays)
        value['additional_base_record_bytes_if_every_chunk_restarts']=sum(
            row['base_record']['record_payload_bytes'] for row in replays)
        if replays:value['native_ld_replays']=replays
        return value

    selected=chunks[(lo-parent_lo)//fine:(hi-parent_lo+fine-1)//fine]
    result=merge(selected,lo,hi)
    result['chunks']=[]
    for start,end in ranges:
        parts=chunks[(start-parent_lo)//fine:(end-parent_lo+fine-1)//fine]
        result['chunks'].append(merge(parts,start,end))
    result['chunk_census_scope']=encoded['chunk_census_scope']
    if explicit_ranges is not None:result['chunk_ranges']=ranges
    return result
