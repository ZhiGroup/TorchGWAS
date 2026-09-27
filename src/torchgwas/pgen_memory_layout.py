"""PGEN reader memory extents from locators, without reading record payloads.

These records deliberately lack decoder-operation counts. They can admit memory
but cannot substitute for an exact work census in the runtime calculator.
"""
from copy import deepcopy
from pathlib import Path
import numpy as np
from .pgen_reader import read_header,PgenFormatError

KIND='torchgwas.pgen_memory_layout.v1'
COMPACT_KIND='torchgwas.pgen_compact_memory_envelope.v1'
_ALLOWED_HARDCALL_TYPES=np.zeros(256,dtype=np.bool_)
_ALLOWED_HARDCALL_TYPES[[0,1,2,3,4,6,7]]=True


def memory_layout(path, chunk_markers, *, header=None, _validated_bases_receiver=None,
                  compact=False, _defer_bases=False):
    if type(chunk_markers) is not int or chunk_markers<1:
        raise ValueError('Positive integer memory-layout chunk size required')
    if type(_defer_bases) is not bool or (_defer_bases and (not compact or _validated_bases_receiver is not None)):
        raise ValueError('Deferred bases require compact admission without an index receiver')
    path=Path(path).resolve(strict=True)
    h=read_header(path) if header is None else header
    if h.variant_ct<1 or h.sample_ct<1:
        raise ValueError('Nonempty genotype dimensions required')
    if not np.all(_ALLOWED_HARDCALL_TYPES[h.vrtypes]):
        raise PgenFormatError('Memory layout requires biallelic unphased hardcalls')
    file_bytes=path.stat().st_size
    final=int(h.record_offsets[-1])+int(h.record_lengths[-1])
    if final>file_bytes or np.any(h.record_lengths==0):
        raise PgenFormatError('PGEN record extents exceed the file or are empty')
    payload=int(h.record_lengths.sum(dtype=np.uint64))
    if file_bytes-payload<12:raise PgenFormatError('Invalid PGEN index extent')
    if int(h.vrtypes[0]) in (2,3):raise PgenFormatError('LD record has no earlier non-LD base')
    bases=(None if _defer_bases else np.flatnonzero((h.vrtypes!=2)&(h.vrtypes!=3)))
    if _validated_bases_receiver is not None:
        if not callable(_validated_bases_receiver):
            raise ValueError('Validated base receiver must be callable')
        _validated_bases_receiver(bases)
    context=dict(kind=KIND,path=str(path),samples=int(h.sample_ct),file_markers=int(h.variant_ct),
        file_bytes=file_bytes,index_bytes=file_bytes-payload,chunk_markers=chunk_markers)
    packed=(h.sample_ct+3)//4
    # Locate all chunk starts in one vector operation. Scalar searchsorted at
    # every LD-compressed start is costly on multi-million-record indexes.
    width=min(chunk_markers,h.variant_ct)
    starts=np.arange(0,h.variant_ct,width,dtype=np.int64)
    ends=np.minimum(starts+width,h.variant_ct)
    offsets=h.record_offsets[starts]
    end_offsets=np.empty(len(starts),dtype=np.uint64)
    end_offsets[:-1]=h.record_offsets[ends[:-1]]
    end_offsets[-1]=final
    amounts=end_offsets-offsets
    prefixes=np.zeros(len(starts),dtype=np.uint64)
    extras=np.zeros(len(starts),dtype=np.uint64)
    ld=(h.vrtypes[starts]==2)|(h.vrtypes[starts]==3)
    if np.any(ld):
        at=np.flatnonzero(ld)
        if bases is None:
            base=starts[at].copy()
            # Only fine-grid starts need predecessors for startup memory bounds.
            # Typical LD replay is short. A long run falls back to the exact
            # full index after a bounded number of vector passes.
            for _ in range(16):
                pending=(h.vrtypes[base]==2)|(h.vrtypes[base]==3)
                if not np.any(pending):break
                base[pending]-=1
            else:
                # Only a few fine starts normally exceed 16 replay records.
                # Bound scalar metadata work before falling back to the full
                # index on an adversarial file with long LD runs.
                remaining=np.flatnonzero((h.vrtypes[base]==2)|(h.vrtypes[base]==3))
                steps=0;exhausted=False
                for index in remaining:
                    position=int(base[index])
                    while int(h.vrtypes[position]) in (2,3):
                        position-=1;steps+=1
                        if steps>1_000_000:
                            exhausted=True;break
                    if exhausted:break
                    base[index]=position
                if exhausted:
                    bases=np.flatnonzero((h.vrtypes!=2)&(h.vrtypes!=3))
                    base=bases[np.searchsorted(bases,starts[at],side='right')-1]
        else:
            base=bases[np.searchsorted(bases,starts[at],side='right')-1]
        prefixes[at]=offsets[at]-h.record_offsets[base]
        extras[at]=packed+8*(starts[at]-base)
    if compact:
        cumulative=np.empty(len(amounts)+1,dtype=np.uint64)
        cumulative[0]=0
        np.cumsum(amounts,dtype=np.uint64,out=cumulative[1:])
        if int(cumulative[-1])!=payload:
            raise PgenFormatError('Compact PGEN payload does not conserve index lengths')
        for array in (amounts,prefixes,extras,cumulative):array.flags.writeable=False
        return dict(context,kind=COMPACT_KIND,markers=int(h.variant_ct),
            variant_range=[0,int(h.variant_ct)],record_payload_bytes=payload,
            fine_payload_bytes=amounts,fine_prefix_bytes=prefixes,
            fine_extra_workspace_bytes=extras,payload_cumulative_bytes=cumulative,
            logical_chunks=len(amounts),
            scope='Full-file validated header vectors for exact fixed and shifted reader-memory maxima; no per-chunk Python records or genotype payload reads.')
    parts=[dict(context,markers=int(hi-lo),variant_range=[int(lo),int(hi)],
                record_payload_bytes=int(amount),read_prefix_bytes=int(prefix),
                reader_extra_workspace_bytes=int(extra))
           for lo,hi,amount,prefix,extra in zip(starts,ends,amounts,prefixes,extras)]
    return dict(context,markers=int(h.variant_ct),variant_range=[0,int(h.variant_ct)],
        record_payload_bytes=payload,chunks=parts,
        scope='Header/record-locator memory extents only; zero genotype payload reads and no decoder work counts.')


def _compact_span(layout,chunk_markers,variant_range):
    if layout.get('kind')!=COMPACT_KIND or type(chunk_markers) is not int or chunk_markers<1:
        raise ValueError('Typed compact layout and positive chunk capacity required')
    fine=layout['chunk_markers'];m=layout['file_markers']
    if type(fine) is not int or fine<1 or chunk_markers%fine:
        raise ValueError('Compact capacity must be a multiple of its fine grid')
    lo,hi=variant_range
    if (any(type(v) is not int for v in (lo,hi)) or not 0<=lo<hi<=m
            or lo%fine or (hi%fine and hi!=m)):
        raise ValueError('Compact range must align to the fine source grid')
    first=lo//fine;last=(hi+fine-1)//fine
    arrays=[layout[key] for key in ('fine_payload_bytes','fine_prefix_bytes',
        'fine_extra_workspace_bytes')]
    cumulative=layout['payload_cumulative_bytes']
    expected=(m+fine-1)//fine
    if (any(not isinstance(a,np.ndarray) or a.ndim!=1 or len(a)!=expected
            or a.flags.writeable for a in arrays)
            or not isinstance(cumulative,np.ndarray) or cumulative.ndim!=1
            or len(cumulative)!=expected+1 or cumulative.flags.writeable
            or int(cumulative[0])!=0 or int(cumulative[-1])!=layout['record_payload_bytes']):
        raise ValueError('Incomplete or mutable compact source vectors')
    return fine,first,last


def compact_rechunk_memory_layout(layout,chunk_markers,variant_range=None):
    """Exact maxima for regular fixed-capacity reads in one aligned shard."""
    span=layout['variant_range'] if variant_range is None else variant_range
    fine,first,last=_compact_span(layout,chunk_markers,span)
    step=chunk_markers//fine
    starts=np.arange(first,last,step,dtype=np.int64)
    stops=np.minimum(starts+step,last)
    cumulative=layout['payload_cumulative_bytes']
    read=cumulative[stops]-cumulative[starts]+layout['fine_prefix_bytes'][starts]
    read_max=int(read.max());extra_max=int(layout['fine_extra_workspace_bytes'][starts].max())
    lo,hi=span
    context={key:layout[key] for key in ('path','samples','file_markers','file_bytes','index_bytes')}
    context.update(kind=KIND,chunk_markers=chunk_markers)
    # A single synthetic row carries the two maxima independently. It is
    # valid only for memory admission; decoder/runtime paths reject KIND.
    row=dict(context,markers=min(chunk_markers,hi-lo),variant_range=[lo,min(lo+chunk_markers,hi)],
        record_payload_bytes=read_max,read_prefix_bytes=0,
        reader_extra_workspace_bytes=extra_max)
    return dict(context,markers=hi-lo,variant_range=[lo,hi],
        record_payload_bytes=int(cumulative[last]-cumulative[first]),chunks=[row],
        logical_chunks=len(starts),
        scope='Exact fixed-grid reader-memory maxima only; synthetic summary row cannot price decoder work.')


def compact_shifted_reader_envelope(layout,variant_range,capacity):
    """Exact maxima over every fine-grid start after a chunk-size change."""
    fine,first,last=_compact_span(layout,capacity,variant_range)
    starts=np.arange(first,last,dtype=np.int64)
    stops=np.minimum(starts+capacity//fine,last)
    cumulative=layout['payload_cumulative_bytes']
    read=cumulative[stops]-cumulative[starts]+layout['fine_prefix_bytes'][starts]
    return int(read.max()),int(layout['fine_extra_workspace_bytes'][first:last].max())


def rechunk_memory_layout(layout, chunk_markers, variant_range=None):
    """Group locators and retain only the first row's LD restart obligation."""
    if layout.get('kind')!=KIND or type(chunk_markers) is not int or chunk_markers<1:
        raise ValueError('Typed memory layout and positive integer capacity required')
    fine=layout['chunk_markers'];start,stop=layout['variant_range']
    if type(fine) is not int or fine<1 or chunk_markers%fine:
        raise ValueError('Memory layout requires a coarser aligned grid')
    lo,hi=(start,stop) if variant_range is None else variant_range
    if any(type(v) is not int for v in (lo,hi)) or not start<=lo<hi<=stop or (lo-start)%fine or (hi!=stop and (hi-start)%fine):
        raise ValueError('Memory layout range must align to its source grid')
    parts=layout['chunks'];result=[]
    context={key:deepcopy(value) for key,value in layout.items()
             if key not in ('chunks','markers','variant_range','record_payload_bytes')}
    context['chunk_markers']=chunk_markers
    for first in range(lo,hi,chunk_markers):
        last=min(first+chunk_markers,hi)
        group=parts[(first-start)//fine:(last-start+fine-1)//fine]
        if not group or group[0]['variant_range'][0]!=first or group[-1]['variant_range'][1]!=last:
            raise ValueError('Incomplete memory layout grouping')
        result.append(dict(context,markers=last-first,variant_range=[first,last],
            record_payload_bytes=sum(row['record_payload_bytes'] for row in group),
            read_prefix_bytes=group[0]['read_prefix_bytes'],
            reader_extra_workspace_bytes=group[0]['reader_extra_workspace_bytes']))
    return dict(context,markers=hi-lo,variant_range=[lo,hi],chunks=result,
                record_payload_bytes=sum(row['record_payload_bytes'] for row in result))


def reader_layout(row):
    values=[row[key] for key in ('record_payload_bytes','read_prefix_bytes','reader_extra_workspace_bytes')]
    if row.get('kind')!=KIND or any(type(v) is not int or v<0 for v in values):
        raise ValueError('Invalid typed PGEN memory extent')
    return dict(read_bytes=values[0]+values[1],extra_workspace_bytes=values[2],
                decode_input_bytes=None)
