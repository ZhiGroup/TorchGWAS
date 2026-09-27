"""Conditional native decoder work intervals from PGEN index metadata only.

These are source counts, never exact censuses or observed capacities. Scope is
the well-formed 32-bit hardcall encoding accepted by decoder_work: each ULEB
occupies 1..5 bytes and the record has no unconsumed trailing payload. Overlong
but <=5-byte integers are allowed. Header inspection cannot establish payload
validity; decoding/census must still reject corrupt or unsupported payloads.
"""
from __future__ import annotations

from collections import Counter, OrderedDict
from copy import deepcopy
from functools import lru_cache
import math
from pathlib import Path

import numpy as np

from .analytical_plan_cache import input_identity
from .decoder_work import _native_header_units, _native_decoder_memory_bytes, price_identified_units
from .first_principles import positive
from .input_read_work import buffered_read_service
from .pgen_reader import PgenFormatError, PgenHeader, _bytes_per_sample_id, ld_safe_start, read_header

KIND = 'torchgwas.pgen_header_work_bounds.v1'
WINDOW_KIND = 'torchgwas.pgen_header_window.v1'
SCHEDULE_KIND = 'torchgwas.pgen_header_schedule_bounds.v1'
_ALLOWED_HARDCALL_TYPES=np.zeros(256,dtype=np.bool_)
_ALLOWED_HARDCALL_TYPES[[0,1,2,3,4,6,7]]=True


def _integer(name, value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def _structure(entries, id_bytes):
    groups = (entries + 63) // 64
    return groups * id_bytes + max(groups - 1, 0) + (entries + 3) // 4


def _lower_bound(stop, predicate):
    """First true on [0, stop), or stop. No enumeration over sample count."""
    lo, hi = 0, stop
    while lo < hi:
        mid = (lo + hi) // 2
        if predicate(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def difflist_bounds(samples, payload_bytes):
    """Bounds for a difflist after removing the one-bit fixed prefix.

    For E entries, G=ceil(E/64), D=E-G, A=G*id+max(G-1,0)+ceil(E/4).
    The byte length L satisfies h_min(E)+A+D <= L <= 5+A+5D.
    Both sides are monotone, so two binary searches bound E. Histogram bounds
    also enforce V=1+D and T=L-A; their endpoints need not be jointly feasible.
    """
    _integer('samples', samples, 1)
    _integer('payload_bytes', payload_bytes, 1)
    if samples > 0xffffffff:
        raise ValueError('PGEN samples must fit uint32')
    width = _bytes_per_sample_id(samples)

    def minimum(e):
        return max(1, (e.bit_length() + 6) // 7) + _structure(e, width) + e - (e + 63) // 64

    def maximum(e):
        return 5 + _structure(e, width) + 5 * (e - (e + 63) // 64)

    low = _lower_bound(samples + 1, lambda e: maximum(e) >= payload_bytes)
    high = _lower_bound(samples + 1, lambda e: minimum(e) > payload_bytes) - 1
    if low > high or high < 0 or low > samples:
        raise PgenFormatError('Record length cannot contain a supported 32-bit difflist')
    groups = [(low + 63) // 64, (high + 63) // 64]
    integers = [1 + low - groups[0], 1 + high - groups[1]]
    variable_bytes = [payload_bytes - _structure(high, width),
                      payload_bytes - _structure(low, width)]
    units = dict(set_category=[low, high], difflist_category_extract=[low, high],
                 difflist_group_absolute_id=groups)
    for length in range(1, 6):
        # Relative to all-one-byte or all-five-byte integers, each occurrence
        # consumes length-1 extra bytes or 5-length missing bytes respectively.
        upper = integers[1]
        if length > 1:
            upper = min(upper, (variable_bytes[1] - integers[0]) // (length - 1))
        if length < 5:
            upper = min(upper, (5 * integers[1] - variable_bytes[0]) // (5 - length))
        lower = (max(0, 2 * integers[0] - variable_bytes[1]) if length == 1
                 else max(0, variable_bytes[0] - 4 * integers[1]) if length == 5 else 0)
        units[f'uleb{length}'] = [lower, max(0, upper)]
    return dict(entries=[low, high], groups=groups, varint_integers=integers,
                varint_bytes=variable_bytes, source_units=units)


def _add(target, units, count=1):
    for name, (lo, hi) in units.items():
        if hi:
            old = target.setdefault(name, [0, 0])
            old[0] += count * lo
            old[1] += count * hi


def _record_units(samples, form, length, *, expand):
    fixed, base_bytes = _native_header_units(samples, {form: 1}, expand=expand)
    units = {name: [value, value] for name, value in fixed.items()}
    varints = dict(varint_integers=[0, 0], varint_bytes=[0, 0])
    if form == 0:
        if length != (samples + 3) // 4:
            raise PgenFormatError('Unexpected plain record length')
    else:
        prefix = 1 + (samples + 7) // 8 if form == 1 else 0
        if length <= prefix:
            raise PgenFormatError('Missing difflist header')
        diff = difflist_bounds(samples, length - prefix)
        _add(units, diff['source_units'])
        varints = {key: diff[key] for key in varints}
        if form == 1 and samples % 8:
            _add(units, {'set_category': [0, samples % 8]})
    return units, base_bytes, varints


class PgenHeaderWork:
    """Read the index once, then bound explicitly budgeted variant ranges.

    Construction is O(file variants) metadata work, not a free startup action.
    Each bounds() call caps both records and distinct (form,length) signatures;
    it can be scheduled as one post-output planning step. No genotype payload
    is read, mmap'ed, or decoded. The object is bound to filesystem identity.
    """
    def __init__(self, path, *, max_cached_signatures=1024, _prepared_header=None,
                 _prepared_index=None, max_cached_bounds=64):
        _integer('max_cached_signatures',max_cached_signatures)
        _integer('max_cached_bounds',max_cached_bounds)
        if max_cached_bounds>256:raise ValueError('max_cached_bounds exceeds 256')
        # Each result depends only on sample count, encoding form, byte length
        # and expansion mode. Overlapping horizons repeatedly use these same
        # structural calculations. Cached rows stay private: _add creates new
        # aggregate lists rather than exposing these mutable source intervals.
        self._record_work=lru_cache(maxsize=max_cached_signatures)(_record_units)
        self._bounds_cache=OrderedDict();self._max_cached_bounds=max_cached_bounds
        self._bounds_hits=self._bounds_misses=0
        self.path = Path(path).resolve(strict=True)
        self.input_identity = input_identity(self.path)
        if _prepared_header is None:
            self._header = read_header(self.path)
        else:
            # The admission path parsed this index before the first useful
            # chunk. Reuse its exact file binding; never accept an unbound or
            # stale array of record locators as fresh source work.
            if (not isinstance(_prepared_header,tuple) or len(_prepared_header)!=2
                    or _prepared_header[0]!=self.input_identity
                    or not isinstance(_prepared_header[1],PgenHeader)):
                raise ValueError('Prepared PGEN header differs from current input')
            self._header = _prepared_header[1]
        h = self._header
        if h.variant_ct < 1 or h.sample_ct < 1:
            raise PgenFormatError('Nonempty genotype dimensions required')
        if _prepared_index is None:
            if not np.all(_ALLOWED_HARDCALL_TYPES[h.vrtypes]):
                raise PgenFormatError('Work bounds require biallelic unphased hardcalls')
            final = int(h.record_offsets[-1]) + int(h.record_lengths[-1])
            if final > self.input_identity['bytes'] or np.any(h.record_lengths == 0):
                raise PgenFormatError('PGEN record extents exceed the file or are empty')
            if int(h.vrtypes[0]) in (2,3):
                raise PgenFormatError('LD record has no earlier non-LD base')
            # Cold windows inspect at most a bounded set of starts; there is
            # no need to build a full-file predecessor array for them.
            self._bases = None
        else:
            # Memory admission already made exactly these full-index checks.
            # Its certificate is private and tied to this same retained header
            # and current file identity, so the first productive window need
            # only inspect the bounded range it actually prices.
            if (not isinstance(_prepared_index,tuple) or len(_prepared_index)!=3
                    or _prepared_header is None or _prepared_index[0]!=self.input_identity
                    or _prepared_index[1] is not h
                    or not isinstance(_prepared_index[2],np.ndarray)
                    or _prepared_index[2].ndim!=1
                    or _prepared_index[2].dtype not in (np.dtype('uint32'),np.dtype('int64'))
                    or _prepared_index[2].flags.writeable or not len(_prepared_index[2])
                    or int(_prepared_index[2][0])!=0
                    or int(_prepared_index[2][-1])>=h.variant_ct):
                raise ValueError('Prepared PGEN index differs from admitted source')
            self._bases = _prepared_index[2]
        self._unchanged()

    def cache_info(self):
        """Per-index structural reuse only; no rates or empirical observations."""
        return self._record_work.cache_info()._asdict()

    def bounds_cache_info(self):
        return dict(hits=self._bounds_hits,misses=self._bounds_misses,
                    entries=len(self._bounds_cache),max_entries=self._max_cached_bounds)

    def _unchanged(self):
        if input_identity(self.path) != self.input_identity:
            raise ValueError('PGEN input changed after header inspection')

    def window(self, start, stop, chunk_markers, *, max_records=65536,
               max_chunks=8, max_signatures=256):
        """Bound one finite pricing window, never expand the remaining job.

        The already-open index still cost O(file variants) to construct. This
        operation reads no payload and refuses excessive work before collecting
        any per-chunk bounds. Original file coordinates and index cost survive.
        """
        for name,value in [('start',start),('stop',stop)]:_integer(name,value)
        for name,value in [('chunk_markers',chunk_markers),('max_records',max_records),
            ('max_chunks',max_chunks),('max_signatures',max_signatures)]:_integer(name,value,1)
        if not 0<=start<stop<=self._header.variant_ct or stop-start>max_records:
            raise ValueError('Header window exceeds its record budget or input range')
        count=(stop-start+chunk_markers-1)//chunk_markers
        if count>max_chunks:raise ValueError('Header window exceeds its chunk budget')
        self._unchanged()
        rows=[self.bounds(lo,min(lo+chunk_markers,stop),max_records=max_records,
            max_signatures=max_signatures) for lo in range(start,stop,chunk_markers)]
        self._unchanged()
        return dict(kind=WINDOW_KIND,path=str(self.path),samples=int(self._header.sample_ct),
            markers=stop-start,file_markers=int(self._header.variant_ct),variant_range=[start,stop],
            chunk_markers=chunk_markers,chunk_ranges=[r['variant_range'] for r in rows],chunks=rows,
            input_identity=deepcopy(self.input_identity),
            structural_work=dict(records=stop-start,chunks=count,payload_bytes_read=0),
            scope='Bounded header-only source window in original file coordinates, not an exact census or a representative-sampling guarantee.')

    def vectorized_window(self, start, stop, chunk_markers, *, max_records=10000000,
                          max_chunks=100000, max_signatures=32768,
                          max_global_signatures=32768, max_signature_pairs=10000000):
        """Materialize exact chunk geometry with one global signature pass.

        Every chunk retains the same conditional per-record intervals and LD
        replay as ``window``. This is a post-output planning path: its bounded
        arrays are proportional to the requested variant and signature-pair
        counts. It does not read genotype payload or price elapsed service.
        """
        for name,value in [('start',start),('stop',stop)]:_integer(name,value)
        for name,value in [('chunk_markers',chunk_markers),('max_records',max_records),
                           ('max_chunks',max_chunks),('max_signatures',max_signatures),
                           ('max_global_signatures',max_global_signatures),
                           ('max_signature_pairs',max_signature_pairs)]:_integer(name,value,1)
        h=self._header
        if not 0<=start<stop<=h.variant_ct or stop-start>max_records:
            raise ValueError('Vectorized window exceeds its record budget or file')
        markers=stop-start;chunks=(markers+chunk_markers-1)//chunk_markers
        if chunks>max_chunks:raise ValueError('Vectorized window exceeds its chunk budget')
        self._unchanged()
        raw=(h.vrtypes[start:stop].astype(np.uint64)<<np.uint64(32))|h.record_lengths[start:stop]
        signatures,inverse=np.unique(raw,return_inverse=True)
        del raw
        if len(signatures)>max_global_signatures:
            raise ValueError('Vectorized window exceeds its global signature budget')
        if chunks*len(signatures)>np.iinfo(np.int64).max:
            raise ValueError('Vectorized signature index overflows int64')
        # The admitted chunk/signature product usually fits uint32. Keep the
        # large per-variant arrays narrow and release each as soon as its one
        # indexing purpose ends; phenotype memory may already occupy the host.
        code_dtype=(np.uint32 if chunks*len(signatures)<=np.iinfo(np.uint32).max
                    else np.int64)
        encoded_inverse=inverse.astype(code_dtype)
        del inverse
        chunk_ids=np.arange(markers,dtype=code_dtype)//chunk_markers
        combined=chunk_ids*len(signatures)+encoded_inverse
        del chunk_ids,encoded_inverse
        pairs,counts=np.unique(combined,return_counts=True)
        del combined
        if len(pairs)>max_signature_pairs:
            raise ValueError('Vectorized window exceeds its signature-pair budget')
        cut=np.searchsorted(pairs,np.arange(chunks,dtype=np.int64)*len(signatures))
        distinct=np.diff(np.append(cut,len(pairs)))
        if np.any(distinct>max_signatures):
            raise ValueError('Requested work exceeds the distinct signature budget')
        codes=pairs%len(signatures)

        def sums(values):
            lookup=np.asarray(values,dtype=np.int64)
            weighted=np.multiply(counts,lookup[codes],dtype=np.int64)
            return np.add.reduceat(weighted,cut,dtype=np.int64)

        primary=[];forms=[]
        n=int(h.sample_ct);packed=(n+3)//4
        for signature in signatures:
            form,length=int(signature)>>32,int(signature)&0xffffffff
            row,base_bytes,integers=self._record_work(n,form,length,expand=True)
            primary.append((row,base_bytes,integers));forms.append(form)
        names=sorted({name for row,_,_ in primary for name in row})
        unit_vectors={name:(sums([row.get(name,[0,0])[0] for row,_,_ in primary]),
                            sums([row.get(name,[0,0])[1] for row,_,_ in primary]))
                      for name in names}
        constraints={name:(sums([integers[name][0] for _,_,integers in primary]),
                           sums([integers[name][1] for _,_,integers in primary]))
                     for name in ('varint_integers','varint_bytes')}
        form_vectors={form:sums([int(value==form) for value in forms]) for form in sorted(set(forms))}
        base_updates=sums([value for _,value,_ in primary])
        del pairs,counts,codes,cut,primary,forms,signatures
        offsets=np.arange(0,markers,chunk_markers,dtype=np.int64)
        payloads=np.add.reduceat(h.record_lengths[start:stop].astype(np.int64),offsets,dtype=np.int64)
        rows=[]
        for i,relative in enumerate(offsets):
            lo=start+int(relative);hi=min(lo+chunk_markers,stop);size=hi-lo
            units={name:[int(pair[0][i]),int(pair[1][i])] for name,pair in unit_vectors.items()
                   if pair[1][i]}
            varints={name:[int(pair[0][i]),int(pair[1][i])] for name,pair in constraints.items()
                     if pair[1][i]}
            payload=int(payloads[i]);replay=None;extra_base=0;prefix=0
            if int(h.vrtypes[lo]) in (2,3):
                base=(ld_safe_start(h.vrtypes,lo) if self._bases is None else
                      int(self._bases[np.searchsorted(self._bases,self._bases.dtype.type(lo),side='right')-1]))
                form,length=int(h.vrtypes[base]),int(h.record_lengths[base])
                replay_units,extra_base,replay_integers=self._record_work(n,form,length,expand=False)
                _add(units,replay_units);_add(varints,replay_integers)
                prefix=int(h.record_offsets[lo])-int(h.record_offsets[base])
                replay=dict(base_variant=base,record_form=form,record_bytes=length,
                            read_prefix_bytes=prefix)
            rows.append(dict(kind=KIND,implementation='torch_native_int8',samples=n,
                markers=size,variant_range=[lo,hi],input_identity=deepcopy(self.input_identity),
                record_form_counts={form:int(values[i]) for form,values in form_vectors.items() if values[i]},
                record_payload_bytes=payload,source_units=units,varint_constraints=varints,
                native_ld_base_update_bytes=int(base_updates[i])+extra_base,
                native_ld_replay_packed_bytes=packed if replay else 0,ld_replay=replay,
                read_bytes=payload+prefix,
                decode_input_bytes=payload+(replay['record_bytes'] if replay else 0),
                logical_final_packed_bytes=packed*size,logical_final_int8_bytes=n*size,
                structural_work=dict(records=size,signatures=int(distinct[i]),
                                     replay_records=int(replay is not None),payload_bytes_read=0),
                scope='Conditional source-count intervals for valid basic hardcall records with 1..5-byte integers and no trailing payload. LD replay decodes only the base, without int8 expansion. Endpoints are not necessarily jointly attainable; no CPU capacity or elapsed time was measured.'))
        self._unchanged()
        return dict(kind=WINDOW_KIND,path=str(self.path),samples=n,markers=markers,
            file_markers=int(h.variant_ct),variant_range=[start,stop],chunk_markers=chunk_markers,
            chunk_ranges=[row['variant_range'] for row in rows],chunks=rows,
            input_identity=deepcopy(self.input_identity),
            structural_work=dict(records=markers,chunks=chunks,payload_bytes_read=0),
            scope='Bounded vectorized header-only source window in original file coordinates, not an exact census, priced service, or a representative sample.')

    def bounds(self, start, stop, *, max_records=65536, max_signatures=256):
        """One independently issued chunk, including its base-only LD replay."""
        for name, value in [('start', start), ('stop', stop)]:
            _integer(name, value)
        _integer('max_records', max_records, 1)
        _integer('max_signatures', max_signatures, 1)
        h = self._header
        if not 0 <= start < stop <= h.variant_ct or stop - start > max_records:
            raise ValueError('Requested variant range exceeds the record budget or file')
        self._unchanged()
        key=(start,stop,max_signatures)
        cached=self._bounds_cache.get(key)
        if cached is not None:
            self._bounds_hits+=1;self._bounds_cache.move_to_end(key)
            self._unchanged()
            return deepcopy(cached)
        self._bounds_misses+=1
        # A uint64 signature avoids a Python loop over every variant.
        keys = (h.vrtypes[start:stop].astype(np.uint64) << np.uint64(32)) | h.record_lengths[start:stop]
        signatures, counts = np.unique(keys, return_counts=True)
        if len(signatures) > max_signatures:
            raise ValueError('Requested work exceeds the distinct signature budget')
        units = {}; forms = Counter(); payload = base_updates = 0
        varints = {}
        n = int(h.sample_ct); packed = (n + 3) // 4
        for signature, count in zip(signatures, counts):
            form, length, count = int(signature) >> 32, int(signature) & 0xffffffff, int(count)
            row, base_bytes, integers = self._record_work(n, form, length, expand=True)
            _add(units, row, count); _add(varints, integers, count)
            forms[form] += count; payload += length * count; base_updates += base_bytes * count
        replay = None
        if int(h.vrtypes[start]) in (2, 3):
            # A Python int can make NumPy promote the entire uint32 base
            # index for each search. Match the retained array dtype exactly.
            base = (ld_safe_start(h.vrtypes,start) if self._bases is None else
                    int(self._bases[np.searchsorted(self._bases,self._bases.dtype.type(start),side='right')-1]))
            form, length = int(h.vrtypes[base]), int(h.record_lengths[base])
            row, base_bytes, integers = self._record_work(n, form, length, expand=False)
            _add(units, row); _add(varints, integers); base_updates += base_bytes
            replay = dict(base_variant=base, record_form=form, record_bytes=length,
                          read_prefix_bytes=int(h.record_offsets[start]) - int(h.record_offsets[base]))
        prefix = replay['read_prefix_bytes'] if replay else 0
        replay_payload = replay['record_bytes'] if replay else 0
        self._unchanged()
        result=dict(kind=KIND, implementation='torch_native_int8', samples=n, markers=stop-start,
            variant_range=[start, stop], input_identity=deepcopy(self.input_identity),
            record_form_counts=dict(forms), record_payload_bytes=payload,
            source_units=units, varint_constraints=varints, native_ld_base_update_bytes=base_updates,
            native_ld_replay_packed_bytes=packed if replay else 0, ld_replay=replay,
            read_bytes=payload+prefix, decode_input_bytes=payload+replay_payload,
            logical_final_packed_bytes=packed*(stop-start), logical_final_int8_bytes=n*(stop-start),
            structural_work=dict(records=stop-start, signatures=len(signatures),
                                 replay_records=int(replay is not None), payload_bytes_read=0),
            scope='Conditional source-count intervals for valid basic hardcall records with 1..5-byte integers and no trailing payload. LD replay decodes only the base, without int8 expansion. Endpoints are not necessarily jointly attainable; no CPU capacity or elapsed time was measured.')
        if self._max_cached_bounds:
            self._bounds_cache[key]=deepcopy(result)
            if len(self._bounds_cache)>self._max_cached_bounds:
                self._bounds_cache.popitem(last=False)
        return result

    def schedule_bounds(self, start, stop, chunk_markers, *, max_records=10000000,
                        max_signatures=32768, max_chunks=100000):
        """Aggregate an entire fixed chunk schedule from index metadata only.

        Each LD chunk start decodes its earlier non-LD base once. The primary
        records are aggregated by (form,length), and only restart positions
        are visited separately. This can be a post-output, bounded source
        continuation; it is not an exact payload census or a pipeline time.
        """
        for name,value in [('chunk_markers',chunk_markers),('max_records',max_records),
                           ('max_signatures',max_signatures),('max_chunks',max_chunks)]:
            _integer(name,value,1)
        for name,value in [('start',start),('stop',stop)]:_integer(name,value)
        h=self._header
        if not 0<=start<stop<=h.variant_ct or stop-start>max_records:
            raise ValueError('Requested schedule exceeds its record budget or file')
        chunks=(stop-start+chunk_markers-1)//chunk_markers
        if chunks>max_chunks:raise ValueError('Requested schedule exceeds its chunk budget')
        result=self.bounds(start,stop,max_records=max_records,max_signatures=max_signatures)
        common=dict(input_identity=deepcopy(result['input_identity']),
            variant_range=list(result['variant_range']),samples=result['samples'],
            markers=result['markers'],record_form_counts=deepcopy(result['record_form_counts']),
            record_payload_bytes=result['record_payload_bytes'],
            source_units=deepcopy(result['source_units']),
            read_bytes=result['read_bytes'],decode_input_bytes=result['decode_input_bytes'],
            native_ld_base_update_bytes=result['native_ld_base_update_bytes'],
            native_ld_replay_packed_bytes=result['native_ld_replay_packed_bytes'])
        starts=np.arange(start+chunk_markers,stop,chunk_markers,dtype=np.int64)
        ld=starts[np.isin(h.vrtypes[starts],(2,3))]
        if self._bases is None:
            bases=(int(ld_safe_start(h.vrtypes,int(at))) for at in ld)
        else:
            locations=np.searchsorted(self._bases,ld,side='right')-1
            bases=(int(value) for value in self._bases[locations])
        replay_rows=[]
        for at,base in zip(ld,bases):
            at=int(at)
            if not 0<=base<at or int(h.vrtypes[base]) in (2,3):
                raise PgenFormatError('LD schedule start has no valid earlier base')
            replay_rows.append([at,base,int(h.record_offsets[at])-int(h.record_offsets[base]),
                                int(h.vrtypes[base]),int(h.record_lengths[base])])
        replay_units={};extra_base_updates=extra_base_bytes=0
        extra_read_prefix=sum(row[2] for row in replay_rows)
        extra_packed=((int(h.sample_ct)+3)//4)*len(ld)
        if replay_rows:
            locations=np.fromiter((row[1] for row in replay_rows),dtype=np.int64,count=len(replay_rows))
            keys=(h.vrtypes[locations].astype(np.uint64)<<np.uint64(32))|h.record_lengths[locations]
            signatures,counts=np.unique(keys,return_counts=True)
            if len(signatures)>max_signatures:
                raise ValueError('Requested replay work exceeds the distinct signature budget')
            units=result['source_units'];varints=result['varint_constraints']
            for signature,count in zip(signatures,counts):
                form,length=int(signature)>>32,int(signature)&0xffffffff
                row,base_bytes,integers=self._record_work(int(h.sample_ct),form,length,expand=False)
                _add(units,row,int(count));_add(varints,integers,int(count))
                _add(replay_units,row,int(count))
                extra_base_updates+=base_bytes*int(count)
                extra_base_bytes+=length*int(count)
            result['native_ld_base_update_bytes']+=extra_base_updates
            result['native_ld_replay_packed_bytes']+=extra_packed
            result['decode_input_bytes']+=extra_base_bytes
            result['read_bytes']+=extra_read_prefix
        result['kind']=SCHEDULE_KIND
        result.pop('ld_replay')
        result['chunk_markers']=chunk_markers
        result['chunk_count']=chunks
        result['ld_replay_count']=result['structural_work']['replay_records']+len(ld)
        result['common_work']=common
        result['additional_replay_work']=dict(source_units=replay_units,
            read_bytes=extra_read_prefix,decode_input_bytes=extra_base_bytes,
            native_ld_base_update_bytes=extra_base_updates,
            native_ld_replay_packed_bytes=extra_packed,replay_count=len(ld),
            entries=replay_rows)
        result['structural_work']['replay_records']=result['ld_replay_count']
        result['structural_work']['chunks']=chunks
        result['scope']='Conditional whole-schedule source-count intervals from the PGEN index, including every fixed chunk-start LD base. Valid supported payload and independent per-chunk reader restart are assumed. No payload census, elapsed service or pipeline completion bound.'
        self._unchanged()
        return result


def _replay_entries_work(entries, samples, span, chunk_markers, *,
                         include_varint_constraints=False):
    if not isinstance(entries,list):raise ValueError('Explicit LD replay entries required')
    signatures=Counter();prefix=0;previous=-1
    for entry in entries:
        if (not isinstance(entry,(list,tuple)) or len(entry)!=5 or
                any(type(value) is not int for value in entry)):
            raise ValueError('Invalid LD replay entry')
        at,base,read_prefix,form,length=entry
        if (not span[0]<at<span[1] or (at-span[0])%chunk_markers or
                at<=previous or not 0<=base<at or form not in (0,1,4,6,7) or
                length<1 or read_prefix<length):
            raise ValueError('LD replay entry differs from its chunk schedule')
        previous=at;prefix+=read_prefix;signatures[(form,length)]+=1
    units={};varints={};base_updates=base_bytes=0
    for (form,length),count in signatures.items():
        row,updates,integers=_record_units(samples,form,length,expand=False)
        _add(units,row,count)
        if include_varint_constraints:_add(varints,integers,count)
        base_updates+=updates*count;base_bytes+=length*count
    result=dict(source_units=units,read_bytes=prefix,decode_input_bytes=base_bytes,
        native_ld_base_update_bytes=base_updates,
        native_ld_replay_packed_bytes=((samples+3)//4)*len(entries),
        replay_count=len(entries))
    if include_varint_constraints:result['varint_constraints']=varints
    return result


def paired_schedule_source_difference(baseline,candidate,cpu_seconds_per_unit):
    """Cancel identical primary source work before bounding a chunk-size move.

    Only additional LD restarts differ. Independent interval endpoints may
    still be loose, but common uncertain difflist work is not subtracted as if
    it were independent in the two layouts. This is source CPU/read work, not
    a paired pipeline makespan or a switch authorization.
    """
    for row in (baseline,candidate):
        if not isinstance(row,dict) or row.get('kind')!=SCHEDULE_KIND:
            raise ValueError('Two typed whole-header schedules required')
        common=row.get('common_work');extra=row.get('additional_replay_work')
        if not isinstance(common,dict) or not isinstance(extra,dict):
            raise ValueError('Complete invariant and replay work required')
        if any(row.get(key)!=common.get(key) for key in
                ('input_identity','variant_range','samples','markers',
                 'record_form_counts','record_payload_bytes')):
            raise ValueError('Header schedule source identity or primary work changed')
        size=row.get('chunk_markers')
        if (type(size) is not int or size<1 or
                row.get('chunk_count')!=(row['markers']+size-1)//size):
            raise ValueError('Header schedule chunk geometry changed')
        summary=_replay_entries_work(extra.get('entries'),row['samples'],
            row['variant_range'],row['chunk_markers'])
        if any(extra.get(key)!=value for key,value in summary.items()):
            raise ValueError('Header schedule replay entries do not conserve source work')
        joined={}
        for part in (common['source_units'],extra['source_units']):
            for name,interval in part.items():
                if (not isinstance(interval,(list,tuple)) or len(interval)!=2 or
                        any(type(value) is not int or value<0 for value in interval) or interval[0]>interval[1]):
                    raise ValueError('Invalid source unit interval')
                _add(joined,{name:interval})
        if joined!=row['source_units'] or any(row[key]!=common[key]+extra[key] for key in
                ('read_bytes','decode_input_bytes','native_ld_base_update_bytes','native_ld_replay_packed_bytes')):
            raise ValueError('Header schedule does not conserve invariant and replay work')
        if row['ld_replay_count']!=int(common['native_ld_replay_packed_bytes']>0)+extra['replay_count']:
            raise ValueError('Header schedule replay count differs from its work')
    if baseline['common_work']!=candidate['common_work']:
        raise ValueError('Paired schedules must share exact source and invariant work')
    a_entries={entry[0]:entry for entry in baseline['additional_replay_work']['entries']}
    b_entries={entry[0]:entry for entry in candidate['additional_replay_work']['entries']}
    if any(a_entries[at]!=b_entries[at] for at in a_entries.keys()&b_entries.keys()):
        raise ValueError('Shared LD replay start differs between paired schedules')
    a_only=[entry for at,entry in a_entries.items() if at not in b_entries]
    b_only=[entry for at,entry in b_entries.items() if at not in a_entries]
    a_work=_replay_entries_work(a_only,baseline['samples'],baseline['variant_range'],
        baseline['chunk_markers'])
    b_work=_replay_entries_work(b_only,candidate['samples'],candidate['variant_range'],
        candidate['chunk_markers'])
    deltas={};missing=[];cpu=[0.,0.]
    first=a_work['source_units'];second=b_work['source_units']
    for name in sorted(set(first)|set(second)):
        a=first.get(name,[0,0]);b=second.get(name,[0,0])
        interval=[b[0]-a[1],b[1]-a[0]]
        if interval==[0,0]:continue
        deltas[name]=interval
        if name not in cpu_seconds_per_unit:missing.append(name);continue
        rate=positive(name,cpu_seconds_per_unit[name],True)
        cpu[0]+=interval[0]*rate;cpu[1]+=interval[1]*rate
    if not all(math.isfinite(value) for value in cpu):
        raise ValueError('Paired source service overflow')
    # A missing nonnegative primitive rate is harmless for one direction if
    # its entire count-difference interval has the opposite sign. For aligned
    # coarser chunks, candidate-only restarts are absent: missing ULEB costs
    # can only make the candidate cheaper, never invalidate an upper bound.
    lower=cpu[0] if all(deltas[name][0]>=0 for name in missing) else None
    upper=cpu[1] if all(deltas[name][1]<=0 for name in missing) else None
    return dict(candidate_minus_baseline_decoder_cpu_seconds=None if missing else cpu,
        identified_decoder_cpu_seconds=cpu,unpriced_source_units=missing,
        decoder_cpu_delta_lower_bound=lower,decoder_cpu_delta_upper_bound=upper,
        source_unit_delta_intervals=deltas,
        read_bytes_delta=candidate['read_bytes']-baseline['read_bytes'],
        decode_input_bytes_delta=candidate['decode_input_bytes']-baseline['decode_input_bytes'],
        native_ld_base_update_bytes_delta=(candidate['native_ld_base_update_bytes']-
                                           baseline['native_ld_base_update_bytes']),
        native_ld_replay_packed_bytes_delta=(candidate['native_ld_replay_packed_bytes']-
                                             baseline['native_ld_replay_packed_bytes']),
        chunk_count_delta=candidate['chunk_count']-baseline['chunk_count'],
        ld_replay_count_delta=candidate['ld_replay_count']-baseline['ld_replay_count'],
        input_identity=deepcopy(baseline['input_identity']),
        variant_range=list(baseline['variant_range']),
        scope='Conditional paired source-work difference with identical primary record work cancelled. Decoder CPU intervals include only additional LD restarts; storage bytes are exact index extents. No GPU, output, contention, pipeline runtime or switch approval.')


def price_header_work(bounds, cpu_seconds_per_unit):
    """Price intervals with the same independent units as the exact calculator.

    Missing rates make cpu_seconds unavailable, including when only the upper
    endpoint can contain that operation. Partial identified costs are retained.
    Empirical provenance, freshness and live-state admission remain the caller's
    responsibility; this arithmetic function does not certify supplied prices.
    """
    if bounds.get('kind') not in (KIND,SCHEDULE_KIND) or bounds.get('implementation') != 'torch_native_int8':
        raise ValueError('Typed native header work bounds required')
    low = {}; high = {}
    for name, interval in bounds['source_units'].items():
        if (not isinstance(interval, (list, tuple)) or len(interval) != 2
                or any(type(v) is not int or v < 0 for v in interval) or interval[0] > interval[1]):
            raise ValueError('Invalid decoder source-count interval')
        if interval[1]:
            low[name], high[name] = interval
    priced = [price_identified_units(dict(source_units=units, uncounted_mechanisms=[]),
                                    cpu_seconds_per_unit) for units in (low, high)]
    missing = sorted(set(priced[0]['unpriced_source_units'] + priced[1]['unpriced_source_units']))
    interval = [row['identified_cpu_seconds'] for row in priced]
    if not missing and bounds.get('varint_constraints'):
        constraints = bounds['varint_constraints']
        for key in ('varint_integers', 'varint_bytes'):
            row = constraints[key]
            if (not isinstance(row, (list, tuple)) or len(row) != 2
                    or any(type(v) is not int or v < 0 for v in row) or row[0] > row[1]):
                raise ValueError('Invalid variable-integer aggregate interval')
        names = [f'uleb{k}' for k in range(1, 6) if high.get(f'uleb{k}', 0)]
        if not names:
            raise ValueError('Variable-integer constraints require priced integer units')
        costs = [cpu_seconds_per_unit[name] for name in names]
        per_byte = [cpu_seconds_per_unit[name]/int(name[4:]) for name in names]
        var_cost = [sum(row['cpu_seconds_by_unit'][name] for name in names) for row in priced]
        # Keep the same additive unit model, but do not allow the five marginal
        # histogram maxima to invent five different full sets of integers.
        lower = max(var_cost[0], constraints['varint_integers'][0]*min(costs),
                    constraints['varint_bytes'][0]*min(per_byte))
        upper = min(var_cost[1], constraints['varint_integers'][1]*max(costs),
                    constraints['varint_bytes'][1]*max(per_byte))
        interval = [sum(value for name, value in row['cpu_seconds_by_unit'].items() if name not in names)+cost
                    for row, cost in zip(priced, (lower, upper))]
    return dict(cpu_seconds=None if missing else interval, identified_cpu_seconds=interval,
                unpriced_source_units=missing,
                scope='Conditional interval of the existing additive primitive CPU service model, not a measured capacity, pipeline makespan bound, or hardware time guarantee.')


def native_decoder_service_bounds(bounds, profile):
    """Read/decode resource-work intervals using the existing scan equations.

    This component has no setup, GPU, output or graph makespan. Work intervals
    can feed monotone conservation bounds; solving endpoint graphs does not
    generally bound a contended pipeline's completion time.
    """
    priced = price_header_work(bounds, profile['decode_units'])
    if priced['unpriced_source_units']:
        raise ValueError('Unpriced possible decoder units: '+str(priced['unpriced_source_units']))
    q = positive('cpu_fraction', profile['cpu_fraction'])
    if q > 1:
        raise ValueError('CPU fraction must be at most one')
    dram = positive('shared_dram_bytes_per_second', profile['shared_dram_bytes_per_second'])
    storage = positive('read_bytes_per_second', profile['read_bytes_per_second'])
    for key in ('samples', 'markers', 'read_bytes', 'decode_input_bytes',
                'native_ld_base_update_bytes', 'native_ld_replay_packed_bytes'):
        _integer(key, bounds[key], 1 if key in ('samples', 'markers') else 0)
    separate = profile.get('input_read_cpu_prices') is not None
    memory = _native_decoder_memory_bytes(bounds['samples'], bounds['markers'], bounds['read_bytes'],
        bounds['decode_input_bytes'], bounds['native_ld_base_update_bytes'],
        bounds['native_ld_replay_packed_bytes'], separate_buffered_read=separate)
    if separate:
        read = buffered_read_service(bounds['read_bytes'], 1, profile['input_read_cpu_prices'],
            storage_bytes_per_second=storage, dram_bytes_per_second=dram, cpu_fraction=q)
    else:
        read = dict(seconds=bounds['read_bytes']/storage, cpu_seconds=0., resources={'input': storage})
    cpu = priced['cpu_seconds']
    return dict(decode_seconds=[max(value/q, memory/dram) for value in cpu],
                decode_resource_work=dict(cpu=cpu, dram=[memory, memory]), read=read,
                scope='Conditional read/decode component in the existing scan model. Logical traffic and primitive prices remain approximations. Reader initialization, tensor work, output, contention and scheduling are outside this component.')


def native_schedule_source_floor(schedule, profile, shared_capacities):
    """Price an entire unissued PGEN read/decode schedule without chunk graphs.

    The result is a necessary source-stage resource/worker floor under the
    declared component prices and availability scenario. The two endpoints
    bound possible *floors*, not the pipeline completion time. Reader setup,
    GPU, transfer, selection, output and final drain remain outside this term.
    """
    if not isinstance(schedule,dict) or schedule.get('kind')!=SCHEDULE_KIND:
        raise ValueError('Typed whole-header schedule required')
    try:
        source=schedule['input_identity']['path']
    except (KeyError,TypeError):
        raise ValueError('Bound PGEN source identity required') from None
    if input_identity(source)!=schedule['input_identity']:
        raise ValueError('PGEN input changed after schedule inspection')
    # This validates the conserved primary/replay work and the exact schedule
    # geometry without reconstructing per-chunk dictionaries.
    paired_schedule_source_difference(schedule,schedule,profile['decode_units'])
    if not isinstance(shared_capacities,dict) or set(shared_capacities)!={'cpu','dram','input'}:
        raise ValueError('Explicit shared source capacities required')
    limits=dict(cpu='cpu_available_cores',dram='shared_dram_bytes_per_second',
                input='read_bytes_per_second')
    capacities={}
    for name,field in limits.items():
        available=positive('shared '+name,shared_capacities[name])
        if available>positive('profile '+field,profile[field]):
            raise ValueError('Shared '+name+' exceeds the bound source profile')
        capacities[name]=available
    q=positive('cpu_fraction',profile['cpu_fraction'])
    if q>1:raise ValueError('CPU fraction must be at most one')
    depth=_integer('depth',profile['depth'],1)
    decoder_workers=_integer('decode_workers',profile['decode_workers'],1)
    chunks=_integer('chunk_count',schedule['chunk_count'],1)
    workers=min(depth,decoder_workers,chunks)
    priced=price_header_work(schedule,profile['decode_units'])
    if priced['unpriced_source_units']:
        raise ValueError('Unpriced possible decoder units: '+str(priced['unpriced_source_units']))
    decoder_cpu=priced['cpu_seconds']
    for field in ('samples','markers','read_bytes','decode_input_bytes',
                  'native_ld_base_update_bytes','native_ld_replay_packed_bytes'):
        _integer(field,schedule[field],1 if field in ('samples','markers') else 0)
    read_bytes=schedule['read_bytes'];separate=profile.get('input_read_cpu_prices') is not None
    decoder_dram=_native_decoder_memory_bytes(schedule['samples'],schedule['markers'],
        read_bytes,schedule['decode_input_bytes'],schedule['native_ld_base_update_bytes'],
        schedule['native_ld_replay_packed_bytes'],separate_buffered_read=separate)
    if separate:
        read=buffered_read_service(read_bytes,chunks,profile['input_read_cpu_prices'],
            storage_bytes_per_second=capacities['input'],
            dram_bytes_per_second=capacities['dram'],cpu_fraction=q)
        read_cpu=read['cpu_seconds'];read_dram=2*read_bytes
    else:
        read_cpu=read_dram=0.
    cpu=[read_cpu+value for value in decoder_cpu]
    dram=decoder_dram+read_dram
    resource_floor=[max(value/capacities['cpu'],dram/capacities['dram'],
                        read_bytes/capacities['input']) for value in cpu]
    # Sum of individual read services is at least the maximum aggregate read
    # demand, and likewise for decode. At most `workers` such chains overlap.
    read_service=max(read_cpu/q,read_dram/capacities['dram'],
                     read_bytes/capacities['input'])
    worker_floor=[(read_service+max(value/q,decoder_dram/capacities['dram']))/workers
                  for value in decoder_cpu]
    floor=[max(a,b) for a,b in zip(resource_floor,worker_floor)]
    if not all(math.isfinite(value) for value in (*cpu,float(dram),*floor)):
        raise ValueError('Whole-source resource work overflow')
    if input_identity(source)!=schedule['input_identity']:
        raise ValueError('PGEN input changed during source-floor pricing')
    return dict(kind='torchgwas.pgen_schedule_source_floor.v1',
        input_identity=deepcopy(schedule['input_identity']),
        variant_range=list(schedule['variant_range']),samples=schedule['samples'],
        chunk_markers=schedule['chunk_markers'],
        chunk_count=chunks,reader_workers=workers,
        resource_work=dict(cpu_seconds=cpu,dram_bytes=dram,input_bytes=read_bytes,
                           read_cpu_seconds=read_cpu,decoder_cpu_seconds=decoder_cpu,
                           read_dram_bytes=read_dram,decoder_dram_bytes=decoder_dram),
        capacities=capacities,resource_floor_seconds=resource_floor,
        reader_worker_floor_seconds=worker_floor,source_stage_floor_seconds=floor,
        prediction_complete=False,selection_validated=False,
        scope='Header-conditional necessary read/decode source-stage floors for the complete fixed chunk schedule. Endpoint floors are not elapsed-time upper/lower bounds for the GPU/output pipeline. Reader initialization, other stages, contention beyond declared capacities and final drain are omitted.')
