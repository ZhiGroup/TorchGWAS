"""Bind bounded model horizons to the exact held productive source frontier.

This adapter changes chunk size only. It does not supply capacity measurements,
live admission, future-source assumptions or in-flight boundary adjustments.
"""
from copy import deepcopy
import math

from .adaptive_chunks import aligned_chunk_sizes
from .calibration_cache import _digest
from .window_forecast import forecast_remaining_windows


def productive_window_proposal(snapshot, comparisons, *, chunk_sizes,
        boundary_adjustments, relative_model_error, max_slope_change,
        max_extrapolation, assumptions, max_partitions=8, max_issued_chunks=128):
    """Return one chunk proposal and its executor/forecast audit.

    Every unfinished partition must appear once in each layout, at its actual
    unissued cursor and native reader ordinal. Fixed device/phenotype ownership
    is preserved. Already-issued work remains in explicit boundary adjustments,
    never silently added to the unissued pair total or declared complete.
    """
    sizes=aligned_chunk_sizes(chunk_sizes)
    for value in (max_partitions,max_issued_chunks):
        if type(value) is not int or value<1:raise ValueError('Positive productive forecast limits required')
    if (not isinstance(snapshot,dict) or snapshot.get('prefix_complete') is not True
            or snapshot.get('finished') is not None or snapshot.get('first_written') is None):
        raise ValueError('Active complete source prefix after useful output required')
    written=snapshot['first_written']
    if isinstance(written,bool) or not isinstance(written,(int,float)) or not math.isfinite(written):
        raise ValueError('Finite first-output timestamp required')
    revision=snapshot.get('issued_revision');current=snapshot.get('current_chunk_size')
    if type(revision) is not int or revision<1 or type(current) is not int or current not in sizes:
        raise ValueError('Bound revision and admitted current chunk size required')
    partitions=snapshot.get('partitions')
    if not isinstance(partitions,list) or not 0<len(partitions)<=max_partitions:
        raise ValueError('Bounded complete partition list required')
    unfinished={};total_issued=0;ids=set();rects=[];remaining=[]
    for row in partitions:
        if not isinstance(row,dict):raise ValueError('Explicit executor partition required')
        key=row.get('id');device=row.get('device')
        if not isinstance(key,str) or not key or key in ids:raise ValueError('Unique executor partition ids required')
        ids.add(key)
        if not isinstance(device,str) or not device.startswith('cuda:') or not device[5:].isdigit():
            raise ValueError('Explicit executor CUDA device required')
        for name in ('variant_range','trait_range'):
            span=row.get(name)
            if (not isinstance(span,(list,tuple)) or len(span)!=2 or
                    any(type(v) is not int for v in span) or not 0<=span[0]<span[1]):
                raise ValueError('Nonempty executor ranges required')
        lo,hi=row['variant_range'];a,b=row['trait_range']
        if any(max(lo,x)<min(hi,y) and max(a,c)<min(b,d) for x,y,c,d in rects):
            raise ValueError('Executor association partitions overlap')
        rects.append((lo,hi,a,b))
        count=row.get('issued_chunks');ranges=row.get('ranges');cursor=lo
        if type(count) is not int or count<0 or not isinstance(ranges,list) or len(ranges)!=count:
            raise ValueError('Complete issued-range history required')
        total_issued+=count
        if total_issued>max_issued_chunks:raise ValueError('Issued history exceeds forecast budget')
        for span in ranges:
            if (not isinstance(span,(list,tuple)) or len(span)!=2 or any(type(v) is not int for v in span)
                    or span[0]!=cursor or not cursor<span[1]<=hi or span[1]-cursor>sizes[-1]):
                raise ValueError('Invalid contiguous admitted source prefix')
            if span[1]-cursor not in sizes and (span[1]!=hi or span[1]-cursor>=sizes[0]):
                raise ValueError('Issued chunk width was not admitted')
            cursor=span[1]
        if row.get('cursor')!=cursor:raise ValueError('Executor cursor differs from its issued prefix')
        if cursor<hi:
            identity=(device,a,b,cursor)
            if identity in unfinished:raise ValueError('Ambiguous unfinished partition')
            unfinished[identity]=row
            remaining.append(dict(id=key,device=device,variant_range=[cursor,hi],trait_range=[a,b],issued_chunks=count))
    if total_issued!=revision:raise ValueError('Executor revision differs from its issued count')
    if not unfinished:raise ValueError('No unissued associations remain')
    if not isinstance(comparisons,(list,tuple)) or len(comparisons)!=3:
        raise ValueError('Exactly three bounded comparisons required')
    contract=comparisons[0].get('comparison_contract')
    if not isinstance(contract,dict):raise ValueError('Bound analytical comparison required')
    span_vectors={};proposed=None;fixed_bindings={}
    for name in ('baseline','candidate'):
        configuration=contract[name];windows=configuration['windows'];seen=set();spans=[]
        if len(windows)!=len(unfinished):raise ValueError('Every unfinished partition must be modeled')
        for w in windows:
            identity=(w['device'],*w['trait_range'],w['variant_start'])
            if identity not in unfinished or identity in seen:
                raise ValueError('Model window differs from held unissued frontier or fixed ownership')
            seen.add(identity);row=unfinished[identity]
            if w['issued_chunks']!=row['issued_chunks']:
                raise ValueError('Model reader ordinal differs from actual issued chunks')
            binding=(w.get('fixed_data_sha256'),w.get('chunk_invariant_profile_sha256'))
            if any(not isinstance(value,str) or not value for value in binding):
                raise ValueError('Fixed data and chunk-invariant profile bindings required')
            if name=='baseline':fixed_bindings[identity]=binding
            elif binding!=fixed_bindings[identity]:
                raise ValueError('Chunk adaptation cannot change fixed data, prices or execution settings')
            size=w['chunk_markers']
            if type(size) is not int or size not in sizes:raise ValueError('Model chunk size was not admitted')
            if name=='baseline' and size!=current:raise ValueError('Baseline must model the actual current chunk size')
            if name=='candidate':
                if proposed is not None and proposed!=size:raise ValueError('One shared candidate chunk size required')
                proposed=size
            spans.append(row['variant_range'][1]-row['cursor'])
        span_vectors[name]=spans
    # This bridge applies one size on fixed partitions; changing a device, tile,
    # decoder allocation, math mode or prices requires a separate admitted move.
    if contract['baseline']['partition_axis']!=contract['candidate']['partition_axis']:
        raise ValueError('Chunk adaptation cannot change partition axis')
    if contract['reduction']=='jagwas' and (contract['baseline']['partition_axis']!='variant' or
            any(tuple(row['trait_range'])!=(0,contract['total_traits']) for row in partitions)):
        raise ValueError('JAGWAS must retain the complete phenotype panel')
    pairs=sum((r['variant_range'][1]-r['variant_range'][0])*(r['trait_range'][1]-r['trait_range'][0]) for r in remaining)
    forecast=forecast_remaining_windows(comparisons,remaining_pairs=pairs,
        remaining_partition_spans=span_vectors,boundary_adjustments=boundary_adjustments,
        relative_model_error=relative_model_error,max_slope_change=max_slope_change,
        max_extrapolation=max_extrapolation,assumptions=assumptions)
    before=forecast['forecasts']['baseline']['lower_seconds'];after=forecast['forecasts']['candidate']['upper_seconds']
    eligible=forecast['status']=='stable_scenario'
    proposal=dict(chunk_size=proposed if eligible else current,
        baseline_seconds=before if eligible else 0.,candidate_seconds=after if eligible else 0.)
    audit=dict(issued_revision=revision,current_size=current,proposed_size=proposed,
        remaining_partitions=remaining,remaining_pairs=pairs,forecast_status=forecast['status'],
        comparison_contract_sha256=_digest(contract),forecasts=deepcopy(forecast['forecasts']),
        partition_extent_scenarios=deepcopy(forecast['partition_extent_scenarios']),
        scenario_gain_floor_seconds=forecast['scenario_gain_floor_seconds'],
        assumptions=deepcopy(assumptions),unpriced_terms=forecast['unpriced_terms'],
        writer_regime_warnings=forecast['writer_regime_warnings'],selection_validated=False,
        scope='Exact held unissued extent with balanced-core/envelope model scenarios for unequal partition progress. Already-issued GPU/queue/writer work remains caller-supplied boundary cost. No hardware-time bound, independent-price qualification or layout reassignment.')
    return proposal,audit
