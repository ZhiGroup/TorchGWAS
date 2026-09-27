"""Bounded early-job checks and refresh of cached component observations.

Comparison is a declared drift heuristic on matching work, not a statistical
confidence interval or conversion of loaded event spans to hardware capacity.
The existing analytical planner remains responsible for candidate scoring and
memory admission. No whole-GWAS timing table is fitted here.
"""
from __future__ import annotations

import copy
import math
from numbers import Integral
from statistics import median
import threading
import time

from .adaptive_chunks import InitialChunkMeasurements, _positive_size


CACHE_NAME = 'initial_chunk_components.v3'


def _finite_nonnegative(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('Finite nonnegative '+name+' required')
    return float(value)


def _row(row):
    if not isinstance(row, dict):
        raise ValueError('Chunk observation must be a dictionary')
    for key in ('start', 'end', 'capacity', 'result_bytes', 'result_blocks'):
        value=row[key]
        if isinstance(value,bool) or not isinstance(value,Integral) or value<0:
            raise ValueError('Invalid observation '+key)
    if row['end']<=row['start'] or row['end']-row['start']>row['capacity'] or not row['result_blocks']:
        raise ValueError('Invalid observation range/capacity')
    device=row['device']
    if not isinstance(device,str) or not device.startswith('cuda:') or not device[5:].isdigit():
        raise ValueError('Explicit observation CUDA device required')
    started=_finite_nonnegative(row['read_started'],'read start')
    finished=_finite_nonnegative(row['read_finished'],'read finish')
    if finished<started:
        raise ValueError('Reader timestamp order is invalid')
    reader_cpu=row.get('read_cpu_seconds')
    consumer_cpu=row.get('consumer_cpu_seconds')
    values={'read_decode':finished-started,
            'read_thread_cpu':None if reader_cpu is None else
                _finite_nonnegative(reader_cpu,'reader thread CPU'),
            'consumer_thread_cpu':None if consumer_cpu is None else
                _finite_nonnegative(consumer_cpu,'consumer thread CPU')}
    for name in ('read_runnable_wait_seconds','read_probe_wall_seconds',
                 'consumer_runnable_wait_seconds','consumer_probe_wall_seconds'):
        if row.get(name) is not None:
            _finite_nonnegative(row[name],name)
    cuda=row.get('cuda')
    names=('h2d','conversion','statistics_and_reduction','result_transfer')
    if cuda is not None:
        if not isinstance(cuda,dict) or set(cuda)!=set(names):
            raise ValueError('Complete CUDA interval declaration required')
        for name in names:
            values[name]=None if cuda[name] is None else _finite_nonnegative(cuda[name],name)
    else:
        values.update(dict.fromkeys(names))
    key=(device,int(row['start']),int(row['end']),int(row['capacity']))
    return key,values


def compare_component_windows(previous, fresh, *, minimum_samples=2, ratio_limit=2., absolute_floor_seconds=2e-6):
    """Compare median spans only for identical device/range/capacity tuples.

    Dependencies (source, data, phenotype/output, library and execution context)
    must already match. Reader spans include I/O and decode; CUDA spans include
    stream scheduling gaps. Consumer CPU is downstream iterator work; its
    output-mode-dependent boundary is not selector-only service or hardware
    availability. Consumer suspension is not interpreted as durable writing.
    """
    minimum_samples=_positive_size(minimum_samples,'minimum_samples')
    if isinstance(ratio_limit,bool) or not isinstance(ratio_limit,(int,float)) or not math.isfinite(ratio_limit) or ratio_limit<=1:
        raise ValueError('Finite ratio_limit greater than one required')
    floor=_finite_nonnegative(absolute_floor_seconds,'absolute floor')
    baseline={}
    for row in previous:
        key,values=_row(row)
        if key in baseline: raise ValueError('Duplicate cached source range')
        baseline[key]=(row,values)
    pairs=[];missing=[];seen=set();changed_output=[]
    for row in fresh:
        key,values=_row(row)
        if key in seen: raise ValueError('Duplicate fresh source range')
        seen.add(key)
        if key not in baseline:
            missing.append(list(key));continue
        old,old_values=baseline[key]
        if (row['result_bytes'],row['result_blocks'])!=(old['result_bytes'],old['result_blocks']):
            changed_output.append(list(key))
        pairs.append((old_values,values))
    result=dict(matched_samples=len(pairs),minimum_samples=minimum_samples,
        missing_ranges=missing,changed_output_ranges=changed_output,metrics={},
        ratio_limit=float(ratio_limit),absolute_floor_seconds=floor,
        scope='Heuristic change detection for measured intervals only; not a hardware-capacity or full-pipeline validity certificate.')
    if missing or len(pairs)<minimum_samples:
        return dict(result,status='incomparable')
    changed=bool(changed_output);unavailable=[]
    for name in ('read_decode','read_thread_cpu','consumer_thread_cpu','h2d','conversion',
                 'statistics_and_reduction','result_transfer'):
        old=[a[name] for a,b in pairs];new=[b[name] for a,b in pairs]
        if all(v is None for v in old+new):
            unavailable.append(name);continue
        if any(v is None for v in old+new):
            return dict(result,status='incomparable',incomparable_metric=name)
        before,after=median(old),median(new)
        # Absolute tolerance prevents tiny near-resolution spans from giving
        # enormous ratios. None represents a zero denominator, not infinity.
        ratio=None if before==0 else after/before
        drift=abs(after-before)>floor and (after>ratio_limit*before or before>ratio_limit*after)
        changed |= drift
        result['metrics'][name]=dict(cached_median_seconds=before,fresh_median_seconds=after,
                                     ratio=ratio,drift=drift)
    return dict(result,status='drift' if changed else 'consistent',unobserved_metrics=unavailable)


class InitialCalibrationController(InitialChunkMeasurements):
    """Check a short per-GPU sample, expand on drift, then stop within budget.

    Pass as _chunk_observer to a native scan. Call finish(successful=True) after
    the complete scan/consumer succeeds. Each GPU has its own immutable record
    and observation timestamp, so validating one does not renew another's age.
    Drift detection is not automatic candidate selection; it supplies evidence
    for the shared calculator and never changes an admitted chunk by itself.
    """
    def __init__(self, devices, *, cache, dependencies, provenance,
                 validation_chunks_per_device=2, max_chunks_per_device=8,
                 max_age_seconds=300., ratio_limit=2., absolute_floor_seconds=2e-6,
                 warmup_chunks=1, stride=4, max_window_seconds=10., cuda_events=True,
                 max_decision_cpu_seconds=.05):
        super().__init__(devices,max_chunks_per_device=max_chunks_per_device,
            warmup_chunks=warmup_chunks,stride=stride,max_window_seconds=max_window_seconds,cuda_events=cuda_events)
        self.validation_chunks=_positive_size(validation_chunks_per_device,'validation_chunks_per_device')
        if self.validation_chunks>self.max_chunks:
            raise ValueError('Validation sample cannot exceed the total sample budget')
        # Validate comparison settings even if this run has no prior record.
        compare_component_windows([],[],minimum_samples=self.validation_chunks,
            ratio_limit=ratio_limit,absolute_floor_seconds=absolute_floor_seconds)
        if not isinstance(provenance,dict) or not provenance:
            raise ValueError('Explicit calibration job provenance required')
        self.cache=cache
        self.dependencies=copy.deepcopy(dependencies)
        self.provenance=copy.deepcopy(provenance)
        self.max_age_seconds=max_age_seconds
        self.ratio_limit=ratio_limit
        self.absolute_floor_seconds=absolute_floor_seconds
        self.max_decision_cpu_seconds=_finite_nonnegative(max_decision_cpu_seconds,'decision CPU budget')
        if self.max_decision_cpu_seconds==0:
            raise ValueError('Positive decision CPU budget required')
        self._decision_cpu_seconds=0.
        self._decision_budget_exhausted=False
        self._controller_lock=threading.Lock()
        self._previous={};self._decisions={};self._finished=None
        for device in self.devices:
            previous=cache.lookup('stage_observations',CACHE_NAME+':'+device,
                                  dependencies=self.dependencies,max_age_seconds=max_age_seconds)
            if previous['hit']:
                try:
                    value=previous['record']['value'];rows=value['observations']
                    if ('observed_unix_seconds' not in previous['record'] or value['device']!=device
                            or value.get('consumer_probe_protocol')!='thread_cpu_schedstat_probe_wall_v2'
                            or not isinstance(rows,list) or len(rows)<self.validation_chunks
                            or any(_row(row)[0][0]!=device for row in rows)):
                        raise ValueError('Incomplete or unbound cached measurements')
                    compare_component_windows(rows,rows,minimum_samples=self.validation_chunks,
                        ratio_limit=ratio_limit,absolute_floor_seconds=absolute_floor_seconds)
                except (KeyError,TypeError,ValueError):
                    previous=dict(hit=False,reason='invalid_component_window')
            self._previous[device]=previous
            self._limits[device]=self.validation_chunks if previous['hit'] else self.max_chunks
            self._decisions[device]=dict(state='validating' if previous['hit'] else 'collecting',
                cache_hit=previous['hit'],cache_reason=previous.get('reason'),
                previous_record_sha256=previous.get('record_sha256'),comparison=None)

    def __call__(self, observation):
        started=time.thread_time()
        try:
            self._observe(observation)
        finally:
            with self._controller_lock:
                self._decision_cpu_seconds+=time.thread_time()-started
                if self._decision_cpu_seconds>=self.max_decision_cpu_seconds:
                    self._decision_budget_exhausted=True
                    self.stop()

    def _observe(self, observation):
        super().__call__(observation)
        with self._controller_lock:
            snapshot=super().snapshot()
            rows=[row for row in snapshot['observations'] if row['device']==observation.device]
            decision=self._decisions[observation.device]
            if decision['state']=='validating' and len(rows)>=self.validation_chunks:
                comparison=compare_component_windows(self._previous[observation.device]['record']['value']['observations'],rows,
                    minimum_samples=self.validation_chunks,ratio_limit=self.ratio_limit,
                    absolute_floor_seconds=self.absolute_floor_seconds)
                decision['comparison']=comparison
                if comparison['status']=='consistent':
                    decision['state']='cached_consistent'
                else:
                    decision['state']='refreshing'
                    with self._lock:
                        self._limits[observation.device]=self.max_chunks
            if decision['state'] in ('collecting','refreshing') and len(rows)>=self.max_chunks:
                decision['state']='collected' if decision['state']=='collecting' else 'refreshed'

    def snapshot(self):
        with self._controller_lock:
            return dict(super().snapshot(),refresh_decisions=copy.deepcopy(self._decisions),
                        validation_chunks_per_device=self.validation_chunks,
                        decision_cpu_seconds=self._decision_cpu_seconds,
                        max_decision_cpu_seconds=self.max_decision_cpu_seconds,
                        decision_budget_exhausted=self._decision_budget_exhausted,
                        publication=copy.deepcopy(self._finished))

    def finish(self, *, successful):
        if type(successful) is not bool:
            raise ValueError('successful must be boolean')
        self.stop()
        with self._controller_lock:
            if self._finished is not None:
                return copy.deepcopy(self._finished)
            snapshot=super().snapshot();publications={};incomplete=[]
            for device in self.devices:
                decision=self._decisions[device]
                rows=[row for row in snapshot['observations'] if row['device']==device]
                pending=any(row['device']==device for row in snapshot['pending'])
                if (not successful or pending or decision['state'] not in
                        ('collected','refreshed','cached_consistent')):
                    incomplete.append(device);continue
                # Use the oldest observation in the saved window. Delayed
                # publication (e.g. an hours-long job) cannot reset freshness.
                observed=snapshot['started_unix_seconds']+min(
                    row['read_started']-snapshot['started_perf_counter'] for row in rows)
                name=CACHE_NAME+':'+device
                if decision['state']=='cached_consistent': name+=':validation'
                value=dict(device=device,observations=rows,decision=copy.deepcopy(decision),
                    consumer_probe_protocol='thread_cpu_schedstat_probe_wall_v2',
                    scope='Matching source-chunk component intervals; not hardware capacity or candidate admission.')
                provenance=dict(self.provenance,previous_record_sha256=decision['previous_record_sha256'],
                    refresh_state=decision['state'])
                try:
                    publications[device]=self.cache.store('stage_observations',name,value,
                        dependencies=self.dependencies,provenance=provenance,
                        max_age_seconds=self.max_age_seconds,observed_unix_seconds=observed)
                except OSError as error:
                    # Optional cache I/O must not discard already-computed GWAS output.
                    publications[device]=dict(saved=False,error=str(error))
            self._finished=dict(successful=successful,publications=publications,incomplete_devices=incomplete,
                decision_cpu_seconds=self._decision_cpu_seconds,
                max_decision_cpu_seconds=self.max_decision_cpu_seconds,
                decision_budget_exhausted=self._decision_budget_exhausted,
                decisions=copy.deepcopy(self._decisions),
                scope='Evidence refresh only; candidate scoring, admission and selection remain with the calculator.')
            return copy.deepcopy(self._finished)
