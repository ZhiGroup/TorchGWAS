"""Bounded reusable calculation work and cooperative incremental planning cost.

Memoization is not calibration evidence or a freshness decision. Callers still
validate immutable parameter records and live admission before using results.
"""
from collections import OrderedDict
from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
import json
import math
import sys
import threading
import time


_ACTIVE=ContextVar('torchgwas_planning_work',default=None)


def _positive_integer(value,name):
    if type(value) is not int or value<1:raise ValueError('Positive integer '+name+' required')
    return value


def _seconds(value,name,zero=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (not zero and value==0):
        raise ValueError('Invalid '+name)
    return float(value)


def _size(value,seen=None):
    # Python-ledger estimate, not a process RSS or allocator-reservation bound.
    seen=set() if seen is None else seen
    if id(value) in seen:return 0
    seen.add(id(value));amount=sys.getsizeof(value)
    if isinstance(value,dict):amount+=sum(_size(k,seen)+_size(v,seen) for k,v in value.items())
    elif isinstance(value,(tuple,list)):amount+=sum(_size(v,seen) for v in value)
    return amount


class PlanningWorkCache:
    """An explicit session may span steps without leaking into another job.

    Both entry count and retained-ledger estimate are bounded. Each consumer
    gets a copy, and changes to binding values or implementation identities
    cause misses. No fitted association durations are stored here.
    """
    def __init__(self,*,max_entries=128,max_bytes=32<<20):
        self.max_entries=_positive_integer(max_entries,'max_entries')
        self.max_bytes=_positive_integer(max_bytes,'max_bytes')
        self._entries=OrderedDict();self._bytes=0;self._lock=threading.RLock()
        self._hits=self._misses=self._bypasses=self._evictions=0
        self._closed=False

    @contextmanager
    def activate(self):
        with self._lock:
            if self._closed:raise ValueError('Planning cache is closed')
        token=_ACTIVE.set(self)
        try:yield self
        finally:_ACTIVE.reset(token)

    def get(self,namespace,bindings,build,implementation):
        with self._lock:
            if self._closed:raise ValueError('Planning cache is closed')
        try:
            encoded=json.dumps(bindings,sort_keys=True,separators=(',',':'),allow_nan=False)
            # Never narrow calculator inputs or merge tuples/non-string keys
            # with superficially similar JSON bindings.
            if json.loads(encoded)!=bindings:raise TypeError('Non-JSON binding')
            key=(namespace,tuple(implementation),encoded);hash(key)
        except (ValueError,TypeError):
            with self._lock:self._bypasses+=1
            return build()
        with self._lock:
            if self._closed:raise ValueError('Planning cache is closed')
            if key in self._entries:
                self._hits+=1;self._entries.move_to_end(key)
                return deepcopy(self._entries[key][0])
            self._misses+=1
        value=build()
        retained=deepcopy(value);amount=_size(key)+_size(retained)
        with self._lock:
            if self._closed:return value
            if amount>self.max_bytes:self._bypasses+=1;return value
            # Concurrent identical misses may finish out of order. Their
            # deterministic ledgers share one entry, never two byte charges.
            if key in self._entries:return value
            while self._entries and (len(self._entries)>=self.max_entries or self._bytes+amount>self.max_bytes):
                _,(_,removed)=self._entries.popitem(last=False);self._bytes-=removed;self._evictions+=1
            self._entries[key]=(retained,amount);self._bytes+=amount
        return value

    def close(self):
        with self._lock:self._closed=True;self._entries.clear();self._bytes=0

    def snapshot(self):
        with self._lock:
            return dict(hits=self._hits,misses=self._misses,bypasses=self._bypasses,evictions=self._evictions,
                entries=len(self._entries),estimated_retained_bytes=self._bytes,
                max_entries=self.max_entries,max_bytes=self.max_bytes,closed=self._closed)


def cached_planning_work(namespace,bindings,build,*,implementation):
    cache=_ACTIVE.get()
    return build() if cache is None else cache.get(namespace,bindings,build,implementation)


@contextmanager
def planning_work_scope():
    """A whole-plan call reuses an enclosing incremental session when present."""
    active=_ACTIVE.get()
    if active is not None:
        yield active
        return
    cache=PlanningWorkCache()
    try:
        with cache.activate():yield cache
    finally:cache.close()


class IncrementalPlanningBudget:
    """Gate individual productive-run planning steps before invoking them.

    Horizons and step costs are explicit forecasts, not statistical bounds.
    An admitted step is cooperative: Python/native calls cannot be preempted
    here. Actual overruns are charged and their late results cannot authorize
    a switch. This class does not claim a Bayesian information-value policy.
    """
    def __init__(self,*,max_steps=4,max_cpu_seconds=.05,max_window_seconds=10.):
        self.max_steps=_positive_integer(max_steps,'max_steps')
        self.max_cpu_seconds=_seconds(max_cpu_seconds,'max_cpu_seconds')
        self.max_window_seconds=_seconds(max_window_seconds,'max_window_seconds')
        self._lock=threading.Lock();self._started=None;self._closed=False;self._running=False
        self._cpu=0.;self._wall=0.;self._rows=[];self._last_skip=None

    def start_after_first_output(self):
        """Called on useful output delivery; repeat calls cannot renew the window."""
        with self._lock:
            if self._closed:raise ValueError('Planning budget is closed')
            if self._started is None:self._started=time.perf_counter()
            return self._started

    def run_step(self,build,*,remaining_seconds,expected_cpu_seconds,expected_wall_seconds,
                 switching_seconds=0.,expected_gain_seconds=None,publication_seconds=0.,reserve_seconds=0.):
        remaining=_seconds(remaining_seconds,'remaining_seconds',True)
        cpu=_seconds(expected_cpu_seconds,'expected_cpu_seconds')
        wall=_seconds(expected_wall_seconds,'expected_wall_seconds')
        switch=_seconds(switching_seconds,'switching_seconds',True)
        publication=_seconds(publication_seconds,'publication_seconds',True)
        reserve=_seconds(reserve_seconds,'reserve_seconds',True)
        try:deferred=math.fsum((switch,publication,reserve))
        except OverflowError:raise ValueError('Tuning cost overflow') from None
        if not math.isfinite(deferred):raise ValueError('Tuning cost overflow')
        gain=expected_gain_seconds
        if gain is not None:
            if isinstance(gain,bool) or not isinstance(gain,(int,float)) or not math.isfinite(gain) or gain>remaining:
                raise ValueError('Finite forecast gain no greater than remaining time required')
            gain=float(gain)
        if not callable(build):raise ValueError('Callable single planning step required')
        with self._lock:
            now=time.perf_counter();reason=None
            if self._closed:reason='closed'
            elif self._started is None:reason='no_useful_output_yet'
            elif self._running:reason='step_in_flight'
            elif len(self._rows)>=self.max_steps:reason='step_budget'
            elif self._cpu+cpu>self.max_cpu_seconds:reason='cpu_budget'
            elif now-self._started+wall+switch>self.max_window_seconds:reason='window_budget'
            elif remaining<=wall+deferred:reason='insufficient_remaining_horizon'
            elif gain is not None and gain<=self._wall+wall+deferred:reason='forecast_gain_does_not_repay_step'
            if reason:
                self._last_skip=reason
                return dict(evaluated=False,reason=reason,value=None)
            self._running=True;began=now;cpu_start=time.thread_time()
        value=None;error=None
        try:
            value=build()
        except BaseException as caught:
            error=type(caught).__name__
            raise
        finally:
            elapsed=time.perf_counter()-began;used=time.thread_time()-cpu_start
            with self._lock:
                self._cpu+=used;self._wall+=elapsed;self._running=False
                on_time=(not self._closed and self._cpu<=self.max_cpu_seconds
                         and time.perf_counter()-self._started+switch<=self.max_window_seconds
                         and elapsed+deferred<remaining)
                cost=math.fsum((self._wall,deferred))
                row=dict(cpu_seconds=used,wall_seconds=elapsed,expected_cpu_seconds=cpu,
                    expected_wall_seconds=wall,remaining_seconds=remaining,switching_seconds=switch,
                    publication_seconds=publication,reserve_seconds=reserve,
                    cumulative_planning_wall_seconds=self._wall,total_tuning_cost_seconds=cost,
                    expected_gain_seconds=gain,error=error,
                    usable_for_decision=on_time and error is None and (gain is None or gain>cost))
                self._rows.append(row)
        return dict(evaluated=True,value=value,**row)

    def finish(self):
        with self._lock:self._closed=True
        return self.snapshot()

    def snapshot(self):
        with self._lock:
            return dict(started=self._started is not None,started_perf_counter=self._started,
                closed=self._closed,step_in_flight=self._running,cpu_seconds=self._cpu,
                wall_seconds=self._wall,max_steps=self.max_steps,max_cpu_seconds=self.max_cpu_seconds,
                max_window_seconds=self.max_window_seconds,last_skip=self._last_skip,steps=deepcopy(self._rows),
                scope='Cooperative planning-thread CPU and cumulative step wall costs after useful output. Each decision must repay all planning steps plus its declared switch, cache publication and reserve costs. Forecast cost gate; not a hard deadline, posterior utility or end-to-end JIT overhead.')
