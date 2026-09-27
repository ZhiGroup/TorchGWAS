"""Optional cross-job reuse of duration-free tensor operation ledgers.

Only meta-operation counts are stored as structural calibration. Prices,
device geometry, loaded timings and live availability do not enter this store.
Missing traces are staged in bounded memory; publication is explicit so fsync
does not interrupt a productive tuning step.
"""
from collections import OrderedDict
from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import platform
import threading
from types import CodeType,FunctionType

from .calibration_cache import CalibrationParameterCache, _digest, _json
from .planning_session import _positive_integer, _size


_ACTIVE=ContextVar('torchgwas_structural_tensor_work',default=None)
_FILES=('structural_tensor_cache.py','tensor_work.py','linear.py','native_scan.py',
        'reduction_tensor_work.py','reduce.py','jagwas_projection.py','jagwas_blocks.py','mechanistic_plan.py')
_NAMES=frozenset({'statistics.v1','jagwas.v3'})


def _sources():
    root=Path(__file__).parent
    return {name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in _FILES}


def _runtime():
    import torch
    return dict(python=platform.python_version(),python_implementation=platform.python_implementation(),
        torch=str(torch.__version__),torch_git=getattr(torch.version,'git_version',None),
        torch_cuda=torch.version.cuda,default_dtype=str(torch.get_default_dtype()))


def _code_identity(code):
    # marshal output can change as Python materializes/references nested code
    # objects. Fingerprint explicit immutable semantics, never refcount state.
    def constant(value):
        if isinstance(value,CodeType):return dict(code=_code_identity(value))
        if isinstance(value,tuple):return dict(tuple=[constant(v) for v in value])
        if isinstance(value,frozenset):return dict(frozenset=sorted((constant(v) for v in value),key=_json))
        if type(value) in (type(None),bool,int,float,complex,str,bytes) or value is Ellipsis:
            return dict(type=type(value).__name__,value=repr(value))
        raise TypeError('Undeclared code constant')
    return dict(bytecode=code.co_code.hex(),constants=[constant(v) for v in code.co_consts],
        names=list(code.co_names),variables=list(code.co_varnames),free=list(code.co_freevars),cells=list(code.co_cellvars),
        args=code.co_argcount,posonly=code.co_posonlyargcount,kwonly=code.co_kwonlyargcount,
        flags=code.co_flags,exceptions=getattr(code,'co_exceptiontable',b'').hex())


def _dependencies(sources,bindings,implementation):
    functions=[]
    for function in implementation:
        # Dynamic wrappers/closures have undeclared state. Keep their ordinary
        # behavior but never load or publish a supposedly permanent ledger.
        if not isinstance(function,FunctionType) or function.__closure__:
            raise TypeError('Dynamic trace implementation')
        functions.append(dict(module=function.__module__,name=function.__qualname__,
            code_sha256=_digest(_code_identity(function.__code__)),
            defaults=function.__defaults__,kwdefaults=function.__kwdefaults__))
    typed={name:dict(type=type(value).__name__,value=value) for name,value in bindings.items()}
    # Restrict the portable request to the wrappers' scalar shape/mode inputs.
    if any(type(value) not in (int,bool,str) for value in bindings.values()):
        raise TypeError('Nonportable tensor shape')
    return json.loads(_json(dict(protocol='torchgwas.structural_tensor_work.v1',
        source_sha256=sources,runtime=_runtime(),shape=typed,implementation=functions)))


class StructuralTensorWorkCache:
    """A scoped, bounded set of reusable source ledgers backed by immutable JSON.

    Construct one for a job, activate it around bounded calculations, and call
    publish(successful=True) after successful work. Loading never republishes
    or refreshes an artifact. Source/library/shape/dtype/implementation changes
    cause misses; this class has no service-rate or memory-admission authority.
    """
    def __init__(self,directory,*,max_entries=16,max_bytes=8<<20):
        self.max_entries=_positive_integer(max_entries,'max_entries')
        self.max_bytes=_positive_integer(max_bytes,'max_bytes')
        self.cache=CalibrationParameterCache(directory)
        self.sources=_sources();self._entries=OrderedDict();self._bytes=0
        self._lock=threading.RLock();self._closed=False
        self._counts=dict(memory_hits=0,disk_hits=0,misses=0,bypasses=0,evictions=0,cache_errors=0)
        self._publication=None

    @contextmanager
    def activate(self):
        with self._lock:
            if self._closed:raise ValueError('Structural tensor cache is closed')
        token=_ACTIVE.set(self)
        try:yield self
        finally:_ACTIVE.reset(token)

    def _get(self,name,bindings,build,implementation):
        if name not in _NAMES:raise ValueError('Unknown structural tensor ledger')
        with self._lock:
            if self._closed:raise ValueError('Structural tensor cache is closed')
        try:
            dependencies=_dependencies(self.sources,bindings,implementation)
            key=_digest(dict(name=name,dependencies=dependencies))
        except (ValueError,TypeError):
            with self._lock:self._counts['bypasses']+=1
            return build()
        with self._lock:
            if key in self._entries:
                self._counts['memory_hits']+=1;self._entries.move_to_end(key)
                return json.loads(self._entries[key]['payload'])
        try:
            hit=self.cache.lookup('source_work','tensor_trace.'+name,dependencies=dependencies)
        except OSError:
            with self._lock:self._counts['cache_errors']+=1
            hit=dict(hit=False)
        if hit['hit']:
            value=hit['record']['value'];pending=False
            with self._lock:self._counts['disk_hits']+=1
        else:
            with self._lock:self._counts['misses']+=1
            value=build();pending=True
        # Keep the immutable serialization, not thousands of Python containers.
        # JSON decoding gives each consumer its own ledger, while byte sizing
        # avoids recursively walking and deep-copying the entire trace again.
        try:
            payload=_json(value)
            if json.loads(payload)!=value:raise TypeError('Nonportable tensor ledger')
        except (ValueError,TypeError):
            with self._lock:self._counts['bypasses']+=1
            return value
        entry=dict(name=name,dependencies=dependencies,payload=payload,pending=pending)
        amount=_size(key)+_size(entry)
        with self._lock:
            if self._closed:return value
            if amount>self.max_bytes:
                self._counts['bypasses']+=1;return value
            if key not in self._entries:
                while self._entries and (len(self._entries)>=self.max_entries or self._bytes+amount>self.max_bytes):
                    _,old=self._entries.popitem(last=False);self._bytes-=old['bytes'];self._counts['evictions']+=1
                entry['bytes']=amount;self._entries[key]=entry;self._bytes+=amount
        return value

    def publish(self,*,successful):
        if type(successful) is not bool:raise ValueError('Boolean successful required')
        with self._lock:
            if self._closed:raise ValueError('Structural tensor cache is closed')
            if not successful:
                self._publication=dict(status='unsuccessful',stored=[]);return deepcopy(self._publication)
            if _sources()!=self.sources:
                self._publication=dict(status='source_changed',stored=[]);return deepcopy(self._publication)
            stored=[]
            for entry in self._entries.values():
                if not entry['pending']:continue
                try:
                    record=self.cache.store('source_work','tensor_trace.'+entry['name'],json.loads(entry['payload']),
                        dependencies=entry['dependencies'],provenance=dict(protocol='meta_trace',
                            scope='Duration-free source operations. No association timings, service prices or live availability.'))
                except OSError:
                    self._counts['cache_errors']+=1;continue
                entry['pending']=False;stored.append(record)
            self._publication=dict(status='partial' if any(entry['pending'] for entry in self._entries.values()) else 'published',stored=stored)
            return deepcopy(self._publication)

    def snapshot(self):
        with self._lock:
            return dict(self._counts,entries=len(self._entries),estimated_retained_bytes=self._bytes,
                pending=sum(entry['pending'] for entry in self._entries.values()),max_entries=self.max_entries,
                max_bytes=self.max_bytes,closed=self._closed,publication=deepcopy(self._publication),
                scope='Structural meta traces only. Retained memory is bounded; cache filesystem lookup/publication costs must also be charged. No timing-price freshness or live admission certificate.')

    def close(self):
        with self._lock:self._closed=True;self._entries.clear();self._bytes=0


def cached_tensor_trace(name,bindings,build,*,implementation):
    active=_ACTIVE.get()
    return build() if active is None else active._get(name,bindings,build,implementation)
