"""Count eager tensor operations from the implementation, without timing it.

Meta tensors carry shape/dtype/storage identities but allocate no genotype data.
The graph includes aliasing, so views do not create fictitious device traffic.
Source-level bytes are deliberately separate from physical cache traffic.
"""
from __future__ import annotations
from collections import Counter, OrderedDict
import hashlib
from pathlib import Path

# Only an explicitly scoped planning invocation retains duration-free traces.
# Context-local storage isolates concurrent planners and is dropped on failure.
from contextvars import ContextVar
from contextlib import contextmanager
from functools import wraps
import copy

_STATISTICS_WORK_CACHE = ContextVar('torchgwas_statistics_work_cache', default=None)


@contextmanager
def tensor_work_cache(max_entries=128):
 """Reuse exact meta-operation ledgers within one bounded planning call."""
 if isinstance(max_entries,bool) or not isinstance(max_entries,int) or max_entries<1:
  raise ValueError('Positive tensor-work cache bound required')
 token=_STATISTICS_WORK_CACHE.set((OrderedDict(),max_entries))
 try:yield
 finally:_STATISTICS_WORK_CACHE.reset(token)


def reuse_tensor_work(function):
 @wraps(function)
 def run(*args,**kwargs):
  with tensor_work_cache():return function(*args,**kwargs)
 return run


def eager_statistics_work(samples, markers, traits=1, covariates=8, validate_range=False):
 state=_STATISTICS_WORK_CACHE.get()
 if state is None:
  return _saved_statistics_work(samples,markers,traits,covariates,validate_range)
 import torch
 from .linear import _dosage_statistics
 cache,maximum=state
 # Include dtype and function identity: neither prices nor device geometry are
 # cached, and a caller's changed trace implementation cannot reuse old work.
 key=(tuple((type(v),v) for v in (samples,markers,traits,covariates,validate_range)),
      torch.get_default_dtype(),_dosage_statistics,_trace_statistics_work)
 if key not in cache:
  value=_saved_statistics_work(samples,markers,traits,covariates,validate_range)
  if len(cache)>=maximum:cache.popitem(last=False)
  cache[key]=value
 cache.move_to_end(key)
 # Consumers receive their own mutable ledger; no previous result can poison
 # a later prediction or alter another candidate's source/memory accounting.
 return copy.deepcopy(cache[key])


def _saved_statistics_work(samples,markers,traits,covariates,validate_range):
 from .structural_tensor_cache import cached_tensor_trace
 from .linear import _dosage_statistics
 return cached_tensor_trace('statistics.v1',dict(samples=samples,markers=markers,traits=traits,
     covariates=covariates,validate_range=validate_range),
     lambda:_trace_statistics_work(samples,markers,traits,covariates,validate_range),
     implementation=(_trace_statistics_work,_dosage_statistics))


def _trace_statistics_work(samples, markers, traits=1, covariates=8, validate_range=False):
 import torch
 from torch.utils._python_dispatch import TorchDispatchMode
 from torch.overrides import TorchFunctionMode
 from torch.utils._pytree import tree_flatten
 from .linear import _dosage_statistics
 if min(samples,markers,traits)<1 or covariates<0:raise ValueError('Invalid dimensions')
 tensors={};storages={};steps=[];keep=[];host_calls=[];active_call=[]
 class HostRecord(TorchFunctionMode):
  def __torch_function__(self,func,types,args=(),kwargs=None):
   call_id=len(host_calls);kwargs=kwargs or {}
   inputs=[v for v in tree_flatten((args,kwargs))[0] if isinstance(v,torch.Tensor)]
   host_calls.append(dict(id=call_id,name=getattr(func,"__name__",str(func)),
    input_dtypes=[str(v.dtype) for v in inputs],input_shapes=[list(v.shape) for v in inputs],
    kwargs={k:str(v) for k,v in kwargs.items() if not isinstance(v,torch.Tensor)},step_indices=[]))
   active_call.append(call_id)
   try:return func(*args,**kwargs)
   finally:active_call.pop()
 def describe(t):
  key=t.untyped_storage()._cdata
  if key not in storages:storages[key]=len(storages)
  identifier=storages[key]
  return dict(storage=identifier,shape=list(t.shape),stride=list(t.stride()),
              dtype=str(t.dtype),device=str(t.device),offset_bytes=t.storage_offset()*t.element_size(),bytes=t.numel()*t.element_size(),storage_bytes=t.untyped_storage().nbytes())
 class Record(TorchDispatchMode):
  def __torch_dispatch__(self,func,types,args=(),kwargs=None):
   kwargs=kwargs or {};inputs=[v for v in tree_flatten((args,kwargs))[0] if isinstance(v,torch.Tensor)]
   inp=[describe(t) for t in inputs];out=func(*args,**kwargs)
   outputs=[v for v in tree_flatten(out)[0] if isinstance(v,torch.Tensor)];op=[describe(t) for t in outputs]
   # Keep meta identities alive, so freed storage handles cannot be recycled.
   keep.extend(outputs)
   name=str(func);view=bool(op) and all(v['storage'] in {a['storage'] for a in inp} for v in op) and not func._schema.is_mutable
   allocation_only=name.startswith(('aten.empty','aten.new_empty'))
   shape_only=name.startswith(('aten.zeros_like','aten.ones_like','aten.full_like','aten.empty_like'))
   scalar_argument=name.startswith('aten.scalar_tensor')
   read={} if shape_only else {v['storage']:v for v in inp if v['device']=='meta'}
   write={v['storage']:v for v in op if v['device']=='meta'}
   # scalar_tensor on meta corresponds to a CUDA scalar fill, as the launch census verifies.
   reads=0 if view or allocation_only else sum(v['bytes'] for v in read.values())
   writes=0 if view or allocation_only else sum(v['bytes'] for v in write.values())
   # A unary input used twice (x*x) is one logical read, not two physical loads.
   call_id=active_call[-1] if active_call else None
   steps.append(dict(op=name,inputs=inp,outputs=op,alias_only=view,allocation_only=allocation_only,host_call_id=call_id,
       read_bytes=reads,write_bytes=writes,logical_bytes=reads+writes,shape_only_inputs=shape_only))
   return out
 x=torch.empty((markers,samples),dtype=torch.int8,device='meta')
 design=torch.empty((samples,traits+covariates+1),device='meta')
 ss=torch.empty(traits,device='meta')
 for t in [x,design,ss]:describe(t)
 with HostRecord(),Record():
  genotype=torch.where(x==-9,torch.nan,x.to(torch.float32))
  conversion_end=len(steps)
  result=_dosage_statistics(genotype,design,ss,traits,samples-covariates-2,validate_range,covariate_rank=covariates)
 for i,step in enumerate(steps):step['phase']='conversion' if i<conversion_end else 'statistics'
 # CUDA sum promotes bool input to a materialized int64 array before reduction.
 # This internal ATen call is hidden from TorchDispatchMode. Its existence and
 # dtype are checked against the duration-free compiled-kernel census.
 expanded=[];next_storage=len(storages)
 for step in steps:
  if step['op']=='aten.sum.dim_IntList' and step['inputs'][0]['dtype']=='torch.bool':
   original=step['inputs'][0];cast=dict(original,storage=next_storage,dtype='torch.int64',
       bytes=8*original['bytes'],storage_bytes=8*original['bytes'],offset_bytes=0)
   next_storage+=1
   expanded.append(dict(op='cuda.sum.bool_to_int64',inputs=[original],outputs=[cast],
       alias_only=False,allocation_only=False,host_call_id=step['host_call_id'],read_bytes=original['bytes'],write_bytes=cast['bytes'],
       logical_bytes=original['bytes']+cast['bytes'],shape_only_inputs=False,phase=step['phase'],
       source='PyTorch CUDA sum dtype promotion, verified in kernel census; no timing coefficient'))
   step=dict(step,inputs=[cast],read_bytes=cast['bytes'],logical_bytes=cast['bytes']+step['write_bytes'])
  expanded.append(step)
 steps=expanded
 for i,step in enumerate(steps):
  if step['host_call_id'] is not None:host_calls[step['host_call_id']]['step_indices'].append(i)
 # Tensor property access is intercepted too, but contributes no ATen work.
 # Preserve it in a separate count; active API calls include views/allocations.
 property_calls=[call for call in host_calls if not call['step_indices']]
 host_calls=[call for call in host_calls if call['step_indices']]
 hashes={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).with_name('linear.py'),Path(__file__).with_name('native_scan.py')]}
 return dict(samples=samples,markers=markers,traits=traits,covariates=covariates,source_sha256=hashes,
   steps=steps,host_calls=host_calls,property_access_count=len(property_calls),logical_bytes=sum(s['logical_bytes'] for s in steps),
   torch_version=torch.__version__,tensor_dispatches=len(steps),active_tensor_operations=sum(not(s['alias_only'] or s['allocation_only']) for s in steps),
   result_bytes=sum(t.numel()*t.element_size() for t in result),
   initial_storages=[describe(t) for t in [x,design,ss]],result_storages=[describe(t) for t in result],
   scope='Eager int8-conversion/FP32 statistics operation graph from meta execution; no elapsed measurements. Logical storage accesses do not assert physical HBM traffic or a one-to-one CUDA-kernel count.')
