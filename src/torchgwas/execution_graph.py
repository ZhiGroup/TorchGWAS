"""Finite source execution schedules, independent of timing coefficients.

Nodes take externally derived service demands. No observations, fits, or
steady-state bottleneck approximation enter these dependency calculations.
"""
from __future__ import annotations
from collections import defaultdict,deque
import math

class ExecutionGraph:
 def __init__(self):self.nodes={};self.demands={};self.capacities={};self.token_capacities={};self.token_actions={};self.fifo_enqueues={};self.fifo_dequeues={};self.resource_waits={};self.conditional_delays={}
 def compose(self,other,prefix,after=(),*,shared_tokens=()):
  """Append an independent graph, preserving its waits and private tokens.

  Resource capacities remain shared. Node names, FIFO queues and token names
  are private except explicitly shared tokens. Every root waits for ``after``; the returned
  completion node includes asynchronous writers and release workers.
  """
  if not prefix or any(prefix+name in self.nodes for name in other.nodes):
   raise ValueError('Unique nonempty graph prefix required')
  if prefix+'complete' in self.nodes or 'complete' in other.nodes:
   raise ValueError('Graph completion name collision')
  for key,value in other.capacities.items():
   if key in self.capacities and self.capacities[key]!=value:
    raise ValueError('Conflicting shared resource capacity: '+key)
  self.capacities.update(other.capacities)
  def token_key(key):return key if key in shared_tokens else prefix+key
  for key,value in other.token_capacities.items():
   mapped=token_key(key)
   if mapped in self.token_capacities and (key not in shared_tokens or self.token_capacities[mapped]!=value):
    raise ValueError('Token name collision or conflicting shared capacity')
   self.token_capacities[mapped]=value
  for name,(seconds,deps) in other.nodes.items():
   self.add(prefix+name,seconds,[prefix+d for d in deps] if deps else after,
            resources=other.demands.get(name))
  for name,actions in other.token_actions.items():
   self.token_actions[prefix+name]={action:{token_key(k):v for k,v in values.items()} for action,values in actions.items()}
  for name,(queue,consumer) in other.fifo_enqueues.items():
   self.fifo_enqueues[prefix+name]=(prefix+queue,prefix+consumer)
  for name,queue in other.fifo_dequeues.items():self.fifo_dequeues[prefix+name]=prefix+queue
  for name,wait in other.resource_waits.items():
   self.resource_waits[prefix+name]=dict(after=prefix+wait['after'],until=prefix+wait['until'],resources=dict(wait['resources']))
  for name,condition in other.conditional_delays.items():
   self.conditional_delays[prefix+name]=dict(attempt=prefix+condition['attempt'],ready=prefix+condition['ready'])
  predecessors={dep for _,deps in other.nodes.values() for dep in deps}
  return self.add(prefix+'complete',after=[prefix+name for name in other.nodes if name not in predecessors] or after)
 def add(self,name,seconds=0.,after=(),resources=None):
  if name in self.nodes:raise ValueError('Duplicate node '+name)
  if not math.isfinite(seconds) or seconds<0:raise ValueError('Invalid service '+name)
  self.nodes[name]=(float(seconds),tuple(after))
  if resources:self.demands[name]=dict(resources)
  return name
 def add_wait_delay(self,name,seconds,attempt,ready):
  """Charge a supplied resume delay only if readiness follows the wait attempt.

  Endpoint timestamps, not a fitted blocking probability, select the branch.
  Supplied delay is elapsed service; wake-path CPU remains separately priced.
  """
  self.add(name,seconds,[attempt,ready])
  self.conditional_delays[name]=dict(attempt=attempt,ready=ready)
  return name
 def with_serial_sections(self,position,resource='host_serial',cpu_resource='cpu'):
  """Represent held CPU as non-preemptible critical sections, not fluid demand.

  Only total held/detached CPU is supplied, so placing each held section first
  or last is an explicit ordering scenario, not a recovered API timeline or
  guaranteed envelope. Other CPU/DRAM demand and node completion endpoints
  are conserved. CPU contention can extend a section while its token is held.
  """
  import copy
  if position not in ('held-first','held-last'):raise ValueError('Unknown serial section order')
  if resource not in self.capacities:raise ValueError('Serial section resource is missing')
  if self.capacities[resource]!=1.:raise ValueError('Serial section capacity must be one')
  if any(wait['resources'].get(resource,0.) for wait in self.resource_waits.values()):
   raise ValueError('Serial resource waits require explicit critical-section ownership')
  token='exclusive:'+resource
  if token in self.token_capacities:raise ValueError('Serial sections already composed')
  # Nodes and demands are rebuilt below. Copy only retained state, preserving
  # its internal aliases, rather than duplicating and discarding large graphs.
  g=copy.copy(self);g.nodes={};g.demands={}
  g.__dict__.update(copy.deepcopy({key:value for key,value in self.__dict__.items()
                                 if key not in ('nodes','demands')}))
  g.token_capacities[token]=1;entries={}
  for name,(seconds,deps) in self.nodes.items():
   demands=dict(self.demands.get(name,{}));held=demands.pop(resource,0.)
   cpu=demands.get(cpu_resource,0.)
   if not math.isfinite(held) or held<0 or held>cpu+1e-12:
    raise ValueError('Serial CPU demand must lie within node CPU demand')
   if not held or not seconds:
    g.add(name,seconds,deps,resources=demands);continue
   if name in self.conditional_delays:raise ValueError('Conditional delay cannot contain serial CPU work')
   held_seconds=seconds*min(1.,held/cpu)
   detached_seconds=seconds-held_seconds
   critical=name
   if detached_seconds:
    first=name+':serial_section_start'
    if first in self.nodes:raise ValueError('Serial section node name collision')
    durations=(held_seconds,detached_seconds) if position=='held-first' else (detached_seconds,held_seconds)
    g.add(first,durations[0],deps,resources=demands)
    g.add(name,durations[1],[first],resources=demands)
    critical=first if position=='held-first' else name
    actions=g.token_actions.pop(name,{})
    early={key:value for key,value in actions.items() if key!='release_finish'}
    late={key:value for key,value in actions.items() if key=='release_finish'}
    if early:g.token_actions[first]=early
    if late:g.token_actions[name]=late
    if name in g.fifo_dequeues:
     g.fifo_dequeues[first]=g.fifo_dequeues.pop(name);entries[name]=first
   else:g.add(name,seconds,deps,resources=demands)
   actions=g.token_actions.setdefault(critical,{})
   actions.setdefault('acquire',{})[token]=1
   actions.setdefault('release_finish',{})[token]=1
  g.fifo_enqueues={name:(key,entries.get(consumer,consumer)) for name,(key,consumer) in g.fifo_enqueues.items()}
  return g
 def solve(self):
  if self.capacities or self.token_capacities or self.resource_waits or self.conditional_delays:return self._solve_shared()
  degree={};children=defaultdict(list);starts={};ends={};ready=deque()
  for name,(seconds,deps) in self.nodes.items():
   deps=set(deps);degree[name]=len(deps);starts[name]=0.
   for dep in deps:
    if dep not in self.nodes:raise ValueError('Unknown predecessor '+dep)
    children[dep].append(name)
   if not deps:ready.append(name)
  while ready:
   name=ready.popleft();ends[name]=starts[name]+self.nodes[name][0]
   for child in children[name]:
    starts[child]=max(starts[child],ends[name]);degree[child]-=1
    if not degree[child]:ready.append(child)
  if len(ends)!=len(self.nodes):raise ValueError('Cyclic execution dependencies')
  return dict(seconds=max(ends.values(),default=0.),start=starts,end=ends)

 def checkpoint(self,at_seconds,*,initial=None):
  """Advance the analytical schedule to a decision time, retaining live state.

  This models a prefix; it does not observe actual CUDA/reader progress. With
  initial, continue that checkpoint instead of charging its completed work.
  Continuation currently uses the Python event solver, never a fresh-run fit.
  """
  if isinstance(at_seconds,bool) or not math.isfinite(at_seconds) or at_seconds<0:
   raise ValueError('Finite nonnegative checkpoint time required')
  return self._solve_shared_python(_initial=initial,_checkpoint_at=at_seconds)['checkpoint']

 def resume(self,checkpoint):
  """Complete a compatible graph from paid work, queues and remaining service.

  Resource capacities and unstarted future nodes may change. Started node and
  token contracts remain fixed. Caller must also preserve issued source ranges
  and physical allocations, which this generic scheduler cannot identify.
  """
  return self._solve_shared_python(_initial=checkpoint)

 def solve_chains(self,chains,*,shared_tokens=(),trace=False,prefix='part'):
  """One event scheduler for independent chains of fully drained graphs.

  Each chain contains (global_order, graph) pairs. Only its current graph is
  materialized; completion joins all leaves before admitting the next graph.
  Capacities and explicitly named tokens stay global. Templates can be reused,
  but their interacting event schedules are never solved independently.
  ``trace`` retains all timestamps for audits; normal operation retains only
  active timestamps and summed handoff service.
  """
  if any(getattr(self,key) for key in ('nodes','demands','token_capacities','token_actions',
         'fifo_enqueues','fifo_dequeues','resource_waits','conditional_delays')):
   raise ValueError('Chain solver requires an empty graph with shared capacities')
  if type(trace) is not bool or not isinstance(prefix,str) or not prefix:
   raise ValueError('Explicit trace boolean and nonempty prefix required')
  chains=[list(chain) for chain in chains];seen=set();capacities=dict(self.capacities)
  for chain in chains:
   previous=-1
   for order,graph in chain:
    if isinstance(order,bool) or not isinstance(order,int) or order<0 or order in seen or order<=previous:
     raise ValueError('Unique increasing graph order required in each chain')
    seen.add(order);previous=order
    for key,value in graph.capacities.items():
     if key in capacities and capacities[key]!=value:raise ValueError('Conflicting shared resource capacity: '+key)
     capacities[key]=value
  active=ExecutionGraph();active.capacities=capacities
  return active._solve_shared(_chains=sorted((c for c in chains if c),key=lambda c:c[0][0]),
                              _shared_tokens=set(shared_tokens),_trace=trace,_prefix=prefix)

 def _solve_shared(self,*,_chains=None,_shared_tokens=(),_trace=False,_prefix='part'):
  from .native_event_solver import solve
  result=solve(self,_chains,_shared_tokens,_trace,_prefix)
  if result is not None:return result
  return self._solve_shared_python(_chains=_chains,_shared_tokens=_shared_tokens,_trace=_trace,_prefix=_prefix)

 def _solve_shared_python(self,*,_chains=None,_shared_tokens=(),_trace=False,_prefix='part',
                          _initial=None,_checkpoint_at=None):
  """Fluid proportional sharing, with identical event rules for full/lazy DAGs.

  Each node has work-equivalent duration at its nominal service capacity and
  resource demand per second at that rate. Contending nodes slow together;
  completion changes the active set and releases capacity immediately.
  This is an explicit scheduler scenario, not a model of OS priorities.
  """
  if (_initial is not None or _checkpoint_at is not None) and _chains is not None:
   raise ValueError('Checkpoint continuation requires a materialized graph')
  children=defaultdict(list);degree={};ready=deque();remaining={};starts={};ends={}
  conditional={};wait_usage=defaultdict(float);waits={}
  wait_after=defaultdict(list);wait_until=defaultdict(list);active_wait_index={}
  fifo=defaultdict(deque);free={};parts={};added=finished=peak_nodes=0
  handoff_totals=dict(possible_waits=0,blocked_waits=0,extra_elapsed_service_seconds=0.)
  trace_starts,trace_ends={},{}
  def register(names,wait_items,token_keys):
   nonlocal added,peak_nodes
   added+=len(names);peak_nodes=max(peak_nodes,len(self.nodes))
   for name in names:
    seconds,deps=self.nodes[name]
    degree[name]=len(set(deps))
    if not deps:ready.append(name)
    for dep in set(deps):
     if dep not in self.nodes:raise ValueError('Unknown predecessor '+dep)
     children[dep].append(name)
    for resource,demand in self.demands.get(name,{}).items():
     if not math.isfinite(demand) or demand<0:raise ValueError('Invalid resource demand')
     if demand and self.capacities.get(resource,0)<=0:raise ValueError('Missing or zero resource '+resource)
    if name in self.conditional_delays:
     condition=self.conditional_delays[name];endpoints=(condition['attempt'],condition['ready'])
     if any(endpoint not in self.nodes for endpoint in endpoints):raise ValueError('Unknown conditional-delay endpoint')
     if not set(endpoints)<=set(deps):raise ValueError('Conditional-delay endpoints must be direct dependencies')
   for index,wait in wait_items:
    waits[index]=wait
    if wait['after'] not in self.nodes or wait['until'] not in self.nodes:raise ValueError('Unknown resource-wait endpoint')
    for resource,demand in wait['resources'].items():
     if not math.isfinite(demand) or demand<0 or (demand and self.capacities.get(resource,0)<=0):
      raise ValueError('Invalid resource-wait demand or missing capacity')
    wait_after[wait['after']].append(index);wait_until[wait['until']].append(index)
   for key in token_keys:
    value=self.token_capacities[key]
    if isinstance(value,bool) or not isinstance(value,int) or value<1:raise ValueError('Token capacities must be positive integers')
    if key not in free:free[key]=value
   for name in names:
    for kind,requests in self.token_actions.get(name,{}).items():
     if kind not in ('acquire','release_start','release_finish'):raise ValueError('Unknown token action')
     for key,amount in requests.items():
      if key not in free or isinstance(amount,bool) or not isinstance(amount,int) or not 0<amount<=self.token_capacities[key]:
       raise ValueError('Invalid token request')
  def admit(chain):
   item=next(chain,None)
   if item is None:return
   order,graph=item;prefix=f'{_prefix}{order}:'
   if set(graph.conditional_delays)-graph.nodes.keys():raise ValueError('Unknown conditional-delay node')
   if set(graph.token_actions)-graph.nodes.keys():raise ValueError('Unknown token node')
   complete=self.compose(graph,prefix,shared_tokens=_shared_tokens)
   names=[prefix+name for name in graph.nodes]+[complete]
   tokens=[key if key in _shared_tokens else prefix+key for key in graph.token_capacities]
   wait_items=[((order,i),self.resource_waits[prefix+name]) for i,name in enumerate(graph.resource_waits)]
   queues={prefix+q for q,_ in graph.fifo_enqueues.values()}|{prefix+q for q in graph.fifo_dequeues.values()}
   parts[complete]=(names,tokens,[key for key,_ in wait_items],queues,chain,
                    [prefix+name for name in graph.resource_waits])
   register(names,wait_items,tokens)
  def completed(name):
   nonlocal finished
   finished+=1
   if _trace:trace_ends[name]=ends[name]
   if name not in parts:return
   names,tokens,wait_keys,queues,chain,wait_names=parts.pop(name)
   if any(key not in ends for key in names):raise ValueError('Cyclic or incomplete graph at chain completion')
   # No other part can depend on these private nodes. Shared token balances
   # persist, including tokens currently held by nodes in another chain.
   for key in names:
    for values in (self.nodes,self.demands,self.token_actions,self.fifo_enqueues,self.fifo_dequeues,
                   self.conditional_delays,degree,children,starts,ends,wait_after,wait_until):values.pop(key,None)
   for key in tokens:
    if key not in _shared_tokens:
     self.token_capacities.pop(key,None);free.pop(key,None)
   for key in wait_keys:
    if key in active_wait_index:raise ValueError('Graph completed with an active resource wait')
    waits.pop(key)
   for key in wait_names:self.resource_waits.pop(key)
   for key in queues:fifo.pop(key,None)
   admit(chain)
  if _chains is None:
   if set(self.conditional_delays)-self.nodes.keys():raise ValueError('Unknown conditional-delay node')
   if set(self.token_actions)-self.nodes.keys():raise ValueError('Unknown token node')
   register(list(self.nodes),list(enumerate(self.resource_waits.values())),list(self.token_capacities))
  else:
   for chain in _chains:admit(iter(chain))
  def update_waits(name):
   for index in wait_until.get(name,()):active_wait_index.pop(index,None)
   for index in wait_after.get(name,()):
    if waits[index]['until'] not in ends:active_wait_index[index]=waits[index]
  def enqueue(name):
   if name in self.fifo_enqueues:
    key,consumer=self.fifo_enqueues[name];fifo[key].append(consumer)
  def release(name,kind):
   for key,amount in self.token_actions.get(name,{}).get(kind,{}).items():
    free[key]+=amount
    if free[key]>self.token_capacities[key]:raise ValueError('Unbalanced token release')
  now=0.;steps=0;origin=0.
  if _initial is not None:
   from .execution_checkpoint import validate_checkpoint
   state=validate_checkpoint(self,_initial)
   now=origin=state['at_seconds'];starts.update(state['start']);ends.update(state['end'])
   if _checkpoint_at is not None and _checkpoint_at<now:raise ValueError('Checkpoint cannot move backward')
   remaining.update(state['remaining_service_seconds']);free.update(state['token_free'])
   fifo.update({key:deque(values) for key,values in state['fifo'].items()})
   wait_usage.update(state['wait_resource_seconds']);conditional.update(state['conditional_delays'])
   steps=state['resource_event_steps'];finished=len(ends)
   for name,(_,deps) in self.nodes.items():degree[name]=sum(dep not in ends for dep in set(deps))
   # Preserve token/FIFO admission order for existing ready nodes. Newly added
   # future work joins behind them; dependencies are rechecked in the new DAG.
   eligible={name for name in self.nodes if name not in starts and not degree[name]}
   ordered=[name for name in state['ready'] if name in eligible]
   ready=deque(ordered+[name for name in self.nodes if name in eligible and name not in ordered])
   active_wait_index.update({index:wait for index,wait in waits.items()
                            if wait['after'] in ends and wait['until'] not in ends})
  def checkpoint_result():
   from .execution_checkpoint import make_checkpoint
   return dict(checkpoint=make_checkpoint(self,now,starts,ends,remaining,ready,free,fifo,
       wait_usage,conditional,steps))
  while ready or remaining:
   while ready:
    progressed=False
    for _ in range(len(ready)):
     name=ready.popleft()
     queue_key=self.fifo_dequeues.get(name)
     if queue_key is not None and (not fifo[queue_key] or fifo[queue_key][0]!=name):
      ready.append(name);continue
     acquire=self.token_actions.get(name,{}).get('acquire',{})
     if any(free[key]<amount for key,amount in acquire.items()):
      ready.append(name);continue
     for key,amount in acquire.items():free[key]-=amount
     if queue_key is not None:fifo[queue_key].popleft()
     release(name,'release_start');progressed=True
     starts[name]=now;duration=self.nodes[name][0]
     if _trace:trace_starts[name]=now
     if name in self.conditional_delays:
      condition=self.conditional_delays[name]
      waited=ends[condition['ready']]-ends[condition['attempt']]
      blocked=waited>1e-12
      if not blocked:duration=0.
      if _chains is None or _trace:
       conditional[name]=dict(blocked=blocked,dependency_wait_seconds=max(0.,waited),extra_elapsed_service_seconds=duration)
      if _chains is not None:
       handoff_totals['possible_waits']+=1;handoff_totals['blocked_waits']+=blocked
       handoff_totals['extra_elapsed_service_seconds']+=duration
     if duration:remaining[name]=duration
     else:
      ends[name]=now;update_waits(name);enqueue(name);release(name,'release_finish')
      for child in children[name]:
       degree[child]-=1
       if not degree[child]:ready.append(child)
      completed(name)
    if not progressed:break
   if _checkpoint_at is not None and now>=_checkpoint_at:return checkpoint_result()
   if not remaining:
    if ready:raise ValueError('Token dependency deadlock')
    continue
   totals=defaultdict(float)
   for name in remaining:
    for resource,demand in self.demands.get(name,{}).items():totals[resource]+=demand
   # Global graph order, not admission order, preserves the summation order
   # when a faster device admits a later tile before another device's tile.
   active_waits=[active_wait_index[index] for index in sorted(active_wait_index)]
   for wait in active_waits:
    for resource,demand in wait['resources'].items():totals[resource]+=demand
   rates={name:min([1.]+[self.capacities[key]/totals[key] for key,value in self.demands.get(name,{}).items() if value and totals[key]]) for name in remaining}
   interval=min(value/rates[name] for name,value in remaining.items())
   cut=_checkpoint_at is not None and now+interval>_checkpoint_at
   if cut:interval=_checkpoint_at-now
   now+=interval;steps+=1
   for wait in active_waits:
    rate=min([1.]+[self.capacities[key]/totals[key] for key,value in wait['resources'].items() if value and totals[key]])
    for resource,demand in wait['resources'].items():wait_usage[resource]+=interval*rate*demand
   done=[]
   for name in remaining:
    remaining[name]-=interval*rates[name]
    if not cut and remaining[name]<=max(1e-14,1e-12*self.nodes[name][0]):done.append(name)
   if cut:return checkpoint_result()
   for name in done:
    remaining.pop(name);ends[name]=now;update_waits(name);enqueue(name);release(name,'release_finish')
    for child in children[name]:
     degree[child]-=1
     if not degree[child]:ready.append(child)
    completed(name)
  if finished!=added:raise ValueError('Cyclic execution dependencies')
  if _checkpoint_at is not None:return checkpoint_result()
  result=dict(seconds=now,start=trace_starts if _trace else starts,end=trace_ends if _trace else ends,
              resource_event_steps=steps,wait_resource_seconds=dict(wait_usage),conditional_delays=conditional,
              resource_policy='Proportional active-demand fluid sharing; capacity released at every completion.')
  if _chains is not None:result.update(conditional_summary=handoff_totals,peak_active_nodes=peak_nodes,scheduled_nodes=added)
  if _initial is not None:result.update(elapsed_before_resume=origin,remaining_seconds=now-origin)
  return result

def plink_block_schedule(blocks):
 """Main reads next, joins previous, spawns next, then formats previous.

 Per-block compute includes its actual worker service and shared-resource
 limits. Read and format/write share the main thread and cannot overlap.
 The final block is formatted only after its compute completes.
 """
 g=ExecutionGraph();main=None;previous=None
 for i,b in enumerate(blocks):
  read=g.add(f'read:{i}',b['read_seconds'],[main] if main else [])
  launch=g.add(f'launch:{i}',b.get('launch_seconds',0),[read]+([previous] if previous else []))
  current=g.add(f'compute:{i}',b['compute_seconds'],[launch])
  main=g.add(f'format:{i-1}',blocks[i-1]['format_write_seconds'],[launch]) if i else launch
  previous=current
 if blocks:g.add(f'format:{len(blocks)-1}',blocks[-1]['format_write_seconds'],[main,previous])
 result=g.solve();result['scope']='Finite read/join/spawn/format schedule; component service must be independently supplied.'
 return result

def torch_scan_schedule(blocks,depth=4,decode_workers=4,consumer=None,shared_capacities=None,return_graph=False,host_serial_policy='fluid',consumer_open_before_scan=False,synchronous_results=False):
 """Pinned decode/three CUDA streams/result ring; no writer approximation.

 Models the source's ordered producer and main-thread yields. Each block
 supplies read+decode duration for a worker, H2D, typed CPU submission/GPU
 operations, D2H, finish and synchronous consumer service. The latter is an
 explicit consumer (e.g. discard or synchronous binary writer); the production
 asynchronous coalescing writer must be composed separately, not hidden here.

 Decoder worker assignment is fixed cyclic; exact for equal full chunks and
 one final short chunk. With arbitrary compressed-record imbalance this is a
 declared scheduling scenario, not an assertion of ThreadPoolExecutor choice.
 """
 if type(synchronous_results) is not bool:raise ValueError('Synchronous result policy must be boolean')
 if any(b.get('result_contract')=='device_significant_status' for b in blocks):
  if not synchronous_results or getattr(consumer,'selection_graphs',None) is None:
   raise ValueError('Device-significant status blocks require synchronous indexed selection graphs')
 if not isinstance(depth,int) or depth<2:raise ValueError('Scan ring depth must be >=2')
 if not isinstance(decode_workers,int) or decode_workers<1:raise ValueError('Invalid decoder count')
 workers=min(depth,decode_workers);g=ExecutionGraph();length=len(blocks);compute_tail=None
 handoffs=any(b.get('handoff_wakeup_seconds') is not None for b in blocks)
 if handoffs:
  for b in blocks:
   values=b.get('handoff_wakeup_seconds')
   if values is None or set(values)!={'queue','future'}:
    raise ValueError('Every block must supply queue and future handoff service')
   if any(isinstance(value,bool) or not math.isfinite(value) or value<0 for value in values.values()):
    raise ValueError('Invalid handoff elapsed service')
  g.add('scan_begin')
 g.capacities=dict(shared_capacities or {})
 if any('host_submit_serial_cpu_seconds' in b for b in blocks):
  g.capacities.setdefault('host_serial',1.)
 for i,b in enumerate(blocks):
  discard=b.get('discard_seconds',0.)
  if isinstance(discard,bool) or not math.isfinite(discard) or discard<0:raise ValueError('Invalid owned-result discard service')
  if discard:
   if consumer is not None:raise ValueError('Owned-result discard requires explicit writer ownership, not immediate-discard service')
   if b.get('host_resources',{}).get('cpu',0.)<=0:raise ValueError('Discard requires explicit host CPU demand')
   g.capacities.setdefault('host_serial',1.)
  deps=[]
  if i:deps.append(f'submit_decode:{i-1}')
  if i>=workers:deps.append(f'publish:{i-workers}')
  if i>=depth:
   if handoffs:
    attempt=g.add(f'attempt_free:{i}',0.,deps or ['scan_begin'])
    deps=[g.add_wait_delay(f'wake_free:{i}',b['handoff_wakeup_seconds']['queue'],attempt,f'release:{i-depth}')]
   else:deps.append(f'release:{i-depth}')
  submit=g.add(f'submit_decode:{i}',b.get('decode_submit_seconds',0.),deps,resources=b.get('host_resources'))
  read_deps=[submit]+([f'decode:{i-workers}'] if i>=workers else [])
  if 'reader_init_seconds' in b:
   initialize=g.add(f'reader_init:{i}',b['reader_init_seconds'],read_deps,
                    resources=b.get('reader_init_resources'))
   read_deps=[initialize]
  read=g.add(f'decode_read:{i}',b.get('decode_read_seconds',0.),read_deps,resources=b.get('read_resources'))
  decode=g.add(f'decode:{i}',b['decode_seconds'],[read],resources=b.get('decode_resources'))
  # Producer submits a full pending window before emitting its oldest item.
  # The final window is emitted after the last submission.
  pubdeps=[decode,f'submit_decode:{min(i+workers-1,length-1)}']
  if i:pubdeps.append(f'publish:{i-1}')
  if handoffs:
   attempt=g.add(f'attempt_publish:{i}',0.,pubdeps[1:])
   pubdeps=[g.add_wait_delay(f'wake_publish:{i}',b['handoff_wakeup_seconds']['future'],attempt,decode)]
  publish=g.add(f'publish:{i}',b.get('publish_seconds',0.),pubdeps,resources=b.get('host_resources'))
  if handoffs:
   attempt=g.add(f'attempt_fetch:{i}',0.,[f'consume:{i-1}' if synchronous_results else f'host_done:{i-1}'] if i else ['scan_begin'])
   fetch=g.add_wait_delay(f'fetch:{i}',b['handoff_wakeup_seconds']['queue'],attempt,publish)
  else:fetch=publish
  main_deps=[fetch]
  if i:main_deps.append(f'host_done:{i-1}')
  if synchronous_results and i:main_deps.append(f'consume:{i-1}')
  elif i>=depth:main_deps.append(f'consume:{i-depth}')
  host=g.add(f'host_start:{i}',b.get('transfer_submit_seconds',0.),main_deps,resources=b.get('host_resources'))
  copy_deps=[host]
  if i:copy_deps.append(f'h2d:{i-1}')
  if i>=depth:copy_deps.append(f'compute_done:{i-depth}')
  copy=g.add(f'h2d:{i}',b['h2d_seconds'],copy_deps,resources=b.get('h2d_resources'))
  early_release=b.get('release_after_transfer',False)
  if type(early_release) is not bool:raise ValueError('Release submission policy must be boolean')
  release_submitted=host if early_release else f'host_done:{i}'
  release_deps=[copy,release_submitted]+([f'release:{i-1}'] if i else [])
  if b.get('event_wait_resources'):
   begin=g.add(f'begin_release_wait:{i}',0.,[release_submitted]+([f'release:{i-1}'] if i else []))
   g.resource_waits[f'release_wait:{i}']=dict(after=begin,until=copy,resources=dict(b['event_wait_resources']))
   release_deps.append(begin)
  g.add(f'release:{i}',b.get('release_seconds',0.),release_deps,resources=b.get('host_resources'))
  previous=compute_tail;last_submit=host;last_offset=0.;last_serial=0.
  serial_total=b.get('host_submit_serial_cpu_seconds')
  if serial_total is not None:
   cpu_rate=b.get('host_resources',{}).get('cpu')
   if cpu_rate is None or cpu_rate<=0 or not math.isfinite(cpu_rate):raise ValueError('Explicit positive host CPU demand required')
   if isinstance(serial_total,bool) or not math.isfinite(serial_total) or not 0<=serial_total<=b['host_submit_seconds']*cpu_rate+1e-12:
    raise ValueError('Serial host total must lie within total CPU work')
  for j,op in enumerate(b['operations']):
   offset=op['host_submit_finish']
   if offset<last_offset:raise ValueError('Host call offsets must be ordered')
   resources=dict(b.get('host_resources',{}));duration=offset-last_offset
   if serial_total is not None:
    serial=op.get('host_serial_cpu_finish')
    if serial is None or isinstance(serial,bool) or not math.isfinite(serial) or not last_serial<=serial<=serial_total+1e-12:
     raise ValueError('Ordered cumulative serial CPU work required')
    delta=serial-last_serial
    if delta>duration*cpu_rate+1e-12:raise ValueError('Serial API CPU exceeds total API CPU')
    resources['host_serial']=delta/duration if duration else 0.
    last_serial=serial
   last_submit=g.add(f'api:{i}:{j}',duration,[last_submit],resources=resources);last_offset=offset
   previous=g.add(f'kernel:{i}:{j}',op['kernel_service_seconds'],[last_submit,copy]+([previous] if previous else []))
  comp=g.add(f'compute_done:{i}',0.,[previous or copy]+([copy] if previous and not b['operations'] else []))
  compute_tail=comp
  tail=b['host_submit_seconds']-last_offset
  if tail< -1e-12:raise ValueError('Host total below last kernel submission')
  result_submit=b.get('result_submit_seconds',0.);duration=max(0.,tail)+result_submit
  resources=dict(b.get('host_resources',{}))
  if serial_total is not None:
   delta=max(0.,serial_total-last_serial)
   if delta>max(0.,tail)*cpu_rate+1e-12:raise ValueError('Serial host tail exceeds total tail CPU')
   resources['host_serial']=(delta+result_submit*resources.get('host_serial',0.))/duration if duration else 0.
  hostdone=g.add(f'host_done:{i}',duration,[last_submit],resources=resources)
  status_submit=hostdone
  if 'status_submit_seconds' in b:
   status_submit=g.add(f'status_submit:{i}',b['status_submit_seconds'],[hostdone],resources=b.get('status_submit_resources'))
  d2h=g.add(f'd2h:{i}',b['d2h_seconds'],[comp,status_submit]+([f'd2h:{i-1}'] if i else []),resources=b.get('d2h_resources'))
  finish_deps=[d2h]+([f'finish:{i-depth}'] if i>=depth else [])
  result_wait_resources=b.get('result_wait_resources',b.get('event_wait_resources'))
  if result_wait_resources:
   begin=g.add(f'begin_finish_wait:{i}',0.,[status_submit]+([f'finish:{i-depth}'] if i>=depth else []))
   g.resource_waits[f'finish_wait:{i}']=dict(after=begin,until=d2h,resources=dict(result_wait_resources))
   finish_deps.append(begin)
  phases=b.get('finish_operations')
  if phases is not None:
   if not phases or any(not math.isfinite(op['seconds']) or op['seconds']<0 for op in phases):
    raise ValueError('Nonempty finite nonnegative finish phases required')
   if not math.isclose(sum(op['seconds'] for op in phases),b['finish_seconds'],rel_tol=1e-10,abs_tol=1e-12):
    raise ValueError('Finish phases must conserve total service')
   for j,op in enumerate(phases):
    resources=dict(b.get('host_resources',{}));resources.update(op.get('resources',{}))
    fraction=op.get('host_serial_fraction')
    if fraction is not None:
     if isinstance(fraction,bool) or not math.isfinite(fraction) or not 0<=fraction<=1:
      raise ValueError('Invalid finish-phase host serial fraction')
     resources['host_serial']=resources.get('cpu',0.)*fraction
     if fraction:g.capacities.setdefault('host_serial',1.)
    phase=g.add(f'finish_phase:{i}:{j}',op['seconds'],finish_deps,resources=resources)
    finish_deps=[phase]
   finish=g.add(f'finish:{i}',0.,finish_deps)
  else:
   finish=g.add(f'finish:{i}',b['finish_seconds'],finish_deps,resources=b.get('host_resources'))
  # Consumer is the same main Python iterator driver. It consumes an old
  # result when requesting the depth-th next block, and drains in order.
  trigger=f'host_done:{min(i+depth-1,length-1)}'
  if synchronous_results:
   # Device selection blocks in the main iterator after compute and status
   # D2H. It cannot submit another chunk until its last output yield resumes.
   # Decode/pinned-input release still run independently in the background.
   deps=[finish]
  else:
   deps=[finish,trigger]+([f'consume:{i-1}'] if i else [])
   if i+depth<length:deps.append(f'fetch:{i+depth}' if handoffs else f'publish:{i+depth}')
  if handoffs and not synchronous_results:
   attempt=g.add(f'attempt_resolve:{i}',0.,deps[1:])
   deps=[g.add_wait_delay(f'wake_resolve:{i}',b['handoff_wakeup_seconds']['future'],attempt,finish)]
  resolve=g.add(f'resolve:{i}',b.get('resolve_seconds',0.),deps,resources=b.get('host_resources'))
  deps=[resolve]
  if consumer is None:
   if discard:
    body=g.add(f'consume:{i}:body',b['consumer_seconds'],deps)
    cpu=b['host_resources']['cpu']
    g.add(f'consume:{i}',discard,[body],resources={'cpu':cpu,'host_serial':cpu})
   else:g.add(f'consume:{i}',b['consumer_seconds'],deps)
  else:
   consumer.append(g,i,deps)
   if synchronous_results and getattr(consumer,'gpu_tail',None) is not None:
    compute_tail=consumer.gpu_tail
 if consumer is not None:consumer.close(g,[f'consume:{length-1}'] if length else [])
 if consumer_open_before_scan:
  # Explicit block_bytes opens streams in BinarySumstatsWriter.__post_init__.
  # Automatic block sizing instead opens them on the first append.
  if consumer is None or 'writer:open' not in g.nodes:raise ValueError('Early open requires a binary writer')
  seconds,old_deps=g.nodes['writer:open'];g.nodes['writer:open']=(seconds,())
  for name,(duration,deps) in list(g.nodes.items()):
   if name=='writer:open':continue
   if 'writer:open' in deps:g.nodes[name]=(duration,tuple(dict.fromkeys((*deps,*old_deps))))
   elif not deps:g.nodes[name]=(duration,('writer:open',))
 if host_serial_policy!='fluid':g=g.with_serial_sections(host_serial_policy)
 if return_graph:return g
 result=g.solve();result.update(depth=depth,decode_workers=workers,blocks=length,
  scope='Finite scan-to-synchronous-consumer schedule only; initialization, asynchronous output staging and shared-resource contention require separate composition. No fitted service constants.')
 return result


def torch_multigpu_schedule(shards, shared_capacities, shared_links=(),
                           ordered=False, result_queue_depth=4, host_serial_fraction=None, queue_service=None,host_serial_policy='fluid',
                           borrow_results=False, acknowledgement_service=None):
 """Compose per-GPU source graphs with shared CPU/DRAM/input/link capacities.

 Shards supply blocks, depth and decode_workers; each CUDA device has its own
 stream dependencies. This models scan to immediate indexed discard only.
 ordered=True adds the legacy shard-order consumer and per-shard queue credits.
 Nonzero output-consumer work must be composed explicitly, not silently omitted.
 """
 if not shards:raise ValueError('At least one shard required')
 if type(borrow_results) is not bool:raise ValueError('Borrowed result policy must be boolean')
 ack=acknowledgement_service
 if borrow_results and len(shards)>1:
  if ack is None:raise ValueError('Borrowed multiGPU scheduling requires acknowledgement service')
  for key in ['create_cpu_seconds','publish_cpu_seconds','receive_cpu_seconds','wakeup_seconds']:
   if isinstance(ack[key],bool) or not math.isfinite(ack[key]) or ack[key]<0:raise ValueError('Invalid acknowledgement service')
 queue_service=dict(queue_service or {'put_cpu_seconds':0.,'get_cpu_seconds':0.,'cpu_fraction':1.})
 q=queue_service['cpu_fraction']
 if not math.isfinite(q) or not 0<q<=1:raise ValueError('Invalid queue CPU fraction')
 for key in ('put_cpu_seconds','get_cpu_seconds'):
  if not math.isfinite(queue_service[key]) or queue_service[key]<0:raise ValueError('Invalid queue CPU service')
 if host_serial_fraction is not None:
  if isinstance(host_serial_fraction,bool) or not math.isfinite(host_serial_fraction) or not 0<=host_serial_fraction<=1:
   raise ValueError('Host serial fraction must be in [0,1]')
 if not isinstance(result_queue_depth,int) or isinstance(result_queue_depth,bool) or result_queue_depth<1:
  raise ValueError('Positive result_queue_depth required')
 g=ExecutionGraph();g.capacities=dict(shared_capacities)
 if host_serial_fraction is not None or any('host_submit_serial_cpu_seconds' in b for shard in shards for b in shard['blocks']):g.capacities['host_serial']=1.
 names=[s['device'] for s in shards]
 if len(set(names))!=len(names):raise ValueError('One active shard graph per device required')
 g.token_capacities['result_consumer']=1
 if not ordered:g.token_capacities['result_queue']=result_queue_depth*len(shards)
 prepared=[]
 for index,shard in enumerate(shards):
  blocks=[dict(b) for b in shard['blocks']]
  if host_serial_fraction is not None:
   for b in blocks:
    host=dict(b.get('host_resources',{}))
    if 'cpu' not in host or not math.isfinite(host['cpu']) or host['cpu']<=0:
     raise ValueError('Explicit positive host CPU demand required for serial-host scenario')
    host['host_serial']=host['cpu']*host_serial_fraction
    b['host_resources']=host
  if not blocks:raise ValueError('Empty device shard')
  if any(b.get('consumer_seconds',0) for b in blocks):
   raise ValueError('Multi-GPU immediate-discard model does not price output consumers')
  if borrow_results and any(b.get('discard_seconds',0) for b in blocks):
   raise ValueError('Borrowed ring views cannot have owned-allocation discard service')
  for direction in ('h2d','d2h'):
   for b in blocks:b[direction+'_resources']=dict(b.get(direction+'_resources',{}))
  prepared.append(blocks)
 for i,link in enumerate(shared_links):
  if not link['devices'] or not set(link['devices'])<=set(names):raise ValueError('Unknown shared-link devices')
  for direction in ('h2d','d2h'):
   key=f'link:{i}:{direction}';capacity=link[direction+'_bytes_per_second']
   if not math.isfinite(capacity) or capacity<=0:raise ValueError('Positive shared-link capacity required')
   g.capacities[key]=capacity
   for shard,blocks in zip(shards,prepared):
    if shard['device'] not in link['devices']:continue
    for b in blocks:
     seconds=b[direction+'_seconds'];amount=b[direction+'_bytes']
     if amount<0 or not math.isfinite(amount) or (amount and seconds<=0):raise ValueError('Invalid transfer service')
     b[direction+'_resources'][key]=amount/seconds if seconds else 0.
 if len(shards)==1:
  shard=shards[0]
  result=torch_scan_schedule(prepared[0],depth=shard['depth'],decode_workers=shard['decode_workers'],shared_capacities=g.capacities,host_serial_policy=host_serial_policy)
  result.update(devices=names,ordered=ordered,result_queue_depth=0,
      resource_policy=dict(queue_service=None,result_queue_capacity=0,host_serial_fraction=host_serial_fraction,host_serial_policy=host_serial_policy,borrow_results=borrow_results),
      scope='Single active device uses the direct scan path; no multiGPU result queue, put/get service or queue credits.')
  return result
 for index,(shard,blocks) in enumerate(zip(shards,prepared)):
  local_blocks=[dict(b,consumer_seconds=(queue_service['put_cpu_seconds']+(ack['create_cpu_seconds'] if borrow_results else 0.))/q,discard_seconds=0.) for b in blocks]
  local=torch_scan_schedule(local_blocks,depth=shard['depth'],decode_workers=shard['decode_workers'],return_graph=True)
  if 'host_serial' in local.capacities:g.capacities.setdefault('host_serial',1.)
  prefix=f'shard{index}:'
  queue_key=prefix+'result_queue' if ordered else 'result_queue'
  if ordered:g.token_capacities[queue_key]=result_queue_depth
  for name,(seconds,deps) in local.nodes.items():
   g.add(prefix+name,seconds,[prefix+d for d in deps],resources=local.demands.get(name))
  for name,condition in local.conditional_delays.items():
   g.conditional_delays[prefix+name]=dict(attempt=prefix+condition['attempt'],ready=prefix+condition['ready'])
  for name,wait in local.resource_waits.items():
   g.resource_waits[prefix+name]=dict(after=prefix+wait['after'],until=prefix+wait['until'],resources=dict(wait['resources']))
  for i in range(len(blocks)):
   consume=prefix+f'consume:{i}';deps=[consume]
   if i:deps.append(prefix+f'deliver:{i-1}')
   elif ordered and index:
    deps.append(f'shard{index-1}:deliver:{len(prepared[index-1])-1}')
   deliver=prefix+f'deliver:{i}'
   resources={'cpu':q} if borrow_results or queue_service['get_cpu_seconds'] or queue_service['put_cpu_seconds'] else None
   if resources and host_serial_fraction is not None:resources['host_serial']=q*host_serial_fraction
   discard=blocks[i].get('discard_seconds',0.)
   if isinstance(discard,bool) or not math.isfinite(discard) or discard<0:raise ValueError('Invalid owned-result discard service')
   get=deliver+':get' if discard or borrow_results else deliver
   g.add(get,queue_service['get_cpu_seconds']/q,deps,resources=resources)
   if discard:
    cpu=blocks[i].get('host_resources',{}).get('cpu',0.)
    if cpu<=0:raise ValueError('Discard requires explicit host CPU demand')
    g.capacities.setdefault('host_serial',1.)
    g.add(deliver,discard,[get],resources={'cpu':cpu,'host_serial':cpu})
   elif borrow_results:
    g.add(deliver,ack['publish_cpu_seconds']/q,[get],resources=resources)
   if resources:g.demands[consume]=resources
   g.fifo_enqueues[consume]=(queue_key,get)
   g.fifo_dequeues[get]=queue_key
   g.token_actions[consume]={'acquire':{queue_key:1}}
   g.token_actions[get]={'acquire':{'result_consumer':1},'release_start':{queue_key:1}}
   g.token_actions.setdefault(deliver,{})['release_finish']={'result_consumer':1}
   if borrow_results:
    wake=g.add_wait_delay(prefix+f'ack_wake:{i}',ack['wakeup_seconds'],consume,deliver)
    resume=g.add(prefix+f'ack_resume:{i}',ack['receive_cpu_seconds']/q,[wake],resources=resources)
    # next(iterator) may immediately resolve a draining result, or submit the
    # next use of this ring slot. Neither may precede acknowledgement.
    successors=([prefix+f'resolve:{i+1}'] if i+1<len(blocks) else [])
    if prefix+f'attempt_resolve:{i+1}' in g.nodes:successors.append(prefix+f'attempt_resolve:{i+1}')
    if i+shard['depth']<len(blocks):successors.append(prefix+f'host_start:{i+shard["depth"]}')
    for successor in successors:
     seconds,deps=g.nodes[successor];g.nodes[successor]=(seconds,(*deps,resume))
 if host_serial_policy!='fluid':g=g.with_serial_sections(host_serial_policy)
 result=g.solve()
 result.update(devices=names,ordered=ordered,result_queue_depth=result_queue_depth,
     resource_policy=dict(queue_service=queue_service,result_queue_capacity=result_queue_depth if ordered else result_queue_depth*len(shards),host_serial_fraction=host_serial_fraction,host_serial_policy=host_serial_policy,borrow_results=borrow_results,acknowledgement_service=ack if borrow_results else None,
                         host_serial_capacity=1. if host_serial_fraction is not None else None,
                         host_serial_scope='Known API/allocator CPU overrides the supplied fraction for other host control. Policy '+host_serial_policy+' is an explicit scheduling scenario, not a recovered GIL timeline or guaranteed runtime bound.'),
     scope='Finite independent CUDA pipelines sharing declared CPU/DRAM/input/PCIe resources, to immediate indexed discard. Setup, reductions, writer and driver contention excluded; host serialization is optional and explicitly supplied.')
 return result
