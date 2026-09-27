"""Conditional remaining-work bounds without guessing live queue progress.

Bounds apply to the declared fluid-resource execution graph, not hardware.
Mandatory nodes are known unstarted. The upper bound includes every node's
full service, including work which may already be complete. This deliberately
overpays an unknown prefix rather than inventing an analytical checkpoint.
"""
from collections import defaultdict, deque
import math

from .execution_graph import ExecutionGraph


def remaining_work_bounds(graph, *, unissued_nodes, wait_chains=(), max_nodes=1000000,
                          max_reachability_visits=1000000):
    """Bound any feasible continuation with these source nodes unstarted.

    A wait chain proves mutually exclusive intervals: the previous wait's
    until endpoint must precede the next after endpoint in the dependency DAG.
    Unassigned waits each receive an independent slot, preserving correctness
    at the expense of a looser upper bound. No event simulation is required.
    """
    if not isinstance(graph, ExecutionGraph):raise ValueError('ExecutionGraph required')
    for value in (max_nodes,max_reachability_visits):
        if type(value) is not int or value<1:raise ValueError('Positive bounded work limits required')
    if len(graph.nodes)>max_nodes:raise ValueError('Remaining bound exceeds node budget')
    if not isinstance(unissued_nodes,(list,tuple)) or len(set(unissued_nodes))!=len(unissued_nodes):
        raise ValueError('Unique explicit unissued node list required')
    if any(name not in graph.nodes for name in unissued_nodes):raise ValueError('Unknown unissued source node')
    capacities=graph.capacities
    if any(isinstance(v,bool) or not math.isfinite(v) or v<=0 for v in capacities.values()):
        raise ValueError('Positive finite resource capacities required')
    def rates(values):
        if not isinstance(values,dict):raise ValueError('Explicit resource-demand map required')
        for key,value in values.items():
            if isinstance(value,bool) or not math.isfinite(value) or value<0 or key not in capacities:
                raise ValueError('Invalid resource demand or missing capacity')
        return {key:value/capacities[key] for key,value in values.items()}
    if set(graph.demands)-set(graph.nodes) or set(graph.conditional_delays)-set(graph.nodes):
        raise ValueError('Unknown service/conditional node')
    children=defaultdict(list);degrees={};ordered=[];ready=deque()
    for name,(duration,deps) in graph.nodes.items():
        if isinstance(duration,bool) or not math.isfinite(duration) or duration<0:
            raise ValueError('Finite nonnegative nominal service required')
        degrees[name]=len(set(deps))
        for dep in set(deps):
            if dep not in graph.nodes:raise ValueError('Unknown dependency')
            children[dep].append(name)
        if not deps:ready.append(name)
    while ready:
        name=ready.popleft();ordered.append(name)
        for child in children[name]:
            degrees[child]-=1
            if not degrees[child]:ready.append(child)
    if len(ordered)!=len(graph.nodes):raise ValueError('Acyclic service dependencies required')
    order={name:index for index,name in enumerate(ordered)}
    mandatory=set(unissued_nodes);pending=list(unissued_nodes)
    while pending:
        for child in children[pending.pop()]:
            if child not in mandatory:mandatory.add(child);pending.append(child)
    normalized={name:rates(graph.demands.get(name,{})) for name in graph.nodes}
    work={key:[] for key in capacities};ends={};potential=[]
    for name in ordered:
        duration,deps=graph.nodes[name];demand=normalized[name]
        # Potential of ALL service, including optional delays at full duration.
        potential.append(duration*(1.+math.fsum(demand.values())))
        seconds=duration if name in mandatory and name not in graph.conditional_delays else 0.
        for key,rate in demand.items():work[key].append(seconds*rate)
        ends[name]=max((ends[dep] for dep in deps),default=0.)+seconds*max([1.,*demand.values()])
    if not isinstance(wait_chains,(list,tuple)):raise ValueError('Explicit bounded wait chains required')
    waits=graph.resource_waits;seen=set();groups=[];visits=0
    for value in waits.values():
        if set(value)!={'after','until','resources'} or value['after'] not in graph.nodes or value['until'] not in graph.nodes:
            raise ValueError('Invalid active-wait endpoints')
        rates(value['resources'])
    for chain in wait_chains:
        if not isinstance(chain,(list,tuple)) or not chain:raise ValueError('Nonempty wait chains required')
        for index,name in enumerate(chain):
            if name not in waits or name in seen:raise ValueError('Unknown or duplicate chained wait')
            seen.add(name)
            if index:
                target=waits[chain[index-1]]['until'];stack=[waits[name]['after']];visited=set();found=False
                while stack:
                    node=stack.pop()
                    if node==target:found=True;break
                    if order[node]<order[target]:continue
                    if node in visited:continue
                    visited.add(node);visits+=1
                    if visits>max_reachability_visits:raise ValueError('Wait-chain proof exceeds traversal budget')
                    stack.extend(graph.nodes[node][1])
                if not found:raise ValueError('Wait-chain intervals are not dependency-ordered')
        groups.append(chain)
    groups.extend([name] for name in waits if name not in seen)
    wait_normalized={key:math.fsum(max((waits[name]['resources'].get(key,0.)/capacity
        for name in group),default=0.) for group in groups) for key,capacity in capacities.items()}
    active_wait_bound=math.fsum(wait_normalized.values())
    resource_floor={key:math.fsum(values) for key,values in work.items()}
    lower=max([0.,*resource_floor.values(),*ends.values()])
    upper=math.fsum(potential)*(1.+active_wait_bound)
    if not math.isfinite(lower) or not math.isfinite(upper):raise ValueError('Remaining work bound overflow')
    if lower>upper+max(1e-10,abs(upper)*1e-12):raise ValueError('Inconsistent remaining-work bounds')
    return dict(lower_seconds=lower,upper_seconds=upper,resource_floor_seconds=resource_floor,
        dependency_floor_seconds=max(ends.values(),default=0.),full_service_potential_seconds=math.fsum(potential),
        active_wait_normalized_bound=active_wait_bound,wait_groups=len(groups),unchained_waits=len(waits)-len(seen),
        nodes=len(graph.nodes),mandatory_nodes=len(mandatory),reachability_visits=visits,
        prediction_complete=False,scope='Conditional fluid-model continuation bounds. Unissued descendants supply the floor; full remaining-service potential and proven wait concurrency supply the ceiling. Feasible reachable execution is assumed. Not hardware-time or empirical uncertainty bounds.')


def source_wait_chains(graph, candidate, *, indexed):
    """Propose source worker chains; remaining_work_bounds proves their order."""
    # A device drains its scan before the next phenotype tile. Reuse its
    # worker slots across those tiles instead of pretending all tile waits
    # can coexist. The DAG proof still checks every proposed link.
    groups={}
    for index,tile in enumerate(candidate['tiles']):
        prefix=f'tile:{index}:' if indexed else f'tile{index}:'
        count=len(tile['data']['encoded']['chunks']);depth=tile['profile']['depth']
        device=tile['device']
        groups.setdefault((device,'release'),[]).extend(f'{prefix}release_wait:{i}' for i in range(count))
        for slot in range(depth):
            groups.setdefault((device,'finish',slot),[]).extend(f'{prefix}finish_wait:{i}' for i in range(slot,count,depth))
    return [[name for name in group if name in graph.resource_waits] for group in groups.values()
            if any(name in graph.resource_waits for name in group)]
