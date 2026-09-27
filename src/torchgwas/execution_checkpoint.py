"""Serializable state of the analytical scheduler, not measured executor state.

Remaining service is in nominal work-equivalent seconds. A loaded wall-clock
interval cannot be substituted for it without an explicit observation model.
"""
from copy import deepcopy
import math


SCHEMA='torchgwas.execution_checkpoint.v1'


def node_contract(graph,name):
    return dict(seconds=graph.nodes[name][0],after=list(graph.nodes[name][1]),
        demands=graph.demands.get(name,{}),tokens=graph.token_actions.get(name,{}),
        enqueue=list(graph.fifo_enqueues[name]) if name in graph.fifo_enqueues else None,
        dequeue=graph.fifo_dequeues.get(name),conditional=graph.conditional_delays.get(name))


def past_waits(graph,started):
    return {name:value for name,value in graph.resource_waits.items()
            if value['after'] in started or value['until'] in started}


def make_checkpoint(graph,now,starts,ends,remaining,ready,free,fifo,wait_usage,conditional,steps):
    return deepcopy(dict(schema=SCHEMA,at_seconds=now,start=starts,end=ends,
        remaining_service_seconds=remaining,ready=list(ready),token_free=free,
        fifo={key:list(values) for key,values in fifo.items() if values},
        wait_resource_seconds=dict(wait_usage),conditional_delays=conditional,
        resource_event_steps=steps,contracts={name:node_contract(graph,name) for name in starts},
        token_capacities=graph.token_capacities,past_resource_waits=past_waits(graph,starts)))


def _time(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
        raise ValueError('Invalid checkpoint '+name)
    return value


def validate_checkpoint(graph,state):
    """Preserve paid work, running services, queue order and acquired tokens.

    The caller may replace unstarted future work. It must separately preserve
    already-issued source ranges and allocation identities, which a generic
    execution DAG does not know. This validator is not memory admission.
    """
    fields={'schema','at_seconds','start','end','remaining_service_seconds','ready','token_free',
            'fifo','wait_resource_seconds','conditional_delays','resource_event_steps',
            'contracts','token_capacities','past_resource_waits'}
    if not isinstance(state,dict) or set(state)!=fields or state['schema']!=SCHEMA:
        raise ValueError('Unknown execution checkpoint schema')
    now=_time(state['at_seconds'],'time')
    starts,ends,remaining=(state[name] for name in ('start','end','remaining_service_seconds'))
    maps=('start','end','remaining_service_seconds','token_free','fifo','wait_resource_seconds',
          'conditional_delays','contracts','token_capacities','past_resource_waits')
    if any(not isinstance(state[name],dict) for name in maps):raise ValueError('Checkpoint maps required')
    if (set(starts)-graph.nodes.keys() or set(ends)-starts.keys() or set(remaining)!=set(starts)-ends.keys()
        or set(state['contracts'])!=set(starts)):
        raise ValueError('Checkpoint completed/running node coverage mismatch')
    ready=state['ready']
    if not isinstance(ready,list) or any(not isinstance(n,str) for n in ready) or len(ready)!=len(set(ready)) or set(ready)&starts.keys():
        raise ValueError('Invalid checkpoint ready order')
    if type(state['resource_event_steps']) is not int or state['resource_event_steps']<0:
        raise ValueError('Invalid checkpoint event count')
    for name,start in starts.items():
        _time(start,'start')
        if start>now or state['contracts'][name]!=node_contract(graph,name):
            raise ValueError('Started node contract changed: '+name)
        for dep in graph.nodes[name][1]:
            if dep not in ends or ends[dep]>start+1e-12:
                raise ValueError('Checkpoint starts before dependency completion')
        if name in ends:
            end=_time(ends[name],'end')
            if not start<=end<=now:raise ValueError('Invalid checkpoint completion order')
        else:
            amount=_time(remaining[name],'remaining service')
            if not 0<amount<=graph.nodes[name][0]:
                raise ValueError('Invalid remaining work for running node')
    if state['token_capacities']!=graph.token_capacities:
        raise ValueError('Checkpoint token capacities cannot change')
    expected=dict(graph.token_capacities)
    for name in starts:
        actions=graph.token_actions.get(name,{})
        for key,value in actions.get('acquire',{}).items():expected[key]-=value
        for key,value in actions.get('release_start',{}).items():expected[key]+=value
        if name in ends:
            for key,value in actions.get('release_finish',{}).items():expected[key]+=value
    if (state['token_free']!=expected or any(type(value) is not int or not 0<=value<=graph.token_capacities[key]
                                         for key,value in state['token_free'].items())):
        raise ValueError('Checkpoint token balance mismatch')
    expected_fifo={};enqueue_times={}
    for name,(queue,consumer) in graph.fifo_enqueues.items():
        if name not in ends or consumer in starts:continue
        if consumer not in graph.nodes or graph.fifo_dequeues.get(consumer)!=queue:
            raise ValueError('Queued checkpoint consumer changed or removed')
        if consumer in enqueue_times:raise ValueError('Duplicate FIFO consumer')
        expected_fifo.setdefault(queue,set()).add(consumer);enqueue_times[consumer]=ends[name]
    if set(state['fifo'])!=set(expected_fifo):raise ValueError('Checkpoint FIFO queues differ')
    for key,values in state['fifo'].items():
        if (not isinstance(values,list) or any(not isinstance(n,str) for n in values)
            or len(values)!=len(set(values)) or set(values)!=expected_fifo[key]):
            raise ValueError('Checkpoint FIFO contents differ')
        times=[enqueue_times[name] for name in values]
        if times!=sorted(times):raise ValueError('Checkpoint FIFO completion order differs')
    if state['past_resource_waits']!=past_waits(graph,starts):
        raise ValueError('Started resource-wait contract changed')
    for value in state['wait_resource_seconds'].values():_time(value,'resource work')
    if any(name not in starts or name not in graph.conditional_delays for name in state['conditional_delays']):
        raise ValueError('Unknown completed/running conditional delay')
    for name in starts:
        if name not in graph.conditional_delays:continue
        condition=graph.conditional_delays[name]
        waited=ends[condition['ready']]-ends[condition['attempt']]
        blocked=waited>1e-12
        expected=dict(blocked=blocked,dependency_wait_seconds=max(0.,waited),
                      extra_elapsed_service_seconds=graph.nodes[name][0] if blocked else 0.)
        if state['conditional_delays'].get(name)!=expected:
            raise ValueError('Checkpoint conditional-delay branch differs')
    return deepcopy(state)
