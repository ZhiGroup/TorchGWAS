"""Capacity accounting for an already solved finite execution graph."""
import math


def resource_balance(graph, solution):
    """Explain shared-resource pressure without changing the chosen schedule.

    Work/capacity is a necessary lower bound under the supplied resource model.
    Include simulated resource-consuming waits and only conditional service
    that actually ran. Dependency stalls use no extra work unless the graph
    explicitly represents them. This is not measured hardware utilization.
    """
    if set(graph.nodes) != set(solution['start']) or set(graph.nodes) != set(solution['end']):
        raise ValueError('Resource accounting requires the complete solved static graph')
    elapsed = solution['seconds']
    conditional = solution.get('conditional_delays', {})
    terms = {key: [] for key in graph.capacities}
    for name, (seconds, _) in graph.nodes.items():
        if name in graph.conditional_delays:
            seconds = conditional[name]['extra_elapsed_service_seconds']
        for key, rate in graph.demands.get(name, {}).items():
            terms[key].append(seconds * rate)
    for key, work in solution.get('wait_resource_seconds', {}).items():
        terms[key].append(work)
    resources = {}
    for key, amounts in terms.items():
        capacity = graph.capacities[key]
        work = math.fsum(amounts)
        seconds = work / capacity
        if seconds > elapsed + max(1e-10, 1e-9 * elapsed):
            raise ValueError('Solved resource work exceeds capacity over elapsed time: ' + key)
        resources[key] = dict(work=work, capacity=capacity,
            capacity_seconds=seconds, mean_capacity_fraction=seconds / elapsed if elapsed else 0.)
    lower = max((row['capacity_seconds'] for row in resources.values()), default=0.)
    limiting = [key for key, row in resources.items()
        if lower > 0 and math.isclose(row['capacity_seconds'], lower, rel_tol=1e-9, abs_tol=1e-12)]
    return dict(resources=resources, resource_lower_bound_seconds=lower,
        largest_resource_bounds=limiting, scheduled_seconds=elapsed,
        scope='Declared-resource work divided by supplied capacity, including simulated active waits. '
              'These are conditional model lower bounds and average capacity fractions, not measured '
              'utilization or a complete explanation of dependency, queue, or unpriced costs.')


def repeated_resource_floor(templates, *, max_templates=1024, max_template_nodes=100000):
    """Necessary work/capacity floor without expanding repeated source graphs.

    Each term is {'graph': ExecutionGraph, 'copies': positive integer}. Graphs
    must represent mandatory work at their declared extents and survivor
    scenario. Full chunks, tails and tile classes can share one template each.
    Optional handoff delays and active waiting are omitted, because their work
    depends on contention; a solved short-run wait cannot simply be multiplied.
    This floor can prune a candidate against a feasible schedule, but cannot
    rank feasible candidates or replace the pipeline dependency calculation.
    """
    from .execution_graph import ExecutionGraph
    from .mechanistic_plan import _integer
    _integer('max_templates',max_templates);_integer('max_template_nodes',max_template_nodes)
    if not isinstance(templates,(list,tuple)) or not templates or len(templates)>max_templates:
        raise ValueError('Nonempty bounded graph-template sequence required')
    capacities={};nodes=0
    for term in templates:
        if not isinstance(term,dict) or set(term)!={'graph','copies'} or not isinstance(term['graph'],ExecutionGraph):
            raise ValueError('Each term requires a source graph and integer copies')
        _integer('copies',term['copies']);g=term['graph'];nodes+=len(g.nodes)
        for key,value in g.capacities.items():
            if isinstance(value,bool) or not math.isfinite(value) or value<=0:
                raise ValueError('Positive finite resource capacity required')
            if key in capacities and capacities[key]!=value:raise ValueError('Conflicting shared resource capacity: '+key)
            capacities[key]=value
    if nodes>max_template_nodes:raise ValueError('Resource floor exceeds max_template_nodes')
    amounts={key:[] for key in capacities};critical=0.;reports=[]
    for term in templates:
        g,copies=term['graph'],term['copies'];path=ExecutionGraph();local={key:[] for key in capacities}
        if set(g.conditional_delays)-set(g.nodes):raise ValueError('Unknown conditional-delay node')
        if set(g.demands)-set(g.nodes):raise ValueError('Unknown resource-demand node')
        for name,(duration,deps) in g.nodes.items():
            if not math.isfinite(duration) or duration<0:raise ValueError('Invalid source service duration')
            seconds=0. if name in g.conditional_delays else duration
            stretch=1.
            for key,rate in g.demands.get(name,{}).items():
                if isinstance(rate,bool) or not math.isfinite(rate) or rate<0:
                    raise ValueError('Invalid resource demand')
                if rate and key not in capacities:raise ValueError('Missing resource capacity: '+key)
                if rate:
                    local[key].append(seconds*rate)
                    stretch=max(stretch,rate/capacities[key])
            path.add(name,seconds*stretch,deps)
        local_work={key:math.fsum(values) for key,values in local.items()}
        for key,value in local_work.items():
            total=value*copies
            if not math.isfinite(total):raise ValueError('Repeated resource work overflow')
            amounts[key].append(total)
        # Copies need not drain serially, so only one template critical path is
        # mandatory here. Cross-copy/token/FIFO ordering can only raise the floor.
        span=path.solve()['seconds'];critical=max(critical,span)
        reports.append(dict(copies=copies,nodes=len(g.nodes),mandatory_work=local_work,
            single_copy_critical_path_seconds=span,omitted_conditional_nodes=len(g.conditional_delays),
            omitted_active_waits=len(g.resource_waits)))
    resources={}
    for key,values in amounts.items():
        work=math.fsum(values)
        if not math.isfinite(work):raise ValueError('Repeated resource work overflow')
        resources[key]=dict(work=work,capacity=capacities[key],capacity_seconds=work/capacities[key])
    resource_floor=max((r['capacity_seconds'] for r in resources.values()),default=0.)
    return dict(resources=resources,resource_lower_bound_seconds=resource_floor,
        single_copy_dependency_lower_bound_seconds=critical,lower_bound_seconds=max(resource_floor,critical),
        largest_resource_bounds=[key for key,row in resources.items() if resource_floor>0 and
            math.isclose(row['capacity_seconds'],resource_floor,rel_tol=1e-9,abs_tol=1e-12)],
        template_count=len(templates),template_nodes=nodes,copies=sum(t['copies'] for t in templates),templates=reports,
        prediction_complete=False,selection_validated=False,
        scope='Necessary lower bound within declared service/capacity scenarios, using mandatory source work only. '
            'No graph replication, observed runtime, steady-state fit, full schedule, or hardware-time guarantee. '
            'Omitted waits, optional handoffs, queue/FIFO/token constraints, and cross-copy dependencies can raise runtime.')
