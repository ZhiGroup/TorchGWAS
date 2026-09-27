"""One source-calculator proposal from an actual reserved read prefix.

No candidate grid, association timing fit or invented live checkpoint is used.
Memory admission, price freshness and the cheap planning-cost gate belong to
the caller. This callback runs only after useful indexed output has begun.
"""
from .adaptive_candidate import future_chunk_candidate
from .adaptive_chunks import aligned_chunk_sizes
from .planning_session import planning_work_scope
from .remaining_work import remaining_work_bounds, source_wait_chains


def _reserved_prefix(candidate, snapshot):
    if (not isinstance(snapshot,dict) or snapshot.get('prefix_complete') is not True
            or snapshot.get('finished') is not None or snapshot.get('first_written') is None):
        raise ValueError('Complete active productive-run prefix after useful output required')
    rows=snapshot.get('partitions')
    if not isinstance(rows,list) or len(rows)!=len(candidate['tiles']):
        raise ValueError('One exact executor partition per candidate tile required')
    def identity(device,variants,traits):return (device,tuple(variants),tuple(traits))
    indexed={}
    for row in rows:
        key=identity(row['device'],row['variant_range'],row['trait_range'])
        if key in indexed:raise ValueError('Ambiguous executor partition')
        ranges=row['ranges'];count=row['issued_chunks']
        if (type(count) is not int or count<0 or not isinstance(ranges,list) or len(ranges)!=count
                or row['cursor']!=(ranges[-1][1] if ranges else row['variant_range'][0])):
            raise ValueError('Incomplete executor prefix or inconsistent cursor')
        indexed[key]=(count,ranges)
    counts=[];prefixes=[]
    for tile in candidate['tiles']:
        key=identity(tile['device'],tile['data']['encoded']['variant_range'],tile['trait_range'])
        if key not in indexed:raise ValueError('Executor partition differs from admitted candidate')
        count,ranges=indexed.pop(key);counts.append(count);prefixes.append(ranges)
    return counts,prefixes


def analytical_chunk_proposal(candidate, *, source_census, snapshot, chunk_sizes,
                              next_size, reduction, graph_factory,
                              max_source_chunks=10000, max_census_chunks=1000000,
                              max_nodes=1000000, max_reachability_visits=1000000):
    """Return a proposal and compact audit using a sufficient model condition.

    baseline_seconds is a necessary floor for the current choice's UNISSUED
    work. candidate_seconds is an upper bound for a feasible continuation of
    the alternative graph, conservatively counting its entire service prefix.
    A positive difference can repay planning even without live queue telemetry,
    conditional on the supplied graph model. Lack of a positive difference does
    not establish that the current size is best. Bounds may be too loose to act.
    """
    if not callable(graph_factory):raise ValueError('Source analytical graph factory required')
    sizes=aligned_chunk_sizes(chunk_sizes)
    current=snapshot.get('current_chunk_size')
    if type(current) is not int or current not in sizes:
        raise ValueError('Current size is outside the admitted choices')
    if type(next_size) is not int or next_size not in sizes:
        raise ValueError('Proposed size is outside the admitted choices')
    counts,prefixes=_reserved_prefix(candidate,snapshot)
    options=dict(source_census=source_census,chunk_sizes=sizes,issued_chunks=counts,
        issued_ranges=prefixes,reduction=reduction,max_source_chunks=max_source_chunks,
        max_census_chunks=max_census_chunks)
    baseline=future_chunk_candidate(candidate,next_size=current,**options)
    alternative=baseline if next_size==current else future_chunk_candidate(candidate,next_size=next_size,**options)
    indexed=reduction is not None
    def evaluate(choice):
        graph=graph_factory(choice)
        unissued=[]
        for index,(tile,count) in enumerate(zip(choice['tiles'],counts)):
            prefix=f'tile:{index}:' if indexed else f'tile{index}:'
            unissued.extend(f'{prefix}submit_decode:{i}' for i in range(count,len(tile['data']['encoded']['chunks'])))
        return remaining_work_bounds(graph,unissued_nodes=unissued,
            wait_chains=source_wait_chains(graph,choice,indexed=indexed),max_nodes=max_nodes,
            max_reachability_visits=max_reachability_visits)
    with planning_work_scope():
        before=evaluate(baseline)
        after=before if next_size==current else evaluate(alternative)
    return dict(chunk_size=next_size,baseline_seconds=before['lower_seconds'],
                candidate_seconds=after['upper_seconds']),dict(
        issued_revision=snapshot['issued_revision'],issued_chunks=counts,
        current_size=current,proposed_size=next_size,baseline=before,candidate=after,
        conditional_gain_floor_seconds=before['lower_seconds']-after['upper_seconds'],
        timing_policy='current unissued-work floor versus candidate full-service continuation ceiling',
        scope='One analytical sufficient-improvement proposal on fixed admitted partitions. Conditional component-service bounds, not hardware guarantees, a Bayesian policy, or a claim of optimality. No full candidate search or graph event simulation.')
