"""Explicit resource scenarios from independent primitive measurements.

Observed probe extrema are scenarios, not confidence intervals or guaranteed
runtime bounds. Association observations are never accepted here.
"""
from collections import defaultdict
import copy
import math
import statistics


def primitive_cpu_scenarios(rows):
    samples=defaultdict(list)
    for row in rows:
        name=row['primitive'];value=row['cpu_seconds_per_call']
        if not isinstance(name,str) or not name or not math.isfinite(value) or value<0:
            raise ValueError('Invalid primitive CPU observation')
        samples[name].append(value)
    if not samples:raise ValueError('Independent primitive measurements required')
    return {
        'nominal':{name:statistics.median(values) for name,values in samples.items()},
        'low_dispatch':{name:min(values) for name,values in samples.items()},
        'high_dispatch':{name:max(values) for name,values in samples.items()},
    }


def torch_runtime_scenarios(data,profile,host_scenarios):
    from .mechanistic_torch import torch_runtime
    if not host_scenarios or 'nominal' not in host_scenarios:
        raise ValueError('Include the nominal resource scenario')
    estimates={}
    for name,prices in host_scenarios.items():
        candidate=copy.deepcopy(profile);candidate['host_primitives']=dict(prices)
        result=torch_runtime(data,candidate)
        if result.get('estimated_seconds') is None:
            raise ValueError('A resource scenario has no finite estimate: '+name)
        estimates[name]=result
    worst=max(estimates,key=lambda name:estimates[name]['estimated_seconds'])
    best=min(estimates,key=lambda name:estimates[name]['estimated_seconds'])
    return dict(status='development_resource_scenarios',prediction_complete=False,
        nominal_seconds=estimates['nominal']['estimated_seconds'],
        scenario_min_seconds=estimates[best]['estimated_seconds'],
        scenario_max_seconds=estimates[worst]['estimated_seconds'],
        worst_scenario=worst,estimates=estimates,
        unpriced_terms=sorted({term for result in estimates.values() for term in result.get('unpriced_terms',[])}),
        scope='Conditional model range across supplied independent CPU dispatch scenarios. '
              'Not a confidence interval, validated runtime bound or complete uncertainty model.')


def compare_runtime_scenarios(candidates):
    """Compare candidates under *matched* resource scenarios, never fitted runs.

    Regret is conditional on the supplied finite scenario set. Unpriced work
    remains unpriced even if one candidate wins every modeled scenario.
    """
    if not candidates:
        raise ValueError('At least one candidate is required')
    scenarios = None
    times = {}
    unpriced = set()
    complete = True
    for name, report in candidates.items():
        estimates = report['estimates']
        keys = set(estimates)
        if not keys or (scenarios is not None and keys != scenarios):
            raise ValueError('Candidates require identical nonempty resource scenarios')
        scenarios = keys
        times[name] = {}
        complete = complete and report.get('prediction_complete', False)
        unpriced.update(report.get('unpriced_terms', []))
        for scenario, estimate in estimates.items():
            value = estimate.get('estimated_seconds')
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError('Candidate scenario time must be finite and positive')
            times[name][scenario] = value
            complete = complete and estimate.get('prediction_complete', False)
            unpriced.update(estimate.get('unpriced_terms', []))
    scenario_names = sorted(scenarios)
    best = {s: min(row[s] for row in times.values()) for s in scenario_names}
    winners = {s: [name for name, row in times.items() if row[s] == best[s]]
               for s in scenario_names}
    regret = {name: {s: row[s] / best[s] - 1 for s in scenario_names}
              for name, row in times.items()}
    worst_regret = {name: max(row.values()) for name, row in regret.items()}
    worst_time = {name: max(row.values()) for name, row in times.items()}
    minimax_regret = min(worst_regret.values())
    minimax_time = min(worst_time.values())
    return dict(
        status='conditional_candidate_comparison',
        prediction_complete=bool(complete and not unpriced),
        automatic_selection_ready=False,
        scenario_winners=winners,
        common_winners=[name for name in times if all(name in winners[s] for s in scenario_names)],
        regret_by_scenario=regret,
        worst_scenario_regret=worst_regret,
        minimax_regret_candidates=[name for name in times if worst_regret[name] == minimax_regret],
        minimax_time_candidates=[name for name in times if worst_time[name] == minimax_time],
        unpriced_terms=sorted(unpriced),
        scope='Matched supplied scenarios only; regret is not an empirical or guaranteed bound. '
              'A common winner does not validate missing work, future resources or automatic selection.')
