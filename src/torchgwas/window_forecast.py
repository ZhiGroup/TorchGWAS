"""Explicit remaining-work scenarios from bounded analytical window horizons.

Slopes come from the source/resource model, never measured GWAS timing fits.
Forecast intervals express supplied scenarios and tolerance, not confidence
intervals, live checkpoints, or hardware-time guarantees.
"""
from copy import deepcopy
import math

from .window_model import _coverage,_integer,_number


def _coverage_units(coverage):
    if not isinstance(coverage,list) or not 0<len(coverage)<=16 or any(
            not isinstance(row['trait_ranges'],list) or not 0<len(row['trait_ranges'])<=8 for row in coverage):
        raise ValueError('Bounded comparison coverage required')
    return sum((row['variant_range'][1]-row['variant_range'][0])*
        sum(hi-lo for lo,hi in row['trait_ranges']) for row in coverage)


def _interval(values,name):
    if not isinstance(values,(list,tuple)) or len(values)!=2:
        raise ValueError('Explicit lower/upper '+name+' required')
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in values) or values[1]<values[0]:
        raise ValueError('Invalid '+name+' interval')
    return list(values)


def forecast_remaining_windows(comparisons, *, remaining_pairs, boundary_adjustments,
        relative_model_error, max_slope_change, max_extrapolation, assumptions,
        remaining_partition_spans=None):
    """Forecast one unchanged pair of layouts from three nested model horizons.

    Work units are variant-phenotype pairs, including for full-panel JAGWAS.
    Every window must grow proportionally across horizons, holding its device,
    trait range, starting variant, chunk size, service prices and output fixed.
    Remaining-work balance/source/selectivity assumptions are explicit caller
    inputs. Boundary adjustments cover in-flight work, existing writer state,
    preparation, tails and final publication absent from prepared windows.

    The caller must separately prove remaining_pairs describes the actual
    unissued workload and that live memory/parameter bindings are still valid.
    Optional exact remaining marker counts per layout/window produce balanced
    core/envelope scenarios when devices have unequal unissued extents. These
    are workload scenarios, not a proof that fluid graph time is monotone.
    No forecast returned here authorizes a public configuration change.
    """
    _integer(remaining_pairs,'remaining_pairs')
    for value,name in [(relative_model_error,'relative_model_error'),(max_slope_change,'max_slope_change')]:
        if _number(value,name,True)>=1:raise ValueError(name+' must be below one')
    if _number(max_extrapolation,'max_extrapolation')<1:raise ValueError('max_extrapolation must be at least one')
    if (not isinstance(assumptions,dict) or set(assumptions)!=
            {'source_work','output_occupancy','resource_capacity','partition_balance'} or
            any(not isinstance(v,str) or not v.strip() for v in assumptions.values())):
        raise ValueError('Explicit source/output/resource/partition forecast assumptions required')
    if not isinstance(boundary_adjustments,dict) or set(boundary_adjustments)!={'baseline','candidate'}:
        raise ValueError('Both continuation boundary adjustments are required')
    adjustments={name:_interval(v,name) for name,v in boundary_adjustments.items()}
    if not isinstance(comparisons,(list,tuple)) or len(comparisons)!=3:
        raise ValueError('Exactly three bounded analytical horizons required')
    contract=comparisons[0].get('comparison_contract')
    if not isinstance(contract,dict) or not isinstance(contract.get('model_identity'),dict) or not contract['model_identity']:
        raise ValueError('Caller-validated source/runtime/model identity required')
    if contract.get('schema')!='torchgwas.prepared_comparison.v1':raise ValueError('Unknown comparison binding')
    if contract['reduction'] is None and contract['output']['block_bytes'] is None:
        raise ValueError('Dense continuation requires explicit existing writer block_bytes')
    units=[];seconds={name:[] for name in ('baseline','candidate')};shapes={name:[] for name in seconds}
    unpriced=set()
    for report in comparisons:
        if report.get('comparison_contract')!=contract:raise ValueError('Model, prices, output or layout changed across horizons')
        work=_coverage_units(report['coverage']);_integer(work,'window coverage')
        if units and work<=units[-1]:raise ValueError('Horizons must strictly increase')
        for name in seconds:
            prediction=report[name];windows=prediction['windows']
            if not isinstance(windows,list) or not 0<len(windows)<=8:raise ValueError('Bounded window reports required')
            rects=[dict(trait_range=w['trait_range'],data=dict(encoded=dict(variant_range=w['variant_range']))) for w in windows]
            if _coverage(rects)!=report['coverage']:raise ValueError('Window report coverage differs from its comparison')
            spans=[]
            for w,bound in zip(windows,contract[name]['windows']):
                if (w['device']!=bound['device'] or w['trait_range']!=bound['trait_range'] or
                        w['variant_range'][0]!=bound['variant_start'] or w['issued_chunks']!=bound['issued_chunks']):
                    raise ValueError('Window differs from its bound layout')
                spans.append(w['variant_range'][1]-w['variant_range'][0])
            if len(windows)!=len(contract[name]['windows']):raise ValueError('Layout window count changed')
            if shapes[name] and any(a*work!=b*units[-1] for a,b in zip(shapes[name][-1],spans)):
                raise ValueError('Window growth must preserve partition balance')
            shapes[name].append(spans)
            seconds[name].append(_number(prediction['estimated_window_seconds'],'modeled window seconds',True))
            unpriced.update(prediction['unpriced_terms'])
        units.append(work)
    if remaining_pairs<units[-1]:raise ValueError('Remaining work is below the modeled horizon')
    try:ratio=remaining_pairs/units[-1]
    except OverflowError:raise ValueError('Remaining forecast extent overflow') from None
    extents={name:dict(lower_pairs=remaining_pairs,upper_pairs=remaining_pairs,
        lower_ratio=ratio,upper_ratio=ratio) for name in seconds}
    if remaining_partition_spans is not None:
        if not isinstance(remaining_partition_spans,dict) or set(remaining_partition_spans)!=set(seconds):
            raise ValueError('Both remaining partition-span vectors required')
        for name,spans in remaining_partition_spans.items():
            if not isinstance(spans,(list,tuple)) or len(spans)!=len(shapes[name][-1]):
                raise ValueError('One remaining span per bound window required')
            for span,modeled in zip(spans,shapes[name][-1]):
                if _integer(span,'remaining partition span')<modeled:
                    raise ValueError('A modeled horizon exceeds its remaining partition')
            actual=sum(span*(w['trait_range'][1]-w['trait_range'][0])
                for span,w in zip(spans,contract[name]['windows']))
            if actual!=remaining_pairs:raise ValueError('Remaining partition spans differ from exact pair total')
            # Integer arithmetic rounds the balanced core down and envelope up.
            # Counts, rather than floating ratios, preserve huge-job coverage.
            low=min(units[-1]*span//modeled for span,modeled in zip(spans,shapes[name][-1]))
            high=max((units[-1]*span+modeled-1)//modeled for span,modeled in zip(spans,shapes[name][-1]))
            try:extents[name]=dict(lower_pairs=low,upper_pairs=high,
                lower_ratio=low/units[-1],upper_ratio=high/units[-1])
            except OverflowError:raise ValueError('Remaining forecast extent overflow') from None
    if any(row['upper_ratio']>max_extrapolation for row in extents.values()):
        raise ValueError('Remaining forecast exceeds the declared extrapolation limit')
    forecasts={};stable=True;regime_warnings=[]
    if contract['reduction'] is None:
        for name in ('baseline','candidate'):
            for index,w in enumerate(comparisons[-1][name]['windows']):
                for stream,settings in w['writer_streams'].items():
                    for previous in comparisons[:-1]:
                        old=previous[name]['windows'][index]['writer_streams'][stream]
                        if any(old[key]!=settings[key] for key in ('block_bytes','writeback_interval_bytes')):
                            raise ValueError('Writer configuration changed across horizons')
                    period=settings['writeback_interval_bytes'];calls=settings['writeback_submit_calls']
                    projected=settings['payload_bytes']*extents[name]['upper_pairs']//units[-1]
                    # Two submissions are the first horizon to include both a
                    # periodic submit and the previous-range wait/drop cycle.
                    if period and calls<2 and projected//period>calls:
                        regime_warnings.append(dict(layout=name,window=index,stream=stream,
                            modeled_writeback_submissions=calls,projected_writeback_submissions=projected//period))
    for name,values in seconds.items():
        slopes=[(values[i+1]-values[i])/(units[i+1]-units[i]) for i in (0,1)]
        if min(slopes)<=0:raise ValueError('Positive marginal service is required for extrapolation')
        change=abs(slopes[1]-slopes[0])/max(slopes)
        steady=change<=max_slope_change;stable &= steady
        extent=extents[name]
        low=values[-1]+(extent['lower_pairs']-units[-1])*min(slopes)
        high=values[-1]+(extent['upper_pairs']-units[-1])*max(slopes)
        low=max(0.,low*(1-relative_model_error)+adjustments[name][0])
        high=max(0.,high*(1+relative_model_error)+adjustments[name][1])
        if not math.isfinite(low) or not math.isfinite(high):raise ValueError('Remaining forecast overflow')
        forecasts[name]=dict(lower_seconds=low,upper_seconds=high,
            modeled_window_seconds=values,marginal_seconds_per_pair=slopes,
            marginal_relative_change=change,stable=steady,boundary_adjustment_seconds=adjustments[name])
    gain=forecasts['baseline']['lower_seconds']-forecasts['candidate']['upper_seconds']
    status='unstable_marginal_cost' if not stable else 'unmodeled_writer_regime' if regime_warnings else 'stable_scenario'
    return dict(status=status,writer_regime_warnings=regime_warnings,
        remaining_pairs=remaining_pairs,horizon_pairs=units,extrapolation_ratio=ratio,
        partition_extent_scenarios=extents,
        remaining_partition_spans=deepcopy(remaining_partition_spans),
        forecasts=forecasts,scenario_gain_floor_seconds=gain,
        calculation_cpu_seconds=sum(_number(r['calculation_cpu_seconds'],'calculation CPU',True) for r in comparisons),
        calculation_wall_seconds=sum(_number(r['calculation_wall_seconds'],'calculation wall',True) for r in comparisons),
        relative_model_error=relative_model_error,max_slope_change=max_slope_change,max_extrapolation=max_extrapolation,
        assumptions=deepcopy(assumptions),comparison_contract=deepcopy(contract),unpriced_terms=sorted(unpriced),
        selection_validated=False,
        scope='Three-horizon extrapolation of analytical model service under declared source, output, resource and balance scenarios. Boundary adjustments are caller supplied. Not measured runtime fitting, statistical confidence, observed continuation state or automatic switch authorization. Calculation costs exclude caller binding/header work and this forecast arithmetic.')


def remaining_forecast_payback(forecast, *, planning_seconds, switching_seconds,
        publication_seconds, reserve_seconds):
    """A forecast must pay for the entire declared tuning cost, including drain.

    This is a numerical decision prerequisite. It does not validate the model,
    future workload assumptions, memory admission or observation freshness.
    """
    costs={name:_number(value,name,True) for name,value in
        [('planning',planning_seconds),('switching',switching_seconds),
         ('publication',publication_seconds),('reserve',reserve_seconds)]}
    try:total=math.fsum(costs.values())
    except OverflowError:raise ValueError('Tuning cost overflow') from None
    if not math.isfinite(total):raise ValueError('Tuning cost overflow')
    gain=forecast['scenario_gain_floor_seconds']
    if isinstance(gain,bool) or not isinstance(gain,(int,float)) or not math.isfinite(gain):raise ValueError('Finite scenario gain required')
    stable=forecast['status']=='stable_scenario'
    return dict(worthwhile=stable and gain>total,scenario_gain_floor_seconds=gain,
        costs_seconds=costs,total_cost_seconds=total,net_scenario_gain_seconds=gain-total,
        reason=forecast['status'] if not stable else 'positive_payback' if gain>total else 'gain_does_not_repay_tuning',
        selection_validated=False)
