"""Validate independent GIL-meter evidence and expose conditional API service."""
import math
import statistics


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('Finite nonnegative ' + name + ' required')
    return value


def gil_service_prices(probe, *, mode, device, torch_version, cpu_affinity):
    if probe.get('context_verified') is not True or probe.get('torch_version') != torch_version or probe.get('affinity') != list(cpu_affinity):
        raise ValueError('GIL probe runtime, affinity or context mismatch')
    workers = [worker for result in probe['results'] if result['mode'] == mode
               for worker in result['workers'] if worker['device'] == device]
    if len(workers) != 1 or workers[0].get('context_verified') is not True:
        raise ValueError('One verified GIL worker context required')
    worker = workers[0]
    for name, intervals in [('native',1), ('nested',3), ('python',0)]:
        control = worker['controls'][name]
        cpu = _number(control['cpu_seconds'], 'control CPU')
        detached = _number(control['detached_cpu_seconds'], 'detached control CPU')
        if control['balance_errors'] != 0 or control['detached_intervals'] != intervals or cpu <= 0 or detached > cpu:
            raise ValueError('Invalid GIL control')
        if (intervals and detached / cpu <= .9) or (not intervals and detached != 0):
            raise ValueError('GIL positive/negative control failed')
    meters = {}
    for row in worker['meter_rows']:
        repeat, recording = row['repeat'], row['recording']
        if isinstance(repeat,bool) or not isinstance(repeat,int) or repeat<0 or not isinstance(recording,bool):
            raise ValueError('Invalid meter repeat')
        if row['intervals'] != 20000 or row['detached_intervals'] != (20000 if recording else 0) or row['balance_errors'] != 0:
            raise ValueError('Meter-only interval accounting failed')
        cpu = _number(row['cpu_seconds'], 'meter CPU')
        detached = _number(row['detached_cpu_seconds'], 'detached meter CPU')
        if detached > cpu or (not recording and detached):raise ValueError('Invalid detached meter CPU')
        pair = meters.setdefault(repeat,{})
        if recording in pair:raise ValueError('Duplicate meter repeat')
        pair[recording] = (cpu,detached)
    if len(meters)<3 or any(set(pair)!={False,True} for pair in meters.values()):
        raise ValueError('Insufficient paired meter repeats')
    overhead = statistics.median((pair[True][0]-pair[False][0])/20000 for pair in meters.values())
    detached_overhead = statistics.median(pair[True][1]/20000 for pair in meters.values())
    if not 0 <= detached_overhead <= overhead:raise ValueError('Inconsistent meter overhead partition')
    pairs = {}
    for row in worker['rows']:
        name,repeat,recording = row['primitive'],row['repeat'],row['recording']
        if not isinstance(name,str) or not name or isinstance(repeat,bool) or not isinstance(repeat,int) or repeat<0 or not isinstance(recording,bool):
            raise ValueError('Invalid API service context')
        count = row['sample_count']
        if isinstance(count,bool) or not isinstance(count,int) or count<1:raise ValueError('API sample count required')
        cpu,raw,total,empty,detached = [_number(row[key],key) for key in
            ['cpu_seconds_per_call','raw_cpu_seconds_per_call','raw_cpu_sum_seconds','empty_cpu_seconds_per_call','detached_cpu_seconds_per_call']]
        if not math.isclose(raw*count,total,rel_tol=1e-10,abs_tol=1e-15) or not math.isclose(cpu,max(0.,raw-empty),rel_tol=1e-10,abs_tol=1e-15):
            raise ValueError('API CPU mean does not conserve raw observations')
        intervals = row['detached_intervals']
        if isinstance(intervals,bool) or not isinstance(intervals,int) or intervals<0 or row['balance_errors']!=0:
            raise ValueError('Unbalanced API GIL transitions')
        if not recording and (intervals or detached):raise ValueError('Disabled meter reports detached service')
        pair = pairs.setdefault(name,{}).setdefault(repeat,{})
        if recording in pair:raise ValueError('Duplicate API repeat')
        pair[recording] = (cpu,detached,intervals/count)
    if not pairs or any(len(repeats)<3 or any(set(pair)!={False,True} for pair in repeats.values()) for repeats in pairs.values()):
        raise ValueError('Insufficient paired API repeats')
    if len({tuple(sorted(repeats)) for repeats in pairs.values()})!=1:raise ValueError('API repeat coverage differs')
    cpu_prices = {};serial_prices = {};comparisons = {}
    for name,repeats in pairs.items():
        costs = [];held = [];deltas = []
        for pair in repeats.values():
            cpu,detached,intervals = pair[True]
            corrected = cpu-intervals*overhead
            detached -= intervals*detached_overhead
            if not 0 <= detached <= corrected:raise ValueError('Meter correction exceeds observed API service: '+name)
            costs.append(corrected);held.append(corrected-detached)
            if pair[False][0]>0:deltas.append(corrected/pair[False][0]-1)
        cpu_prices[name] = statistics.median(costs)
        serial_prices[name] = statistics.median(held)
        comparisons[name] = dict(min_relative_delta=min(deltas,default=0.),max_relative_delta=max(deltas,default=0.))
    return dict(cpu_primitives=cpu_prices,serial_primitives=serial_prices,
        meter_cpu_seconds_per_interval=overhead,meter_detached_cpu_seconds_per_interval=detached_overhead,
        paired_recording_deltas=comparisons,
        measurement_controls_passed=True,total_cpu_transfer_qualified=False,
        serial_transfer_qualified=False,prediction_complete=False,
        scope='Independent fixed-API CPU and held-or-unknown CPU, meter-only overhead removed per observed interval. Median repeat means; unobserved paths remain held/unknown. Probe perturbation and context transfer remain assumptions, not certified bounds.')
