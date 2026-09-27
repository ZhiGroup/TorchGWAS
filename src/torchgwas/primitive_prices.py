"""CPU service statistics for independent component experiments."""
import math
import statistics


def cpu_service_samples(values, timer_samples):
    """Keep expensive calls in additive CPU work; median is descriptive only.

    Subtract one independently measured empty-timer interval per observation.
    This does not remove all probe perturbation or estimate another context.
    """
    values, timer_samples = list(values), list(timer_samples)
    if not values or not timer_samples:
        raise ValueError('Nonempty call and timer samples required')
    for sample in values + timer_samples:
        if isinstance(sample, bool) or not isinstance(sample, (float, int)) or not math.isfinite(sample) or sample < 0:
            raise ValueError('Finite nonnegative CPU samples required')
    raw = math.fsum(values)
    overhead = statistics.fmean(timer_samples)
    return dict(cpu_seconds_per_call=max(0., raw / len(values) - overhead),
                median_cpu_seconds_per_call=max(0., statistics.median(values) - statistics.median(timer_samples)),
                raw_cpu_seconds_per_call=raw / len(values),
                raw_cpu_sum_seconds=raw, sample_count=len(values),
                timer_overhead_seconds=overhead,
                cpu_cost_statistic='arithmetic mean within repeat, empty timer mean removed')


def mean_cpu_service_prices(rows, *, device, minimum_repeats=3):
    """Validate an additive service ledger before taking median repeat means."""
    selected = [row for row in rows if row.get('device') == device]
    groups = {}
    for row in selected:
        name, repeat = row.get('primitive'), row.get('repeat')
        if not isinstance(name, str) or not name or isinstance(repeat, bool) or not isinstance(repeat, int) or repeat < 0:
            raise ValueError('Primitive and repeat identifiers required')
        count = row.get('sample_count')
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError('Raw service sample count required')
        names = ['cpu_seconds_per_call', 'raw_cpu_seconds_per_call', 'raw_cpu_sum_seconds', 'timer_overhead_seconds']
        values = [row.get(key) for key in names]
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in values):
            raise ValueError('Finite nonnegative service ledger required')
        price, raw_mean, raw_sum, overhead = values
        if not math.isclose(raw_mean * count, raw_sum, rel_tol=1e-10, abs_tol=1e-15):
            raise ValueError('Raw mean must conserve observed CPU sum')
        if not math.isclose(price, max(0., raw_mean - overhead), rel_tol=1e-10, abs_tol=1e-15):
            raise ValueError('CPU price must equal timer-corrected arithmetic mean')
        group = groups.setdefault(name, {})
        if repeat in group:
            raise ValueError('Duplicate primitive repeat')
        group[repeat] = price
    if not groups or any(len(group) < minimum_repeats for group in groups.values()):
        raise ValueError('Insufficient independent service repeats')
    if len({tuple(sorted(group)) for group in groups.values()}) != 1:
        raise ValueError('Primitive repeat coverage must agree')
    return {name: statistics.median(group.values()) for name, group in groups.items()}
