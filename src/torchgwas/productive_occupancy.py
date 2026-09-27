"""Bounded reduced-output scenarios on explicit association grids."""
from copy import deepcopy
from math import gcd


def occupancy_kind(value):
    if isinstance(value, str):
        if value in ('empty', 'dense'):
            return value
        raise ValueError('Reduced output requires empty, dense or explicit sparse occupancy')
    if not isinstance(value, dict) or set(value) != {'retained_fraction', 'placement'}:
        raise ValueError('Reduced output requires empty, dense or explicit sparse occupancy')
    fraction = value['retained_fraction']
    placement = value['placement']
    # Spread uses a logarithmic floor-sum. Clustered placement can require
    # visiting one full coordinate period, so keep its tighter bound.
    limit = 1_000_000_000_000 if placement == 'spread' else 1_000_000
    if (not isinstance(fraction, (list, tuple)) or len(fraction) != 2 or
            any(type(part) is not int for part in fraction) or
            not 0 < fraction[0] < fraction[1] <= limit or
            placement not in ('spread', 'clustered')):
        raise ValueError('Sparse occupancy requires a bounded rational fraction and placement')
    return 'sparse'


def validate_occupancy_scenarios(scenarios):
    if (not isinstance(scenarios, dict) or not scenarios or
            any(not isinstance(name, str) or not name.strip()
                for name in scenarios)):
        raise ValueError('Named reduced-output occupancy scenarios required')
    for value in scenarios.values():
        occupancy_kind(value)


def survivor_bins(partitions, horizon, *, fine_chunk, reduction,
                  total_traits, scenario, max_bins=256,
                  max_cluster_visits=1_000_000):
    """Keep partial counts whole when either chunk size rebins the same work.

    Sparse counts are declared conditional scenarios on the same global
    association ledger as whole-layout continuation. They are not inferred
    from the threshold or first chunks.
    """
    kind = occupancy_kind(scenario)
    if (reduction not in ('significant', 'jagwas') or
            type(total_traits) is not int or total_traits < 1 or
            type(fine_chunk) is not int or fine_chunk < 1 or
            type(horizon) is not int or horizon < 1 or horizon % fine_chunk or
            type(max_bins) is not int or max_bins < 1 or
            type(max_cluster_visits) is not int or max_cluster_visits < 1 or
            not isinstance(partitions, list) or not partitions):
        raise ValueError('Bounded reduced-output source geometry required')
    count = len(partitions) * (horizon // fine_chunk if kind == 'sparse' else 1)
    if count > max_bins:
        raise ValueError('Sparse survivor evidence exceeds the productive bin budget')
    bins = []
    for part in partitions:
        if not isinstance(part, dict):
            raise ValueError('Explicit productive partition required')
        lo = part.get('cursor')
        end = part.get('variant_range')
        trait = part.get('trait_range')
        if (type(lo) is not int or not isinstance(end, (list, tuple)) or
                len(end) != 2 or any(type(x) is not int for x in end) or
                not end[0] <= lo < lo + horizon <= end[1] or
                not isinstance(trait, (list, tuple)) or len(trait) != 2 or
                any(type(x) is not int for x in trait) or
                not 0 <= trait[0] < trait[1] <= total_traits or
                (reduction == 'jagwas' and tuple(trait) != (0, total_traits))):
            raise ValueError('Reduced-output scenario differs from unissued partition')
        step = fine_chunk if kind == 'sparse' else horizon
        possible_per_marker = trait[1] - trait[0] if reduction == 'significant' else 1
        for offset in range(0, horizon, step):
            cells = step * possible_per_marker
            retained = 0 if kind in ('empty', 'sparse') else cells
            bins.append(dict(variant_range=[lo + offset, lo + offset + step],
                             trait_range=list(trait), retained=retained))
    if kind == 'sparse':
        layout = dict(kind='torchgwas.pgen_layout_source_floor.v1',
                      reduction=reduction, total_traits=total_traits,
                      partitions=[
                          dict(id=str(index),
                               variant_range=row['variant_range'],
                               trait_range=row['trait_range'])
                          for index, row in enumerate(bins)])
        counts = whole_layout_survivors(
            layout, scenario, max_partitions=max_bins,
            max_cluster_visits=max_cluster_visits)['retained_ranges']
        for index, row in enumerate(bins):
            row['retained'] = counts[str(index)][0]
    return bins


def _floor_sum(count, modulus, slope, offset):
    """Sum floor((slope*i+offset)/modulus), i in [0,count), in log time."""
    total = 0
    while True:
        quotient, slope = divmod(slope, modulus)
        total += count * (count - 1) // 2 * quotient
        quotient, offset = divmod(offset, modulus)
        total += count * quotient
        upper = slope * count + offset
        if upper < modulus:
            return total
        count, offset = divmod(upper, modulus)
        modulus, slope = slope, modulus


def whole_layout_survivors(layout, scenario, *, max_partitions=10_000,
                           max_cluster_visits=1_000_000):
    """Count one declared sparse pattern on a global association coordinate.

    A significant pair has coordinate variant*total_traits+trait; a JAGWAS
    retained row has coordinate variant. Thus candidate tiles and GPU shards
    see the same pattern even when ownership changes. This is conditional
    output work, not an estimate of future selectivity from observed chunks.
    """
    kind = occupancy_kind(scenario)
    if (not isinstance(layout, dict) or
            layout.get('kind') != 'torchgwas.pgen_layout_source_floor.v1' or
            layout.get('reduction') not in ('significant', 'jagwas') or
            type(layout.get('total_traits')) is not int or
            layout['total_traits'] < 1 or
            not isinstance(layout.get('partitions'), list) or
            type(max_partitions) is not int or max_partitions < 1 or
            not 0 < len(layout['partitions']) <= max_partitions or
            type(max_cluster_visits) is not int or max_cluster_visits < 1):
        raise ValueError('Bounded reduced-output layout required')
    reduction = layout['reduction']
    traits = layout['total_traits']
    if kind == 'sparse':
        numerator, denominator = scenario['retained_fraction']
        placement = scenario['placement']

        def cumulative(cells):
            if placement == 'spread':
                return (cells * numerator + denominator - 1) // denominator
            periods, remainder = divmod(cells, denominator)
            return periods * numerator + min(remainder, numerator)

    retained = {}
    possible_total = retained_total = visits = 0
    for row in layout['partitions']:
        if not isinstance(row, dict) or not isinstance(row.get('id'), str) or not row['id']:
            raise ValueError('Explicit reduced-output partition required')
        key = row['id']
        if key in retained:
            raise ValueError('Unique reduced-output partition ids required')
        variant = row.get('variant_range')
        trait = row.get('trait_range')
        if (not isinstance(variant, (list, tuple)) or len(variant) != 2 or
                not isinstance(trait, (list, tuple)) or len(trait) != 2 or
                any(type(value) is not int for value in (*variant, *trait)) or
                not 0 <= variant[0] < variant[1] or
                not 0 <= trait[0] < trait[1] <= traits):
            raise ValueError('Invalid reduced-output partition geometry')
        if reduction == 'jagwas' and list(trait) != [0, traits]:
            raise ValueError('JAGWAS requires the complete phenotype panel')
        lo, hi = variant
        a, b = trait
        possible = (hi - lo) * (b - a if reduction == 'significant' else 1)
        if kind == 'empty':
            count = 0
        elif kind == 'dense':
            count = possible
        elif reduction == 'jagwas':
            count = cumulative(hi) - cumulative(lo)
        elif placement == 'spread':
            # F(x)=ceil(x*n/d). Sum F(v*K+b)-F(v*K+a) without
            # visiting the potentially millions of variants.
            span = hi - lo
            slope = traits * numerator
            upper = _floor_sum(span, denominator, slope,
                               (lo * traits + b) * numerator + denominator - 1)
            lower = _floor_sum(span, denominator, slope,
                               (lo * traits + a) * numerator + denominator - 1)
            count = upper - lower
        else:
            period = denominator // gcd(traits, denominator)
            if visits + min(period, hi - lo) > max_cluster_visits:
                raise ValueError('Clustered occupancy exceeds the whole-layout visit budget')
            cycles, tail = divmod(hi - lo, period)
            period_total = tail_total = 0
            for index in range(min(period, hi - lo)):
                marker = lo + index
                selected = (cumulative(marker * traits + b) -
                            cumulative(marker * traits + a))
                period_total += selected
                if index < tail:
                    tail_total += selected
            visits += min(period, hi - lo)
            count = cycles * period_total + tail_total
        if not 0 <= count <= possible:
            raise ValueError('Reduced-output scenario exceeds association capacity')
        retained[key] = [count, count]
        retained_total += count
        possible_total += possible
    return dict(kind='torchgwas.whole_layout_survivors.v1',
                input_identity=deepcopy(layout.get('input_identity')),
                reduction=reduction, scenario=deepcopy(scenario),
                retained_ranges=retained, possible=possible_total,
                retained=retained_total, clustered_variant_visits=visits,
                scope='Exact conditional retained counts on a global pair or variant coordinate, invariant to candidate tile/GPU ownership. No future-selectivity inference, runtime service or switch authorization.')
