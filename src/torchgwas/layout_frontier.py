"""Exact association coverage at a productive JIT issue frontier.

The executor may have many prefetched or in-flight chunks. Only each
partition's contiguous *unissued* suffix may be reassigned. This module
proves rectangle equality for a proposed trait/variant/GPU layout; it does
not transfer an in-flight buffer or change writer ownership.
"""
from copy import deepcopy

from .analytical_plan_cache import input_identity


KIND = 'torchgwas.unissued_association_frontier.v1'


def _span(value, name):
    if (not isinstance(value, (list, tuple)) or len(value) != 2 or
            any(type(x) is not int for x in value) or
            not 0 <= value[0] < value[1]):
        raise ValueError('Nonempty ' + name + ' range required')
    return tuple(value)


def _rectangles_match(required, candidate, max_visits):
    """Sweep variant boundaries and compare phenotype coverage on each strip."""
    events = {}
    for side, rows in enumerate((required, candidate)):
        for index, row in enumerate(rows):
            lo, hi = _span(row['variant_range'], 'variant')
            trait = _span(row['trait_range'], 'phenotype')
            events.setdefault(lo, []).append((side, index, True, trait))
            events.setdefault(hi, []).append((side, index, False, trait))
    boundaries = sorted(events)
    active = ({}, {})
    visits = strips = 0
    for left, right in zip(boundaries, boundaries[1:]):
        for side, index, entering, trait in events[left]:
            if entering:
                active[side][index] = trait
            else:
                active[side].pop(index)
        if not any(active):
            continue
        strips += 1
        visits += len(active[0]) + len(active[1])
        if visits > max_visits:
            raise ValueError('Unissued coverage proof exceeds its work budget')
        y_events = {}
        for side in (0, 1):
            for lo, hi in active[side].values():
                change = y_events.setdefault(lo, [0, 0])
                change[side] += 1
                change = y_events.setdefault(hi, [0, 0])
                change[side] -= 1
        counts = [0, 0]
        y = sorted(y_events)
        for low, high in zip(y, y[1:]):
            for side in (0, 1):
                counts[side] += y_events[low][side]
            if low < high and (counts[0] != counts[1] or counts[0] > 1):
                raise ValueError('Candidate changes or overlaps unissued association coverage')
    return dict(variant_strips=strips, active_rectangle_visits=visits)


def unissued_frontier(snapshot, *, source_identity, reduction, total_traits,
                      job_variant_range,
                      max_partitions=10000, max_issued_ranges=100000,
                      max_coverage_visits=1000000):
    """Extract only unreserved suffixes from one active productive snapshot."""
    for name, value in (('max_partitions', max_partitions),
                        ('max_issued_ranges', max_issued_ranges),
                        ('max_coverage_visits', max_coverage_visits)):
        if type(value) is not int or value < 1:
            raise ValueError('Positive bounded ' + name + ' required')
    if (not isinstance(snapshot, dict) or snapshot.get('prefix_complete') is not True or
            snapshot.get('finished') is not None or snapshot.get('first_written') is None or
            snapshot.get('stop_reason') is not None):
        raise ValueError('Active productive snapshot after written output required')
    for name in ('issued_revision', 'written_events'):
        if type(snapshot.get(name)) is not int or snapshot[name] < 0:
            raise ValueError('Exact productive revision and writer progress required')
    if reduction not in (None, 'significant', 'jagwas') or type(total_traits) is not int or total_traits < 1:
        raise ValueError('Explicit native output mode and phenotype count required')
    job_variant = _span(job_variant_range, 'job variant')
    if (not isinstance(source_identity, dict) or 'path' not in source_identity or
            input_identity(source_identity['path']) != source_identity):
        raise ValueError('Bound live PGEN identity required')
    rows = snapshot.get('partitions')
    if not isinstance(rows, list) or not 0 < len(rows) <= max_partitions:
        raise ValueError('Bounded productive source partitions required')
    rectangles = []
    full_rectangles = []
    identifiers = set()
    range_count = 0
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError('Productive source partition required')
        key = row.get('id')
        if not isinstance(key, str) or not key or key in identifiers:
            raise ValueError('Unique productive source partition ids required')
        identifiers.add(key)
        lo, hi = _span(row.get('variant_range'), 'partition variant')
        trait = _span(row.get('trait_range'), 'partition phenotype')
        if trait[1] > total_traits or (reduction == 'jagwas' and trait != (0, total_traits)):
            raise ValueError('Productive phenotype panel differs from output mode')
        full_rectangles.append(dict(variant_range=[lo, hi], trait_range=list(trait)))
        cursor = row.get('cursor')
        ranges = row.get('ranges')
        count = row.get('issued_chunks')
        if (type(cursor) is not int or not lo <= cursor <= hi or
                type(count) is not int or count < 0 or
                not isinstance(ranges, list) or len(ranges) != count):
            raise ValueError('Incomplete productive source reservation prefix')
        range_count += count
        if range_count > max_issued_ranges:
            raise ValueError('Productive source prefix exceeds retained range budget')
        expected = lo
        for span in ranges:
            first, last = _span(span, 'issued source')
            if first != expected or last > hi:
                raise ValueError('Issued source ranges are not a contiguous prefix')
            expected = last
        if expected != cursor:
            raise ValueError('Source cursor differs from issued reservations')
        if cursor < hi:
            rectangles.append(dict(id=key, device=row.get('device'),
                                   variant_range=[cursor, hi],
                                   trait_range=list(trait)))
    if not rectangles:
        raise ValueError('No unissued association work remains')
    initial_coverage = _rectangles_match([
        dict(variant_range=list(job_variant), trait_range=[0, total_traits])],
        full_rectangles, max_coverage_visits)
    proof = _rectangles_match(rectangles, rectangles, max_coverage_visits)
    if input_identity(source_identity['path']) != source_identity:
        raise ValueError('PGEN input changed during frontier extraction')
    return dict(kind=KIND, input_identity=deepcopy(source_identity),
                reduction=reduction, total_traits=total_traits,
                job_variant_range=list(job_variant),
                issued_revision=snapshot['issued_revision'],
                written_events=snapshot['written_events'],
                rectangles=rectangles, initial_coverage=initial_coverage,
                required_cells=sum(
                    (row['variant_range'][1] - row['variant_range'][0]) *
                    (row['trait_range'][1] - row['trait_range'][0])
                    for row in rectangles), proof=proof,
                scope='Exact contiguous unissued suffixes after actual written output. Issued, prefetched and in-flight work stays on the original partitions. Not a live GPU/queue checkpoint or authorization to move writers.')


def bind_layout_to_frontier(frontier, layout, *, max_coverage_visits=1000000):
    """Prove a fixed calculator layout covers exactly the unissued pairs."""
    if (not isinstance(frontier, dict) or frontier.get('kind') != KIND or
            not isinstance(layout, dict) or
            layout.get('kind') != 'torchgwas.pgen_layout_source_floor.v1'):
        raise ValueError('Typed unissued frontier and source layout required')
    if type(max_coverage_visits) is not int or max_coverage_visits < 1:
        raise ValueError('Positive bounded coverage work required')
    for name in ('input_identity', 'reduction', 'total_traits'):
        if frontier.get(name) != layout.get(name):
            raise ValueError('Calculator layout differs from productive frontier: ' + name)
    source = frontier['input_identity']
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed before layout coverage proof')
    candidate = layout.get('partitions')
    required = frontier.get('rectangles')
    if not isinstance(candidate, list) or not candidate or not isinstance(required, list) or not required:
        raise ValueError('Nonempty calculator and productive partitions required')
    proof = _rectangles_match(required, candidate, max_coverage_visits)
    if input_identity(source['path']) != source:
        raise ValueError('PGEN input changed during layout coverage proof')
    return dict(kind='torchgwas.layout_frontier_coverage.v1',
                input_identity=deepcopy(source),
                issued_revision=frontier['issued_revision'],
                written_events=frontier['written_events'],
                required_cells=frontier['required_cells'],
                required_rectangles=len(required), candidate_partitions=len(candidate),
                proof=proof, exact_coverage=True,
                selection_validated=False,
                scope='Exact unissued variant-phenotype rectangle equality under a bounded sweep. Issued work, result queues, writer ownership, memory admission, service prices and live snapshot revalidation are separate obligations.')
