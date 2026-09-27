"""Bounded observations from the first completed reduced-output chunks.

These are loaded writer observations. They can falsify an empty/full-only
scenario set after partial output, but cannot establish future selectivity or an
independent writer capacity.
"""
from copy import deepcopy
import math

class ProductiveOutputSample:
    def __init__(self, partitions, reduction, *, max_events=32):
        if reduction not in (None, 'significant', 'jagwas'):
            raise ValueError('Unknown productive reduction')
        if type(max_events) is not int or not 1 <= max_events <= 128:
            raise ValueError('Bounded positive output sample required')
        self.reduction = reduction
        self.max_events = max_events
        self.partitions = {p['id']: p for p in partitions}
        self.events = []
        self.counts = dict(indexed_events=0, empty=0, partial=0, full=0,
                           unbound=0, invalid=0, sampled=0)
        self._seen = set()

    def observe(self, event):
        if self.reduction is None:
            return
        from .sumstats_indexed import IndexedChunkWrite
        if not isinstance(event, IndexedChunkWrite):
            return
        self.counts['indexed_events'] += 1
        try:
            partition = event.partition
            source = event.source_variant_range
            if partition is None or source is None:
                self.counts['unbound'] += 1
                return
            matches = [p for p in self.partitions.values()
                       if p['device'] == partition.device
                       and tuple(p['variant_range']) == partition.variant_range
                       and tuple(p['trait_range']) == partition.trait_range]
            if len(matches) != 1 or event.kind != self.reduction:
                raise ValueError('Output producer or reduction differs from admission')
            lo, hi = source
            if (type(lo) is not int or type(hi) is not int or
                    not partition.variant_range[0] <= lo < hi <= partition.variant_range[1]):
                raise ValueError('Output source interval differs from partition')
            width = partition.trait_range[1] - partition.trait_range[0]
            possible = (hi - lo) * (width if self.reduction == 'significant' else 1)
            if type(event.rows) is not int or not 0 <= event.rows <= possible:
                raise ValueError('Output survivor count exceeds source work')
            if (type(event.part_bytes) is not int or event.part_bytes < 0 or
                    not math.isfinite(event.started) or not math.isfinite(event.completed) or
                    event.completed < event.started):
                raise ValueError('Invalid completed writer observation')
            key = (matches[0]['id'], lo, hi)
            # Duplicate detection covers the retained first-window sample.
            # Later events keep only aggregate counters, never a job-sized set.
            if len(self._seen) < self.max_events:
                if key in self._seen:
                    raise ValueError('Duplicate completed output source interval')
                self._seen.add(key)
            regime = 'empty' if event.rows == 0 else 'full' if event.rows == possible else 'partial'
            self.counts[regime] += 1
            if len(self.events) < self.max_events:
                self.events.append(dict(partition_id=matches[0]['id'], device=partition.device,
                                        variant_range=[lo, hi], trait_range=list(partition.trait_range),
                                        rows=event.rows, possible=possible,
                                        part_bytes=event.part_bytes,
                                        part_file_fsynced=event.part_file_fsynced,
                                        writer_wall_seconds=event.completed-event.started))
                self.counts['sampled'] += 1
        except (ValueError, TypeError, KeyError):
            self.counts['invalid'] += 1

    def check_scenarios(self, occupancies):
        """Fail optional planning if a completed chunk refutes every scenario."""
        if self.reduction is None:
            return
        if self.counts['invalid'] or self.counts['unbound']:
            raise ValueError('Productive output occupancy cannot be bound to the source')
        from .productive_occupancy import occupancy_kind
        scenarios = {occupancy_kind(value) for value in occupancies.values()}
        if scenarios == {'empty'} and (self.counts['partial'] or self.counts['full']):
            raise ValueError('Observed reduced output refutes the empty-only forecast')
        if scenarios == {'dense'} and (self.counts['empty'] or self.counts['partial']):
            raise ValueError('Observed reduced output refutes the full-only forecast')
        if self.counts['partial'] and 'sparse' not in scenarios:
            raise ValueError('Observed partial output requires a sparse occupancy scenario')

    def snapshot(self):
        return dict(counts=deepcopy(self.counts), events=deepcopy(self.events),
                    max_events=self.max_events, reduction=self.reduction,
                    scope='First completed indexed chunks only. Writer intervals include loaded formatting and optional part fsync; they are not independent capacity prices. Observed survivor counts do not establish future occupancy or renew any cached price.')
