"""Exact whole-header PGEN schedule assembled in bounded productive steps.

This work is duration-free metadata accounting. A caller can advance it after
successive completed output chunks rather than paying the whole source walk in
one cold callback. It never extrapolates a sampled source window.
"""

from collections import Counter
from copy import deepcopy
import math
import time

import numpy as np

from .analytical_plan_cache import input_identity
from .pgen_work_bounds import (PgenHeaderWork, SCHEDULE_KIND,
                               _replay_entries_work)
from .pgen_reader import ld_safe_start, PgenFormatError


_SCALARS = ('record_payload_bytes', 'read_bytes', 'decode_input_bytes',
            'native_ld_base_update_bytes', 'native_ld_replay_packed_bytes',
            'logical_final_packed_bytes', 'logical_final_int8_bytes')


def _sum_intervals(target, source):
    for key, pair in source.items():
        if (not isinstance(pair, (list, tuple)) or len(pair) != 2 or
                any(type(value) is not int or value < 0 for value in pair) or
                pair[0] > pair[1]):
            raise ValueError('Invalid incremental source interval')
        if pair[1] == 0:
            continue
        current = target.setdefault(key, [0, 0])
        current[0] += pair[0]
        current[1] += pair[1]


def _subtract_intervals(total, removed):
    result = {}
    for name in set(total) | set(removed):
        a, b = total.get(name, [0, 0])
        c, d = removed.get(name, [0, 0])
        pair = [a - c, b - d]
        if pair[0] < 0 or pair[0] > pair[1]:
            raise ValueError('Incremental source intervals do not conserve work')
        if pair != [0, 0]:
            result[name] = pair
    return result


def _initial_replay(header, at):
    bound = header.bounds(at, at + 1, max_records=1)
    replay = bound['ld_replay']
    if replay is None:
        return None
    units, updates, varints = header._record_work(
        bound['samples'], replay['record_form'],
        replay['record_bytes'], expand=False)
    return dict(entry=[at, replay['base_variant'],
                       replay['read_prefix_bytes'], replay['record_form'],
                       replay['record_bytes']],
                source_units=units, varint_constraints=varints,
                base_update_bytes=updates)


def _primary_segment(header, first, last, signatures):
    bound = header.bounds(first, last, max_records=last - first,
                          max_signatures=32768)
    replay = _initial_replay(header, first)
    units = _subtract_intervals(bound['source_units'],
                                {} if replay is None else replay['source_units'])
    varints = _subtract_intervals(bound['varint_constraints'],
                                  {} if replay is None else replay['varint_constraints'])
    base = bound['native_ld_base_update_bytes'] - (
        0 if replay is None else replay['base_update_bytes'])
    if base < 0:
        raise ValueError('Initial LD replay exceeds primary base updates')
    return dict(start=first, stop=last, source_units=units,
                varint_constraints=varints,
                base_update_bytes=base,
                record_form_counts=bound['record_form_counts'],
                record_payload_bytes=bound['record_payload_bytes'],
                signatures=signatures)


class IncrementalPgenSchedule:
    """One fixed chunk grid; `advance` inspects at most one bounded segment."""

    def __init__(self, header, start, stop, chunk_markers, *,
                 records_per_step=65536, max_chunks_per_step=65536,
                 max_total_records=20_000_000, max_total_chunks=200_000,
                 max_segments=512, max_signatures=32768):
        if not isinstance(header, PgenHeaderWork):
            raise ValueError('Bound PGEN header required')
        for name, value in (('start', start), ('stop', stop),
                            ('chunk_markers', chunk_markers),
                            ('records_per_step', records_per_step),
                            ('max_chunks_per_step', max_chunks_per_step),
                            ('max_total_records', max_total_records),
                            ('max_total_chunks', max_total_chunks),
                            ('max_segments', max_segments),
                            ('max_signatures', max_signatures)):
            if type(value) is not int or value < (0 if name == 'start' else 1):
                raise ValueError('Bounded positive ' + name + ' required')
        if (not 0 <= start < stop <= header._header.variant_ct or
                stop - start > max_total_records):
            raise ValueError('Incremental source span exceeds the file or budget')
        chunks = (stop - start + chunk_markers - 1) // chunk_markers
        if chunks > max_total_chunks:
            raise ValueError('Incremental chunk budget exceeded')
        step_chunks = min(records_per_step // chunk_markers,
                          max_chunks_per_step)
        if step_chunks < 1 or math.ceil(chunks / step_chunks) > max_segments:
            raise ValueError('Incremental segment budget exceeded')
        source = header.input_identity
        if input_identity(source['path']) != source:
            raise ValueError('PGEN input changed before incremental schedule')
        self.header = header
        self.input_identity = deepcopy(source)
        self.start, self.stop = start, stop
        self.chunk_markers = chunk_markers
        self.step_records = step_chunks * chunk_markers
        self.max_chunks_per_step = max_chunks_per_step
        self.max_signatures = max_signatures
        self.max_segments = max_segments
        self.max_total_chunks = max_total_chunks
        self.chunk_count = chunks
        self.cursor = start
        self.segments = 0
        self.wall_seconds = 0.
        self.cpu_seconds = 0.
        self._first = None
        self._totals = {name: 0 for name in _SCALARS}
        self._units = {}
        self._varints = {}
        self._forms = Counter()
        self._signatures = set()
        self._entries = []
        self._ld_replays = 0
        self._primary_segments = []

    def advance(self):
        """Price one aligned segment, or return completed progress unchanged."""
        if self.cursor == self.stop:
            return self.snapshot()
        if self.segments >= self.max_segments:
            raise ValueError('Incremental segment budget exhausted')
        began = time.perf_counter()
        cpu_began = time.thread_time()
        first = self.cursor
        last = min(first + self.step_records, self.stop)
        schedule = self.header.schedule_bounds(
            first, last, self.chunk_markers,
            max_records=self.step_records,
            max_chunks=self.max_chunks_per_step,
            max_signatures=self.max_signatures)
        if (schedule['input_identity'] != self.input_identity or
                schedule['variant_range'] != [first, last]):
            raise ValueError('Incremental PGEN source changed')
        entries = []
        if first != self.start:
            boundary = self.header.bounds(first, first + 1,
                max_records=1, max_signatures=self.max_signatures)
            replay = boundary['ld_replay']
            if replay is not None:
                entries.append([first, replay['base_variant'],
                                replay['read_prefix_bytes'],
                                replay['record_form'], replay['record_bytes']])
        entries.extend(schedule['additional_replay_work']['entries'])
        # The signature set tracks primary records only. It is accumulated
        # here so finalization does not rescan the entire 8M-record header.
        h = self.header._header
        keys = ((h.vrtypes[first:last].astype(np.uint64) << np.uint64(32)) |
                h.record_lengths[first:last])
        signatures = set(int(value) for value in np.unique(keys))
        if len(self._signatures | signatures) > self.max_signatures:
            raise ValueError('Incremental global signature budget exceeded')
        primary = _primary_segment(self.header, first, last, signatures)
        if self._first is None:
            self._first = deepcopy(schedule)
        for name in _SCALARS:
            self._totals[name] += schedule[name]
        _sum_intervals(self._units, schedule['source_units'])
        _sum_intervals(self._varints, schedule['varint_constraints'])
        self._forms.update(schedule['record_form_counts'])
        self._signatures.update(signatures)
        self._entries.extend(entries)
        self._ld_replays += schedule['ld_replay_count']
        self._primary_segments.append(primary)
        self.cursor = last
        self.segments += 1
        self.wall_seconds += time.perf_counter() - began
        self.cpu_seconds += time.thread_time() - cpu_began
        return self.snapshot()

    def snapshot(self):
        return dict(kind='torchgwas.incremental_pgen_schedule.v1',
                    input_identity=deepcopy(self.input_identity),
                    variant_range=[self.start, self.stop],
                    chunk_markers=self.chunk_markers,
                    cursor=self.cursor, complete=self.cursor == self.stop,
                    completed_segments=self.segments,
                    remaining_records=self.stop - self.cursor,
                    total_chunks=self.chunk_count,
                    observed_replay_entries=len(self._entries),
                    calculation_wall_seconds=self.wall_seconds,
                    calculation_cpu_seconds=self.cpu_seconds,
                    scope='Exact bounded header segments on one fixed chunk grid. No payload reads, service-time completion estimate or JIT switch authorization.')

    def finish(self):
        """Return the same conserved whole-schedule work after all segments."""
        if self.cursor != self.stop or self._first is None:
            raise ValueError('Incremental whole-source schedule is incomplete')
        if input_identity(self.input_identity['path']) != self.input_identity:
            raise ValueError('PGEN input changed before schedule finalization')
        extra = _replay_entries_work(self._entries,
            self._first['samples'], [self.start, self.stop],
            self.chunk_markers)
        common_units = {}
        for name in set(self._units) | set(extra['source_units']):
            total = self._units.get(name, [0, 0])
            replay = extra['source_units'].get(name, [0, 0])
            pair = [total[i] - replay[i] for i in (0, 1)]
            if pair[0] < 0 or pair[0] > pair[1]:
                raise ValueError('Incremental invariant source work differs')
            if pair != [0, 0]:
                common_units[name] = pair
        common = dict(input_identity=deepcopy(self.input_identity),
            variant_range=[self.start, self.stop],
            samples=self._first['samples'], markers=self.stop - self.start,
            record_form_counts=dict(self._forms),
            record_payload_bytes=self._totals['record_payload_bytes'],
            source_units=common_units)
        for name in ('read_bytes', 'decode_input_bytes',
                     'native_ld_base_update_bytes',
                     'native_ld_replay_packed_bytes'):
            value = self._totals[name] - extra[name]
            if value < 0:
                raise ValueError('Incremental replay work exceeds total')
            common[name] = value
        result = deepcopy(self._first)
        result.update(kind=SCHEDULE_KIND, input_identity=deepcopy(self.input_identity),
            variant_range=[self.start, self.stop],
            markers=self.stop - self.start,
            chunk_count=self.chunk_count,
            ld_replay_count=self._ld_replays,
            record_form_counts=dict(self._forms),
            source_units=deepcopy(self._units),
            varint_constraints=deepcopy(self._varints),
            common_work=common,
            additional_replay_work=dict(extra, entries=deepcopy(self._entries)),
            structural_work=dict(records=self.stop - self.start,
                signatures=len(self._signatures),
                replay_records=self._ld_replays,
                payload_bytes_read=0, chunks=self.chunk_count),
            scope='Exact whole-schedule source-count intervals assembled from bounded post-output header steps. No payload census, elapsed completion or switch authorization.')
        result.update(self._totals)
        if (self._ld_replays !=
                int(common['native_ld_replay_packed_bytes'] > 0) +
                extra['replay_count']):
            raise ValueError('Incremental LD replay count differs')
        if input_identity(self.input_identity['path']) != self.input_identity:
            raise ValueError('PGEN input changed during schedule finalization')
        return result

    def rebase(self, cursor, chunk_markers, *, stop=None):
        """Reuse staged primary work for any later cursor and admitted B.

        Only the partial leading segment and actual LD chunk starts are
        inspected. This remains a source-work report, not a runtime forecast.
        """
        if self.cursor != self.stop or self._first is None:
            raise ValueError('Incremental primary source ledger is incomplete')
        stop = self.stop if stop is None else stop
        if (type(cursor) is not int or type(stop) is not int or
                not self.start <= cursor < stop <= self.stop or
                type(chunk_markers) is not int or chunk_markers < 1):
            raise ValueError('Bound later source range and chunk size required')
        chunks = (stop - cursor + chunk_markers - 1) // chunk_markers
        if chunks > self.max_total_chunks:
            raise ValueError('Rebased source chunk budget exceeded')
        if input_identity(self.input_identity['path']) != self.input_identity:
            raise ValueError('PGEN input changed before source rebase')
        h = self.header._header
        primary_units = {}
        primary_varints = {}
        forms = Counter()
        signatures = set()
        payload = base_updates = 0
        for segment in self._primary_segments:
            if segment['stop'] <= cursor or segment['start'] >= stop:
                continue
            first = max(cursor, segment['start'])
            last = min(stop, segment['stop'])
            if first != segment['start'] or last != segment['stop']:
                keys = ((h.vrtypes[first:last].astype(np.uint64)
                         << np.uint64(32)) |
                        h.record_lengths[first:last])
                current = _primary_segment(self.header, first, last,
                                           set(int(value) for value in np.unique(keys)))
            else:
                current = segment
            _sum_intervals(primary_units, current['source_units'])
            _sum_intervals(primary_varints, current['varint_constraints'])
            forms.update(current['record_form_counts'])
            signatures.update(current['signatures'])
            payload += current['record_payload_bytes']
            base_updates += current['base_update_bytes']
        if len(signatures) > self.max_signatures:
            raise ValueError('Rebased source signature budget exceeded')
        initial = _initial_replay(self.header, cursor)
        common_units = deepcopy(primary_units)
        common_varints = deepcopy(primary_varints)
        if initial is not None:
            _sum_intervals(common_units, initial['source_units'])
            _sum_intervals(common_varints, initial['varint_constraints'])
        starts = np.arange(cursor + chunk_markers, stop,
                           chunk_markers, dtype=np.int64)
        ld = starts[np.isin(h.vrtypes[starts], (2, 3))]
        if self.header._bases is None:
            bases = (int(ld_safe_start(h.vrtypes, int(at))) for at in ld)
        else:
            locations = np.searchsorted(self.header._bases, ld,
                                        side='right') - 1
            bases = (int(value) for value in self.header._bases[locations])
        entries = []
        for at, base in zip(ld, bases):
            at = int(at)
            if not 0 <= base < at or int(h.vrtypes[base]) in (2, 3):
                raise PgenFormatError('LD schedule start has no valid earlier base')
            entries.append([at, base,
                int(h.record_offsets[at]) - int(h.record_offsets[base]),
                int(h.vrtypes[base]), int(h.record_lengths[base])])
        extra = _replay_entries_work(entries, self._first['samples'],
                                     [cursor, stop], chunk_markers,
                                     include_varint_constraints=True)
        extra_varints = extra.pop('varint_constraints')
        units = deepcopy(common_units)
        varints = deepcopy(common_varints)
        _sum_intervals(units, extra['source_units'])
        _sum_intervals(varints, extra_varints)
        packed = (self._first['samples'] + 3) // 4
        common = dict(input_identity=deepcopy(self.input_identity),
            variant_range=[cursor, stop],
            samples=self._first['samples'], markers=stop - cursor,
            record_form_counts=dict(forms), record_payload_bytes=payload,
            source_units=common_units,
            read_bytes=payload + (0 if initial is None else initial['entry'][2]),
            decode_input_bytes=payload + (0 if initial is None else initial['entry'][4]),
            native_ld_base_update_bytes=base_updates + (
                0 if initial is None else initial['base_update_bytes']),
            native_ld_replay_packed_bytes=packed * int(initial is not None))
        result = deepcopy(self._first)
        result.update(kind=SCHEDULE_KIND,
            input_identity=deepcopy(self.input_identity),
            variant_range=[cursor, stop],
            samples=self._first['samples'], markers=stop - cursor,
            record_form_counts=dict(forms),
            record_payload_bytes=payload,
            source_units=units, varint_constraints=varints,
            read_bytes=common['read_bytes'] + extra['read_bytes'],
            decode_input_bytes=(common['decode_input_bytes'] +
                                extra['decode_input_bytes']),
            native_ld_base_update_bytes=(
                common['native_ld_base_update_bytes'] +
                extra['native_ld_base_update_bytes']),
            native_ld_replay_packed_bytes=(
                common['native_ld_replay_packed_bytes'] +
                extra['native_ld_replay_packed_bytes']),
            logical_final_packed_bytes=packed * (stop - cursor),
            logical_final_int8_bytes=self._first['samples'] * (stop - cursor),
            chunk_markers=chunk_markers, chunk_count=chunks,
            ld_replay_count=int(initial is not None) + len(entries),
            common_work=common,
            additional_replay_work=dict(extra, entries=entries),
            structural_work=dict(records=stop - cursor,
                signatures=len(signatures),
                replay_records=int(initial is not None) + len(entries),
                payload_bytes_read=0, chunks=chunks),
            scope='Exact future PGEN source schedule rebased from staged primary work to a later cursor and chunk grid. No payload census, elapsed completion or JIT switch authorization.')
        if input_identity(self.input_identity['path']) != self.input_identity:
            raise ValueError('PGEN input changed during source rebase')
        return result
