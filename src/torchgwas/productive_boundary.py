"""Bounded written-output progress at a productive JIT issue frontier.

Issued source work and completed output are distinct. This observer retains
only early per-partition prefixes/ranges needed to describe their boundary;
it does not inspect private GPU streams or poll writer queues and never changes
the scientific scan on a malformed optional observation.
"""
from copy import deepcopy


class ProductiveBoundaryProgress:
    def __init__(self, partitions, reduction, *, max_events=1024):
        if reduction not in (None, 'significant', 'jagwas'):
            raise ValueError('Unknown productive output mode')
        if type(max_events) is not int or not 1 <= max_events <= 4096:
            raise ValueError('Bounded positive output-event budget required')
        if not isinstance(partitions, (list, tuple)) or not partitions:
            raise ValueError('Explicit productive partitions required')
        self.reduction = reduction
        self.max_events = max_events
        self.partitions = {}
        for row in partitions:
            if (not isinstance(row, dict) or set(row) !=
                    {'id', 'device', 'variant_range', 'trait_range'} or
                    not isinstance(row['id'], str) or not row['id'] or
                    row['id'] in self.partitions):
                raise ValueError('Unique productive output partition required')
            self.partitions[row['id']] = deepcopy(row)
        self.job_variant_start = min(
            row['variant_range'][0] for row in self.partitions.values())
        self.events = 0
        self.invalid = 0
        self.truncated = False
        self.latest_writer_queue = None
        self.writer_queues = {}
        self.rows = {key: dict(indexed_ranges=[], rows=0, part_bytes=0,
                               fsynced_parts=0,
                               matrix_to=row['variant_range'][0],
                               df_to=row['variant_range'][0],
                               matrix_statistic_bytes=0, events=0)
                     for key, row in self.partitions.items()}

    def _match_indexed(self, event):
        partition = event.partition
        source = event.source_variant_range
        if partition is None or source is None or event.kind != self.reduction:
            return None
        matches = [key for key, row in self.partitions.items()
                   if (row['device'] == partition.device and
                       tuple(row['variant_range']) == partition.variant_range and
                       tuple(row['trait_range']) == partition.trait_range and
                       row['variant_range'][0] <= source[0] < source[1] <=
                       row['variant_range'][1])]
        return matches[0] if len(matches) == 1 else None

    def _match_dense(self, event):
        lo = event.start + self.job_variant_start
        hi = event.end + self.job_variant_start
        matches = [key for key, row in self.partitions.items()
                   if ((event.device is None or row['device'] == event.device) and
                       tuple(row['trait_range']) == event.trait_range and
                       row['variant_range'][0] <= lo < hi <=
                       row['variant_range'][1])]
        # The native binary writer emits a global completed row prefix without
        # a producer device. Source/trait ownership is still exact when only
        # one admitted partition contains the whole written interval. A
        # cross-shard write, or overlapping candidate owners, stays invalid.
        return (matches[0], lo, hi) if len(matches) == 1 else (None, lo, hi)

    def observe(self, event):
        from .sumstats import DenseWriteProgress
        from .sumstats_indexed import IndexedChunkWrite
        if not isinstance(event, (DenseWriteProgress, IndexedChunkWrite)):
            self.invalid += 1
            return
        self.events += 1
        if self.events > self.max_events:
            self.truncated = True
            return
        try:
            if isinstance(event, IndexedChunkWrite):
                if self.reduction is None:
                    raise ValueError('Indexed output differs from dense admission')
                key = self._match_indexed(event)
                if key is None:
                    raise ValueError('Unbound indexed output')
                lo, hi = event.source_variant_range
                state = self.rows[key]
                if (type(lo) is not int or type(hi) is not int or
                        type(event.rows) is not int or event.rows < 0 or
                        type(event.part_bytes) is not int or event.part_bytes < 0 or
                        any(lo < end and first < hi
                            for first, end in state['indexed_ranges'])):
                    raise ValueError('Invalid or duplicate indexed completion')
                state['indexed_ranges'].append([lo, hi])
                state['rows'] += event.rows
                state['part_bytes'] += event.part_bytes
                state['fsynced_parts'] += bool(event.part_file_fsynced)
            else:
                if self.reduction is not None:
                    raise ValueError('Dense output differs from indexed admission')
                queue = event.writer_queue
                if queue is not None:
                    if (not isinstance(queue, dict) or
                            queue.get('kind') !=
                            'torchgwas.dense_writer_queue_observation.v1' or
                            queue.get('valid') is not True or
                            not isinstance(event.directory, str) or
                            not event.directory):
                        raise ValueError('Invalid dense writer queue observation')
                key, lo, hi = self._match_dense(event)
                if key is None:
                    raise ValueError('Unbound dense output')
                state = self.rows[key]
                df = (None if event.variant_df_complete_to is None else
                      event.variant_df_complete_to + self.job_variant_start)
                if (type(event.start) is not int or type(event.end) is not int or
                        lo != state['matrix_to'] or
                        type(event.statistic_bytes) is not int or
                        event.statistic_bytes < 0 or
                        (df is not None and
                         (type(df) is not int or not state['df_to'] <= df <=
                          self.partitions[key]['variant_range'][1]))):
                    raise ValueError('Noncontiguous dense write progress')
                state['matrix_to'] = hi
                if df is not None:
                    state['df_to'] = df
                state['matrix_statistic_bytes'] += event.statistic_bytes
                if queue is not None:
                    observed=dict(written_event=self.events,
                                  event_partition_id=key,
                                  directory=event.directory,
                                  observation=deepcopy(queue),
                                  scope='Queue belongs to the writer directory; event_partition_id identifies only this completed prefix. A shared writer may also hold another partition.')
                    self.latest_writer_queue=observed
                    self.writer_queues[event.directory]=observed
            state['events'] += 1
        except (AttributeError, IndexError, KeyError, TypeError, ValueError):
            self.invalid += 1

    def bind(self, issued_snapshot):
        """Reconcile output progress with a held productive issue snapshot."""
        if (not isinstance(issued_snapshot, dict) or
                issued_snapshot.get('prefix_complete') is not True or
                issued_snapshot.get('written_events') != self.events or
                not isinstance(issued_snapshot.get('partitions'), list)):
            raise ValueError('Complete synchronized productive snapshot required')
        issued = {row['id']: row for row in issued_snapshot['partitions']}
        if set(issued) != set(self.partitions):
            raise ValueError('Issued and output partitions differ')
        rows = []
        valid = self.invalid == 0 and not self.truncated
        pending_pairs = 0
        for key, spec in self.partitions.items():
            source = issued[key]
            state = self.rows[key]
            if any(source.get(name) != spec[name]
                   for name in ('device', 'variant_range', 'trait_range')):
                valid = False
            ranges = source.get('ranges')
            cursor = source.get('cursor')
            if (not isinstance(ranges, list) or type(cursor) is not int or
                    not spec['variant_range'][0] <= cursor <=
                    spec['variant_range'][1]):
                valid = False
                ranges = []
            if self.reduction is None:
                matrix_to = state['matrix_to']
                if matrix_to > cursor or state['df_to'] > cursor:
                    valid = False
                width = spec['trait_range'][1] - spec['trait_range'][0]
                pending = max(0, cursor - matrix_to) * width
                pending_pairs += pending
                rows.append(dict(id=key, device=spec['device'],
                    variant_range=list(spec['variant_range']),
                    trait_range=list(spec['trait_range']),
                    issued_to=cursor, issued_ranges=deepcopy(ranges),
                    matrix_written_to=matrix_to, df_written_to=state['df_to'],
                    issued_not_matrix_written=[matrix_to, cursor],
                    issued_not_matrix_written_pairs=pending,
                    issued_not_df_written_markers=max(0, cursor-state['df_to']),
                    matrix_statistic_bytes=state['matrix_statistic_bytes'],
                    output_events=state['events']))
            else:
                completed = sorted(map(tuple, state['indexed_ranges']))
                reserved = {tuple(span) for span in ranges}
                if any(span not in reserved for span in completed):
                    valid = False
                completed_set = set(completed)
                unwritten = [list(span) for span in ranges
                             if tuple(span) not in completed_set]
                width = spec['trait_range'][1] - spec['trait_range'][0]
                pending = sum(hi-lo for lo,hi in unwritten)*width
                pending_pairs += pending
                rows.append(dict(id=key, device=spec['device'],
                    variant_range=list(spec['variant_range']),
                    trait_range=list(spec['trait_range']),
                    issued_to=cursor, issued_ranges=deepcopy(ranges),
                    indexed_written_ranges=[list(span) for span in completed],
                    issued_not_indexed_written=unwritten,
                    issued_not_indexed_written_pairs=pending,
                    indexed_rows=state['rows'], part_bytes=state['part_bytes'],
                    fsynced_parts=state['fsynced_parts'],
                    output_events=state['events']))
        return dict(kind='torchgwas.productive_output_boundary.v1',
                    issued_revision=issued_snapshot['issued_revision'],
                    written_events=self.events, reduction=self.reduction,
                    valid=valid, invalid_events=self.invalid,
                    truncated=self.truncated, partitions=rows,
                    issued_not_written_pairs=pending_pairs,
                    writer_queue_observation=deepcopy(self.latest_writer_queue),
                    writer_queues=deepcopy(self.writer_queues),
                    scope='Bounded issued versus completed output progress. Dense writer queue observations atomically cover one writer at callbacks, not a synchronized execution checkpoint. Dense beta/t prefixes may precede df and fsync; indexed part fsync precedes manifest/directory durability. GPU streams, producer queues and remaining service within issued chunks are unknown.')

    def snapshot(self):
        return dict(kind='torchgwas.productive_output_observations.v1',
                    reduction=self.reduction, written_events=self.events,
                    invalid_events=self.invalid, truncated=self.truncated,
                    max_events=self.max_events, partitions=deepcopy(self.rows),
                    writer_queue_observation=deepcopy(self.latest_writer_queue),
                    writer_queues=deepcopy(self.writer_queues))
