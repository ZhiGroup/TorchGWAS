"""Compact read-only state for a productive indexed result consumer."""
from copy import deepcopy
import threading
import time


class IndexedWriterProgress:
    """Observe one synchronous indexed writer without retaining result arrays."""

    def __init__(self, kind):
        if kind not in ('jagwas', 'significant'):
            raise ValueError('Observed indexed writer mode required')
        self.kind = kind
        self._lock = threading.Lock()
        self._phase = 'waiting_for_result'
        self._active = None
        self._revision = 0
        self._completed_chunks = 0
        self._fsynced_parts = 0
        self._part_bytes = 0

    def begin_chunk(self, source_range, partition):
        from .sumstats_indexed import IndexedOutputPartition
        if (not isinstance(partition, IndexedOutputPartition) or
                not isinstance(source_range, tuple) or len(source_range) != 2 or
                any(type(value) is not int for value in source_range) or
                not partition.variant_range[0] <= source_range[0] <
                    source_range[1] <= partition.variant_range[1]):
            raise ValueError('Observed indexed result requires one owned source range')
        with self._lock:
            if self._phase != 'waiting_for_result' or self._active is not None:
                raise ValueError('Indexed writer already has an active result')
            self._active = dict(device=partition.device,
                variant_range=list(source_range),
                trait_range=list(partition.trait_range))
            self._phase = 'selecting'
            self._revision += 1

    def begin_emit(self):
        with self._lock:
            if self._phase != 'selecting' or self._active is None:
                raise ValueError('Indexed writer cannot emit without an active result')
            self._phase = 'emitting_part'
            self._revision += 1

    def complete_chunk(self, rows, part_bytes, fsynced):
        if (type(rows) is not int or rows < 0 or
                type(part_bytes) is not int or part_bytes < 0 or
                type(fsynced) is not bool or (rows == 0 and (part_bytes or fsynced))):
            raise ValueError('Observed indexed part completion required')
        with self._lock:
            if self._phase != 'emitting_part' or self._active is None:
                raise ValueError('Indexed writer has no emitted result to complete')
            self._active = None
            self._phase = 'waiting_for_result'
            self._revision += 1
            self._completed_chunks += 1
            self._fsynced_parts += int(fsynced)
            self._part_bytes += part_bytes

    def begin_finalization(self):
        with self._lock:
            if self._phase != 'waiting_for_result' or self._active is not None:
                raise ValueError('Indexed writer cannot finalize with an active result')
            self._phase = 'writing_metadata'
            self._revision += 1

    def begin_publish(self):
        with self._lock:
            if self._phase != 'writing_metadata':
                raise ValueError('Indexed writer cannot publish before metadata')
            self._phase = 'publishing_manifest'
            self._revision += 1

    def published(self):
        with self._lock:
            if self._phase != 'publishing_manifest':
                raise ValueError('Indexed writer cannot publish twice')
            self._phase = 'published'
            self._revision += 1

    def failed(self):
        with self._lock:
            if self._phase != 'published':
                self._phase = 'failed'
                self._revision += 1

    def snapshot(self):
        with self._lock:
            anchor = time.perf_counter()
            return dict(kind='torchgwas.indexed_writer_progress.v1',
                mode=self.kind, capture_anchor_seconds=anchor,
                revision=self._revision, phase=self._phase,
                active=deepcopy(self._active),
                completed_chunks=self._completed_chunks,
                fsynced_parts=self._fsynced_parts,
                part_bytes=self._part_bytes,
                scope='One synchronous indexed consumer. Active result arrays, selection and part service are not timed or priced; completed parts exclude final metadata and manifest durability.')
