"""Early productive-run control with exact reserved source prefixes.

This internal bridge does not perform memory admission or price calibration.
Its caller must admit all sizes on the fixed partitions before starting work.
Planning callbacks use the analytical calculator; they may not wait for the
executor. A synchronous step holds the issue frontier while its cost is charged.
"""
from copy import deepcopy
from contextlib import nullcontext
import math
from pathlib import Path
import threading
import time

from .adaptive_chunks import AlignedChunkSizeControl, _positive_size
from .planning_session import IncrementalPlanningBudget, PlanningWorkCache
from .structural_tensor_cache import StructuralTensorWorkCache


class _PartitionControl:
    def __init__(self, run, key):
        self.run, self.key, self.capacity = run, key, run.capacity

    def __call__(self, start, stop, capacity):
        return self.run._reserve(self.key, start, stop, capacity)


class ProductiveTuningRun:
    """Connect real source reservations, written output and one-step planning.

    Partitions have unique ids, explicit devices, variant and trait ranges.
    Reserved reads include prefetched work and reads about to be submitted.
    Failed submissions fail the scan; reservations are never silently removed.
    After tuning ends only bounded counters advance, not a full-job range log.
    A written empty significant chunk is productive work but no material part.
    """
    def __init__(self, partitions, *, chunk_sizes, initial, max_issued_chunks=128,
                 budget=None, max_cache_entries=128, max_cache_bytes=32 << 20,
                 structural_cache_dir=None):
        self._control = AlignedChunkSizeControl(chunk_sizes, initial=initial)
        self.capacity = self._control.capacity
        self.max_issued_chunks = _positive_size(max_issued_chunks, 'max_issued_chunks')
        self.budget = budget if budget is not None else IncrementalPlanningBudget()
        if not isinstance(self.budget, IncrementalPlanningBudget):
            raise ValueError('An IncrementalPlanningBudget is required')
        self.cache = PlanningWorkCache(max_entries=max_cache_entries, max_bytes=max_cache_bytes)
        if structural_cache_dir is not None and (not isinstance(structural_cache_dir,(str,Path)) or not str(structural_cache_dir).strip()):
            raise ValueError('Nonempty structural cache directory required')
        self._structural_cache_dir=structural_cache_dir
        self._structural_cache=None
        self._structural_publication=None
        self._lock = threading.RLock()
        self._partitions = {}
        if not isinstance(partitions, (list, tuple)) or not partitions:
            raise ValueError('Explicit nonempty partitions required')
        for spec in partitions:
            if not isinstance(spec, dict) or set(spec) != {'id', 'device', 'variant_range', 'trait_range'}:
                raise ValueError('Partition id, device, variant_range and trait_range required')
            key, device = spec['id'], spec['device']
            if not isinstance(key, str) or not key or key in self._partitions:
                raise ValueError('Unique nonempty partition ids required')
            if not isinstance(device, str) or not device.startswith('cuda:') or not device[5:].isdigit():
                raise ValueError('Explicit CUDA partition device required')
            for name in ('variant_range', 'trait_range'):
                span = spec[name]
                if (not isinstance(span, (list, tuple)) or len(span) != 2
                        or any(type(x) is not int for x in span) or not 0 <= span[0] < span[1]):
                    raise ValueError('Nonempty integer partition ranges required')
            self._partitions[key] = dict(deepcopy(spec), cursor=spec['variant_range'][0],
                issued_chunks=0, ranges=[])
        self._created = time.perf_counter()
        self._first_written = self._first_material = self._first_fsynced_part = None
        self._first_material_output = None
        self._written_events = self._dense_statistic_bytes = 0
        self._revision = self._written_chunks = self._written_rows = self._part_bytes = 0
        self._retained_ranges = 0
        self._current_size = int(initial)
        self._stop_reason = None
        self._finished = None
        self._decisions = []

    def __call__(self, start, stop, capacity):
        raise ValueError('Bind the productive run to a specific scan partition')

    def for_partition(self, key):
        with self._lock:
            if key not in self._partitions:
                raise ValueError('Unknown productive-run partition')
        return _PartitionControl(self, key)

    def for_scan(self, device, variant_range, n_traits):
        """Automatic binding for unique variant shards; tile callers bind ids."""
        with self._lock:
            matches = [key for key, row in self._partitions.items()
                if row['device'] == str(device) and list(row['variant_range']) == list(variant_range)
                and row['trait_range'][1] - row['trait_range'][0] == n_traits]
            if len(matches) != 1:
                raise ValueError('Scan must identify exactly one productive-run partition; bind trait tiles explicitly')
        return self.for_partition(matches[0])

    def _stop(self, reason):
        if self._stop_reason is None:
            self._stop_reason = reason
            self.budget.finish()
            self.cache.close()

    def stop_planning(self, reason):
        """End optional work while source reservations and the scan continue."""
        if reason not in ('insufficient_unissued_work', 'partition_budget', 'source_fully_issued',
                          'extrapolation_limit', 'background_error',
                          'source_staging_evidence_only'):
            raise ValueError('Unknown productive planning stop reason')
        with self._lock:
            self._stop(reason)

    def _expire(self):
        if self._first_written is not None and time.perf_counter() - self._first_written >= self.budget.max_window_seconds:
            self._stop('window_budget')

    def _reserve(self, key, start, stop, capacity):
        with self._lock:
            if self._finished is not None:
                raise ValueError('Productive run is already finished')
            row = self._partitions[key]
            if (start != row['cursor'] or stop != row['variant_range'][1]
                    or capacity != self.capacity or start >= stop):
                raise ValueError('Source reservation differs from its exact contiguous partition')
            self._expire()
            count = self._control(start, stop, capacity)
            if self._retained_ranges >= self.max_issued_chunks:
                self._stop('issued_prefix_budget')
            if self._stop_reason is None:
                row['ranges'].append([start, start + count])
                self._retained_ranges += 1
            row['cursor'] += count
            row['issued_chunks'] += 1
            self._revision += 1
            return count

    def output_written(self, event):
        """Receive completed indexed writes or dense beta/t prefixes, not queues.

        Dense progress can split/combine scan chunks and precede df sidecar
        flushing. It starts useful-output planning but does not establish
        complete store durability or a measured writer capacity.
        """
        from .sumstats_indexed import IndexedChunkWrite
        from .sumstats import DenseWriteProgress
        if not isinstance(event, (IndexedChunkWrite,DenseWriteProgress)):
            raise ValueError('Typed writer completion event required')
        dense=isinstance(event,DenseWriteProgress)
        with self._lock:
            if self._finished is not None:
                raise ValueError('Output arrived after productive run completion')
            self._written_events += 1
            if not dense:self._written_chunks += 1
            self._written_rows += event.rows
            if dense:self._dense_statistic_bytes += event.statistic_bytes
            else:self._part_bytes += event.part_bytes
            if self._first_written is None:
                self._first_written = event.completed
                if self._stop_reason is None:
                    self.budget.start_after_first_output()
            else:
                self._first_written = min(self._first_written,event.completed)
            if event.rows:
                self._first_material_output = (event.completed if self._first_material_output is None
                    else min(self._first_material_output,event.completed))
            if not dense and event.rows:
                self._first_material = event.completed if self._first_material is None else min(self._first_material,event.completed)
            if not dense and event.part_file_fsynced:
                self._first_fsynced_part = (event.completed if self._first_fsynced_part is None
                    else min(self._first_fsynced_part,event.completed))
            self._expire()

    def planning_step(self, build, *, publication_seconds=None, reserve_seconds=0.,
                      release_issue_frontier=False, **forecasts):
        """Evaluate at most one proposal against a held source issue frontier.

        A proposal contains chunk_size, baseline_seconds and candidate_seconds
        for equivalent remaining work. It is applied only if its estimated gain
        repays cumulative planning wall time plus declared switching,
        publication and reserve costs. Publication cost must be explicit when
        structural persistence is enabled, including an explicitly chosen zero.
        These are forecasts, not a posterior or a guaranteed runtime bound.
        By default the source issue frontier is held during the step. The
        optional released mode captures one exact frontier, lets reservations
        and writer events continue, and rejects a result if either advanced.
        It is intended for a separate planner thread; calling it directly on
        a writer still blocks that writer until calculation completes.
        Executor progress inside the issued prefix continues in either mode.
        Actual end-of-job cache publication cost is reported separately.
        """
        if type(release_issue_frontier) is not bool:
            raise ValueError('release_issue_frontier must be boolean')
        guard = nullcontext() if release_issue_frontier else self._lock
        with guard:
            with self._lock:
                self._expire()
                if self._stop_reason is None and all(row['cursor'] == row['variant_range'][1] for row in self._partitions.values()):
                    self._stop('source_fully_issued')
                if self._stop_reason is not None:
                    return dict(evaluated=False, applied=False, reason=self._stop_reason)
                if self._first_written is not None and self._structural_cache_dir is not None and publication_seconds is None:
                    return dict(evaluated=False, applied=False, reason='publication_cost_required')
                held=self.snapshot()
            def evaluate():
                nonlocal held
                # Construction and loading occur only inside an admitted step,
                # after written output, and count against its actual budget.
                if self._structural_cache_dir is not None and self._structural_cache is None:
                    self._structural_cache=StructuralTensorWorkCache(self._structural_cache_dir)
                scope=nullcontext() if self._structural_cache is None else self._structural_cache.activate()
                with self.cache.activate(),scope:
                    # Include the lazily constructed cache in the bound state.
                    # In released mode this is the last short lock before the
                    # model runs; any later issue/output advance is rejected.
                    with self._lock:held=self.snapshot()
                    proposal = build(held)
                if not isinstance(proposal, dict) or set(proposal) != {'chunk_size', 'baseline_seconds', 'candidate_seconds'}:
                    raise ValueError('One analytical remaining-work proposal required')
                size = _positive_size(proposal['chunk_size'], 'proposed chunk size')
                if size not in self._control.sizes:
                    raise ValueError('Proposed size was not admitted')
                for key in ('baseline_seconds', 'candidate_seconds'):
                    value = proposal[key]
                    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                        raise ValueError('Finite nonnegative completion forecasts required')
                if proposal['baseline_seconds'] > forecasts['remaining_seconds']:
                    raise ValueError('Proposal exceeds the declared remaining horizon')
                return deepcopy(proposal)
            try:
                result = self.budget.run_step(evaluate,publication_seconds=0. if publication_seconds is None else publication_seconds,
                    reserve_seconds=reserve_seconds,**forecasts)
            except Exception as error:
                with self._lock:
                    self._stop('planning_error')
                    result = dict(evaluated=True, applied=False, error=type(error).__name__, message=str(error))
                    self._decisions.append(deepcopy(result))
                    return result
            with self._lock:
                result['applied'] = False
                stale=(release_issue_frontier and result.get('evaluated',False) and
                    (self._revision!=held['issued_revision'] or
                     self._written_events!=held['written_events'] or
                     self._current_size!=held['current_chunk_size'] or
                     self._stop_reason is not None or self._finished is not None))
                if stale:
                    result['usable_for_decision']=False
                    result['stale_frontier']=True
                    result['reason']='source_or_output_advanced'
                if result['evaluated'] and result['usable_for_decision']:
                    proposal = result['value']
                    gain = proposal['baseline_seconds'] - proposal['candidate_seconds']
                    if gain > result['total_tuning_cost_seconds']:
                        self._control.set_size(proposal['chunk_size'])
                        self._current_size = proposal['chunk_size']
                        result['applied'] = True
                    result['forecast_gain_seconds'] = gain
                result['planned_issued_revision']=held['issued_revision']
                result['issued_revision'] = self._revision
                # Rejected attempts do not grow an unbounded trace of caller polls.
                if result['evaluated']:
                    self._decisions.append(deepcopy(result))
                if result.get('evaluated') and not result.get('usable_for_decision', False) and not stale:
                    self._stop('planning_budget_or_horizon')
                return result

    def forecast_step(self, build_comparisons, *, forecast_options, price_profile=None, **cost_forecasts):
        """Construct bounded window comparisons inside the held issue frontier.

        The same productive budget charges header/model/forecast work together.
        Callers still supply valid prices, admitted shapes and explicit missing
        continuation costs. A stale, incomplete or incompatible model stops
        optional tuning while the scientific scan continues.
        An explicit price profile is checked before calculation and after the
        proposal, including exact comparison profile binding and original ages.
        These checks are charged; they certify only declared measurement fields.
        """
        audit=None
        def build(snapshot):
            nonlocal audit
            from .productive_forecast import productive_window_proposal
            bound=None
            if price_profile is not None:
                from .price_binding import validate_price_bindings,validate_comparison_prices
                bound=deepcopy(price_profile)
                if validate_price_bindings(bound)['status']!='declared_targets_verified':
                    raise ValueError('Productive price profile has no declared measurement evidence')
            comparisons=build_comparisons(snapshot)
            proposal,audit=productive_window_proposal(snapshot,comparisons,
                chunk_sizes=self._control.sizes,max_issued_chunks=self.max_issued_chunks,**forecast_options)
            if bound is not None:audit['price_evidence']=validate_comparison_prices(bound,comparisons)
            return proposal
        with self._lock:
            result=self.planning_step(build,**cost_forecasts)
            if audit is not None:
                result['forecast_audit']=deepcopy(audit)
                if result.get('evaluated') and self._decisions:
                    self._decisions[-1]['forecast_audit']=deepcopy(audit)
            return result

    def revision_token(self):
        """Compact issue/output token for nonblocking observation joins."""
        with self._lock:
            return dict(issued_revision=self._revision,
                        written_events=self._written_events,
                        current_chunk_size=self._current_size,
                        stop_reason=self._stop_reason,
                        finished=deepcopy(self._finished))

    def snapshot(self):
        with self._lock:
            return dict(partitions=deepcopy(list(self._partitions.values())), issued_revision=self._revision,
                prefix_complete=self._stop_reason is None, current_chunk_size=self._current_size,
                stop_reason=self._stop_reason, created=self._created, first_written=self._first_written,
                first_material_part=self._first_material, first_fsynced_part=self._first_fsynced_part,
                first_material_output=self._first_material_output,written_events=self._written_events,
                dense_statistic_bytes=self._dense_statistic_bytes,
                written_chunks=self._written_chunks, written_rows=self._written_rows, part_bytes=self._part_bytes,
                planning=self.budget.snapshot(), cache=self.cache.snapshot(), decisions=deepcopy(self._decisions),
                structural_cache=dict(enabled=self._structural_cache_dir is not None,
                    state=None if self._structural_cache is None else self._structural_cache.snapshot(),
                    publication=deepcopy(self._structural_publication)),
                finished=self._finished,
                scope='Reserved source ranges, completed indexed chunks and written dense beta/t prefixes. Dense progress need not follow df flushing or fsync and does not identify scan chunks. Not live GPU/queue service state. Caller owns memory admission and calibrated analytical forecasts.')

    def finish(self, *, successful):
        if type(successful) is not bool:
            raise ValueError('successful must be boolean')
        with self._lock:
            if self._finished is None:
                if successful and any(row['cursor'] != row['variant_range'][1] for row in self._partitions.values()):
                    raise ValueError('Successful run requires complete reserved source coverage')
                self._finished = dict(successful=successful, completed=time.perf_counter())
                self._stop('completed' if successful else 'failed')
                if self._structural_cache is not None:
                    started=time.perf_counter();cpu=time.thread_time()
                    try:
                        report=self._structural_cache.publish(successful=successful)
                        self._structural_publication=dict(report,wall_seconds=time.perf_counter()-started,
                            cpu_seconds=time.thread_time()-cpu)
                    except Exception as error:
                        # An optional cache failure must not invalidate already
                        # completed association output or leave resources open.
                        self._structural_publication=dict(status='cache_error',error=type(error).__name__,
                            wall_seconds=time.perf_counter()-started,cpu_seconds=time.thread_time()-cpu)
                    finally:self._structural_cache.close()
            return self.snapshot()
