"""Post-output, bounded PGEN metadata staging for future JIT decisions.

One useful output event earns at most one source-accounting step.  This is a
read-only evidence collector; it never changes an issued or future layout.
"""

import threading
import time
import sys

def _retained_size(value, seen=None):
    """Count Python containers, including the per-segment signature sets."""
    seen = set() if seen is None else seen
    if id(value) in seen:
        return 0
    seen.add(id(value))
    size = sys.getsizeof(value)
    if isinstance(value, dict):
        size += sum(_retained_size(key, seen) + _retained_size(item, seen)
                    for key, item in value.items())
    elif isinstance(value, (list, tuple, set, frozenset)):
        size += sum(_retained_size(item, seen) for item in value)
    return size


def validate_source_stage_config(config, initial_chunk):
    if type(initial_chunk) is not int or initial_chunk < 1:
        raise ValueError('Positive initial source chunk required')
    required = {'records_per_step', 'max_steps', 'max_cpu_seconds',
                'max_window_seconds', 'max_retained_bytes',
                'extra_host_reserve_bytes'}
    optional = {'max_cached_signatures', 'max_cached_bounds'}
    if not isinstance(config, dict) or set(config) - optional != required:
        raise ValueError('Explicit bounded source staging settings required')
    for key in ('records_per_step', 'max_steps', 'max_retained_bytes',
                'extra_host_reserve_bytes'):
        value = config[key]
        if type(value) is not int or value < 1:
            raise ValueError('Positive source staging ' + key + ' required')
    if (config['records_per_step'] < initial_chunk or
            config['records_per_step'] > 1_048_576 or
            config['max_steps'] > 32):
        raise ValueError('Source staging step geometry exceeds the admitted bounds')
    if config['max_retained_bytes'] > config['extra_host_reserve_bytes']:
        raise ValueError('Source staging retained ledger exceeds its host reserve')
    signatures=config.get('max_cached_signatures',1024)
    bounds=config.get('max_cached_bounds',64)
    if (type(signatures) is not int or not 1<=signatures<=32768 or
            type(bounds) is not int or not 1<=bounds<=256):
        raise ValueError('Bounded source staging header cache entries required')
    # Extra retained cache entries are not included in the source-segment
    # ledger cap. Reserve a deliberately loose Python-object envelope at
    # admission; this is not an operating-system peak-RSS guarantee.
    additional_cache_reserve=(max(0,signatures-1024)+max(0,bounds-64))*(16<<10)
    if config['max_retained_bytes']+additional_cache_reserve>config['extra_host_reserve_bytes']:
        raise ValueError('Source staging header cache exceeds its host reserve')
    for key in ('max_cpu_seconds', 'max_window_seconds'):
        value = config[key]
        if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                not 0 < value < float('inf')):
            raise ValueError('Positive finite source staging ' + key + ' required')
    return config


class ProductiveSourceStage:
    """Spend at most one bounded source step per completed output event."""

    def __init__(self, header_factory, source_identity, chunk_markers, config, *,
                 on_step_observation=None, on_complete=None):
        validate_source_stage_config(config, chunk_markers)
        if on_step_observation is not None and not callable(on_step_observation):
            raise ValueError('Source-stage observation callback must be callable')
        if on_complete is not None and not callable(on_complete):
            raise ValueError('Source-stage completion callback must be callable')
        self.on_step_observation = on_step_observation
        self.on_complete = on_complete
        self.header_factory = header_factory
        self.source_identity = source_identity
        self.chunk_markers = chunk_markers
        self.config = dict(config)
        self.lock = threading.Lock()
        self.worker = None
        self.pending = 0
        self.started = None
        self.steps = []
        self.cpu_seconds = 0.
        self.stage = None
        self.result = None
        self.stop_reason = None
        self.start_error = None
        self.completion_callback_error = None
        self.closed = False

    def output_written(self):
        with self.lock:
            if self.closed or self.stop_reason is not None or self.result is not None:
                return
            if self.started is None:
                self.started = time.perf_counter()
            self.pending += 1
            if self.worker is None or not self.worker.is_alive():
                try:
                    self.worker = threading.Thread(target=self._work,
                        name='torchgwas-source-stage', daemon=True)
                    self.worker.start()
                except (OSError, RuntimeError) as error:
                    self.start_error = type(error).__name__ + ': ' + str(error)
                    self.stop_reason = 'worker_start_error'
                    self.worker = None

    def _work(self):
        while True:
            with self.lock:
                if self.closed or self.stop_reason is not None or self.pending == 0:
                    self.worker = None
                    return
                if len(self.steps) >= self.config['max_steps']:
                    self.stop_reason = 'step_budget'
                    self.worker = None
                    return
                if self.cpu_seconds >= self.config['max_cpu_seconds']:
                    self.stop_reason = 'cpu_budget'
                    self.worker = None
                    return
                if time.perf_counter() - self.started >= self.config['max_window_seconds']:
                    self.stop_reason = 'window_budget'
                    self.worker = None
                    return
                self.pending -= 1
            wall_began = time.perf_counter()
            cpu_began = time.thread_time()
            error = None
            result = None
            retained = None
            observation = None
            try:
                if self.stage is None:
                    from .incremental_pgen_schedule import IncrementalPgenSchedule
                    header = self.header_factory()
                    if header.input_identity != self.source_identity:
                        raise ValueError('Staged PGEN identity differs from admission')
                    self.stage = IncrementalPgenSchedule(header, 0,
                        header._header.variant_ct, self.chunk_markers,
                        records_per_step=self.config['records_per_step'],
                        max_segments=self.config['max_steps'])
                progress = self.stage.advance()
                retained = (_retained_size(self.stage._primary_segments) +
                            _retained_size(self.steps))
                if retained > self.config['max_retained_bytes']:
                    raise ValueError('Source staging retained-ledger budget exceeded')
                if progress['complete']:
                    result = self.stage.finish()
                if self.on_step_observation is not None:
                    try:
                        observation = self.on_step_observation()
                    except Exception as caught:
                        observation = dict(status='observation_error',
                                           error=type(caught).__name__+': '+str(caught))
                    if retained + _retained_size(observation) > self.config['max_retained_bytes']:
                        observation = dict(status='omitted_retained_budget')
                    retained += _retained_size(observation)
            except Exception as caught:
                error = type(caught).__name__ + ': ' + str(caught)
            if retained is None and self.stage is not None:
                retained = _retained_size(self.stage._primary_segments)
            row = dict(wall_seconds=time.perf_counter() - wall_began,
                       cpu_seconds=time.thread_time() - cpu_began,
                       cursor=None if self.stage is None else self.stage.cursor,
                       retained_ledger_bytes=retained, error=error)
            if self.stage is not None:
                row['header_cache']=dict(signatures=self.stage.header.cache_info(),
                    bounds=self.stage.header.bounds_cache_info())
            if observation is not None:row['writer_observation']=observation
            with self.lock:
                self.steps.append(row)
                self.cpu_seconds += row['cpu_seconds']
                if result is not None:
                    self.result = result
                if error is not None:
                    self.stop_reason = 'source_error'
                elif self.result is not None:
                    self.stop_reason = 'complete'
                elif self.cpu_seconds >= self.config['max_cpu_seconds']:
                    self.stop_reason = 'cpu_budget'
                elif time.perf_counter() - self.started >= self.config['max_window_seconds']:
                    self.stop_reason = 'window_budget'
            if result is not None and self.on_complete is not None:
                # Run outside the stage lock. A callback may inspect the
                # completed ledger or start a separate evidence worker.
                try:
                    self.on_complete()
                except Exception as caught:
                    with self.lock:
                        self.completion_callback_error = (
                            type(caught).__name__ + ': ' + str(caught))

    def completed_ledger(self):
        """Return exact source work with the full observed worker cost.

        The incremental schedule's own timer excludes header setup, the
        optional live output observation and finalization. A later JIT
        decision must charge these measured worker steps once per job.
        """
        with self.lock:
            if self.result is None or self.stage is None:
                raise ValueError('Completed productive source stage required')
            return self.stage, dict(
                input_identity=self.source_identity,
                cpu_seconds=self.cpu_seconds,
                wall_seconds=sum(row['wall_seconds'] for row in self.steps),
                steps=len(self.steps),
                scope='Complete post-output source-stage worker cost, including header setup, metadata calculation, live observation and finalization. Charge once per job before any JIT comparison.')

    def finish(self):
        with self.lock:
            self.closed = True
            self.pending = 0
            worker = self.worker
        if worker is not None and worker is not threading.current_thread():
            worker.join()
        with self.lock:
            if self.stop_reason is None and self.result is None:
                self.stop_reason = 'job_finished'
        return self.snapshot()

    def snapshot(self):
        with self.lock:
            return dict(enabled=True, closed=self.closed,
                complete=self.result is not None,
                stop_reason=self.stop_reason,
                start_error=self.start_error,
                completion_callback_error=self.completion_callback_error,
                steps=list(self.steps), pending_events=self.pending,
                cpu_seconds=self.cpu_seconds,
                step_wall_seconds=sum(row['wall_seconds'] for row in self.steps),
                elapsed_window_seconds=None if self.started is None else
                    time.perf_counter() - self.started,
                input_identity=self.source_identity,
                header_cache=(None if self.stage is None else dict(
                    signatures=self.stage.header.cache_info(),
                    bounds=self.stage.header.bounds_cache_info())),
                scope='Read-only PGEN source accounting after useful output. Earlier CPU and wall costs must be charged to any later JIT decision; no runtime completion or switch is authorized.')
