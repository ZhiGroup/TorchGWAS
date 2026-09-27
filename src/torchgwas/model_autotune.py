"""Chunk-size tuning from per-chunk measurements, adjusted for load, re-planned on drift.

Replaces the segment trials of EmpiricalChunkTuner (forward/reverse order,
~25% of the job) with per-chunk samples:

* Every chunk that completes on a device right after another chunk of the
  same size, once the switch has settled, gives one rate sample
  (rows / completion interval). Settling skips `depth` completions per
  device: the first new-size chunks were decoded while the GPU still ran the
  old size and complete back to back. A few samples per device per size suffice,
  so probing costs a few percent of the job instead of a quarter.
* Load is measured, not assumed to drift smoothly. Native reads report the
  reader's CPU time and scheduler wait; the wait fraction w is a per-chunk
  covariate. Sizes are compared at the same load through a fitted
  log(period) = alpha_size + lambda * w (lambda = 0 without data).
* Per-variant stage costs are linear in chunk size, but how the stages
  overlap is not predictable from them (docs/autotune_design_20260924.md,
  2.2). So the model is the observed period itself, period(c) = a + b*c,
  fitted on the measured sizes. It is reported, and used to explain a
  decision, not to pick a size that was never measured.
* The first probe visits only sizes above the start. With period(c) =
  a + b*c and a per-chunk cost a >= 0, per-row time a/c + b cannot fall as
  the chunk shrinks. The start size is then revisited: its warmup samples
  include the pipeline's start-up, which would bias the comparison against
  it (H100 JAGWAS: 54-63k rows/s per GPU at start-up, ~94k in steady state).
* After committing, the committed size's rate is tracked against its
  estimate. A sustained deviation re-probes a neighbouring size only when the
  fitted model says it could gain more than the margin. For a larger size the
  bound credits the whole observed slowdown to per-chunk cost. A smaller size
  qualifies only when the fit's intercept is negative beyond its error, i.e.
  per-row cost measurably grows with chunk size. A slowdown that scales every
  size alike (host load) changes no ranking and re-probes nothing; on the
  H100, load drift triggered the maximum three re-probes in 4 of 5 runs.

Every candidate's memory was checked before the run (plan_layout) and the
ring is allocated for the largest, so switching never allocates. Changing the
size never changes results (AlignedChunkSizeControl).
"""
from __future__ import annotations

import math
import random
from collections import deque
import statistics
import threading
import time

from .adaptive_chunks import AlignedChunkSizeControl, aligned_chunk_sizes
from .empirical_autotune import PREFERRED_CHUNK_SIZE, _process_cpu


def _median_error(values):
    """Robust standard error of a median: 1.253 * 1.4826 * MAD / sqrt(n); 0 for fewer than 3 values."""
    values = list(values)
    if len(values) < 3:
        return 0.0
    middle = statistics.median(values)
    mad = statistics.median(abs(v - middle) for v in values)
    return 1.253 * 1.4826 * mad / math.sqrt(len(values))


def _rank_z(sample, reference):
    """Mann-Whitney z of `sample` against `reference` (average ranks for ties); positive = sample larger."""
    values = sorted([(v, 0) for v in sample] + [(v, 1) for v in reference])
    ranks, i = {}, 0
    rank_sum = 0.0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[j + 1][0] == values[i][0]:
            j += 1
        average = (i + j) / 2 + 1
        rank_sum += average * sum(1 for k in range(i, j + 1) if values[k][1] == 0)
        i = j + 1
    n1, n2 = len(sample), len(reference)
    u = rank_sum - n1 * (n1 + 1) / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    return (u - n1 * n2 / 2) / sd if sd else 0.0


def _load_fraction(observation):
    cpu = getattr(observation, 'read_cpu_seconds', None)
    wait = getattr(observation, 'read_runnable_wait_seconds', None)
    if cpu is None or wait is None or cpu + wait <= 0:
        return None
    return wait / (cpu + wait)


def _fit_load(samples_by_size):
    """lambda in log(period) = alpha_size + lambda * w, from samples that carry w."""
    rows = [(size, math.log(size / rate), w) for size, samples in samples_by_size.items()
            for rate, w in samples if w is not None and rate > 0]
    if len(rows) < 6:
        return 0.0
    # Within-size centring removes alpha_size.
    by_size = {}
    for size, y, w in rows:
        by_size.setdefault(size, []).append((y, w))
    num = den = 0.0
    for values in by_size.values():
        my = statistics.fmean(y for y, _ in values)
        mw = statistics.fmean(w for _, w in values)
        num += sum((y - my) * (w - mw) for y, w in values)
        den += sum((w - mw) ** 2 for _, w in values)
    return num / den if den > 1e-6 else 0.0


class ModelChunkTuner:
    """Chunk-size observer for one job; `control` is the matching selector.

    Same interface as EmpiricalChunkTuner: pass `tuner.control` as every scan's
    `_chunk_size_selector` and the tuner as its `_chunk_observer`.
    """
    accepts_minimal_observations = True

    def __init__(self, sizes, *, total_rows, initial=None, depth=4, concurrent=1,
                 warmup_fraction=0.02, probe_chunks=4, settle_chunks=1, margin=0.03,
                 max_probe_share=0.3, residual_threshold=0.15, residual_chunks=8,
                 max_reprobes=3, min_job_seconds=20.0, seed=20260924,
                 clock=time.perf_counter, cpu_clock=_process_cpu):
        self.sizes = aligned_chunk_sizes(sizes)
        if isinstance(total_rows, bool) or not isinstance(total_rows, int) or total_rows < 1:
            raise ValueError('Positive total source rows required')
        if not 0 <= warmup_fraction < 1 or not 0 < max_probe_share < 1 or probe_chunks < 2:
            raise ValueError('Invalid tuning fractions or probe length')
        if initial is None:
            initial = PREFERRED_CHUNK_SIZE if PREFERRED_CHUNK_SIZE in self.sizes else self.sizes[(len(self.sizes)-1)//2]
        if initial not in self.sizes:
            raise ValueError('Initial chunk size must be a candidate')
        self.total_rows, self.initial = total_rows, initial
        self.depth, self.concurrent = int(depth), int(concurrent)
        self.probe_chunks, self.settle_chunks = int(probe_chunks), int(settle_chunks)
        self.margin, self.max_probe_share = float(margin), float(max_probe_share)
        self.residual_threshold, self.residual_chunks = float(residual_threshold), int(residual_chunks)
        self.max_reprobes, self.min_job_seconds = int(max_reprobes), float(min_job_seconds)
        self.warmup_rows = max(int(warmup_fraction * total_rows), 2 * self.depth * self.concurrent * initial)
        self.control = AlignedChunkSizeControl(self.sizes, initial=initial)
        self.capacity = self.control.capacity
        self._clock, self._cpu, self._rng = clock, cpu_clock, random.Random(seed)
        self._lock = threading.Lock()
        self.state = 'fixed' if len(self.sizes) == 1 else 'warmup'
        self.reason = 'single_candidate' if len(self.sizes) == 1 else None
        self.current = initial
        self.choice = initial if self.state == 'fixed' else None
        self.completed_rows = 0
        self.devices = set()
        self.first_completed = self.last_completed = None
        self.estimated_job_seconds = None
        self._next_check = 0
        self._last = {}          # device -> (rows, completed)
        self._settle = {}        # device -> same-size completions still to skip
        self._visit = None       # the probe visit in progress
        self._queue = []         # sizes still to probe
        self._visits = {}        # size -> samples [(rate, w)] of its latest visit
        self._load_samples = {}  # size -> every sample, for the within-size load fit
        self.segments = []       # every probe visit, for the audit
        self.decisions = []
        self.reprobes = 0
        self.skipped_reprobes = []  # drift events the model gate declined
        self._baseline = None    # committed size's expected per-device rate at the reference load
        self._reference_load = None
        self._recent = deque(maxlen=4 * residual_chunks)  # committed-size samples (rate, w)
        self._ewma = None
        self._strikes = 0
        self.load_coefficient = 0.0
        self._errors = {}
        self._warm, self._reference = [], []
        self.model = None
        self.started = clock()

    # -- observer protocol -------------------------------------------------
    def for_scan(self, device, variant_range, n_traits, **_):
        return self

    def observer(self, device, variant_range=None, trait_range=None, **_):
        return self

    def __call__(self, observation):
        rows = int(observation.end - observation.start)
        now = observation.completed
        device = observation.device
        load = _load_fraction(observation)
        with self._lock:
            self.completed_rows += rows
            self.devices.add(device)
            self.first_completed = now if self.first_completed is None else self.first_completed
            self.last_completed = now
            previous = self._last.get(device)
            self._last[device] = (rows, now)
            sample = None
            if rows == self.current:
                if self._settle.get(device, 0) > 0:
                    self._settle[device] -= 1
                elif previous is not None and previous[0] == rows and now > previous[1]:
                    sample = (rows / (now - previous[1]), load)
            if self.state in ('warmup', 'deferred'):
                if sample is not None:
                    self._visits.setdefault(self.current, []).append(sample)
                    self._load_samples.setdefault(self.current, []).append(sample)
                self._maybe_start()
            elif self.state == 'probe' and sample is not None:
                self._visit['samples'].append(sample)
                if self._clearly_worse(self._visit):
                    # Time, not rows, is what a probe costs: a size several
                    # times slower (chunk 512 at full scale on a loaded A100:
                    # 2.5k vs 23-28k rows/s per GPU) spends ~20 s on its visit.
                    self._visit['stopped_early'] = True
                    self._finish_visit()
                elif len(self._visit['samples']) >= self.probe_chunks * max(1, len(self.devices)):
                    self._finish_visit()
            elif self.state == 'committed' and sample is not None:
                self._monitor(sample)

    # -- internals ---------------------------------------------------------
    def _remaining(self):
        return self.total_rows - self.completed_rows

    def _probe_cost(self, sizes, start):
        """Rows spent probing `sizes` in order from `start`: settle + samples + drain."""
        cost, previous = 0, start
        for size in sizes:
            cost += (self._settle_count() + self.probe_chunks + 1) * size * self.concurrent
            cost += self.depth * self.concurrent * previous
            previous = size
        return cost

    def _maybe_start(self):
        if self.completed_rows < max(self.warmup_rows, self._next_check):
            return
        elapsed = self.last_completed - self.first_completed
        rate = self.completed_rows / elapsed if elapsed > 0 else None
        if self.state == 'warmup':
            self.estimated_job_seconds = None if rate is None else self.total_rows / rate
        remaining_seconds = None if rate is None else self._remaining() / rate
        # Larger sizes only (module docstring); the start size is revisited last.
        others = [size for size in self.sizes if size > self.current]
        if not others:
            self._commit(self.current, 'no_larger_size_to_probe')
            self.state = 'skipped'
            return
        # Keep the budget: drop the largest sizes first.
        while others and (self._probe_cost(others + [self.current], self.current)
                          > self.max_probe_share * self._remaining()):
            others.pop()
        if not others:
            self._commit(self.current, 'job_too_short_to_probe')
            self.state = 'skipped'
            return
        if remaining_seconds is None or remaining_seconds < self.min_job_seconds:
            # Early chunks can run far faster than later ones; keep checking.
            self.state, self.reason = 'deferred', 'remaining job shorter than min_job_seconds'
            self._next_check = self.completed_rows + max(1, self.total_rows // 100)
            return
        self._rng.shuffle(others)  # random order: no assumption about how load moves
        self._queue = others + [self.current]
        self.state = 'probe'
        # The warmup at the initial size is that size's first visit.
        self.segments.append(dict(size=self.current, kind='warmup', rows_at_start=0,
                                  samples=list(self._visits.get(self.current, []))))
        self._next_visit()

    def _clearly_worse(self, visit):
        """True once the visited size is slower than the incumbent beyond three standard errors and the margin."""
        samples = [rate for rate, _ in visit['samples']]
        incumbent = self.initial if self.choice is None else self.choice
        reference = [rate for rate, _ in self._visits.get(incumbent, [])]
        if visit['size'] == incumbent or len(samples) < max(3, len(self.devices)) or len(reference) < 3:
            return False
        error = math.hypot(_median_error(samples), _median_error(reference))
        return statistics.median(samples) + 3 * error < statistics.median(reference) * (1 - self.margin)

    def _next_visit(self):
        size = self._queue.pop(0)
        self._switch(size)
        incumbent = self.initial if self.choice is None else self.choice
        kind = 'revisit' if size == incumbent else 'probe' if self.reprobes == 0 else 'reprobe'
        self._visit = dict(size=size, kind=kind, rows_at_start=self.completed_rows, samples=[])

    def _switch(self, size):
        self.current = size
        self.control.set_size(size)
        # The first new-size chunks were decoded while the GPU still ran the
        # old size, so up to `depth` of them complete back to back. H100
        # JAGWAS: the first size probed after warmup read 178-248k rows/s per
        # GPU against ~92k at its revisit and in steady state. Skip them.
        self._settle = {device: self._settle_count() for device in self.devices}

    def _settle_count(self):
        return max(self.settle_chunks, self.depth)

    def _finish_visit(self):
        visit = self._visit
        # Estimates use a size's latest visit (the start size: its revisit,
        # not its start-up); the load fit keeps every sample of every size.
        self._visits[visit['size']] = list(visit['samples'])
        self._load_samples.setdefault(visit['size'], []).extend(visit['samples'])
        self.segments.append(visit)
        self._visit = None
        if self._queue:
            self._next_visit()
        else:
            self._decide()

    def _adjusted(self, rate, w):
        """A rate moved to the reference load, with the fitted load coefficient."""
        if w is None or self._reference_load is None or self.load_coefficient == 0.0:
            return rate
        return rate * math.exp(self.load_coefficient * (w - self._reference_load))

    def _estimates(self):
        """Per-device rate per size, at the median observed load."""
        fitted = _fit_load(self._load_samples)
        # Too few loaded samples in a re-probe: keep the earlier coefficient.
        self.load_coefficient = fitted if fitted != 0.0 else self.load_coefficient
        loads = [w for samples in self._visits.values() for _, w in samples if w is not None]
        reference = statistics.median(loads) if loads else None
        self._reference_load = reference
        estimates = {}
        self._errors = {}
        for size, samples in self._visits.items():
            if not samples:
                continue
            if reference is None or self.load_coefficient == 0.0:
                adjusted = [rate for rate, _ in samples]
            else:
                adjusted = [rate * math.exp(self.load_coefficient * ((w if w is not None else reference) - reference))
                            for rate, w in samples]
            estimates[size] = statistics.median(adjusted)
            self._errors[size] = _median_error(adjusted)
        # The observed period model. Decisions use measured sizes only; the
        # fit bounds what a drift re-probe could gain (_reprobe_candidates).
        points = sorted((size, size / rate) for size, rate in estimates.items() if rate > 0)
        if len(points) >= 2:
            xs, ys = zip(*points)
            mx, my = statistics.fmean(xs), statistics.fmean(ys)
            sxx = sum((x - mx) ** 2 for x in xs)
            b = sum((x - mx) * (y - my) for x, y in points) / sxx if sxx else 0.0
            # Intercept error from each period's error (period = size / rate):
            # a = sum(l_i * y_i) with l_i = 1/n - mx * (x_i - mx) / sxx.
            n = len(points)
            a_error = math.sqrt(sum(
                ((1 / n - (mx * (x - mx) / sxx if sxx else 0.0))
                 * x * self._errors.get(x, 0.0) / estimates[x] ** 2) ** 2
                for x, _ in points))
            self.model = dict(a_seconds=my - b * mx, b_seconds_per_row=b, a_error_seconds=a_error,
                              points=[list(p) for p in points])
        return estimates

    def _decide(self):
        estimates = self._estimates()
        # First decision: the starting size; after a drift re-probe: the committed one.
        incumbent = self.initial if self.choice is None else self.choice
        best = max(estimates, key=estimates.get)
        noise = 0.0
        if incumbent in estimates:
            # Switch only on a gain the samples resolve: above the margin and
            # one standard error of the difference of medians (~84% one-sided).
            # Full-scale A100 under load: per-size estimates moved 3x between visits.
            noise = math.hypot(self._errors.get(best, 0.0), self._errors.get(incumbent, 0.0))
        if incumbent in estimates and estimates[best] - estimates[incumbent] <= max(
                self.margin * estimates[incumbent], noise):
            reason = ('within_noise_of_incumbent' if best != incumbent and
                      estimates[best] >= estimates[incumbent] * (1 + self.margin) else 'within_margin_of_incumbent')
            best = incumbent
        else:
            reason = 'highest_measured_throughput'
        self.decisions.append(dict(rows=self.completed_rows, estimates={str(k): v for k, v in estimates.items()},
                                   errors={str(k): v for k, v in self._errors.items()},
                                   load_coefficient=self.load_coefficient, choice=best, reason=reason))
        self._commit(best, reason)
        # The drift baseline comes from the committed size's own first
        # samples (_monitor), not the few probe samples behind the decision.

    def _commit(self, size, reason):
        self.choice, self.reason = size, reason
        if size != self.current:
            self._switch(size)
        self.state = 'committed'
        self._ewma, self._strikes, self._baseline = None, 0, None
        self._warm = []
        self._recent.clear()
        self.decided_at_rows = self.completed_rows

    def _monitor(self, sample):
        self._recent.append(sample)
        # Load the model already explains is not drift.
        rate = self._adjusted(*sample)
        self._warm.append(rate)
        span = 2 * self.residual_chunks
        if self._baseline is None:
            if len(self._warm) >= span:
                self._reference = list(self._warm)
                self._baseline = self._ewma = statistics.median(self._reference)
                self._warm = []
            return
        # Non-overlapping windows of 2*residual_chunks samples, each compared
        # with the committed size's own first samples by a rank test: per-chunk
        # rates are skewed and heavy-tailed on a busy host, so a mean or a
        # normal-theory error of the median misjudges them. A window strikes
        # when the shift is both significant (|z| > 3) and larger than the
        # fixed threshold; two strikes in a row are drift.
        if len(self._warm) < span:
            return
        window, self._warm = self._warm, []
        self._ewma = statistics.median(window)
        drift = self._ewma / self._baseline - 1
        z = _rank_z(window, self._reference)
        self._strikes = self._strikes + 1 if abs(z) > 3 and abs(drift) > self.residual_threshold else 0
        if self._strikes < 2 or self.reprobes >= self.max_reprobes:
            return
        neighbours, bounds = self._reprobe_candidates()
        elapsed = self.last_completed - self.first_completed
        remaining_seconds = self._remaining() * elapsed / self.completed_rows if elapsed > 0 else None
        declined = ('no neighbour can gain the margin' if not neighbours else
                    'probe budget' if self._probe_cost(neighbours, self.current) > self.max_probe_share * self._remaining()
                    else 'remaining job shorter than min_job_seconds'
                    if remaining_seconds is None or remaining_seconds < self.min_job_seconds else None)
        if declined is not None:
            self._strikes = 0
            self.skipped_reprobes.append(dict(rows=self.completed_rows, size=self.current, drift=drift,
                                              reason=declined, gain_bounds=bounds))
            return
        # Conditions changed: the committed size's recent samples are its new estimate.
        self.reprobes += 1
        recent = list(self._recent)[-2 * self.residual_chunks:]
        self._visits = {self.current: recent}
        self.segments.append(dict(size=self.current, kind='drift', rows_at_start=self.completed_rows,
                                  samples=recent, drift=drift))
        self._queue = neighbours
        self.state = 'probe'
        self._next_visit()

    def _reprobe_candidates(self):
        """Neighbouring sizes the model says could beat the committed one by the margin.

        Per-row time is r(c) = a/c + b. The fit (a, b) is anchored at the
        committed size's baseline, and the observed slowdown dP is credited to
        a for a larger neighbour (a per-chunk cost grew) and to b for a
        smaller one. Both are the most favourable readings for switching, so
        a neighbour that cannot clear the margin even then is not probed.
        """
        size = self.current
        index = self.sizes.index(size)
        if not self._baseline or not self._ewma:
            return [], {}
        before, now = size / self._baseline, size / self._ewma  # seconds per chunk per device
        change = now - before
        model = self.model or {}
        a = float(model.get('a_seconds', 0.0))
        error = float(model.get('a_error_seconds', 0.0))
        per_row = now / size
        chosen, bounds = [], {}
        for neighbour in (self.sizes[i] for i in (index - 1, index + 1) if 0 <= i < len(self.sizes)):
            if neighbour > size:
                intercept = a + max(change, 0.0)
            elif a + error < 0:
                intercept = a + min(change, 0.0)
            else:
                bounds[str(neighbour)] = None  # a >= 0 within error: smaller cannot win
                continue
            gain = intercept * (1 / size - 1 / neighbour) / per_row
            bounds[str(neighbour)] = gain
            if gain > self.margin:
                chosen.append(neighbour)
        return chosen, bounds

    def audit(self):
        with self._lock:
            segments = []
            for visit in self.segments:
                rates = [rate for rate, _ in visit['samples']]
                segments.append(dict(size=visit['size'], kind=visit['kind'], rows_at_start=visit['rows_at_start'],
                                     samples=len(rates), rows_per_second=statistics.median(rates) if rates else None,
                                     loads=[w for _, w in visit['samples']], other_size_rows=0))
            return dict(method='model_per_chunk', sizes=list(self.sizes), initial=self.initial,
                        capacity=self.capacity, state=self.state, reason=self.reason, choice=self.choice,
                        total_rows=self.total_rows, warmup_rows=self.warmup_rows,
                        probe_chunks=self.probe_chunks, margin=self.margin,
                        estimated_job_seconds=self.estimated_job_seconds, min_job_seconds=self.min_job_seconds,
                        decided_fraction=(None if getattr(self, 'decided_at_rows', None) is None
                                          else self.decided_at_rows / self.total_rows),
                        load_coefficient=self.load_coefficient, model=self.model, reprobes=self.reprobes,
                        skipped_reprobes=list(self.skipped_reprobes),
                        decisions=list(self.decisions), segments=segments, devices=sorted(self.devices),
                        first_chunk_seconds=None if self.first_completed is None else self.first_completed - self.started,
                        last_chunk_seconds=None if self.last_completed is None else self.last_completed - self.started)
