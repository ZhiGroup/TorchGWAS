"""Bounded independent CPU probes distributed over productive job callbacks.

The caller runs each advance inside its charged productive planning step and
owns the probe buffers. A probe reports CPU service for fixed work, never a
loaded pipeline interval. Drift thresholds are heuristics, not confidence
bounds. A completed record still requires price binding before model use.
"""
from copy import deepcopy
import math
from statistics import median
import time

from .calibration_cache import CalibrationParameterCache, _age, _key, read_calibration_record


SCHEMA = 'torchgwas.cpu_service.v1'


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError('Positive finite '+name+' required')
    return float(value)


class CpuServiceRefresh:
    """Reuse an original record after a short check, or collect a new window.

    Construction performs no cache I/O or measurement. Each advance invokes
    sample() at most samples_per_step times. sample() must itself be bounded;
    this controller cannot preempt it. No partial or unstable window is saved.
    Matching checks do not renew the original record's age. Check samples can
    contribute to a replacement window, without repeating their work.
    """
    def __init__(self, directory, name, *, dependencies, work_units, max_age_seconds,
                 samples_per_step=2, check_samples=2, refresh_samples=7,
                 drift_ratio=2., max_sample_ratio=4., absolute_cpu_seconds=2e-6):
        self._key = _key('cpu_capacity', name, dependencies)
        if not isinstance(self._key['dependencies'].get('measurement_protocol'), dict) or not self._key['dependencies']['measurement_protocol']:
            raise ValueError('Explicit measurement protocol required')
        self.work_units = _positive(work_units, 'fixed work units')
        self.max_age_seconds = _age(max_age_seconds)
        for value in (samples_per_step, check_samples, refresh_samples):
            if type(value) is not int or not 1 <= value <= 128:
                raise ValueError('Sample counts must be integers in [1,128]')
        if not 2 <= check_samples <= refresh_samples or refresh_samples < 3:
            raise ValueError('At least two check samples and three refresh samples required')
        self.samples_per_step, self.check_samples, self.refresh_samples = samples_per_step, check_samples, refresh_samples
        self.drift_ratio = _positive(drift_ratio, 'drift ratio')
        self.max_sample_ratio = _positive(max_sample_ratio, 'sample ratio')
        if min(self.drift_ratio, self.max_sample_ratio) <= 1:
            raise ValueError('Ratio thresholds must exceed one')
        self.absolute_cpu_seconds = _positive(absolute_cpu_seconds, 'absolute CPU tolerance')
        self.cache = CalibrationParameterCache(directory)
        self._state = 'pending'
        self._rows = []
        self._cached = self._result = self._check = self._lookup_reason = None
        self._steps = 0

    def _row(self, row):
        if not isinstance(row, dict) or set(row) != {'cpu_seconds', 'wall_seconds', 'work_units', 'observed_unix_seconds'}:
            raise ValueError('Complete independent CPU sample required')
        row = deepcopy(row)
        for key in row:
            _positive(row[key], key)
        if row['work_units'] != self.work_units:
            raise ValueError('Probe work differs from the fixed measurement protocol')
        return row

    def _spread(self, rows):
        values = [row['cpu_seconds'] for row in rows]
        lo, hi = min(values), max(values)
        return dict(ratio=hi/lo, minimum_cpu_seconds=lo, maximum_cpu_seconds=hi,
                    stable=hi-lo <= self.absolute_cpu_seconds or hi <= self.max_sample_ratio*lo)

    def _window_stability(self, rows):
        """Check dispersion and drift across nonoverlapping window halves.

        A bounded overall spread does not imply a stable measurement: costs
        can move from one level to another while remaining inside that bound.
        Use the existing drift tolerance for early/recent medians too. This is
        a heuristic, not a confidence interval or a future-capacity guarantee.
        """
        if len(rows)<2:
            return dict(assessed=False,stable=False,samples=len(rows))
        half=len(rows)//2
        early=median(row['cpu_seconds'] for row in rows[:half])
        recent=median(row['cpu_seconds'] for row in rows[-half:])
        drift=(abs(recent-early)>self.absolute_cpu_seconds and
               (recent>self.drift_ratio*early or early>self.drift_ratio*recent))
        spread=self._spread(rows)
        return dict(assessed=True,stable=spread['stable'] and not drift,
            samples=len(rows),samples_per_half=half,spread=spread,
            early_median_cpu_seconds=early,recent_median_cpu_seconds=recent,
            recent_over_early=recent/early,temporal_drift=drift,
            method='nonoverlapping_first_last_half_medians')

    def _value(self, rows):
        return dict(schema=SCHEMA, work_units=self.work_units, samples=deepcopy(rows),
                    cpu_seconds_per_unit=median(row['cpu_seconds']/self.work_units for row in rows))

    def _read(self, path):
        return read_calibration_record(path, **self._key, max_age_seconds=self.max_age_seconds)

    def _load(self):
        hit = self.cache.lookup(**self._key, max_age_seconds=self.max_age_seconds)
        self._lookup_reason = hit.get('reason')
        if hit['hit']:
            try:
                exact = self._read(hit['path'])
                record = exact['record']; value = record['value']
                rows = [self._row(row) for row in value['samples']]
                if (len(rows) != self.refresh_samples or value != self._value(rows)
                        or record['observed_unix_seconds'] != min(row['observed_unix_seconds'] for row in rows)
                        or any(row['observed_unix_seconds'] > record['created_unix_seconds'] for row in rows)
                        or not self._window_stability(rows)['stable']):
                    raise ValueError('Incompatible or unstable independent CPU window')
                self._cached = dict(path=hit['path'], **exact)
            except (ValueError, KeyError, TypeError, OSError):
                self._lookup_reason = 'invalid_or_expired_cpu_window'
        self._state = 'checking' if self._cached else 'measuring'

    def advance(self, sample, *, validate=None):
        """One charged batch. Failures are terminal and propagate to the budget."""
        if self._state in ('ready', 'unstable_measurement', 'expired_measurement', 'failed'):
            return self.snapshot()
        try:
            self._steps += 1
            if self._state == 'pending':
                self._load()
            target = self.check_samples if self._state == 'checking' else self.refresh_samples
            for _ in range(min(self.samples_per_step, target-len(self._rows))):
                began = time.time()
                row = self._row(sample())
                ended = time.time()
                if not began <= row['observed_unix_seconds'] <= ended:
                    raise ValueError('Probe must report its actual observation time within this call')
                self._rows.append(row)
            # The caller can reject a context change during sampling before
            # either reuse or publication. A failed check saves no new record.
            if validate is not None:validate()
            if self._state == 'checking' and len(self._rows) == self.check_samples:
                old = self._cached['record']['value']['cpu_seconds_per_unit']
                fresh = self._value(self._rows)['cpu_seconds_per_unit']
                drift = (abs(fresh-old)*self.work_units > self.absolute_cpu_seconds
                         and (fresh > self.drift_ratio*old or old > self.drift_ratio*fresh))
                self._check = dict(cached_cpu_seconds_per_unit=old, fresh_cpu_seconds_per_unit=fresh,
                    ratio=fresh/old, drift=drift, spread=self._spread(self._rows),
                    window_stability=self._window_stability(self._rows),
                    original_record_sha256=self._cached['record_sha256'])
                try:
                    exact = self._read(self._cached['path'])
                except (ValueError, OSError):
                    self._check['status'] = 'expired_or_invalid_during_check'
                else:
                    if not drift and self._check['window_stability']['stable']:
                        self._check['status'] = 'consistent'
                        self._result = dict(path=self._cached['path'], **exact, status='reused_original')
                        self._state = 'ready'
                        return self.snapshot()
                    self._check['status'] = 'drift' if drift else 'unstable_check'
                self._state = 'measuring'
                # Deliberately return: the replacement continues after a later
                # useful chunk, even if this call had spare sample capacity.
                return self.snapshot()
            if self._state == 'measuring' and len(self._rows) == self.refresh_samples:
                if not self._window_stability(self._rows)['stable']:
                    self._state = 'unstable_measurement'
                    return self.snapshot()
                observed = min(row['observed_unix_seconds'] for row in self._rows)
                if time.time()-observed >= self.max_age_seconds:
                    self._state = 'expired_measurement'
                    return self.snapshot()
                value = self._value(self._rows)
                published = self.cache.store(**self._key, value=value,
                    provenance=dict(protocol=self._key['dependencies']['measurement_protocol'],
                        samples_per_step=self.samples_per_step, check=deepcopy(self._check),
                        window_stability=self._window_stability(self._rows),
                        scope='Fixed-work independent CPU service; no loaded-stage capacity or prediction guarantee.'),
                    max_age_seconds=self.max_age_seconds, observed_unix_seconds=observed)
                # Publication can take long enough to expire a short lifetime.
                try:
                    exact = self._read(published['path'])
                except ValueError:
                    self._state = 'expired_measurement'
                    return self.snapshot()
                self._result = dict(path=published['path'], **exact, status='measured_and_published')
                self._state = 'ready'
            return self.snapshot()
        except Exception:
            self._state = 'failed'
            raise

    def snapshot(self):
        return deepcopy(dict(state=self._state, steps=self._steps, samples=self._rows,
            lookup_reason=self._lookup_reason, check=self._check, result=self._result,
            window_stability=self._window_stability(self._rows),
            policy=dict(samples_per_step=self.samples_per_step, check_samples=self.check_samples,
                refresh_samples=self.refresh_samples, drift_ratio=self.drift_ratio,
                max_sample_ratio=self.max_sample_ratio, absolute_cpu_seconds=self.absolute_cpu_seconds),
            scope='Heuristic early check of one fixed-work CPU coefficient. Original age is never renewed.'))
