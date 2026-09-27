"""Immutable calibration records with dependency and freshness checks.

Structural values have no age expiry, but only match identical dependencies.
Empirical values require an explicit maximum age. Availability observations
are audit-only: every new job must query live memory and contention again.
This cache stores evidence; a hit does not certify a model or a candidate.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import time


SCHEMA = 'torchgwas.calibration_parameters.v1'
STRUCTURAL = frozenset({'device_properties', 'source_work', 'kernel_geometry', 'workspace'})
EMPIRICAL = frozenset({'cpu_capacity', 'gpu_capacity', 'transfer_capacity',
                       'storage_capacity', 'stage_observations'})
LIVE = frozenset({'available_memory', 'contention'})


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _digest(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _age(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError('Positive finite maximum age required for empirical calibration')
    return float(value)


def _key(kind, name, dependencies):
    if kind not in STRUCTURAL | EMPIRICAL | LIVE:
        raise ValueError('Unknown calibration parameter kind')
    if not isinstance(name, str) or not name.strip():
        raise ValueError('Explicit calibration parameter name required')
    if not isinstance(dependencies, dict) or not dependencies:
        raise ValueError('Explicit nonempty calibration dependencies required')
    # Round trip also rejects non-JSON values and keeps later caller mutation
    # from changing the published binding. Callers must include measurement
    # protocol plus relevant source, library, GPU, shape, math/thread/I/O context.
    return json.loads(_json(dict(kind=kind, name=name, dependencies=dependencies)))


class CalibrationParameterCache:
    def __init__(self, directory):
        self.directory = Path(directory)

    def store(self, kind, name, value, *, dependencies, provenance, max_age_seconds=None,
              observed_unix_seconds=None):
        key = _key(kind, name, dependencies)
        if not isinstance(provenance, dict) or not provenance:
            raise ValueError('Explicit measurement provenance required')
        if kind in EMPIRICAL:
            lifetime = _age(max_age_seconds)
        else:
            if max_age_seconds is not None:
                raise ValueError('Only empirical calibration accepts a maximum age')
            lifetime = None
        record = dict(schema=SCHEMA, **key, value=value, provenance=provenance,
                      created_unix_seconds=time.time(), max_age_seconds=lifetime)
        if observed_unix_seconds is not None:
            if (kind not in EMPIRICAL or isinstance(observed_unix_seconds, bool)
                    or not isinstance(observed_unix_seconds, (int, float))
                    or not math.isfinite(observed_unix_seconds)
                    or observed_unix_seconds > record['created_unix_seconds']):
                raise ValueError('Empirical observation time must be finite and no later than publication')
            record['observed_unix_seconds'] = float(observed_unix_seconds)
        body = _json(record)
        identity = hashlib.sha256(body.encode()).hexdigest()
        root = self.directory / _digest(key)
        root.mkdir(parents=True, exist_ok=True)
        path = root / (identity + '.json')
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=root, delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(body+'\n'); stream.flush(); os.fsync(stream.fileno())
            try:
                os.link(temporary, path)  # Publication cannot overwrite another record.
            except FileExistsError:
                if path.read_text(encoding='utf-8') != body+'\n':
                    raise ValueError('Calibration cache identity collision')
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return dict(record_sha256=identity, path=str(path), kind=kind,
                    reusable=kind not in LIVE)

    def lookup(self, kind, name, *, dependencies, max_age_seconds=None):
        """Find the newest compatible observation still fresh after cache I/O.

        A directory can contain records from many jobs. Reading it is not
        instantaneous, and none of that time renews an observation's age.
        Consumers must still recheck evidence immediately before a decision.
        """
        key = _key(kind, name, dependencies)
        if kind in LIVE:
            return dict(hit=False, reason='requires_live_observation')
        # A caller can demand a shorter freshness interval than the producer,
        # but cannot extend the lifetime written into an immutable record.
        requested_age = _age(max_age_seconds) if max_age_seconds is not None else None
        if kind in STRUCTURAL and requested_age is not None:
            raise ValueError('Structural calibration uses dependency validity, not an age override')
        root = self.directory / _digest(key)
        started = time.time(); candidates = []; invalid = expired = future = 0
        for path in root.glob('*.json'):
            try:
                row = json.loads(path.read_text(encoding='utf-8'))
                if (set(row)-{'observed_unix_seconds'} != {'schema', 'kind', 'name', 'dependencies', 'value', 'provenance',
                                 'created_unix_seconds', 'max_age_seconds'}
                        or row['schema'] != SCHEMA or any(row[k] != v for k, v in key.items())
                        or _digest(row) != path.stem or not isinstance(row['provenance'], dict)
                        or not row['provenance']):
                    raise ValueError('Invalid calibration record')
                created = row['created_unix_seconds']
                if isinstance(created, bool) or not isinstance(created, (int, float)) or not math.isfinite(created):
                    raise ValueError('Invalid calibration timestamp')
                observed = row.get('observed_unix_seconds', created)
                if (isinstance(observed, bool) or not isinstance(observed, (int, float))
                        or not math.isfinite(observed) or observed > created
                        or ('observed_unix_seconds' in row and kind not in EMPIRICAL)):
                    raise ValueError('Invalid observation timestamp')
                lifetime = _age(row['max_age_seconds']) if kind in EMPIRICAL else None
                if kind in STRUCTURAL and row['max_age_seconds'] is not None:
                    raise ValueError('Invalid structural lifetime')
                if requested_age is not None:
                    lifetime = min(lifetime, requested_age)
                # Do not retain already unusable payloads from a long cache
                # history. Survivors still need the post-read check below.
                if created > started:
                    future += 1; continue
                if lifetime is not None and started-observed >= lifetime:
                    expired += 1; continue
                # Jobs can finish out of order. Publishing an old measurement
                # later must not displace a more recently observed capacity.
                candidates.append((observed, created, path.stem, row, lifetime))
            except (ValueError, TypeError, KeyError, OSError):
                invalid += 1
        # Use the post-read time for every candidate, including an older one
        # whose longer lifetime might outlast a more recent observation.
        now = time.time()
        if now < started:
            return dict(hit=False, reason='clock_moved_backwards', invalid_records=invalid,
                        expired_records=expired, future_records=future)
        fresh = []
        for observed, created, identity, row, lifetime in candidates:
            if created > now:
                future += 1
            elif lifetime is not None and now-observed >= lifetime:
                expired += 1
            else:
                fresh.append((observed, created, identity, row))
        if not fresh:
            return dict(hit=False, reason=('expired' if expired else 'future_timestamp' if future
                        else 'invalid_record' if invalid else 'missing_or_dependencies_changed'),
                        invalid_records=invalid, expired_records=expired, future_records=future)
        _, _, identity, record = max(fresh, key=lambda item: item[:3])
        return dict(hit=True, record_sha256=identity, path=str(root/(identity+'.json')), record=record,
                    age_seconds=now-record.get('observed_unix_seconds',record['created_unix_seconds']), invalid_records=invalid,
                    scope='Matching evidence only; live admission and early-job rate checks remain required.')


def inspect_calibration_record(path, *, kind, name, dependencies, max_age_seconds=None):
    """Validate immutable identity and report freshness, including expired data.

    Inspection is not permission to use a price. Call read_calibration_record
    for empirical model inputs; admission can inspect named pending refreshes.
    """
    key=_key(kind,name,dependencies)
    if kind not in EMPIRICAL:
        raise ValueError('Bound service prices require empirical calibration')
    requested_age=_age(max_age_seconds) if max_age_seconds is not None else None
    def unique(pairs):
        result={}
        for key,value in pairs:
            if key in result:raise ValueError('Duplicate calibration record key')
            result[key]=value
        return result
    path=Path(path)
    row=json.loads(path.read_text(encoding='utf-8'),object_pairs_hook=unique)
    if (not isinstance(row,dict) or set(row)!={'schema','kind','name','dependencies','value','provenance',
            'created_unix_seconds','observed_unix_seconds','max_age_seconds'}
            or row['schema']!=SCHEMA or any(row[k]!=v for k,v in key.items())
            or _digest(row)!=path.stem or not isinstance(row['provenance'],dict) or not row['provenance']):
        raise ValueError('Invalid or incompatible bound calibration record')
    now=time.time();created=row['created_unix_seconds'];observed=row['observed_unix_seconds']
    for value in (created,observed):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise ValueError('Invalid bound calibration observation time')
    if observed>created or created>now:
        raise ValueError('Future bound calibration observation/publication time')
    lifetime=_age(row['max_age_seconds'])
    if requested_age is not None:lifetime=min(lifetime,requested_age)
    return dict(record=row,record_sha256=path.stem,age_seconds=now-observed,
                fresh=now-observed<lifetime)


def read_calibration_record(path, *, kind, name, dependencies, max_age_seconds=None):
    """Read exactly one fresh immutable artifact; never substitute a new one."""
    result=inspect_calibration_record(path,kind=kind,name=name,dependencies=dependencies,
                                     max_age_seconds=max_age_seconds)
    if not result.pop('fresh'):
        raise ValueError('Expired bound calibration measurement; collect and bind a fresh record')
    return result
