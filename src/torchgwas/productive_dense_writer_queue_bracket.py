"""Common-instant bounds from two passes over live dense writer counters."""
import math


def _count(value, name):
    if type(value) is not int or value < 0:
        raise ValueError('Nonnegative writer ' + name + ' required')
    return value


def _time(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError('Finite writer ' + name + ' required')
    return value


def bracket_dense_writer_queues(first, second, anchor_seconds):
    """Bound logical accepted-minus-written bytes at one time between passes.

    Accepted and written counters never decrease. Every first observation
    must finish before the anchor, and every second observation must begin
    after it. Writer callbacks and GPU work need not stop during sampling.
    """
    if (not isinstance(first, dict) or not first or
            not isinstance(second, dict) or set(first) != set(second) or
            type(anchor_seconds) not in (int, float) or
            not math.isfinite(anchor_seconds)):
        raise ValueError('Two bounded passes over the same writers required')
    writers = {}
    totals = [0, 0]
    for directory, earlier in first.items():
        later = second[directory]
        if (not isinstance(directory, str) or not directory or
                not isinstance(earlier, dict) or not isinstance(later, dict) or
                earlier.get('device') != later.get('device') or
                not isinstance(earlier.get('device'), str)):
            raise ValueError('Bound writer directory and device required')
        a, b = earlier.get('observation'), later.get('observation')
        if (not isinstance(a, dict) or not isinstance(b, dict) or
                any(row.get('kind') != 'torchgwas.dense_writer_queue_observation.v1' or
                    row.get('valid') is not True or
                    row.get('atomic_writer_streams') is not True or
                    not isinstance(row.get('streams'), dict)
                    for row in (a, b)) or
                not a['streams'] or 't_stat' not in a['streams'] or
                set(a['streams']) != set(b['streams'])):
            raise ValueError('Atomic writer snapshots must cover matching streams')
        if (_time(a.get('capture_started_seconds'), 'first start') >
                _time(a.get('capture_finished_seconds'), 'first finish') or
                _time(b.get('capture_started_seconds'), 'second start') >
                _time(b.get('capture_finished_seconds'), 'second finish') or
                a['capture_finished_seconds'] > anchor_seconds or
                b['capture_started_seconds'] < anchor_seconds):
            raise ValueError('Atomic writer snapshots must bracket the anchor')
        stream_rows = {}
        for name, old in a['streams'].items():
            new = b['streams'][name]
            if not isinstance(old, dict) or not isinstance(new, dict):
                raise ValueError('Writer stream state required')
            accepted0 = _count(old.get('accepted_bytes'), name + ' accepted')
            accepted1 = _count(new.get('accepted_bytes'), name + ' accepted')
            written0 = _count(old.get('written_bytes'), name + ' written')
            written1 = _count(new.get('written_bytes'), name + ' written')
            for state, accepted, written in ((old, accepted0, written0),
                                             (new, accepted1, written1)):
                staged = _count(state.get('staging_bytes'), name + ' staging')
                queued = _count(state.get('queued_bytes'), name + ' queued')
                active = _count(state.get('active_bytes'), name + ' active')
                if (state.get('error') is not False or
                        accepted-written != staged+queued+active or
                        state.get('pending_write_bytes_interval') !=
                        [staged+queued, accepted-written]):
                    raise ValueError('Writer stream byte invariant failed')
            if (accepted0 > accepted1 or written0 > written1 or
                    written0 > accepted0 or written1 > accepted1):
                raise ValueError('Writer byte counters moved backward')
            lower = max(0, accepted0 - written1)
            upper = accepted1 - written0
            stream_rows[name] = [lower, upper]
            totals[0] += lower
            totals[1] += upper
        writers[directory] = dict(device=earlier['device'],
                                  logical_pending_bytes_interval=[
                                      sum(row[0] for row in stream_rows.values()),
                                      sum(row[1] for row in stream_rows.values())],
                                  streams=stream_rows)
    return dict(kind='torchgwas.bracketed_dense_writer_queues.v1',
                anchor_seconds=anchor_seconds,
                writers=writers,
                logical_pending_bytes_interval=totals,
                os_write_bytes_upper=totals[1],
                scope='All writer logical accepted-minus-completely-written bytes at one common time between two sampling passes. Monotone counters bound this interval even while writers progress. An active os.write may already have partly completed, so the lower logical endpoint is not a lower bound on physical bytes still to write. Producer/GPU state, dirty writeback and final fsync are outside the bracket.')
