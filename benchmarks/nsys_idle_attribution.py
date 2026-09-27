"""Attribute each GPU's idle time to what its launching thread was doing.

    python benchmarks/nsys_idle_attribution.py run.sqlite

Per device: busy = union of its kernels and copies inside the scan window
(first to last GEMM on that device). Idle intervals are then intersected
with the launching thread's CUDA runtime calls and OS runtime calls
(sem_wait = Python lock/queue waits, pthread_cond_timedwait = GIL waits,
read/poll = I/O). Time in neither is the thread running Python or native
code without a traced call. CUDA calls take precedence over OS calls.
"""
import argparse
import bisect
import collections
import sqlite3


def merge(intervals):
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out


def subtract(intervals, holes):
    """intervals minus holes; both merged and sorted."""
    out, j = [], 0
    for start, end in intervals:
        cursor = start
        while j < len(holes) and holes[j][1] <= cursor:
            j += 1
        k = j
        while k < len(holes) and holes[k][0] < end:
            if holes[k][0] > cursor:
                out.append([cursor, holes[k][0]])
            cursor = max(cursor, holes[k][1])
            k += 1
        if cursor < end:
            out.append([cursor, end])
    return out


def overlap(intervals, start, end):
    """Total overlap of [start, end) with merged sorted intervals."""
    i = bisect.bisect_right(intervals, [start, float('inf')]) - 1
    i = max(i, 0)
    total = 0
    while i < len(intervals) and intervals[i][0] < end:
        total += max(0, min(end, intervals[i][1]) - max(start, intervals[i][0]))
        i += 1
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sqlite')
    parser.add_argument('--top', type=int, default=8)
    args = parser.parse_args()
    db = sqlite3.connect(args.sqlite)
    names = dict(db.execute('SELECT id, value FROM StringIds'))
    runtime = list(db.execute('SELECT start, end, nameId, globalTid, correlationId FROM CUPTI_ACTIVITY_KIND_RUNTIME'))
    launcher_of = {corr: tid for _, _, _, tid, corr in runtime}
    kernels = list(db.execute('SELECT start, end, deviceId, correlationId, shortName FROM CUPTI_ACTIVITY_KIND_KERNEL'))
    copies = list(db.execute('SELECT start, end, deviceId, correlationId, copyKind, bytes FROM CUPTI_ACTIVITY_KIND_MEMCPY'))
    try:
        osrt = list(db.execute('SELECT start, end, nameId, globalTid FROM OSRT_API'))
    except sqlite3.OperationalError:
        osrt = []
    by_thread_cuda = collections.defaultdict(list)
    for s, e, n, tid, _ in runtime:
        by_thread_cuda[tid].append((s, e, names.get(n, str(n))))
    by_thread_os = collections.defaultdict(list)
    for s, e, n, tid in osrt:
        by_thread_os[tid].append((s, e, names.get(n, str(n))))
    for device in sorted({k[2] for k in kernels}):
        mine = [k for k in kernels if k[2] == device]
        gemm = [k for k in mine if 'gemm' in names.get(k[4], '').lower()]
        if not gemm:
            continue
        first, last = min(k[0] for k in gemm), max(k[1] for k in gemm)
        launchers = collections.Counter(launcher_of.get(k[3]) for k in mine)
        tid = launchers.most_common(1)[0][0]
        busy = merge([(s, e) for s, e, *_ in mine if first <= s and e <= last] +
                     [(s, e) for s, e, d, *_ in copies if d == device and first <= s and e <= last])
        window = last - first
        idle = subtract([[first, last]], busy)
        idle_total = sum(e - s for s, e in idle)
        chunks = len(gemm)
        print(f'device {device}: window {window / 1e9:.3f} s, {chunks} GEMMs, busy {(window - idle_total) / 1e9:.3f} s, '
              f'idle {idle_total / 1e9:.3f} s ({idle_total / chunks / 1e6:.2f} ms/chunk); launcher tid {tid % 2**24}')
        cuda_calls = [(s, e, n) for s, e, n in by_thread_cuda[tid] if e > first and s < last]
        cuda_iv = merge([(s, e) for s, e, _ in cuda_calls])
        attribution = collections.Counter()
        for s, e, name in cuda_calls:
            attribution['cuda:' + name] += overlap(idle, s, e)
        remaining = subtract(idle, cuda_iv)
        os_calls = [(s, e, n) for s, e, n in by_thread_os[tid] if e > first and s < last]
        for s, e, name in os_calls:
            attribution['os:' + name] += overlap(remaining, s, e)
        accounted = sum(attribution.values())
        attribution['running (no traced call)'] = idle_total - accounted
        for name, value in attribution.most_common(args.top):
            print(f'    {value / 1e9:7.3f} s  {100 * value / idle_total:5.1f}%  {value / chunks / 1e6:6.2f} ms/chunk  {name}')
        calls = collections.Counter(n for _, _, n in cuda_calls)
        per_chunk = {name: round(count / chunks, 1) for name, count in calls.most_common(4)}
        print(f'    CUDA calls per chunk: {per_chunk}')


if __name__ == '__main__':
    main()
