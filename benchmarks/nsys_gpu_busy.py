"""GPU busy fraction, kernel/copy totals and host API time from an nsys sqlite export.

    nsys export --type sqlite -o run.sqlite run.nsys-rep
    python benchmarks/nsys_gpu_busy.py run.sqlite [--top 15]

The window runs from the first to the last GEMM kernel (the scan), so
start-up and output publication do not dilute the busy fraction.
"""
import argparse
import collections
import re
import sqlite3


def union_seconds(intervals):
    total, current_start, current_end = 0, None, None
    for start, end in sorted(intervals):
        if current_end is None or start > current_end:
            if current_end is not None:
                total += current_end - current_start
            current_start, current_end = start, end
        else:
            current_end = max(current_end, end)
    if current_end is not None:
        total += current_end - current_start
    return total / 1e9


def gaps(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [b[0] - a[1] for a, b in zip(merged, merged[1:])]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sqlite')
    parser.add_argument('--top', type=int, default=15)
    parser.add_argument('--gap-ms', type=float, default=1.0, help='list idle gaps longer than this')
    args = parser.parse_args()
    db = sqlite3.connect(args.sqlite)
    names = dict(db.execute('SELECT id, value FROM StringIds'))
    kernels = [(s, e, names.get(n, str(n)), stream, device, names.get(d, '')) for s, e, n, stream, device, d in db.execute(
        'SELECT start, end, shortName, streamId, deviceId, demangledName FROM CUPTI_ACTIVITY_KIND_KERNEL')]
    gemm = [k for k in kernels if 'gemm' in k[2].lower() or 'cutlass' in k[2].lower()]
    first, last = min(k[0] for k in gemm), max(k[1] for k in gemm)
    inside = lambda s, e: s >= first and e <= last
    window = (last - first) / 1e9
    kernels = [k for k in kernels if inside(k[0], k[1])]
    copies = [row for row in db.execute('SELECT start, end, bytes, copyKind, streamId FROM CUPTI_ACTIVITY_KIND_MEMCPY')
              if inside(row[0], row[1])]
    osrt = []
    try:
        osrt = [(s, e, names.get(n, str(n)), tid) for s, e, n, tid in db.execute(
            'SELECT start, end, nameId, globalTid FROM OSRT_API')]
    except sqlite3.OperationalError:
        pass
    try:
        sets = [row for row in db.execute('SELECT start, end FROM CUPTI_ACTIVITY_KIND_MEMSET') if inside(*row)]
    except sqlite3.OperationalError:
        sets = []
    kernel_iv = [(s, e) for s, e, *_ in kernels]
    copy_iv = [(s, e) for s, e, *_ in copies]
    print(f'scan window {window:.3f} s   GEMM kernels {len(gemm)}')
    print(f'busy (kernels) {union_seconds(kernel_iv):.3f} s  '
          f'busy (kernels+copies+memsets) {union_seconds(kernel_iv + copy_iv + sets):.3f} s')
    print(f'GEMM time {sum(e - s for s, e, *_ in gemm) / 1e9:.3f} s')
    idle = gaps(kernel_iv + copy_iv + sets)
    for limit in (10_000, 100_000, 1_000_000):
        chosen = [g for g in idle if g > limit]
        print(f'  idle gaps > {limit / 1000:g} us: {len(chosen)} totalling {sum(chosen) / 1e9:.3f} s')
    devices = sorted({k[4] for k in kernels})
    if len(devices) > 1:
        for device in devices:
            mine = [(s, e) for s, e, _, _, d, _ in kernels if d == device]
            gemm_s = sum(e - s for s, e, name, _, d, _ in kernels if d == device and ('gemm' in name.lower()))
            print(f'  device {device}: kernel busy {union_seconds(mine):.3f} s, GEMM {gemm_s / 1e9:.3f} s')
    by_stream = collections.defaultdict(lambda: [0, 0])
    for s, e, name, stream, device, _ in kernels:
        by_stream[(device, stream)][0] += 1
        by_stream[(device, stream)][1] += e - s
    print('kernels per (device, stream) (count, total s):')
    for key, (count, total) in sorted(by_stream.items()):
        print(f'  {key}: {count} {total / 1e9:.3f}')
    runtime = [(s, e, names.get(n, str(n)), tid) for s, e, n, tid in db.execute(
        'SELECT start, end, nameId, globalTid FROM CUPTI_ACTIVITY_KIND_RUNTIME')]
    merged = []
    for start, end in sorted(kernel_iv + copy_iv + sets):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    long_gaps = [(a[1], b[0]) for a, b in zip(merged, merged[1:]) if b[0] - a[1] > args.gap_ms * 1e6]
    print(f'idle gaps > {args.gap_ms:g} ms (offset ms, length ms; host calls overlapping >= 20% of the gap):')
    for start, end in long_gaps[:40]:
        busy = collections.defaultdict(int)
        for s, e, name, tid in runtime + osrt:
            overlap = min(e, end) - max(s, start)
            if overlap > 0.2 * (end - start):
                busy[(tid % 2**24, name)] += overlap
        calls = ', '.join(f'{tid}:{name} {overlap / 1e6:.1f}' for (tid, name), overlap in sorted(busy.items(), key=lambda item: -item[1])[:6])
        print(f'  {(start - first) / 1e6:9.1f} {(end - start) / 1e6:7.1f}  {calls}')
    by_kernel = collections.defaultdict(lambda: [0, 0])
    for s, e, name, stream, device, demangled in kernels:
        functor = re.findall(r'(?:at::native::|\(anonymous namespace\)::)(\w+)', demangled)
        key = (name, stream, '/'.join(dict.fromkeys(functor[:4])))
        by_kernel[key][0] += 1
        by_kernel[key][1] += e - s
    print('kernels (count, total s, mean us, stream, name / functors):')
    for (name, stream, functor), (count, total) in sorted(by_kernel.items(), key=lambda item: -item[1][1])[:args.top]:
        print(f'  {count:7d} {total / 1e9:8.3f} {total / count / 1e3:9.1f}  s{stream}  {name[:40]}  {functor[:110]}')
    kinds = {1: 'HtoD', 2: 'DtoH', 8: 'DtoD', 10: 'PtoP'}
    by_copy = collections.defaultdict(lambda: [0, 0, 0])
    for s, e, size, kind, _ in copies:
        entry = by_copy[kinds.get(kind, kind)]
        entry[0] += 1; entry[1] += size; entry[2] += e - s
    print('copies (count, MiB, total s):')
    for kind, (count, size, total) in by_copy.items():
        print(f'  {kind:5} {count:7d} {size / 2**20:10.1f} {total / 1e9:8.3f}')
    api = collections.defaultdict(lambda: [0, 0])
    threads = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    for s, e, n, tid in db.execute('SELECT start, end, nameId, globalTid FROM CUPTI_ACTIVITY_KIND_RUNTIME'):
        if not inside(s, e):
            continue
        name = names.get(n, str(n))
        api[name][0] += 1; api[name][1] += e - s
        threads[tid][name][0] += 1; threads[tid][name][1] += e - s
    print('runtime API (count, total s, mean us):')
    for name, (count, total) in sorted(api.items(), key=lambda item: -item[1][1])[:args.top]:
        print(f'  {count:8d} {total / 1e9:8.3f} {total / count / 1e3:9.1f}  {name}')
    print('runtime API per thread (total s):')
    for tid, calls in sorted(threads.items(), key=lambda item: -sum(v[1] for v in item[1].values())):
        top = sorted(calls.items(), key=lambda item: -item[1][1])[:5]
        print(f'  tid {tid % 2**24}: ' + ', '.join(f'{name} {count}x {total / 1e9:.2f}s' for name, (count, total) in top))


if __name__ == '__main__':
    main()
