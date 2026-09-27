"""Where a host-bound scan's time goes: readers, ring depth, pinning and the tuner.

Full scale (22,250 x 8.09M hard calls), reduce='min-p' or dense output,
fresh interleaved processes, executor seconds. Configs vary one factor at a
time from autotune's own choice (2 shards, 16 readers per GPU, depth 32,
chunk 4096) toward the fixed layout that beat it (depth 4):

- r4d4, r8d8, r16d16, r16d32: one GPU, chunk 4096, readers x depth;
- s2_r16d32, s2_r16d16, s2_r8d8, s2_r4d4: two shards at the given per-GPU
  readers and depth, fixed chunk 4096 (autotune's layout without its tuner);
- pin: cudaHostAlloc seconds for one depth-32 ring of 4096-row slots, alone.

    python benchmarks/host_stage_attribution_20260927.py observe --root results/host_stage_attribution_20260927 \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --devices cuda:4 cuda:5 \\
        --output-root /data/zxie3/torchgwas_bench/min_p_out --reduce min-p
"""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def configs(devices, reduce):
    base = dict(chunk_size=4096, **({} if reduce == 'none' else dict(reduce=reduce)))
    rows = [dict(name=f'r{r}d{d}', kwargs=dict(base, device=devices[0], reader_workers=r, prefetch_chunks=d))
            for r, d in ((4, 4), (8, 8), (16, 16), (16, 32))]
    rows += [dict(name=f's2_r{r}d{d}', kwargs=dict(base, variant_devices=devices[:2], reader_workers=2 * r,
                                                   prefetch_chunks=d))
             for r, d in ((16, 32), (16, 16), (8, 8), (4, 4))]
    return rows


def pin_seconds(device, slots=32, rows=4096, row_bytes=22_250):
    import torch
    torch.cuda.init()
    with torch.cuda.device(device):
        started = time.perf_counter()
        ring = [torch.empty(rows * row_bytes, dtype=torch.uint8, pin_memory=True) for _ in range(slots)]
        seconds = time.perf_counter() - started
    del ring
    return dict(slots=slots, bytes=slots * rows * row_bytes, seconds=round(seconds, 3),
                gb_per_second=round(slots * rows * row_bytes / seconds / 1e9, 2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', nargs='?', default='observe', choices=['observe', 'child', 'pin'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--reduce', default='min-p', choices=['min-p', 'none'])
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    if args.stage == 'pin':
        row = pin_seconds(args.devices[0])
        bench.save(args.root/f'pin_r{args.repeat or 0}.json', row)
        print(json.dumps(row), flush=True)
        return
    rows = configs(args.devices, args.reduce)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        out = args.output_root/args.root.name/f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, config)
        summary = dict(name=args.name, repeat=args.repeat, reduce=args.reduce,
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       api_seconds=round(row['api_seconds'], 2), variant_devices=row.get('variant_devices'),
                       reader_workers=row.get('reader_workers'), prefetch_chunks=row.get('prefetch_chunks'),
                       process_cpu=row.get('process_cpu'))
        bench.save(args.root/f'{args.name}_r{args.repeat}.json', dict(row, summary=summary))
        shutil.rmtree(out)
        print(json.dumps(summary), flush=True)
        return
    schedule = []
    rng = random.Random(20260927)
    for repeat in range(args.repeats):
        names = [row['name'] for row in rows] + ['pin']
        rng.shuffle(names)
        schedule += [(name, repeat) for name in names]
    for name, repeat in schedule:
        if (args.root/f'{name}_r{repeat}.json').exists():
            continue
        stage = ['pin'] if name == 'pin' else ['child', '--name', name]
        subprocess.run([sys.executable, __file__, *stage, '--root', str(args.root), '--data', str(args.data),
                        '--devices', *args.devices, '--output-root', str(args.output_root),
                        '--reduce', args.reduce, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
