"""The Triton statistics at full scale: Torch, Triton on int8 rows, Triton on packed rows.

22,250 x 8.09M hard calls on four GPUs (fixed variant shards, chunk 4096,
16 readers, depth 4), `reduce='min-p'` and full output, K = 512 and the
dataset given. Three statistics paths, in fresh interleaved processes:
- torch: TORCHGWAS_STATS_BACKEND=torch (the old default), int8 transport;
- triton_int8: the Triton kernels on the same int8 transport;
- triton_packed: the Triton kernels on packed two-bit rows (the new default:
  a quarter of the host-to-device bytes).

    python benchmarks/triton_full_scale_20260928.py --root results/triton_full_scale_20260928 \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
        --output-root /data/zxie3/torchgwas_bench/triton_out
"""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402

PATHS = dict(torch=dict(TORCHGWAS_STATS_BACKEND='torch'),
             triton_int8=dict(TORCHGWAS_STATS_BACKEND='triton', TORCHGWAS_PGEN_PACKED='0'),
             # '' is neither 0 nor 1: the transport is chosen automatically.
             triton_packed=dict(TORCHGWAS_STATS_BACKEND='triton', TORCHGWAS_PGEN_PACKED=''))


def configs(args):
    rows = []
    for mode in args.modes:
        kwargs = dict(prefetch_chunks=4, chunk_size=4096, reader_workers=16, variant_devices=args.devices[:4])
        if mode != 'full':
            kwargs['reduce'] = mode
        for path, env in PATHS.items():
            rows.append(dict(name=f'{mode}_{path}', kwargs=kwargs, env=env))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', nargs='?', default='observe', choices=['observe', 'child'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--modes', nargs='+', default=['min-p', 'full'])
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    rows = configs(args)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        out = args.output_root/args.root.name/f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, config)
        summary = dict(name=args.name, repeat=args.repeat, api_seconds=round(row['api_seconds'], 2),
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       rows=row.get('rows'))
        bench.save(args.root/f'{args.name}_r{args.repeat}.json', dict(row, summary=summary))
        shutil.rmtree(out)
        print(json.dumps(summary), flush=True)
        return
    schedule = []
    rng = random.Random(20260928)
    for repeat in range(args.repeats):
        names = [row['name'] for row in rows]
        rng.shuffle(names)
        schedule += [(name, repeat) for name in names]
    for name, repeat in schedule:
        if (args.root/f'{name}_r{repeat}.json').exists():
            continue
        subprocess.run([sys.executable, __file__, 'child', '--root', str(args.root), '--data', str(args.data),
                        '--devices', *args.devices, '--output-root', str(args.output_root),
                        '--modes', *args.modes, '--name', name, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
