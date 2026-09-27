"""Autotune against fixed variant shards, after pricing the scan's device work.

Full scale (22,250 x 8.09M hard calls), fresh interleaved processes,
executor seconds; rows keep autotune's layout and chunk decisions. `--reduce`
is min-p, significant, jagwas or none (dense output, written under
--output-root).

    python benchmarks/autotune_vs_fixed_20260927.py --root results/autotune_vs_fixed_20260927/minp_k512 \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
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

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def configs(devices, reduce):
    base = {} if reduce == 'none' else dict(reduce=reduce)
    return [dict(name='fixed4', kwargs=dict(base, chunk_size=4096, prefetch_chunks=4, reader_workers=16,
                                            variant_devices=devices[:4])),
            dict(name='fixed4_r8d8', kwargs=dict(base, chunk_size=4096, prefetch_chunks=8, reader_workers=32,
                                                 variant_devices=devices[:4])),
            dict(name='auto', kwargs=dict(base, autotune=True, autotune_options=dict(devices=devices)))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', nargs='?', default='observe', choices=['observe', 'child'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--reduce', default='min-p', choices=['min-p', 'significant', 'jagwas', 'none'])
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    rows = configs(args.devices, args.reduce)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        out = args.output_root/args.root.name/f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, config)
        autotune = row.get('autotune') or {}
        layout, chunk = autotune.get('layout') or {}, autotune.get('chunk') or {}
        summary = dict(name=args.name, repeat=args.repeat, reduce=args.reduce,
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       api_seconds=round(row['api_seconds'], 2), variant_devices=row.get('variant_devices'),
                       reader_workers=row.get('reader_workers'), prefetch_chunks=row.get('prefetch_chunks'),
                       why=layout.get('why'), cpu_demand=layout.get('cpu_demand'),
                       shard_model=layout.get('shard_model'), chunk=chunk.get('choice'),
                       chunk_state=chunk.get('state'), chunk_reason=chunk.get('reason'))
        bench.save(args.root/f'{args.name}_r{args.repeat}.json', dict(row, summary=summary))
        shutil.rmtree(out)
        print(json.dumps(summary), flush=True)
        return
    schedule = []
    rng = random.Random(20260927)
    for repeat in range(args.repeats):
        names = [row['name'] for row in rows]
        rng.shuffle(names)
        schedule += [(name, repeat) for name in names]
    for name, repeat in schedule:
        if (args.root/f'{name}_r{repeat}.json').exists():
            continue
        subprocess.run([sys.executable, __file__, 'child', '--root', str(args.root), '--data', str(args.data),
                        '--devices', *args.devices, '--output-root', str(args.output_root),
                        '--reduce', args.reduce, '--name', name, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
