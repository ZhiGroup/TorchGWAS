"""min-p against significant pairs on the same layout, same window.

Both reduce a (chunk x K) statistic on the device and move almost nothing to
the host, so they should cost about the same. Fixed 4 variant shards, chunk
4096, fresh interleaved processes, executor seconds.

    python benchmarks/min_p_vs_significant_20260927.py --root results/min_p_vs_significant_20260927 \\
        --data /data/zxie3/torchgwas_bench/full_scale_20260925 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
        --output-root /data/zxie3/torchgwas_bench/min_p_out
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', nargs='?', default='observe', choices=['observe', 'child'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    common = dict(chunk_size=4096, prefetch_chunks=4, reader_workers=16, variant_devices=args.devices)
    rows = {'significant': dict(common, reduce='significant'), 'min-p': dict(common, reduce='min-p')}
    if args.stage == 'child':
        out = args.output_root/args.root.name/f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, dict(name=args.name, kwargs=rows[args.name]))
        summary = dict(name=args.name, repeat=args.repeat,
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       api_seconds=round(row['api_seconds'], 2), rows=row.get('rows'),
                       process_cpu=row.get('process_cpu'))
        bench.save(args.root/f'{args.name}_r{args.repeat}.json', dict(row, summary=summary))
        shutil.rmtree(out)
        print(json.dumps(summary), flush=True)
        return
    schedule = []
    rng = random.Random(20260927)
    for repeat in range(args.repeats):
        names = list(rows)
        rng.shuffle(names)
        schedule += [(name, repeat) for name in names]
    for name, repeat in schedule:
        if (args.root/f'{name}_r{repeat}.json').exists():
            continue
        subprocess.run([sys.executable, __file__, 'child', '--root', str(args.root), '--data', str(args.data),
                        '--devices', *args.devices, '--output-root', str(args.output_root),
                        '--name', name, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
