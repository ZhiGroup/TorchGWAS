"""The post-scan tail: embedded variant IDs with NumPy's hugepage advice on or off.

Significant pairs on four variant shards (the fixed layout of
autotune_overhead_20260926), in fresh interleaved processes:
- embed_advice_on: sumstats_variant_ids=True, TORCHGWAS_NUMPY_HUGEPAGE=1 (NumPy's default);
- embed_advice_off: sumstats_variant_ids=True, advice off (the run default);
- index_only: the default store, advice off.
Each row keeps the writer's publication timings (sumstats_write.publication).

    python benchmarks/publication_tail_20260927.py --root results/publication_tail_20260927 \
        --data /data/zxie3/torchgwas_bench/full_scale_20260925 --devices cuda:4 cuda:5 cuda:6 cuda:7 \
        --output-root /dev/shm/zxie3_torchgwas_bench
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


def configs(devices):
    base = dict(reduce='significant', significance_threshold=1e-5, chunk_size=1024, prefetch_chunks=4,
                reader_workers=4 * len(devices), variant_devices=devices)
    return [dict(name='embed_advice_on', kwargs=dict(base, sumstats_variant_ids=True),
                 env={'TORCHGWAS_NUMPY_HUGEPAGE': '1'}),
            dict(name='embed_advice_off', kwargs=dict(base, sumstats_variant_ids=True), env={}),
            dict(name='index_only', kwargs=dict(base), env={})]


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
    rows = configs(args.devices)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        out = args.output_root / args.root.name / f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, config)
        with open('/proc/pressure/memory') as stream:
            pressure = stream.readline().strip()
        summary = dict(name=args.name, repeat=args.repeat, api_seconds=round(row['api_seconds'], 2),
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       publication=row['sumstats_write'].get('publication'),
                       ids_embedded=row['sumstats_write'].get('variant_ids_embedded'),
                       advice=row.get('numpy_hugepage_advice'), memory_pressure=pressure)
        bench.save(args.root / f'{args.name}_r{args.repeat}.json', dict(row, summary=summary, output=str(out)))
        shutil.rmtree(out)
        print(json.dumps(summary), flush=True)
        return
    args.root.mkdir(parents=True, exist_ok=True)
    schedule = []
    rng = random.Random(20260927)
    for repeat in range(args.repeats):
        names = [row['name'] for row in rows]
        rng.shuffle(names)
        schedule += [(name, repeat) for name in names]
    for name, repeat in schedule:
        if (args.root / f'{name}_r{repeat}.json').exists():
            continue
        subprocess.run([sys.executable, __file__, 'child', '--root', str(args.root), '--data', str(args.data),
                        '--devices', *args.devices, '--output-root', str(args.output_root),
                        '--name', name, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
