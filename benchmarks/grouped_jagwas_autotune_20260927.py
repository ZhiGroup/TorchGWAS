"""Grouped JAGWAS (22 groups over the K = 8,192 panel): autotune priced by groups vs as one panel.

Configs, fresh interleaved processes:
- fixed4: 4 variant shards, chunk 1024;
- auto_grouped: autotune (groups priced per group, empirical_autotune.jagwas_factor_bytes /
  jagwas_projection_flops);
- auto_panel: autotune with group_sizes withheld from the planner, i.e. the
  earlier pricing of one panel of the total width.
Rows keep the layout decision and the executor seconds.

    python benchmarks/grouped_jagwas_autotune_20260927.py --root results/grouped_jagwas_autotune_20260927 \\
        --data /data/zxie3/torchgwas_bench/full_scale_20260925 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
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

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402

GROUPS = 22


def groups(traits):
    edges = np.linspace(0, traits, GROUPS + 1).astype(int)
    return [(f'g{i:02d}', list(range(int(edges[i]), int(edges[i + 1])))) for i in range(GROUPS)]


def configs(devices, traits):
    base = dict(reduce='jagwas', jagwas_groups=groups(traits))
    return [dict(name='fixed4', kwargs=dict(base, chunk_size=1024, prefetch_chunks=4,
                                            reader_workers=4 * len(devices), variant_devices=devices)),
            dict(name='auto_grouped', kwargs=dict(base, autotune=True, autotune_options=dict(devices=devices))),
            dict(name='auto_panel', kwargs=dict(base, autotune=True, autotune_options=dict(devices=devices)),
                 panel_pricing=True)]


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
    traits = int(json.loads((args.data / 'manifest.json').read_text())['traits'])
    rows = configs(args.devices, traits)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        if config.get('panel_pricing'):
            # The earlier pricing: factor and projection priced as one panel of
            # the total width (the planner looks these up at call time).
            import torchgwas.empirical_autotune as planner
            factor, projection = planner.jagwas_factor_bytes, planner.jagwas_projection_flops
            planner.jagwas_factor_bytes = lambda n_traits, group_sizes=None: factor(n_traits)
            planner.jagwas_projection_flops = lambda n_traits, group_sizes=None: projection(n_traits)
        out = args.output_root / args.root.name / f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, config)
        layout = (row.get('autotune') or {}).get('layout') or {}
        summary = dict(name=args.name, repeat=args.repeat, executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       api_seconds=round(row['api_seconds'], 2), variant_devices=row.get('variant_devices'),
                       reader_workers=row.get('reader_workers'), why=layout.get('why'),
                       cpu_demand=layout.get('cpu_demand'), rows=row.get('rows'))
        bench.save(args.root / f'{args.name}_r{args.repeat}.json', dict(row, summary=summary))
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
