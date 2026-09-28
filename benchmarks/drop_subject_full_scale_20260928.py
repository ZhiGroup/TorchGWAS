"""missing_phenotype='drop_subject' at full scale: what dropping subjects costs.

The same scan (22,250 x 8.09M hard calls, the default int8 transport and
Torch statistics, reduce='min-p', four GPUs) three ways, in fresh
interleaved processes:
- complete: the panel as stored;
- drop_whole: `--drop` of the subjects given a missing value, so they are
  dropped; rows travel whole and the GPU masks them (the default);
- drop_host: the same drop decoded on the host by the reader's selection
  (TORCHGWAS_PGEN_WHOLE_ROWS=0).

`prepare` writes the dropped panel beside the dataset's own files (links to
the genotype, covariates and metadata cache; a phenotype with NaN rows).

    python benchmarks/drop_subject_full_scale_20260928.py prepare --root results/drop_subject_full_scale_20260928 \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
        --output-root /data/zxie3/torchgwas_bench/drop_subject_out
    python benchmarks/drop_subject_full_scale_20260928.py observe ...same arguments...
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


def dropped_data(args):
    return args.output_root/f'{args.data.name}_drop{args.drop}'


def configs(args):
    kwargs = dict(reduce='min-p', prefetch_chunks=4, chunk_size=4096, reader_workers=16,
                  variant_devices=args.devices[:4])
    return [dict(name='complete', data=args.data, kwargs=kwargs),
            dict(name='drop_whole', data=dropped_data(args), kwargs=kwargs),
            dict(name='drop_host', data=dropped_data(args), kwargs=kwargs,
                 env=dict(TORCHGWAS_PGEN_WHOLE_ROWS='0'))]


def prepare(args):
    target = dropped_data(args)
    target.mkdir(parents=True, exist_ok=True)
    for name in ('input.pgen', 'input.pvar', 'input.psam', 'covariates.npy', 'manifest.json', 'metadata_cache'):
        link = target/name
        if not link.exists():
            link.symlink_to(args.data/name)
    phenotype = np.load(args.data/'phenotype.npy')
    rng = np.random.default_rng(20260928)
    rows = np.sort(rng.choice(phenotype.shape[0], int(round(args.drop * phenotype.shape[0])), replace=False))
    columns = rng.integers(0, phenotype.shape[1], size=rows.size)
    phenotype = phenotype.copy()
    phenotype[rows, columns] = np.nan
    np.save(target/'phenotype.npy', phenotype)
    summary = dict(dropped=int(rows.size), samples=int(phenotype.shape[0]), traits=int(phenotype.shape[1]))
    bench.save(args.root/'prepare.json', summary)
    print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'observe', 'child'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--drop', type=float, default=0.01)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    if args.stage == 'prepare':
        prepare(args)
        return
    rows = configs(args)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        out = args.output_root/args.root.name/f'{args.name}_r{args.repeat}'
        row = bench.child(config['data'], out, config)
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
                        '--drop', str(args.drop), '--name', name, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
