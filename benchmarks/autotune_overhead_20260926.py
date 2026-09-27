"""Where JAGWAS autotune loses time on a ~30 s job: layout setup or chunk probing.

Three configurations on the same GPUs, interleaved in fresh processes:
- fixed4: 4 variant shards at chunk 1024, no autotune;
- auto_chunk1024: autotune with one chunk size (layout decisions, no probing);
- auto: full autotune.
auto - auto_chunk1024 is the chunk-probing cost; auto_chunk1024 - fixed4 is
autotune's layout/pre-scan cost. Uses the layout bench's child(), so the
timings match its table (executor = setup, scan and write).

    python benchmarks/autotune_overhead_20260926.py --root results/autotune_overhead_20260926 \\
        --data /data/zxie3/torchgwas_bench/full_scale_20260925 --devices cuda:1 cuda:3 cuda:4 cuda:5
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


def configs(devices, reduce):
    if reduce == 'full':
        # Dense binary output (beta + t). fixed1 uses the 09-15 benchmark flags.
        return [dict(name='fixed1', kwargs=dict(device=devices[0], chunk_size=4096, reader_workers=16,
                                                prefetch_chunks=32)),
                dict(name='fixed4', kwargs=dict(variant_devices=devices, chunk_size=4096,
                                                reader_workers=8 * len(devices), prefetch_chunks=16)),
                # The shard count dense autotune chose on 2026-09-27, without its planning.
                dict(name='fixed2', kwargs=dict(variant_devices=devices[:2], chunk_size=4096,
                                                reader_workers=16, prefetch_chunks=16)),
                dict(name='auto', kwargs=dict(autotune=True, autotune_options=dict(devices=devices)))]
    base = dict(reduce=reduce)
    if reduce == 'significant':
        base.update(significance_threshold=1e-5)
    return [dict(name='fixed4', kwargs=dict(base, chunk_size=1024, prefetch_chunks=4,
                                            reader_workers=4 * len(devices), variant_devices=devices)),
            dict(name='auto_chunk1024', kwargs=dict(base, autotune=True,
                                                    autotune_options=dict(devices=devices, chunk_sizes=[1024]))),
            dict(name='auto', kwargs=dict(base, autotune=True, autotune_options=dict(devices=devices)))]


def dense_sample(out, rows=512):
    """t at evenly spaced variants: a small cross-layout check kept after the output is deleted."""
    from torchgwas.sumstats import open_binary_sumstats
    _, t_stat, _logp, _ = open_binary_sumstats(out / 'sumstats')
    n = t_stat.shape[0]
    return np.stack([np.asarray(t_stat[int(i)], dtype=np.float32)
                     for i in np.linspace(0, n - 1, min(rows, n)).astype(np.int64)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', nargs='?', default='observe', choices=['observe', 'child'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', default=['cuda:0', 'cuda:1', 'cuda:2', 'cuda:3'])
    parser.add_argument('--reduce', choices=['jagwas', 'significant', 'full'], default='jagwas')
    parser.add_argument('--keep-output', action='store_true', help='keep each run output (deleted by default)')
    parser.add_argument('--only', nargs='+', help='run only these config names')
    parser.add_argument('--output-root', type=Path, default=None,
                        help='write run outputs here (e.g. tmpfs) instead of <data>/outputs')
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    rows = configs(args.devices, args.reduce)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        env = {'TORCHGWAS_SIGNIFICANCE_BACKEND': 'device'} if args.reduce == 'significant' else {}
        out = (args.output_root or args.data / 'outputs') / args.root.name / f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, dict(config, env=env))
        chunk = (row.get('autotune') or {}).get('chunk') or {}
        summary = dict(name=args.name, repeat=args.repeat, api_seconds=round(row['api_seconds'], 2),
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       chunk=chunk.get('choice'), state=chunk.get('state'), reason=chunk.get('reason'),
                       reprobes=chunk.get('reprobes'), first_chunk_seconds=chunk.get('first_chunk_seconds'),
                       last_chunk_seconds=chunk.get('last_chunk_seconds'),
                       variant_devices=row.get('variant_devices'), reader_workers=row.get('reader_workers'),
                       trait_devices=row.get('trait_devices'), prefetch_chunks=row.get('prefetch_chunks'))
        if args.reduce == 'full':
            np.save(args.root / f'{args.name}_r{args.repeat}_tsample.npy', dense_sample(out))
        bench.save(args.root / f'{args.name}_r{args.repeat}.json', dict(row, summary=summary, output=str(out)))
        if not args.keep_output:
            shutil.rmtree(out)
        print(json.dumps(summary), flush=True)
        return
    args.root.mkdir(parents=True, exist_ok=True)
    # Warm the page cache once, as the layout bench does, so the first
    # scheduled run is not the only one reading cold input.
    manifest = json.loads((args.data / 'manifest.json').read_text())
    for name in manifest.get('genotype_files', [str(args.data / 'input.pgen')]) + [
            args.data / 'phenotype.npy', args.data / 'covariates.npy']:
        with open(name, 'rb') as stream:
            while stream.read(1 << 24):
                pass
    schedule = []
    rng = random.Random(20260926)
    for repeat in range(args.repeats):
        names = [row['name'] for row in rows]
        rng.shuffle(names)
        schedule += [(name, repeat) for name in names]
    for name, repeat in schedule:
        if (args.root / f'{name}_r{repeat}.json').exists() or (args.only and name not in args.only):
            continue
        subprocess.run([sys.executable, __file__, 'child', '--root', str(args.root), '--data', str(args.data),
                        '--devices', *args.devices, '--reduce', args.reduce, '--name', name,
                        '--repeat', str(repeat)] + (['--keep-output'] if args.keep_output else [])
                       + (['--output-root', str(args.output_root)] if args.output_root else []),
                       check=True, env={**os.environ, **bench.ENVIRONMENT})
    if args.reduce == 'full':
        samples = {path.stem[:-len('_tsample')]: np.load(path) for path in sorted(args.root.glob('*_tsample.npy'))}
        if samples:
            reference_name, reference = next(iter(samples.items()))
            for name, sample in samples.items():
                finite = np.isfinite(sample) & np.isfinite(reference)
                print(f'{name}: max |t - {reference_name}| on sampled rows',
                      float(np.max(np.abs(sample[finite] - reference[finite]), initial=0.0)),
                      'same NaN pattern', bool(np.array_equal(np.isnan(sample), np.isnan(reference))))


if __name__ == '__main__':
    main()
