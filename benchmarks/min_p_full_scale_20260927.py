"""reduce='min-p' at full scale (22,250 x 8.09M hard calls): layouts and a check.

Stages:
- check: on the first `--check-variants` variants, the min-p store equals the
  dense store's per-variant minimum p (trait, t, -log10 P), and one GPU
  equals two variant shards row for row;
- observe: fresh interleaved processes per config, `--repeats` rounds;
  rows keep the layout decision and the executor seconds.

    python benchmarks/min_p_full_scale_20260927.py observe --root results/min_p_full_scale_20260927 \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
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

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def configs(devices):
    base = dict(reduce='min-p', prefetch_chunks=4)
    rows = [dict(name='fixed1', kwargs=dict(base, device=devices[0], chunk_size=4096, reader_workers=8)),
            dict(name='fixed2', kwargs=dict(base, chunk_size=4096, reader_workers=16, variant_devices=devices[:2]))]
    if len(devices) >= 4:
        rows.append(dict(name='fixed4', kwargs=dict(base, chunk_size=4096, reader_workers=16,
                                                    variant_devices=devices[:4])))
    rows.append(dict(name='auto', kwargs=dict(reduce='min-p', autotune=True,
                                              autotune_options=dict(devices=devices))))
    return rows


def read_min_p(directory):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    manifest, parts = open_indexed_sumstats(directory/'sumstats')
    parts = list(parts)
    values = {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}
    order = np.argsort(values['variant_index'], kind='stable')
    return manifest, {key: value[order] for key, value in values.items()}


def check(args):
    """Dense versus min-p on a variant range, and one GPU versus two shards."""
    import torch  # noqa: F401
    from torchgwas.sumstats import open_binary_sumstats
    out = args.output_root/args.root.name/'check'
    span = (0, int(args.check_variants))
    common = dict(chunk_size=4096, prefetch_chunks=4, variant_range=span)
    summary = {}
    for name, kwargs in (('dense', dict(common, device=args.devices[0], reader_workers=8)),
                         ('minp1', dict(common, device=args.devices[0], reader_workers=8, reduce='min-p')),
                         ('minp2', dict(common, reader_workers=16, reduce='min-p', variant_devices=args.devices[:2]))):
        row = bench.child(args.data, out/name, dict(name=name, kwargs=kwargs))
        summary[name] = dict(api_seconds=round(row['api_seconds'], 2), executor_seconds=row['executor_seconds'])
    _beta, t, logp, _ = open_binary_sumstats(out/'dense/sumstats')
    t, logp = np.asarray(t), np.asarray(logp, dtype=np.float64)
    scores = np.where(np.isfinite(logp), logp, -np.inf)
    winner = scores.argmax(axis=1)
    rows = np.arange(t.shape[0])
    manifest, one = read_min_p(out/'minp1')
    _, two = read_min_p(out/'minp2')
    valid = np.isfinite(scores.max(axis=1))
    np.testing.assert_array_equal(one['variant_index'], rows[valid])
    same = one['trait_index'] == winner[valid]
    # Where the dense store's float32 -log10 P ties or nearly ties, float32
    # rounding decides the dense argmax; compare the value there instead.
    np.testing.assert_allclose(one['neg_log10_p'], scores[rows[valid], winner[valid]], rtol=2e-6, atol=1e-6)
    np.testing.assert_allclose(one['t_stat'][same], t[rows[valid][same], winner[valid][same]], rtol=0, atol=0)
    for key in one:
        np.testing.assert_array_equal(two[key], one[key])
    summary.update(variants=int(t.shape[0]), valid=int(valid.sum()), trait_agreement=int(same.sum()),
                   traits=int(t.shape[1]), reduction=manifest.get('reduction'), df_layout=manifest.get('df'),
                   max_logp=float(one['neg_log10_p'].max()))
    bench.save(args.root/'check.json', summary)
    shutil.rmtree(out)
    print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', nargs='?', default='observe', choices=['observe', 'child', 'check'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--check-variants', type=int, default=200_000)
    parser.add_argument('--name'); parser.add_argument('--repeat', type=int)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    if args.stage == 'check':
        os.environ.update(bench.ENVIRONMENT)
        check(args)
        return
    rows = configs(args.devices)
    if args.stage == 'child':
        config = next(row for row in rows if row['name'] == args.name)
        out = args.output_root/args.root.name/f'{args.name}_r{args.repeat}'
        row = bench.child(args.data, out, config)
        layout = (row.get('autotune') or {}).get('layout') or {}
        summary = dict(name=args.name, repeat=args.repeat,
                       executor_seconds=row['executor_seconds'] and round(row['executor_seconds'], 2),
                       api_seconds=round(row['api_seconds'], 2), variant_devices=row.get('variant_devices'),
                       trait_block=row.get('trait_block'), reader_workers=row.get('reader_workers'),
                       prefetch_chunks=row.get('prefetch_chunks'), why=layout.get('why'),
                       chunk=((row.get('autotune') or {}).get('chunk') or {}).get('choice'), rows=row.get('rows'))
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
                        '--name', name, '--repeat', str(repeat)],
                       check=True, env={**os.environ, **bench.ENVIRONMENT})


if __name__ == '__main__':
    main()
