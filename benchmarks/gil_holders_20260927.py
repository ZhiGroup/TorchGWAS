"""Who holds the GIL in a host-bound scan? One fixed layout, run under py-spy --gil.

Run this script under `py-spy record --gil -f raw -o out.txt -- python ...`
and pass the raw file to `--summarize`: samples are taken only from the
thread holding the GIL, so the summary is the GIL's time by function (the
innermost torchgwas frame, and by thread role).

    py-spy record --gil -r 200 -f raw -o /tmp/gil.txt -- python benchmarks/gil_holders_20260927.py run \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 --devices cuda:4 cuda:5 cuda:6 cuda:7 \\
        --out /data/zxie3/torchgwas_bench/gil_out --reduce min-p
    python benchmarks/gil_holders_20260927.py summarize /tmp/gil.txt
"""
import argparse
import collections
import json
import os
from pathlib import Path
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).parent))


def run(args):
    import empirical_layout_bench_20260923 as bench
    os.environ.update(bench.ENVIRONMENT)
    from torchgwas.api import run_linear_gwas
    options = {} if args.reduce == 'none' else dict(reduce=args.reduce)
    started = time.perf_counter()
    result = run_linear_gwas(genotype=str(args.data/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                             genotype_cache_dir=str(args.data/'metadata_cache'), phenotype=args.data/'phenotype.npy',
                             covariates=args.data/'covariates.npy', output_dir=args.out, compute_dtype='float32',
                             chunk_size=args.chunk, prefetch_chunks=args.depth,
                             reader_workers=args.readers * len(args.devices), variant_devices=args.devices, **options)
    timing = result.run_metadata.get('sumstats_write') or {}
    print(json.dumps(dict(api_seconds=round(time.perf_counter() - started, 2),
                          executor_seconds=timing.get('setup_scan_and_write_seconds'))), flush=True)
    shutil.rmtree(args.out)


def summarize(path, top=25):
    by_frame, by_thread, total = collections.Counter(), collections.Counter(), 0
    for line in open(path):
        stack, _, count = line.rstrip().rpartition(' ')
        count = int(count)
        frames = stack.split(';')
        total += count
        by_thread[frames[0].split(' ')[0] if frames else '?'] += count
        mine = [f for f in frames if 'torchgwas' in f]
        by_frame[mine[-1] if mine else frames[-1]] += count
    print(json.dumps(dict(total_samples=total)))
    for frame, count in by_frame.most_common(top):
        print(f'{100 * count / total:5.1f}%  {frame[:160]}')
    print('-- by thread')
    for thread, count in by_thread.most_common(12):
        print(f'{100 * count / total:5.1f}%  {thread[:120]}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['run', 'summarize'])
    parser.add_argument('path', nargs='?')
    parser.add_argument('--data', type=Path)
    parser.add_argument('--devices', nargs='+')
    parser.add_argument('--out', type=Path)
    parser.add_argument('--reduce', default='min-p')
    parser.add_argument('--chunk', type=int, default=4096)
    parser.add_argument('--depth', type=int, default=4)
    parser.add_argument('--readers', type=int, default=4)
    args = parser.parse_args()
    if args.stage == 'run':
        run(args)
    else:
        summarize(args.path)


if __name__ == '__main__':
    main()
