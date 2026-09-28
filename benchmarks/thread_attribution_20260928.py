"""Which host thread is the serial stage of a multi-GPU scan? py-spy over every thread.

At K = 512 four GPUs finish no sooner than one: with the Triton statistics a
GPU needs ~1.1 us per variant, and four took as long as one would (7.7 s
for 8.09M variants), so some stage serializes at ~3.9 ms per 4,096-variant
chunk. Every thread is sampled (not only the GIL holder) while the default
scan runs; a thread busy for the whole run is that stage.

    py-spy record --threads --idle -r 100 -f raw -o /tmp/threads.txt -- \\
        python benchmarks/thread_attribution_20260928.py run --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 \\
        --devices cuda:4 cuda:5 cuda:6 cuda:7 --out /data/zxie3/torchgwas_bench/thread_out
    python benchmarks/thread_attribution_20260928.py summarize /tmp/threads.txt
"""
import argparse
import collections
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).parent))

IDLE = ('wait', 'acquire', 'sleep', 'select', 'poll', 'epoll', '_wait_for_tstate_lock', 'get (queue')


def run(args):
    import empirical_layout_bench_20260923 as bench
    # The benchmark environment, but the default statistics and transport.
    os.environ.update({k: v for k, v in bench.ENVIRONMENT.items()
                       if k not in ('TORCHGWAS_NATIVE_STATS', 'TORCHGWAS_PGEN_PACKED')})
    for key in ('TORCHGWAS_NATIVE_STATS', 'TORCHGWAS_PGEN_PACKED'):
        os.environ.pop(key, None)
    from torchgwas.api import run_linear_gwas
    started = time.perf_counter()
    result = run_linear_gwas(genotype=str(args.data/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                             genotype_cache_dir=str(args.data/'metadata_cache'), phenotype=args.data/'phenotype.npy',
                             covariates=args.data/'covariates.npy', output_dir=args.out, compute_dtype='float32',
                             chunk_size=4096, prefetch_chunks=4, reader_workers=16, variant_devices=args.devices,
                             reduce=args.reduce)
    timing = result.run_metadata.get('sumstats_write') or {}
    print(json.dumps(dict(api_seconds=round(time.perf_counter() - started, 2),
                          executor_seconds=timing.get('setup_scan_and_write_seconds'))), flush=True)
    shutil.rmtree(args.out)


def summarize(path, top=8):
    """Per thread: share of samples busy (not in a wait), and its busiest frames."""
    busy, total, frames = collections.Counter(), collections.Counter(), collections.defaultdict(collections.Counter)
    for line in open(path):
        stack, _, count = line.rstrip().rpartition(' ')
        count = int(count)
        parts = stack.split(';')
        thread = re.sub(r'0x[0-9a-f]+', '', parts[0]).strip()
        total[thread] += count
        leaf = parts[-1] if len(parts) > 1 else ''
        if any(word in leaf for word in IDLE):
            continue
        busy[thread] += count
        mine = [f for f in parts[1:] if 'torchgwas' in f]
        frames[thread][(mine[-1] if mine else leaf)[:150]] += count
    samples = max(total.values())
    for thread, count in busy.most_common(12):
        print(f'{100 * count / samples:5.1f}% busy  {thread[:90]}  (of {total[thread]} samples)')
        for frame, n in frames[thread].most_common(top):
            print(f'        {100 * n / samples:5.1f}%  {frame}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['run', 'summarize'])
    parser.add_argument('path', nargs='?')
    parser.add_argument('--data', type=Path)
    parser.add_argument('--devices', nargs='+')
    parser.add_argument('--out', type=Path)
    parser.add_argument('--reduce', default='min-p')
    args = parser.parse_args()
    if args.stage == 'run':
        run(args)
    else:
        summarize(args.path)


if __name__ == '__main__':
    main()
