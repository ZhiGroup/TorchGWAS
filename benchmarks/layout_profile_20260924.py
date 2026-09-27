"""Where a significant-pair run spends its time, per layout, selection backend and fsync.

Usage: python layout_profile_20260924.py DATA OUT.jsonl [--layouts 1gpu vshards4 tiles4]
       [--backends host device] [--fsync 1 0] [--repeats 1]

Each case runs in a fresh process with TORCHGWAS_SCAN_PROFILE=1 on the
benchmark data (empirical_layout_bench_20260923.py make) and appends one JSON
line: API and executor seconds, process CPU, the writer's time, and each
shard's or tile's scan profile (wait for the next decoded chunk, copy-slot and
result waits, GPU compute, chunk count).
"""
import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ENV = dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0',
           TORCHGWAS_SCAN_PROFILE='1', NUMPY_MADVISE_HUGEPAGE='0', OMP_NUM_THREADS='4', MKL_NUM_THREADS='1',
           OPENBLAS_NUM_THREADS='1', OMP_WAIT_POLICY='PASSIVE', TORCHGWAS_SHARED_DECODE='1',
           TORCHGWAS_GPU_FANOUT='pcie')
FIELDS = ('setup_seconds', 'fetch_seconds', 'copy_wait_seconds', 'result_wait_seconds',
          'gpu_compute_milliseconds', 'chunks')


def child(data, out, layout, fsync):
    import numpy as np
    from torchgwas.api import run_linear_gwas
    from torchgwas.io import load_genotype
    traits = np.load(data/'phenotype.npy', mmap_mode='r').shape[1]
    count = int(''.join(ch for ch in layout if ch.isdigit()) or 1)
    devices = [f'cuda:{i}' for i in range(count)]
    kwargs = dict(reduce='significant', significance_threshold=1e-5, chunk_size=1024, prefetch_chunks=4,
                  reader_workers=4*count, sumstats_fsync=bool(fsync))
    if layout.startswith('vshards'):
        kwargs['variant_devices'] = devices
    elif layout.startswith('tiles'):
        kwargs.update(trait_block=-(-traits//count), trait_devices=devices)
    else:
        kwargs['device'] = devices[0]
    genotype = load_genotype(str(data/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                             reader_workers=4*count)[0]
    started = time.perf_counter()
    result = run_linear_gwas(genotype=genotype, phenotype=data/'phenotype.npy', covariates=data/'covariates.npy',
                             compute_dtype='float32', output_dir=out, **kwargs)
    api = time.perf_counter()-started
    meta = result.run_metadata
    write = meta.get('sumstats_write') or {}
    profile = getattr(genotype, '_last_scan_profile', {}) or {}
    parts = profile.get('shards') or profile.get('tiles') or [dict(profile=profile)]
    rows = [dict(device=p.get('device'), **{k: round(float((p.get('profile') or {}).get(k) or 0), 3) for k in FIELDS})
            for p in parts]
    print(json.dumps(dict(api_seconds=round(api, 2), executor_seconds=write.get('setup_scan_and_write_seconds'),
                          write_seconds=write.get('write_seconds'), ordering=write.get('ordering'),
                          rows=meta.get('n_result_rows'), process_cpu=[round(x, 1) for x in os.times()[:2]],
                          scans=rows)), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--layouts', nargs='+', default=['1gpu', 'vshards4', 'tiles4'])
    parser.add_argument('--backends', nargs='+', default=['host', 'device'])
    parser.add_argument('--fsync', type=int, nargs='+', default=[1, 0])
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--child', nargs=3)
    args = parser.parse_args()
    if args.child:
        layout, backend, fsync = args.child
        child(args.data, args.out, layout, int(fsync))
        return
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cases = list(itertools.product(range(args.repeats), args.layouts, args.backends, args.fsync))
    for repeat, layout, backend, fsync in cases:
        work = args.data/'outputs'/'layout_profile'/f'{layout}_{backend}_f{fsync}_r{repeat}'
        done = subprocess.run([sys.executable, __file__, str(args.data), str(work), '--child', layout, backend, str(fsync)],
                              env={**os.environ, **ENV, 'TORCHGWAS_SIGNIFICANCE_BACKEND': backend},
                              check=True, capture_output=True, text=True)
        line = [l for l in done.stdout.splitlines() if l.startswith('{')][-1]
        record = dict(json.loads(line), layout=layout, backend=backend, fsync=fsync, repeat=repeat)
        with args.out.open('a') as handle:
            handle.write(json.dumps(record)+'\n')
        scans = record['scans']
        print(f"{layout:9s} {backend:6s} fsync={fsync} r{repeat}: api {record['api_seconds']:6.1f}  exec "
              f"{record['executor_seconds'] or -1:6.1f}  write {record['write_seconds'] or -1:6.1f}  cpu {record['process_cpu']}  "
              f"fetch/GPU-ms per scan {[(s['fetch_seconds'], round(s['gpu_compute_milliseconds'])) for s in scans]}",
              flush=True)


if __name__ == '__main__':
    main()
