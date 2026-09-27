"""Per-tile scan profile of shared decode with per-GPU PCIe copies versus NVLink fan-out.

Usage: python fanout_profile_20260924.py DATA_DIR OUT_DIR [--tiles 4] [--repeats 2]

Runs the benchmark's significant-pair job (DATA_DIR from
empirical_layout_bench_20260923.py make) with TORCHGWAS_SCAN_PROFILE=1 in a
fresh process per run, alternating modes, and prints each tile's wait for the
next chunk (fetch), copy-slot wait, result wait and GPU compute time.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ENV = dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0',
           TORCHGWAS_SCAN_PROFILE='1', TORCHGWAS_BLOCKING_EVENTS='1', NUMPY_MADVISE_HUGEPAGE='0',
           OMP_NUM_THREADS='4', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', TORCHGWAS_SHARED_DECODE='1')
FIELDS = ('fetch_seconds', 'copy_wait_seconds', 'result_wait_seconds', 'gpu_compute_milliseconds', 'chunks')


def child(data, out, tiles, mode):
    import numpy as np
    from torchgwas.api import run_linear_gwas
    traits = np.load(data/'phenotype.npy', mmap_mode='r').shape[1]
    started = time.perf_counter()
    result = run_linear_gwas(genotype=str(data/'input.pgen'), phenotype=data/'phenotype.npy',
                             covariates=data/'covariates.npy', pgen_mode='hardcall', compute_dtype='float32',
                             output_dir=out, reduce='significant', significance_threshold=1e-5,
                             trait_block=-(-traits//tiles), trait_devices=[f'cuda:{i}' for i in range(tiles)],
                             chunk_size=1024, reader_workers=4*tiles, prefetch_chunks=4)
    api = time.perf_counter()-started
    meta = result.run_metadata
    profile = None
    for value in meta.values():
        if isinstance(value, dict) and value.get('mode') == 'significant_trait_tiles':
            profile = value
    shared = meta.get('shared_decode') or {}
    rows = [dict(device=t['device'], **{k: round(t['profile'].get(k, 0), 3) for k in FIELDS})
            for t in (profile or {}).get('tiles', [])]
    print(json.dumps(dict(mode=mode, api_seconds=round(api, 2), transfer=shared.get('transfer'),
                          process_cpu=[round(x, 1) for x in os.times()[:2]], tiles=rows)), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--tiles', type=int, default=4)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--child', choices=['pcie', 'nvlink'])
    args = parser.parse_args()
    if args.child:
        child(args.data, args.out, args.tiles, args.child)
        return
    for repeat in range(args.repeats):
        for mode in (('pcie', 'nvlink') if repeat % 2 == 0 else ('nvlink', 'pcie')):
            out = args.out/f'{mode}_r{repeat}'
            subprocess.run([sys.executable, __file__, str(args.data), str(out), '--tiles', str(args.tiles),
                            '--child', mode], check=True, env={**os.environ, **ENV, 'TORCHGWAS_GPU_FANOUT': mode})


if __name__ == '__main__':
    main()
