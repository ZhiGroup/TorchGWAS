"""cProfile one significant-pair run: where the non-GPU time goes.

Usage: python run_cprofile_20260924.py DATA OUT.prof [--layout 1gpu|tiles4|vshards4] [--backend device|host]

Prints the top functions by cumulative and by own time. Threads other than
the caller (readers, writers, shard producers) are not profiled; their
effect shows as waits in the caller.
"""
import argparse
import cProfile
import os
import pstats
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--layout', default='1gpu')
    parser.add_argument('--backend', default='device')
    args = parser.parse_args()
    os.environ.update(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0',
                      NUMPY_MADVISE_HUGEPAGE='0', TORCHGWAS_SHARED_DECODE='1', TORCHGWAS_GPU_FANOUT='pcie',
                      TORCHGWAS_SIGNIFICANCE_BACKEND=args.backend)
    import numpy as np
    from torchgwas.api import run_linear_gwas
    traits = np.load(args.data/'phenotype.npy', mmap_mode='r').shape[1]
    count = int(''.join(c for c in args.layout if c.isdigit()) or 1)
    devices = [f'cuda:{i}' for i in range(count)]
    kwargs = dict(reduce='significant', significance_threshold=1e-5, chunk_size=1024, prefetch_chunks=4,
                  reader_workers=4*count)
    if args.layout.startswith('tiles'):
        kwargs.update(trait_block=-(-traits//count), trait_devices=devices)
    elif args.layout.startswith('vshards'):
        kwargs['variant_devices'] = devices
    else:
        kwargs['device'] = devices[0]
    out = args.data/'outputs'/'cprofile'/f'{args.layout}_{args.backend}'
    profiler = cProfile.Profile()
    started = time.perf_counter()
    profiler.enable()
    run_linear_gwas(genotype=str(args.data/'input.pgen'), phenotype=args.data/'phenotype.npy',
                    covariates=args.data/'covariates.npy', pgen_mode='hardcall', compute_dtype='float32',
                    output_dir=out, **kwargs)
    profiler.disable()
    print(f'wall {time.perf_counter()-started:.1f} s')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    profiler.dump_stats(str(args.out))
    stats = pstats.Stats(profiler)
    stats.sort_stats('cumulative').print_stats(35)
    stats.sort_stats('tottime').print_stats(20)


if __name__ == '__main__':
    main()
