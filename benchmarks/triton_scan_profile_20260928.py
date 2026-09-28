"""Where one GPU's scan time goes, per statistics backend (TORCHGWAS_SCAN_PROFILE).

A slice of the full-scale store (22,250 samples, `--variants` from 1M),
min-p with -log10 P, chunk 4096, on one device; run once per backend in a
fresh process (the backend is fixed at import):

    TORCHGWAS_STATS_BACKEND=torch python benchmarks/triton_scan_profile_20260928.py --device cuda:5 \\
        --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926
"""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--device', required=True)
    parser.add_argument('--variants', type=int, default=200_000)
    parser.add_argument('--readers', type=int, default=8)
    parser.add_argument('--mode', default='min-p', choices=['min-p', 'full'])
    args = parser.parse_args()
    os.environ['TORCHGWAS_SCAN_PROFILE'] = '1'
    os.environ.setdefault('TORCHGWAS_PGEN_BACKEND', 'native')
    from torchgwas.linear import linear_scan_streaming_chunks
    from torchgwas.min_p import MinPReduction
    from torchgwas.pgen import PgenGenotype
    source = PgenGenotype(args.data/'input.pgen', mode='hardcall', reader_workers=args.readers,
                          metadata_cache_dir=args.data/'metadata_cache')
    phenotype = np.load(args.data/'phenotype.npy').astype(np.float32)
    covariates = np.load(args.data/'covariates.npy').astype(np.float32)
    span = (1_000_000, 1_000_000 + args.variants)
    extra = dict(reduction=MinPReduction()) if args.mode == 'min-p' else {}
    for attempt in ('warm', 'timed'):
        started = time.perf_counter()
        iterator, _ = linear_scan_streaming_chunks(
            source, phenotype, covariates, chunk_size=4096, device=args.device, compute_dtype='float32',
            reader_workers=args.readers, prefetch_chunks=4, compute_p_values=False, compute_log10_p=True,
            log10_p_dtype='float32', variant_range=span if attempt == 'timed' else (span[0], span[0] + 20_000),
            borrow_results=True, **extra)
        for _ in iterator:
            pass
        seconds = time.perf_counter() - started
    profile = {key: (round(value, 4) if isinstance(value, float) else value)
               for key, value in getattr(source, '_last_scan_profile', {}).items()
               if isinstance(value, (int, float, str, bool))}
    print(json.dumps(dict(backend=profile.get('statistics_backend'), encoding=profile.get('native_encoding'),
                          seconds=round(seconds, 3), us_per_variant=round(1e6 * seconds / args.variants, 3),
                          profile=profile)), flush=True)


if __name__ == '__main__':
    main()
