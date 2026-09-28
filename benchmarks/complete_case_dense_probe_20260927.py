"""Dense output with missing phenotypes: one GPU, a variant range, per-stage profile.

Runs dense output on `--variants` variants of the complete and the missing
K = 512 panels (TORCHGWAS_SCAN_PROFILE=1) and prints the API seconds and the
scan's own timings, to locate host work the complete-case path adds.

    python benchmarks/complete_case_dense_probe_20260927.py --root /data/zxie3/torchgwas_bench \\
        --device cuda:4 --out /data/zxie3/torchgwas_bench/cc_probe
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--device', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--variants', type=int, default=200_000)
    parser.add_argument('--reduce', default=None)
    parser.add_argument('--shards', nargs='*', default=None, help='variant_devices instead of --device')
    args = parser.parse_args()
    os.environ.update(bench.ENVIRONMENT)
    os.environ['TORCHGWAS_SCAN_PROFILE'] = '1'
    from torchgwas.api import run_linear_gwas
    from torchgwas.pgen import PgenGenotype
    for name in ('full_scale_k512_20260926', 'full_scale_k512_missing_20260926'):
        data = args.root/name
        source = PgenGenotype(data/'input.pgen', mode='hardcall', reader_workers=8,
                              metadata_cache_dir=data/'metadata_cache')
        started = time.perf_counter()
        placement = dict(variant_devices=args.shards) if args.shards else dict(device=args.device)
        result = run_linear_gwas(source, data/'phenotype.npy', data/'covariates.npy', output_dir=args.out/name,
                                 compute_dtype='float32', chunk_size=4096, prefetch_chunks=8,
                                 reader_workers=8 * max(1, len(args.shards or [])), variant_range=(0, args.variants),
                                 **placement,
                                 **({} if args.reduce is None else dict(reduce=args.reduce)))
        seconds = time.perf_counter() - started
        profile = {key: (round(value, 3) if isinstance(value, float) else value)
                   for key, value in (getattr(source, '_last_scan_profile', None) or {}).items()
                   if isinstance(value, (int, float, str, bool))}
        timing = {key: round(value, 3) for key, value in (result.run_metadata.get('sumstats_write') or {}).items()
                  if isinstance(value, float)}
        print(json.dumps(dict(panel=name, api_seconds=round(seconds, 2), write=timing, profile=profile)), flush=True)
        shutil.rmtree(args.out/name)


if __name__ == '__main__':
    main()
