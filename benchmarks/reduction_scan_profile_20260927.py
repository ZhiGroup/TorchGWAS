"""Per-stage scan profile of min-p against significant pairs, one GPU, K = 8,192.

TORCHGWAS_SCAN_PROFILE=1 on the first `--variants` variants of the
full-scale store; prints the scan's own timings (fetch wait, result wait,
GPU compute and result milliseconds) for each mode, plus min-p ranked by |t|
without its device tail (VariantReduction('min-p'), internal) to isolate it.

    TORCHGWAS_SCAN_PROFILE=1 python benchmarks/reduction_scan_profile_20260927.py \\
        --data /data/zxie3/torchgwas_bench/full_scale_20260925 --device cuda:4 --out /data/zxie3/torchgwas_bench/min_p_out/profile
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--device', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--variants', type=int, default=400_000)
    args = parser.parse_args()
    os.environ.update(bench.ENVIRONMENT)
    os.environ['TORCHGWAS_SCAN_PROFILE'] = '1'
    from torchgwas.api import run_linear_gwas
    from torchgwas.reduce import VariantReduction
    manifest = json.loads((args.data/'manifest.json').read_text())
    source = dict(genotype=manifest.get('genotype', str(args.data/'input.pgen')), genotype_format='pgen',
                  pgen_mode='hardcall', genotype_cache_dir=str(args.data/'metadata_cache'))
    common = dict(device=args.device, chunk_size=4096, prefetch_chunks=4, reader_workers=4,
                  variant_range=(0, args.variants), compute_dtype='float32')
    for name, extra in (('significant', dict(reduce='significant')), ('min-p', dict(reduce='min-p')),
                        ('max-abs-t', dict(_internal_reduction=VariantReduction('min-p')))):
        out = args.out/name
        started = time.perf_counter()
        result = run_linear_gwas(**source, phenotype=args.data/'phenotype.npy', covariates=args.data/'covariates.npy',
                                 output_dir=out, **common, **extra)
        seconds = time.perf_counter() - started
        meta = result.run_metadata
        profile = {key: value for key, value in (meta.get('scan_profile') or meta.get('profile') or {}).items()
                   if isinstance(value, (int, float, str, bool))}
        timing = meta.get('sumstats_write') or {}
        print(json.dumps(dict(name=name, api_seconds=round(seconds, 2),
                              executor_seconds=timing.get('setup_scan_and_write_seconds'),
                              profile=profile, keys=sorted(k for k in meta if 'prof' in k or 'scan' in k))), flush=True)
        shutil.rmtree(out)


if __name__ == '__main__':
    main()
