"""One GPU's scan pipeline at K = 512: where 2.6 us per variant goes.

TORCHGWAS_SCAN_PROFILE=1 on the first `--variants` variants of the
full-scale store, one GPU, chunk 4096. Prints the scan's own timings: the
main loop's wait for decoded chunks (fetch), its wait for results, GPU
compute / input conversion / result milliseconds, per chunk and per variant,
for min-p and dense output at a few reader counts.

    python benchmarks/gpu_pipeline_profile_20260927.py --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 \\
        --device cuda:4 --out /data/zxie3/torchgwas_bench/min_p_out/pipeline
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
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--device', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--variants', type=int, default=800_000)
    args = parser.parse_args()
    os.environ.update(bench.ENVIRONMENT)
    os.environ['TORCHGWAS_SCAN_PROFILE'] = '1'
    from torchgwas.api import run_linear_gwas
    from torchgwas.pgen import PgenGenotype
    for name, readers, extra in (('min-p', 4, dict(reduce='min-p')), ('min-p', 8, dict(reduce='min-p')),
                                 ('min-p', 16, dict(reduce='min-p')), ('dense', 8, {})):
        source = PgenGenotype(args.data/'input.pgen', mode='hardcall', reader_workers=readers,
                              metadata_cache_dir=args.data/'metadata_cache')
        out = args.out/f'{name}_r{readers}'
        started = time.perf_counter()
        run_linear_gwas(source, args.data/'phenotype.npy', args.data/'covariates.npy', output_dir=out,
                        device=args.device, compute_dtype='float32', chunk_size=4096, prefetch_chunks=readers,
                        reader_workers=readers, variant_range=(0, args.variants), **extra)
        seconds = time.perf_counter() - started
        profile = {key: value for key, value in (getattr(source, '_last_scan_profile', None) or {}).items()
                   if isinstance(value, (int, float, str, bool))}
        chunks = max(1, int(profile.get('chunks', 0) or 0))
        per_variant = {key: round(1e6 * profile[key] / args.variants, 3) for key in
                       ('fetch_seconds', 'result_wait_seconds', 'setup_seconds') if key in profile}
        per_variant.update({key.replace('milliseconds', 'us'): round(1e3 * profile[key] / args.variants, 3)
                            for key in profile if key.endswith('milliseconds')})
        print(json.dumps(dict(name=name, readers=readers, api_seconds=round(seconds, 2),
                              us_per_variant_total=round(1e6 * seconds / args.variants, 3), chunks=chunks,
                              per_variant_us=per_variant, profile=profile)), flush=True)
        source.close() if hasattr(source, 'close') else None
        shutil.rmtree(out)


if __name__ == '__main__':
    main()
