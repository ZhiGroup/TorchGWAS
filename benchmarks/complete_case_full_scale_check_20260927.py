"""Full-scale dense output with missing phenotypes against exact complete-case OLS.

Runs the dense scan on `--variants` variants of the store (one GPU), then
recomputes, for every trait with missing values, FP64 least squares on the
samples with the phenotype and the call observed (intercept, covariates,
genotype) and compares the stored t and -log10 P.

    python benchmarks/complete_case_full_scale_check_20260927.py --data /data/zxie3/torchgwas_bench/full_scale_k512_missing_20260926 \\
        --device cuda:4 --out /data/zxie3/torchgwas_bench/cc_check
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from scipy import special

sys.path.insert(0, str(Path(__file__).parent))
import empirical_layout_bench_20260923 as bench  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--device', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--first', type=int, default=1_000_000)
    parser.add_argument('--variants', type=int, default=300)
    args = parser.parse_args()
    os.environ.update(bench.ENVIRONMENT)
    from torchgwas.api import run_linear_gwas
    from torchgwas.pgen import PgenGenotype
    from torchgwas.sumstats import open_binary_sumstats
    span = (args.first, args.first + args.variants)
    run_linear_gwas(genotype=str(args.data/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                    genotype_cache_dir=str(args.data/'metadata_cache'), phenotype=args.data/'phenotype.npy',
                    covariates=args.data/'covariates.npy', output_dir=args.out, device=args.device,
                    compute_dtype='float32', chunk_size=128, variant_range=span)
    _, t, logp, _ = open_binary_sumstats(args.out/'sumstats')
    t, logp = np.asarray(t), np.asarray(logp)
    source = PgenGenotype(args.data/'input.pgen', mode='hardcall', metadata_cache_dir=args.data/'metadata_cache')
    chunk = np.asarray(source.read_chunk(*span, dtype=np.float64))
    calls = chunk if chunk.shape[0] == args.variants else chunk.T
    y = np.load(args.data/'phenotype.npy').astype(np.float64)
    covariates = np.load(args.data/'covariates.npy').astype(np.float64)
    missing = np.flatnonzero(np.isnan(y).any(axis=0))
    errors_t, errors_logp = [], []
    for j in missing:
        keep_y = ~np.isnan(y[:, j])
        for v in range(args.variants):
            keep = keep_y & ~np.isnan(calls[v])
            x = np.column_stack([np.ones(keep.sum()), covariates[keep], calls[v, keep]])
            coef, *_ = np.linalg.lstsq(x, y[keep, j], rcond=None)
            resid = y[keep, j] - x @ coef
            df = keep.sum() - x.shape[1]
            exact_t = coef[-1] / np.sqrt(np.linalg.inv(x.T @ x)[-1, -1] * (resid @ resid) / df)
            exact_logp = -np.log10(2 * special.stdtr(df, -abs(exact_t)))
            errors_t.append((abs(t[v, j] - exact_t), abs(exact_t)))
            errors_logp.append((abs(logp[v, j] - exact_logp), exact_logp))
    et, el = np.array(errors_t), np.array(errors_logp)
    large_t, large_logp = et[:, 1] > 1, el[:, 1] > 1
    summary = dict(variants=args.variants, missing_traits=int(missing.size), pairs=len(et),
                   missing_calls=int(np.isnan(calls).sum()),
                   t_max_abs=float(et[:, 0].max()), t_median_abs=float(np.median(et[:, 0])),
                   t_max_rel_where_abs_t_over_1=float((et[large_t, 0] / et[large_t, 1]).max()),
                   logp_max_abs=float(el[:, 0].max()),
                   logp_max_rel_where_logp_over_1=float((el[large_logp, 0] / el[large_logp, 1]).max()),
                   pairs_logp_over_1=int(large_logp.sum()))
    print(json.dumps(summary), flush=True)
    shutil.rmtree(args.out)


if __name__ == '__main__':
    main()
