"""Missing phenotypes: which scan statistic is the complete-case OLS t, and how closely?

On `--variants` variants of the full-scale store, for every trait with
missing values, compares against exact FP64 complete-case OLS (samples with
the phenotype and the call observed; intercept, covariates, genotype):

- current: the release convention, mean-imputed t times sqrt(trait_df / df),
  at pair df variant_df * trait_df / df;
- imputed: the mean-imputed t unscaled, at the same pair df;
- imputed_exact_df: the mean-imputed t at the pair's own count minus rank - 2;
- corrected: the imputed statistic with the genotype's sums over the pair's
  own samples (sum of g and g^2 of the centred call), i.e. the exact
  denominator except for the covariates' subset refit.

The scan's statistic is reproduced in FP64 (linear._dosage_statistics on
preprocess.residualize_and_standardize). Two panels: the dataset's own
(64 traits with up to 2% missing) and a synthetic one on the same genotypes
and covariates, where missingness depends on a covariate (30% of samples,
those with the largest first covariate plus noise) and traits carry effects.

    python benchmarks/missing_phenotype_conventions_20260927.py --data /data/zxie3/torchgwas_bench/full_scale_k512_missing_20260926
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import special

from torchgwas.pgen import PgenGenotype
from torchgwas.preprocess import residualize_and_standardize


def neg_log10_p(t, df):
    return -np.log10(2 * special.stdtr(df, -np.abs(t)))


def scan_statistics(genotype, y_proc, q, rank):
    """_dosage_statistics in FP64: mean-imputed calls, full-sample covariate projection."""
    observed = ~np.isnan(genotype)
    present = observed.sum(1)
    mean = np.where(observed, genotype, 0).sum(1) / present
    g = np.where(observed, genotype - mean[:, None], 0.0)
    gy = g @ y_proc
    gq = g @ q
    rss = (g * g).sum(1) - (gq * gq).sum(1)
    ss = (y_proc * y_proc).sum(0)
    vdf = present - rank - 2
    beta = gy / rss[:, None]
    se = np.sqrt((ss[None, :] - gy * gy / rss[:, None]) / vdf[:, None] / rss[:, None])
    return beta / se, g, observed, vdf


def complete_case_t(call, y, covariates):
    keep = ~np.isnan(call) & ~np.isnan(y)
    x = np.column_stack([np.ones(keep.sum()), covariates[keep], call[keep]])
    coef, *_ = np.linalg.lstsq(x, y[keep], rcond=None)
    resid = y[keep] - x @ coef
    df = keep.sum() - x.shape[1]
    cov = np.linalg.pinv(x.T @ x) * (resid @ resid / df)
    return coef[-1] / np.sqrt(cov[-1, -1]), df


def compare(name, genotype, phenotype, covariates, max_pairs):
    y_proc, q, counts = residualize_and_standardize(phenotype, covariates, return_observed_counts=True)
    n = phenotype.shape[0]
    rank = q.shape[1]
    df = n - rank - 2
    t_imp, g, observed, vdf = scan_statistics(genotype, y_proc, q, rank)
    missing = np.flatnonzero(counts < n)
    trait_df = counts - rank - 2
    # Each missing trait residualized on [1, covariates] over its own observed rows.
    z = np.column_stack([np.full(n, 1 / np.sqrt(n)), q])
    y_complete = {}
    for j in missing:
        rows = ~np.isnan(phenotype[:, j])
        coef, *_ = np.linalg.lstsq(z[rows], phenotype[rows, j], rcond=None)
        column = np.full(n, np.nan)
        column[rows] = phenotype[rows, j] - z[rows] @ coef
        y_complete[j] = column
    rows = []
    rng = np.random.default_rng(0)
    pairs = [(v, j) for v in range(genotype.shape[0]) for j in missing]
    rng.shuffle(pairs)
    for v, j in pairs[:max_pairs]:
        call = genotype[v]
        exact_t, exact_df = complete_case_t(call, phenotype[:, j], covariates)
        pair_df = vdf[v] * trait_df[j] / df
        keep = observed[v] & ~np.isnan(phenotype[:, j])
        own_df = keep.sum() - rank - 2
        current_t = t_imp[v, j] * np.sqrt(trait_df[j] / df)
        # The denominator over the pair's own samples: centred calls restricted
        # to them and recentred there (the covariates are not refitted).
        gv = g[v, keep]
        gv = gv - gv.mean()
        y_obs = y_proc[keep, j]
        gy = gv @ y_obs
        rss = gv @ gv
        ss = y_obs @ y_obs - y_obs.sum() ** 2 / keep.sum()
        corrected_t = (gy / rss) / np.sqrt((ss - gy * gy / rss) / own_df / rss)
        # Exact complete case from full-sample quantities and the missing rows:
        # y residualized on its own observed rows (y_cc), and the call's
        # residual sum over the pair's rows by a downdate of the full-sample
        # basis Z = [1/sqrt(n), Q] over the rows it lacks.
        y_cc = y_complete[j]
        drop = ~keep
        z_drop = z[drop]
        g_full = g[v]
        zg = z.T @ g_full - z_drop.T @ g_full[drop]
        gram = np.eye(z.shape[1]) - z_drop.T @ z_drop
        rss_exact = g_full @ g_full - g_full[drop] @ g_full[drop] - zg @ np.linalg.solve(gram, zg)
        gy_exact = g_full @ np.nan_to_num(y_cc) - g_full[drop] @ np.nan_to_num(y_cc)[drop]
        ss_exact = np.nansum(y_cc[keep] ** 2)
        exact_formula_t = (gy_exact / rss_exact) / np.sqrt((ss_exact - gy_exact ** 2 / rss_exact) / own_df / rss_exact)
        # Only y residualized on its own rows; the call's denominator as in `corrected`.
        y_obs2 = y_cc[keep]
        gy2 = gv @ y_obs2
        ss2 = y_obs2 @ y_obs2
        y_only_t = (gy2 / rss) / np.sqrt((ss2 - gy2 * gy2 / rss) / own_df / rss)
        rows.append(dict(exact=neg_log10_p(exact_t, exact_df),
                         current=neg_log10_p(current_t, pair_df),
                         imputed=neg_log10_p(t_imp[v, j], pair_df),
                         imputed_exact_df=neg_log10_p(t_imp[v, j], own_df),
                         corrected=neg_log10_p(corrected_t, own_df),
                         y_only=neg_log10_p(y_only_t, own_df),
                         exact_formula=neg_log10_p(exact_formula_t, own_df),
                         exact_t=exact_t, imputed_t=t_imp[v, j], corrected_t=corrected_t))
    exact = np.array([r['exact'] for r in rows])
    summary = dict(panel=name, pairs=len(rows), missing_traits=int(missing.size),
                   max_missing_fraction=round(float(1 - counts.min() / n), 4))
    strong = exact > 2
    for key in ('current', 'imputed', 'imputed_exact_df', 'corrected', 'y_only', 'exact_formula'):
        value = np.array([r[key] for r in rows])
        relative = (value - exact) / exact
        summary[key] = dict(median_rel=float(np.median(relative)),
                            p90_abs_rel=float(np.quantile(np.abs(relative), 0.9)),
                            max_abs_rel_strong=float(np.abs(relative[strong]).max()) if strong.any() else None,
                            mean_rel_strong=float(relative[strong].mean()) if strong.any() else None)
    summary['strong_pairs'] = int(strong.sum())
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--variants', type=int, default=300)
    parser.add_argument('--pairs', type=int, default=6000)
    args = parser.parse_args()
    source = PgenGenotype(args.data/'input.pgen', mode='hardcall', metadata_cache_dir=args.data/'metadata_cache')
    first = 1_000_000
    chunk = np.asarray(source.read_chunk(first, first + args.variants, dtype=np.float64))
    genotype = chunk if chunk.shape[0] == args.variants else chunk.T  # variants x samples
    phenotype = np.load(args.data/'phenotype.npy').astype(np.float64)
    covariates = np.load(args.data/'covariates.npy').astype(np.float64)
    print(json.dumps(compare('dataset', genotype, phenotype, covariates, args.pairs)), flush=True)
    # Synthetic: covariate-dependent missingness, effects at a third of the variants.
    rng = np.random.default_rng(20260927)
    n = phenotype.shape[0]
    traits = 16
    calls = np.nan_to_num(genotype - np.nanmean(genotype, axis=1, keepdims=True))
    effects = rng.normal(scale=0.03, size=(traits, genotype.shape[0])) * (rng.random((traits, genotype.shape[0])) < 0.3)
    synthetic = rng.normal(size=(n, traits)) + calls.T @ effects.T + 0.3 * covariates[:, :1]
    score = covariates[:, 0] + rng.normal(scale=0.5, size=n)
    synthetic[score > np.quantile(score, 0.7)] = np.nan
    print(json.dumps(compare('synthetic_covariate_missing_30pct', genotype, synthetic, covariates, args.pairs)),
          flush=True)


if __name__ == '__main__':
    main()
