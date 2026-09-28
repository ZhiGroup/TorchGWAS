"""Missing phenotypes and calls: which denominator terms to correct, against exact complete-case OLS.

The scan's statistic (linear._dosage_statistics on the mean-imputed panel)
takes every sum over all n samples. For a pair with S = the trait's observed
rows O intersect the variant's observed calls A, exact OLS on [1, Q, g] over S
needs, with Z = [1/sqrt(n), Q] and H_S = Z_S'Z_S:

    rss = g_S'g_S - c'H_S^-1 c,   c = Z_S'g_S = Z'g - Z_M'g_M          (M: trait's missing rows)
    gy  = g_S'y_S - c'H_S^-1 d,   d = Z_S'y_S
    ss  = y_S'y_S - d'H_S^-1 d

g is zero at its own missing calls, so every trait-side term gathers at the
trait's (static) missing rows; the call-side terms (Y_over, the call rows in
H_S, d) gather at each variant's missing calls. Residualising each trait on
its own observed rows, with exact zeros elsewhere, makes Z_O'y = 0, so gy and
ss need no trait-side term. Compared, in relative error of -log10 P against
FP64 lstsq over S (intercept, covariates, genotype), and in absolute t error
(the null-scale view: relative error blows up where the exact -log10 P is ~0):

- release: the imputed t times sqrt(trait_df / df), at df_pair = variant_df x trait_df / df;
- expected: rss x n_trait / n and ss x n_call / n (zero cost), at df_pair;
- sums: the release statistic minus G_over (g^2 at M) and Y_over (y^2 at the missing calls), at df_pair;
- trait_exact: per-trait residualisation on O, G_over, c with the trait's H_O, and Y_over, at df_pair;
- call_gram: trait_exact with the pair's own H_S (call rows included), at df_pair;
- exact: call_gram plus d, at the pair's own count (an algebra check: should be ~0).

Panels: the full-scale store's own (64 traits up to 2% missing); the same with
5% of calls missing at random, and by a covariate; and a synthetic panel on
the same calls with 30% of samples missing by a covariate.

    python benchmarks/missing_value_scaling_20260927.py --data /data/zxie3/torchgwas_bench/full_scale_k512_missing_20260926
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import special

from torchgwas.pgen import PgenGenotype
from torchgwas.preprocess import residualize_and_standardize

METHODS = ('release', 'expected', 'sums', 'trait_exact', 'call_gram', 'exact')


def neg_log10_p(t, df):
    return -np.log10(2 * special.stdtr(df, -np.abs(t)))


def t_stat(gy, rss, ss, df):
    return gy / np.sqrt(rss * (ss - gy * gy / rss) / df)


def complete_case_residual(phenotype, z):
    """Each trait residualised on its own observed rows with exact zeros elsewhere, and each trait's Gram Z_O'Z_O."""
    y, grams = np.zeros(phenotype.shape), []
    for j in range(phenotype.shape[1]):
        o = ~np.isnan(phenotype[:, j])
        y0 = phenotype[o, j] - phenotype[o, j].mean()
        h = z[o].T @ z[o]
        y[o, j] = y0 - z[o] @ np.linalg.solve(h, z[o].T @ y0)
        grams.append(h)
    return y, grams


def compare(name, calls, phenotype, covariates, pairs, rng):
    n = phenotype.shape[0]
    y, q, counts = residualize_and_standardize(phenotype, covariates, return_observed_counts=True)
    rank = q.shape[1]
    df = n - rank - 2
    z = np.column_stack([np.full(n, n ** -0.5), q])
    y_cc, grams = complete_case_residual(phenotype, z)
    observed = ~np.isnan(calls)
    present = observed.sum(1)
    mean = np.where(observed, calls, 0).sum(1) / present
    g = np.where(observed, calls - mean[:, None], 0.0)
    gy, gz = g @ y, g @ z
    gg = (g * g).sum(1)
    rss = gg - (gz * gz).sum(1)
    ss, ss_cc = (y * y).sum(0), (y_cc * y_cc).sum(0)
    vdf = present - rank - 2
    trait_df = counts - rank - 2
    missing = np.flatnonzero(counts < n) if (counts < n).any() else np.arange(y.shape[1])
    chosen = [(v, j) for v in range(calls.shape[0]) for j in missing]
    rng.shuffle(chosen)
    rows, ts = [], []
    for v, j in chosen[:pairs]:
        m, b = np.isnan(phenotype[:, j]), ~observed[v]
        s = ~m & ~b
        x = np.column_stack([np.ones(s.sum()), covariates[s], calls[v, s]])
        coef, *_ = np.linalg.lstsq(x, phenotype[s, j], rcond=None)
        resid = phenotype[s, j] - x @ coef
        exact_df = s.sum() - x.shape[1]
        exact_t = coef[-1] / np.sqrt(np.linalg.inv(x.T @ x)[-1, -1] * (resid @ resid) / exact_df)
        pair_df = vdf[v] * trait_df[j] / df
        gm = g[v, m]
        g_over = gm @ gm
        t_imp = t_stat(gy[v, j], rss[v], ss[j], vdf[v])
        t_expected = t_stat(gy[v, j], rss[v] * counts[j] / n, ss[j] * present[v] / n, pair_df)
        t_sums = t_stat(gy[v, j], rss[v] - g_over, ss[j] - y[b, j] @ y[b, j], pair_df)
        c = gz[v] - z[m].T @ gm
        gy_cc = g[v] @ y_cc[:, j]
        ss_b = ss_cc[j] - y_cc[b, j] @ y_cc[b, j]
        t_trait = t_stat(gy_cc, gg[v] - g_over - c @ np.linalg.solve(grams[j], c), ss_b, pair_df)
        h_s = z[s].T @ z[s]
        rss_s = gg[v] - g_over - c @ np.linalg.solve(h_s, c)
        t_gram = t_stat(gy_cc, rss_s, ss_b, pair_df)
        d = -(z[b].T @ y_cc[b, j])
        t_exact = t_stat(gy_cc - c @ np.linalg.solve(h_s, d), rss_s, ss_b - d @ np.linalg.solve(h_s, d),
                         s.sum() - z.shape[1] - 1)
        t_release = t_imp * np.sqrt(trait_df[j] / df)
        rows.append((neg_log10_p(exact_t, exact_df), neg_log10_p(t_release, pair_df),
                     neg_log10_p(t_expected, pair_df), neg_log10_p(t_sums, pair_df),
                     neg_log10_p(t_trait, pair_df), neg_log10_p(t_gram, pair_df),
                     neg_log10_p(t_exact, s.sum() - z.shape[1] - 1)))
        ts.append((exact_t, t_release, t_expected, t_sums, t_trait, t_gram, t_exact))
    rows, ts = np.array(rows), np.array(ts)
    exact, strong = rows[:, 0], rows[:, 0] > 2
    summary = dict(panel=name, pairs=len(rows), strong_pairs=int(strong.sum()),
                   max_missing_trait=round(float(1 - counts.min() / n), 3),
                   missing_calls=round(float(1 - observed.mean()), 3))
    for column, key in enumerate(METHODS, start=1):
        relative = (rows[:, column] - exact) / exact
        t_error = np.abs(ts[:, column] - ts[:, 0])
        summary[key] = dict(median=round(float(np.median(relative)), 5),
                            t_abs_p99=round(float(np.quantile(t_error, 0.99)), 5),
                            t_abs_max=round(float(t_error.max()), 5),
                            strong_mean=round(float(relative[strong].mean()), 5) if strong.any() else None,
                            strong_max_abs=round(float(np.abs(relative[strong]).max()), 5) if strong.any() else None)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--variants', type=int, default=300)
    parser.add_argument('--pairs', type=int, default=4000)
    args = parser.parse_args()
    source = PgenGenotype(args.data/'input.pgen', mode='hardcall', metadata_cache_dir=args.data/'metadata_cache')
    chunk = np.asarray(source.read_chunk(1_000_000, 1_000_000 + args.variants, dtype=np.float64))
    calls = chunk if chunk.shape[0] == args.variants else chunk.T
    phenotype = np.load(args.data/'phenotype.npy').astype(np.float64)
    covariates = np.load(args.data/'covariates.npy').astype(np.float64)
    rng = np.random.default_rng(20260927)
    n = phenotype.shape[0]
    print(json.dumps(compare('dataset', calls, phenotype, covariates, args.pairs, rng)), flush=True)
    holes = calls.copy()
    holes[rng.random(holes.shape) < 0.05] = np.nan
    print(json.dumps(compare('dataset_5pct_calls_random', holes, phenotype, covariates, args.pairs, rng)), flush=True)
    score = covariates[:, 0] + rng.normal(scale=0.5, size=n)
    holes = calls.copy()
    holes[rng.random(holes.shape) < 0.1 * (score > np.median(score))[None, :]] = np.nan
    print(json.dumps(compare('dataset_5pct_calls_by_covariate', holes, phenotype, covariates, args.pairs, rng)),
          flush=True)
    centred = np.nan_to_num(calls - np.nanmean(calls, axis=1, keepdims=True))
    effects = rng.normal(scale=0.03, size=(16, calls.shape[0])) * (rng.random((16, calls.shape[0])) < 0.3)
    synthetic = rng.normal(size=(n, 16)) + centred.T @ effects.T + 0.3 * covariates[:, :1]
    synthetic[score > np.quantile(score, 0.7)] = np.nan
    print(json.dumps(compare('synthetic_30pct_by_covariate', calls, synthetic, covariates, args.pairs, rng)),
          flush=True)


if __name__ == '__main__':
    main()
