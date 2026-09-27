"""Spurious joint signal in near-collinear panels: t versus the linear score.

t_j = r_j sqrt(df / (1 - r_j^2)) is not linear in the phenotype, so an exact
or near-exact linear dependence among traits (r_c = a r_1 + b r_2) does not
hold for t when the variant has a real effect: the collinear direction v gets
v't ~ z^3 / (2N), and dividing by its small pivot inflates T. This injects
synthetic genotype effects into real residualised panels and compares, per
pivot tolerance (the production rank selection picks the kept set S):

  T_z   production: FP32 t (the scan's formula) through the reduction, which
        forms the score z = t / sqrt(1 + t^2 / df)
  T_t   a quadratic form of t itself over the same kept factor (the old statistic)
  T_ref exact score form over S: FP64 z = sqrt(df) r, FP64 R_SS (same values)
  T_all information in the full panel: FP64 score form with R's pseudo-inverse

Effects: g = a s + sqrt(1 - a^2) e with s = Y u (residualised and
standardised), a set so the strongest single trait has |z| ~ target. u is a
random direction in trait space, a single trait among the 8 most collinear
ones (smallest 1 - R^2 given the others), or a random single trait.
    jagwas_nonlinearity_20260926.py --covariates C.npy P1.npy ...
"""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import torch

from torchgwas.jagwas_projection import JagwasRankSelection, JagwasReduction
from torchgwas.preprocess import residualize_and_standardize

Z_LEVELS = (0, 5, 10, 20, 30)
# Rounding targets for T (the selection's knob since 2026-09-26; the results
# in results/jagwas_nonlinearity_20260926.jsonl swept pivot tolerances).
TOLERANCES = (10.0, 1.0, 0.1, None, 1e-3, 1e-4, 1e-5)


def residualise(values, q):
    values = values - values.mean(0)
    if q is not None:
        values = values - q @ (q.T @ values)
    return values / values.std(0)


def scan_t(g32, p32, df):
    """The scan's FP32 t (linear._linear_stats) for residualised, centred genotypes."""
    gy = g32.T @ p32
    ss = (g32 * g32).sum(0)
    beta = gy / ss[:, None]
    explained = gy * gy / ss[:, None]
    residual = torch.clamp((p32 * p32).sum(0)[None, :] - explained, min=1e-12)
    se = torch.sqrt(residual / df / ss[:, None])
    return beta / se, gy, ss


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--covariates', required=True)
    parser.add_argument('--variants', type=int, default=200)
    parser.add_argument('panels', nargs='+')
    args = parser.parse_args()
    covariates = np.load(args.covariates).astype(np.float32)
    rng = np.random.default_rng(20260926)
    torch.set_num_threads(16)
    for path in args.panels:
        raw = np.load(path).astype(np.float32)
        processed, q = residualize_and_standardize(raw, covariates)
        p32 = torch.as_tensor(np.ascontiguousarray(processed, dtype=np.float32))
        n, k = p32.shape
        q64 = None if q is None else np.asarray(q, np.float64)
        df = float(n - (0 if q is None else q.shape[1]) - 2)
        p64 = p32.double()
        gram = p64.T @ p64 / n                       # the reduction's R (FP64 Gram of the FP32 panel)
        eigen, vectors = torch.linalg.eigh(gram)
        keep = eigen > 1e-12 * eigen[-1]
        pinv = (vectors[:, keep] / eigen[keep]) @ vectors[:, keep].T
        residual_variance = 1.0 / torch.diagonal(torch.linalg.pinv(gram, hermitian=True))
        collinear = torch.argsort(residual_variance)[:8].tolist()
        reductions = {}
        for tolerance in TOLERANCES:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                reductions[tolerance] = JagwasReduction(JagwasRankSelection(tolerance)).prepare(p32)
        for kind in ('random_direction', 'collinear_trait', 'random_trait'):
            for z in Z_LEVELS:
                m = args.variants
                if kind == 'random_direction':
                    u = rng.standard_normal((k, m))
                elif kind == 'collinear_trait':
                    u = np.zeros((k, m)); u[rng.choice(collinear, m), np.arange(m)] = 1.0
                else:
                    u = np.zeros((k, m)); u[rng.integers(0, k, m), np.arange(m)] = 1.0
                u = torch.as_tensor(u)
                ru = gram @ u
                spread = ru.abs().max(0).values / torch.sqrt((u * ru).sum(0))
                alpha = torch.clamp(z / (np.sqrt(df) * spread), max=0.95).numpy()
                s = residualise((p64 @ u).numpy(), q64)
                e = residualise(rng.standard_normal((n, m)), q64)
                g = residualise(alpha * s + np.sqrt(1 - alpha ** 2) * e, q64)
                g32 = torch.as_tensor(g, dtype=torch.float32)
                t32, gy32, ss32 = scan_t(g32, p32, df)
                variant_df = torch.full((m,), df, dtype=torch.float32)
                g64 = g32.double()
                gy64 = g64.T @ p64
                zref = np.sqrt(df) * gy64 / torch.sqrt((g64 * g64).sum(0)[:, None] * (p64 * p64).sum(0)[None, :])
                t_all = ((zref @ pinv) * zref).sum(1)
                status = torch.zeros(m, dtype=torch.uint8)
                for tolerance, reduction in reductions.items():
                    kept = torch.arange(k) if reduction._kept is None else reduction._kept
                    beta = torch.zeros_like(t32)
                    # Production (the score form from the scan's t) and, over the
                    # same factor, a quadratic form of t itself for comparison.
                    tz = reduction.reduce(beta, t32, status, variant_df, 1)[1][:, 0].double()
                    tt = (t32[:, kept].double() @ reduction._inverse_cholesky.T).square().sum(1)
                    zs = zref[:, kept]
                    factor, info = torch.linalg.cholesky_ex(gram[kept][:, kept])
                    if int(info):  # the kept block is not positive definite even in FP64
                        print(json.dumps(dict(panel=Path(path).name, kind=kind, z=z, tolerance=tolerance,
                                              rank=int(kept.numel()), error='kept R not positive definite')))
                        continue
                    tref = torch.linalg.solve_triangular(factor, zs.T, upper=False).square().sum(0)
                    dt, dz = (tt - tref).numpy(), (tz - tref).numpy()
                    report = reduction.rank_report()
                    print(json.dumps(dict(
                        panel=Path(path).name, traits=k, kind=kind, z=z, tolerance=tolerance,
                        rounding_error=report['rounding_error'], rank=report['rank'],
                        max_abs_z=float(zref.abs().max(1).values.median()),
                        t_ref_median=float(tref.median()),
                        t_excess_median=float(np.median(dt)), t_excess_p99=float(np.quantile(dt, .99)),
                        t_excess_max=float(dt.max()), t_ratio_max=float((tt / tref).max()),
                        score_error_max=float(np.abs(dz).max()),
                        # T_all - T_S is chi-square on (rank_all - r) df under the null; the
                        # excess over that is signal the dropped traits carried.
                        dropped_signal_median=float((t_all - tref).median() - (int(keep.sum()) - report['rank'])),
                        dropped_signal_max=float((t_all - tref).max() - (int(keep.sum()) - report['rank'])))),
                        flush=True)


if __name__ == '__main__':
    main()
