"""How much of T is FP32 rounding of R, as a function of the pivot tolerance.

For each real panel: the scan's R (FP32 residualisation and FP32 Gram, as a
float32 scan forms it) against an FP64 reference R (FP64 residualisation and
Gram). For each relative tolerance the kept set is chosen on the scan's R
(as the reduction does); T over that set is evaluated for null z ~ N(0, R64)
with the FP32-formed factor and with the FP64 one. Reports max |R32 - R64|,
the rank, and the relative T error (median / 99.9th percentile / max).
    jagwas_rank_error_20260925.py --covariates C.npy P1.npy ...
"""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import torch

from torchgwas.jagwas_projection import JagwasRankSelection, JagwasReduction
from torchgwas.preprocess import residualize_and_standardize

# Rounding targets for T (the selection's knob since 2026-09-26; the results
# in results/jagwas_rank_error_2026092[56].jsonl swept pivot tolerances).
TOLERANCES = (None, 10.0, 1.0, 0.1, 0.01, 1e-3, 1e-4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--covariates', required=True)
    parser.add_argument('--draws', type=int, default=20000)
    parser.add_argument('panels', nargs='+')
    args = parser.parse_args()
    covariates = np.load(args.covariates)
    for path in args.panels:
        raw = np.load(path)
        p32, _ = residualize_and_standardize(raw.astype(np.float32), covariates.astype(np.float32))
        p64, _ = residualize_and_standardize(raw.astype(np.float64), covariates.astype(np.float64))
        p32 = torch.as_tensor(np.ascontiguousarray(p32, dtype=np.float32))
        p64 = torch.as_tensor(np.ascontiguousarray(p64, dtype=np.float64))
        n, k = p32.shape
        r32 = ((p32.T @ p32) / float(n)).double()
        r64 = (p64.T @ p64) / float(n)
        r32 = torch.tril(r32) + torch.tril(r32, -1).T  # the triangle the factor reads
        # The t-statistics come from the FP32-residualised panel, so its exact
        # correlation (FP64 Gram of the same FP32 values) isolates the Gram rounding.
        exact = (p32.double().T @ p32.double()) / float(n)
        eigen, vectors = torch.linalg.eigh(r64)
        rng = np.random.default_rng(1)
        z = (torch.as_tensor(rng.standard_normal((args.draws, k))) * eigen.clamp(min=0).sqrt()) @ vectors.T
        rows = []
        for tolerance in TOLERANCES:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    reduction = JagwasReduction(JagwasRankSelection(tolerance)).prepare(p32)
                except ValueError as error:
                    rows.append(dict(tolerance=tolerance, error=str(error)[:60]))
                    continue
            report = reduction.rank_report()
            kept = torch.arange(k) if reduction._kept is None else reduction._kept
            zk = z[:, kept]
            row = dict(tolerance=tolerance, rounding_error=report['rounding_error'], rank=report['rank'])
            try:
                t32 = torch.linalg.solve_triangular(torch.linalg.cholesky(r32[kept][:, kept]), zk.T, upper=False).square().sum(0)
                t64 = torch.linalg.solve_triangular(torch.linalg.cholesky(r64[kept][:, kept]), zk.T, upper=False).square().sum(0)
                error = ((t32 - t64).abs() / t64).numpy()
                row.update(t_error_median=float(np.median(error)), t_error_p999=float(np.quantile(error, .999)),
                           t_error_max=float(error.max()), mean_t=float(t64.mean()))
                tx = torch.linalg.solve_triangular(torch.linalg.cholesky(exact[kept][:, kept]), zk.T, upper=False).square().sum(0)
                gemm = ((t32 - tx).abs() / tx).numpy()
                row.update(gemm_error_median=float(np.median(gemm)), gemm_error_max=float(gemm.max()))
                # With an exact (FP64 Gram) R, what is left is the FP32 rounding of z itself.
                zr = zk.float().double()
                tz = torch.linalg.solve_triangular(torch.linalg.cholesky(exact[kept][:, kept]), zr.T, upper=False).square().sum(0)
                rounding = ((tz - tx).abs() / tx).numpy()
                row.update(z_rounding_error_median=float(np.median(rounding)), z_rounding_error_max=float(rounding.max()))
                # The reduction's own factor (whatever Gram it forms) on the FP32-rounded z.
                tr = (reduction._inverse_cholesky.cpu() @ zr.T).square().sum(0)
                source = ((tr - tx).abs() / tx).numpy()
                row.update(reduction_error_median=float(np.median(source)), reduction_error_max=float(source.max()))
            except torch.linalg.LinAlgError:
                row.update(t_error='cholesky failed')
            rows.append(row)
        print(json.dumps(dict(panel=Path(path).name, samples=n, traits=k,
                              max_abs_r_error=float((r32 - r64).abs().max()),
                              max_abs_gemm_error=float((r32 - exact).abs().max()),
                              min_eigenvalue_r64=float(eigen[0]), rows=rows)), flush=True)


if __name__ == '__main__':
    main()
