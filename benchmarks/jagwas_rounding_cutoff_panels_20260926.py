"""The accuracy-target JAGWAS cutoff on real panels.

Keep the longest prefix of the greedy (pivoted Cholesky) order whose null
rounding error in T, 2 eps_z sqrt(tr R_S^-1), is at most the target. For each
panel: its FP64 numerical rank, tr(R^-1), and the kept rank for eps_z = u
sqrt(N) (the random-rounding bound) and for the largest measured GPU eps_z,
against the previous K eps32 pivot rule.
    jagwas_rounding_cutoff_panels_20260926.py --covariates C.npy P1.npy ...
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy.linalg import lapack, solve_triangular

from torchgwas.preprocess import residualize_and_standardize

U32 = 2.0 ** -24


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--covariates', required=True)
    parser.add_argument('--target', type=float, default=0.01)
    parser.add_argument('--measured', type=float, default=2.4e-6, help='largest measured eps_z at this N')
    parser.add_argument('panels', nargs='+')
    args = parser.parse_args()
    covariates = np.load(args.covariates).astype(np.float32)
    for path in args.panels:
        processed, _ = residualize_and_standardize(np.load(path).astype(np.float32), covariates)
        p = torch.as_tensor(np.ascontiguousarray(processed, dtype=np.float32)).double()
        n, k = p.shape
        r = (p.T @ p / n).numpy()
        lower, pivots, rank, _ = lapack.dpstrf(r.copy(), tol=-1.0, lower=1)
        factor = np.tril(lower)[:rank, :rank]
        inverse = solve_triangular(factor, np.eye(rank), lower=True)
        prefix = np.cumsum(np.square(inverse).sum(axis=1))  # tr of each leading block's inverse
        bound = U32 * np.sqrt(n)
        kept = {name: int(np.searchsorted(prefix, (args.target / (2 * eps)) ** 2, side='right'))
                for name, eps in (('bound', bound), ('measured', args.measured))}
        old_tol = k * np.finfo(np.float32).eps
        old = int(np.sum(np.square(np.diag(factor)) > old_tol))
        print(json.dumps(dict(panel=Path(path).name, traits=k, fp64_rank=int(rank),
                              trace_full=float(prefix[-1]), eps_bound=bound,
                              kept_bound=kept['bound'], kept_measured=kept['measured'], kept_k_eps32=old,
                              rms_error_bound_at_kept=float(2 * bound * np.sqrt(prefix[kept['bound'] - 1])),
                              rms_error_bound_full=float(2 * bound * np.sqrt(prefix[-1])))), flush=True)


if __name__ == '__main__':
    main()
