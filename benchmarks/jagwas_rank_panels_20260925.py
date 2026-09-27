"""Rank decision of the default JAGWAS factor on real phenotype panels (no genotypes).

Residualises each panel with the scan's own preprocessing, prepares the
default (rank-revealing) reduction and reports, per panel: K, kept rank,
dropped traits, the rounding target and estimated error, the smallest pivot any order reaches
(1/max diag R^-1 = min 1 - R^2_j|others), and the eigenvalue range of R.
    jagwas_rank_panels_20260925.py --covariates C.npy [--target T] P1.npy P2.npy ...
"""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import torch

from torchgwas.jagwas_projection import JagwasRankSelection, JagwasReduction
from torchgwas.preprocess import residualize_and_standardize


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--covariates', required=True)
    parser.add_argument('--target', type=float, default=None)
    parser.add_argument('--device', default='cuda:0' if torch.cuda.is_available() else 'cpu')
    parser.add_argument('panels', nargs='+')
    args = parser.parse_args()
    covariates = np.load(args.covariates).astype(np.float64)
    for path in args.panels:
        phenotype = np.load(path, mmap_mode='r')
        if phenotype.shape[0] != covariates.shape[0] or not np.isfinite(phenotype).all():
            print(json.dumps(dict(panel=Path(path).name, skipped='rows differ from covariates or missing values')))
            continue
        processed, _ = residualize_and_standardize(np.asarray(phenotype, np.float32), covariates)
        processed = torch.as_tensor(np.ascontiguousarray(processed, dtype=np.float32))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            reduction = JagwasReduction(JagwasRankSelection(args.target)).prepare(processed, device=args.device)
        report = reduction.rank_report()
        correlation = ((processed.T @ processed) / float(processed.shape[0])).double()
        eigen = torch.linalg.eigvalsh(correlation)
        inverse = torch.linalg.pinv(correlation, hermitian=True)
        print(json.dumps(dict(panel=Path(path).name, samples=int(processed.shape[0]), traits=report['traits'],
            rank=report['rank'], method=report['method'], target=report['target'], rounding_error=report['rounding_error'],
            min_residual_variance=float(1 / torch.diagonal(inverse).max()), dropped=len(report['dropped']),
            dropped_residual_variance=[round(item['residual_variance'], 9) for item in report['dropped']][:8],
            eigenvalue_min=float(eigen[0]), eigenvalue_max=float(eigen[-1]),
            condition=float(eigen[-1] / eigen[0]) if eigen[0] > 0 else None)), flush=True)


if __name__ == '__main__':
    main()
