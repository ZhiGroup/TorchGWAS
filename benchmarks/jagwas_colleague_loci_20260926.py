"""Distance-clumped loci for the stored colleague-panel scans (jagwas_colleague_panels_20260926.py).

Recomputes the old (t over every numerical direction), current (production)
and rcond=1e-3 (score form) statistics from the stored full-t scans, clumps
P<5e-8 hits greedily by P into +-500 kb loci, and counts how many current loci
have no rcond=1e-3 locus within 500 kb.
    jagwas_colleague_loci_20260926.py --out DIR GROUP ...
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy import stats

from torchgwas.preprocess import residualize_and_standardize
from torchgwas.sumstats import open_binary_df, open_binary_sumstats
from torchgwas.sumstats_indexed import open_indexed_sumstats

ROOT = Path('/data484_4/txia2/gwas_practice')
PHENOS = ROOT / 'torchgwas/batch_seven_plus_torchgwas/phenos'
COVARIATES = ROOT / 'torchgwas/torchgwas_jagwas_nceq_35k_cnnonly/covar_design_35k_fusionN.npy'


def clump(p, position, window=500_000, threshold=5e-8):
    hits = np.flatnonzero(p < threshold)
    leads = []
    for index in hits[np.argsort(p[hits])]:
        if all(abs(int(position[index]) - int(position[lead])) > window for lead in leads):
            leads.append(index)
    return np.asarray(leads, dtype=np.int64)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    parser.add_argument('groups', nargs='+')
    args = parser.parse_args()
    out = Path(args.out)
    for group in args.groups:
        manifest, parts = open_indexed_sumstats(out / group / 'jagwas' / 'sumstats')
        with np.load(out / group / 'jagwas' / 'sumstats' / 'variant_metadata.npz') as meta:
            position = np.asarray(meta['position'], np.int64)
        kept = sorted(set(range(manifest['jagwas_rank']['traits'])) - {d['index'] for d in manifest['jagwas_rank']['dropped']})
        _, t, _logp, _ = open_binary_sumstats(out / group / 'full' / 'sumstats')
        variant_df = np.asarray(open_binary_df(out / group / 'full' / 'sumstats'), np.float64).reshape(-1)
        processed, _ = residualize_and_standardize(np.load(PHENOS / f'{group}.npy').astype(np.float32),
                                                   np.load(COVARIATES).astype(np.float32))
        panel = torch.as_tensor(np.ascontiguousarray(processed, dtype=np.float32)).double()
        r = panel.T @ panel / panel.shape[0]
        eigen, vectors = torch.linalg.eigh(r)
        full, keep = eigen > 1e-14 * eigen[-1], eigen > 1e-3 * eigen[-1]
        old_w = vectors[:, full] / torch.sqrt(eigen[full])
        rcond_w = vectors[:, keep] / torch.sqrt(eigen[keep])
        factor = torch.linalg.solve_triangular(torch.linalg.cholesky(r[kept][:, kept]),
                                               torch.eye(len(kept), dtype=torch.float64), upper=False)
        values = {'old': [], 'current': [], 'rcond': []}
        for start in range(0, t.shape[0], 50_000):
            block = torch.as_tensor(np.asarray(t[start:start + 50_000]), dtype=torch.float64)
            df = torch.as_tensor(variant_df[start:start + 50_000])[:, None]
            z = block / torch.sqrt(1 + block * block / df)
            values['old'].append((block @ old_w).square().sum(1))
            values['current'].append((z[:, kept] @ factor.T).square().sum(1))
            values['rcond'].append((z @ rcond_w).square().sum(1))
        dfs = {'old': int(full.sum()), 'current': len(kept), 'rcond': int(keep.sum())}
        loci, row = {}, dict(group=group)
        for name, parts_ in values.items():
            statistic = torch.cat(parts_).numpy()
            p = np.where(np.isfinite(statistic), stats.chi2.sf(np.nan_to_num(statistic), dfs[name]), 1.0)
            loci[name] = clump(p, position)
            row[name] = dict(df=dfs[name], hits=int((p < 5e-8).sum()), loci=int(len(loci[name])))
        near = lambda a, b: np.array([np.any(np.abs(position[b] - position[i]) <= 500_000) for i in a], dtype=bool)
        row['current_loci_without_rcond_locus'] = int((~near(loci['current'], loci['rcond'])).sum()) if len(loci['current']) else 0
        row['rcond_loci_without_current_locus'] = int((~near(loci['rcond'], loci['current'])).sum()) if len(loci['rcond']) else 0
        row['span_mb'] = float((position.max() - position.min()) / 1e6)
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
