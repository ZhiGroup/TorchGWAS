"""The inflated collinear JAGWAS groups (graphunet, fourier_PE_*): old versus current statistic.

For each phenotype group, over the first V variants of the real BGEN:
  production  run_linear_gwas(reduce='jagwas') with the current code (score form,
              FP64 Gram, rounding-target cutoff): kept rank, lambda_GC, P<5e-8.
  full t      the same scan with full t output, from which, offline:
    old       a t-form quadratic form over every numerical direction of the full
              panel (what produced the inflated loci counts, without FP32 rounding luck)
    current   the score form over the production kept set (checks the offline R)
    rcond     the score form over the eigen-directions with eigenvalue > 1e-3 x max
              (the truncation under which torchGWAS and fastGWA agreed)
    jagwas_colleague_panels_20260926.py --out DIR [--variants 400000] [--device cuda:0] GROUP ...
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy import stats

from torchgwas.api import run_linear_gwas
from torchgwas.preprocess import residualize_and_standardize
from torchgwas.sumstats import open_binary_df, open_binary_sumstats
from torchgwas.sumstats_indexed import open_indexed_sumstats

ROOT = Path('/data484_4/txia2/gwas_practice')
PREP = ROOT / 'torchgwas/torchgwas_jagwas_nceq_35k_cnnonly'
PHENOS = ROOT / 'torchgwas/batch_seven_plus_torchgwas/phenos'
BGEN = ROOT / 'UKB_bgen/step4_hetqc_fusionN_rsid.bgen'
COVARIATES = PREP / 'covar_design_35k_fusionN.npy'
SAMPLES = PREP / 'sample_order_35k_fusionN_eid_eid.npy'


def summary(statistic, df):
    statistic = np.asarray(statistic, np.float64)
    statistic = statistic[np.isfinite(statistic)]
    p = stats.chi2.sf(statistic, df)
    return dict(df=int(df), variants=int(statistic.size), hits_5e8=int((p < 5e-8).sum()),
                lambda_gc=float(stats.chi2.isf(np.median(p), 1) / stats.chi2.ppf(0.5, 1)),
                max_T=float(statistic.max()), min_p=float(p.min()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    parser.add_argument('--variants', type=int, default=400_000)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('groups', nargs='+')
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    common = dict(genotype_format='bgen', bgen_decode_backend='gpu', covariates=str(COVARIATES),
                  sample_ids=str(SAMPLES), genotype_cache_dir=str(out / 'bgen_cache'), device=args.device,
                  variant_range=(0, args.variants), reader_workers=8, prefetch_chunks=16)
    for group in args.groups:
        phenotype = PHENOS / f'{group}.npy'
        row = dict(group=group)
        if not (out / group / 'jagwas' / 'sumstats' / 'manifest.json').exists():
            run_linear_gwas(str(BGEN), str(phenotype), reduce='jagwas', sumstats_fields='t',
                            output_dir=str(out / group / 'jagwas'), **common)
        manifest, parts = open_indexed_sumstats(out / group / 'jagwas' / 'sumstats')
        production = np.full(manifest['shape'][0], np.nan)
        for part in parts:
            production[np.asarray(part['variant_index'])] = np.asarray(part['chi2'], np.float64)
        rank = manifest['jagwas_rank']
        row['production'] = dict(summary(production, manifest['df']), traits=rank['traits'],
                                 estimated_rounding_error=rank['rounding_error'])
        if not (out / group / 'full' / 'sumstats' / 'manifest.json').exists():
            run_linear_gwas(str(BGEN), str(phenotype), sumstats_fields='t', sumstats_format='binary',
                            output_dir=str(out / group / 'full'), **common)
        _, t, _logp, full_manifest = open_binary_sumstats(out / group / 'full' / 'sumstats')
        variant_df = np.asarray(open_binary_df(out / group / 'full' / 'sumstats'), np.float64).reshape(-1)
        # The scan's processed panel: complete phenotypes, the same covariates and sample order.
        processed, _ = residualize_and_standardize(np.load(phenotype).astype(np.float32),
                                                   np.load(COVARIATES).astype(np.float32))
        panel = torch.as_tensor(np.ascontiguousarray(processed, dtype=np.float32))
        n = panel.shape[0]
        r64 = (panel.double().T @ panel.double() / n)
        kept = sorted(set(range(panel.shape[1])) - {item['index'] for item in rank['dropped']})
        eigen, vectors = torch.linalg.eigh(r64)
        keep = eigen > 1e-3 * eigen[-1]
        # The old statistic: t itself against the full panel. The FP32 Gram of some
        # panels is not positive definite, so use the exact R over every direction it
        # numerically has (eigenvalue > 1e-14 x max): the t-form effect, without luck.
        full = eigen > 1e-14 * eigen[-1]
        old_whitened = vectors[:, full] / torch.sqrt(eigen[full])
        new_factor = torch.linalg.solve_triangular(torch.linalg.cholesky(r64[kept][:, kept]),
                                                   torch.eye(len(kept), dtype=torch.float64), upper=False)
        whitened = vectors[:, keep] / torch.sqrt(eigen[keep])
        old, current, rcond = [], [], []
        for start in range(0, t.shape[0], 50_000):
            block = torch.as_tensor(np.asarray(t[start:start + 50_000]), dtype=torch.float64)
            df = torch.as_tensor(variant_df[start:start + 50_000])[:, None]
            z = block / torch.sqrt(1 + block * block / df)
            old.append((block @ old_whitened).square().sum(1))
            current.append((z[:, kept] @ new_factor.T).square().sum(1))
            rcond.append((z @ whitened).square().sum(1))
        old, current, rcond = (torch.cat(x).numpy() for x in (old, current, rcond))
        finite = np.isfinite(production)
        row['offline_check_max_relative'] = float(np.nanmax(np.abs(current[finite] - production[finite]) / production[finite]))
        row['old_t_form_all_directions'] = summary(old, int(full.sum()))
        row['current_offline'] = summary(current, len(kept))
        row['rcond_1e-3_score'] = summary(rcond, int(keep.sum()))
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
