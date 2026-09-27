"""Full output's t and df for a partly missing trait against observed-sample OLS.

Runs the full CPU scan on a small PGEN with trait 1 missing for the first
`missing` samples and compares, per variant, the stored t (and the -log10 P
it implies at the scan's pair df) with ordinary least squares on the
observed samples only (intercept, two covariates, genotype).

    python benchmarks/missing_phenotype_t_check_20260927.py /tmp/scratch_dir
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, 'tests')
from test_pgen_native_reader import write_pgen  # noqa: E402
from torchgwas.api import run_linear_gwas  # noqa: E402
from torchgwas.sumstats import open_binary_sumstats  # noqa: E402


def ols_t(g, y, covariates):
    design = np.column_stack([np.ones(len(y)), covariates, g])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    df = len(y) - design.shape[1]
    cov = np.linalg.inv(design.T @ design) * (resid @ resid / df)
    return coef[-1] / np.sqrt(cov[-1, -1]), df


def main():
    root = Path(sys.argv[1]); root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    n, m, missing = 400, 12, 150
    calls = rng.integers(0, 3, size=(m, n), dtype=np.uint8)
    path = root/'input.pgen'
    write_pgen(path, calls)
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    y = rng.normal(size=(n, 2))
    y += 0.25 * calls[:2].T  # an effect at variants 0 and 1 in each trait
    y[:missing, 1] = np.nan
    covariates = rng.normal(size=(n, 2))
    run_linear_gwas(path, y, covariates, output_dir=root/'full', genotype_format='pgen', pgen_mode='hardcall',
                    device='cpu', compute_dtype='float64', chunk_size=4)
    _, t, logp, _ = open_binary_sumstats(root/'full/sumstats')
    t, logp = np.asarray(t), np.asarray(logp)
    observed = np.isfinite(y[:, 1])
    from torchgwas.tails import upper_tail_log10_from_t
    for trait, rows in ((0, np.ones(n, bool)), (1, observed)):
        for v in range(4):
            ref_t, ref_df = ols_t(calls[v, rows].astype(float), y[rows, trait], covariates[rows])
            print(dict(trait=trait, variant=v, scan_t=round(float(t[v, trait]), 4), ols_t=round(float(ref_t), 4),
                       scan_logp=round(float(logp[v, trait]), 4),
                       ols_logp=round(float(upper_tail_log10_from_t(np.array(ref_t), ref_df)), 4)), flush=True)


if __name__ == '__main__':
    main()
