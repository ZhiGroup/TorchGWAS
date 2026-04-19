# Data Format

## Genotype

- supported inputs:
  `npy`
  `PLINK bed/bim/fam`
  `BGEN + sample`
- internal shape after loading: `(n_samples, n_markers)`
- encoding: numeric dosage or hard-call counts

## Phenotype

- `.npy`, `.csv`, `.tsv`, or `.txt`
- shape: `(n_samples, n_traits)`
- quantitative traits only in v0.1
- tabular mode should include `IID`; `FID` is allowed and ignored for alignment

## Covariates

- `.npy`, `.csv`, `.tsv`, or `.txt`
- shape: `(n_samples, n_covariates)`
- columns with zero variance are dropped and reported in `qc.json`
- user should provide raw covariates, not a precomputed `covarQ`

## Alignment

When genotype is loaded from PLINK/BGEN, TorchGWAS uses genotype sample order as the source of truth.
Phenotype and covariate tables are reordered to that sample order via `IID`.

## Outputs

Linear GWAS writes:

- `results.tsv.gz`
- `run.json`
- `qc.json`

The experimental multivariate mode also writes:

- `phenotype_correlation.json`
