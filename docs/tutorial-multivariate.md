# Multivariate GWAS Tutorial

## Overview

TorchGWAS provides a multivariate association mode for phenotype panels composed of correlated quantitative traits.
Instead of testing each trait independently, this workflow constructs a joint test statistic for each marker across the full phenotype panel.
This mode is currently treated as an experimental extension to the main linear GWAS workflow, not the primary public entry point.

The current implementation is intended for dense, continuous phenotype matrices such as imaging-derived latent traits.
As in the linear workflow, TorchGWAS aligns phenotype and covariate tables to genotype sample order and computes the covariate projection basis internally.

## When to Use This Mode

The multivariate mode is appropriate when:

- the phenotype set is composed of related quantitative traits
- a joint marker-level test is preferred to many separate univariate tests
- the primary scientific question concerns shared association across a phenotype panel

The multivariate mode is less suitable when:

- phenotypes are weakly related or belong to unrelated domains
- binary outcomes are required
- the trait count is so large that phenotype correlation inversion becomes unstable or computationally expensive

In most routine uses of TorchGWAS, start with `torchgwas linear` and treat `torchgwas multi` as a follow-on analysis once the single-trait workflow is validated.

## Minimal Reproducible Example

```bash
torchgwas multi \
  --genotype examples/toy/genotype.npy \
  --phenotype-table examples/toy/pheno.tsv \
  --covariates-table examples/toy/covar.tsv \
  --trait-columns trait_0,trait_1,trait_2 \
  --covariate-columns covar_0,covar_1,covar_2,covar_3 \
  --marker-ids examples/toy/markers.tsv \
  --output-dir multi_out
```

## Statistical Output

The multivariate workflow produces:

- `results.tsv.gz`
- `run.json`
- `qc.json`
- `phenotype_correlation.json`

Each row in `results.tsv.gz` reports:

- `marker_id`
- `n_traits`
- `chi2`
- `df`
- `p_value`
- `-log10_p`

The `phenotype_correlation.json` file records the phenotype correlation matrix used in the multivariate test.

## Recommended Workflow

The preferred sequence is:

1. validate the phenotype panel definition on the toy or linear workflow first
2. align all phenotype and covariate records by `IID`
3. run `torchgwas prep`
4. run `torchgwas multi`
5. inspect the phenotype correlation matrix and the `ridge` parameter recorded in `run.json`

Example:

```bash
torchgwas prep \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns QT0,QT1,QT2,QT3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --output-dir prep_out

torchgwas multi \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns QT0,QT1,QT2,QT3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --output-dir multi_out
```

For `BGEN` input, add `--sample-file`.

## Notes for Reporting

For publication, the multivariate analysis description should report:

- the phenotype panel definition
- the number of traits included in the joint test
- the covariates used for residualization
- the final aligned sample size
- the regularization setting used in correlation inversion
- the TorchGWAS version and exact command line
