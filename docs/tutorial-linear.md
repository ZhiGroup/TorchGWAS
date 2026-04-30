# Linear GWAS Tutorial

## Overview

This tutorial describes the recommended TorchGWAS workflow for single-marker linear association testing across one or more quantitative phenotypes.
The procedure is designed for publication-quality analyses in which genotype data are supplied in matrix, PLINK, or BGEN form, and phenotype/covariate information are provided as raw tables rather than precomputed projection matrices.

TorchGWAS performs the following steps internally:

1. determine genotype sample order from the input genotype resource
2. align phenotype and covariate records to that sample order by `IID`
3. remove zero-variance columns and report these decisions in `qc.json`
4. compute the covariate QR basis internally
5. residualize and standardize the phenotype matrix
6. perform marker-wise linear association testing

Accordingly, `covarQ.npy` is not a user-facing input in the standardized workflow.

## Required Inputs

TorchGWAS accepts one genotype source and one phenotype source.
Covariates are optional but strongly recommended for real analyses.

### Genotype input

One of the following formats may be used:

- `NumPy` matrix: `.npy`, shape `(n_samples, n_markers)`
- `PLINK 1` files: `.bed/.bim/.fam`
- `BGEN` genotype file with accompanying `.sample`

For large `BGEN` inputs, TorchGWAS converts the genotype resource through `plink2` into a disk-backed `memmap` cache and then reads the genotype matrix in marker batches during linear GWAS.

### Phenotype and covariate input

Phenotype and covariate data may be provided either as matrices or as tabular files.
For publication and cohort-scale use, the tabular form is preferred.

Tabular files should contain:

- an `IID` column used for sample alignment
- optional `FID`
- one or more quantitative phenotype columns
- one or more covariate columns

The final aligned analysis matrices must not contain missing values.

## Minimal Reproducible Example

The repository includes a small toy dataset that reproduces the full workflow.

```bash
torchgwas linear \
  --genotype examples/toy/genotype.npy \
  --phenotype-table examples/toy/pheno.tsv \
  --covariates-table examples/toy/covar.tsv \
  --trait-columns trait_0,trait_1,trait_2 \
  --covariate-columns covar_0,covar_1,covar_2,covar_3 \
  --marker-ids examples/toy/markers.tsv \
  --sample-ids examples/toy/samples.tsv \
  --output-dir linear_out
```

## Interpreting the Outputs

The linear workflow writes three core files:

- `results.tsv.gz`
  Marker-by-trait association statistics
- `run.json`
  Runtime metadata, input dimensions, requested columns, and analysis settings
- `qc.json`
  Quality-control summary, including dropped zero-variance columns

Each row in `results.tsv.gz` contains:

- `marker_id`
- `trait`
- `n`
- `beta`
- `se`
- `t_stat`
- `p_value`
- `-log10_p`

## Recommended Workflow for Real Cohorts

For cohort-scale analyses, the recommended sequence is:

1. prepare raw phenotype and covariate tables
2. run `torchgwas prep` to verify sample alignment and preprocessing
3. run `torchgwas linear` on the validated inputs
4. archive `run.json`, `qc.json`, and the exact command line with the results

Example:

```bash
torchgwas prep \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns trait1,trait2,trait3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --output-dir prep_out

torchgwas linear \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns trait1,trait2,trait3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --output-dir linear_out
```

The same interface applies to `BGEN` input, with the addition of `--sample-file`.

For large `BGEN` cohorts, it is recommended to keep a persistent cache directory:

```bash
torchgwas linear \
  --genotype /path/to/study.bgen \
  --genotype-format bgen \
  --sample-file /path/to/study.sample \
  --genotype-cache-dir /path/to/torchgwas_cache \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns trait1,trait2,trait3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --compute-dtype float32 \
  --output-dir linear_out
```

When a disk-backed genotype is used together with `--output-dir`, TorchGWAS streams genotype chunks from disk and writes `results.tsv.gz` incrementally.
This avoids constructing the full marker-by-trait results table in memory before writing.
For this large-scale path, `compute-dtype=auto` resolves to `float32` by default.

## Large-Scale Output Control

For very large GWAS runs, writing the complete marker-by-trait long table can become a major bottleneck.
TorchGWAS therefore supports two output-control options for the streaming linear path:

- `--topk-per-trait K`
  retain only the strongest `K` associations for each phenotype
- `--p-value-threshold X`
  retain only associations with `p <= X`

These options are intended for large cohort runs in which the primary goal is to report top signals or thresholded findings rather than to materialize the entire result matrix as a TSV long table.

Example:

```bash
torchgwas linear \
  --genotype /path/to/study.bgen \
  --genotype-format bgen \
  --sample-file /path/to/study.sample \
  --genotype-cache-dir /path/to/torchgwas_cache \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns trait1,trait2,trait3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --compute-dtype float32 \
  --topk-per-trait 100 \
  --output-dir linear_out
```

## Tracked Toy Paths

The tracked toy dataset is the preferred validation target for this tutorial.
The corresponding inputs live under:

- `examples/toy/genotype.npy`
- `examples/toy/pheno.tsv`
- `examples/toy/covar.tsv`
- `examples/toy/markers.tsv`
- `examples/toy/samples.tsv`

For a quick repository-level smoke test, use:

```bash
torchgwas demo --output-dir demo_run
```

Local cohort-specific helpers may exist under `real_case/`, but they are optional and not required for the documented public workflow.

## Unrelated-Sample Workflow

If the linear GWAS should be restricted to unrelated individuals, use:

- `scripts/run_linear_unrelated.py`

This script uses `plink2` with a `KING` kinship cutoff, filters phenotype and covariate tables to the retained `IID`s, and then runs linear TorchGWAS on the filtered cohort.

The default cutoff is:

- `0.0884`

This is the standard threshold for excluding first- and second-degree relatives.
Reference: Manichaikul A, et al. *Bioinformatics*. 2010;26(22):2867-2873.

Example:

```bash
PYTHONPATH=src python scripts/run_linear_unrelated.py \
  --genotype /path/to/study.bgen \
  --genotype-format bgen \
  --sample-file /path/to/study.sample \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --sample-id-column IID \
  --trait-columns trait1,trait2,trait3 \
  --covariate-columns SEX,AGE,PC1,PC2,PC3 \
  --output-dir unrelated_linear_run
```

The script writes:

- an unrelated genotype subset in `PLINK` format
- filtered phenotype and covariate TSV files
- a keep file listing retained samples
- the final TorchGWAS linear results under `linear/`

## Current Large-Scale Scope

The out-of-core path is currently implemented for `linear` GWAS.
The experimental `multi` workflow does not yet support disk-backed genotype streaming and should not be treated as the large-cohort path.

## Notes for Reporting

For manuscript preparation, it is advisable to report:

- genotype source and encoding format
- phenotype and covariate column definitions
- final sample size after alignment
- number of tested markers and traits
- covariates included in the linear model
- TorchGWAS version and command line
- any exclusion criteria, including unrelated-subject filtering
