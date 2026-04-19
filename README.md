# TorchGWAS

TorchGWAS is a software companion repository and reproducible project with stable inputs, command-line workflows, tutorial material, and benchmark assets.

## Availability And Implementation

TorchGWAS provides:

- linear GWAS for one or more quantitative phenotypes
- preprocessing directly from raw phenotype and covariate tables
- support for `NumPy`, `PLINK`, and `BGEN` genotype sources
- an optional unrelated-sample linear workflow built on `KING/plink2` filtering
- bundled toy examples and reproducible smoke-test workflows
- an experimental multivariate mode for future phenotype-panel extensions
- benchmark summaries and manuscript-facing figure assets

The public interface is intentionally centered on raw analysis inputs.
Users provide genotype data together with phenotype and covariate tables; TorchGWAS performs sample alignment, covariate projection, and phenotype preprocessing internally.
Accordingly, precomputed `covarQ` matrices are not part of the standardized user workflow.

## Installation

```bash
cd TorchGWAS
python -m pip install -e .
```

If editable installation is restricted by the local environment, the command-line interface can be invoked directly:

```bash
PYTHONPATH=src python -m torchgwas ...
```

## Quick Start

The repository includes a toy dataset that reproduces the complete preprocessing and association workflow:

```bash
torchgwas demo --output-dir demo_run
```

This command creates:

- `demo_run/linear/results.tsv.gz`
- `demo_run/multi/results.tsv.gz`
- `demo_run/run_summary.json`

The default user-facing workflow centers on the toy example and the linear GWAS interface.
The `demo` command also emits a multivariate result directory for exploratory comparison, but `multi` remains an extension rather than the primary entry point.

## Supported Inputs

TorchGWAS accepts the following genotype sources:

- `.npy` genotype matrices
- `PLINK 1` `.bed/.bim/.fam`
- `BGEN` genotype files with an accompanying `.sample` file

Phenotype and covariate information may be supplied as:

- aligned matrices in `.npy`, `.csv`, or `.tsv` form
- tabular files containing an `IID` column for sample alignment

For tabular inputs, TorchGWAS uses genotype sample order as the source of truth and aligns phenotype and covariate records by `IID`.
The final aligned analysis matrices are expected to be complete and quantitative.

## Recommended Workflow

The recommended analysis sequence is:

1. provide genotype data in `NumPy`, `PLINK`, or `BGEN` form
2. provide raw phenotype and covariate tables
3. run `torchgwas prep` to validate alignment and preprocessing
4. run `torchgwas linear`
5. optionally evaluate `torchgwas multi` for exploratory phenotype-panel analyses
6. archive `run.json`, `qc.json`, and the command line together with the result tables

Illustrative command pattern:

```bash
torchgwas prep \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --output-dir prep_out

torchgwas linear \
  --genotype /path/to/study.bed \
  --genotype-format plink \
  --phenotype-table /path/to/pheno.tsv \
  --covariates-table /path/to/covar.tsv \
  --output-dir linear_out
```

## Documentation

Primary user documentation is available in:

- [docs/data-format.md](docs/data-format.md)
- [docs/tutorial-linear.md](docs/tutorial-linear.md)
- [docs/tutorial-multivariate.md](docs/tutorial-multivariate.md)
- [docs/reproduce-benchmarks.md](docs/reproduce-benchmarks.md)

## Toy Workflow

The repository is organized around a tracked toy dataset under `examples/toy/`.
This path is the recommended starting point for documentation, smoke tests, and command-line examples.

Run the end-to-end demo with:

```bash
torchgwas demo --output-dir demo_run
```

Or invoke the linear workflow directly with tracked inputs:

```bash
torchgwas linear \
  --genotype examples/toy/genotype.npy \
  --phenotype-table examples/toy/pheno.tsv \
  --covariates-table examples/toy/covar.tsv \
  --marker-ids examples/toy/markers.tsv \
  --sample-ids examples/toy/samples.tsv \
  --output-dir linear_out
```

## Unrelated-Sample Linear GWAS

For studies that require unrelated samples before linear GWAS, the repository includes:

- `scripts/run_linear_unrelated.py`

This helper script:

1. runs `plink2 --king-cutoff` on the input genotype resource
2. exports an unrelated subset as `PLINK bed/bim/fam`
3. filters phenotype and covariate tables to the retained `IID`s
4. runs linear TorchGWAS on the filtered cohort

The default relatedness threshold is `0.0884`, which is the standard `KING` cutoff for removing first- and second-degree relatives.
Reference: Manichaikul A, et al. *Bioinformatics*. 2010;26(22):2867-2873.

Illustrative command pattern:

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

Local cohort-specific helpers remain under `real_case/`, but they are not part of the default public workflow.

## Verification And Benchmarks

The repository includes two complementary reproducibility layers:

- toy example outputs in `examples/`
- toy-oriented smoke-test outputs in `examples/` and `real_case/`

Benchmark summaries and manuscript-facing assets are stored in:

- `results/benchmarks/`
- `results/figures/`

These assets include PLINK concordance summaries, runtime tables, multivariate case summaries, and figure-generation scripts.
The benchmark and figure materials are primarily reported around the linear workflow; multivariate summaries are retained as exploratory assets.

## Repository Layout

- `src/torchgwas/`
  library and CLI implementation
- `docs/`
  user documentation and tutorials
- `examples/`
  toy data and expected outputs
- `scripts/`
  data-preparation, smoke-test, and figure-generation utilities
- `real_case/`
  optional local workflow assets that are not required for the toy-first public interface
- `results/`
  benchmark tables and figures
- `tests/`
  regression tests

## Scope

TorchGWAS currently focuses on a stable, documented public interface for linear association testing workflows.
The multivariate mode is retained as an experimental extension path rather than the primary advertised analysis mode.
Historical notebooks and exploratory backends from the original research environment are not treated as part of the supported public API unless explicitly migrated into this repository.

## Citation

Citation metadata is provided in [CITATION.cff](CITATION.cff).
If TorchGWAS is used in a manuscript, cite the software release together with the accompanying application note when available.
