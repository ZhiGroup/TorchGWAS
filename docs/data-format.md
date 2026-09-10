# Data Format

## Genotype

- supported inputs:
  `npy`
  `PLINK bed/bim/fam`
  `BGEN` with embedded sample IDs or `BGEN + sample`
- internal shape after loading: `(n_samples, n_markers)`
- encoding: numeric dosage or hard-call counts

PLINK BED is read directly in variant-major order. The effect allele is BIM A2.

BGEN is converted once into the established torchGWAS1.1 zstore representation:

- expected dosage of BGEN allele 2: `P(A1/A2) + 2P(A2/A2)`
- quantized byte: `round(clip(dosage, 0, 2) × 127.5)`
- level-15 Zstandard compression
- one independent frame per 2,500 retained variants (`sub=1`)
- multiallelic, phased, non-diploid, incomplete, and invalid-probability variants skipped during conversion

Skipping an unsupported site does not abort the conversion. Cohort-level MAF, HLA/MHC, chromosome, marker-ID, and allele-alphabet filters are left to the user.

The `.complete.json` sidecar records the dosage scale and exclusion counts. Values are dequantized when beta and standard error are reported; scale-invariant test statistics can operate directly on the byte codes.

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
