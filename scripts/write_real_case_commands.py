from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="real_case/t1_128")
    parser.add_argument("--genotype", default="/data4012/zxie3/all_filtered.bgen")
    parser.add_argument("--sample-file", default="/data4012/zxie3/MRI_samples_chr1.sample")
    parser.add_argument("--python-prefix", default="PYTHONPATH=src python -m torchgwas")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    prep_dir = out / "prep"
    linear_dir = out / "linear"
    multi_dir = out / "multi"
    pheno = out / "t1_pheno.tsv"
    covar = out / "t1_covar.tsv"
    trait_cols = ",".join([f"QT{i}" for i in range(128)])
    covar_cols = "SEX,GE,PC0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,AGE,25735,54"

    script = f"""#!/usr/bin/env bash
set -euo pipefail

{args.python_prefix} prep \\
  --genotype {args.genotype} \\
  --genotype-format bgen \\
  --sample-file {args.sample_file} \\
  --phenotype-table {pheno} \\
  --covariates-table {covar} \\
  --sample-id-column IID \\
  --trait-columns {trait_cols} \\
  --covariate-columns {covar_cols} \\
  --output-dir {prep_dir}

{args.python_prefix} linear \\
  --genotype {args.genotype} \\
  --genotype-format bgen \\
  --sample-file {args.sample_file} \\
  --phenotype-table {pheno} \\
  --covariates-table {covar} \\
  --sample-id-column IID \\
  --trait-columns {trait_cols} \\
  --covariate-columns {covar_cols} \\
  --output-dir {linear_dir}

{args.python_prefix} multi \\
  --genotype {args.genotype} \\
  --genotype-format bgen \\
  --sample-file {args.sample_file} \\
  --phenotype-table {pheno} \\
  --covariates-table {covar} \\
  --sample-id-column IID \\
  --trait-columns {trait_cols} \\
  --covariate-columns {covar_cols} \\
  --output-dir {multi_dir}
"""
    path = out / "run_t1_real_case.sh"
    path.write_text(script)
    path.chmod(0o755)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
