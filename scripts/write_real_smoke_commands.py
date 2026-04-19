from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="real_case/t1_128_smoke")
    parser.add_argument("--python-prefix", default="PYTHONPATH=src python -m torchgwas")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    script = f"""#!/usr/bin/env bash
set -euo pipefail

{args.python_prefix} prep \\
  --genotype {out / 'genotype.npy'} \\
  --phenotype-table {out / 'pheno.tsv'} \\
  --covariates-table {out / 'covar.tsv'} \\
  --sample-ids {out / 'samples.tsv'} \\
  --trait-columns $(paste -sd, {out / 'trait_columns.txt'}) \\
  --covariate-columns $(paste -sd, {out / 'covariate_columns.txt'}) \\
  --output-dir {out / 'prep'}

{args.python_prefix} linear \\
  --genotype {out / 'genotype.npy'} \\
  --phenotype-table {out / 'pheno.tsv'} \\
  --covariates-table {out / 'covar.tsv'} \\
  --sample-ids {out / 'samples.tsv'} \\
  --marker-ids {out / 'markers.tsv'} \\
  --trait-columns $(paste -sd, {out / 'trait_columns.txt'}) \\
  --covariate-columns $(paste -sd, {out / 'covariate_columns.txt'}) \\
  --output-dir {out / 'linear'}

{args.python_prefix} multi \\
  --genotype {out / 'genotype.npy'} \\
  --phenotype-table {out / 'pheno.tsv'} \\
  --covariates-table {out / 'covar.tsv'} \\
  --sample-ids {out / 'samples.tsv'} \\
  --marker-ids {out / 'markers.tsv'} \\
  --trait-columns $(paste -sd, {out / 'trait_columns.txt'}) \\
  --covariate-columns $(paste -sd, {out / 'covariate_columns.txt'}) \\
  --output-dir {out / 'multi'}
"""
    path = out / "run_smoke_test.sh"
    path.write_text(script)
    path.chmod(0o755)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
