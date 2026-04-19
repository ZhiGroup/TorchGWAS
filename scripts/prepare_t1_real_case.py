from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_COVARIATE_COLUMNS = [
    "SEX",
    "GE",
    "PC0",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "AGE",
    "25735",
    "54",
]


def read_space_table(path: str | Path) -> pd.DataFrame:
    return pd.read_table(path, sep=" ")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--covar-discovery", default="/data4012/zxie3/MRI_AE_small_training_set_GWAS/T1_covar_discovery")
    parser.add_argument("--covar-replication", default="/data4012/zxie3/MRI_AE_small_training_set_GWAS/T1_covar_replication")
    parser.add_argument("--pheno-discovery", default="/data4012/zxie3/MRI_AE_small_training_set_GWAS/T1_pheno_discovery")
    parser.add_argument("--pheno-replication", default="/data4012/zxie3/MRI_AE_small_training_set_GWAS/T1_pheno_replication")
    parser.add_argument("--unrelated-ids", default="/data4012/zxie3/fast_GWAS/unrelated_UKB_MRI_subjects")
    parser.add_argument("--output-dir", default="real_case/t1_128")
    parser.add_argument("--keep-related", action="store_true")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    covar = pd.concat(
        [read_space_table(args.covar_discovery), read_space_table(args.covar_replication)],
        ignore_index=True,
    ).drop_duplicates(subset=["FID", "IID"])
    pheno = pd.concat(
        [read_space_table(args.pheno_discovery), read_space_table(args.pheno_replication)],
        ignore_index=True,
    ).drop_duplicates(subset=["FID", "IID"])

    covar["IID"] = covar["IID"].astype(str)
    covar["FID"] = covar["FID"].astype(str)
    pheno["IID"] = pheno["IID"].astype(str)
    pheno["FID"] = pheno["FID"].astype(str)

    if not args.keep_related:
        unrelated = {line.strip() for line in Path(args.unrelated_ids).read_text().splitlines() if line.strip()}
        covar = covar[covar["IID"].isin(unrelated)].copy()
        pheno = pheno[pheno["IID"].isin(unrelated)].copy()

    trait_columns = [col for col in pheno.columns if col.startswith("QT")]
    covariate_columns = [col for col in DEFAULT_COVARIATE_COLUMNS if col in covar.columns]

    covar = covar[["FID", "IID", *covariate_columns]].dropna()
    pheno = pheno[["FID", "IID", *trait_columns]].dropna()

    covar.to_csv(out / "t1_covar.tsv", sep="\t", index=False)
    pheno.to_csv(out / "t1_pheno.tsv", sep="\t", index=False)
    (out / "trait_columns.txt").write_text("\n".join(trait_columns) + "\n")
    (out / "covariate_columns.txt").write_text("\n".join(covariate_columns) + "\n")

    summary = {
        "n_covar_rows": int(covar.shape[0]),
        "n_pheno_rows": int(pheno.shape[0]),
        "n_traits": int(len(trait_columns)),
        "n_covariates": int(len(covariate_columns)),
        "filtered_to_unrelated": not args.keep_related,
        "covariate_columns": covariate_columns,
        "trait_columns_first10": trait_columns[:10],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

