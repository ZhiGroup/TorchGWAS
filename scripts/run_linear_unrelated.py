from __future__ import annotations

"""
Filter related samples with KING/plink2, then run linear TorchGWAS on the unrelated subset.

Reference:
- Manichaikul A, et al. Robust relationship inference in genome-wide association studies.
  Bioinformatics. 2010;26(22):2867-2873.

Common KING kinship thresholds used by plink2/KING:
- > 0.354: duplicate / monozygotic twin
- > 0.177: up to 1st-degree relatives
- > 0.0884: up to 2nd-degree relatives
- > 0.0442: up to 3rd-degree relatives
"""

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from torchgwas.api import run_linear_gwas
from torchgwas.io import _find_plink2_binary, infer_genotype_format

DEFAULT_KING_CUTOFF = 0.0884


def _read_id_file(path: str | Path) -> list[str]:
    values = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line:
            values.append(line.split()[0])
    return values


def _filter_table_to_ids(
    table_path: str | Path,
    keep_ids: list[str],
    output_path: str | Path,
    sample_id_column: str = "IID",
    require_all: bool = True,
) -> pd.DataFrame:
    table_path = Path(table_path)
    table = pd.read_csv(table_path) if table_path.suffix.lower() == ".csv" else pd.read_table(table_path)
    if sample_id_column not in table.columns:
        raise ValueError(f"{table_path} is missing sample ID column {sample_id_column}")
    table[sample_id_column] = table[sample_id_column].astype(str)
    keep_index = pd.Index(keep_ids, dtype=object)
    filtered = table[table[sample_id_column].isin(set(keep_ids))].copy()
    filtered = filtered.drop_duplicates(subset=[sample_id_column]).set_index(sample_id_column)
    missing = [sample_id for sample_id in keep_ids if sample_id not in filtered.index]
    if missing and require_all:
        raise ValueError(f"{table_path} is missing {len(missing)} kept samples after relatedness filtering")
    filtered = filtered.loc[keep_index].reset_index()
    filtered = filtered.rename(columns={"index": sample_id_column})
    filtered.to_csv(output_path, sep="\t", index=False)
    return filtered


def _table_available_ids(table_path: str | Path, sample_id_column: str = "IID") -> list[str]:
    table_path = Path(table_path)
    table = pd.read_csv(table_path) if table_path.suffix.lower() == ".csv" else pd.read_table(table_path)
    if sample_id_column not in table.columns:
        raise ValueError(f"{table_path} is missing sample ID column {sample_id_column}")
    return table[sample_id_column].astype(str).drop_duplicates().tolist()


def _write_keep_file(sample_ids: list[str], path: str | Path) -> Path:
    path = Path(path)
    path.write_text("".join(f"{sample_id}\t{sample_id}\n" for sample_id in sample_ids))
    return path


def _write_keep_file_from_fam(fam_path: str | Path, sample_ids: list[str], path: str | Path) -> Path:
    fam = pd.read_table(fam_path, sep=r"\s+", header=None, names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
    fam["IID"] = fam["IID"].astype(str)
    fam["FID"] = fam["FID"].astype(str)
    keep_set = set(sample_ids)
    filtered = fam[fam["IID"].isin(keep_set)].drop_duplicates(subset=["IID"]).set_index("IID")
    missing = [sample_id for sample_id in sample_ids if sample_id not in filtered.index]
    if missing:
        raise ValueError(f"{fam_path} is missing {len(missing)} requested sample IDs for keep export")
    path = Path(path)
    path.write_text("".join(f"{filtered.loc[sample_id, 'FID']}\t{sample_id}\n" for sample_id in sample_ids))
    return path


def _run_plink2_unrelated_export(
    genotype: str | Path,
    genotype_format: str,
    output_prefix: str | Path,
    plink2_binary: str | Path | None = None,
    sample_file: str | Path | None = None,
    king_cutoff: float = 0.0884,
) -> tuple[Path, list[str]]:
    plink2 = _find_plink2_binary(plink2_binary)
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [plink2]
    if genotype_format == "bgen":
        if sample_file is None:
            raise ValueError("BGEN relatedness filtering requires --sample-file")
        cmd.extend(["--bgen", str(genotype), "ref-first", "--sample", str(sample_file)])
    elif genotype_format == "plink":
        cmd.extend(["--bfile", str(Path(genotype).with_suffix("")) if str(genotype).endswith(".bed") else str(genotype)])
    else:
        raise ValueError(f"relatedness filtering only supports PLINK/BGEN input, got {genotype_format}")

    cmd.extend(["--king-cutoff", str(king_cutoff), "--make-bed", "--out", str(output_prefix)])
    subprocess.run(cmd, check=True)

    fam_path = output_prefix.with_suffix(".fam")
    if not fam_path.exists():
        raise FileNotFoundError(f"expected filtered FAM file at {fam_path}")
    fam = pd.read_table(fam_path, sep=r"\s+", header=None)
    keep_ids = fam.iloc[:, 1].astype(str).tolist()
    return output_prefix, keep_ids


def _run_plink2_keep_export(
    genotype_prefix: str | Path,
    keep_file: str | Path,
    output_prefix: str | Path,
    plink2_binary: str | Path | None = None,
) -> Path:
    plink2 = _find_plink2_binary(plink2_binary)
    genotype_prefix = Path(genotype_prefix)
    output_prefix = Path(output_prefix)
    cmd = [
        plink2,
        "--bfile",
        str(genotype_prefix),
        "--keep",
        str(keep_file),
        "--make-bed",
        "--out",
        str(output_prefix),
    ]
    subprocess.run(cmd, check=True)
    return output_prefix


def _run_plink2_complete_case_export(
    genotype_prefix: str | Path,
    output_prefix: str | Path,
    plink2_binary: str | Path | None = None,
) -> Path:
    plink2 = _find_plink2_binary(plink2_binary)
    genotype_prefix = Path(genotype_prefix)
    output_prefix = Path(output_prefix)
    cmd = [
        plink2,
        "--bfile",
        str(genotype_prefix),
        "--geno",
        "0",
        "--mind",
        "0",
        "--make-bed",
        "--out",
        str(output_prefix),
    ]
    subprocess.run(cmd, check=True)
    return output_prefix


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Filter related samples from genotype input with plink2, then run linear TorchGWAS on the unrelated subset."
    )
    parser.add_argument("--genotype", required=True)
    parser.add_argument("--genotype-format", default="auto", choices=["auto", "plink", "bgen"])
    parser.add_argument("--sample-file", default=None)
    parser.add_argument("--phenotype-table", required=True)
    parser.add_argument("--covariates-table", required=True)
    parser.add_argument("--sample-id-column", default="IID")
    parser.add_argument("--trait-columns", default=None)
    parser.add_argument("--covariate-columns", default=None)
    parser.add_argument(
        "--king-cutoff",
        type=float,
        default=DEFAULT_KING_CUTOFF,
        help="KING kinship cutoff. 0.0884 is the standard threshold for removing 1st- and 2nd-degree relatives; 0.177 removes only 1st-degree relatives.",
    )
    parser.add_argument("--chunk-size", type=int, default=None)
    parser.add_argument("--plink2-binary", default=None)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    genotype_format = infer_genotype_format(args.genotype, genotype_format=args.genotype_format)
    filtered_prefix, keep_ids = _run_plink2_unrelated_export(
        genotype=args.genotype,
        genotype_format=genotype_format,
        sample_file=args.sample_file,
        output_prefix=output_dir / "unrelated_genotype",
        plink2_binary=args.plink2_binary,
        king_cutoff=args.king_cutoff,
    )

    pheno_ids = set(_table_available_ids(args.phenotype_table, sample_id_column=args.sample_id_column))
    covar_ids = set(_table_available_ids(args.covariates_table, sample_id_column=args.sample_id_column))
    common_keep_ids = [sample_id for sample_id in keep_ids if sample_id in pheno_ids and sample_id in covar_ids]
    if not common_keep_ids:
        raise ValueError("no unrelated samples remain after intersecting genotype, phenotype, and covariate tables")

    keep_path = _write_keep_file_from_fam(filtered_prefix.with_suffix(".fam"), common_keep_ids, output_dir / "unrelated_samples.keep")
    aligned_prefix = _run_plink2_keep_export(
        filtered_prefix,
        keep_path,
        output_dir / "unrelated_genotype_aligned",
        plink2_binary=args.plink2_binary,
    )
    complete_prefix = _run_plink2_complete_case_export(
        aligned_prefix,
        output_dir / "unrelated_genotype_complete",
        plink2_binary=args.plink2_binary,
    )
    complete_fam = pd.read_table(complete_prefix.with_suffix(".fam"), sep=r"\s+", header=None)
    complete_ids = complete_fam.iloc[:, 1].astype(str).tolist()
    final_keep_path = _write_keep_file_from_fam(
        complete_prefix.with_suffix(".fam"),
        complete_ids,
        output_dir / "unrelated_samples_complete.keep",
    )
    filtered_pheno = _filter_table_to_ids(
        args.phenotype_table,
        complete_ids,
        output_dir / "phenotype_unrelated.tsv",
        sample_id_column=args.sample_id_column,
        require_all=False,
    )
    filtered_covar = _filter_table_to_ids(
        args.covariates_table,
        complete_ids,
        output_dir / "covariates_unrelated.tsv",
        sample_id_column=args.sample_id_column,
        require_all=False,
    )

    trait_columns = None if args.trait_columns is None else [token.strip() for token in args.trait_columns.split(",") if token.strip()]
    covariate_columns = None if args.covariate_columns is None else [token.strip() for token in args.covariate_columns.split(",") if token.strip()]

    run_linear_gwas(
        genotype=complete_prefix.with_suffix(".bed"),
        genotype_format="plink",
        phenotype=None,
        covariates=None,
        phenotype_table=output_dir / "phenotype_unrelated.tsv",
        covariates_table=output_dir / "covariates_unrelated.tsv",
        trait_columns=trait_columns,
        covariate_columns=covariate_columns,
        sample_id_column=args.sample_id_column,
        chunk_size=args.chunk_size,
        output_dir=output_dir / "linear",
    )

    summary = {
        "input_genotype": str(args.genotype),
        "input_genotype_format": genotype_format,
        "king_cutoff": args.king_cutoff,
        "king_cutoff_default_recommendation": DEFAULT_KING_CUTOFF,
        "relatedness_reference": "Manichaikul A, et al. Bioinformatics (2010) 26(22):2867-2873",
        "n_unrelated_samples_before_table_intersection": len(keep_ids),
        "n_unrelated_samples": len(complete_ids),
        "keep_file": str(final_keep_path),
        "filtered_genotype_prefix": str(complete_prefix),
        "filtered_phenotype_rows": int(filtered_pheno.shape[0]),
        "filtered_covariate_rows": int(filtered_covar.shape[0]),
        "linear_output_dir": str(output_dir / "linear"),
    }
    (output_dir / "run_linear_unrelated.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
