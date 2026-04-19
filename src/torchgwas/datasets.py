from __future__ import annotations

from pathlib import Path


def get_toy_dataset_paths(base_dir: str | Path | None = None) -> dict[str, Path]:
    root = Path(base_dir) if base_dir is not None else Path(__file__).resolve().parents[2] / "examples" / "toy"
    return {
        "genotype": root / "genotype.npy",
        "phenotype": root / "pheno.npy",
        "phenotype_table": root / "pheno.tsv",
        "covariates": root / "covar.npy",
        "covariates_table": root / "covar.tsv",
        "marker_ids": root / "markers.tsv",
        "sample_ids": root / "samples.tsv",
    }
