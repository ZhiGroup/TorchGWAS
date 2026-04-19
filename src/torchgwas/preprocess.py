from __future__ import annotations

import numpy as np

from .utils import check_aligned_rows, chunk_bounds, column_std_mask, ensure_2d, validate_no_missing


def prepare_inputs(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, dict]:
    genotype = ensure_2d(np.asarray(genotype, dtype=np.float64), "genotype")
    phenotype = ensure_2d(np.asarray(phenotype, dtype=np.float64), "phenotype")
    covariates = None if covariates is None else ensure_2d(np.asarray(covariates, dtype=np.float64), "covariates")

    validate_no_missing(genotype, "genotype")
    validate_no_missing(phenotype, "phenotype")
    if covariates is not None:
        validate_no_missing(covariates, "covariates")
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype), ("covariates", covariates))
    else:
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype))

    geno_mask = column_std_mask(genotype)
    pheno_mask = column_std_mask(phenotype)
    if covariates is not None:
        covar_mask = column_std_mask(covariates)
        covariates = covariates[:, covar_mask]
    else:
        covar_mask = np.array([], dtype=bool)

    qc = {
        "genotype_columns_input": int(genotype.shape[1]),
        "genotype_columns_kept": int(geno_mask.sum()),
        "phenotype_columns_input": int(phenotype.shape[1]),
        "phenotype_columns_kept": int(pheno_mask.sum()),
        "covariate_columns_input": int(0 if covariates is None else len(covar_mask)),
        "covariate_columns_kept": int(0 if covariates is None else covariates.shape[1]),
        "dropped_genotype_columns": int((~geno_mask).sum()),
        "dropped_phenotype_columns": int((~pheno_mask).sum()),
        "dropped_covariate_columns": int(0 if covariates is None else (~covar_mask).sum()),
        "n_samples": int(genotype.shape[0]),
    }
    if not geno_mask.all():
        genotype = genotype[:, geno_mask]
    if not pheno_mask.all():
        phenotype = phenotype[:, pheno_mask]
    if phenotype.shape[1] == 0:
        raise ValueError("all phenotype columns were dropped due to zero variance")
    if genotype.shape[1] == 0:
        raise ValueError("all genotype columns were dropped due to zero variance")
    return genotype, phenotype, covariates, qc


def prepare_inputs_for_prep(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None = None,
    genotype_chunk_size: int | None = None,
) -> tuple[np.ndarray, np.ndarray | None, dict]:
    genotype = ensure_2d(np.asarray(genotype), "genotype")
    phenotype = ensure_2d(np.asarray(phenotype, dtype=np.float64), "phenotype")
    covariates = None if covariates is None else ensure_2d(np.asarray(covariates, dtype=np.float64), "covariates")

    validate_no_missing(phenotype, "phenotype")
    if covariates is not None:
        validate_no_missing(covariates, "covariates")
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype), ("covariates", covariates))
    else:
        check_aligned_rows(("genotype", genotype), ("phenotype", phenotype))

    geno_mask = _chunked_genotype_std_mask(genotype, chunk_size=genotype_chunk_size)
    pheno_mask = column_std_mask(phenotype)
    if covariates is not None:
        covar_mask = column_std_mask(covariates)
        covariates = covariates[:, covar_mask]
    else:
        covar_mask = np.array([], dtype=bool)

    qc = {
        "genotype_columns_input": int(genotype.shape[1]),
        "genotype_columns_kept": int(geno_mask.sum()),
        "phenotype_columns_input": int(phenotype.shape[1]),
        "phenotype_columns_kept": int(pheno_mask.sum()),
        "covariate_columns_input": int(0 if covariates is None else len(covar_mask)),
        "covariate_columns_kept": int(0 if covariates is None else covariates.shape[1]),
        "dropped_genotype_columns": int((~geno_mask).sum()),
        "dropped_phenotype_columns": int((~pheno_mask).sum()),
        "dropped_covariate_columns": int(0 if covariates is None else (~covar_mask).sum()),
        "n_samples": int(genotype.shape[0]),
        "genotype_qc_mode": "chunked",
        "genotype_qc_chunk_size": int(genotype_chunk_size or min(genotype.shape[1], 4096) or 1),
    }

    if not pheno_mask.all():
        phenotype = phenotype[:, pheno_mask]
    if phenotype.shape[1] == 0:
        raise ValueError("all phenotype columns were dropped due to zero variance")
    return phenotype, covariates, qc


def _chunked_genotype_std_mask(genotype: np.ndarray, chunk_size: int | None = None) -> np.ndarray:
    n_markers = genotype.shape[1]
    mask = np.ones(n_markers, dtype=bool)
    for start, end in chunk_bounds(n_markers, chunk_size):
        geno_chunk = np.asarray(genotype[:, start:end], dtype=np.float64)
        if np.isnan(geno_chunk).any():
            raise ValueError("genotype contains missing values; v0.1 requires complete matrices")
        mask[start:end] = np.nanstd(geno_chunk, axis=0) > 0
    return mask


def residualize_and_standardize(phenotype: np.ndarray, covariates: np.ndarray | None) -> tuple[np.ndarray, np.ndarray | None]:
    phenotype = phenotype - phenotype.mean(axis=0, keepdims=True)
    q_matrix = None
    if covariates is not None and covariates.shape[1] > 0:
        cov_centered = covariates - covariates.mean(axis=0, keepdims=True)
        cov_std = covariates.std(axis=0, keepdims=True)
        cov_std[cov_std == 0] = 1.0
        cov_scaled = cov_centered / cov_std
        q_matrix, _ = np.linalg.qr(cov_scaled, mode="reduced")
        phenotype = phenotype - q_matrix @ (q_matrix.T @ phenotype)
    std = phenotype.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    phenotype = phenotype / std
    return phenotype, q_matrix


def standardize_genotype(genotype_chunk: np.ndarray) -> np.ndarray:
    centered = genotype_chunk - genotype_chunk.mean(axis=0, keepdims=True)
    std = centered.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    return centered / std
