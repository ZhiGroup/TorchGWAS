from __future__ import annotations

from pathlib import Path

import numpy as np

from .io import align_table_to_samples, load_array, load_genotype, load_vector, write_table
from .linear import linear_scan
from .multivariate import multivariate_scan
from .preprocess import prepare_inputs
from .types import GWASResult, MultiGWASResult
from .utils import choose_device, elapsed, mkdir, timestamp, upper_tail_log10, write_json


def _coerce_array_or_path(value):
    if value is None:
        return None
    if isinstance(value, (str, Path)):
        return load_array(value)
    return np.asarray(value)


def _coerce_vector_or_path(value):
    if value is None:
        return None
    if isinstance(value, (str, Path)):
        return load_vector(value)
    return np.asarray(value)


def _prefer_vector(primary, fallback):
    return primary if primary is not None else fallback


def run_linear_gwas(
    genotype,
    phenotype,
    covariates=None,
    phenotype_table=None,
    covariates_table=None,
    trait_columns: list[str] | None = None,
    covariate_columns: list[str] | None = None,
    genotype_format: str = "auto",
    sample_file: str | Path | None = None,
    genotype_cache_dir: str | Path | None = None,
    bim: str | Path | None = None,
    fam: str | Path | None = None,
    plink2_binary: str | Path | None = None,
    sample_id_column: str = "IID",
    sample_ids=None,
    marker_ids=None,
    device: str = "auto",
    chunk_size: int | None = None,
    return_beta: bool = True,
    return_se: bool = True,
    return_t: bool = True,
    output_dir: str | Path | None = None,
) -> GWASResult:
    start = timestamp()
    resolved_device = choose_device(device)
    genotype_meta = {}
    if isinstance(genotype, (str, Path)):
        genotype, geno_sample_ids, geno_marker_ids, genotype_meta = load_genotype(
            genotype,
            genotype_format=genotype_format,
            bim=bim,
            fam=fam,
            sample_file=sample_file,
            genotype_cache_dir=genotype_cache_dir,
            plink2_binary=plink2_binary,
        )
    else:
        genotype = np.asarray(genotype)
        geno_sample_ids, geno_marker_ids = None, None
    marker_ids = _prefer_vector(_coerce_vector_or_path(marker_ids), geno_marker_ids)
    sample_ids = _prefer_vector(_coerce_vector_or_path(sample_ids), geno_sample_ids)
    if phenotype_table is not None:
        if sample_ids is None:
            raise ValueError("tabular phenotype input requires genotype sample IDs")
        phenotype, trait_columns = align_table_to_samples(
            phenotype_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=trait_columns, sample_id_column=sample_id_column
        )
    else:
        phenotype = _coerce_array_or_path(phenotype)
    if covariates_table is not None:
        if sample_ids is None:
            raise ValueError("tabular covariate input requires genotype sample IDs")
        covariates, covariate_columns = align_table_to_samples(
            covariates_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=covariate_columns, sample_id_column=sample_id_column
        )
    else:
        covariates = _coerce_array_or_path(covariates)
    genotype, phenotype, covariates, qc = prepare_inputs(genotype, phenotype, covariates)
    beta, t_stat, p_value, q_matrix = linear_scan(
        genotype,
        phenotype,
        covariates,
        chunk_size=chunk_size,
        device=str(resolved_device),
    )
    logp = upper_tail_log10(p_value)

    marker_names = [f"marker_{i}" for i in range(beta.shape[0])] if marker_ids is None else [str(v) for v in marker_ids[: beta.shape[0]]]
    trait_names = trait_columns or [f"trait_{i}" for i in range(beta.shape[1])]
    table = []
    for marker_index, marker_name in enumerate(marker_names):
        for trait_index, trait_name in enumerate(trait_names):
            row = {
                "marker_id": marker_name,
                "trait": trait_name,
                "n": int(genotype.shape[0]),
                "p_value": float(p_value[marker_index, trait_index]),
                "-log10_p": float(logp[marker_index, trait_index]),
            }
            if return_beta:
                row["beta"] = float(beta[marker_index, trait_index])
            if return_t:
                row["t_stat"] = float(t_stat[marker_index, trait_index])
            if return_se:
                denom = abs(t_stat[marker_index, trait_index])
                row["se"] = float(abs(beta[marker_index, trait_index]) / denom) if denom > 0 else float("nan")
            table.append(row)

    run_metadata = {
        "analysis": "linear",
        "chunk_size": int(chunk_size or min(beta.shape[0], 4096)),
        "device_requested": device,
        "device_used": str(resolved_device),
        "genotype_shape": list(genotype.shape),
        "phenotype_shape": list(phenotype.shape),
        "covariate_shape": None if covariates is None else list(covariates.shape),
        "has_sample_ids": bool(sample_ids is not None),
        "has_marker_ids": bool(marker_ids is not None),
        "sample_id_column": sample_id_column,
        "trait_columns": trait_names,
        "covariate_columns": covariate_columns,
        "q_matrix_shape": None if q_matrix is None else list(q_matrix.shape),
        "runtime_seconds": elapsed(start),
        "version": "0.1.0",
    }
    run_metadata.update(genotype_meta)
    result = GWASResult(table=table, run_metadata=run_metadata, qc_summary=qc)
    if output_dir is not None:
        out = mkdir(output_dir)
        write_table(table, out / "results.tsv.gz")
        write_json(run_metadata, out / "run.json")
        write_json(qc, out / "qc.json")
    return result


def run_multivariate_gwas(
    genotype,
    phenotype,
    covariates=None,
    phenotype_table=None,
    covariates_table=None,
    trait_columns: list[str] | None = None,
    covariate_columns: list[str] | None = None,
    genotype_format: str = "auto",
    sample_file: str | Path | None = None,
    genotype_cache_dir: str | Path | None = None,
    bim: str | Path | None = None,
    fam: str | Path | None = None,
    plink2_binary: str | Path | None = None,
    sample_id_column: str = "IID",
    sample_ids=None,
    marker_ids=None,
    device: str = "auto",
    chunk_size: int | None = None,
    output_dir: str | Path | None = None,
    ridge: float = 1e-6,
) -> MultiGWASResult:
    start = timestamp()
    resolved_device = choose_device(device)
    genotype_meta = {}
    if isinstance(genotype, (str, Path)):
        genotype, geno_sample_ids, geno_marker_ids, genotype_meta = load_genotype(
            genotype,
            genotype_format=genotype_format,
            bim=bim,
            fam=fam,
            sample_file=sample_file,
            genotype_cache_dir=genotype_cache_dir,
            plink2_binary=plink2_binary,
        )
    else:
        genotype = np.asarray(genotype)
        geno_sample_ids, geno_marker_ids = None, None
    marker_ids = _prefer_vector(_coerce_vector_or_path(marker_ids), geno_marker_ids)
    sample_ids = _prefer_vector(_coerce_vector_or_path(sample_ids), geno_sample_ids)
    if phenotype_table is not None:
        if sample_ids is None:
            raise ValueError("tabular phenotype input requires genotype sample IDs")
        phenotype, trait_columns = align_table_to_samples(
            phenotype_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=trait_columns, sample_id_column=sample_id_column
        )
    else:
        phenotype = _coerce_array_or_path(phenotype)
    if covariates_table is not None:
        if sample_ids is None:
            raise ValueError("tabular covariate input requires genotype sample IDs")
        covariates, covariate_columns = align_table_to_samples(
            covariates_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=covariate_columns, sample_id_column=sample_id_column
        )
    else:
        covariates = _coerce_array_or_path(covariates)
    genotype, phenotype, covariates, qc = prepare_inputs(genotype, phenotype, covariates)
    chi2, p_value, corr, corr_inv, q_matrix = multivariate_scan(
        genotype,
        phenotype,
        covariates,
        chunk_size=chunk_size,
        ridge=ridge,
        device=str(resolved_device),
    )
    logp = upper_tail_log10(p_value)
    marker_names = [f"marker_{i}" for i in range(chi2.shape[0])] if marker_ids is None else [str(v) for v in marker_ids[: chi2.shape[0]]]
    table = [
        {
            "marker_id": marker_names[idx],
            "n_traits": int(phenotype.shape[1]),
            "chi2": float(chi2[idx]),
            "df": int(phenotype.shape[1]),
            "p_value": float(p_value[idx]),
            "-log10_p": float(logp[idx]),
        }
        for idx in range(chi2.shape[0])
    ]
    run_metadata = {
        "analysis": "multivariate",
        "chunk_size": int(chunk_size or min(chi2.shape[0], 4096)),
        "device_requested": device,
        "device_used": str(resolved_device),
        "genotype_shape": list(genotype.shape),
        "phenotype_shape": list(phenotype.shape),
        "covariate_shape": None if covariates is None else list(covariates.shape),
        "has_sample_ids": bool(sample_ids is not None),
        "has_marker_ids": bool(marker_ids is not None),
        "sample_id_column": sample_id_column,
        "trait_columns": trait_columns or [f"trait_{i}" for i in range(phenotype.shape[1])],
        "covariate_columns": covariate_columns,
        "q_matrix_shape": None if q_matrix is None else list(q_matrix.shape),
        "corr_inverse_shape": list(corr_inv.shape),
        "ridge": ridge,
        "runtime_seconds": elapsed(start),
        "version": "0.1.0",
    }
    run_metadata.update(genotype_meta)
    result = MultiGWASResult(
        table=table,
        run_metadata=run_metadata,
        qc_summary=qc,
        phenotype_correlation=corr.tolist(),
    )
    if output_dir is not None:
        out = mkdir(output_dir)
        write_table(table, out / "results.tsv.gz")
        write_json(run_metadata, out / "run.json")
        write_json(qc, out / "qc.json")
        write_json({"phenotype_correlation": corr.tolist()}, out / "phenotype_correlation.json")
    return result
