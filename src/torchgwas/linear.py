from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import torch
from scipy import special

from .io import DiskBackedGenotype
from .kernels import linear_chunk_kernel
from .preprocess import residualize_and_standardize
from .utils import choose_device, chunk_bounds


def _resolve_compute_dtypes(compute_dtype: str) -> tuple[np.dtype, torch.dtype]:
    if compute_dtype == "float32":
        return np.float32, torch.float32
    if compute_dtype == "float64":
        return np.float64, torch.float64
    raise ValueError(f"unsupported compute dtype: {compute_dtype}")


def _two_sided_t_pvalue(t_stat: np.ndarray, df: int) -> np.ndarray:
    return 2.0 * special.stdtr(df, -np.abs(t_stat))


def linear_scan(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
    compute_dtype: str = "float64",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    n_traits = pheno_proc.shape[1]
    torch_device = choose_device(device)
    np_dtype, torch_dtype = _resolve_compute_dtypes(compute_dtype)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch_dtype, device=torch_device)

    beta = np.empty((n_markers, n_traits), dtype=np_dtype)
    t_stat = np.empty((n_markers, n_traits), dtype=np_dtype)
    p_value = np.empty((n_markers, n_traits), dtype=np.float64)

    for start, end in chunk_bounds(n_markers, chunk_size):
        geno_chunk = np.asarray(genotype[:, start:end], dtype=np_dtype)
        geno_t = torch.as_tensor(geno_chunk, dtype=torch_dtype, device=torch_device)
        corr_t, t_chunk_t = linear_chunk_kernel(geno_t, pheno_t)
        corr = corr_t.cpu().numpy()
        t_chunk = t_chunk_t.cpu().numpy()
        p_chunk = _two_sided_t_pvalue(t_chunk, df=n_samples - 2)
        beta[start:end] = corr
        t_stat[start:end] = t_chunk
        p_value[start:end] = p_chunk
    return beta, t_stat, p_value, q_matrix


def linear_scan_streaming(
    genotype: DiskBackedGenotype,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
    compute_dtype: str = "float32",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    n_traits = pheno_proc.shape[1]
    chunk = chunk_size or min(n_markers, 4096) or 1
    torch_device = choose_device(device)
    np_dtype, torch_dtype = _resolve_compute_dtypes(compute_dtype)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch_dtype, device=torch_device)

    beta = np.empty((n_markers, n_traits), dtype=np_dtype)
    t_stat = np.empty((n_markers, n_traits), dtype=np_dtype)
    p_value = np.empty((n_markers, n_traits), dtype=np.float64)

    for start, end, geno_chunk in genotype.iter_chunks(chunk_size=chunk, dtype=np_dtype, prefetch_chunks=2):
        geno_t = torch.as_tensor(geno_chunk, dtype=torch_dtype, device=torch_device)
        corr_t, t_chunk_t = linear_chunk_kernel(geno_t, pheno_t)
        corr = corr_t.cpu().numpy()
        t_chunk = t_chunk_t.cpu().numpy()
        p_chunk = _two_sided_t_pvalue(t_chunk, df=n_samples - 2)
        beta[start:end] = corr
        t_stat[start:end] = t_chunk
        p_value[start:end] = p_chunk
    return beta, t_stat, p_value, q_matrix


def linear_scan_streaming_chunks(
    genotype: DiskBackedGenotype,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
    compute_dtype: str = "float32",
) -> tuple[Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray]], np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    chunk = chunk_size or min(n_markers, 4096) or 1
    torch_device = choose_device(device)
    np_dtype, torch_dtype = _resolve_compute_dtypes(compute_dtype)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch_dtype, device=torch_device)

    def _iterator() -> Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray]]:
        for start, end, geno_chunk in genotype.iter_chunks(chunk_size=chunk, dtype=np_dtype, prefetch_chunks=2):
            geno_t = torch.as_tensor(geno_chunk, dtype=torch_dtype, device=torch_device)
            corr_t, t_chunk_t = linear_chunk_kernel(geno_t, pheno_t)
            corr = corr_t.cpu().numpy()
            t_chunk = t_chunk_t.cpu().numpy()
            p_chunk = _two_sided_t_pvalue(t_chunk, df=n_samples - 2)
            yield start, end, corr, t_chunk, p_chunk

    return _iterator(), q_matrix
