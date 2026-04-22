from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import torch
from scipy import stats

from .io import DiskBackedGenotype
from .kernels import linear_chunk_kernel
from .preprocess import residualize_and_standardize
from .utils import choose_device, chunk_bounds


def linear_scan(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    device: str = "auto",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    n_traits = pheno_proc.shape[1]
    torch_device = choose_device(device)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch.float64, device=torch_device)

    beta = np.empty((n_markers, n_traits), dtype=np.float64)
    t_stat = np.empty((n_markers, n_traits), dtype=np.float64)
    p_value = np.empty((n_markers, n_traits), dtype=np.float64)

    for start, end in chunk_bounds(n_markers, chunk_size):
        geno_chunk = np.asarray(genotype[:, start:end], dtype=np.float64)
        geno_t = torch.as_tensor(geno_chunk, dtype=torch.float64, device=torch_device)
        corr_t, t_chunk_t = linear_chunk_kernel(geno_t, pheno_t)
        corr = corr_t.cpu().numpy()
        t_chunk = t_chunk_t.cpu().numpy()
        p_chunk = 2.0 * stats.t.sf(np.abs(t_chunk), df=n_samples - 2)
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    n_traits = pheno_proc.shape[1]
    chunk = chunk_size or min(n_markers, 4096) or 1
    torch_device = choose_device(device)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch.float64, device=torch_device)

    beta = np.empty((n_markers, n_traits), dtype=np.float64)
    t_stat = np.empty((n_markers, n_traits), dtype=np.float64)
    p_value = np.empty((n_markers, n_traits), dtype=np.float64)

    for start, end, geno_chunk in genotype.iter_chunks(chunk_size=chunk, dtype=np.float64, prefetch_chunks=2):
        geno_t = torch.as_tensor(geno_chunk, dtype=torch.float64, device=torch_device)
        corr_t, t_chunk_t = linear_chunk_kernel(geno_t, pheno_t)
        corr = corr_t.cpu().numpy()
        t_chunk = t_chunk_t.cpu().numpy()
        p_chunk = 2.0 * stats.t.sf(np.abs(t_chunk), df=n_samples - 2)
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
) -> tuple[Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray]], np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples = pheno_proc.shape[0]
    n_markers = genotype.shape[1]
    chunk = chunk_size or min(n_markers, 4096) or 1
    torch_device = choose_device(device)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch.float64, device=torch_device)

    def _iterator() -> Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray]]:
        for start, end, geno_chunk in genotype.iter_chunks(chunk_size=chunk, dtype=np.float64, prefetch_chunks=2):
            geno_t = torch.as_tensor(geno_chunk, dtype=torch.float64, device=torch_device)
            corr_t, t_chunk_t = linear_chunk_kernel(geno_t, pheno_t)
            corr = corr_t.cpu().numpy()
            t_chunk = t_chunk_t.cpu().numpy()
            p_chunk = 2.0 * stats.t.sf(np.abs(t_chunk), df=n_samples - 2)
            yield start, end, corr, t_chunk, p_chunk

    return _iterator(), q_matrix
