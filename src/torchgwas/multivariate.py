from __future__ import annotations

import numpy as np
import torch
from scipy import stats

from .kernels import multivariate_chunk_kernel
from .preprocess import residualize_and_standardize
from .utils import choose_device, chunk_bounds


def multivariate_scan(
    genotype: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    chunk_size: int | None = None,
    ridge: float = 1e-6,
    device: str = "auto",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    n_samples, n_traits = pheno_proc.shape
    corr = np.corrcoef(pheno_proc, rowvar=False)
    corr = np.asarray(corr, dtype=np.float64)
    corr += ridge * np.eye(n_traits)
    corr_inv = np.linalg.pinv(corr)
    torch_device = choose_device(device)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch.float64, device=torch_device)
    corr_inv_t = torch.as_tensor(corr_inv, dtype=torch.float64, device=torch_device)

    n_markers = genotype.shape[1]
    chi2 = np.empty(n_markers, dtype=np.float64)
    p_value = np.empty(n_markers, dtype=np.float64)

    for start, end in chunk_bounds(n_markers, chunk_size):
        geno_chunk = np.asarray(genotype[:, start:end], dtype=np.float64)
        geno_t = torch.as_tensor(geno_chunk, dtype=torch.float64, device=torch_device)
        chi2_chunk = multivariate_chunk_kernel(geno_t, pheno_t, corr_inv_t).cpu().numpy()
        chi2[start:end] = chi2_chunk
        p_value[start:end] = stats.chi2.sf(chi2_chunk, df=n_traits)
    return chi2, p_value, corr, corr_inv, q_matrix
