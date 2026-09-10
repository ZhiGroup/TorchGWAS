from __future__ import annotations

from typing import Optional, Tuple

import torch


@torch.jit.script
def linear_chunk_kernel(
    genotype_chunk: torch.Tensor,
    phenotype: torch.Tensor,
    q_matrix: Optional[torch.Tensor],
    df: int,
    eps: float = 1e-12,
) -> Tuple[torch.Tensor, torch.Tensor]:
    centered = genotype_chunk - torch.mean(genotype_chunk, dim=0, keepdim=True)
    gy = torch.matmul(torch.transpose(centered, 0, 1), phenotype)
    residual_ss = torch.sum(centered * centered, dim=0)
    if q_matrix is not None:
        gq = torch.matmul(torch.transpose(centered, 0, 1), q_matrix)
        residual_ss = residual_ss - torch.sum(gq * gq, dim=1)
    valid = residual_ss > eps
    safe_ss = torch.clamp(residual_ss, min=eps)
    beta = gy / safe_ss[:, None]
    phenotype_ss = torch.sum(phenotype * phenotype, dim=0)
    explained_ss = gy * gy / safe_ss[:, None]
    residual_y_ss = torch.clamp(phenotype_ss[None, :] - explained_ss, min=eps)
    standard_error = torch.sqrt(residual_y_ss / float(df) / safe_ss[:, None])
    t_stat = beta / standard_error
    beta = torch.where(valid[:, None], beta, torch.zeros_like(beta))
    t_stat = torch.where(valid[:, None], t_stat, torch.zeros_like(t_stat))
    return beta, t_stat


@torch.jit.script
def multivariate_chunk_kernel(
    genotype_chunk: torch.Tensor,
    phenotype: torch.Tensor,
    corr_inv: torch.Tensor,
    q_matrix: Optional[torch.Tensor],
    df: int,
    eps: float = 1e-12,
) -> torch.Tensor:
    _, t_stat = linear_chunk_kernel(genotype_chunk, phenotype, q_matrix, df, eps)
    return torch.sum(torch.matmul(t_stat, corr_inv) * t_stat, dim=1)
