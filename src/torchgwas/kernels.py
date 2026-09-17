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
    covariate_rank: int = -1,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    # Missing calls arrive as NaN; mask them so they take no part in the mean
    # and contribute exactly zero to every product below. See _dosage_statistics.
    observed = ~torch.isnan(genotype_chunk)
    present = observed.sum(dim=0, keepdim=True).clamp(min=1).to(genotype_chunk.dtype)
    filled = torch.where(observed, genotype_chunk,
                         torch.zeros_like(genotype_chunk))
    mean = filled.sum(dim=0, keepdim=True) / present
    centered = torch.where(observed, genotype_chunk - mean,
                           torch.zeros_like(genotype_chunk))
    gy = torch.matmul(torch.transpose(centered, 0, 1), phenotype)
    residual_ss = torch.sum(centered * centered, dim=0)
    if q_matrix is not None:
        gq = torch.matmul(torch.transpose(centered, 0, 1), q_matrix)
        residual_ss = residual_ss - torch.sum(gq * gq, dim=1)
    valid = residual_ss > eps
    # A variant spends only the samples it observed, so it keeps its own
    # residual degrees of freedom. Passing no rank keeps the scalar df.
    if covariate_rank < 0:
        variant_df = torch.full_like(residual_ss, float(df))
    else:
        variant_df = (present.squeeze(0).to(residual_ss.dtype)
                      - float(covariate_rank) - 2.0)
        valid = valid & (variant_df > 0)
    safe_df = torch.clamp(variant_df, min=1.0)
    safe_ss = torch.clamp(residual_ss, min=eps)
    beta = gy / safe_ss[:, None]
    phenotype_ss = torch.sum(phenotype * phenotype, dim=0)
    explained_ss = gy * gy / safe_ss[:, None]
    residual_y_ss = torch.clamp(phenotype_ss[None, :] - explained_ss, min=eps)
    standard_error = torch.sqrt(residual_y_ss / safe_df[:, None] / safe_ss[:, None])
    t_stat = beta / standard_error
    beta = torch.where(valid[:, None], beta, torch.zeros_like(beta))
    t_stat = torch.where(valid[:, None], t_stat, torch.zeros_like(t_stat))
    return beta, t_stat, variant_df
