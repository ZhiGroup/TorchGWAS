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


def linear_chunk_statistics(
    genotype_chunk: torch.Tensor,
    phenotype: torch.Tensor,
    q_matrix: Optional[torch.Tensor],
    df: int,
    covariate_rank: int = -1,
    complete_case=None,
):
    """linear_chunk_kernel, with complete-case OLS for traits with missing phenotype values.

    complete_case (complete_case.CompleteCasePlan): those traits' beta and t
    are replaced and a fourth result is the (chunk, K) pair df. The kernel is
    TorchScript, so the correction runs here, after it: the centred calls are
    formed again and only the missing traits' products are recomputed.
    """
    beta, t_stat, variant_df = linear_chunk_kernel(
        genotype_chunk, phenotype, q_matrix, df, covariate_rank=covariate_rank)
    if complete_case is None:
        return beta, t_stat, variant_df
    if not getattr(complete_case, "needs_calls", True):
        return beta, t_stat, variant_df, complete_case.correct(
            None, None, None, None, None, variant_df, beta, t_stat)
    if covariate_rank < 0:
        raise ValueError("complete-case statistics need the covariate rank")
    observed = ~torch.isnan(genotype_chunk)
    present = observed.sum(dim=0, keepdim=True).clamp(min=1).to(genotype_chunk.dtype)
    filled = torch.where(observed, genotype_chunk, torch.zeros_like(genotype_chunk))
    centered = torch.where(observed, genotype_chunk - filled.sum(dim=0, keepdim=True) / present,
                           torch.zeros_like(genotype_chunk))
    rows = centered.transpose(0, 1)  # (chunk, n)
    # Z^T g for Z = [1/sqrt(n), Q], the plan's basis.
    intercept = rows.sum(dim=1, keepdim=True) / float(rows.shape[1]) ** 0.5
    projections = intercept if q_matrix is None else torch.cat((intercept, rows @ q_matrix), dim=1)
    traits = torch.as_tensor(complete_case.traits, device=phenotype.device)
    pair_df = complete_case.correct(
        rows, torch.sum(rows * rows, dim=1), projections, None, torch.sum(phenotype * phenotype, dim=0),
        variant_df, beta, t_stat, calls_observed=lambda index: observed.transpose(0, 1).index_select(1, index),
        missing_products=rows @ phenotype.index_select(1, traits))
    return beta, t_stat, variant_df, pair_df
