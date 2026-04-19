from __future__ import annotations

from typing import Tuple

import torch


@torch.jit.script
def _standardize_genotype_chunk(genotype_chunk: torch.Tensor, eps: float = 1e-12) -> torch.Tensor:
    centered = genotype_chunk - torch.mean(genotype_chunk, dim=0, keepdim=True)
    std = torch.std(centered, dim=0, unbiased=False, keepdim=True)
    std = torch.clamp(std, min=eps)
    return centered / std


@torch.jit.script
def linear_chunk_kernel(
    genotype_chunk: torch.Tensor,
    phenotype: torch.Tensor,
    eps: float = 1e-12,
) -> Tuple[torch.Tensor, torch.Tensor]:
    n_samples = genotype_chunk.size(0)
    geno = _standardize_genotype_chunk(genotype_chunk, eps)
    corr = torch.matmul(torch.transpose(geno, 0, 1), phenotype) / float(n_samples)
    corr = torch.clamp(corr, -0.999999999, 0.999999999)
    denom = torch.clamp(1.0 - corr * corr, min=eps)
    t_stat = corr * torch.sqrt(float(n_samples - 2) / denom)
    return corr, t_stat


@torch.jit.script
def multivariate_chunk_kernel(
    genotype_chunk: torch.Tensor,
    phenotype: torch.Tensor,
    corr_inv: torch.Tensor,
    eps: float = 1e-12,
) -> torch.Tensor:
    n_samples = genotype_chunk.size(0)
    geno = _standardize_genotype_chunk(genotype_chunk, eps)
    beta = torch.matmul(torch.transpose(geno, 0, 1), phenotype) / float(n_samples)
    beta = torch.clamp(beta, -0.999999999, 0.999999999)
    denom = torch.clamp(1.0 - beta * beta, min=eps)
    z = beta * torch.sqrt(float(n_samples - 2) / denom)
    return torch.sum(torch.matmul(z, corr_inv) * z, dim=1)
