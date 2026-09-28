"""reduce='min-p': the smallest-p trait of each variant, with its exact -log10 P.

One row per variant: the winning trait, its beta, t, df and -log10 P, all
from the device. The 8.09M winners of a full-scale scan would cost 11 s of
host tail at 1.35 us a value; tails.neg_log10_p_device prices them in
milliseconds.

**For a complete phenotype panel min-P is max |t|.** The residual df is then
per *variant* (it depends on that variant's genotype missingness), so every
trait of a variant shares it and the two-sided p-value is strictly
decreasing in |t|. The winner is the largest |t| (VariantReduction's
ranking) and the tail is computed for the winners only.

**With missing phenotypes it is not.** Each trait is then tested on its own
observed samples (complete-case OLS, complete_case.py), so a pair's df is its
own sample count less rank and genotype, and p is no longer monotone in |t|
across a variant's traits. min-p then ranks by the exact tail of every cell
of the chunk at its pair df: 0.4-0.8 ns a cell on an H100 (the two tail forms,
tails._choose_form), against 0.9 ns of scan GEMM per cell at 22,250 samples
(K = 512).

This is a separate module, not a change to reduce.py: that file is the device
selector's source, and its bytes identify the recorded selector census
(device_significance_work.source_sha256).
"""

from __future__ import annotations

import torch

from .reduce import VariantReduction


class MinPReduction(VariantReduction):
    """VariantReduction('min-p') that also stages each winner's -log10 P and df.

    `stages_log10_p` tells the scan to pass `log10_p` and stage the two extra
    columns: chunks are `(start, end, beta, t, p, trait_index, -log10 P, df)`.
    """

    stages_log10_p = True

    def __init__(self):
        super().__init__("min-p")

    def host_buffers(self, chunk_size: int, width: int, pin_memory: bool = True,
                     log10_p_dtype=None):
        """The narrow ring slot, plus the kept pairs' -log10 P (at `log10_p_dtype`) and float32 df."""
        buffers = super().host_buffers(chunk_size, width, pin_memory)
        if log10_p_dtype is None:
            return buffers
        return buffers + (
            torch.empty((chunk_size, width), dtype=log10_p_dtype, pin_memory=pin_memory),
            torch.empty((chunk_size, width), dtype=torch.float32, pin_memory=pin_memory),
        )

    def reduce(self, beta, t, status, variant_df, width: int, log10_p=None):
        """`(chunk, K)` device tensors in, `(chunk, 1)` out.

        `log10_p=(pair_df, dtype)` appends the kept pairs' exact -log10 P
        (at `dtype`) and float32 df: `(beta, t, index, status, variant_df,
        logp, pair_df)`. pair_df is None for a complete panel; for missing
        phenotypes it is the scan's (chunk, K) complete-case df, t is already
        each pair's complete-case t, and the ranking is by the tail itself
        (module docstring).
        """
        if log10_p is None:
            return super().reduce(beta, t, status, variant_df, width)
        from .tails import neg_log10_p_device
        pair_df, dtype = log10_p
        if pair_df is None:
            kept = super().reduce(beta, t, status, variant_df, width)
            # The winners only: one df per variant, shared by its traits.
            pair_df = variant_df.reshape(-1, 1).expand(kept[1].shape)
            logp = neg_log10_p_device(kept[1], pair_df,
                                      out=torch.empty(kept[1].shape, dtype=dtype, device=t.device))
            return (*kept, logp, pair_df.float().contiguous())
        scores = neg_log10_p_device(t, pair_df, out=torch.empty(t.shape, dtype=torch.float64, device=t.device))
        # Same NaN discipline as VariantReduction.reduce: a NaN would sort
        # above every finite value, so it is replaced rather than compared.
        scores.nan_to_num_(nan=float("-inf"))
        scores.masked_fill_((status != 0).unsqueeze(1), float("-inf"))
        if width == 1:
            _, index = scores.max(dim=1, keepdim=True)
        else:
            _, index = scores.topk(width, dim=1, largest=True, sorted=True)
        return (beta.gather(1, index), t.gather(1, index), index.to(torch.int32), status, variant_df,
                scores.gather(1, index).to(dtype), pair_df.gather(1, index).float())

    def from_winners(self, beta, t, index, status, variant_df, dtype):
        """reduce(..., log10_p=(None, dtype)) from a complete panel's winners, already ranked.

        The Triton finish (triton_scan.finish_min_p) keeps each variant's
        largest |t| with VariantReduction's rules; this adds the winners'
        -log10 P at the variant df, so the result is reduce()'s tuple.
        """
        from .tails import neg_log10_p_device
        pair_df = variant_df.reshape(-1, 1).expand(t.shape)
        logp = neg_log10_p_device(t, pair_df, out=torch.empty(t.shape, dtype=dtype, device=t.device))
        return (beta, t, index, status, variant_df, logp, pair_df.float().contiguous())

    def merge(self, running, incoming, trait_offset: int):
        """Combine two trait blocks' winners, ranking by their -log10 P.

        Tuples are VariantReduction.merge's `(beta, t, p, index, status,
        variant_df)` plus `(logp, pair_df)` when the blocks staged them (all
        blocks or none). For a complete panel that ranking is the |t| order,
        so blocks should carry float64 -log10 P: float32 would tie pairs
        whose |t| still differ.
        """
        beta, t, p, index, status, variant_df, *pair = incoming
        if running is None:
            return (beta, t, p, index + int(trait_offset), status, variant_df, *pair)
        prior_pair = running[6:]
        if len(pair) != len(prior_pair):
            raise ValueError("merged blocks must all carry -log10 P and pair df, or none")
        if not pair:
            return super().merge(running, incoming, trait_offset)
        prior_beta, prior_t, prior_p, prior_index, prior_status, prior_df = running[:6]
        combined = [torch.cat((before, after), dim=1) for before, after in
                    ((prior_beta, beta), (prior_t, t), (prior_index, index + int(trait_offset)),
                     *zip(prior_pair, pair))]
        combined_p = None if prior_p is None or p is None else torch.cat((prior_p, p), dim=1)
        width = min(self.width, combined[1].shape[1])
        scores = combined[3].clone()
        scores.nan_to_num_(nan=float("-inf"))
        scores.masked_fill_((prior_status != 0).unsqueeze(1), float("-inf"))
        if width == 1:
            _, pick = scores.max(dim=1, keepdim=True)
        else:
            _, pick = scores.topk(width, dim=1, largest=True, sorted=True)
        kept_beta, kept_t, kept_index, kept_logp, kept_df = (value.gather(1, pick) for value in combined)
        return (kept_beta, kept_t, None if combined_p is None else combined_p.gather(1, pick),
                kept_index, prior_status, prior_df, kept_logp, kept_df)
