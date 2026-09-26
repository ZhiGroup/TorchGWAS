"""Device-side per-variant reduction across traits.

A scan produces one statistic per (variant, trait) pair, and past a few thousand
traits that product *is* the run: 2M variants by 10^5 traits is 800 GB of
float32 t-statistics, which will not fit in the pinned result ring, will not fit
on the output disk, and is not what anyone asked for.  What a high-trait scan
wants is the best trait per variant, or the best few.

That reduction is one pass over a tensor the statistics kernel has already
produced on the device, so it belongs there: the host then receives a
`chunk x k` result instead of a `chunk x K` one, and at K = 10^5 with k = 1 that
is a 100,000-fold narrowing of everything downstream -- the pinned ring, the
D2H copy, the Student-t evaluation, and the file.

**min-P and max-t^2 are the same reduction here, and that is worth stating
rather than implementing twice.** The residual degrees of freedom are per
*variant* (they depend on that variant's missingness), not per trait, so within
one variant every trait shares a df and the two-sided p-value is strictly
decreasing in |t|.  Ranking by |t| and ranking by p therefore give the identical
ordering and the identical winner.  `min-p` is kept as a spelling because it is
what the literature calls the quantity, but it selects by |t| and the p-value is
computed afterwards on the k survivors only.
"""

from __future__ import annotations

import numpy as np
import torch


_SINGLE = ("max-abs-t", "max-t2", "min-p")


class VariantReduction:
    """Select the k strongest traits per variant, on the device.

    `mode` is one of `max-abs-t`, `max-t2`, `min-p` (all k = 1 and all the same
    ranking -- see the module docstring) or `top-k`, which needs `top_k`.
    """

    def __init__(self, mode: str = "max-abs-t", top_k: int | None = None):
        key = str(mode).replace("_", "-").lower()
        if key in _SINGLE:
            if top_k not in (None, 1):
                raise ValueError(
                    f"reduction mode {mode!r} selects one trait per variant; "
                    "use mode 'top-k' to keep more than one")
            self.width = 1
        elif key == "top-k":
            if top_k is None or int(top_k) < 1:
                raise ValueError("mode 'top-k' requires a positive top_k")
            self.width = int(top_k)
        else:
            raise ValueError(
                f"unknown reduction mode {mode!r}: "
                f"expected one of {', '.join((*_SINGLE, 'top-k'))}")
        self.mode = key

    def resolved_width(self, n_traits: int) -> int:
        """k, capped at the number of traits actually present."""
        if n_traits < 1:
            raise ValueError("a reduction needs at least one trait")
        return min(self.width, int(n_traits))

    def host_buffers(self, chunk_size: int, width: int, pin_memory: bool = True):
        """Pinned staging for one ring slot: the narrow results, not the wide ones."""
        return (
            torch.empty((chunk_size, width), dtype=torch.float32, pin_memory=pin_memory),
            torch.empty((chunk_size, width), dtype=torch.float32, pin_memory=pin_memory),
            torch.empty((chunk_size, width), dtype=torch.int32, pin_memory=pin_memory),
            torch.empty(chunk_size, dtype=torch.uint8, pin_memory=pin_memory),
            torch.empty(chunk_size, pin_memory=pin_memory),
        )

    def merge(self, running, incoming, trait_offset: int):
        """Combine two reduced results for the same variants, keeping the best k.

        This is what makes trait tiling possible. A scan that cannot hold all K
        traits on the device processes them in blocks, and each block returns
        its own best k with trait indices numbered *within the block*; merging
        rebases those onto the full trait axis and keeps the k strongest across
        everything seen so far.

        Both arguments are `(beta, t, p, index, status, variant_df)` tuples, or
        `running` may be None for the first block. `p` may be None when the
        caller defers p-values; when present it is carried through the same
        selection rather than recomputed, which is what lets a blocked scan keep
        **per-variant** degrees of freedom. Recomputing p downstream from a
        single scalar df would be wrong wherever missingness varies by variant,
        and would be wrong silently.

        The merge is exact, not approximate: the top k of a union is contained
        in the union of the two top-k sets, so nothing that should have been
        kept can have been discarded by an earlier block.

        `status` and `variant_df` are per variant and identical across blocks --
        the same genotypes decide both -- so the running copy is kept.
        """
        beta, t, p, index, status, variant_df = incoming
        index = index + int(trait_offset)
        if running is None:
            return (beta, t, p, index, status, variant_df)
        prior_beta, prior_t, prior_p, prior_index, prior_status, prior_df = running
        combined_beta = torch.cat((prior_beta, beta), dim=1)
        combined_t = torch.cat((prior_t, t), dim=1)
        combined_index = torch.cat((prior_index, index), dim=1)
        combined_p = (None if (prior_p is None or p is None)
                      else torch.cat((prior_p, p), dim=1))
        width = min(self.width, combined_t.shape[1])
        # Same NaN discipline as `reduce`: a NaN would sort above every finite
        # value, so it is replaced rather than compared against.
        scores = combined_t.abs()
        scores.nan_to_num_(nan=float("-inf"))
        scores.masked_fill_((prior_status != 0).unsqueeze(1), float("-inf"))
        if width == 1:
            _, pick = scores.max(dim=1, keepdim=True)
        else:
            _, pick = scores.topk(width, dim=1, largest=True, sorted=True)
        return (combined_beta.gather(1, pick), combined_t.gather(1, pick),
                None if combined_p is None else combined_p.gather(1, pick),
                combined_index.gather(1, pick), prior_status, prior_df)

    def reduce(self, beta, t, status, variant_df, width: int):
        """`(chunk, K)` device tensors in, `(chunk, k)` device tensors out.

        `status` and `variant_df` pass through untouched: they are per variant
        already, and the caller applies the same NaN-on-invalid rule it applies
        to an unreduced chunk, so a reduced scan and a full scan disagree
        nowhere.
        """
        # torch orders NaN *above* every finite value in `topk` and `max`, so an
        # invalid variant would otherwise win its own row and be reported as the
        # top hit. Guard the row by status and scrub any remaining NaN, rather
        # than trusting a comparison against NaN to do the right thing.
        #
        # The mask is applied unconditionally, never behind an `invalid.any()`
        # test: reading that predicate on the host would synchronise the compute
        # stream on every chunk, which is exactly the overlap this scan is built
        # around. An unconditional elementwise pass is far cheaper than a stall.
        # `t.abs()` already allocates, so the rest is in place on that copy.
        scores = t.abs()
        scores.nan_to_num_(nan=float("-inf"))
        scores.masked_fill_((status != 0).unsqueeze(1), float("-inf"))
        if width == 1:
            _, index = scores.max(dim=1, keepdim=True)
        else:
            _, index = scores.topk(width, dim=1, largest=True, sorted=True)
        return (beta.gather(1, index), t.gather(1, index),
                index.to(torch.int32), status, variant_df)


GENOME_WIDE_ALPHA = 5e-8
"""The conventional genome-wide significance threshold for a single trait.

It is not `0.05 / M`. It is the long-standing 5e-8, which already carries a
Bonferroni correction for the *effective* number of independent common variants
in a European-ancestry genome -- roughly a million, not the several million
actually tested. Dividing 0.05 by the number of variants in the file would
therefore correct twice and would move with the file rather than with the
genome.
"""


def bonferroni_threshold(n_traits: int, alpha: float = GENOME_WIDE_ALPHA) -> float:
    """The per-test threshold for a `n_traits`-trait scan: `5e-8 / K`.

    The genome axis is already paid for by `alpha`; what a multi-trait scan
    adds is the trait axis, so the correction is over K and only K.
    """
    if int(n_traits) < 1:
        raise ValueError("a threshold needs at least one trait")
    if not (0.0 < float(alpha) <= 1.0):
        raise ValueError("alpha must be in (0, 1]")
    return float(alpha) / int(n_traits)


class SignificantPairs:
    """Keep only (variant, trait) pairs passing a significance threshold.

    This is the other thing a high-trait scan might want, and it is not a
    special case of `VariantReduction`: that one keeps a *fixed* k per variant,
    so its output is `chunk x k` and fits a preallocated ring. This one keeps
    however many pairs happen to pass, which is data-dependent and usually
    almost none -- at K = 128 the threshold is 5e-8/128 = 3.9e-10, so the
    expected count under the null across 8.93M variants is **0.0035 pairs**.

    That asymmetry is the whole point. The full result is `M x K` -- 1.14
    billion rows at this shape, 9.2 GB even in binary -- and the answer is a
    handful of rows. Selecting on the device means the discarded 1.14 billion
    never cross PCIe.

    **The threshold is compared on |t|, not on p, and that is exact rather than
    an approximation.** Residual degrees of freedom are per *variant*, so
    within one variant the two-sided p-value is strictly decreasing in |t|;
    converting the p threshold to a critical |t| once per variant and comparing
    against that gives precisely the same set as computing every p and
    comparing. It costs one `stdtrit` call per chunk of variants instead of one
    per (variant, trait) cell.

    **NaN needs no special handling here, unlike in `VariantReduction`.** There
    the danger was that `topk` sorts NaN *above* every finite value, so an
    invalid variant would win its own row. A comparison is the opposite: every
    comparison with NaN is False, so an invalid statistic simply fails the
    threshold. The status mask is still applied, because a variant can be
    invalid while its statistics are finite.
    """

    def __init__(self, threshold: float | None = None,
                 alpha: float = GENOME_WIDE_ALPHA):
        if threshold is not None and not (0.0 < float(threshold) <= 1.0):
            raise ValueError("threshold must be in (0, 1]")
        self.threshold = None if threshold is None else float(threshold)
        self.alpha = float(alpha)

    def resolved_threshold(self, n_traits: int) -> float:
        """The explicit threshold if given, else `alpha / K`."""
        if self.threshold is not None:
            return self.threshold
        return bonferroni_threshold(n_traits, self.alpha)

    def critical_abs_t(self, variant_df, n_traits: int):
        """Per-variant |t| at which the two-sided p-value equals the threshold.

        `variant_df` may be a numpy array (one df per variant) or a scalar.
        Returned as a numpy float64 array shaped to broadcast down a column.
        """
        from scipy import special

        threshold = self.resolved_threshold(n_traits)
        df = np.asarray(variant_df, dtype=np.float64)
        critical = np.abs(special.stdtrit(df, threshold / 2.0))
        return critical

    def select(self, beta, t, status, critical_t):
        """`(chunk, K)` device tensors in, four 1-D device tensors out.

        Returns `(variant_offset, trait_index, beta, t)`, all of length
        `n_surviving`, with `variant_offset` relative to the start of the
        chunk. The caller rebases it onto absolute variant numbering.

        One device-to-host synchronisation is unavoidable here and it is worth
        being explicit about why, because the rest of this scan is built to
        avoid exactly that: the number of survivors is a property of the data,
        so `nonzero` cannot know its own output shape without the device
        telling the host. It is one sync per chunk rather than per operation,
        and the payload that follows is a few rows rather than `chunk x K`, so
        it buys far more than it costs.
        """
        import torch

        if critical_t.dim() == 1:
            critical_t = critical_t.unsqueeze(1)
        keep = t.abs() >= critical_t
        keep &= (status == 0).unsqueeze(1)
        pairs = keep.nonzero(as_tuple=False)
        rows, cols = pairs[:, 0], pairs[:, 1]
        return rows, cols, beta[rows, cols], t[rows, cols]


# The one JAGWAS reduction (score statistic, FP64 correlation, rounding-target
# rank cutoff, block-triangular projection) lives in jagwas_projection.
from .jagwas_projection import JagwasReduction  # noqa: E402,F401
