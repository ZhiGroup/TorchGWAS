"""Decide before scanning whether a plan fits, and say what would.

**The failure this removes.** Supplementary Methods S3.1 lists as a limitation
that "users must currently split phenotype groups manually if the processed
matrix does not fit". That is what happens today on the unreduced path: there
is no feasibility check at all, so an over-large trait panel runs until CUDA
raises an out-of-memory error -- after the genotype read has already begun,
with a message that names neither the cause nor a workable setting.

Automatic trait blocking already exists, but only fires when a reduction or a
significance filter is requested, because blocking an UNREDUCED scan does not
save work: every (variant, trait) cell is emitted regardless, so the blocks
just re-read the genotypes. That reasoning is right about cost and wrong about
feasibility -- when the design does not fit, extra passes are not a waste, they
are the only way to run at all.

So this module answers three questions before anything is read:

    does the plan fit?
    if not, can it be made to fit automatically?
    if it cannot, what exactly would the user have to change?

The third is the point. "No manual tuning" does not mean the tool silently
guesses; it means the tool computes the answer instead of making the user
bisect for it. A refusal that names the largest workable trait count is
actionable in a way that `torch.cuda.OutOfMemoryError` is not.
"""
from __future__ import annotations

from .pipeline_model import auto_trait_block, predicted_peak_bytes


class PlanTooLarge(RuntimeError):
    """The requested scan cannot fit, and blocking cannot rescue it.

    Carries the numbers a caller needs to act, so the message is a plan rather
    than a complaint.
    """

    def __init__(self, message: str, *, predicted_bytes: float,
                 available_bytes: float, largest_traits: int | None):
        super().__init__(message)
        self.predicted_bytes = predicted_bytes
        self.available_bytes = available_bytes
        self.largest_traits = largest_traits


def _gigabytes(value: float) -> str:
    return f"{value / 1e9:,.1f} GB"


def largest_fitting_traits(*, chunk_variants: int, depth: int, n_samples: int,
                           covariate_rank: int,
                           transfer_bytes_per_variant: float,
                           device_memory_bytes: float,
                           n_traits: int,
                           headroom: float = 0.85,
                           decode_on_gpu: bool = False,
                           compute_log10_p: bool = False) -> int:
    """Largest trait count whose device rings fit, by bisection on the model.

    Bisects the same `predicted_peak_bytes` the planner uses, so the answer a
    refusal quotes is the answer the planner would accept -- rather than a
    rule of thumb that disagrees with it.
    """
    budget = device_memory_bytes * headroom

    def fits(traits: int) -> bool:
        peak = predicted_peak_bytes(
            chunk_variants=chunk_variants, depth=depth, n_samples=n_samples,
            n_traits=traits, covariate_rank=covariate_rank,
            transfer_bytes_per_variant=transfer_bytes_per_variant,
            decode_on_gpu=decode_on_gpu,
            compute_log10_p=compute_log10_p)
        return peak["gpu_bytes_per_device"] <= budget

    if fits(n_traits):
        return int(n_traits)
    if not fits(1):
        return 0
    low, high = 1, int(n_traits)
    while low < high:
        middle = (low + high + 1) // 2
        if fits(middle):
            low = middle
        else:
            high = middle - 1
    return low


def check_plan(*, chunk_variants: int, depth: int, n_samples: int,
               n_traits: int, covariate_rank: int,
               transfer_bytes_per_variant: float,
               device_memory_bytes: float,
               reduced: bool,
               trait_devices: int = 1,
               headroom: float = 0.85,
               decode_on_gpu: bool = False,
               compute_log10_p: bool = False) -> dict:
    """Report whether this plan fits, and what to do when it does not.

    `reduced` says whether a reduction or significance filter is active, which
    is what makes automatic trait blocking legal: blocked results can only be
    merged when they are reduced across traits. An unreduced scan emits every
    cell, so its blocks would each have to write their own slice of the output
    -- possible, but not something to do silently behind the caller's back.

    Returns a dict rather than raising, so a caller can decide whether an
    over-large plan is fatal or merely worth a warning. `require_fit` raises.
    """
    peak = predicted_peak_bytes(
        chunk_variants=chunk_variants, depth=depth, n_samples=n_samples,
        n_traits=n_traits, covariate_rank=covariate_rank,
        transfer_bytes_per_variant=transfer_bytes_per_variant,
        decode_on_gpu=decode_on_gpu, trait_devices=trait_devices,
        compute_log10_p=compute_log10_p)
    budget = device_memory_bytes * headroom
    fits = peak["gpu_bytes_per_device"] <= budget

    block = None
    if not fits and reduced:
        block = auto_trait_block(
            n_samples=n_samples, n_traits=n_traits,
            covariate_rank=covariate_rank, chunk_variants=chunk_variants,
            depth=depth, transfer_bytes_per_variant=transfer_bytes_per_variant,
            device_memory_bytes=int(device_memory_bytes), headroom=headroom)

    largest = None
    if not fits and not reduced:
        largest = largest_fitting_traits(
            chunk_variants=chunk_variants, depth=depth, n_samples=n_samples,
            covariate_rank=covariate_rank,
            transfer_bytes_per_variant=transfer_bytes_per_variant,
            device_memory_bytes=device_memory_bytes, n_traits=n_traits,
            headroom=headroom, decode_on_gpu=decode_on_gpu,
            compute_log10_p=compute_log10_p)

    return {
        "fits": fits,
        "predicted_device_bytes": peak["gpu_bytes_per_device"],
        "predicted_host_bytes": peak["host_pinned_bytes"],
        "budget_bytes": budget,
        "device_memory_bytes": device_memory_bytes,
        "auto_trait_block": block,
        "largest_fitting_traits": largest,
        "reduced": reduced,
    }


def require_fit(**kwargs) -> dict:
    """`check_plan`, but refuse an impossible plan BEFORE reading anything.

    The message names the shortfall and a setting that works. Failing here
    costs milliseconds; failing inside the scan costs however long the read
    had been running, and says only that some allocation did not succeed.
    """
    report = check_plan(**kwargs)
    if report["fits"] or report["auto_trait_block"]:
        return report

    predicted = report["predicted_device_bytes"]
    available = report["device_memory_bytes"]
    largest = report["largest_fitting_traits"]
    traits = kwargs["n_traits"]

    if largest and largest >= 1:
        advice = (
            f"reduce the trait count to {largest:,} or fewer, or pass "
            f"reduce= (top-k, max-abs-t or significant), which lets the scan "
            f"block the traits automatically and merge the blocked results")
    else:
        advice = (
            f"not even a single trait fits at chunk_variants="
            f"{kwargs['chunk_variants']:,}; lower the chunk size or use a "
            f"device with more memory")

    raise PlanTooLarge(
        f"this scan needs {_gigabytes(predicted)} on the device but only "
        f"{_gigabytes(available)} is present "
        f"({_gigabytes(report['budget_bytes'])} after headroom). "
        f"{traits:,} traits at {kwargs['n_samples']:,} samples: {advice}.",
        predicted_bytes=predicted, available_bytes=available,
        largest_traits=largest)
