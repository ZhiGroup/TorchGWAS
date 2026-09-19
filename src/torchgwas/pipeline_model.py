"""One resource-accounting model for all genotype input paths.

All bandwidths are bytes/second, rates are FLOP/second, and decode costs
are seconds/variant. Supply independent component measurements or operation
counts; never infer storage bandwidth or decode cost from a complete scan.

**CALIBRATE COMPONENTS INDEPENDENTLY; VALIDATE TERMS AND HELD-OUT TOTALS.**

Read this before "improving" the model's accuracy, because the structure has
essentially never been the problem -- the inputs have, and end-to-end checking
is what hid that. Every serious error found so far:

    symptom                              actual cause
    "the pipeline does not overlap,      the disk probe read 1.5 GB and
     1.45-1.84x is available"            reported 9.00 GB/s where the
                                         sustained rate is 6.4
    "output-bound above K = 179"         the write was priced at the
                                         single-stream fsync rate, 1.46 GB/s,
                                         against ~4.9 GB/s in the scan
    "the store is worth 7.63x"           there was no decode term at all
    "end to end is accurate, 0.92-0.97"  a disk rate 1.4x too high and a write
                                         rate 3.3x too low, cancelling at
                                         exactly the K that was checked

That last line is the trap. A single total can be right for two wrong reasons
and stay right for weeks; four separate terms cannot. Checked term by term on
corrected rates, the same model runs 1.00-1.01 across K = 1 to 512, and the old
end-to-end-only numbers run 0.51-1.12 -- a 96% worst error that one wall clock
never showed.

The scan publishes its own phase clocks, so this needs no new instrumentation.
Set `TORCHGWAS_SCAN_PROFILE=1` and read `genotype._last_scan_profile`:
`fetch_seconds` (read plus host decode), `gpu_compute_milliseconds` and
`gpu_result_milliseconds` (CUDA events, so they time the kernel and not the
launch), `setup_seconds`, `result_copy_seconds`, `result_wait_seconds`.
`benchmarks/direct_term_by_term.py` puts each model term beside the clock that
measures that term and nothing else, and prints "no clock" for a term nothing
measures rather than folding it into a residual.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import itertools
import json
import math


# Binary sumstats always hold float32 t_stat and neg_log10_p, with optional
# float32 beta, so output volume is a property of the requested fields.
SUMSTATS_BYTES_PER_TEST = {"beta+t": 12.0, "t": 8.0}


def sumstats_output_bytes_per_test(fields: str = "beta+t") -> float:
    """Output bytes per marker-trait test for a binary sumstats field set."""
    try:
        return SUMSTATS_BYTES_PER_TEST[fields]
    except KeyError:
        raise ValueError(
            f"unknown sumstats field set {fields!r}; expected one of "
            f"{sorted(SUMSTATS_BYTES_PER_TEST)}"
        ) from None


MIN_AUTO_CHUNK_VARIANTS = 128

# Measured, not assumed, and it is the opposite of what this module used to
# believe -- though not in the simple way first recorded. The curve below only
# reaches down to 2,048 and so looked monotone; a later sweep on a quiet box,
# two rounds within 2%, reaching down to 256, shows it is U-SHAPED, with the
# optimum at 1,024 (see the second table further down). "Monotonically worse as
# the chunk grows" is true only above the optimum, and is retracted as a
# general claim.
#
# Cold, M = 200,000, pgen hard calls at 35,365 samples and depth 32:
#
#   chunk    chunks (200k)   scan     peak GPU
#    2,048        97         2.00 s     1.4 GB
#    4,096        48         2.33 s     2.8 GB
#    8,192        24         3.12 s     5.1 GB
#   13,267        15         3.20 s     8.0 GB
#   26,624         7         4.00 s    15.6 GB
#   65,536         3        23.28 s    37.6 GB
#
# Same shape at 1M variants. So memory is a CEILING, not the objective.
# `auto_chunk_variants` still bisects `device_ring_bytes` to guarantee the ring
# fits, but the answer is capped here, because past this point larger is simply
# worse. The cap was 65,536, which the corrected memory model was happy to hand
# out -- and doing so made scans 3.5x slower.
#
# **Why 4096 and not another point on that curve.** `benchmarks/
# direct_chunk_curve_check.py` compares each row above against `explain_time`
# plus the separately measured 2.10 s pgen startup. Up to 4096 the model is the
# whole story; past it a residual appears that the model has no term for:
#
#   chunk    measured   model+startup   residual
#    2,048     2.00 s       2.30 s       -0.30 s
#    4,096     2.33 s       2.31 s       +0.02 s
#    8,192     3.12 s       2.31 s       +0.81 s
#   26,624     4.00 s       2.35 s       +1.65 s
#   65,536    23.28 s       2.44 s      +20.84 s
#
# The cap is therefore near where the model stops being valid, which is a
# better reason to put it here than "larger measured slower". Earlier scan-only
# checks were called unresolved because they omitted the startup intercept and
# so came out uniformly ~10x low.
#
# The QUIET-BOX curve is the better dataset and is less flattering. Two rounds
# within 2%, M = 1,000,000, K = 128, warm cache (so the model is run with the
# disk term dropped, making it a lower bound):
#
#   chunk   chunks   measured    model   residual
#     256    3,906     3.32 s    2.32 s   +1.00 s
#     512    1,953     2.50 s    2.32 s   +0.18 s
#   1,024      977     2.10 s    2.32 s   -0.22 s   <-- measured optimum
#   2,048      488     2.23 s    2.32 s   -0.09 s
#   4,096      244     2.67 s    2.32 s   +0.35 s   <-- the cap
#   8,192      122     3.29 s    2.32 s   +0.97 s
#
# Every modelled term except fill/drain is flat in chunk, and fill/drain is
# milliseconds at these counts, so the model predicts a HORIZONTAL LINE against
# a measured U. It is within ~10% near the bottom (1,024-2,048) and misses both
# arms. The low arm is not one fixed cost per chunk either -- residual over
# chunk count gives 0.256 ms then 0.092 ms, not a constant.
#
# So the cap sits 1.26x off the measured optimum and uses 2.9x the memory.
# Moving it to 1,024 is the indicated change and is NOT made here: every
# published timing was taken at 4,096, so it invalidates the comparison matrix
# and is a deliberate call, not a side effect of a calculator fix.
#
# This module used to explain the penalty as lost overlap ("three chunks, so
# fill/drain is the whole run"). That is arithmetically wrong and is retracted:
# fill/drain is `total / chunks`, and at chunk 65,536 that is 0.14 s against a
# 20.84 s residual. Whatever the mechanism is, it is not the pipeline running
# out of chunks to overlap.
#
# The residual is NOT identified. Pinning the host ring at the measured
# 1.91 GB/s is the right order of magnitude but the wrong shape -- it
# over-predicts 2.2-3.5x between 8,192 and 26,624 and under-predicts at
# 65,536 -- so it is not one clean mechanism, and the 5x jump in seconds per GB
# at the top row says two regimes rather than one curve. Anything above this
# cap is extrapolation, and `explain_time` says so in
# `chunk_within_validated_range` rather than predicting silently.
MAX_AUTO_CHUNK_VARIANTS = 4096

# Above this chunk the measured scan departs from every term the model has.
#
# No longer equal to the cap. It used to be, with the note "the cap is placed
# at the edge of validity" -- but the edge has since been measured further out.
# Four round-robin rounds on an idle H100 (M=8,931,083, N=35,365, K=128,
# sumstats=none) put chunks 8,192 and 16,384 at 4.63 s and 4.75 s against
# 4.68 s at 4,096, all three inside a within-chunk spread of 1.38x, and the
# model predicts them 1.15x and 1.13x -- the same accuracy it has at the cap.
# So the model is checked to 16,384; the unmodelled residual that motivated the
# ceiling was measured at 26,624 and 65,536, both still above this.
#
# The CAP stays at 4,096 regardless, and not because larger is worse: it is
# not, it is flat. Larger simply buys nothing measurable while costing device
# memory in proportion, and the ring is the scarcer resource.
VALIDATED_MAX_CHUNK_VARIANTS = 16_384


def device_ring_bytes(*, chunk_variants: int, depth: int, n_samples: int,
                      n_traits: int, covariate_rank: int,
                      transfer_bytes_per_variant: float,
                      decode_tile: int | None = None,
                      decode_on_gpu: bool = False,
                      fused_packed_statistics: bool = False,
                      gpu_decoded_bytes_per_variant: float = 0.0) -> float:
    """Device memory the scan's rings occupy, term by term.

    ONE formula, used both to *predict* a plan and to *choose* a chunk. There
    used to be two: this accounting, and a separate per-variant estimate that
    explicitly omitted "padding and workspace" and so had to be paired with a
    flat quarter-of-memory budget. Two models of the same quantity drift, and
    the cheaper one was the one actually wired into the scan.

    **Calibrated against measured peaks, not assumed.** The first version
    over-predicted by 2x to 9x. It charged a `depth`-deep ring of *decoded
    float32* genotype, which the native fused path does not have -- `prepare`
    centres one chunk at a time -- and for a float32 transport that
    double-counted the staging ring outright. It also counted the per-variant
    result cells on the device, when `native_scan` allocates them with
    `pin_memory=True`, on the host. Against the K sweep the old form predicted
    71 GB for every format while the real peaks were 8.0 GB (pgen hard calls),
    16.2 GB (zstd) and 35.2 GB (pgen dosage) -- so the planner was choosing
    chunks up to 9x smaller than the card could hold. This form reproduces
    those peaks to within about 5%.
    """
    tile = chunk_variants if decode_tile is None else max(int(decode_tile), 1)
    n, k, c = n_samples, n_traits, covariate_rank

    # The staging ring: `depth` buffers of the TRANSFERRED representation.
    # This is the dominant term and the one that differs by 16x between a
    # packed two-bit transport and float32 dosage.
    total = depth * tile * transfer_bytes_per_variant

    # The centred float32 genotype the GEMM reads. Produced one chunk at a
    # time by `prepare`, not held `depth` deep -- allow two for the overlap
    # between the chunk being centred and the one being multiplied. A
    # device-native decoder inflates into its own buffers, so it carries an
    # extra slot.
    # A FUSED PACKED kernel never materialises this at all: it reads two-bit
    # codes straight out of the staging buffer and centres on the fly, so there
    # is no float32 genotype chunk to hold. Charging it anyway is what made the
    # model predict 2.41 GB for a BED scan that measured 0.803 GB -- a 3x
    # over-prediction, and the wrong direction, because it would refuse chunk
    # sizes the card can actually take.
    #
    # Measured, BED at chunk 4096 / depth 18 / 8,842 B per variant: staging
    # 651.9 MB + design and products 24 MB = 676 MB predicted against 803 MB
    # measured. The residual ~127 MB is statistics-kernel workspace that has
    # NOT been itemised yet -- read it out of the kernel rather than folding it
    # into a coefficient here.
    # Three cases, not two:
    #   * fused packed kernel -- reads two-bit codes directly, so no float32
    #     genotype chunk exists at all (BED).
    #   * device decoder with an output ring -- the decoded float32 output IS
    #     the genotype the GEMM reads, already charged as
    #     `depth * chunk * decoded_bytes`. Charging a centred copy too
    #     double-counts: BGEN came out at 6.379 G against a measured 4.635 G,
    #     over by exactly the 1.744 G of a three-slot centred term.
    #   * host-staged -- `prepare` centres into its own buffers, two slots to
    #     overlap the chunk being centred with the one being multiplied.
    if fused_packed_statistics or gpu_decoded_bytes_per_variant:
        centred_slots = 0.0
    else:
        centred_slots = 3.0 if decode_on_gpu else 2.0
    total += centred_slots * max(chunk_variants, tile) * n * 4.0

    # Row means and the per-chunk products, then the design, which is the only
    # term that does not scale with the chunk.
    total += chunk_variants * 4.0
    total += chunk_variants * (k + c) * 4.0
    total += n * (k + c) * 4.0

    # A DEVICE-NATIVE decoder has no HOST staging ring -- it reports
    # `transfer_bytes_per_variant: 0` -- but it still has a ring: the decoded
    # float32 OUTPUT, `depth` deep, which the consumer drains. That term was
    # missing entirely, which is why the model predicted 1.70 GB against a
    # measured 5.35 GB.
    #
    # Measured exactly, at N=35,365, chunk 4096, prefetch 8:
    #     8 * 4096 * 35365 * 4 = 4,634,951,680 B = 4.635 G
    # against a measured peak of 4.635 G -- to the byte. The ring depth is what
    # sets BGEN's device memory, not the decoder's own buffers.
    #
    # An earlier guess blamed nvCOMP's scratch. Instrumenting the decoder to
    # report `nvcompBatchedDeflateDecompressGetTempSizeAsync` showed it is
    # **zero** for this configuration, so that hypothesis is dead and this term
    # is what replaced it.
    #
    # BGEN Layout-2 at 16 bits stores two probabilities per sample per
    # biallelic variant, so the decompressed width is n * 2 values * 2 bytes =
    # 4n = 141,460 B/variant at N=35,365 -- the SAME width as float32 dosage,
    # just inflated on the device instead of the host.
    #
    # NOT INCLUDED, deliberately: nvCOMP's own scratch, sized at runtime by
    # `nvcompBatchedDeflateDecompressGetTempSizeAsync`. It is opaque and must
    # be measured, not guessed; leaving it out means this term is a floor for
    # a device-decoded format and the model will under-predict by whatever
    # nvCOMP reserves.
    if gpu_decoded_bytes_per_variant:
        total += depth * chunk_variants * gpu_decoded_bytes_per_variant
    return total


def host_pinned_bytes(*, chunk_variants: int, depth: int, n_traits: int,
                      transfer_bytes_per_variant: float,
                      reduction_width: int | None = None,
                      stage_on_host: bool = True,
                      compute_log10_p: bool = False) -> float:
    """Pinned host memory the scan holds, term by term.

    The device ring is not the only ring. Every slot is pinned on BOTH sides:
    `streaming` pins a staging buffer per slot, and `native_scan` pins a
    result buffer per slot. Pinned pages are not swappable, so this is a hard
    floor on resident memory, and at voxel scale it is the term that dies
    first -- the result ring carries `chunk * K` cells, so at depth 32,
    chunk 4096 and K = 2,085,000 it asks for 2.2 PB. That is the arithmetic
    making trait blocking mandatory rather than an optimisation, and it was
    previously nowhere in the model: `device_ring_bytes` deliberately excludes
    these buffers because they are not on the card, which left the host side
    unaccounted for entirely.
    """
    chunk = max(int(chunk_variants), 1)
    total = 0.0

    # The staging ring, mirrored on the host by `streaming.PinnedLoader`. A
    # device-native source stages nothing here.
    if stage_on_host:
        total += depth * chunk * transfer_bytes_per_variant

    # The result ring. Unreduced: float32 beta and t per marker-trait cell,
    # optional float64 -log10(P), plus a uint8 status and a float32 residual
    # df per VARIANT. Reduced: the narrow width carries an extra int32
    # trait-index column, and the whole point of the mode is that `width`
    # replaces `K` here.
    if reduction_width is None:
        result_bytes_per_test = 16.0 if compute_log10_p else 8.0
        total += depth * chunk * (result_bytes_per_test * n_traits + 5.0)
    else:
        total += depth * chunk * (12.0 * max(int(reduction_width), 1) + 5.0)
    return total


def predicted_peak_bytes(*, chunk_variants: int, depth: int, n_samples: int,
                         n_traits: int, covariate_rank: int,
                         transfer_bytes_per_variant: float,
                         decode_tile: int | None = None,
                         decode_on_gpu: bool = False,
                         reduction_width: int | None = None,
                         trait_block: int | None = None,
                         trait_devices: int = 1,
                         compute_log10_p: bool = False) -> dict:
    """Both peaks for one resolved plan, as the benchmark tables report them.

    `trait_block` is what actually sits on a device at once; `n_traits` is the
    full trait axis. When the traits are blocked, the device and result rings
    are sized by the block, and sharding the blocks across `trait_devices`
    cards divides the per-card device peak but NOT the host peak, which is
    shared -- each shard pins its own result ring in the same address space.
    """
    block = int(trait_block) if trait_block else int(n_traits)
    block = max(min(block, int(n_traits)), 1)
    devices = max(int(trait_devices), 1)
    gpu = device_ring_bytes(
        chunk_variants=chunk_variants, depth=depth, n_samples=n_samples,
        n_traits=block, covariate_rank=covariate_rank,
        transfer_bytes_per_variant=transfer_bytes_per_variant,
        decode_tile=decode_tile, decode_on_gpu=decode_on_gpu)
    host = host_pinned_bytes(
        chunk_variants=chunk_variants, depth=depth, n_traits=block,
        transfer_bytes_per_variant=transfer_bytes_per_variant,
        reduction_width=reduction_width,
        stage_on_host=not decode_on_gpu,
        compute_log10_p=compute_log10_p) * devices
    return {
        "gpu_bytes_per_device": gpu,
        "host_pinned_bytes": host,
        "trait_block": block,
        "trait_devices": devices,
        "trait_passes": math.ceil(int(n_traits) / block),
    }


def explain_memory(*, chunk_variants: int, depth: int, n_samples: int,
                   n_traits: int, covariate_rank: int,
                   transfer_bytes_per_variant: float,
                   decode_tile: int | None = None,
                   decode_on_gpu: bool = False,
                   fused_packed_statistics: bool = False,
                   gpu_decoded_bytes_per_variant: float = 0.0,
                   device_memory_bytes: int | None = None) -> dict:
    """Device memory broken into its terms, with the binding one named.

    `device_ring_bytes` returns one number, which is what `auto_chunk_variants`
    needs. This returns the same arithmetic DECOMPOSED, because a single number
    cannot answer the question that actually gets asked: *what do I change?*

    The terms differ by orders of magnitude between formats and the largest one
    is not the same in each, so 'reduce memory' has a different answer per
    format. Measured at N=35,365, chunk 4096, depth 32: pgen dosage is 94%
    staging ring (float32 dosage at 141,460 B/variant), while a fused packed
    BED scan has no centred term at all. Telling a caller only the total hides
    exactly the fact that would let them act.

    `headroom_ratio` is against the card when `device_memory_bytes` is given:
    below 1.0 the plan does not fit.
    """
    tile = chunk_variants if decode_tile is None else max(int(decode_tile), 1)
    n, k, c = n_samples, n_traits, covariate_rank
    # Three cases, not two:
    #   * fused packed kernel -- reads two-bit codes directly, so no float32
    #     genotype chunk exists at all (BED).
    #   * device decoder with an output ring -- the decoded float32 output IS
    #     the genotype the GEMM reads, already charged as
    #     `depth * chunk * decoded_bytes`. Charging a centred copy too
    #     double-counts: BGEN came out at 6.379 G against a measured 4.635 G,
    #     over by exactly the 1.744 G of a three-slot centred term.
    #   * host-staged -- `prepare` centres into its own buffers, two slots to
    #     overlap the chunk being centred with the one being multiplied.
    if fused_packed_statistics or gpu_decoded_bytes_per_variant:
        centred_slots = 0.0
    else:
        centred_slots = 3.0 if decode_on_gpu else 2.0
    terms = {
        # `depth` buffers of the TRANSFERRED representation. The term that
        # differs 16x between packed two-bit and float32 dosage.
        'staging_ring': depth * tile * transfer_bytes_per_variant,
        # float32 genotype the GEMM reads; absent for a fused packed kernel.
        'centred_genotype': centred_slots * max(chunk_variants, tile) * n * 4.0,
        # Row means, one per variant in the chunk.
        'row_means': chunk_variants * 4.0,
        # Per-chunk products against the design.
        'chunk_products': chunk_variants * (k + c) * 4.0,
        # The design itself -- the only term that does not scale with chunk,
        # and the one that decides trait blocking at voxel-scale K.
        'design': n * (k + c) * 4.0,
        # A device decoder's ring of decoded float32 output. Zero for
        # host-staged formats; for BGEN it is the whole story.
        'decoded_output_ring': depth * chunk_variants * gpu_decoded_bytes_per_variant,
    }
    total = sum(terms.values())
    binding = max(terms, key=terms.get)
    out = {
        'terms_bytes': terms,
        'total_bytes': total,
        'binding_term': binding,
        'binding_fraction': terms[binding] / total if total else 0.0,
        'share': {name: (value / total if total else 0.0)
                  for name, value in terms.items()},
    }
    if device_memory_bytes:
        out['headroom_ratio'] = device_memory_bytes / total
    return out


def memory_sensitivity(*, chunk_variants: int, depth: int, n_samples: int,
                       n_traits: int, covariate_rank: int,
                       transfer_bytes_per_variant: float,
                       decode_on_gpu: bool = False,
                       fused_packed_statistics: bool = False) -> dict:
    """What halving each knob actually buys, in bytes.

    The point of a decision tool is to rank the levers, and the ranking is not
    obvious: halving `depth` and halving `chunk` both halve the staging ring,
    but only `chunk` also halves the centred genotype, while NEITHER touches
    the design -- so at voxel-scale K, where the design dominates, both are
    nearly useless and only trait blocking helps. This makes that explicit
    rather than leaving it to be rediscovered.
    """
    base = dict(chunk_variants=chunk_variants, depth=depth, n_samples=n_samples,
                n_traits=n_traits, covariate_rank=covariate_rank,
                transfer_bytes_per_variant=transfer_bytes_per_variant,
                decode_on_gpu=decode_on_gpu,
                fused_packed_statistics=fused_packed_statistics)
    reference = device_ring_bytes(**base)
    out = {}
    for knob, changed in (
            ('halve_chunk', dict(chunk_variants=max(chunk_variants // 2, 1))),
            ('halve_depth', dict(depth=max(depth // 2, 1))),
            ('halve_traits', dict(n_traits=max(n_traits // 2, 1))),
            ('halve_transport_width',
             dict(transfer_bytes_per_variant=transfer_bytes_per_variant / 2.0)),
    ):
        trial = dict(base); trial.update(changed)
        after = device_ring_bytes(**trial)
        out[knob] = {'bytes_after': after, 'bytes_saved': reference - after,
                     'fraction_saved': (reference - after) / reference
                     if reference else 0.0}
    out['reference_bytes'] = reference
    return out

# Achieved FP32 is NOT a constant -- it is a strong function of the design
# width, and treating it as a scalar is a real source of error. Measured on the
# H100 with CUDA events around the statistics kernel, 400,000 variants:
#
#     K      design width   achieved TFLOPS
#     8           36              8.35
#    32           60              8.71
#   128          156             20.57
#   512          540             20.28
#
# A 2.5x spread. The roofline says exactly this: arithmetic intensity is about
# (K + C) / 2, so a narrow design is memory-bound and a wide one approaches the
# FLOP limit. Using one number -- the model used 44.88 TFLOPS, calibrated at
# width 540 -- over-estimates the GPU roughly 2x at K = 128 and 5x at K = 8.
# It was invisible end to end because the GEMM rarely binds.
#
# Two caveats on the values above, both pointing the same way: they were taken
# at host load 155 with a competitor benchmark running, and the CUDA-event
# bracket includes any wait on the input copy, so they are LOWER bounds. The
# SHAPE is what matters here and it matches an independent earlier note
# (~5.6 TFLOPS at 155 columns, 29.3 at 540). Re-take on a quiet host before
# quoting the absolute numbers.
# RE-TAKEN QUIET, 2026-09-15, which the note above had been asking for. Idle
# H100, TF32 explicitly DISABLED (a silent TF32 promotion would be compared
# against the wrong roofline), chunk 4,096, samples 35,365, min of 7, by
# `benchmarks/direct_gemm_realized_fraction.py`:
#
#     width     36    60   156    540   1052   2076   4124
#     TFLOP/s 25.6  42.4  38.2   26.1   29.7   26.5   22.7
#
# Three things this changes. First, the old values were 2-4x LOW, exactly as
# their own note predicted ("taken at host load 155 ... LOWER bounds"). Second,
# the table now reaches width 4,124 instead of clamping at 540, so K=1024 is no
# longer priced at the K=512 rate. Third, the shape is NOT monotonic, and the
# reason is measurable rather than mysterious:
#
#   THE DESIGN MATRIX IS THE REUSED OPERAND, AND IT EITHER FITS IN L2 OR DOES
#   NOT. At 35,365 samples the design is `samples * width * 4` bytes, against a
#   driver-reported 52.4 MB L2:
#
#       width 36 -> 5.1 MB, 60 -> 8.5, 156 -> 22.1   RESIDENT, 25.6-42.4 TFLOP/s
#       width 540 -> 76.4 MB, 1052 -> 148.8, 4124 -> 583.4   NOT, 22.7-29.7
#
#   Every chunk row multiplies against that same design, so when it is resident
#   the reads come from L2 and when it is not they come from HBM. Confirmed
#   independently on square products, where the achieved rate climbs to
#   51.4 TFLOP/s at the last size whose working set fits (50.3 MB) and drops to
#   41.4 at the first that does not (63.7 MB).
#
# ALSO MEASURED, and not expressible as a rate at all: a fixed per-call cost of
# ~0.31 ms. Time is flat from width 2 to 28 and flat again from 36 to 60 while
# the FLOPs grow 14x, so at low K this product is LATENCY-bound. The old
# table's 8.35 TFLOP/s at width 36 was that floor expressed as a fake rate.
# `benchmarks/direct_gemm_l2_and_floor.py` has both sweeps.
#
# Prefer `gemm_model.DeviceSpec` over this table where portability matters:
# this is one card, and the derived bound answers for cards not present.
MEASURED_GEMM_TFLOPS_BY_WIDTH = {36: 25.6, 60: 42.4, 156: 38.2, 540: 26.1,
                                 1052: 29.7, 2076: 26.5, 4124: 22.7}

# Superseded, kept so the retraction is legible rather than a silent edit.
CONTENDED_GEMM_TFLOPS_BY_WIDTH_RETRACTED = {36: 8.35, 60: 8.71,
                                            156: 20.57, 540: 20.28}


def gemm_rate_at_width(rate, design_width: float,
                       chunk_variants: int | None = None,
                       realized_fraction: float | None = None) -> float:
    """Achieved FLOP/s for this design width.

    `rate` is either a scalar -- used as-is, which is what every caller did
    before this existed -- or a mapping from design width to FLOP/s, which is
    interpolated. A mapping is strongly preferred: see the table above.

    Linear interpolation between the two nearest measured widths, and the
    nearest endpoint outside the measured range. Not a fit: no curve is
    assumed, and outside the measured span it refuses to extrapolate a trend
    it has no evidence for.
    """
    # A DeviceSpec means "derive it" -- see `gemm_model`. This is the
    # first-principles path and the one that does not clamp: SM count, FP32
    # lanes and clock give the peak, per-tile arithmetic intensity caps it, and
    # tile/wave quantization reduces it. Preferred over the measured table
    # below, which holds four widths and returns the width-540 rate for
    # anything wider -- a clamp worth 1.81x at K=1024 on the K ladder.
    from .gemm_model import DeviceSpec, gemm_flops_per_second

    if isinstance(rate, DeviceSpec):
        if chunk_variants is None:
            raise ValueError(
                "deriving a GEMM rate needs the chunk: the bound depends on "
                "both output dimensions, since tile and wave quantization are "
                "what separate it from peak")
        return gemm_flops_per_second(rate, chunk_variants, int(design_width),
                                     realized_fraction=realized_fraction)
    if not isinstance(rate, dict):
        return float(rate)
    if not rate:
        raise ValueError("empty gemm rate table")
    widths = sorted(float(w) for w in rate)
    lookup = {float(w): float(v) for w, v in rate.items()}
    width = float(design_width)
    if width <= widths[0]:
        return lookup[widths[0]]
    if width >= widths[-1]:
        return lookup[widths[-1]]
    for low, high in zip(widths, widths[1:]):
        if low <= width <= high:
            if high == low:
                return lookup[low]
            share = (width - low) / (high - low)
            return lookup[low] + share * (lookup[high] - lookup[low])
    return lookup[widths[-1]]


def finite_pipeline_seconds(resource_seconds, variants, chunk_variants):
    """Flow-shop makespan for uniform per-row service and a partial last chunk.

    Each chunk traverses the stages in order; each stage handles one chunk
    at a time. For q equal chunks, completion at stage j is the prefix sum
    of services plus (q-1) times the prefix maximum. Append the partial chunk
    with C[i,j] = max(C[i-1,j], C[i,j-1]) + service[i,j]. No fitted margin.
    """
    q, tail = divmod(variants, chunk_variants)
    per_full = [t * chunk_variants / variants for t in resource_seconds]
    complete, prefix, bottleneck = [], 0.0, 0.0
    for service in per_full:
        prefix += service
        bottleneck = max(bottleneck, service)
        complete.append(prefix + (q-1)*bottleneck if q else 0.0)
    if tail:
        previous_stage = 0.0
        for j, service in enumerate(per_full):
            previous_stage = max(previous_stage, complete[j]) + service*tail/chunk_variants
            complete[j] = previous_stage
    return complete[-1] if complete else 0.0


def explain_time(*, variants: int, samples: int, traits: int,
                 covariate_rank: int, chunk_variants: int,
                 stored_bytes: float, transfer_bytes_per_variant: float,
                 output_bytes_per_test: float,
                 disk_bytes_per_second: float, h2d_bytes_per_second: float,
                 d2h_bytes_per_second: float,
                 # float, or {design_width: FLOP/s} -- see
                 # `gemm_rate_at_width`. A scalar is wrong at every
                 # width but the one it was measured at.
                 gemm_flops_per_second,
                 write_bytes_per_second: float,
                 metadata_bytes: float = 0.0,
                 text_parse_bytes_per_second: float = 0.0,
                 decoded_bytes_per_variant: float = 0.0,
                 decode_bytes_per_second: float = 0.0,
                 overlap: float = 1.0,
                 write_overlap: float = 0.0,
                 # DEFAULTS TO OFF, and that is deliberate. Derived at ~453 us
                 # on an idle H100, but several callers calibrate with ACHIEVED
                 # rates -- the anchor test sets `disk_bytes_per_second` to
                 # `stored_bytes / 14.21` from a measured scan -- and an
                 # achieved rate already contains this overhead, so switching
                 # it on there would double-count. Pass it explicitly alongside
                 # micro-benchmarked rates.
                 #
                 # Treat the VALUE as provisional. It was derived from two
                 # points that disagreed with each other (479 us at chunk 512
                 # against 661 us at 4,096) and averaged, which is fitting. The
                 # non-GEMM GPU time it was meant to capture is 8.35, 6.42,
                 # 1.53, 1.44, 1.45 s at chunks 512..8,192 -- flat from 2,048
                 # up and then a cliff, which no `a * chunks` term reproduces.
                 per_chunk_seconds: float = 0.0,
                 # Variants per frame in the source store, when it is framed.
                 # See the straddling term below: this is what makes a chunk
                 # of 1,536 cost 4.6x a chunk of 2,048.
                 frame_variants: int = 0,
                 # Only meaningful when `gemm_flops_per_second` is a DeviceSpec
                 # and the rate is being DERIVED. It is the device's realized
                 # share of its own roofline -- one machine number, measured
                 # once with a large square GEMM. None returns the bound, which
                 # is optimistic by about 1.4x on an H100 but is honest about
                 # being a bound. There is deliberately no default.
                 gemm_realized_fraction: float | None = None) -> dict:
    """Scan and end-to-end seconds, decomposed, with the binding term named.

    Rates come from `calibrate.py` -- measured on the host, never looked up
    from a spec sheet and never fitted to a scan. Feeding peak FP32 in is what
    made a hand-estimate predict 60-100 s for a run that measured 24.26 s.

    Three parts, because a single number hides which one to attack:

      steady state -- the slowest resource once the ring is full. This is the
        term that shrinks when you fix the bottleneck, and nothing else.
      fill and drain -- service not already counted in the bottleneck total.
        For equal chunks it is `(sum(resources)-max(resources))/chunks`;
        a partial final chunk uses the flow-shop recurrence. Leaving it out is why scan-only predictions came out ~1.5x
        low while end-to-end came out right.
      fixed -- the variant-metadata parse, which is flat in M. At M=1,000,000
        it is 76% of an end-to-end run and at full genome under 10%, so a model
        validated only at full scale never sees it.

    `overlap` is 1.0 for a perfectly overlapped pipeline (cost is the slowest
    stage) and 0.0 for a fully serialised one (cost is the sum). It is a
    STRUCTURAL parameter, not a fudge.

    **THE VALUE 0.93 IS A FIT, AND IT IS WORSE THAN NOT FITTING.** Read the
    ladder below as what it is: `measured`, `max` and `sum` are used to SOLVE
    for the overlap that reproduces the clock. That is inverting the model
    against wall time, which is the one thing a cost model here may not do.

    Re-checked after the decode rate was measured rather than inferred, the
    unfitted structural value 1.0 predicts BETTER at every point tested, on
    both regimes (`_scratch/_fit_audit.sh`, `_scratch/_overlap_regimes.sh`):

        zstd store, aligned chunks   1.0 -> 0.99-1.04x    0.93 -> 1.10-1.15x
        raw .bed K ladder            1.0 -> 1.00/0.98/1.13/1.81
                                     0.93 -> 1.02/1.01/1.21/1.88

    So 0.93 was never measuring overlap. It was absorbing the error in the
    inferred 11.3 GB/s decode rate -- both push the prediction the same way --
    and once that rate was measured at 20.01 GB/s the fit became dead weight.
    Use 1.0. What the ladder below DOES establish, and what survives, is the
    structural claim that the pipeline overlaps at all.

    (The K=1024 row is 1.81x off under either value. That is not overlap: it is
    the GEMM width table clamping above width 540 -- see
    `MEASURED_GEMM_TFLOPS_BY_WIDTH`, whose own note says it was taken at host
    load 155 and is a lower bound awaiting a quiet re-take.)

    The original derivation, kept because the structural argument is sound:
    a K ladder at full genome with the
    sumstats write disabled, so the only terms growing with K are the GEMM and
    the D2H copy while the 79.0 GB read stays fixed
    (`benchmarks/direct_overlap_from_kladder.py`):

        K      measured    max     sum    implied overlap
        1        14.21s  14.21s  16.05s        1.00
        128      14.52s  14.21s  18.18s        0.92
        512      14.90s  14.21s  24.64s        0.93
        1024     18.09s  14.81s  33.26s        0.82

    From K=1 to K=512 the compute grows by 8.60 s and the scan by 0.69 s --
    **8% of it reaches the clock**, which serialisation cannot produce. At
    K=1024 the compute finally overtakes the read and the curve turns up,
    exactly where the overlapped model says it should. So the default of 1.0 is
    a mild optimism, not an unverified guess.

    **This overturns an earlier conclusion.** The M sweep was read as showing
    the pipeline did NOT overlap, with 1.45-1.84x available from fixing it.
    That compared the measured per-variant slope against a `max` bound computed
    with the disk at 9.00 GB/s -- the 16-reader microbenchmark peak, which the
    scan does not reach; it achieves 5.56 GB/s. Against the achieved rate the
    same measurement sits at the max, not the sum. sum-at-peak-rate happens to
    land near max-at-achieved-rate, and reading that coincidence as evidence of
    serialisation invented a headroom that is not there.
    """
    m, n, k, c = variants, samples, traits, covariate_rank
    for name, value in (("variants", m), ("samples", n), ("traits", k),
                        ("chunk_variants", chunk_variants)):
        if isinstance(value, bool) or int(value) != value or value <= 0:
            raise ValueError(f"{name} must be a positive integer")
    if not 0 <= overlap <= 1 or not 0 <= write_overlap <= 1:
        raise ValueError("overlap values must be in [0, 1]")
    chunks = math.ceil(m / chunk_variants)
    result_bytes_per_test = 16.0 if output_bytes_per_test else 8.0
    resources = {
        'disk': stored_bytes / disk_bytes_per_second,
        'h2d': m * transfer_bytes_per_variant / h2d_bytes_per_second,
        'gemm': (m * 2.0 * n * (k + c + 1))
                / gemm_rate_at_width(gemm_flops_per_second, k + c + 1,
                                     chunk_variants=chunk_variants,
                                     realized_fraction=gemm_realized_fraction),
        # Dense output also returns float64 -log10(P) from the GPU. It is
        # stored as float32, but the device-to-host transfer preserves the
        # calculation before the writer performs that final cast.
        'd2h': m * (result_bytes_per_test * k + 5.0) / d2h_bytes_per_second,
    }
    # Host-side decompression, for a compressed store. Priced on the bytes it
    # PRODUCES, because that is what the codec's rate is quoted against and
    # what actually has to be written somewhere.
    #
    # This term did not exist, and its absence is not academic: without it the
    # model predicted a 7.63x speedup for the zstd hard-call store where the
    # measurement gives 1.84x. Compression trades disk bytes for CPU work and
    # a model with no CPU term can only ever see the saving.
    # FRAME STRADDLING. A framed store decompresses whole frames, so a chunk
    # that does not line up with the frame grid decodes rows it discards. The
    # amplification is exact arithmetic, not a fitted factor: over one
    # `lcm(chunk, frame)` period there are `period/chunk` chunks and
    # `period/frame` frame boundaries, each boundary landing inside exactly one
    # chunk except the single coincident one, so
    #
    #     frames touched = period/chunk + period/frame - 1
    #     amplification  = frame * that / period = 1 + (frame - gcd) / chunk
    #
    # Measured on an idle H100 against a frame of 2,048 (four round-robin
    # rounds, M=8,931,083, N=35,365, K=128, sumstats=none), predicted vs the
    # penalty relative to the best aligned chunk:
    #
    #     chunk    1024   1536   2048   3072   4096   8192  16384
    #     predicted 2.00   2.00   1.00   1.33   1.00   1.00   1.00
    #     measured  2.48   4.64   1.03   1.45   1.00   0.99   1.01
    #
    # It classifies all seven correctly and gets the magnitude right from
    # 2,048 up. It UNDER-predicts the two chunks smaller than one frame --
    # 1,536 costs 4.64x where this says 2.00x. That residual is real and is
    # NOT modelled here; inventing a coefficient to close it would be fitting.
    # What is known is that the excess appears in GPU time as well as decode
    # (gpu_compute 9.60 s at chunk 1,536 against 4.37 s at 2,048), which
    # straddling alone cannot explain, so it is a second mechanism and wants
    # its own measurement. In practice `auto_chunk_variants` now aligns the
    # chunk it hands out, so a caller never reaches the under-predicted region
    # by accident -- the residual is bounded to chunks a user forces by hand.
    amplification = 1.0
    if frame_variants and frame_variants > 0 and chunk_variants > 0:
        from math import gcd
        amplification = 1.0 + (frame_variants
                               - gcd(int(chunk_variants), int(frame_variants))
                               ) / float(chunk_variants)
    if decoded_bytes_per_variant and decode_bytes_per_second:
        resources['decode'] = (m * decoded_bytes_per_variant * amplification
                               / decode_bytes_per_second)
    binding = max(resources, key=resources.get)
    steady = resources[binding]
    total = sum(resources.values())
    # Between the two bounds: fully overlapped is the max, serialised is the
    # sum. `overlap` interpolates, and 1.0 recovers the makespan form.
    steady_state = steady + (1.0 - overlap) * (total - steady)
    # The slowest stage's first chunk is ALREADY included in steady.
    # Adding total/chunks used to count it twice, even for a single chunk.
    # Order is disk -> host decode -> H2D -> GEMM -> D2H.
    ordered = [resources[key] for key in ('disk', 'decode', 'h2d', 'gemm', 'd2h')
               if key in resources]
    pipeline = finite_pipeline_seconds(ordered, m, chunk_variants)
    fill_drain = overlap * (pipeline - steady)
    # PER-CHUNK HOST ISSUE COST, paid once per chunk, so it FALLS as chunks
    # grow. It is the only term here that does, and without it the model has no
    # minimum in chunk size at all -- `fill_drain` rises with chunk, so the
    # model named the smallest chunk optimal and predicted a 0.3% spread across
    # a 64x range where the measurement spreads 4.5x.
    #
    # Measured, not fitted to a wall clock. The scan's own `gpu_compute` phase
    # at M=8,931,083, N=35,365, K=128 on an idle H100 runs 11.29 s over 17,444
    # chunks (chunk 512) and 4.02 s over 2,181 (chunk 4,096). Subtracting the
    # separately measured GEMM -- 0.168 ms and 1.182 ms per chunk, from
    # `benchmarks/direct_gemm_chunk_curve.py` -- leaves 479 us and 660 us per
    # chunk, i.e. a fixed ~453 us plus ~0.05 us per variant. That fixed part
    # times the chunk-count difference is 6.91 s against a measured 7.27 s.
    #
    # It is NOT kernel launch latency, which is ~5 us for the ~14 CUDA calls a
    # chunk makes. At 17,444 iterations the host cannot issue work fast enough
    # and the GPU idles *inside* the timing window, so this is host issue cost
    # seen through the device clock.
    #
    # The GEMM was the other candidate and measurement ruled it out: its
    # achieved rate is flat in chunk (0.91-1.00x of peak from 1,024 up at width
    # 156; 0.93-1.07x throughout at width 540), so chunk-dependent GEMM
    # efficiency cannot explain a 2.8x change in `gpu_compute`.
    per_chunk = per_chunk_seconds * chunks
    fixed = (metadata_bytes / text_parse_bytes_per_second
             if text_parse_bytes_per_second else 0.0)
    scan = steady_state + fill_drain + per_chunk
    write = (m * k * output_bytes_per_test / write_bytes_per_second
             if write_bytes_per_second else 0.0)
    # **The write does NOT overlap the scan.** This used to be `max(scan,
    # write)`, on the assumption that the writer runs alongside. Measured, it
    # does not: at K=512 on the full genome, enabling binary sumstats takes the
    # scan from 14.90 s to 22.42 s -- **+7.52 s**, where `max` predicts +0
    # because the 7.55 s write fits entirely inside a 14.90 s scan. And
    # 14.90 + 7.52 = 22.42 to the hundredth.
    #
    # `write_overlap` interpolates for an implementation where it does overlap:
    # 0.0 is serial (measured, and the default), 1.0 recovers the old `max`.
    write_cost = write - write_overlap * min(scan, write)
    return {
        'resource_seconds': resources,
        'binding_term': binding,
        'steady_state_seconds': steady_state,
        'fill_drain_seconds': fill_drain,
        'per_chunk_seconds_total': per_chunk,
        'fixed_seconds': fixed,
        'write_seconds': write,
        # `scan_seconds` EXCLUDES the metadata parse, because that is what the
        # instrumentation calls `open` and reports separately. Folding it in
        # made every scan prediction look 3x too slow at small M for a purely
        # bookkeeping reason.
        'scan_seconds': scan,
        'timing_device_count': 1,
        'runtime_prediction_validated': False,
        'unresolved_timing_terms': ['non-GEMM GPU statistics service',
            'process startup and phenotype/covariate preprocessing',
            'variable-record and finite-buffer scheduling'],
        'open_seconds': fixed,
        'end_to_end_seconds': scan + write_cost + fixed,
        'write_overlap_assumed': write_overlap,
        'chunks': chunks,
        'overlap_assumed': overlap,
        # The model is validated against measured scans only up to
        # `VALIDATED_MAX_CHUNK_VARIANTS`; above it a residual appears with no
        # term here, reaching +20.84 s at chunk 65,536. A calculator meant to
        # guide decisions has to say where it stops knowing, rather than
        # returning a confident number outside the range it was checked in.
        'chunk_within_validated_range':
            chunk_variants <= VALIDATED_MAX_CHUNK_VARIANTS,
        'frame_decode_amplification': amplification,
        'caveats': [c for c in (
            # Below one frame there is no aligned chunk to choose -- every size
            # straddles, and the amplification is ~2x whatever you pick. That
            # is not a planning failure, it is the STORE being built with a
            # frame this device cannot afford, and the fix is on the encoder
            # side. Worth saying explicitly because the alternative is a user
            # tuning `chunk_size` forever against a floor they cannot move.
            ('chunk %d is smaller than the store frame of %d, so every read '
             'straddles and decode does ~%.1fx the necessary work; no chunk '
             'size fixes this -- re-encode the store with frame_variants <= %d'
             % (chunk_variants, frame_variants, amplification, chunk_variants))
            if frame_variants and 0 < chunk_variants < frame_variants else None,
            ('chunk %d straddles the store frame of %d; decode does %.2fx the '
             'necessary work. The nearest aligned chunk is %d'
             % (chunk_variants, frame_variants, amplification,
                (chunk_variants // frame_variants) * frame_variants))
            if (frame_variants and chunk_variants >= frame_variants
                and chunk_variants % frame_variants) else None,
            ('chunk %d exceeds the validated ceiling of %d; the measured '
             'residual above it is unmodelled and reached +20.8 s at 65,536'
             % (chunk_variants, VALIDATED_MAX_CHUNK_VARIANTS))
            if chunk_variants > VALIDATED_MAX_CHUNK_VARIANTS else None,
            # RETRACTED: this used to warn that the 1.0 default was "mildly
            # optimistic" against a measured 0.93. The opposite is true. 0.93
            # was obtained by INVERTING this model against measured wall
            # clocks -- a fit -- and re-checked against both regimes it is
            # worse than the structural 1.0 at every point tested: on the zstd
            # store 0.99-1.04x (1.0) against 1.10-1.15x (0.93), and on the raw
            # .bed K ladder it was fitted on, 1.00/0.98/1.13/1.81 against
            # 1.02/1.01/1.21/1.88. A fitted parameter that predicts worse than
            # the unfitted structural value is not describing the machine; it
            # was compensating for the decode rate that has since been
            # measured. Passing 0.93 is what now needs the warning.
            ('overlap=%.2f was FITTED by inverting this model against measured '
             'scans and predicts worse than the structural 1.0 on both regimes '
             'checked; prefer 1.0 unless you have re-derived it' % overlap)
            if 0.0 < overlap < 1.0 else None,
        ) if c],
    }

def _align_down(value: int, multiple: int) -> int:
    """Largest multiple of `multiple` not exceeding `value`, or 0 if none.

    Down, never up: every caller here is working against a memory ceiling that
    the unaligned value already satisfies, and rounding up can breach it.
    """
    if multiple <= 0:
        return value
    return (value // multiple) * multiple


def auto_chunk_variants(n_samples: int, n_traits: int, covariate_rank: int,
                        transfer_bytes_per_variant: float,
                        device_memory_bytes: int, depth: int,
                        headroom: float = 0.85,
                        decode_on_gpu: bool = False,
                        frame_variants: int | None = None) -> int:
    """Largest chunk whose rings fit in the given device memory.

    Chunk size cannot be a constant. Every per-chunk buffer scales with the
    sample count, so a chunk tuned at 22,250 samples asks for a hundred times
    the memory at 1,000,000. A fixed 5000-variant chunk decoded to float32 is
    445 MB at N=22,250 and 20 GB at N=1,000,000, times the ring depth.

    Variants are the axis to cut: every sample of a variant is needed for its
    regression, so samples cannot be split, and cutting traits only shrinks
    the result and projection terms, which stay small until K is very large.

    Solved by bisection against `device_ring_bytes` rather than by dividing a
    fixed fraction: the ring is not proportional to the chunk (the design term
    does not scale with it), and inverting the real formula is what makes this
    the same model the planner reports. `headroom` is the only free parameter
    and it covers allocator fragmentation and kernel workspace -- not
    unmodelled buffers, which is what the old 0.25 was standing in for.

    **Memory is a ceiling here, not a target.** "The largest chunk that fits"
    is a feasibility answer, and feasibility is not optimality.

    This docstring used to justify the cap with "scan time is monotonically
    worse as the chunk grows". That is FALSE and is retracted. Measured on an
    idle H100 at M=8,931,083, N=35,365, K=128, four round-robin rounds:

        chunk    1024    1536    2048    3072    4096    8192   16384
        median  11.61   21.71    4.83    6.80    4.68    4.63    4.75

    Larger is not worse -- it is FLAT from 2,048 up (4.63-4.83 s over an 8x
    range, against a within-chunk spread of up to 1.38x, so the differences up
    there are not even resolved). Smaller is dramatically worse. The cap is
    retained for the separate reason recorded at `MAX_AUTO_CHUNK_VARIANTS`
    (unmodelled residual above it), not for a monotonicity that does not exist.

    `frame_variants`, when the source reads from a framed store, is the far
    bigger effect and the reason this function cannot simply return a
    bisection result. A store frame is the unit of decompression, so a chunk
    that straddles frames pays for every frame it touches: at frame 2,048 the
    unaligned chunks above cost 2.5x (1,024), 4.6x (1,536) and 1.5x (3,072)
    while every multiple of 2,048 lands within noise of the best. The
    bisection returns arbitrary integers -- 1,537 is as reachable as 2,048 --
    so without this the chunk lands in the penalised class by luck of the
    memory arithmetic. Aligned DOWN, since the unaligned value is the one
    already known to fit.
    """
    if min(n_samples, n_traits, depth) <= 0 or device_memory_bytes <= 0:
        raise ValueError("positive samples, traits, depth and memory required")
    budget = device_memory_bytes * headroom

    def fits(chunk: int) -> bool:
        return device_ring_bytes(
            chunk_variants=chunk, depth=depth, n_samples=n_samples,
            n_traits=n_traits, covariate_rank=covariate_rank,
            transfer_bytes_per_variant=transfer_bytes_per_variant,
            decode_on_gpu=decode_on_gpu) <= budget

    def snap(chunk: int) -> int:
        """Align to the store's frame, or leave it alone if that is worse.

        Below one whole frame there is nothing to align to -- a chunk smaller
        than a frame straddles regardless -- so the unaligned value stands and
        the memory ceiling wins. Alignment is an optimisation; feasibility is
        not negotiable.
        """
        if not frame_variants or frame_variants <= 0:
            return chunk
        aligned = _align_down(chunk, frame_variants)
        return aligned if aligned >= frame_variants else chunk

    if not fits(MIN_AUTO_CHUNK_VARIANTS):
        return MIN_AUTO_CHUNK_VARIANTS
    if fits(MAX_AUTO_CHUNK_VARIANTS):
        return snap(MAX_AUTO_CHUNK_VARIANTS)
    low, high = MIN_AUTO_CHUNK_VARIANTS, MAX_AUTO_CHUNK_VARIANTS
    while low < high:
        middle = (low + high + 1) // 2
        if fits(middle):
            low = middle
        else:
            high = middle - 1
    return max(MIN_AUTO_CHUNK_VARIANTS, snap(low))


def choose_chunk_variants(*, n_samples: int, n_traits: int,
                          covariate_rank: int,
                          transfer_bytes_per_variant: float,
                          device_memory_bytes: int, depth: int,
                          time_model_rates: dict,
                          variants: int, stored_bytes: float,
                          headroom: float = 0.85,
                          decode_on_gpu: bool = False,
                          frame_variants: int | None = None,
                          candidates=None) -> dict:
    """Pick a chunk by PREDICTED TIME among the sizes that fit.

    `auto_chunk_variants` answers a different question -- "the largest chunk
    whose rings fit" -- and never consults the time model at all. That was a
    real gap: this module can rank chunks (the frame-straddling term puts every
    aligned size ahead of every straddling one, which is the measurement's own
    verdict) and nothing asked it to.

    The split of responsibilities is deliberate:

      * MEMORY IS A CONSTRAINT, not a score. A chunk that does not fit is not a
        slow chunk, it is an impossible one, so `device_ring_bytes` filters the
        candidates and never contributes to the ranking.
      * TIME IS THE SCORE, over what survives.

    Validated against the measured curve on an idle H100 (M=8,931,083,
    N=35,365, K=128, store frame 2,048, four round-robin rounds):

        chunk    1024    1536    2048    3072    4096    8192   16384
        median  11.61   21.71    4.83    6.80    4.68    4.63    4.75

    The measured region above the frame is FLAT -- 4.63 to 4.83 s across an 8x
    span, against a within-chunk spread of up to 1.38x, so the differences up
    there are not resolved and any of them is a correct answer. What a chooser
    has to get right is avoiding 1,024 and 1,536, and the model does: it ranks
    every aligned size ahead of every straddling one.

    Returns the choice AND the ranking, because a planner that cannot say why
    it chose is not much better than a constant.
    """
    if candidates is None:
        # Powers of two up to the cap, plus the frame and its multiples, which
        # are the sizes that can actually be optimal. Sweeping every integer
        # would cost nothing in accuracy and a great deal in time.
        span = {1 << k for k in range(7, 15)}
        if frame_variants:
            span |= {frame_variants * m for m in range(1, 9)}
        candidates = sorted(c for c in span
                            if MIN_AUTO_CHUNK_VARIANTS <= c <= MAX_AUTO_CHUNK_VARIANTS)
    budget = device_memory_bytes * headroom

    feasible = []
    for chunk in candidates:
        rings = device_ring_bytes(
            chunk_variants=chunk, depth=depth, n_samples=n_samples,
            n_traits=n_traits, covariate_rank=covariate_rank,
            transfer_bytes_per_variant=transfer_bytes_per_variant,
            decode_on_gpu=decode_on_gpu)
        if rings <= budget:
            feasible.append(chunk)

    if not feasible:
        # Nothing fits the scored set; fall back to the memory answer, which
        # bisects below the smallest candidate rather than giving up.
        return {'chunk_variants': auto_chunk_variants(
                    n_samples=n_samples, n_traits=n_traits,
                    covariate_rank=covariate_rank,
                    transfer_bytes_per_variant=transfer_bytes_per_variant,
                    device_memory_bytes=device_memory_bytes, depth=depth,
                    headroom=headroom, decode_on_gpu=decode_on_gpu,
                    frame_variants=frame_variants),
                'ranking': [], 'scored': False,
                'reason': 'no candidate fits; memory bisection used'}

    ranking = []
    for chunk in feasible:
        predicted = explain_time(
            variants=variants, samples=n_samples, traits=n_traits,
            covariate_rank=covariate_rank, chunk_variants=chunk,
            stored_bytes=stored_bytes,
            transfer_bytes_per_variant=transfer_bytes_per_variant,
            frame_variants=frame_variants or 0,
            **time_model_rates)
        ranking.append({
            'chunk_variants': chunk,
            'seconds': predicted['end_to_end_seconds'],
            'binding': predicted['binding_term'],
            'frame_decode_amplification':
                predicted.get('frame_decode_amplification', 1.0),
        })
    ranking.sort(key=lambda row: row['seconds'])
    return {'chunk_variants': ranking[0]['chunk_variants'],
            'ranking': ranking, 'scored': True,
            'reason': 'lowest predicted end-to-end among the chunks that fit'}


MIN_AUTO_TRAIT_BLOCK = 64


def auto_trait_block(n_samples: int, n_traits: int, covariate_rank: int,
                     chunk_variants: int, depth: int,
                     transfer_bytes_per_variant: float,
                     device_memory_bytes: int, headroom: float = 0.85,
                     decode_on_gpu: bool = False,
                     host_memory_bytes: int | None = None,
                     reduction_width: int | None = None,
                     trait_devices: int = 1) -> int:
    """Largest trait block whose rings fit, or all of them if they already do.

    The same inversion as `auto_chunk_variants`, on the other axis. At voxel
    scale K is the wall rather than the variant count: the design matrix alone
    is `n_samples * K * 4`, which at 33,417 subjects and 2,085,000 voxels is
    279 GB -- no device holds it, whatever the chunk size. Blocking the traits
    cuts that term, and `device_ring_bytes` already carries it as `n * (k + c)
    * 4`, so the block width is found by bisecting the same formula.

    **The device is not the only ceiling, and at high K it is not the binding
    one.** Every ring slot is pinned on the host as well, and the result ring
    carries `chunk * block` cells against the device ring's `n * block`. At
    33,417 samples, chunk 4096 and depth 32 the host ring is 6.8x the device
    ring, so bisecting the device alone chose a 471,971-trait block that needs
    1,984 GB of pinned host memory across four shards -- twice what the host
    has. Pinned pages cannot be swapped, so that is not a slow run, it is a
    failed allocation. Pass `host_memory_bytes` and the block satisfies both.

    `reduction_width` is what collapses the host side: a reduced scan stages
    `chunk * width` instead of `chunk * K`, which at width 100 turns those
    1,984 GB into 5 GB. That is why the stress test is feasible reduced and
    not feasible unreduced, and the calculator now says so instead of finding
    out at allocation time.

    Returns `n_traits` when the whole matrix fits, so a small K pays nothing.
    """
    if min(n_samples, n_traits, chunk_variants, depth) <= 0:
        raise ValueError("positive samples, traits, chunk and depth required")
    budget = device_memory_bytes * headroom
    host_budget = (None if host_memory_bytes is None
                   else host_memory_bytes * headroom)
    shards = max(int(trait_devices), 1)

    def fits(block: int) -> bool:
        if device_ring_bytes(
                chunk_variants=chunk_variants, depth=depth, n_samples=n_samples,
                n_traits=block, covariate_rank=covariate_rank,
                transfer_bytes_per_variant=transfer_bytes_per_variant,
                decode_on_gpu=decode_on_gpu) > budget:
            return False
        if host_budget is None:
            return True
        # Each shard pins its own staging and result rings in the one host
        # address space, so the host cost multiplies by the device count
        # exactly where the device cost divides by it.
        return host_pinned_bytes(
            chunk_variants=chunk_variants, depth=depth, n_traits=block,
            transfer_bytes_per_variant=transfer_bytes_per_variant,
            reduction_width=reduction_width,
            stage_on_host=not decode_on_gpu) * shards <= host_budget

    if fits(n_traits):
        return int(n_traits)
    if not fits(MIN_AUTO_TRAIT_BLOCK):
        return MIN_AUTO_TRAIT_BLOCK
    low, high = MIN_AUTO_TRAIT_BLOCK, int(n_traits)
    while low < high:
        middle = (low + high + 1) // 2
        if fits(middle):
            low = middle
        else:
            high = middle - 1
    return max(MIN_AUTO_TRAIT_BLOCK, low)


@dataclass(frozen=True)
class Workload:
    variants: int
    samples: int
    traits: int
    covariates: int  # rank excluding the intercept; design adds one intercept
    output_bytes_per_test: float = 0.0
    statistics_kernel: str = 'torch'


@dataclass(frozen=True)
class InputProfile:
    # Labels are descriptive: GPU labels do not assert an implemented decoder.
    format: str
    stored_bytes: int
    transfer_bytes_per_variant: float
    decoded_bytes_per_value: float
    cpu_decode_core_seconds_per_variant: float = 0.0
    gpu_decode_seconds_per_variant: float = 0.0
    # Additional decode traffic beyond stored read, decoded write and DMA.
    extra_host_bytes_per_variant: float = 0.0
    max_stored_bytes_per_variant: int = 0
    decode_on_gpu: bool = False
    decode_tile_multiple_of_chunk: bool = False
    direct_native_fill: bool = False
    decode_tile_matches_chunk: bool = False
    # Decoded/compressed host buffer -> pinned buffer.
    #
    # **This default is a trap for hand-written profiles.** A reader that fills
    # pinned buffers directly makes NO staging copy, so `direct_native_fill`
    # usually wants 0 here, and leaving the default charges host memory for
    # traffic that does not happen. Measured while predicting the 2026-09-13
    # format matrix: with the default left in place the model returned 137%,
    # 111%, 200% and 171% of measured wall time -- a documented *lower bound*
    # sitting above the measurement for four of five formats. Set to 0 for the
    # direct-fill readers it returns 66-105% and behaves as a bound again.
    #
    # It is deliberately NOT implied by `direct_native_fill` and not a hard
    # error: a direct-fill source still stages once when the sample order is
    # not the identity (see `native_host_staging_copies`). Profiles derived
    # from a live source have this cross-checked against the source's actual
    # staging path; hand-written profiles do not, so state it explicitly.
    host_staging_copies: int = 1
    native_encoding: str = 'dosage'
    native_row_width: int = 0  # physical transfer elements, including padding
    native_input_bytes_per_value: float = 0.0
    gpu_input_conversion_seconds_per_variant: float = 0.0  # zero is an explicit omission/absence assumption
    # Direct service of the complete resident statistics kernel (prep/GEMM/finish),
    # independently measured. Avoid disguising non-GEMM work as an effective FLOP rate.
    gpu_statistics_seconds_per_variant: float | None = None
    gpu_statistics_component_shape: tuple | None = None  # (samples, traits, covariates, chunk)
    gpu_statistics_component_source: str | None = None
    gpu_statistics_component_kernel: str | None = None
    calibration_context: dict | None = None

    # Optional exact-plan work aggregation, including LD-prefix replays.
    cpu_decode_core_seconds_total: float | None = None
    genotype_record_read_bytes_total: int | None = None
    decode_work_shape: tuple | None = None  # M,N,read,decode,chunk
    decode_work_source: str | None = None


@dataclass(frozen=True)
class Hardware:
    disk_bytes_per_second: float
    h2d_bytes_per_second: float
    host_bytes_per_second: float
    gpu_flops_per_second: float
    gpu_bytes_per_second: float
    host_memory_bytes: int
    device_memory_bytes: int
    cpu_workers: int
    d2h_bytes_per_second: float
    # Zero means output disabled; output-enabled workloads require a rate.
    output_bytes_per_second: float = 0.0
    read_latency_seconds: float = 0.0
    decode_launch_seconds: float = 0.0
    compute_launch_seconds: float = 0.0
    device_count: int = 1  # timing supports one GPU only
    calibration_context: dict | None = None
    strict_calibration: bool = False


@dataclass(frozen=True)
class ExecutionProfile:
    """Independent execution components, not a fit to association elapsed.

    Setup excludes variant-metadata parsing, which scales with file rows.
    Consumer service includes H2D, statistics, D2H and owned result delivery;
    it replaces those resource timings rather than being added to them.
    All times are wall service under the supplied placement/runtime context.
    """
    setup_seconds: float
    metadata_seconds_per_variant: float
    consumer_seconds_per_chunk: float
    shape: tuple  # (samples, traits, covariate rank, chunk, depth, workers)
    source: str
    calibration_context: dict
    statistics_kernel: str = 'torch'
    output_mode: str = 'none'
    setup_low_seconds: float | None = None
    setup_high_seconds: float | None = None
    producer_seconds_per_chunk: float | None = None
    first_use_seconds: float = 0.0


@dataclass(frozen=True)
class PipelinePlan:
    read_variants: int
    decode_variants: int
    chunk_variants: int
    workers: int
    depth: int


def _validate(w, p, h, plan):
    from .model_provenance import check_calibration_context
    check_calibration_context(p.calibration_context, h.calibration_context,
                              strict=h.strict_calibration)
    for obj in (w, p, h, plan):
        for key, value in asdict(obj).items():
            if isinstance(value, (int, float)) and (not math.isfinite(value) or value < 0):
                raise ValueError(f'{key} must be finite and nonnegative')
    for key in ('variants', 'samples', 'traits'):
        if getattr(w, key) <= 0:
            raise ValueError(f'{key} must be positive')
    if h.device_count != 1:
        raise ValueError('the current timing calculator supports exactly one GPU')
    if p.gpu_statistics_seconds_per_variant is not None:
        expected = (w.samples, w.traits, w.covariates, plan.chunk_variants)
        if tuple(p.gpu_statistics_component_shape or ()) != expected:
            raise ValueError('resident statistics component shape does not match the workload/plan')
        if not p.gpu_statistics_component_source:
            raise ValueError('resident statistics component provenance is required')
        if p.gpu_statistics_component_kernel != w.statistics_kernel:
            raise ValueError('resident statistics component kernel does not match execution')
        if h.compute_launch_seconds:
            raise ValueError('resident statistics service already includes submission overhead')
    if p.native_encoding == 'pgen_2bit':
        expected_row = ((w.samples + 3) // 4 + 63) // 64 * 64
        if p.native_row_width != expected_row or p.native_input_bytes_per_value != 1:
            raise ValueError('packed PGEN requires 64-byte-padded uint8 physical rows')
        if p.transfer_bytes_per_variant != expected_row or p.decode_on_gpu:
            raise ValueError('packed PGEN transfer width or decode placement mismatch')
        if w.statistics_kernel != 'native_fused':
            raise ValueError('packed PGEN requires native fused statistics')
    if w.statistics_kernel not in {'torch', 'native_fused'}:
        raise ValueError('unknown statistics_kernel')
    if w.statistics_kernel == 'native_fused':
        if p.native_input_bytes_per_value <= 0:
            raise ValueError('native fused statistics require positive native input width')
        if p.gpu_input_conversion_seconds_per_variant:
            raise ValueError('native fused statistics already include input conversion')
    if w.covariates + 2 >= w.samples:
        raise ValueError('insufficient residual degrees of freedom')
    for key in ('disk_bytes_per_second', 'h2d_bytes_per_second', 'host_bytes_per_second',
                'gpu_flops_per_second', 'gpu_bytes_per_second', 'd2h_bytes_per_second', 'cpu_workers'):
        if getattr(h, key) <= 0:
            raise ValueError(f'{key} must be positive')
    if p.stored_bytes <= 0 or p.transfer_bytes_per_variant <= 0 or p.decoded_bytes_per_value <= 0:
        raise ValueError('input byte sizes must be positive')
    if not isinstance(p.host_staging_copies, int) or isinstance(p.host_staging_copies, bool) or p.host_staging_copies < 0:
        raise ValueError('host_staging_copies must be a nonnegative integer')
    if min(asdict(plan).values()) <= 0 or plan.workers > h.cpu_workers:
        raise ValueError('invalid pipeline dimensions or worker count')
    if plan.depth < 2:
        raise ValueError('overlapped pipeline requires depth >= 2')
    if w.output_bytes_per_test and h.output_bytes_per_second <= 0:
        raise ValueError('output bandwidth required when output is enabled')


def estimate(w: Workload, p: InputProfile, h: Hardware, plan: PipelinePlan,
             *, execution: ExecutionProfile | None = None, operating_conditions=None) -> dict:
    """Return a resource lower bound and an explicit finite-pipeline estimate.

    GPU decode and statistical kernels are conservatively serialized on the
    same single GPU. CPU decode and DMA share host bandwidth. Output is
    assumed to use the same storage as input. Multi-GPU timing is unsupported.
    """
    _validate(w, p, h, plan)
    if execution is not None and operating_conditions is not None:
        raise ValueError('execution wall-service already includes load; do not apply operating conditions twice')
    if w.output_bytes_per_test and not h.output_bytes_per_second:
        raise ValueError(
            "workload writes output_bytes_per_test but hardware gives no "
            "output_bytes_per_second; a missing write rate silently prices "
            "serialization at zero"
        )
    m, n, k, c = w.variants, w.samples, w.traits, w.covariates + 1
    decode_tile = plan.decode_variants
    read_tile = plan.read_variants
    if p.decode_on_gpu and p.decode_tile_multiple_of_chunk:
        decode_tile = math.ceil(max(plan.chunk_variants, decode_tile) / plan.chunk_variants) * plan.chunk_variants
        read_tile = decode_tile  # BGEN reads precisely one rounded decode batch
    if p.decode_tile_matches_chunk:
        read_tile = decode_tile = plan.chunk_variants
    has_decode_work = p.cpu_decode_core_seconds_total is not None or p.genotype_record_read_bytes_total is not None
    if has_decode_work:
        if tuple(p.decode_work_shape or ()) != (m,n,read_tile,decode_tile,plan.chunk_variants):
            raise ValueError('decode work statistics do not match workload and plan')
        if not p.decode_work_source:
            raise ValueError('decode work statistics require a source')
        if execution is not None:
            raise ValueError('execution producer wall-service already includes input-specific decoding')
    reads = math.ceil(m / read_tile)
    tiles = math.ceil(m / decode_tile)
    from .model_provenance import check_calibration_context
    chunks = math.ceil(m / plan.chunk_variants)
    transfer = m * p.transfer_bytes_per_variant
    decoded_row_bytes = (p.native_row_width if p.native_encoding == 'pgen_2bit'
                         else n * p.decoded_bytes_per_value)
    decoded = m * decoded_row_bytes
    output = m * k * w.output_bytes_per_test
    # Dense output also returns float64 -log10(P) before the writer casts it
    # to float32. Scan-only paths retain the older beta+t result contract.
    result_bytes_per_test = 16 if w.output_bytes_per_test else 8
    result_bytes = m * (result_bytes_per_test * k + 5)
    record_read_bytes = p.stored_bytes if p.genotype_record_read_bytes_total is None else p.genotype_record_read_bytes_total
    storage = record_read_bytes / h.disk_bytes_per_second
    output_time = output / h.output_bytes_per_second if output else 0.0
    read_service = reads * h.read_latency_seconds
    cpu_work = m * p.cpu_decode_core_seconds_per_variant if p.cpu_decode_core_seconds_total is None else p.cpu_decode_core_seconds_total
    cpu_capacity = min(plan.workers, plan.depth)
    if operating_conditions is not None:
        cpu_capacity = operating_conditions.capacity(cpu_capacity)
    cpu = (cpu_work + (0 if p.decode_on_gpu else tiles * h.decode_launch_seconds)) / cpu_capacity
    # fp32 genotype materialization, G@B, sum-of-squares and statistic arrays.
    # This is linear GWAS only: no JAGWAS, exact p-values or clumping.
    flops = m * (2 * n * (k + c) + 4 * n + 2 * c + 12 * k)
    genotype_bytes = (2 * p.native_input_bytes_per_value + 8
                      if w.statistics_kernel == 'native_fused' else 32)
    genotype_row_bytes = n * genotype_bytes
    if p.native_encoding == 'pgen_2bit':
        genotype_row_bytes = 2 * p.native_row_width + 8 * n
    gpu_bytes = m * (genotype_row_bytes + 16 * (k + c))
    roofline_compute = max(flops / h.gpu_flops_per_second, gpu_bytes / h.gpu_bytes_per_second)
    compute = (m * p.gpu_statistics_seconds_per_variant
               if p.gpu_statistics_seconds_per_variant is not None else roofline_compute)
    compute += chunks * h.compute_launch_seconds
    gpu_conversion = m * p.gpu_input_conversion_seconds_per_variant
    gpu_decode = m * p.gpu_decode_seconds_per_variant
    if p.decode_on_gpu:
        gpu_decode += tiles * h.decode_launch_seconds
    host_traffic = record_read_bytes + (1 + 2 * p.host_staging_copies) * transfer + 3 * result_bytes + m * p.extra_host_bytes_per_variant
    if not p.decode_on_gpu:
        host_traffic += decoded
    resource = {
        'storage': storage + output_time + read_service,
        'cpu_decode': cpu,
        'host_memory': host_traffic / h.host_bytes_per_second,
        'h2d': transfer / h.h2d_bytes_per_second,
        'd2h': result_bytes / h.d2h_bytes_per_second,
        'gpu': compute + gpu_decode + gpu_conversion,
    }
    resource_detail = {
        'input_read_seconds': storage,
        'input_record_read_bytes': record_read_bytes,
        'cpu_decode_work_core_seconds': cpu_work,
        'cpu_decode_capacity_cores': cpu_capacity,
        'operating_conditions': asdict(operating_conditions) if operating_conditions is not None else None,
        'output_write_seconds': output_time,
        'read_service_seconds': read_service,
        'output_bytes': output,
    }
    bottleneck = max(resource, key=resource.get)
    lower = resource[bottleneck]
    # Explicit coarse batch model: serial startup/drain of one largest batch.
    batch_fraction = min(m, max(read_tile, decode_tile, plan.chunk_variants)) / m
    fill_drain = batch_fraction * (sum(resource.values()) - lower)
    average_stored = p.stored_bytes / m
    max_stored = p.max_stored_bytes_per_variant or average_stored
    if max_stored < average_stored:
        raise ValueError('maximum record size smaller than average stored size')
    host_ring = plan.depth * ((0 if p.direct_native_fill else read_tile * max_stored) +
                              plan.chunk_variants * p.transfer_bytes_per_variant)
    if not p.decode_on_gpu:
        if not p.direct_native_fill:
            host_ring += plan.depth * plan.chunk_variants * n * 4
        if not p.direct_native_fill or p.host_staging_copies:
            host_ring += min(plan.workers, plan.depth) * decode_tile * decoded_row_bytes
    host_ring += plan.depth * plan.chunk_variants * (
        result_bytes_per_test * k + 5)
    # The same accounting the scan uses to choose its chunk, so a predicted
    # plan and a live scan cannot disagree about what fits.
    device_ring = device_ring_bytes(
        chunk_variants=plan.chunk_variants, depth=plan.depth, n_samples=n,
        n_traits=k, covariate_rank=c,
        transfer_bytes_per_variant=p.transfer_bytes_per_variant,
        decode_tile=decode_tile, decode_on_gpu=p.decode_on_gpu)
    feasible = host_ring <= h.host_memory_bytes and device_ring <= h.device_memory_bytes
    execution_result = {}
    if execution is not None:
        if not feasible:
            raise ValueError('execution plan exceeds supplied host/device memory budgets')
        if tuple(execution.shape) != (n, k, w.covariates, plan.chunk_variants, plan.depth, plan.workers):
            raise ValueError('execution components do not match workload and queue geometry')
        check_calibration_context(execution.calibration_context, h.calibration_context, strict=True)
        if not execution.source or execution.statistics_kernel != w.statistics_kernel:
            raise ValueError('execution component source/kernel mismatch')
        if execution.output_mode != 'none' or output:
            raise ValueError('execution components currently cover discarded results only')
        for key in ('setup_seconds', 'metadata_seconds_per_variant', 'consumer_seconds_per_chunk', 'first_use_seconds'):
            value = getattr(execution, key)
            if not math.isfinite(value) or value < 0:
                raise ValueError(key + ' must be finite and nonnegative')
        if execution.consumer_seconds_per_chunk <= 0:
            raise ValueError('consumer service must be positive')
        # No GPU/PCIe/result-copy repricing: these operations are already in
        # consumer service. Source read/decode and the consumer overlap.
        producer = max(storage + read_service, cpu)
        if execution.producer_seconds_per_chunk is not None:
            rate = execution.producer_seconds_per_chunk
            if not math.isfinite(rate) or rate <= 0:
                raise ValueError('producer service must be finite and positive')
            # Measured producer includes native decode and pinned staging;
            # storage is a capacity constraint, not another serial charge.
            producer = max(storage + read_service, m / plan.chunk_variants * rate)
        consumer = m / plan.chunk_variants * execution.consumer_seconds_per_chunk
        scan = finite_pipeline_seconds([producer, consumer], m, plan.chunk_variants)
        # Shared host traffic is a capacity constraint, not another serial stage.
        scan = max(scan, resource['host_memory']) + execution.first_use_seconds
        fixed = execution.setup_seconds + m * execution.metadata_seconds_per_variant
        execution_result = dict(prediction_seconds=fixed + scan,
            execution_setup_seconds=fixed, execution_scan_seconds=scan,
            execution_producer_seconds=producer, execution_consumer_seconds=consumer,
            execution_first_use_seconds=execution.first_use_seconds,
            prediction_scope='process setup + metadata + finite producer/consumer model',
            execution_component_source=execution.source,
            execution_assumptions=['homogeneous full chunks and a proportional final chunk',
                'producer and consumer component rates transfer under the recorded placement/load',
                'consumer measurement includes transfers and result delivery; no double counting',
                'metadata parsing scales with input rows; changed record widths require recalibration'],
            unresolved_timing_terms=['variable-record imbalance and shared-host scheduling changes'])
        if execution.setup_low_seconds is not None and execution.setup_high_seconds is not None:
            lo, hi = execution.setup_low_seconds, execution.setup_high_seconds
            if not (math.isfinite(lo) and math.isfinite(hi) and 0 <= lo <= execution.setup_seconds <= hi):
                raise ValueError('invalid observed setup range')
            execution_result['setup_sensitivity_seconds'] = [lo + fixed - execution.setup_seconds + scan,
                                                            hi + fixed - execution.setup_seconds + scan]
    return {
        'plan': asdict(plan), 'effective_geometry': {'read_variants': read_tile, 'decode_variants': decode_tile, 'chunk_variants': plan.chunk_variants},
        'resource_seconds': resource, 'storage_detail': resource_detail,
        'bottleneck': bottleneck,
        'io_lower_bound_seconds': storage, 'resource_lower_bound_seconds': lower,
        'planning_seconds': lower + fill_drain, 'fill_drain_seconds': fill_drain,
        'resource_service_seconds': lower,
        'is_guaranteed_lower_bound': False,
        'prediction_seconds': None,  # incomplete end-to-end model; see unresolved terms
        'io_bound_possible': lower <= storage * (1 + 1e-12),
        'host_buffer_bytes': math.ceil(host_ring), 'device_buffer_bytes': math.ceil(device_ring),
        'memory_feasible': feasible, 'memory_is_upper_bound': False,
        'memory_estimate_only': True, 'result_d2h_bytes': result_bytes,
        'gpu_compute_seconds': compute, 'gpu_decode_seconds': gpu_decode,
        'gpu_roofline_seconds_reference': roofline_compute,
        'gpu_service_basis': ('independent_resident_component'
            if p.gpu_statistics_seconds_per_variant is not None else 'operation_traffic_roofline'),
        'timing_device_count': 1,
        'runtime_prediction_validated': False,
        'unresolved_timing_terms': check_calibration_context(p.calibration_context, h.calibration_context) + ['process startup and metadata/preprocessing',
            'finite-buffer and variable-record scheduling'] +
            ([] if p.gpu_statistics_seconds_per_variant is not None else
             ['non-GEMM statistics service beyond roofline']),
        'native_encoding': p.native_encoding, 'native_row_width': p.native_row_width,
        'input_transfer_bytes': transfer,
        'statistics_kernel': w.statistics_kernel,
        'conversion_fused_into_statistics': w.statistics_kernel == 'native_fused',
        'gpu_input_conversion_seconds': gpu_conversion,
        'gpu_input_conversion_status': ('fused_into_statistics' if w.statistics_kernel == 'native_fused' else
                                        'supplied' if gpu_conversion else 'zero_assumed_absent_or_unmeasured'),
        'gpu_flops': flops, 'gpu_memory_bytes': gpu_bytes,
        'assumptions': [
            'Linear GWAS with genotype mean-centering; exact p-values, JAGWAS and clumping excluded.',
            'Sustained capacities are supplied; shape efficiency is not fitted.',
            'Same-GPU decode, input conversion and compute service demands add.',
            'Zero conversion cost is an absence/omission assumption, not evidence that conversion is free.',
            'Torch statistical genotype traffic is logical32*N; native fused uses(2*input_width+8)*N with conversion included. Small-array traffic is approximate.',
            'Effective GPU compute rates must match the selected statistics kernel and precision; they are not hardware peak rates.',
            'CPU concurrency is limited by both worker count and prefetch depth.',
            'Input and output storage demands add; variable records need measured maximum sizes.',
            'Finite-pipeline estimate assumes sufficient overlap and no unmodeled dependencies.',
            'Buffer estimate includes float32 genotype and pinned result rings but is not an allocation guarantee; reserve decoder workspace, allocator overhead and retained results separately.',
        ],
        **execution_result,
    }


def choose_plan(w: Workload, p: InputProfile, h: Hardware, *,
                read_variants=(1024, 4096, 16384), decode_variants=(256, 1024, 4096),
                chunk_variants=(256, 1024, 4096), workers=None, depths=None) -> dict:
    """Choose among explicit candidates; no empirical saturation constants.

    Worker speedup is ideal until host bandwidth binds. This is a starting plan,
    requiring end-to-end validation. Depth must accommodate active CPU workers;
    beyond that, ties choose the smallest memory footprint, then fewest workers.
    """
    if workers is None:
        workers = tuple(sorted({1 << exponent for exponent in range(max(1, int(h.cpu_workers).bit_length()))} | {h.cpu_workers}))
    if depths is None:
        depths = tuple(sorted({max(2, min(value, h.cpu_workers))
                               for value in (2, 3, 4, 8, 16, h.cpu_workers)}))
    best = None
    for values in itertools.product(read_variants, decode_variants, chunk_variants, workers, depths):
        plan = PipelinePlan(*values)
        result = estimate(w, p, h, plan)
        if not result['memory_feasible']:
            continue
        score = (result['planning_seconds'], result['host_buffer_bytes'] + result['device_buffer_bytes'], plan.workers)
        if best is None or score < best[0]:
            best = score, result
    if best is None:
        raise ValueError('no candidate fits the supplied memory budgets')
    return best[1]


def plan_parameters(w: Workload, p: InputProfile, h: Hardware, **candidates) -> dict:
    """Return common scan kwargs plus separately named adapter recommendations.

    Only ``scan_kwargs`` can be passed directly to the shared scan iterator.
    Source-specific controls must be translated by its adapter; unsupported
    read/decode batch controls must never be silently claimed as applied.
    """
    result = choose_plan(w, p, h, **candidates)
    plan = result['plan']
    return {
        'scan_kwargs': {'chunk_size': plan['chunk_variants'],
                        'reader_workers': plan['workers'],
                        'prefetch_chunks': plan['depth']},
        'source_recommendations': {'read_variants': plan['read_variants'],
                                   'decode_variants': plan['decode_variants']},
        'estimate': result,
    }


def apply_pipeline_profile(source, phenotype, covariates, profile: dict) -> dict:
    """Validate an explicit profile against this scan, then apply supported knobs.

    ``profile`` contains ``hardware`` and ``input`` dictionaries, optionally
    ``workload`` assertions/output volume, ``candidates`` and ``source_controls``.
    Explicit scan arguments must constrain the corresponding candidates BEFORE
    calling this helper: chunk_size -> chunk_variants, reader_workers -> workers,
    prefetch_chunks -> depths. Never overwrite returned kwargs while reporting
    the original estimate. No genotypes are read by this helper.
    """
    import numpy as np
    from pathlib import Path
    from .preprocess import residualize_and_standardize

    allowed = {'hardware', 'input', 'workload', 'candidates', 'source_controls'}
    if not isinstance(profile, dict) or set(profile) - allowed:
        raise ValueError('invalid pipeline profile fields')
    if 'hardware' not in profile or 'input' not in profile:
        raise ValueError('pipeline profile requires explicit hardware and input records')
    required_costs = {'cpu_decode_core_seconds_per_variant', 'gpu_decode_seconds_per_variant'}
    if not required_costs.issubset(profile['input']):
        raise ValueError('pipeline profile must explicitly supply both CPU and GPU decode costs (zero if absent)')
    h, p = Hardware(**profile['hardware']), InputProfile(**profile['input'])
    shape = tuple(source.shape)
    y = np.asarray(phenotype)
    if y.ndim != 2 or len(shape) != 2 or min(shape) <= 0 or y.shape[0] != shape[0] or y.shape[1] < 1:
        raise ValueError('phenotype must be a sample-aligned two-dimensional matrix')
    if not np.isfinite(y).all():
        raise ValueError('phenotype must be finite before pipeline planning')
    rank = 0
    if covariates is not None:
        cov = np.asarray(covariates)
        if cov.ndim != 2 or cov.shape[0] != shape[0] or not np.isfinite(cov).all():
            raise ValueError('covariates must be a finite sample-aligned matrix')
        # Use the engine's centered/scaled SVD tolerance, including its dtype.
        _, basis = residualize_and_standardize(np.zeros((shape[0], 1)), cov)
        rank = 0 if basis is None else basis.shape[1]
    supplied_workload = dict(profile.get('workload', {}))
    output_bytes = supplied_workload.pop('output_bytes_per_test', 0.0)
    from .scan_gpu import resolve_statistics_backend
    statistics_kernel = resolve_statistics_backend()
    w = Workload(int(shape[1]), int(shape[0]), int(y.shape[1]), int(rank), output_bytes, statistics_kernel)
    encoding = getattr(source, 'native_encoding', 'dosage')
    native_width = 4 if p.decode_on_gpu else np.dtype(getattr(source, 'native_transfer_dtype', getattr(source, 'native_dtype', np.float32))).itemsize
    physical_row = int(getattr(source, 'native_row_width', w.samples))
    if 'native_encoding' in profile['input'] and p.native_encoding != encoding:
        raise ValueError('profile native encoding differs from actual source')
    if p.native_row_width not in (0, physical_row):
        raise ValueError('profile physical row width differs from actual source')
    if ('native_input_bytes_per_value' in profile['input'] and p.native_input_bytes_per_value not in (0, native_width)):
        raise ValueError('profile native input width differs from actual tensor dtype')
    p = replace(p, native_input_bytes_per_value=native_width, native_encoding=encoding, native_row_width=physical_row)
    actual_workload = asdict(w)
    for name, value in supplied_workload.items():
        if name not in actual_workload or value != actual_workload[name]:
            raise ValueError(f'pipeline workload {name} does not match actual scan')

    actual_bytes = getattr(source, 'input_bytes', None)
    byte_evidence = 'source.input_bytes' if actual_bytes is not None else None
    if actual_bytes is None:
        for name in ('zst_path', 'bed_path', 'genotype_path'):
            path = getattr(source, name, None)
            if path is not None and Path(path).is_file():
                actual_bytes = Path(path).stat().st_size
                byte_evidence = name + '.stat().st_size'
                break
    if actual_bytes is not None and p.stored_bytes != int(actual_bytes):
        raise ValueError(f'profile stored_bytes={p.stored_bytes} differs from actual source bytes={actual_bytes}')
    direct_fill = bool(getattr(source, 'allows_direct_native_fill', False)) and not p.decode_on_gpu
    expected_staging = (getattr(source, 'native_host_staging_copies', 1) if direct_fill
                        else getattr(source, 'host_staging_copies', 1))
    if 'host_staging_copies' not in profile['input']:
        p = replace(p, host_staging_copies=expected_staging)
    if p.host_staging_copies != expected_staging:
        raise ValueError('profile host_staging_copies does not match source staging path')
    known_backend = getattr(source, 'backend_used', None) or getattr(source, 'decode_backend', None)
    if known_backend in {'cpu', 'gpu'} and p.decode_on_gpu != (known_backend == 'gpu'):
        raise ValueError('profile decode placement does not match resolved source backend')
    if not p.decode_on_gpu and hasattr(source, 'native_dtype'):
        expected_transfer = physical_row * native_width
        if p.transfer_bytes_per_variant != expected_transfer:
            raise ValueError('profile H2D representation does not match native source dtype/sample count')

    controls = dict(profile.get('source_controls', {}))
    if set(controls) - {'read_workers', 'read_ahead_batches'}:
        raise ValueError('source_controls supports only read_workers and read_ahead_batches')
    for name, value in controls.items():
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f'{name} must be a positive integer')
    p = replace(p, direct_native_fill=direct_fill,
                decode_tile_matches_chunk=direct_fill and bool(getattr(source, 'decode_tile_matches_chunk', False)),
                decode_tile_multiple_of_chunk=bool(getattr(source, 'decode_tile_multiple_of_chunk', False)))
    candidates = dict(profile.get('candidates', {}))
    model_constraints = {}
    parse_parameter = getattr(source, 'cpu_parse_worker_parameter', None)
    if parse_parameter is None and hasattr(source, 'read_workers'):
        parse_parameter = 'read_workers'
    if parse_parameter is not None and not hasattr(source, parse_parameter):
        raise ValueError('source declares an unavailable CPU parsing worker parameter')
    if p.decode_on_gpu:
        if parse_parameter is None:
            parse_workers = 1
            if 'workers' in candidates and any(value != 1 for value in candidates['workers']):
                raise ValueError('GPU profile workers must be 1 without a declared read_workers/CPU parsing pool')
            candidates['workers'] = [1]
        elif 'read_workers' in controls:
            parse_workers = controls['read_workers']
            if 'workers' in candidates and any(value != parse_workers for value in candidates['workers']):
                raise ValueError('GPU profile workers must match source_controls.read_workers')
            candidates['workers'] = [parse_workers]
        elif 'workers' not in candidates:
            candidates['workers'] = [int(getattr(source, parse_parameter))]
        model_constraints['cpu_parse_worker_parameter'] = parse_parameter
    result = plan_parameters(w, p, h, **candidates)
    plan = result['estimate']['plan']
    geometry = result['estimate']['effective_geometry']
    result['source_recommendations'] = {'read_variants': geometry['read_variants'], 'decode_variants': geometry['decode_variants']}
    if p.decode_on_gpu:
        model_constraints['cpu_parse_workers'] = plan['workers']
    applied, advisory = {}, {}
    proposed = {
        'decode_batch_size': geometry['decode_variants'],
        'read_batch_bytes': math.ceil(plan['read_variants'] * p.stored_bytes / w.variants),
    }
    # CPU adapters expose these controls; fixed GPU decoders currently do not
    # promise independently configurable CPU/decode worker pools.
    if p.decode_on_gpu and p.decode_tile_multiple_of_chunk:
        proposed.pop('read_batch_bytes')
        model_constraints['read_batch_control'] = 'decode_batch_size'
    if p.decode_tile_matches_chunk:
        proposed.pop('read_batch_bytes', None)
        proposed.pop('decode_batch_size', None)
        model_constraints.update(read_batch_control='chunk_size', decode_batch_control='chunk_size',
                                 host_staging_copies=p.host_staging_copies, direct_native_fill=True)
    if hasattr(source, 'prefetch_chunks'):
        proposed['prefetch_chunks'] = plan['depth']
    if not p.decode_on_gpu:
        for worker_parameter in ('decode_workers', 'reader_workers'):
            if hasattr(source, worker_parameter):
                proposed[worker_parameter] = plan['workers']
    elif parse_parameter is not None:
        proposed[parse_parameter] = plan['workers']
    else:
        advisory['decode_workers'] = {'value': plan['workers'], 'reason': 'GPU adapter worker mapping is not declared'}
    for name, value in controls.items():
        target = parse_parameter if p.decode_on_gpu and name == 'read_workers' and parse_parameter else name
        proposed[target] = value
    for name, value in proposed.items():
        if name == 'decode_batch_size' and p.decode_tile_multiple_of_chunk and not p.decode_on_gpu:
            advisory[name] = {'value': value, 'reason': 'CPU backend decodes compute chunks; independent GPU decode tile is inactive'}
        elif hasattr(source, name):
            applied[name] = value
        else:
            advisory[name] = {'value': value, 'reason': 'source does not expose this control'}
    previous = {name: getattr(source, name) for name in applied}
    try:
        for name, value in applied.items():
            setattr(source, name, value)
    except BaseException:
        for name, value in previous.items():
            setattr(source, name, value)
        raise
    result.update(workload=actual_workload, decode_input_samples=int(getattr(source, '_n_bgen_samples', w.samples)), model_constraints=model_constraints, applied_source_settings=applied,
                  advisory_source_settings=advisory, stored_bytes_evidence=byte_evidence,
                  all_source_controls_applied=not bool(advisory))
    result['estimate']['application_status'] = ('conditional_on_advisory_controls' if advisory else 'controls_applied')
    result['estimate']['gpu_input_conversion_status'] = (
        'fused_into_statistics' if w.statistics_kernel == 'native_fused' else
        'supplied' if 'gpu_input_conversion_seconds_per_variant' in profile['input']
        else 'unmeasured_default_zero')
    result['estimate']['advisory_conditions'] = list(advisory)
    result['estimate']['applied_io_bound_possible'] = (None if advisory else result['estimate']['io_bound_possible'])
    result['estimate']['assumptions'].append(
        'If source controls are advisory, timing and memory results are conditional on implementing them and do not describe the actual configured runtime.')
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profile', help='JSON with workload, input, hardware; optional candidates or a measured execution profile and exact plan')
    args = parser.parse_args()
    with open(args.profile, encoding='utf-8') as handle:
        data = json.load(handle)
    if 'execution' in data or 'operating_conditions' in data or 'plan' in data:
        if 'candidates' in data or 'plan' not in data:
            raise ValueError('execution components require an exact plan, not candidate search')
        from .work_statistics import CpuOperatingConditions
        result = estimate(Workload(**data['workload']), InputProfile(**data['input']),
                          Hardware(**data['hardware']), PipelinePlan(**data['plan']),
                          execution=ExecutionProfile(**data['execution']) if 'execution' in data else None,
                          operating_conditions=CpuOperatingConditions(**data['operating_conditions']) if 'operating_conditions' in data else None)
    else:
        result = choose_plan(Workload(**data['workload']), InputProfile(**data['input']),
                             Hardware(**data['hardware']), **data.get('candidates', {}))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
















def binary_output_pipeline_seconds(producer_chunk, consumer_chunk, variants,
                                   traits, chunk_variants, append_seconds,
                                   storage_seconds, final_seconds,
                                   block_bytes=1 << 20, queue_depth=3):
    """Finite producer/main/background-write schedule for dense binary output.

    append_seconds excludes queue stalls; storage_seconds is the critical
    service of the concurrent array writers for the whole payload.
    final_seconds contains measured fsync and metadata publication service.
    Blocks become writable only after their bytes have been produced. A bounded
    staging pool propagates storage backpressure to the consumer. Rates are
    independent component measurements, never fitted association elapsed.
    """
    if variants <= 0 or traits <= 0 or chunk_variants <= 0 or block_bytes <= 0 or queue_depth < 1:
        raise ValueError("positive geometry and queue depth required")
    if any(not math.isfinite(v) or v < 0 for v in
           (producer_chunk, consumer_chunk, append_seconds, storage_seconds, final_seconds)):
        raise ValueError("finite nonnegative component service required")
    per_array_bytes = variants * traits * 4
    main = producer = storage = 0.0
    pending = []
    filled = 0
    for start in range(0, variants, chunk_variants):
        count = min(chunk_variants, variants - start)
        producer += producer_chunk * count / chunk_variants
        main = max(main, producer) + consumer_chunk * count / chunk_variants
        main += append_seconds * count / variants
        filled += count * traits * 4
        while filled >= block_bytes:
            storage = max(storage, main) + storage_seconds * block_bytes / per_array_bytes
            pending.append(storage)
            filled -= block_bytes
            # One active staging buffer plus queue_depth queued buffers.
            if len(pending) > queue_depth:
                main = max(main, pending.pop(0))
    if filled:
        storage = max(storage, main) + storage_seconds * filled / per_array_bytes
    return max(main, storage) + final_seconds
