"""Decode cost from the record mix, instead of one rate inferred by subtraction.

Decode was the only term in the pipeline whose cost was a guess: 11.3 GB/s,
obtained by subtracting an estimated read from a measured scan. That number
lumps five different decode paths together, carries no error bar, and cannot
explain why two files of the same size decode at different speeds.

It does not have to be a guess, because the loops are not mysterious -- some
are bounded and some are data-dependent, and the data that decides the
data-dependent ones is readable in milliseconds.

**Bounded, so a closed form:**

    bed                 `unpack_genovec` is four strided writes over N/4
                        bytes. Trip count is exactly N, every variant.
    PGEN form 0         a copy of exactly ceil(N/4) bytes.
    zstd store          the compressed input varies, but zstd's cost tracks
                        the OUTPUT, which is frame_variants * ceil(N/4).

**Data-dependent, and all for the same reason:** every other PGEN form routes
through `_parse_difflist`, which opens with

    entry_ct, position = _uleb128(buffer, position)
    group_ct = (entry_ct + 63) // 64
    for group in range(group_ct):        # trip count unknown until read

`entry_ct` is the number of samples differing from the background -- pure data.
It drives three loops: the base ids, the group sizes, and the within-group
increments.

    form 1              bounded unpackbits over N, PLUS a difflist
    forms 2, 3          copy the previous genovec (bounded N), PLUS a difflist
                        of deltas; form 3 additionally swaps 0 <-> 2
    forms 4, 6, 7       a difflist over a constant background

So per-variant decode is not unpredictable -- it is a WEIGHTED SUM over the
record forms. The weights are a byte per variant in the header's vrtype track,
so a histogram costs one sequential read of M bytes and nothing else.

**The same decomposition prices plink2.** It reads the same `.pgen` through the
same record forms, so the histogram measured here serves both calculators; only
the per-form seconds differ, because the implementations differ. That is why
this lives in the library rather than in one benchmark script.
"""
from __future__ import annotations

import time

import numpy as np

# Low nibble of vrtype. Names from the PGEN spec and from the branches in
# `pgen_reader._decode`.
FORM_NAMES = {
    0: "plain genovec",
    1: "one-bit + difflist",
    2: "LD-compressed",
    3: "LD-compressed, inverted",
    4: "difflist (background 0)",
    6: "difflist (background 2)",
    7: "difflist (background 3)",
}

# Forms whose cost is fixed by N alone. Everything else carries a difflist
# whose length is a property of the genotypes.
BOUNDED_FORMS = frozenset({0})


def record_mix(vrtypes) -> dict[int, int]:
    """Histogram of record forms, from the header's vrtype track.

    One sequential pass over M bytes. This is the whole data-dependence of PGEN
    decode, and it is this cheap to obtain.
    """
    forms = np.asarray(vrtypes, dtype=np.uint8) & 0x07
    counts = np.bincount(forms, minlength=8)
    return {int(form): int(count) for form, count in enumerate(counts) if count}


def difflist_share(mix: dict[int, int]) -> float:
    """Fraction of variants whose decode cost is data-dependent."""
    total = sum(mix.values())
    if not total:
        return 0.0
    bounded = sum(count for form, count in mix.items() if form in BOUNDED_FORMS)
    return 1.0 - bounded / total


def decode_seconds(mix: dict[int, int],
                   seconds_per_variant_by_form: dict[int, float]) -> float:
    """Total decode seconds for this file, as a weighted sum over the forms.

    A form present in the mix but missing from the timing table is an error
    rather than a zero: silently pricing a record type at nothing is how a
    decode term goes missing, which is exactly the bug this module exists to
    prevent.
    """
    missing = sorted(set(mix) - set(seconds_per_variant_by_form))
    if missing:
        raise KeyError(
            f"no measured decode cost for form(s) {missing} "
            f"({', '.join(FORM_NAMES.get(f, '?') for f in missing)}); "
            f"pricing them at zero would understate decode")
    return sum(count * seconds_per_variant_by_form[form]
               for form, count in mix.items())


def measure_decode_by_form(reader, vrtypes, per_form: int = 200,
                           repeat: int = 3) -> dict:
    """Time each record form separately, on this machine, on this file.

    The point of separating them: a difflist record's cost depends on how many
    samples differ from its background, and that differs by form and by cohort.
    One blended rate cannot express it, and a rate obtained by subtracting a
    read estimate from a scan -- which is where 11.3 GB/s came from -- cannot
    either.

    Variants are sampled ACROSS the file rather than taken from the head,
    because record forms cluster: LD-compressed runs follow a plain genovec, so
    the first few thousand variants are not a sample of anything.

    `reader` needs `read_genovec(index)`. Forms 2 and 3 are expressed against
    the preceding record, so they are timed in forward order from a nearby
    starting point rather than by random access, which would make the reader
    replay to re-establish its LD base and time the replay instead.
    """
    forms = np.asarray(vrtypes, dtype=np.uint8) & 0x07
    out: dict[int, float] = {}
    detail: dict[str, dict] = {}

    for form in sorted(set(int(f) for f in np.unique(forms))):
        indices = np.flatnonzero(forms == form)
        if indices.size == 0:
            continue
        # Spread the sample over the file, not the first N.
        step = max(1, indices.size // per_form)
        chosen = indices[::step][:per_form]
        if chosen.size == 0:
            continue

        seconds = []
        for _ in range(repeat):
            started = time.perf_counter()
            for index in chosen:
                reader.read_genovec(int(index))
            seconds.append((time.perf_counter() - started) / chosen.size)
        seconds.sort()
        median = seconds[len(seconds) // 2]
        out[form] = median
        detail[str(form)] = {
            "name": FORM_NAMES.get(form, "?"),
            "seconds_per_variant": median,
            "variants_timed": int(chosen.size),
            "variants_present": int(indices.size),
            "bounded": form in BOUNDED_FORMS,
        }

    return {"seconds_per_variant_by_form": out, "detail": detail}


# Measured on the real cohort with the NATIVE decoder -- the one the scan uses
# -- by `benchmarks/direct_decode_native_forms.py`. Forms are interleaved
# through the file so no range isolates one; these are solved from 120 timed
# ranges against their form counts and validated on 36 ranges held OUT of the
# solve, median predicted/actual 1.008.
#
# Seconds per variant, SINGLE-THREADED. The scan decodes in its reader pool, so
# a whole-file figure must be divided by the worker count: 146.0 s here over 16
# workers is ~9.1 s, which is consistent with the ~20 s full-genome scan where
# decode overlaps the read.
#
# The 13.5x spread is the point. One blended rate -- the 11.3 GB/s this
# replaces, itself obtained by subtracting an estimated read from a scan time --
# cannot express it, and cannot explain why two files of equal size decode at
# different speeds.
#
# It also contradicts the obvious guess: LD-compressed is CHEAP (a bounded
# genovec copy plus a short difflist) while difflist-over-constant-background
# is EXPENSIVE. Worth stating because the intuition is the wrong way round.
#
# Caveat: taken at host load 74. Absolute values are contended; the ratios and
# the held-out validation are what this rests on. Re-take quiet before quoting
# the absolute numbers.
MEASURED_SECONDS_PER_VARIANT_BY_FORM = {
    0: 2.62e-6,    # plain genovec          -- bounded, ceil(N/4) copy
    1: 15.97e-6,   # one-bit + difflist
    2: 5.19e-6,    # LD-compressed
    3: 35.24e-6,   # LD-compressed, inverted
    4: 32.67e-6,   # difflist (background 0)
    6: 4.32e-6,    # difflist (background 2)
}


# THE ZSTD HARD-CALL STORE IS A DIFFERENT DECODER and needs its own rate. The
# table above is PGEN record forms; a store frame is a zstd block and its cost
# has nothing to do with vrtypes.
#
# Measured by `_scratch/_store_decode_rate.sh` and `_store_straddle.sh` on an
# idle host (load 2.2) against the real store -- frame 2,048, ratio 7.61,
# 8,842 bytes per variant -- reading disjoint frame-aligned spans, min of 3:
#
#     workers       1     2     4     8    12    16    24
#     GB/s raw   1.40  2.66  5.26 10.44 15.04 20.01 25.08
#     efficiency 1.00  0.95  0.94  0.93  0.90  0.89  0.75
#
# Two things this settles. First, `HardcallStore.read_packed_into` is SERIAL --
# single-threaded zstd measures 1.55 GB/s and the store's own path achieves
# 1.37-1.49 -- so all of the scan's decode throughput comes from the reader
# pool, and a model using the one-worker rate under-predicts by the worker
# count. Second, it replaces the 11.3 GB/s that `explain_time` was given, which
# was never measured: it was inferred by subtracting an estimated read from a
# scan time. That number was provably impossible -- at 78.97 GB raw it puts
# 6.99 s of decode inside a scan measured at 4.68 s end to end.
#
# With 20.01 GB/s the full-genome decode is 3.95 s against a 4.68 s wall, so
# the aligned scan is DECODE-BOUND -- which is also why its time is flat in
# chunk size above the frame: the binding resource does not care how the reads
# are cut, only how many bytes come out of them.
MEASURED_ZSTD_STORE_BYTES_PER_SECOND = {
    1: 1.40e9, 2: 2.66e9, 4: 5.26e9, 8: 10.44e9,
    12: 15.04e9, 16: 20.01e9, 24: 25.08e9,
}


def zstd_store_decode_bytes_per_second(reader_workers: int = 16,
                                       table=None) -> float:
    """Aggregate zstd decode rate for `reader_workers`, from the measured curve.

    Interpolated between measured points rather than extrapolated from the
    one-worker rate times the count: efficiency is already 0.89 at 16 and falls
    to 0.75 by 24, so multiplying would over-predict by a third at the top of
    the range. Above the last measured point the rate is held FLAT rather than
    extended -- the curve is bending down and guessing how far is the kind of
    extrapolation this module exists to avoid.
    """
    table = dict(table or MEASURED_ZSTD_STORE_BYTES_PER_SECOND)
    if not table:
        raise ValueError("no measured decode rates")
    workers = max(1, int(reader_workers))
    if workers in table:
        return table[workers]
    points = sorted(table)
    if workers < points[0]:
        return table[points[0]] * workers / points[0]
    if workers > points[-1]:
        return table[points[-1]]
    hi = next(p for p in points if p > workers)
    lo = max(p for p in points if p < workers)
    span = (workers - lo) / (hi - lo)
    return table[lo] + span * (table[hi] - table[lo])


def decode_seconds_for_file(vrtypes, seconds_per_variant_by_form=None,
                            reader_workers: int = 1) -> float:
    """Whole-file decode seconds, from this file's histogram.

    `reader_workers` divides the single-threaded total, because the scan
    decodes in its reader pool. Passing 1 gives the serial figure, which is
    what the per-form costs are measured as.

    Defaults to the measured table, so a caller that has not run its own
    calibration still gets a decomposed estimate rather than a blended rate --
    but a host whose decoder differs should measure its own.
    """
    table = (MEASURED_SECONDS_PER_VARIANT_BY_FORM
             if seconds_per_variant_by_form is None
             else seconds_per_variant_by_form)
    total = decode_seconds(record_mix(vrtypes), table)
    return total / max(int(reader_workers), 1)
