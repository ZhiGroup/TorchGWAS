"""A chunk that straddles store frames pays for what it throws away.

`hardcall_store` compresses variants into independently decodable frames of
`frame_variants` rows, and a read decompresses every frame it overlaps. So the
cost of a chunk depends on whether it lines up with that grid, and the penalty
is not subtle. Measured on an idle H100, M=8,931,083, N=35,365, K=128, four
round-robin rounds, `sumstats=none`, against a store with frame 2,048:

    chunk    1024    1536    2048    3072    4096    8192   16384
    median  11.61   21.71    4.83    6.80    4.68    4.63    4.75
    frame     no      no     yes      no     yes     yes     yes
    penalty 2.48x   4.64x   1.03x   1.45x   1.00x   0.99x   1.01x

Seven for seven: every multiple of the frame size lands within noise of the
best, every straddling size is penalised. `auto_chunk_variants` bisects the
memory formula and can return ANY integer -- 1,537 is as reachable as 2,048 --
so before this it picked a penalised size whenever the memory arithmetic
happened to land on one, silently and with no way for a caller to tell.

These tests pin the alignment, not the seconds: the measurement above is one
machine and one store, but "the returned chunk is a whole number of frames"
is a property, and it is the property that makes the measurement moot.
"""
import unittest

from torchgwas.pipeline_model import (MAX_AUTO_CHUNK_VARIANTS,
                                      MIN_AUTO_CHUNK_VARIANTS,
                                      auto_chunk_variants)


# A regime where the BISECTION actually runs. At the measured cohort
# (N=35,365, two-bit on the wire) the ring is small enough that every budget
# from 2 GB up returns the cap, so a sweep there exercises only the early
# return and proves nothing about alignment. The first version of this file
# swept exactly that and its own premise test failed, which is the only reason
# the gap was noticed. N=250,000 sent as float32 puts the answer in the range
# that matters, and the unaligned values it produces are not hypothetical:
#
#   4 GB -> 286   8 GB -> 625   12 GB -> 965   16 GB -> 1305
#  24 GB -> 1985  32 GB -> 2665  48 GB -> 4025  64 GB -> 4096
#
# 4,025 against a 2,048 frame is the field version of the measured 1,536 case.
BISECTING_SAMPLES = 250_000
BUDGETS_GB = (4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64)


def choose(memory_bytes, frame_variants=None, samples=BISECTING_SAMPLES,
           traits=512, depth=8):
    """The planner's own call, with only the memory ceiling varied."""
    return auto_chunk_variants(
        n_samples=samples, n_traits=traits, covariate_rank=27,
        transfer_bytes_per_variant=4.0 * samples,
        device_memory_bytes=int(memory_bytes), depth=depth,
        frame_variants=frame_variants)


class FrameAlignmentTests(unittest.TestCase):
    def test_an_unframed_source_is_unconstrained(self):
        """No frame, no alignment: a flat `.bed` costs what it reads."""
        seen = set()
        for gigabytes in BUDGETS_GB:
            seen.add(choose(gigabytes * 1e9))
        # The point is that arbitrary sizes ARE reachable without a frame --
        # otherwise the aligned test below proves nothing.
        unaligned = [c for c in seen if c % 2048 and c != MAX_AUTO_CHUNK_VARIANTS]
        self.assertTrue(unaligned,
                        f"expected some unaligned chunk across the sweep, got {sorted(seen)}")

    def test_every_returned_chunk_is_a_whole_number_of_frames(self):
        """The property the measurement makes matter, swept over memory."""
        for frame in (1024, 2048, 4096):
            for gigabytes in BUDGETS_GB:
                chunk = choose(gigabytes * 1e9, frame_variants=frame)
                if chunk <= MIN_AUTO_CHUNK_VARIANTS or chunk < frame:
                    continue          # below one frame there is nothing to align
                self.assertEqual(
                    chunk % frame, 0,
                    f"frame {frame}, {gigabytes} GB -> chunk {chunk}, "
                    f"which straddles by {chunk % frame} variants")

    def test_the_measured_bad_sizes_are_never_returned(self):
        """1,536 and 3,072 cost 4.6x and 1.5x; they must be unreachable."""
        penalised = {1536, 3072}
        for gigabytes in BUDGETS_GB:
            chunk = choose(gigabytes * 1e9, frame_variants=2048)
            self.assertNotIn(
                chunk, penalised,
                f"{gigabytes} GB returned {chunk}, a size measured "
                f"{'4.6x' if chunk == 1536 else '1.5x'} slower than its "
                f"aligned neighbour")

    def test_alignment_never_raises_the_memory_the_chunk_needs(self):
        """Aligned DOWN. An aligned chunk must not exceed the unaligned one.

        Alignment is an optimisation and feasibility is not negotiable: the
        unaligned value is the one the bisection proved fits, so rounding up
        could breach the ceiling the caller asked to respect.
        """
        for gigabytes in BUDGETS_GB:
            loose = choose(gigabytes * 1e9)
            tight = choose(gigabytes * 1e9, frame_variants=2048)
            self.assertLessEqual(
                tight, loose,
                f"{gigabytes} GB: aligned {tight} exceeds unaligned {loose}")

    def test_a_frame_larger_than_the_budget_does_not_wedge_the_chunk_to_zero(self):
        """Below one frame, alignment is impossible, so it must not be tried.

        Every size straddles when the chunk is smaller than a frame, so there
        is nothing to gain -- and aligning down would floor the chunk to zero
        and take the scan with it.
        """
        chunk = choose(4e9, frame_variants=65_536)
        self.assertGreater(chunk, 0)
        self.assertGreaterEqual(chunk, MIN_AUTO_CHUNK_VARIANTS)

    def test_the_cap_is_aligned_too(self):
        """The fast path returns the cap directly; it gets the same treatment.

        4,096 is already a multiple of 2,048, so this would pass by accident.
        A frame of 3,000 is the case that catches a missed `snap` on the
        early return.
        """
        chunk = auto_chunk_variants(
            n_samples=1000, n_traits=4, covariate_rank=2,
            transfer_bytes_per_variant=250.0,
            device_memory_bytes=int(200e9), depth=4, frame_variants=3000)
        self.assertLessEqual(chunk, MAX_AUTO_CHUNK_VARIANTS)
        self.assertEqual(chunk % 3000, 0, f"cap returned unaligned {chunk}")

    def test_zero_or_negative_frames_are_ignored_not_crashed_on(self):
        """A malformed manifest must not divide by zero inside the planner."""
        for frame in (0, -1, None):
            chunk = choose(24e9, frame_variants=frame)
            self.assertGreater(chunk, 0)


class AlignmentIsReportedByTheSourceTests(unittest.TestCase):
    def test_a_plain_bed_reports_no_alignment(self):
        """Attribute must exist on the real class, not just in the planner.

        `linear.py` reads it with `getattr(..., None)`, which would silently
        accept a typo forever, so the name is pinned here.
        """
        from torchgwas.bed import PlinkBedGenotype

        self.assertTrue(hasattr(PlinkBedGenotype, "chunk_alignment_variants"))
        source = PlinkBedGenotype.__new__(PlinkBedGenotype)
        source._store = None
        self.assertIsNone(source.chunk_alignment_variants)

    def test_a_stored_bed_reports_its_frame(self):
        from torchgwas.bed import PlinkBedGenotype

        class FakeStore:
            frame_variants = 2048

        source = PlinkBedGenotype.__new__(PlinkBedGenotype)
        source._store = FakeStore()
        self.assertEqual(source.chunk_alignment_variants, 2048)


if __name__ == "__main__":
    unittest.main()
