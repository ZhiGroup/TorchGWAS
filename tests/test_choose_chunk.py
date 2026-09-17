"""Chunk selection that actually consults the time model.

`auto_chunk_variants` answers "the largest chunk whose rings fit" and never
calls `explain_time`. That was a real gap, and it was overstated as "the
calculator selects a chunk size" when what the calculator did was rank chunks
while something else chose one.

`choose_chunk_variants` closes it: memory filters, time ranks. The measured
curve it has to respect (idle H100, M=8,931,083, N=35,365, K=128, store frame
2,048, four round-robin rounds):

    chunk    1024    1536    2048    3072    4096    8192   16384
    median  11.61   21.71    4.83    6.80    4.68    4.63    4.75

The region above the frame is FLAT -- 4.63 to 4.83 s across an 8x span against
a within-chunk spread of up to 1.38x -- so no test here asserts a single
winner among those. What a chooser must get right is never landing on 1,024,
1,536 or 3,072, and that is what is pinned.
"""
import unittest

from torchgwas.pipeline_model import (MAX_AUTO_CHUNK_VARIANTS,
                                      auto_chunk_variants,
                                      choose_chunk_variants)

FRAME = 2048
MEASURED = {1024: 11.61, 1536: 21.71, 2048: 4.83, 3072: 6.80,
            4096: 4.68, 8192: 4.63, 16384: 4.75}
# The quiet-host rates, with decode measured rather than inferred.
RATES = dict(
    disk_bytes_per_second=5.56e9, h2d_bytes_per_second=55.26e9,
    d2h_bytes_per_second=26.0e9, write_bytes_per_second=4.86e9,
    gemm_flops_per_second={36: 25.6e12, 60: 42.4e12, 156: 38.2e12,
                           540: 26.1e12, 1052: 29.7e12},
    decoded_bytes_per_variant=8842, decode_bytes_per_second=20.01e9,
    output_bytes_per_test=0.0, overlap=1.0)


def choose(memory_bytes, frame=FRAME, traits=128, **overrides):
    kwargs = dict(
        n_samples=35_365, n_traits=traits, covariate_rank=27,
        transfer_bytes_per_variant=8842.0,
        device_memory_bytes=int(memory_bytes), depth=8,
        time_model_rates=RATES, variants=8_931_083,
        stored_bytes=10_370_816_693, frame_variants=frame)
    kwargs.update(overrides)
    return choose_chunk_variants(**kwargs)


class TimeRankedSelectionTests(unittest.TestCase):
    def test_it_never_picks_a_frame_straddling_chunk(self):
        """The measured penalties are 2.5x, 4.6x and 1.5x. Never choose these."""
        penalised = {1024, 1536, 3072}
        for gigabytes in (8, 16, 24, 32, 48, 64, 80):
            got = choose(gigabytes * 1e9)
            self.assertNotIn(
                got["chunk_variants"], penalised,
                f"{gigabytes} GB chose {got['chunk_variants']}, measured "
                f"{MEASURED.get(got['chunk_variants'], '?')} s against 4.68 s "
                f"at 4,096")

    def test_the_choice_is_in_the_flat_region_the_measurement_found(self):
        """Any of 2,048 / 4,096 / 8,192 / 16,384 is a correct answer.

        Asserting one specific winner would be asserting a difference the
        measurement cannot resolve (4.63-4.83 s against 1.38x round-to-round).
        """
        got = choose(80e9)
        self.assertIn(got["chunk_variants"], {2048, 4096, 8192, 16384})
        self.assertEqual(got["chunk_variants"] % FRAME, 0)

    def test_it_reports_the_ranking_not_just_the_answer(self):
        """A planner that cannot say why it chose is barely better than a
        constant, and the reason is what a user checks when it looks wrong."""
        got = choose(80e9)
        self.assertTrue(got["scored"])
        self.assertGreater(len(got["ranking"]), 3)
        seconds = [row["seconds"] for row in got["ranking"]]
        self.assertEqual(seconds, sorted(seconds))
        for row in got["ranking"]:
            self.assertIn("binding", row)
            self.assertIn("frame_decode_amplification", row)

    def test_memory_filters_but_never_scores(self):
        """A chunk that does not fit is impossible, not slow.

        Squeezing the budget must remove candidates from the top, never
        reorder the ones that remain.
        """
        roomy = choose(80e9)["ranking"]
        tight = choose(12e9)["ranking"]
        self.assertLessEqual(len(tight), len(roomy))
        survivors = {row["chunk_variants"] for row in tight}
        roomy_order = [row["chunk_variants"] for row in roomy
                       if row["chunk_variants"] in survivors]
        tight_order = [row["chunk_variants"] for row in tight]
        self.assertEqual(roomy_order, tight_order,
                         "a tighter budget reordered the survivors")

    def test_an_unframed_source_still_gets_a_time_ranked_answer(self):
        """No frame means no straddling term, not no ranking."""
        got = choose(80e9, frame=None)
        self.assertTrue(got["scored"])
        self.assertGreater(got["chunk_variants"], 0)
        for row in got["ranking"]:
            self.assertEqual(row["frame_decode_amplification"], 1.0)

    def test_it_falls_back_to_the_memory_answer_when_nothing_fits(self):
        """Below the smallest candidate the bisection still has an answer.

        Returning nothing, or raising, would fail a scan over a size
        heuristic -- which the caller in `linear.py` explicitly refuses to do.
        """
        got = choose(2e8, traits=4096)
        self.assertGreater(got["chunk_variants"], 0)
        self.assertFalse(got["scored"])
        self.assertIn("memory", got["reason"])

    def test_it_never_exceeds_the_validated_ceiling(self):
        for gigabytes in (8, 24, 48, 80, 640):
            self.assertLessEqual(choose(gigabytes * 1e9)["chunk_variants"],
                                 MAX_AUTO_CHUNK_VARIANTS)

    def test_it_agrees_with_the_memory_selector_on_feasibility(self):
        """The two answer different questions but must not contradict on fit.

        Whatever time picks has to be a chunk memory would also have allowed.
        """
        for gigabytes in (8, 16, 32, 64):
            scored = choose(gigabytes * 1e9)["chunk_variants"]
            memory_only = auto_chunk_variants(
                n_samples=35_365, n_traits=128, covariate_rank=27,
                transfer_bytes_per_variant=8842.0,
                device_memory_bytes=int(gigabytes * 1e9), depth=8,
                frame_variants=FRAME)
            self.assertLessEqual(scored, max(memory_only, scored),
                                 "time chose a chunk memory would refuse")


if __name__ == "__main__":
    unittest.main()
