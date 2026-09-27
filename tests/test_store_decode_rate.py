"""The store's decode rate, and the straddling term built on top of it.

Both numbers here were MEASURED on an idle host against the real store, and
both replaced a number that was inferred. The inferred one was not merely
imprecise -- `decode_bytes_per_second = 11.3e9` put 6.99 s of decode inside a
scan that measures 4.68 s end to end, which is impossible on its face. These
tests exist so that the next person to touch the rate has to confront the
measurement rather than the guess.
"""
import unittest
from math import gcd

from torchgwas.decode_model import (MEASURED_ZSTD_STORE_BYTES_PER_SECOND,
                                    zstd_store_decode_bytes_per_second)
from torchgwas.pipeline_model import explain_time


class ZstdStoreRateTests(unittest.TestCase):
    def test_the_default_matches_the_reader_pool_the_scan_actually_runs(self):
        """16 workers is the scan's default, and 20.01 GB/s is what it got."""
        self.assertAlmostEqual(
            zstd_store_decode_bytes_per_second(16) / 1e9, 20.01, places=2)

    def test_it_is_not_the_impossible_inferred_rate(self):
        """The whole point: 11.3 GB/s cannot be right.

        78.97 GB of raw hard calls at 11.3 GB/s is 6.99 s of decode, and the
        scan containing it measures 4.68 s. A part cannot exceed its whole.
        """
        raw_bytes = 8_931_083 * 8842
        measured_scan_seconds = 4.68
        rate = zstd_store_decode_bytes_per_second(16)
        self.assertLess(
            raw_bytes / rate, measured_scan_seconds,
            "decode alone cannot take longer than the scan that contains it")

    def test_decode_is_the_binding_resource_at_this_scale(self):
        """Why the curve is flat in chunk: the bottleneck ignores chunking.

        At 84% of the wall this is not a detail -- it is the reason a chunk
        size above the frame makes no difference, and the reason one below it
        makes such a large one.
        """
        raw_bytes = 8_931_083 * 8842
        share = (raw_bytes / zstd_store_decode_bytes_per_second(16)) / 4.68
        self.assertGreater(share, 0.75)
        self.assertLess(share, 1.0)

    def test_scaling_is_sublinear_and_is_not_extrapolated_upward(self):
        """Efficiency falls to 0.75 by 24 workers, so multiplying is wrong."""
        single = zstd_store_decode_bytes_per_second(1)
        self.assertLess(zstd_store_decode_bytes_per_second(24), 24 * single)
        # Beyond the measured range the rate is HELD, not extended.
        top = max(MEASURED_ZSTD_STORE_BYTES_PER_SECOND)
        self.assertEqual(zstd_store_decode_bytes_per_second(top + 40),
                         zstd_store_decode_bytes_per_second(top))

    def test_interpolation_stays_between_its_neighbours(self):
        for workers in (3, 6, 10, 14, 20):
            rate = zstd_store_decode_bytes_per_second(workers)
            points = sorted(MEASURED_ZSTD_STORE_BYTES_PER_SECOND)
            lo = max(p for p in points if p < workers)
            hi = min(p for p in points if p > workers)
            self.assertGreaterEqual(rate, MEASURED_ZSTD_STORE_BYTES_PER_SECOND[lo])
            self.assertLessEqual(rate, MEASURED_ZSTD_STORE_BYTES_PER_SECOND[hi])

    def test_zero_or_negative_workers_do_not_produce_a_zero_rate(self):
        """A zero rate divides into infinite seconds and poisons the model."""
        for workers in (0, -4):
            self.assertGreater(zstd_store_decode_bytes_per_second(workers), 0)


class StraddleAmplificationTests(unittest.TestCase):
    """`1 + (frame - gcd(chunk, frame)) / chunk`, against measurement.

    Measured in the decode path alone -- no GPU, no pipeline, no reader pool --
    by reading the same spans through `HardcallStore` at each chunk size:

        chunk      1024  1536  2048  3072  4096  8192  16384
        predicted  2.00  2.00  1.00  1.33  1.00  1.00   1.00
        measured   2.08  2.85  0.99  1.40  1.00  1.03   1.05
    """

    FRAME = 2048
    MEASURED_IN_DECODE = {1024: 2.08, 1536: 2.85, 2048: 0.99, 3072: 1.40,
                          4096: 1.00, 8192: 1.03, 16384: 1.05}

    @staticmethod
    def amplification(chunk, frame):
        return 1.0 + (frame - gcd(chunk, frame)) / chunk

    def test_it_is_exact_for_every_chunk_except_the_one_below_a_frame(self):
        """Six of seven within 8%. 1,536 is the stated exception."""
        for chunk, measured in sorted(self.MEASURED_IN_DECODE.items()):
            if chunk == 1536:
                continue
            predicted = self.amplification(chunk, self.FRAME)
            self.assertLess(
                abs(predicted - measured) / measured, 0.08,
                f"chunk {chunk}: predicted {predicted:.2f} against measured "
                f"{measured:.2f}")

    def test_the_known_shortfall_at_1536_is_recorded_not_papered_over(self):
        """A failing prediction stays visible, with its size written down.

        If someone closes this gap with a coefficient, this test tells them
        what they are actually fitting. The residual is real -- 1,536 costs
        2.85x in decode where the frame arithmetic allows 2.00x -- and the
        honest state of the model is that it does not know why.
        """
        predicted = self.amplification(1536, self.FRAME)
        self.assertAlmostEqual(predicted, 2.00, places=2)
        self.assertGreater(self.MEASURED_IN_DECODE[1536] / predicted, 1.4)

    def test_any_multiple_of_the_frame_is_free(self):
        for multiple in range(1, 12):
            self.assertEqual(
                self.amplification(self.FRAME * multiple, self.FRAME), 1.0)

    def test_the_model_ranks_aligned_above_straddling(self):
        """The property that makes the term useful for CHOOSING a chunk.

        Absolute seconds carry a uniform ~1.13x over-prediction on this
        machine, which cancels in a comparison; the ranking does not.
        """
        rates = dict(
            disk_bytes_per_second=5.56e9, h2d_bytes_per_second=55.26e9,
            d2h_bytes_per_second=26.0e9, write_bytes_per_second=4.86e9,
            gemm_flops_per_second={36: 8.35e12, 60: 8.71e12,
                                   156: 20.57e12, 540: 20.28e12})

        def predict(chunk):
            return explain_time(
                variants=8_931_083, samples=35_365, traits=128,
                covariate_rank=27, chunk_variants=chunk,
                stored_bytes=10_370_816_693, transfer_bytes_per_variant=8842,
                output_bytes_per_test=0.0, decoded_bytes_per_variant=8842,
                decode_bytes_per_second=zstd_store_decode_bytes_per_second(16),
                overlap=0.93, frame_variants=self.FRAME,
                **rates)["end_to_end_seconds"]

        aligned = [predict(c) for c in (2048, 4096, 8192, 16384)]
        straddling = [predict(c) for c in (1024, 1536, 3072)]
        self.assertLess(
            max(aligned), min(straddling),
            "every frame-aligned chunk must be predicted faster than every "
            "straddling one; that is the verdict the measurement gives")

    def test_without_a_frame_the_term_is_inert(self):
        """An unframed source must not pay a penalty it cannot incur."""
        base = dict(
            variants=1_000_000, samples=10_000, traits=8, covariate_rank=5,
            stored_bytes=int(1e9), transfer_bytes_per_variant=2500,
            output_bytes_per_test=8.0, decoded_bytes_per_variant=2500,
            decode_bytes_per_second=5e9, disk_bytes_per_second=2e9,
            h2d_bytes_per_second=20e9, d2h_bytes_per_second=15e9,
            write_bytes_per_second=1e9,
            gemm_flops_per_second={36: 8e12, 60: 8e12, 156: 20e12, 540: 20e12})
        without = explain_time(chunk_variants=1536, **base)
        zeroed = explain_time(chunk_variants=1536, frame_variants=0, **base)
        self.assertEqual(without["end_to_end_seconds"],
                         zeroed["end_to_end_seconds"])


if __name__ == "__main__":
    unittest.main()


class StraddleCaveatTests(unittest.TestCase):
    """The model must SAY when the chunk cannot be fixed by choosing better."""

    BASE = dict(
        variants=8_931_083, samples=35_365, traits=128, covariate_rank=27,
        stored_bytes=10_370_816_693, transfer_bytes_per_variant=8842,
        output_bytes_per_test=0.0, decoded_bytes_per_variant=8842,
        decode_bytes_per_second=20.01e9, overlap=0.93,
        disk_bytes_per_second=5.56e9, h2d_bytes_per_second=55.26e9,
        d2h_bytes_per_second=26.0e9, write_bytes_per_second=4.86e9,
        gemm_flops_per_second={36: 8.35e12, 60: 8.71e12,
                               156: 20.57e12, 540: 20.28e12})

    def test_a_chunk_below_one_frame_names_the_encoder_not_the_chunk(self):
        """No chunk size fixes it, so the advice must not be "tune the chunk".

        Every size below a frame straddles. A user told only that their chunk
        is bad will tune it forever against a floor they cannot move; the
        actionable fact is that the STORE was built with a frame their device
        cannot afford.
        """
        got = explain_time(chunk_variants=1024, frame_variants=2048,
                           **self.BASE)
        text = " ".join(got["caveats"])
        self.assertIn("re-encode", text)
        self.assertIn("frame_variants", text)

    def test_a_straddling_chunk_above_a_frame_names_the_aligned_neighbour(self):
        got = explain_time(chunk_variants=3072, frame_variants=2048,
                           **self.BASE)
        text = " ".join(got["caveats"])
        self.assertIn("nearest aligned chunk is 2048", text)

    def test_an_aligned_chunk_is_not_warned_about(self):
        for chunk in (2048, 4096):
            got = explain_time(chunk_variants=chunk, frame_variants=2048,
                               **self.BASE)
            self.assertEqual(got["frame_decode_amplification"], 1.0)
            self.assertNotIn("straddle", " ".join(got["caveats"]))

    def test_the_amplification_is_reported_so_a_caller_can_act_on_it(self):
        got = explain_time(chunk_variants=1536, frame_variants=2048,
                           **self.BASE)
        self.assertAlmostEqual(got["frame_decode_amplification"], 2.0, places=2)
