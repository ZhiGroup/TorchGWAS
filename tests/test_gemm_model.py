"""The GEMM bound, checked against published specs and measured shapes.

Two independent kinds of check, and both matter:

  * the PEAK must reproduce the vendor's headline number for three different
    architectures from nothing but SM count, lanes and clock. If it does not,
    the lanes-per-SM table is wrong and every prediction built on it is wrong
    by a factor of two;
  * the BOUND must not be violated by any measured shape, and must not be so
    loose as to be useless. A roofline that measurement exceeds is not a
    roofline.
"""
import unittest

from torchgwas.gemm_model import (KNOWN_DEVICES, DeviceSpec,
                                  gemm_bound_flops_per_second,
                                  gemm_flops_per_second,
                                  tile_intensity_flops_per_byte)


class PublishedSpecTests(unittest.TestCase):
    """Peak from first principles, against what the vendor advertises."""

    EXPECTED_TFLOPS = {"H100": 66.9, "A100": 19.5, "2080Ti": 13.4}

    def test_peak_matches_the_published_figure_on_three_architectures(self):
        for key, expected in self.EXPECTED_TFLOPS.items():
            got = KNOWN_DEVICES[key].peak_flops_per_second / 1e12
            self.assertAlmostEqual(
                got, expected, delta=0.2,
                msg=f"{key}: derived {got:.1f} TFLOP/s against published "
                    f"{expected} -- check FP32_LANES_PER_SM")

    def test_lanes_differ_by_architecture_and_that_is_the_whole_point(self):
        """A100 and H100 differ by more than their SM counts suggest.

        108 -> 132 SMs is 1.22x, but the peak ratio is 3.4x, because Hopper
        doubled the FP32 lanes per SM. A model that scaled by SM count alone
        would be wrong by a factor of two, which is exactly the error the
        lookup table it replaces could not even express.
        """
        self.assertEqual(KNOWN_DEVICES["A100"].fp32_lanes_per_sm, 64)
        self.assertEqual(KNOWN_DEVICES["H100"].fp32_lanes_per_sm, 128)

    def test_an_unknown_capability_refuses_rather_than_guessing(self):
        with self.assertRaises(ValueError):
            DeviceSpec("invented", 100, (99, 9), 1.5e9, 1e12)
        # ...but an explicit value is accepted, because a caller who knows the
        # architecture should not be blocked by this module's table.
        spec = DeviceSpec("invented", 100, (99, 9), 1.5e9, 1e12,
                          fp32_lanes_per_sm=128)
        self.assertGreater(spec.peak_flops_per_second, 0)


class TileIntensityTests(unittest.TestCase):
    def test_intensity_is_independent_of_the_reduction_length(self):
        """K cancels, which is why the bound is a function of two dimensions."""
        self.assertAlmostEqual(tile_intensity_flops_per_byte(32, 32), 8.0)
        self.assertAlmostEqual(tile_intensity_flops_per_byte(128, 128), 32.0)

    def test_a_small_tile_is_memory_bound_on_an_h100(self):
        """8 FLOP/byte against a machine balance of 20: starved, and by how much.

        This is the term a global roofline misses. The whole product at chunk
        4,096 and width 156 has an intensity of 74.8 FLOP/byte and looks
        comfortably compute-bound; the 32x32 tiles computing it are not.
        """
        h100 = KNOWN_DEVICES["H100"]
        self.assertLess(tile_intensity_flops_per_byte(32, 32),
                        h100.machine_balance_flops_per_byte)
        self.assertGreater(tile_intensity_flops_per_byte(128, 128),
                           h100.machine_balance_flops_per_byte)


class BoundAgainstMeasurementTests(unittest.TestCase):
    """Measured `torch.mm`, width 156, samples 35,365, min of 7 per point.

    From `benchmarks/direct_gemm_chunk_curve.py` on an H100.
    """

    MEASURED_TFLOPS = {256: 24.72, 512: 33.66, 1024: 37.34, 2048: 35.25,
                       4096: 38.24, 8192: 37.23, 16384: 38.55}
    WIDTH = 156

    def bound(self, chunk):
        return gemm_bound_flops_per_second(
            KNOWN_DEVICES["H100"], chunk, self.WIDTH
        )['bound_flops_per_second'] / 1e12

    def test_no_measured_shape_exceeds_the_bound(self):
        """The defining property. A violated bound is a broken model."""
        for chunk, measured in sorted(self.MEASURED_TFLOPS.items()):
            self.assertGreater(
                self.bound(chunk), measured,
                f"chunk {chunk}: measured {measured} TFLOP/s exceeds the "
                f"bound {self.bound(chunk):.1f}, so a term is missing")

    def test_the_bound_is_tight_enough_to_be_useful(self):
        """Within 1.6x everywhere, and flat to 1.27x across a 64x chunk span.

        Looseness is the price of a bound. What would make it useless is
        looseness that VARIES with the shape, because then it cannot rank two
        shapes -- and ranking is what the calculator needs.
        """
        ratios = [self.bound(c) / m for c, m in self.MEASURED_TFLOPS.items()]
        self.assertLess(max(ratios), 1.7)
        self.assertLess(max(ratios) / min(ratios), 1.35,
                        "the bound's looseness must not vary much with shape")

    def test_it_does_not_clamp_beyond_the_measured_widths(self):
        """The failure that motivated all of this.

        The table it replaces holds four widths and returns the width-540 rate
        for anything wider, so K=1024 (width 1052) was priced at the K=512
        rate -- a 1.81x error on the ladder. A derived bound has no edge to
        fall off: the rate keeps responding to the shape.
        """
        device = KNOWN_DEVICES["H100"]
        at_540 = gemm_flops_per_second(device, 4096, 540)
        at_1052 = gemm_flops_per_second(device, 4096, 1052)
        at_4096 = gemm_flops_per_second(device, 4096, 4096)
        self.assertNotEqual(at_540, at_1052)
        self.assertNotEqual(at_1052, at_4096)
        # Wider designs have better reuse, so the rate must not fall.
        self.assertGreaterEqual(at_1052, at_540 * 0.95)

    def test_it_answers_for_a_device_that_is_not_present(self):
        """Planning a run on a card you do not have is most of the point."""
        for key in ("A100", "2080Ti"):
            rate = gemm_flops_per_second(KNOWN_DEVICES[key], 4096, 156)
            self.assertGreater(rate, 0)
            self.assertLess(rate, KNOWN_DEVICES[key].peak_flops_per_second)

    def test_the_realized_fraction_is_explicit_and_has_no_default(self):
        """A bound must not silently masquerade as a prediction.

        Callers get the bound unless they pass a measured fraction, and there
        is deliberately no default value: inventing one is how the old scalar
        44.88 TFLOPS came to over-estimate the GPU by 2x at K=128.
        """
        device = KNOWN_DEVICES["H100"]
        bound = gemm_flops_per_second(device, 4096, 156)
        realized = gemm_flops_per_second(device, 4096, 156,
                                         realized_fraction=0.70)
        self.assertAlmostEqual(realized, bound * 0.70, places=3)
        for bad in (0.0, -0.5, 1.5):
            with self.assertRaises(ValueError):
                gemm_flops_per_second(device, 4096, 156, realized_fraction=bad)


class WiredIntoTheTimeModelTests(unittest.TestCase):
    """A DeviceSpec must reach `explain_time`, not just sit beside it."""

    BASE = dict(
        variants=8_931_083, samples=35_365, covariate_rank=27,
        chunk_variants=4096, stored_bytes=78_968_635_889,
        transfer_bytes_per_variant=8842, output_bytes_per_test=0.0,
        decoded_bytes_per_variant=0.0, decode_bytes_per_second=0.0,
        overlap=1.0, disk_bytes_per_second=5.56e9,
        h2d_bytes_per_second=55.26e9, d2h_bytes_per_second=26.0e9,
        write_bytes_per_second=4.86e9)

    def test_the_clamp_that_cost_1_81x_at_k1024_is_gone(self):
        """The measured ladder rung the lookup table could not reach.

        Table: 32.79 s predicted against 18.09 s measured, because width 1052
        is priced at the width-540 rate. Derived: 14.22 s, which is wrong in
        the other direction and by far less -- and wrong for a different
        reason, since something else grows with K that is not the GEMM.
        """
        from torchgwas.pipeline_model import explain_time

        table = {36: 8.35e12, 60: 8.71e12, 156: 20.57e12, 540: 20.28e12}
        measured = 18.095
        with_table = explain_time(traits=1024, gemm_flops_per_second=table,
                                  **self.BASE)["end_to_end_seconds"]
        with_spec = explain_time(traits=1024,
                                 gemm_flops_per_second=KNOWN_DEVICES["H100"],
                                 **self.BASE)["end_to_end_seconds"]
        table_error = max(with_table / measured, measured / with_table)
        spec_error = max(with_spec / measured, measured / with_spec)
        self.assertGreater(table_error, 1.7)
        self.assertLess(spec_error, 1.35)
        self.assertLess(spec_error, table_error)

    def test_deriving_without_a_chunk_refuses_rather_than_assuming_one(self):
        """The bound needs both output dimensions; guessing one is silent harm."""
        from torchgwas.pipeline_model import gemm_rate_at_width

        with self.assertRaises(ValueError):
            gemm_rate_at_width(KNOWN_DEVICES["H100"], 156)
        self.assertGreater(
            gemm_rate_at_width(KNOWN_DEVICES["H100"], 156, chunk_variants=4096),
            0)

    def test_a_scalar_and_a_table_still_work(self):
        """The derived path is added, not substituted -- old callers unchanged."""
        from torchgwas.pipeline_model import gemm_rate_at_width

        self.assertEqual(gemm_rate_at_width(44.88e12, 156), 44.88e12)
        self.assertEqual(
            gemm_rate_at_width({36: 1e12, 540: 2e12}, 540), 2e12)


if __name__ == "__main__":
    unittest.main()
