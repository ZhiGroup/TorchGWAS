"""An impossible plan must be refused before reading, with a workable setting.

Supplementary Methods S3.1 lists as a limitation that "users must currently
split phenotype groups manually if the processed matrix does not fit". Today
the unreduced path has no feasibility check at all: it runs until CUDA raises
an out-of-memory error, after the genotype read has begun, naming neither the
cause nor a setting that would work.
"""
from __future__ import annotations

import unittest

from torchgwas.preflight import (PlanTooLarge, check_plan,
                                 largest_fitting_traits, require_fit)

# The real cohort, and the voxel panel the stress test targets.
COHORT = dict(chunk_variants=4096, depth=32, n_samples=35_365,
              covariate_rank=27, transfer_bytes_per_variant=8842.0)
H100_BYTES = 85.0e9


class FitsTests(unittest.TestCase):
    def test_an_ordinary_panel_fits(self):
        report = check_plan(n_traits=512, device_memory_bytes=H100_BYTES,
                            reduced=False, **COHORT)
        self.assertTrue(report["fits"])
        self.assertIsNone(report["largest_fitting_traits"])

    def test_a_fitting_plan_is_not_refused(self):
        report = require_fit(n_traits=512, device_memory_bytes=H100_BYTES,
                             reduced=False, **COHORT)
        self.assertTrue(report["fits"])


class RefusalTests(unittest.TestCase):
    """The voxel case: 2,085,000 traits cannot fit on any single card."""

    def test_a_voxel_panel_is_refused_when_unreduced(self):
        with self.assertRaises(PlanTooLarge) as caught:
            require_fit(n_traits=2_085_000, device_memory_bytes=H100_BYTES,
                        reduced=False, **COHORT)
        message = str(caught.exception)
        # The message must name the shortfall AND a way forward, because
        # "out of memory" alone is what this replaces.
        self.assertIn("GB", message)
        self.assertIn("reduce", message)

    def test_the_refusal_carries_the_numbers_not_just_prose(self):
        with self.assertRaises(PlanTooLarge) as caught:
            require_fit(n_traits=2_085_000, device_memory_bytes=H100_BYTES,
                        reduced=False, **COHORT)
        error = caught.exception
        self.assertGreater(error.predicted_bytes, error.available_bytes)
        self.assertEqual(error.available_bytes, H100_BYTES)

    def test_it_names_the_largest_trait_count_that_would_work(self):
        """A number the caller can act on, not 'try something smaller'."""
        with self.assertRaises(PlanTooLarge) as caught:
            require_fit(n_traits=2_085_000, device_memory_bytes=H100_BYTES,
                        reduced=False, **COHORT)
        largest = caught.exception.largest_traits
        self.assertIsNotNone(largest)
        self.assertGreater(largest, 0)
        self.assertLess(largest, 2_085_000)
        # And that number must itself actually fit.
        self.assertTrue(
            check_plan(n_traits=largest, device_memory_bytes=H100_BYTES,
                       reduced=False, **COHORT)["fits"])
        # One more trait than it should not.
        self.assertFalse(
            check_plan(n_traits=largest + 1, device_memory_bytes=H100_BYTES,
                       reduced=False, **COHORT)["fits"])


class ReducedPathTests(unittest.TestCase):
    """With a reduction, blocking is legal, so the same panel proceeds."""

    def test_a_reduced_voxel_panel_is_blocked_rather_than_refused(self):
        report = check_plan(n_traits=2_085_000, device_memory_bytes=H100_BYTES,
                            reduced=True, **COHORT)
        self.assertFalse(report["fits"])
        self.assertIsNotNone(report["auto_trait_block"])
        self.assertLess(report["auto_trait_block"], 2_085_000)

    def test_require_fit_lets_a_blockable_plan_through(self):
        report = require_fit(n_traits=2_085_000,
                             device_memory_bytes=H100_BYTES,
                             reduced=True, **COHORT)
        self.assertIsNotNone(report["auto_trait_block"])

    def test_blocking_is_why_reduce_is_suggested_in_the_refusal(self):
        """The advice has to be advice that works."""
        with self.assertRaises(PlanTooLarge):
            require_fit(n_traits=2_085_000, device_memory_bytes=H100_BYTES,
                        reduced=False, **COHORT)
        # Same panel, reduction on: no refusal.
        require_fit(n_traits=2_085_000, device_memory_bytes=H100_BYTES,
                    reduced=True, **COHORT)


class SmallDeviceTests(unittest.TestCase):
    def test_a_tiny_card_that_cannot_hold_one_trait_says_so(self):
        """The advice must change when no trait count rescues it."""
        with self.assertRaises(PlanTooLarge) as caught:
            require_fit(n_traits=128, device_memory_bytes=1.0e8,
                        reduced=False, **COHORT)
        self.assertIn("chunk", str(caught.exception).lower())

    def test_largest_fitting_returns_zero_when_nothing_fits(self):
        got = largest_fitting_traits(
            n_traits=128, device_memory_bytes=1.0e8, **COHORT)
        self.assertEqual(got, 0)

    def test_the_2080ti_holds_far_fewer_traits_than_the_h100(self):
        """11.4 GB against 85 GB, which the planner must reflect."""
        small = largest_fitting_traits(n_traits=1_000_000,
                                       device_memory_bytes=11.4e9, **COHORT)
        large = largest_fitting_traits(n_traits=1_000_000,
                                       device_memory_bytes=85.0e9, **COHORT)
        self.assertLess(small, large)


if __name__ == "__main__":
    unittest.main()
