"""Decode must be a weighted sum over record forms, not one blended rate.

The term it replaces was 11.3 GB/s, obtained by subtracting an estimated read
from a measured scan. That lumps five decode paths together, carries no error
bar, and cannot explain why two files of the same size decode at different
speeds. These tests pin the properties that make the replacement honest.
"""
from __future__ import annotations

import unittest

import numpy as np

from torchgwas.decode_model import (BOUNDED_FORMS, decode_seconds,
                                    difflist_share, record_mix)


class RecordMixTests(unittest.TestCase):
    def test_it_counts_the_low_nibble_only(self):
        """vrtype carries other tracks in the high bits; the form is 0x07."""
        vrtypes = np.array([0x00, 0x08, 0x10, 0x02, 0x42], dtype=np.uint8)
        # 0x08 and 0x10 are form 0 with other track bits set; 0x42 is form 2.
        self.assertEqual(record_mix(vrtypes), {0: 3, 2: 2})

    def test_absent_forms_are_not_reported_as_zero(self):
        """A key with count 0 invites pricing a form that is not there."""
        mix = record_mix(np.array([0, 0, 0], dtype=np.uint8))
        self.assertEqual(mix, {0: 3})
        self.assertNotIn(2, mix)

    def test_an_empty_file_has_an_empty_mix(self):
        self.assertEqual(record_mix(np.array([], dtype=np.uint8)), {})


class DifflistShareTests(unittest.TestCase):
    """How much of decode is data-dependent at all."""

    def test_all_plain_genovec_is_fully_bounded(self):
        self.assertEqual(difflist_share({0: 1000}), 0.0)

    def test_all_difflist_is_fully_data_dependent(self):
        self.assertEqual(difflist_share({4: 1000}), 1.0)

    def test_ld_compressed_counts_as_data_dependent(self):
        """Forms 2 and 3 copy a genovec AND apply a difflist of deltas."""
        self.assertNotIn(2, BOUNDED_FORMS)
        self.assertAlmostEqual(difflist_share({0: 50, 2: 50}), 0.5)

    def test_the_measured_cohort(self):
        """20.72% plain genovec on this cohort, so ~79% is data-dependent."""
        mix = {0: 2072, 2: 2128, 4: 5800}
        self.assertAlmostEqual(difflist_share(mix), 1 - 0.2072, places=3)

    def test_an_empty_mix_is_zero_not_a_division_error(self):
        self.assertEqual(difflist_share({}), 0.0)


class DecodeSecondsTests(unittest.TestCase):
    def test_it_is_a_weighted_sum(self):
        mix = {0: 100, 4: 50}
        per_form = {0: 1e-6, 4: 4e-6}
        self.assertAlmostEqual(decode_seconds(mix, per_form),
                               100 * 1e-6 + 50 * 4e-6)

    def test_a_form_with_no_measured_cost_RAISES(self):
        """Pricing an unmeasured record type at zero is how decode vanished.

        The model had no decode term at all and predicted the zstd store at
        7.63x where it measures 2.15x. Silence is the failure mode to prevent.
        """
        with self.assertRaises(KeyError) as caught:
            decode_seconds({0: 10, 2: 10}, {0: 1e-6})
        self.assertIn("2", str(caught.exception))

    def test_extra_timings_for_absent_forms_are_harmless(self):
        """A machine profile may time forms this particular file lacks."""
        self.assertAlmostEqual(
            decode_seconds({0: 10}, {0: 1e-6, 2: 5e-6, 4: 9e-6}), 10e-6)

    def test_a_cheap_form_and_an_expensive_one_do_not_average(self):
        """The whole point: one blended rate cannot express this.

        Same variant count, different mixes, different totals -- which a single
        bytes-per-second figure keyed on file size would call identical.
        """
        per_form = {0: 1e-6, 4: 10e-6}
        mostly_plain = decode_seconds({0: 900, 4: 100}, per_form)
        mostly_difflist = decode_seconds({0: 100, 4: 900}, per_form)
        self.assertGreater(mostly_difflist / mostly_plain, 4.0)


class MeasuredTableTests(unittest.TestCase):
    """The shipped per-form costs, and the whole-file helper that uses them."""

    def test_the_table_covers_every_form_the_cohort_contains(self):
        """The real file has forms 0, 1, 2, 3, 4 and 6.

        A missing entry would raise at scan time rather than silently price a
        record type at zero -- but it should not be missing in the first place.
        """
        from torchgwas.decode_model import MEASURED_SECONDS_PER_VARIANT_BY_FORM
        for form in (0, 1, 2, 3, 4, 6):
            self.assertIn(form, MEASURED_SECONDS_PER_VARIANT_BY_FORM)

    def test_plain_genovec_is_the_cheapest(self):
        """It is the only bounded form, so it must be."""
        from torchgwas.decode_model import MEASURED_SECONDS_PER_VARIANT_BY_FORM
        table = MEASURED_SECONDS_PER_VARIANT_BY_FORM
        self.assertEqual(min(table, key=table.get), 0)

    def test_the_spread_is_large_enough_to_matter(self):
        """13.5x measured. A single blended rate cannot express that.

        If this ever collapses toward 1, the decomposition has stopped earning
        its complexity and a scalar would do.
        """
        from torchgwas.decode_model import MEASURED_SECONDS_PER_VARIANT_BY_FORM
        values = MEASURED_SECONDS_PER_VARIANT_BY_FORM.values()
        self.assertGreater(max(values) / min(values), 5.0)

    def test_ld_compressed_is_cheaper_than_difflist_over_background(self):
        """The counter-intuitive part, pinned so it is not 'corrected' away.

        LD-compressed copies the previous genovec -- bounded -- plus a short
        difflist of deltas. A difflist over a constant background can carry
        far more entries. Measured 5.19 us against 32.67 us.
        """
        from torchgwas.decode_model import MEASURED_SECONDS_PER_VARIANT_BY_FORM
        table = MEASURED_SECONDS_PER_VARIANT_BY_FORM
        self.assertLess(table[2], table[4])

    def test_whole_file_seconds_divide_by_the_reader_pool(self):
        """The scan decodes in its reader threads, so a serial total misleads."""
        from torchgwas.decode_model import decode_seconds_for_file
        vrtypes = np.zeros(1000, dtype=np.uint8)      # all form 0
        serial = decode_seconds_for_file(vrtypes, reader_workers=1)
        pooled = decode_seconds_for_file(vrtypes, reader_workers=16)
        self.assertAlmostEqual(serial / 16, pooled, places=9)

    def test_a_worker_count_of_zero_does_not_divide_by_zero(self):
        from torchgwas.decode_model import decode_seconds_for_file
        vrtypes = np.zeros(10, dtype=np.uint8)
        self.assertGreater(decode_seconds_for_file(vrtypes, reader_workers=0),
                           0.0)

    def test_a_caller_can_supply_its_own_table(self):
        """A host whose decoder differs must not be stuck with ours."""
        from torchgwas.decode_model import decode_seconds_for_file
        vrtypes = np.zeros(100, dtype=np.uint8)
        got = decode_seconds_for_file(vrtypes, {0: 1.0}, reader_workers=1)
        self.assertAlmostEqual(got, 100.0)


class Plink2FormTests(unittest.TestCase):
    """The competitor model must use the same histogram, not its own constant."""

    def setUp(self):
        import importlib.util
        import pathlib
        path = pathlib.Path("benchmarks/direct_plink2_cost_model.py")
        if not path.exists():
            self.skipTest("competitor model not present")
        import sys
        spec = importlib.util.spec_from_file_location("plink2_model", path)
        self.model = importlib.util.module_from_spec(spec)
        # Register BEFORE executing: `@dataclass` resolves annotations through
        # `sys.modules[cls.__module__]`, so a module loaded by spec alone dies
        # with `'NoneType' object has no attribute '__dict__'` the moment it
        # defines one. The competitor model defines several.
        sys.modules[spec.name] = self.model
        self.addCleanup(sys.modules.pop, spec.name, None)
        spec.loader.exec_module(self.model)

    def test_the_record_mix_no_longer_sets_the_dense_fraction(self):
        """RETRACTED: `dense_fraction_from_mix` and DEFAULT_DENSE_PATH_FRACTION.

        Two tests here used to assert that plain-genovec plus LD-compressed
        records (20.72% + 21.28%) reproduce a 0.42 dense fraction. plink2
        consults the record form only when there are no covariates
        (plink2_glm_linear.cc:2762), and the branch itself is decided by
        `missing_ct == 0 && prev_nm` (2948, 2972). The helper and the constant
        are gone so nothing can quietly derive a regression-path fraction from
        a decode-path histogram again; `direct_plink2_sparse_predicate.py`
        computes the real one.
        """
        self.assertFalse(hasattr(self.model, "dense_fraction_from_mix"))
        self.assertFalse(hasattr(self.model, "DEFAULT_DENSE_PATH_FRACTION"))

    def test_a_cohort_cannot_inherit_a_dense_fraction_by_omission(self):
        """The 0.077 and the 0.42 were both inherited defaults. No default."""
        with self.assertRaises(TypeError):
            self.model.Cohort(variants=10, samples=10, covariates=2,
                              stored_bytes_per_variant=1.0, pvar_bytes=1.0)

    def test_dense_fraction_from_census_requires_provenance(self):
        census = {"dense_fraction": 0.99, "plink_ng_commit": "abc", "covariate_ct": 27}
        self.assertAlmostEqual(self.model.dense_fraction_from_census(census), 0.99)
        with self.assertRaises(KeyError):
            self.model.dense_fraction_from_census({"dense_fraction": 0.99})

    def test_decode_is_paid_once_per_subbatch(self):
        """plink2 re-reads and re-decodes per 240-trait subbatch."""
        mix = {0: 100}
        per_form = {0: 1e-6}
        one = self.model.decode_seconds_by_form(mix, per_form, passes_=1)
        three = self.model.decode_seconds_by_form(mix, per_form, passes_=3)
        self.assertAlmostEqual(three, 3 * one)

    def test_an_unmeasured_form_raises_here_too(self):
        with self.assertRaises(KeyError):
            self.model.decode_seconds_by_form({0: 5, 4: 5}, {0: 1e-6})


if __name__ == "__main__":
    unittest.main()
