"""plink2's dense/sparse decision, as the source states it, on known inputs.

These pin the predicate transcribed in
`benchmarks/direct_plink2_sparse_predicate.py` from plink2_glm_linear.cc
(commit ca0f464): with covariates, a variant is sparse iff it has no missing
call and the previous processed variant had none either; without covariates
every autosomal variant is sparse; chrX is never sparse; constant-allele
variants are skipped and pay neither path.

The earlier model derived a "dense fraction" from the PGEN record form. One
test here writes a file whose record forms are all difflist -- the form the
old rule called sparse -- with a missing call in every record, and checks the
predicate calls every one of them dense. That is the case that separated the
two rules by a factor of two on the real cohort.
"""
from __future__ import annotations

import importlib.util
import pathlib
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from test_pgen_native_reader import write_pgen, write_pgen_mixed  # noqa: E402

from torchgwas import pgen_native  # noqa: E402
from torchgwas.pgen_reader import pack_genovec  # noqa: E402


def _load():
    path = pathlib.Path(__file__).resolve().parents[1] / "benchmarks" / "direct_plink2_sparse_predicate.py"
    spec = importlib.util.spec_from_file_location("plink2_sparse_predicate", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class GenocountTests(unittest.TestCase):
    def setUp(self):
        self.mod = _load()

    def test_counts_match_bincount_including_padding(self):
        rng = np.random.default_rng(7)
        for sample_ct in (1, 2, 3, 4, 5, 37, 64, 101):
            cats = rng.integers(0, 4, size=(9, sample_ct)).astype(np.uint8)
            rows = np.stack([pack_genovec(r, sample_ct) for r in cats])
            # Wider buffer with garbage past the genovec and in the padding.
            wide = np.full((9, rows.shape[1] + 3), 0xFF, dtype=np.uint8)
            wide[:, :rows.shape[1]] = rows
            rem = sample_ct % 4
            if rem:
                wide[:, rows.shape[1] - 1] |= (0xFF << (2 * rem)) & 0xFF
            got = self.mod.genocounts_from_packed(wide, sample_ct)
            want = np.stack([np.bincount(r, minlength=4) for r in cats])
            np.testing.assert_array_equal(got, want, err_msg=f"sample_ct={sample_ct}")

    def test_constant_allele_mask_matches_source_clauses(self):
        counts = np.array([
            [10, 0, 0, 0],   # all hom-ref: skipped
            [0, 10, 0, 0],   # all het: skipped
            [0, 0, 10, 0],   # all hom-alt: skipped
            [0, 0, 0, 10],   # all missing: skipped (0 het, 0 hom)
            [9, 1, 0, 0],    # polymorphic: regressed
            [5, 0, 5, 0],    # regressed
            [0, 5, 5, 0],    # regressed
        ])
        np.testing.assert_array_equal(
            self.mod.constant_allele_mask(counts),
            [True, True, True, True, False, False, False])


class PredicateTests(unittest.TestCase):
    def setUp(self):
        self.mod = _load()

    def test_with_covariates_missing_and_first_after_missing_are_dense(self):
        missing = np.array([0, 0, 1, 0, 0, 1, 1, 0])
        dense = self.mod.dense_path_mask(missing, covariate_ct=27)
        # index 0: slice start, prev_nm == 0 -> dense
        # index 1: prev_nm == 1, no missing -> sparse
        # index 2: missing -> dense, prev_nm := 0
        # index 3: no missing but prev_nm == 0 -> dense, prev_nm := 1
        # index 4: sparse
        # index 5, 6: missing -> dense
        # index 7: first after missing -> dense
        np.testing.assert_array_equal(dense, [1, 0, 1, 1, 0, 1, 1, 1])

    def test_without_covariates_nothing_is_dense(self):
        missing = np.array([0, 3, 1, 0, 0, 7])
        np.testing.assert_array_equal(
            self.mod.dense_path_mask(missing, covariate_ct=0), np.zeros(6, bool))

    def test_chrx_is_always_dense_whatever_missingness(self):
        missing = np.array([0, 0, 0, 0])
        is_x = np.array([False, True, True, False])
        dense = self.mod.dense_path_mask(missing, 27, is_x=is_x)
        np.testing.assert_array_equal(dense, [1, 1, 1, 0])
        # chrX blocks the sparse path even with no covariates (2761 && 2948).
        np.testing.assert_array_equal(
            self.mod.dense_path_mask(missing, 0, is_x=is_x), [0, 1, 1, 0])

    def test_skipped_variants_pay_neither_path_but_missing_ones_clear_prev_nm(self):
        missing = np.array([0, 0, 0, 0, 0, 5, 0])
        skipped = np.array([0, 0, 1, 0, 0, 1, 0], dtype=bool)
        dense = self.mod.dense_path_mask(missing, 27, skipped=skipped)
        # 0 dense (start), 1 sparse, 2 skipped (no missing: prev_nm kept),
        # 3 sparse, 4 sparse, 5 skipped WITH missing: prev_nm := 0,
        # 6 dense.
        np.testing.assert_array_equal(dense, [1, 0, 0, 0, 0, 0, 1])

    def test_slice_starts_restart_prev_nm(self):
        missing = np.zeros(6, dtype=int)
        dense = self.mod.dense_path_mask(missing, 27, slice_starts=(0, 3))
        np.testing.assert_array_equal(dense, [1, 0, 0, 1, 0, 0])

    def test_census_fractions_count_regressed_variants_only(self):
        counts = np.array([
            [10, 0, 0, 0],   # skipped
            [8, 1, 1, 0],    # dense (opens the run)
            [8, 2, 0, 0],    # sparse
            [7, 2, 0, 1],    # dense (missing)
            [8, 2, 0, 0],    # dense (first after missing)
            [8, 1, 1, 0],    # sparse
        ])
        c = self.mod.census_from_counts(counts, sample_ct=10, covariate_ct=27)
        self.assertEqual(c.constant_allele_skipped, 1)
        self.assertEqual(c.dense_variants, 3)
        self.assertEqual(c.sparse_variants, 2)
        self.assertAlmostEqual(c.dense_fraction, 3 / 5)
        self.assertAlmostEqual(c.dense_fraction_no_covariates, 0.0)
        # Only the dense variant WITH the missing call reaches the Gram.
        self.assertEqual(c.gram_variants, 1)
        self.assertAlmostEqual(c.gram_fraction, 1 / 5)
        self.assertEqual(c.variants_with_missing, 1)
        # Minor carriers: ALT is minor everywhere here -> het + hom-alt.
        self.assertAlmostEqual(c.mean_minor_carriers_sparse, 2.0)

    def test_minor_carriers_follow_the_omitted_major_allele(self):
        # ALT is the MAJOR allele: carriers of the minor (REF) allele are
        # het + hom-ref, per plink2_glm.cc:2622 and the genovec inversion.
        counts = np.array([[1, 2, 7, 0], [1, 2, 7, 0]])
        c = self.mod.census_from_counts(counts, sample_ct=10, covariate_ct=27)
        self.assertAlmostEqual(c.mean_minor_carriers_all, 3.0)


@unittest.skipUnless(pgen_native.available(), "native PGEN decoder not built")
class FileCensusTests(unittest.TestCase):
    """End to end on real PGEN bytes through the native decoder."""

    def setUp(self):
        self.mod = _load()

    def test_all_difflist_records_with_a_missing_call_each_are_all_dense(self):
        """The case that separates the source predicate from the vrtype rule.

        Form 4 is a difflist -- 'sparse-eligible' under the retired rule. With
        one missing call per record and 27 covariates, plink2 takes the dense
        path on every one of them (missing_ct != 0 -> 2972 -> dense).
        """
        rng = np.random.default_rng(11)
        variants, samples = 30, 53
        cats = (rng.random((variants, samples)) < 0.05).astype(np.uint8)  # mostly hom-ref
        cats[:, 0] = 3  # one missing call per variant
        cats[:, 1] = 1  # keep every variant polymorphic
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "difflist.pgen"
            write_pgen_mixed(path, cats, [4] * variants)
            counts = self.mod.count_file(str(path), chunk=7, workers=2)
            want = np.stack([np.bincount(r, minlength=4) for r in cats])
            np.testing.assert_array_equal(counts, want)
            c = self.mod.census_from_counts(counts, samples, 27)
            self.assertEqual(c.dense_fraction, 1.0)
            self.assertEqual(c.variants_with_missing, variants)

    def test_mixed_forms_count_exactly_and_ld_replay_is_handled(self):
        rng = np.random.default_rng(5)
        variants, samples = 24, 41
        cats = rng.integers(0, 3, size=(variants, samples)).astype(np.uint8)
        cats[::5, 3] = 3
        forms = [0, 2, 3, 1, 4, 2, 0, 3] * 3
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "mixed.pgen"
            write_pgen_mixed(path, cats, forms)
            # Chunks that start on LD-compressed records force the replay path.
            counts = self.mod.count_file(str(path), chunk=5, workers=3)
            want = np.stack([np.bincount(r, minlength=4) for r in cats])
            np.testing.assert_array_equal(counts, want)

    def test_missing_free_plain_records_with_covariates_are_sparse_after_the_first(self):
        rng = np.random.default_rng(2)
        cats = rng.integers(0, 3, size=(12, 20)).astype(np.uint8)
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "plain.pgen"
            write_pgen(path, cats)
            counts = self.mod.count_file(str(path), chunk=4, workers=1)
            c = self.mod.census_from_counts(counts, 20, 27)
            # Plain genovec -- 'dense' under the retired rule -- is sparse here
            # except for the slice-opening variant.
            self.assertEqual(c.dense_variants, 1)
            self.assertEqual(c.sparse_variants, 12 - 1 - c.constant_allele_skipped)


if __name__ == "__main__":
    unittest.main()
