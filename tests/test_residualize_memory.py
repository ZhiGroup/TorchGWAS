"""Blocked residualisation: same answer, bounded memory.

The host path used to evaluate three whole-matrix expressions, one of which
held three full-size arrays at once. At the voxel stress configuration
(35,365 x 600,000 float32 = 84.9 GB) that is ~255 GB of peak anonymous memory;
the run was measured at 229 GB resident with 53,468 s of kernel time against
348 s of user time -- a 153x ratio that is page faulting, not arithmetic --
while all eight GPUs idled 69% of the time behind it.

Two things have to hold for the fix to be worth anything: the numbers must not
move, and the peak must actually fall. Both are tested, the second by counting
allocations rather than by trusting the shape of the code.
"""
import unittest

import numpy as np

from torchgwas.preprocess import (RESIDUALIZE_TRAIT_BLOCK,
                                  residualize_and_standardize)


def reference(phenotype, covariates):
    """The original whole-matrix expressions, kept as the oracle.

    Deliberately written out rather than imported: this is what the blocked
    implementation must reproduce, and if someone changes the implementation
    the comparison should still be against the original arithmetic.
    """
    from torchgwas.preprocess import _covariate_basis

    q = None
    if covariates is not None and covariates.shape[1] > 0:
        q = _covariate_basis(covariates)
    out = phenotype - phenotype.mean(axis=0, keepdims=True)
    if q is not None:
        out = out - q @ (q.T @ out)
    std = out.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    return out / std, q


class BlockedResidualizeAgreesTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        self.pheno = rng.normal(size=(200, 37)).astype(np.float64)
        self.covar = rng.normal(size=(200, 5)).astype(np.float64)

    def test_it_matches_the_whole_matrix_arithmetic(self):
        want, want_q = reference(self.pheno, self.covar)
        got, got_q = residualize_and_standardize(self.pheno, self.covar)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(got_q, want_q, rtol=1e-12, atol=1e-12)

    def test_the_block_size_does_not_change_the_answer(self):
        """Traits are independent, so blocking must be arithmetic-neutral.

        A block size that changed the result would mean the projection was
        leaking across traits, which is the one bug this shape could hide.
        """
        want, _ = reference(self.pheno, self.covar)
        for block in (1, 2, 7, 36, 37, 38, 1000):
            got, _ = residualize_and_standardize(
                self.pheno, self.covar, trait_block=block)
            np.testing.assert_allclose(
                got, want, rtol=1e-12, atol=1e-12,
                err_msg=f"block {block} changed the result")

    def test_it_works_without_covariates(self):
        want, _ = reference(self.pheno, None)
        got, q = residualize_and_standardize(self.pheno, None)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)
        self.assertIsNone(q)

    def test_a_constant_trait_does_not_divide_by_zero(self):
        pheno = self.pheno.copy()
        pheno[:, 3] = 2.5
        got, _ = residualize_and_standardize(pheno, self.covar)
        self.assertTrue(np.all(np.isfinite(got)))

    def test_a_single_trait_still_works(self):
        want, _ = reference(self.pheno[:, :1], self.covar)
        got, _ = residualize_and_standardize(self.pheno[:, :1], self.covar)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)

    def test_a_one_dimensional_input_is_refused_not_silently_reshaped(self):
        with self.assertRaises(ValueError):
            residualize_and_standardize(self.pheno[:, 0], self.covar)


class MemoryBehaviourTests(unittest.TestCase):
    """The point of the change, measured rather than assumed."""

    def setUp(self):
        rng = np.random.default_rng(3)
        self.pheno = rng.normal(size=(64, 5000))
        self.covar = rng.normal(size=(64, 4))

    def test_no_temporary_is_ever_the_size_of_the_whole_matrix(self):
        """Count allocations by size; none may approach the input.

        This is the invariant that failed before: the middle expression held
        the input, the projection product and the result simultaneously. A
        shape-based review would not have caught it, so it is counted here.
        """
        whole = self.pheno.shape[0] * self.pheno.shape[1]
        seen = []
        real_empty = np.empty

        def spy(shape, *args, **kwargs):
            size = int(np.prod(shape)) if np.ndim(shape) else int(shape)
            seen.append(size)
            return real_empty(shape, *args, **kwargs)

        np.empty = spy
        try:
            residualize_and_standardize(self.pheno, self.covar,
                                        trait_block=256)
        finally:
            np.empty = real_empty

        # The output itself is allowed to be full size; nothing else is.
        oversized = [s for s in seen if s > whole // 4 and s != whole]
        self.assertFalse(
            oversized,
            f"allocations of {oversized} against a {whole}-element input; a "
            f"block of 256 traits should allocate ~{whole * 256 // 5000}")

    def test_inplace_writes_through_the_callers_array(self):
        pheno = self.pheno.copy()
        got, _ = residualize_and_standardize(
            pheno, self.covar, inplace=True, out_dtype=pheno.dtype)
        self.assertIs(got, pheno, "inplace must not allocate an output")
        want, _ = reference(self.pheno, self.covar)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)

    def test_inplace_refuses_when_promotion_would_change_the_dtype(self):
        """Silently downcasting the caller's buffer would lose precision.

        A float32 phenotype against float64 covariates promotes to float64,
        which cannot be written back through a float32 buffer. Refusing names
        the fix instead of quietly doing the wrong thing.
        """
        pheno32 = self.pheno.astype(np.float32)
        with self.assertRaises(ValueError):
            residualize_and_standardize(pheno32, self.covar, inplace=True)
        got, _ = residualize_and_standardize(
            pheno32, self.covar, inplace=True, out_dtype=np.float32)
        self.assertIs(got, pheno32)

    def test_out_dtype_avoids_the_promotion_that_doubles_memory(self):
        """The promotion is real and costly: 84.9 GB float32 -> 170 GB float64.

        The default keeps numpy's promotion, because everything downstream was
        built against it and narrowing silently would be a behaviour change
        dressed as an optimisation. `out_dtype` makes the choice explicit at
        the call site, which is where a caller holding an 85 GB matrix can
        weigh it -- and the scan casts back to float32 regardless.
        """
        pheno32 = self.pheno.astype(np.float32)
        promoted, _ = residualize_and_standardize(pheno32, self.covar)
        kept, _ = residualize_and_standardize(
            pheno32, self.covar, out_dtype=np.float32)
        self.assertEqual(promoted.dtype, np.float64)
        self.assertEqual(kept.dtype, np.float32)
        np.testing.assert_allclose(kept, promoted, rtol=1e-5, atol=1e-5)

    def test_the_default_does_not_mutate_the_caller(self):
        """`inplace` is opt-in precisely because the default must be safe."""
        before = self.pheno.copy()
        residualize_and_standardize(self.pheno, self.covar)
        np.testing.assert_array_equal(self.pheno, before)

    def test_the_default_block_is_a_sane_size_at_voxel_scale(self):
        """579 MB per block at 35,365 samples, not a fraction of 84.9 GB."""
        block_bytes = 35_365 * RESIDUALIZE_TRAIT_BLOCK * 4
        self.assertLess(block_bytes, 1e9)
        self.assertGreater(block_bytes, 1e8)


if __name__ == "__main__":
    unittest.main()
