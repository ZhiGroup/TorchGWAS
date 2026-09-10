from __future__ import annotations

import unittest

import numpy as np
from scipy import stats

from torchgwas.linear import linear_scan
from torchgwas.preprocess import residualize_and_standardize


class ExactLinearStatisticsTestCase(unittest.TestCase):
    def test_matches_explicit_ols_with_rank_deficient_covariates(self):
        rng = np.random.default_rng(20260910)
        n_samples = 96
        genotype = rng.integers(0, 3, size=(n_samples, 7)).astype(np.float64)
        covariate = rng.normal(size=n_samples)
        covariates = np.column_stack([covariate, 2.0 * covariate])
        phenotype = np.column_stack(
            [
                0.4 * genotype[:, 1] + 0.8 * covariate + rng.normal(size=n_samples),
                -0.2 * genotype[:, 4] - 0.5 * covariate + rng.normal(size=n_samples),
            ]
        )

        beta, t_stat, p_value, q_matrix = linear_scan(
            genotype,
            phenotype,
            covariates,
            chunk_size=3,
            device="cpu",
            compute_dtype="float64",
        )
        phenotype_processed, q_reference = residualize_and_standardize(phenotype, covariates)
        self.assertIsNotNone(q_matrix)
        self.assertEqual(q_matrix.shape[1], 1)
        np.testing.assert_allclose(np.abs(q_matrix), np.abs(q_reference), atol=1e-12)

        beta_reference = np.empty_like(beta)
        t_reference = np.empty_like(t_stat)
        p_reference = np.empty_like(p_value)
        for marker in range(genotype.shape[1]):
            design = np.column_stack([np.ones(n_samples), q_reference, genotype[:, marker]])
            coefficients, _, rank, _ = np.linalg.lstsq(design, phenotype_processed, rcond=None)
            residual = phenotype_processed - design @ coefficients
            df = n_samples - rank
            sigma2 = np.sum(residual * residual, axis=0) / df
            covariance_last = np.linalg.pinv(design.T @ design)[-1, -1]
            standard_error = np.sqrt(sigma2 * covariance_last)
            beta_reference[marker] = coefficients[-1]
            t_reference[marker] = coefficients[-1] / standard_error
            p_reference[marker] = 2.0 * stats.t.sf(np.abs(t_reference[marker]), df=df)

        np.testing.assert_allclose(beta, beta_reference, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(t_stat, t_reference, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(p_value, p_reference, rtol=1e-10, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
