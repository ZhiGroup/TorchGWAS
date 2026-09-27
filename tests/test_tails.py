"""The Student-t tail must stay exact where scipy stops working.

`2 * scipy.special.stdtr(df, -|t|)` underflows to exactly zero from |t| = 38.354
at df = 20,000, and `upper_tail_log10` then clips, so every association past
that point reports the identical -log10 p = 307.6527. In a min-P scan the
ranking disappears exactly where the answer is.
"""

from __future__ import annotations

import unittest

import numpy as np
from scipy import special, stats

from torchgwas.tails import log_two_sided_t_sf, upper_tail_log10_from_t


class StudentTailTestCase(unittest.TestCase):
    def test_matches_scipy_wherever_scipy_has_not_underflowed(self):
        # The accuracy claim: agreement to near machine precision over the
        # whole range scipy can still represent, across five orders of
        # magnitude of df.
        for df in (2.0, 5.0, 30.0, 500.0, 20000.0, 499998.0):
            t = np.linspace(0.0, 30.0, 400)
            p = 2.0 * special.stdtr(df, -t)
            usable = p > 1e-290
            want = -np.log10(p[usable])
            got = upper_tail_log10_from_t(t[usable], df)
            with self.subTest(df=df):
                np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-9)

    def test_it_keeps_going_where_scipy_returns_zero(self):
        df = 20000.0
        # scipy's cliff.
        self.assertEqual(2.0 * special.stdtr(df, -39.0), 0.0)
        values = [float(upper_tail_log10_from_t(t, df))
                  for t in (39.0, 50.0, 100.0, 300.0, 1000.0)]
        for value in values:
            self.assertTrue(np.isfinite(value))
        # Strictly increasing, so the ranking a min-P scan needs survives.
        self.assertEqual(values, sorted(values))
        self.assertGreater(values[0], 307.6527, "must exceed the old clip")
        # Sanity against the asymptotic normal bound: the t tail is heavier
        # than the normal, so -log10 p must be smaller than the normal's.
        for t, value in zip((50.0, 100.0, 300.0), values[1:4]):
            normal = -(special.log_ndtr(-t) + np.log(2.0)) / np.log(10.0)
            self.assertLess(value, normal)

    def test_no_cliff_and_no_flat_region_across_the_old_boundary(self):
        df = 20000.0
        t = np.linspace(30.0, 200.0, 4000)
        values = upper_tail_log10_from_t(t, df)
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertTrue(np.all(np.diff(values) > 0),
                        "a flat region is the defect this exists to remove")

    def test_degrees_of_freedom_broadcast_per_variant(self):
        t = np.array([[10.0, 40.0], [80.0, 5.0]])
        df = np.array([1000.0, 50000.0])[:, None]
        got = upper_tail_log10_from_t(t, df)
        self.assertEqual(got.shape, (2, 2))
        # Each cell must equal the scalar-df computation for its own row.
        for row in range(2):
            np.testing.assert_allclose(
                got[row], upper_tail_log10_from_t(t[row], float(df[row, 0])),
                rtol=1e-12)

    def test_the_zero_statistic_is_a_tail_of_one(self):
        for df in (3.0, 100.0, 20000.0):
            with self.subTest(df=df):
                self.assertAlmostEqual(float(log_two_sided_t_sf(0.0, df)), 0.0,
                                       places=12)

    def test_it_is_symmetric_in_the_sign_of_t(self):
        df = 750.0
        t = np.array([0.5, 3.0, 25.0, 60.0])
        np.testing.assert_allclose(upper_tail_log10_from_t(t, df),
                                   upper_tail_log10_from_t(-t, df), rtol=1e-14)

    def test_small_degrees_of_freedom_still_agree(self):
        # df = 1 is Cauchy, the heaviest case and the one most likely to expose
        # a continued fraction that was tuned for large df.
        t = np.linspace(0.1, 50.0, 200)
        for df in (1.0, 2.0, 3.0):
            want = -np.log10(2.0 * stats.t.sf(t, df))
            got = upper_tail_log10_from_t(t, df)
            with self.subTest(df=df):
                np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-9)


class DeviceStudentTailTestCase(unittest.TestCase):
    """The device implementation must agree with the host one, bit for bit enough.

    It exists because the unreduced scan evaluates a tail per (variant, trait)
    cell, and on 2M values the torch path measured 42.8 ms on an A100 against
    987 ms for `scipy.special.stdtr` -- which is 23x faster *and*, unlike stdtr,
    does not return zero past |t| = 38.
    """

    def _both(self, t, df):
        import torch

        from torchgwas.tails import upper_tail_log10_from_t_torch

        host = upper_tail_log10_from_t(t, df)
        device = upper_tail_log10_from_t_torch(
            torch.as_tensor(t, dtype=torch.float64),
            torch.as_tensor(df, dtype=torch.float64)).numpy()
        return host, device

    def test_device_matches_host_where_scipy_still_works(self):
        for df in (1.0, 2.0, 5.0, 30.0, 500.0, 20000.0):
            t = np.linspace(0.0, 30.0, 200)
            host, device = self._both(t, df)
            with self.subTest(df=df):
                np.testing.assert_allclose(device, host, rtol=1e-8, atol=1e-8)

    def test_device_matches_host_far_past_the_cliff(self):
        df = 20000.0
        t = np.array([39.0, 50.0, 100.0, 300.0, 1000.0])
        self.assertEqual(2.0 * special.stdtr(df, -39.0), 0.0)
        host, device = self._both(t, df)
        np.testing.assert_allclose(device, host, rtol=1e-12)
        self.assertTrue(np.all(np.isfinite(device)))
        self.assertTrue(np.all(np.diff(device) > 0))

    def test_device_broadcasts_degrees_of_freedom_per_variant(self):
        t = np.array([[10.0, 40.0], [80.0, 5.0]])
        df = np.array([[1000.0], [50000.0]])
        host, device = self._both(t, df)
        self.assertEqual(device.shape, (2, 2))
        np.testing.assert_allclose(device, host, rtol=1e-10)

    def test_forty_iterations_is_enough(self):
        # The iteration count is fixed on the device (no early exit), so this
        # pins the claim that 40 suffices rather than leaving it a comment.
        from torchgwas.tails import _betacf_torch, _TORCH_ITERATIONS

        import torch

        self.assertEqual(_TORCH_ITERATIONS, 40)
        a = torch.full((64,), 10000.0, dtype=torch.float64)
        x = torch.linspace(1e-4, 0.9, 64, dtype=torch.float64)
        few = _betacf_torch(a, 0.5, x, iterations=40)
        many = _betacf_torch(a, 0.5, x, iterations=300)
        np.testing.assert_allclose(few.numpy(), many.numpy(), rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
