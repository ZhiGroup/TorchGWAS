"""Contention belongs to hardware, not to every number in the dict.

`apply_to_rates` divides every rate by one host-CPU factor, and justified
including the GPU's GEMM rate on the grounds that skipping it "would leave the
one term that matters most at high K uncontended". That is backwards: a busy
host cannot slow a GEMM that is executing on the GPU.

It has not shown up yet because the scan is decode-bound -- 78.97 GB at
20.01 GB/s is 3.95 s of a 4.68 s wall -- so contending everything and
contending only decode give nearly the same number, and the one validated
measurement (9.01 s quiet, 35.14 s at load 160) cannot separate them. At high
K the GEMM binds and the difference becomes the whole answer.
"""
import unittest

from torchgwas.contention import (RATE_RESOURCE, apply_by_resource,
                                  apply_to_rates, contention_factor)

RATES = {
    'disk_bytes_per_second': 5.56e9,
    'h2d_bytes_per_second': 55.26e9,
    'd2h_bytes_per_second': 26.0e9,
    'write_bytes_per_second': 4.86e9,
    'decode_bytes_per_second': 20.01e9,
    'gemm_flops_per_second': {36: 25.6e12, 156: 38.2e12, 1052: 29.7e12},
}


class ResourceAwareContentionTests(unittest.TestCase):
    def test_host_load_does_not_slow_the_gpu(self):
        """The defect, stated as a test.

        A GEMM running on the GPU is not slowed by 160 foreign threads on the
        host. `apply_to_rates` says it is, by a factor of 3.667.
        """
        factor = contention_factor(16, 160, 48)
        by_resource = apply_by_resource(RATES, host_cpu_factor=factor)
        everything = apply_to_rates(RATES, factor)

        self.assertEqual(by_resource['gemm_flops_per_second'],
                         RATES['gemm_flops_per_second'])
        self.assertLess(everything['gemm_flops_per_second'][156],
                        RATES['gemm_flops_per_second'][156] * 0.9)

    def test_host_load_does_slow_the_decode(self):
        """Decode is the reader pool's zstd work and is genuinely contended."""
        factor = contention_factor(16, 160, 48)
        got = apply_by_resource(RATES, host_cpu_factor=factor)
        self.assertAlmostEqual(
            got['decode_bytes_per_second'],
            RATES['decode_bytes_per_second'] / factor, places=3)

    def test_disk_contention_is_its_own_mechanism(self):
        """CPU load and a busy disk are different things with different factors.

        Conflating them was already recorded as an error elsewhere in this
        module; this pins that disk rates follow the disk factor and host rates
        do not.
        """
        got = apply_by_resource(RATES, host_cpu_factor=4.0, disk_factor=2.0)
        self.assertAlmostEqual(got['disk_bytes_per_second'],
                               RATES['disk_bytes_per_second'] / 2.0, places=3)
        self.assertAlmostEqual(got['write_bytes_per_second'],
                               RATES['write_bytes_per_second'] / 2.0, places=3)
        self.assertAlmostEqual(got['decode_bytes_per_second'],
                               RATES['decode_bytes_per_second'] / 4.0, places=3)

    def test_an_unrecognised_rate_is_left_alone_not_contended(self):
        """Refusing to guess. Contending a term whose hardware was never
        identified is exactly how the GEMM came to be divided by a CPU factor.
        """
        extra = dict(RATES, mystery_bytes_per_second=1.0e9)
        got = apply_by_resource(extra, host_cpu_factor=4.0)
        self.assertEqual(got['mystery_bytes_per_second'], 1.0e9)

    def test_an_idle_host_changes_nothing(self):
        got = apply_by_resource(RATES, host_cpu_factor=1.0, disk_factor=1.0)
        self.assertEqual(got, RATES)

    def test_every_rate_the_model_takes_has_a_declared_resource(self):
        """A rate with no entry is silently uncontended, so the list must be
        complete for the rates `explain_time` actually accepts."""
        for key in RATES:
            self.assertIn(key, RATE_RESOURCE,
                          f"{key} has no declared hardware, so it would be "
                          f"left uncontended by default without anyone saying so")

    def test_factors_must_be_positive(self):
        for bad in (0.0, -1.0):
            with self.assertRaises(ValueError):
                apply_by_resource(RATES, host_cpu_factor=bad)
            with self.assertRaises(ValueError):
                apply_by_resource(RATES, disk_factor=bad)

    def test_it_reproduces_the_measured_slowdown_where_decode_binds(self):
        """The validated case still validates, for the right reason now.

        9.01 s quiet against 35.14 s at load ~160 is 3.90x measured. The scan
        is decode-bound, so contending decode alone reproduces it; the old
        function got the same answer by also slowing four terms that were not
        binding and one that cannot be slowed at all.
        """
        factor = contention_factor(16, 160, 48)
        self.assertAlmostEqual(factor, 3.667, places=2)
        measured = 35.14 / 9.01
        self.assertLess(abs(factor - measured) / measured, 0.10)


if __name__ == "__main__":
    unittest.main()
