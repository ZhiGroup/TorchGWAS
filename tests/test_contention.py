"""The contention rule must reproduce the measurements it was derived from.

These are not invented cases. Every number below was measured on this cohort
and is recorded in the session log, so a change that breaks one of them has
broken agreement with something real rather than with a unit test's opinion.
"""
import sys
import unittest

sys.path.insert(0, "src")

from torchgwas.contention import (apply_to_rates, contention_factor,
                                  disk_contention_factor, effective_cores,
                                  runnable_threads_from_measurement)


class EffectiveCoresTestCase(unittest.TestCase):
    def test_idle_host_gives_every_thread_a_core(self):
        self.assertEqual(effective_cores(16, 0, 48), 16)
        self.assertEqual(effective_cores(48, 0, 48), 48)

    def test_undersubscribed_host_is_not_shared(self):
        # 16 of ours plus 10 foreign is 26 on 48 cores: nobody waits.
        self.assertEqual(effective_cores(16, 10, 48), 16)

    def test_oversubscribed_host_shares_proportionally(self):
        # Our reader pool at the load the benchmarks actually ran under.
        self.assertAlmostEqual(effective_cores(16, 160, 48), 48 * 16 / 176)

    def test_never_exceeds_the_machine(self):
        self.assertLessEqual(effective_cores(1000, 0, 48), 1000)
        self.assertLessEqual(effective_cores(48, 10_000, 48), 48)

    def test_zero_threads_get_nothing(self):
        self.assertEqual(effective_cores(0, 100, 48), 0.0)

    def test_negative_load_is_clamped_not_rewarded(self):
        self.assertEqual(effective_cores(16, -5, 48), 16)


class ContentionFactorTestCase(unittest.TestCase):
    def test_our_measured_slowdown_is_predicted_within_ten_percent(self):
        # Our K=1 binary scan: 9.01 s quiet, 35.14 s at load ~160, so 3.90x.
        # 16 reader workers, 48 cores. Nothing here is fitted.
        predicted = contention_factor(16, 160, 48)
        self.assertAlmostEqual(predicted, 3.667, places=2)
        measured = 35.14 / 9.01
        self.assertLess(abs(predicted - measured) / measured, 0.10)

    def test_idle_host_has_no_penalty(self):
        self.assertEqual(contention_factor(16, 0, 48), 1.0)

    def test_factor_grows_with_load(self):
        light = contention_factor(16, 50, 48)
        heavy = contention_factor(16, 300, 48)
        self.assertGreater(heavy, light)


class RunnableThreadsTestCase(unittest.TestCase):
    def test_recovers_plink2s_effective_thread_count(self):
        # plink2 asked for 48 threads and got 10.2 cores at load ~128. The rule
        # says that is ~34.5 continuously runnable threads, not 48 -- which is
        # a fact about plink2's sync points, not about the scheduler.
        threads = runnable_threads_from_measurement(10.2, 128, 48)
        self.assertGreater(threads, 30)
        self.assertLess(threads, 40)

    def test_round_trips_through_effective_cores(self):
        for observed, load in ((10.2, 128), (2.7, 150), (9.5, 140)):
            threads = runnable_threads_from_measurement(observed, load, 48)
            self.assertAlmostEqual(effective_cores(threads, load, 48),
                                   observed, places=6)

    def test_small_work_shows_the_tool_limiting_not_the_scheduler(self):
        # plink2 at N=200 got 2.7 cores where the scheduler alone would have
        # allowed ~13. The inverted thread count is therefore far below 48, and
        # that gap is the tool's own parallel limit.
        threads = runnable_threads_from_measurement(2.7, 150, 48)
        self.assertLess(threads, 12)

    def test_saturated_process_clamps_to_the_machine(self):
        self.assertEqual(runnable_threads_from_measurement(48, 10, 48), 48)
        self.assertEqual(runnable_threads_from_measurement(60, 10, 48), 48)


class DiskContentionTestCase(unittest.TestCase):
    """Disk is a SEPARATE parameter from CPU, with different physics."""

    def test_a_measured_probe_beats_any_model(self):
        # The disk probe measured 9.00 GB/s on a quiet host and 6.4 GB/s
        # against contention; that ratio IS the factor, no structure needed.
        self.assertAlmostEqual(
            disk_contention_factor(16, 0, achieved_bytes_per_second=6.4e9,
                                   quiet_bytes_per_second=9.0e9),
            9.0 / 6.4, places=6)

    def test_falls_back_to_queued_streams_when_no_probe(self):
        self.assertAlmostEqual(disk_contention_factor(16, 16), 2.0)
        self.assertAlmostEqual(disk_contention_factor(16, 0), 1.0)

    def test_blocked_count_not_load_average_is_the_signal(self):
        # A host whose load is entirely D-state has an idle CPU and a saturated
        # disk. The CPU rule and the disk rule must disagree there, which is
        # the whole reason they are two parameters rather than one.
        cpu = contention_factor(16, 0, 48)        # nothing RUNNABLE competing
        disk = disk_contention_factor(16, 120)    # 120 blocked on I/O
        self.assertEqual(cpu, 1.0)
        self.assertGreater(disk, 8.0)

    def test_rejects_a_nonsense_probe(self):
        with self.assertRaises(ValueError):
            disk_contention_factor(16, 0, achieved_bytes_per_second=0.0,
                                   quiet_bytes_per_second=9.0e9)


class ApplyToRatesTestCase(unittest.TestCase):
    def test_scales_scalars_and_curves_alike(self):
        rates = {"disk_bytes_per_second": 5.56e9,
                 "gemm_flops_per_second": {36: 8.35e12, 156: 20.57e12}}
        got = apply_to_rates(rates, 2.0)
        self.assertAlmostEqual(got["disk_bytes_per_second"], 2.78e9)
        # The width curve must be scaled too: skipping it would leave the term
        # that dominates at high K uncontended.
        self.assertAlmostEqual(got["gemm_flops_per_second"][36], 4.175e12)
        self.assertAlmostEqual(got["gemm_flops_per_second"][156], 10.285e12)

    def test_leaves_non_numeric_entries_alone(self):
        got = apply_to_rates({"label": "h100", "rate": 10.0}, 2.0)
        self.assertEqual(got["label"], "h100")
        self.assertEqual(got["rate"], 5.0)

    def test_refuses_a_nonsense_factor(self):
        with self.assertRaises(ValueError):
            apply_to_rates({"rate": 1.0}, 0.0)


if __name__ == "__main__":
    unittest.main()
