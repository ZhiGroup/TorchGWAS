"""The run's wall clock must add up, and a phase that never happened must be
absent rather than zero.

Written because the 600,000-trait stress run reported 4,183 s total and
2,379 s of scan-and-write, and there was no way to ask where the other
1,804 seconds went -- 43% of the run, invisible in the output it produced.
A gap that large should be nameable from `run.json` alone.
"""
import sys
import unittest

sys.path.insert(0, "src")

from torchgwas.api import _phase_breakdown


class PhaseBreakdownTestCase(unittest.TestCase):
    def test_phases_sum_to_the_total(self):
        phases = _phase_breakdown(100.0, 102.0, 110.0, 140.0, 200.0)
        self.assertAlmostEqual(phases["open_and_resolve"], 2.0)
        self.assertAlmostEqual(phases["input_qc"], 8.0)
        self.assertAlmostEqual(phases["scan_setup"], 30.0)
        self.assertAlmostEqual(phases["scan_and_write"], 60.0)
        self.assertAlmostEqual(phases["total"], 100.0)
        named = sum(v for k, v in phases.items() if k != "total")
        self.assertAlmostEqual(named, phases["total"])

    def test_the_stress_run_gap_is_named(self):
        # The real numbers. `scan_setup` is the phase that was invisible.
        phases = _phase_breakdown(0.0, 1.0, 20.0, 1804.3, 4183.3)
        self.assertAlmostEqual(phases["scan_and_write"], 2379.0, places=1)
        self.assertGreater(phases["scan_setup"], 1700.0)
        self.assertAlmostEqual(
            phases["scan_setup"] + phases["scan_and_write"]
            + phases["input_qc"] + phases["open_and_resolve"],
            phases["total"], places=6)

    def test_a_phase_that_never_ran_is_absent_not_zero(self):
        # An in-memory result never reaches the writer. Reporting zero seconds
        # of scan-and-write would read as "instant" rather than "not measured".
        phases = _phase_breakdown(0.0, 1.0, 5.0, None, 9.0)
        self.assertNotIn("scan_and_write", phases)
        self.assertNotIn("scan_setup", phases)
        self.assertIn("input_qc", phases)
        self.assertAlmostEqual(phases["total"], 9.0)

    def test_nothing_measured_still_reports_the_total(self):
        phases = _phase_breakdown(0.0, None, None, None, 3.0)
        self.assertEqual(set(phases), {"total"})
        self.assertAlmostEqual(phases["total"], 3.0)


if __name__ == "__main__":
    unittest.main()
