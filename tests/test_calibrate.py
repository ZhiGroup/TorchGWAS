"""The calibration must not report a rate it did not measure.

Every test here corresponds to a measurement that came back *plausible and
wrong*. That is the failure mode worth guarding: an obviously broken number
gets noticed, a merely optimistic one becomes a profile and then a bad plan.
"""
from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from torchgwas.calibrate import (host_contention, is_quiet, measure_disk_read,
                                 measure_text_parse, measure_write)


class DiskConcurrencyTests(unittest.TestCase):
    """One stream is the wrong number, and reporting it alone misleads.

    Measured on the H100 host: 1.48 GB/s at one reader against 9.00 GB/s at
    sixteen -- a 6x spread. A profile carrying only the single-stream rate
    under-prices the read term by that factor, and the model then names the
    wrong binding resource entirely.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "payload.bin"
        # Small but multi-block, so the worker split is exercised.
        self.path.write_bytes(os.urandom(8 << 20))

    def test_reports_every_concurrency_not_just_one(self):
        result = measure_disk_read(self.path, block_bytes=1 << 20, blocks=8,
                                   repeat=1, workers=(1, 2, 4))
        self.assertEqual(set(result["bytes_per_second_by_workers"]),
                         {"1", "2", "4"})
        self.assertIn("single_stream_bytes_per_second", result)
        # The headline rate is the best achievable, not the single stream --
        # the scan reads with many workers and should be priced that way.
        self.assertGreaterEqual(
            result["bytes_per_second"],
            result["single_stream_bytes_per_second"] - 1.0)

    def test_records_whether_the_cache_was_actually_dropped(self):
        """A warm read measures RAM; the profile must say which it got."""
        result = measure_disk_read(self.path, block_bytes=1 << 20, blocks=4,
                                   repeat=1, workers=(1,))
        self.assertIn("cache_dropped", result)
        self.assertIsInstance(result["cache_dropped"], bool)
        # Residency is verified, not assumed: `posix_fadvise` is advisory.
        self.assertIn("resident_fraction_before", result)

    def test_reports_the_span_it_actually_read(self):
        """The rate is meaningless without it.

        The same host reported 9.00 GB/s over a 1.5 GB span and ~6.1 GB/s over
        20 GB. The small span was flattered by caching, and that optimistic
        number propagated into a "the pipeline does not overlap, 1.45-1.84x
        available" conclusion that was purely an artifact of it. A consumer
        cannot tell those two measurements apart unless the span is reported.
        """
        result = measure_disk_read(self.path, block_bytes=1 << 20, blocks=4,
                                   repeat=1, workers=(1,))
        self.assertEqual(result["span_bytes"], result["bytes"])
        self.assertGreater(result["span_bytes"], 0)

    def test_never_claims_to_read_past_the_end_of_the_file(self):
        """A short file would otherwise report the nominal span at a huge rate."""
        result = measure_disk_read(self.path, block_bytes=1 << 20,
                                   blocks=10_000, repeat=1, workers=(1,))
        self.assertLessEqual(result["span_bytes"], self.path.stat().st_size)


class ContentionTests(unittest.TestCase):
    """A profile taken on a busy box must not pass itself off as the machine.

    Observed on one host within the same hour: HBM 2.52 TB/s quiet against
    1.05 TB/s while another user held ~35 of 48 cores, moving the computed
    roofline ridge from 17.8 to 43.1 FLOP/byte -- which changes which resource
    the model says binds. Contention is not noise around a true value; it is a
    different machine, so more repeats cannot fix it.
    """

    def test_reports_load_and_foreign_cpu(self):
        state = host_contention()
        for key in ("load_1min", "foreign_cpu_percent", "cpu_count", "quiet",
                    "blocked_processes"):
            self.assertIn(key, state)
        self.assertIsInstance(state["quiet"], bool)

    def test_a_saturated_disk_disqualifies_even_when_the_cpu_looks_fine(self):
        """The case the CPU signals cannot see, and the one that matters most.

        Observed: load 17.04 and 429% foreign on 48 cores -- borderline by the
        CPU measures -- while twelve of another user's processes sat in D state
        with the filesystem pinned, enough that a test suite on the same
        storage ran at 4.8% CPU waiting on I/O. This module's headline product
        is a disk rate, so calibrating there would be exactly wrong.
        """
        self.assertFalse(is_quiet(load_1min=8.0, foreign_cpu_percent=50.0,
                                  cpu_count=48, blocked_processes=12))
        # And the same box with the disk free is fine.
        self.assertTrue(is_quiet(load_1min=8.0, foreign_cpu_percent=50.0,
                                 cpu_count=48, blocked_processes=0))

    def test_a_couple_of_blocked_processes_is_normal(self):
        """Some D-state processes are ordinary; only saturation matters."""
        self.assertTrue(is_quiet(load_1min=2.0, foreign_cpu_percent=5.0,
                                 cpu_count=48, blocked_processes=2))

    def test_unknown_io_state_does_not_veto(self):
        """`ps` is absent on some hosts; unknown must not mean contended."""
        self.assertTrue(is_quiet(load_1min=2.0, foreign_cpu_percent=5.0,
                                 cpu_count=48, blocked_processes=None))

    def test_the_box_this_session_ran_on_is_not_called_quiet(self):
        """The conditions every contaminated measurement was taken under.

        3,641% foreign CPU at load 65.67 on 48 cores -- the state the host sat
        in for over an hour while a concurrency sweep waited it out. A profile
        taken here must be marked untrustworthy, not published as the machine.
        """
        self.assertFalse(is_quiet(load_1min=65.67, foreign_cpu_percent=3641.0,
                                  cpu_count=48))

    def test_high_load_alone_is_enough_to_disqualify(self):
        """Load is a one-minute mean, so it catches what `ps` has already missed.

        A batch that has just finished leaves load high and foreign CPU near
        zero; its page-cache and thermal effects are still present.
        """
        self.assertFalse(is_quiet(load_1min=40.0, foreign_cpu_percent=0.0,
                                  cpu_count=48))

    def test_high_foreign_cpu_alone_is_enough_to_disqualify(self):
        """And `ps` catches the burst that load has not yet decayed to show.

        Three ripgreps starting this second read as ~300% CPU while the load
        average is still reporting the idle minute before them.
        """
        self.assertFalse(is_quiet(load_1min=0.5, foreign_cpu_percent=300.0,
                                  cpu_count=48))

    def test_an_idle_box_is_quiet(self):
        self.assertTrue(is_quiet(load_1min=0.4, foreign_cpu_percent=1.2,
                                 cpu_count=48))

    def test_the_load_threshold_scales_with_core_count(self):
        """Load 8 is idle on 48 cores and oversubscribed on 4."""
        self.assertTrue(is_quiet(load_1min=8.0, foreign_cpu_percent=0.0,
                                 cpu_count=48))
        self.assertFalse(is_quiet(load_1min=8.0, foreign_cpu_percent=0.0,
                                  cpu_count=4))

    def test_a_missing_signal_does_not_veto(self):
        """`/proc/loadavg` and `ps` are both absent on some hosts.

        Treating unknown as contended would mean never calibrating there.
        """
        self.assertTrue(is_quiet(load_1min=None, foreign_cpu_percent=None,
                                 cpu_count=None))
        self.assertTrue(is_quiet(load_1min=None, foreign_cpu_percent=2.0,
                                 cpu_count=48))
        # But a signal that IS present still vetoes.
        self.assertFalse(is_quiet(load_1min=None, foreign_cpu_percent=3641.0,
                                  cpu_count=48))


class TextParseTests(unittest.TestCase):
    """The metadata parse is a real term, and it is parsing, not reading.

    It sets `open`, which is flat in M and is 76% of an end-to-end run at
    M = 1,000,000. Timing a bare read instead would report storage bandwidth
    and under-price it by more than an order of magnitude.
    """

    def test_measures_a_rate_and_counts_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.tsv"
            row = "\t".join(str(i) for i in range(6)) + "\n"
            path.write_text(row * 20000)
            result = measure_text_parse(path, max_bytes=1 << 20)
        self.assertGreater(result["bytes_per_second"], 0.0)
        self.assertGreater(result["bytes"], 0)
        # Fields actually counted, i.e. the content was parsed rather than
        # merely streamed past.
        self.assertGreater(result["fields_seen"], 0)


class ZstdDecodeTests(unittest.TestCase):
    """Compression's price has to be measurable, or it cannot be modelled.

    The model predicted the zstd hard-call store at 7.63x against a measured
    1.84x, entirely because it had no decode term. Adding the term to
    `explain_time` is only half the fix -- the rate has to come from the host,
    like every other rate in the profile.
    """

    def test_reports_an_aggregate_rate_across_threads(self):
        from torchgwas.calibrate import measure_zstd_decode
        result = measure_zstd_decode(payload_bytes=1 << 20, threads=4,
                                     repeat=1)
        if not result.get("available"):
            self.skipTest("zstandard is not installed")
        self.assertGreater(result["bytes_per_second"], 0.0)
        self.assertEqual(result["threads"], 4)

    def test_the_default_payload_is_compressible(self):
        """Random bytes would be the pathological case, not a genotype proxy.

        They barely compress, so they also decompress almost trivially, and
        the measured rate would flatter any compressed store.
        """
        from torchgwas.calibrate import measure_zstd_decode
        result = measure_zstd_decode(payload_bytes=1 << 20, threads=2,
                                     repeat=1)
        if not result.get("available"):
            self.skipTest("zstandard is not installed")
        self.assertGreater(result["ratio"], 2.0)

    def test_real_bytes_can_be_supplied(self):
        """The rate depends on the content, so the caller should pass its own."""
        from torchgwas.calibrate import measure_zstd_decode
        payload = bytes([0] * (1 << 20))
        result = measure_zstd_decode(sample_bytes=payload, threads=2, repeat=1)
        if not result.get("available"):
            self.skipTest("zstandard is not installed")
        self.assertEqual(result["output_bytes"], 1 << 20)


class GemmCurveTests(unittest.TestCase):
    """The profile must be able to carry a curve, not just a scalar.

    `explain_time` accepts `{design_width: flops}`, but that is only half the
    fix: if calibration still produces one number, the model ships with a rate
    that is wrong at every width except the one measured. Measured spread on
    the H100 is 2.5x (8.35 TFLOPS at width 36 against 20.57 at 156).
    """

    def test_it_returns_a_width_keyed_table(self):
        from torchgwas.calibrate import measure_gemm_curve
        try:
            import torch
            if not torch.cuda.is_available():
                self.skipTest("no CUDA device")
        except ImportError:
            self.skipTest("torch is not installed")
        result = measure_gemm_curve(2048, covariates=7,
                                    trait_counts=(4, 16),
                                    chunk_variants=256, repeat=1)
        curve = result["gemm_flops_by_width"]
        # Widths are traits + covariates + 1, so the keys must be shaped by
        # the design and not by the trait count alone.
        self.assertEqual(sorted(curve), [12, 24])
        for width, rate in curve.items():
            self.assertGreater(rate, 0.0, f"width {width}")

    def test_the_table_plugs_straight_into_the_model(self):
        """The two halves have to fit together without a conversion step."""
        from torchgwas.pipeline_model import gemm_rate_at_width
        curve = {12: 1.0e12, 24: 2.0e12}
        self.assertEqual(gemm_rate_at_width(curve, 12), 1.0e12)
        self.assertEqual(gemm_rate_at_width(curve, 24), 2.0e12)
        self.assertAlmostEqual(gemm_rate_at_width(curve, 18), 1.5e12, delta=1e9)


class WriteTests(unittest.TestCase):
    def test_write_rate_includes_the_fsync(self):
        """Without the fsync this measures the page cache.

        A profile that prices the sumstats write at cache speed under-predicts
        every end-to-end number, and the write is the binding term at K = 512.
        """
        with tempfile.TemporaryDirectory() as tmp:
            result = measure_write(tmp, block_bytes=1 << 20, blocks=2)
        self.assertTrue(result["fsynced"])
        self.assertGreater(result["bytes_per_second"], 0.0)


if __name__ == "__main__":
    unittest.main()
