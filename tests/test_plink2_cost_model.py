"""The plink2 cost model's two former constants, now derived, against the source
and the recorded observations.

`Machine.cores_used = 13.6` and `DEFAULT_DENSE_PATH_FRACTION = 0.42` are gone.
The first is derived from plink2's thread-pool sizing (plink2_glm_linear.cc:
4097-4100, commit ca0f464) through the CFS proportional-share rule; the second
is computed from the .pgen by `direct_plink2_sparse_predicate.py`. The tests
below pin the source transcription, check the derived share against every
`Percent of CPU` observation that was recorded with its load, and read the
real cohort's census JSON so the number the model uses arrives with its
provenance.

Numbers in the observation tests are measurements from this session's bench
JSONs (`direct_plink2_cores_validation.py` prints the full table); a change
that breaks one of them has broken agreement with something real.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load():
    path = ROOT / "benchmarks" / "direct_plink2_cost_model.py"
    if not path.exists():
        raise unittest.SkipTest(f"benchmarks/{path.name} is not in this repository")
    spec =importlib.util.spec_from_file_location("plink2_cost_model_under_test", path)
    module = importlib.util.module_from_spec(spec)
    # Register before executing: @dataclass resolves annotations through
    # sys.modules[cls.__module__].
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ThreadPoolTests(unittest.TestCase):
    """plink2_glm_linear.cc:4097-4100 and the kMaxThreads clip, transcribed."""

    def setUp(self):
        self.mod = _load()

    def test_default_thread_count_is_the_cpu_count_clipped_at_kmaxthreads(self):
        # include/plink2_thread.cc:67-95 NumCpu; plink2_thread.h:96 kMaxThreads=512
        self.assertEqual(self.mod.plink2_default_thread_ct(48), 48)
        self.assertEqual(self.mod.plink2_default_thread_ct(600), 512)
        with self.assertRaises(ValueError):
            self.mod.plink2_default_thread_ct(0)

    def test_one_thread_is_kept_for_the_main_thread_above_eight(self):
        # calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct
        self.assertEqual(self.mod.plink2_calc_thread_ct(48), 47)
        self.assertEqual(self.mod.plink2_calc_thread_ct(9), 8)
        self.assertEqual(self.mod.plink2_calc_thread_ct(8), 8)
        self.assertEqual(self.mod.plink2_calc_thread_ct(1), 1)

    def test_calc_threads_never_exceed_the_variant_count(self):
        # if (calc_thread_ct > variant_ct) calc_thread_ct = variant_ct;
        self.assertEqual(self.mod.plink2_calc_thread_ct(48, variant_ct=10), 10)
        self.assertEqual(self.mod.plink2_calc_thread_ct(48, variant_ct=10_000), 47)

    def test_threads_flag_above_kmaxthreads_is_reduced(self):
        # plink2.cc:12906-12908
        self.assertEqual(self.mod.plink2_calc_thread_ct(1000), 511)


class DerivedCoresTests(unittest.TestCase):
    def setUp(self):
        self.mod = _load()

    def test_quiet_host_gives_the_whole_calc_pool(self):
        # 47 threads + 0 foreign on 48 cores: nobody waits.
        self.assertEqual(self.mod.plink2_cores_used(0.0, 48), 47.0)

    def test_loaded_host_shares_proportionally_with_the_source_thread_count(self):
        # C*T/(L+T) with T=47, not the flag's 48: 48*47/175.
        self.assertAlmostEqual(self.mod.plink2_cores_used(128.0, 48), 48 * 47 / 175, places=9)
        # The flag value would have given 13.09; the source count gives 12.89.
        self.assertLess(self.mod.plink2_cores_used(128.0, 48), 48 * 48 / 176)

    def test_explicit_threads_and_variant_cap_flow_through(self):
        self.assertEqual(self.mod.plink2_cores_used(0.0, 48, max_thread_ct=8), 8.0)
        self.assertEqual(self.mod.plink2_cores_used(0.0, 48, max_thread_ct=48, variant_ct=5), 5.0)

    def test_derived_share_bounds_every_recorded_observation_from_above(self):
        """Observed `Percent of CPU` / 100 next to the foreign runnable count.

        Sources: plink2_fullm/cell.json (busy 71.5% -> 34.3), tie_point/
        M1000000_N35365.json (saturated, loadavg 140.8), plink2_k_sweep/
        k1.json (loadavg 127.8), tie_real/M1000000_N35365.json (busy 1.4%
        -> 0.7), the paper-table orphan (1131% at load 125.1), and the quiet
        K=128 sampling (40.9). Only observations whose load was RECORDED are
        here; the K-sweep orphan (1226%) has its own test below because its
        load is known only as a range.
        """
        observed = [
            (22.82, 0.715 * 48),
            (10.85, 140.8), (10.58, 140.8),
            (11.00, 127.8), (10.16, 127.8),
            (41.12, 0.0144 * 48), (40.61, 0.0144 * 48),
            (11.31, 125.1),
            (40.9, 0.0),
        ]
        for cores, foreign in observed:
            derived = self.mod.plink2_cores_used(foreign, 48)
            self.assertGreaterEqual(derived, cores, msg=f"observed {cores} at L={foreign}")

    def test_the_k_sweep_orphan_is_consistent_within_its_load_range(self):
        # 1226% was read off `ps` after the K sweep was cancelled; the load at
        # that moment was not recorded, only the sweep's cells' host_before
        # load averages, 127.8-160.3. The derived share across that range
        # brackets the observation. At the midpoint it would read 1.04, which
        # is why this is not in the upper-bound test above: L is undetermined
        # to within the range, not the rule.
        low = self.mod.plink2_cores_used(160.3, 48)
        high = self.mod.plink2_cores_used(127.8, 48)
        self.assertLess(low, 12.26)
        self.assertGreater(high, 12.26)

    def test_the_shortfall_below_the_share_is_load_independent(self):
        """Quiet and loaded cells sit at the same fraction of the derived share.

        That is what makes it the tool's duty cycle rather than a scheduler
        error: a wrong L would move with load. The band here is the recorded
        range on the calc-dominated cells (N=35,365, M >= 200,000).
        """
        cells = [(41.12, 0.0144 * 48), (40.98, 0.0167 * 48), (22.82, 0.715 * 48),
                 (10.85, 140.8), (10.76, 139.7), (11.00, 127.8)]
        ratios = [c / self.mod.plink2_cores_used(l, 48) for c, l in cells]
        self.assertGreater(min(ratios), 0.78)
        self.assertLess(max(ratios), 0.95)

    def test_retracted_constant_is_not_what_the_rule_gives_on_a_quiet_host(self):
        # 13.6 came from a host at load ~125 (48*47/172 = 13.1) and was carried
        # to a quiet host, where the same rule gives 47: the 2.69x the
        # re-taken whole-cohort cell exposed.
        self.assertAlmostEqual(self.mod.plink2_cores_used(125.1, 48), 13.1, delta=0.1)
        self.assertGreater(self.mod.plink2_cores_used(0.0, 48) / 13.6, 3.0)


class MachineTests(unittest.TestCase):
    def setUp(self):
        self.mod = _load()

    def _machine(self, **overrides):
        kwargs = dict(disk_bytes_per_second=6.0e9, memory_bytes_per_second=28.0e9,
                      flops_per_second=550.0e9, write_bytes_per_second=1.6e9,
                      cores=48, cores_used=47.0)
        kwargs.update(overrides)
        return self.mod.Machine(**kwargs)

    def test_cores_used_has_no_default(self):
        with self.assertRaises(TypeError):
            self.mod.Machine(disk_bytes_per_second=6.0e9, memory_bytes_per_second=28.0e9,
                             flops_per_second=550.0e9, write_bytes_per_second=1.6e9, cores=48)

    def test_no_attribute_carries_the_retracted_constants(self):
        from dataclasses import MISSING
        self.assertFalse(hasattr(self.mod, "DEFAULT_DENSE_PATH_FRACTION"))
        self.assertFalse(hasattr(self.mod, "dense_fraction_from_mix"))
        self.assertIs(self.mod.Machine.__dataclass_fields__["cores_used"].default, MISSING)
        self.assertIs(self.mod.Cohort.__dataclass_fields__["dense_path_fraction"].default, MISSING)

    def test_calc_term_scales_inversely_with_cores_used(self):
        """The property that survived the model's restructuring.

        This test used to read `seconds["memory"]` and call `plink2_seconds`
        without `rates`. Both are gone: profiling showed the per-variant steps
        (covariate fill, syrk, X'y) are the SAME thread's consecutive
        instructions, so the model sums them instead of taking
        `max(compute, memory)`, and the step rates are now a required argument
        rather than an aggregate constant. Updated rather than deleted, because
        the invariant it guards is unchanged and is the whole point of having
        derived `cores_used`: the calc pool's time must fall as the share rises.
        """
        cohort = self.mod.Cohort(variants=1_000_000, samples=35_365, covariates=27,
                                 stored_bytes_per_variant=4000.0, pvar_bytes=0.0,
                                 dense_path_fraction=0.9882, mean_sparse_carriers=7744.0,
                                 gram_path_fraction=0.923)
        rates = self._rates()
        at_47 = self.mod.plink2_seconds(
            cohort, 1, self._machine(cores_used=47.0), rates)["seconds"]
        at_13 = self.mod.plink2_seconds(
            cohort, 1, self._machine(cores_used=13.6), rates)["seconds"]
        key = "calc" if "calc" in at_47 else "compute"
        self.assertAlmostEqual(at_13[key] / at_47[key], 47.0 / 13.6, places=6)

    def _rates(self):
        """Measured Gram-path step rates, or skip if the harness JSON is absent.

        These come from `benchmarks/direct_plink2_gram_path_rates.c`, which
        transcribes the real path and times each step at concurrency. There is
        deliberately no default: a step rate measured at one thread count does
        not describe another (fill 0.73 -> 1.31 ns/element from 1 to 47).
        """
        import glob

        candidates = sorted(glob.glob(
            str(ROOT / "_scratch" / "results" / "gram_path_rates_t*.json")))
        if not candidates:
            self.skipTest("gram path rates JSON not present")
        return self.mod.load_gram_path_rates(candidates[-1])

    def test_fractional_cpu_share_is_not_rounded_up(self):
        cohort = self.mod.Cohort(1000, 35365, 27, 4000, 0, 1, gram_path_fraction=1)
        a = self.mod.plink2_seconds(cohort, 1, self._machine(cores_used=1), self._rates())
        b = self.mod.plink2_seconds(cohort, 1, self._machine(cores_used=0.5), self._rates())
        self.assertAlmostEqual(b['seconds']['calc'], 2*a['seconds']['calc'])

    def test_passes_turn_over_at_the_subbatch_size(self):
        self.assertEqual(self.mod.passes(1), 1)
        self.assertEqual(self.mod.passes(240), 1)
        self.assertEqual(self.mod.passes(241), 2)


class TraitAccountingTests(MachineTests):
    def prediction(self, k, **kw):
        cohort = self.mod.Cohort(variants=1000, samples=35365, covariates=27,
            stored_bytes_per_variant=4000, pvar_bytes=0,
            dense_path_fraction=1, gram_path_fraction=1)
        return self.mod.plink2_seconds(cohort, k, self._machine(), self._rates(), **kw)

    def test_each_trait_repeats_the_full_xty_product(self):
        a, b = self.prediction(1), self.prediction(8)
        self.assertAlmostEqual(b['calc_breakdown']['xty'], 8*a['calc_breakdown']['xty'])
        self.assertAlmostEqual(b['calc_breakdown']['syrk'], a['calc_breakdown']['syrk'])

    def test_partial_last_batch_counts_only_remaining_traits(self):
        one, full, tail = self.prediction(1), self.prediction(240), self.prediction(241)
        self.assertEqual(tail['subbatch_sizes'], [240, 1])
        for step in full['calc_breakdown']:
            self.assertAlmostEqual(tail['calc_breakdown'][step],
                full['calc_breakdown'][step] + one['calc_breakdown'][step], msg=step)

    def test_formatting_counts_each_test_once_across_passes(self):
        result = self.prediction(481, format_seconds_per_row=0.001)
        self.assertAlmostEqual(result['seconds']['format'], 481)

    def test_read_and_format_share_the_main_thread(self):
        result = self.prediction(8, format_seconds_per_row=0.001)
        self.assertAlmostEqual(result['seconds']['main_read_format'],
            result['seconds']['read'] + result['seconds']['format'])

    def test_missing_decoder_measurement_is_reported(self):
        result = self.prediction(1)
        self.assertIn('PGEN decode component rate', result['unmeasured_or_approximated_terms'])
        self.assertIsNone(result['prediction_seconds'])

    def test_invalid_trait_counts_are_rejected(self):
        for k in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                self.prediction(k)


class CensusProvenanceTests(unittest.TestCase):
    """The dense fraction the model uses is the real file's census, with its commit."""

    CENSUS = ROOT / "_scratch" / "results" / "sparse_census_full_hardcall.json"

    def setUp(self):
        self.mod = _load()
        if not self.CENSUS.exists():
            self.skipTest("whole-file census JSON not present")
        with open(self.CENSUS) as fh:
            self.census = json.load(fh)

    def test_census_is_for_the_cited_commit_and_the_benchmark_configuration(self):
        self.assertEqual(self.census["plink_ng_commit"], self.mod.PLINK_NG_COMMIT)
        self.assertEqual(self.census["covariate_ct"], 27)
        self.assertEqual(self.census["variants"], 8_931_083)
        self.assertEqual(self.census["samples"], 35_365)

    def test_computed_dense_fraction_and_its_disagreement_with_the_retracted_values(self):
        dense = self.mod.dense_fraction_from_census(self.census)
        self.assertAlmostEqual(dense, 0.988206, places=5)
        # 92.3% of variants carry a missing call from --hard-call-threshold 0.1;
        # that, not the record form, is what forces the dense path.
        self.assertAlmostEqual(self.census["variants_with_missing"] / self.census["variants"],
                               0.923, places=3)
        self.assertGreater(dense / 0.42, 2.3)      # the retracted 0.42
        self.assertGreater(dense / 0.077, 12.0)    # the retracted 0.077
        self.assertLessEqual(self.census["gram_fraction"], dense)

    def test_a_bare_number_is_refused(self):
        with self.assertRaises(KeyError):
            self.mod.dense_fraction_from_census({"dense_fraction": 0.42})


class SyrkRateTests(unittest.TestCase):
    """The THIRD asserted constant: the aggregate FLOP rate.

    Once `cores_used` was derived, `compute` became the binding term, and it
    was charged at an aggregate 550.67 GFLOP/s -- one number that did not move
    when the load changed or when the machine had a different core count,
    which is the defining property of something that is not a machine rate.
    It is replaced by a measured per-core `cblas_dsyrk` rate times the derived
    share (`benchmarks/direct_plink2_syrk_rate.c`, 15.23 GFLOP/s best,
    single-threaded, verified at 100% CPU, N=35,365 P=29, one triangle).

    Note this made AGREEMENT WORSE, from 1.85x to 2.35x below the measured
    854.62 s, and that is the correct trade: the model is a floor over work
    the source names, and the remaining ~490 s is unattributed plink2 work
    rather than something to absorb into a coefficient.
    """

    def setUp(self):
        # `benchmarks/` is not an importable package; the rest of this file
        # loads the module by path and so must this.
        self.mod = _load()

    def machine(self, cores_used, per_core=15.23e9):
        return self.mod.Machine(
            disk_bytes_per_second=5.56e9, memory_bytes_per_second=40e9,
            flops_per_second=550.67e9, write_bytes_per_second=4.86e9,
            cores=48, cores_used=cores_used,
            syrk_flops_per_second_per_core=per_core)

    def test_compute_rate_scales_with_the_cores_plink2_gets(self):
        quiet = self.machine(46.0).compute_flops_per_second
        loaded = self.machine(11.3).compute_flops_per_second
        self.assertAlmostEqual(quiet / loaded, 46.0 / 11.3, places=6)

    def test_machine_does_not_inherit_a_different_hosts_syrk_rate(self):
        self.assertIsNone(self.mod.Machine.__dataclass_fields__[
            "syrk_flops_per_second_per_core"].default)
        with self.assertRaises(ValueError):
            _ = self.machine(12, per_core=None).compute_flops_per_second

    def test_zero_and_fractional_core_shares_do_not_invent_capacity(self):
        self.assertEqual(self.machine(0.0).compute_flops_per_second, 0)
        self.assertAlmostEqual(self.machine(0.25).compute_flops_per_second, 0.25*15.23e9)


if __name__ == "__main__":
    unittest.main()
