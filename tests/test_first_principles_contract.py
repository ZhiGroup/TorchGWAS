"""Independent invariants for source-derived service models."""
from dataclasses import replace
from pathlib import Path
import sys
import unittest
BENCHMARKS = Path(__file__).resolve().parents[1] / 'benchmarks'
sys.path.insert(0, str(BENCHMARKS))
# The competitor and scheduler models are benchmark scripts that are not in
# this repository; their tests skip, the pipeline-model tests still run.
try:
    import direct_parallel_scaling_model as scaling
except ImportError:
    scaling = None
try:
    import direct_fastgwa_cost_model as fast
except ImportError:
    fast = None
needs_scaling = unittest.skipIf(scaling is None, 'benchmarks/direct_parallel_scaling_model.py is not in this repository')
needs_fast = unittest.skipIf(fast is None, 'benchmarks/direct_fastgwa_cost_model.py is not in this repository')
needs_plink = unittest.skipUnless((BENCHMARKS / 'direct_plink2_cost_model.py').exists(), 'benchmarks/direct_plink2_cost_model.py is not in this repository')
from torchgwas.pipeline_model import Workload, InputProfile, Hardware, PipelinePlan, estimate

@needs_scaling
class SchedulerTests(unittest.TestCase):
    def test_no_unjustified_point_or_upper_bound_under_oversubscription(self):
        for foreign in (1, 100, 10000):
            r = scaling.predict_block_seconds(48, foreign_cores=foreign)
            self.assertIsNone(r['block_seconds'])
            self.assertIsNone(r['block_high'])
            self.assertIsNone(r['oversubscription_stall'])

    def test_uncontended_service_has_shared_format_and_parallel_work(self):
        r = scaling.predict_block_seconds(12, foreign_cores=0)
        self.assertAlmostEqual(r['block_seconds'], r['format_serial'] + r['region_overhead'] + r['compute_ideal'])

@needs_fast
class FastGwaTests(unittest.TestCase):
    def cohort(self):
        return fast.Cohort(10000,35365,4,250,200000)

    def test_missing_cycle_is_unknown_not_zero(self):
        m=fast.Machine(1e9,1e10,1e11,1e9,48)
        r=fast.fastgwa_seconds(self.cohort(),1,m)
        self.assertIsNone(r['prediction_seconds'])
        self.assertIn('decode',r['unresolved_terms'])
        self.assertGreater(r['service_seconds'],0)

    def test_traits_repeat_the_entire_invocation(self):
        a=fast.fastgwa_seconds(self.cohort(),1,fast.H100_QUIET_MIN_CYCLE)
        b=fast.fastgwa_seconds(self.cohort(),8,fast.H100_QUIET_MIN_CYCLE)
        self.assertAlmostEqual(b['service_seconds'],8*a['service_seconds'])

    def test_unknown_scheduler_cycle_propagates(self):
        r=fast.fastgwa_seconds(self.cohort(),1,fast.H100_QUIET)
        self.assertIsNone(r['prediction_seconds'])
        self.assertIn('block_cycle_at_48_threads',r['unresolved_terms'])

    def test_sample_shape_cannot_be_silently_extrapolated(self):
        with self.assertRaisesRegex(ValueError,'rates measured'):
            fast.fastgwa_seconds(replace(self.cohort(), samples=1000),1,fast.H100_QUIET)

    def test_chromosome_boundaries_and_partial_blocks(self):
        c=replace(self.cohort(),variants=2000,chromosome_variant_counts=(1000,1000))
        self.assertEqual(fast.block_sizes(c),[1000,1000])
        self.assertEqual(fast.block_sizes(replace(c,chromosome_variant_counts=(2000,))),[1024,976])
        with self.assertRaises(ValueError): fast.block_sizes(replace(c,chromosome_variant_counts=(1000,)))

    def test_invalid_trait_count_is_rejected(self):
        for k in (-1,0,0.5,True):
            with self.assertRaises(ValueError): fast.passes(k)

class ResidentGpuTests(unittest.TestCase):
    def setUp(self):
        self.w=Workload(10000,1000,8,2)
        self.p=InputProfile('pgen',2500000,1000,1,
            gpu_statistics_seconds_per_variant=1e-4,
            gpu_statistics_component_shape=(1000,8,2,1000),
            gpu_statistics_component_source='independent-component.json',
            gpu_statistics_component_kernel='torch')
        self.h=Hardware(1e9,1e10,1e11,1e12,1e11,2**30,2**30,8,1e10)
        self.plan=PipelinePlan(1000,1000,1000,2,2)

    def test_full_resident_service_is_not_repriced_as_flops(self):
        r=estimate(self.w,self.p,self.h,self.plan)
        self.assertAlmostEqual(r['gpu_compute_seconds'],1)
        faster_flops=estimate(self.w,self.p,replace(self.h,gpu_flops_per_second=1e15),self.plan)
        self.assertEqual(r['gpu_compute_seconds'],faster_flops['gpu_compute_seconds'])

    def test_wrong_shape_and_extra_launch_cost_are_rejected(self):
        with self.assertRaises(ValueError): estimate(replace(self.w,traits=16),self.p,self.h,self.plan)
        with self.assertRaises(ValueError): estimate(self.w,self.p,replace(self.h,compute_launch_seconds=1e-6),self.plan)

    def test_only_single_gpu_is_supported(self):
        with self.assertRaisesRegex(ValueError,'one GPU'):
            estimate(self.w,self.p,replace(self.h,device_count=2),self.plan)


class PortabilityTests(unittest.TestCase):
    def context(self):
        return dict(machine='box-a', processor='cpu/topology', accelerator='gpu',
                    software='tool/blas/cuda/compiler', placement='affinity/numa/threads')

    def test_changed_execution_context_requires_recalibration(self):
        from torchgwas.model_provenance import check_calibration_context
        c = self.context()
        for key in c:
            with self.assertRaisesRegex(ValueError, 'recalibrate'):
                check_calibration_context(c, {**c, key: 'different'})
        self.assertEqual(check_calibration_context(c,c,strict=True), [])

    def test_legacy_profiles_are_unverified_and_strict_mode_refuses_them(self):
        from torchgwas.model_provenance import check_calibration_context
        self.assertEqual(len(check_calibration_context(None,None)),2)
        with self.assertRaises(ValueError):
            check_calibration_context(None,None,strict=True)
        with self.assertRaises(ValueError):
            check_calibration_context({},self.context())

    def test_gpu_rejects_rate_from_another_box(self):
        case=ResidentGpuTests(); case.setUp()
        c=self.context()
        with self.assertRaisesRegex(ValueError,'recalibrate'):
            estimate(case.w,replace(case.p,calibration_context=c),
                     replace(case.h,calibration_context={**c,'machine':'box-b'}),case.plan)

    @needs_fast
    def test_fastgwa_does_not_inherit_a_host_memory_rate(self):
        m=fast.Machine(1e9,1e10,1e11,1e9,48)
        c=fast.Cohort(10000,35365,4,250,200000)
        result=fast.fastgwa_seconds(c,1,m)
        self.assertIsNone(result['seconds_per_pass']['main.compute_memory_bound(reference)'])
        with self.assertRaises(ValueError):
            fast.fastgwa_seconds(c,1,replace(m,strict_calibration=True))


class ExecutionComponentTests(unittest.TestCase):
    def setUp(self):
        from torchgwas.pipeline_model import ExecutionProfile
        self.w=Workload(10000,1000,8,2)
        context=PortabilityTests().context()
        self.p=InputProfile('pgen',2500000,1000,1,calibration_context=context)
        self.h=Hardware(1e9,1e10,1e11,1e12,1e11,2**30,2**30,8,1e10,calibration_context=context)
        self.plan=PipelinePlan(1000,1000,1000,2,2)
        self.e=ExecutionProfile(2,1e-5,.1,(1000,8,2,1000,2,2),'prefilled.json',context)

    def test_component_replaces_gpu_roofline_and_adds_setup_once(self):
        r=estimate(self.w,self.p,self.h,self.plan,execution=self.e)
        self.assertAlmostEqual(r['execution_setup_seconds'],2.1)
        self.assertAlmostEqual(r['execution_consumer_seconds'],1)
        self.assertAlmostEqual(r['prediction_seconds'],2.1+r['execution_scan_seconds'])
        slow=estimate(self.w,self.p,replace(self.h,gpu_flops_per_second=1),self.plan,execution=self.e)
        self.assertEqual(r['prediction_seconds'],slow['prediction_seconds'])

    def test_metadata_and_consumer_scale_but_startup_does_not(self):
        r=estimate(replace(self.w,variants=20000),self.p,self.h,self.plan,execution=self.e)
        self.assertAlmostEqual(r['execution_setup_seconds'],2.2)
        self.assertAlmostEqual(r['execution_consumer_seconds'],2)

    def test_wrong_geometry_context_and_invalid_rate_rejected(self):
        for e in (replace(self.e,shape=(1000,8,2,2000,2,2)),
                  replace(self.e,consumer_seconds_per_chunk=float('nan')),
                  replace(self.e,calibration_context={**self.e.calibration_context,'machine':'another'})):
            with self.assertRaises(ValueError):estimate(self.w,self.p,self.h,self.plan,execution=e)

    def test_measured_producer_replaces_ideal_worker_scaling(self):
        e=replace(self.e,producer_seconds_per_chunk=.2,first_use_seconds=.3)
        r=estimate(self.w,self.p,self.h,self.plan,execution=e)
        slow=estimate(self.w,replace(self.p,cpu_decode_core_seconds_per_variant=1),self.h,self.plan,execution=e)
        self.assertEqual(r['prediction_seconds'],slow['prediction_seconds'])
        self.assertAlmostEqual(r['execution_producer_seconds'],2)
        no_init=estimate(self.w,self.p,self.h,self.plan,execution=replace(e,first_use_seconds=0))
        self.assertAlmostEqual(r['prediction_seconds']-no_init['prediction_seconds'],.3)

    def test_setup_sensitivity_is_observed_range(self):
        r=estimate(self.w,self.p,self.h,self.plan,execution=replace(self.e,setup_low_seconds=1,setup_high_seconds=3))
        self.assertAlmostEqual(r['setup_sensitivity_seconds'][1]-r['setup_sensitivity_seconds'][0],2)


@needs_plink
class FiniteCpuScheduleTests(unittest.TestCase):
    def test_one_plink_block_must_drain_all_three_phases(self):
        import direct_plink2_cost_model as plink
        self.assertEqual(plink.finite_main_calc_schedule(1,4,2,100),7)

    def test_plink_format_overlaps_next_compute_but_not_next_read(self):
        import direct_plink2_cost_model as plink
        # r0=1,c0 ends5; r1 ends2, join5, c1 ends9,f0 ends7;
        # r2 ends8, join9,c2 ends13,f1 ends11, final f2 ends15.
        self.assertEqual(plink.finite_main_calc_schedule(3,12,6,3,1),15)
        self.assertEqual(plink.finite_main_calc_schedule(3,0,6,3,1),9)

    def test_partial_plink_block_and_invalid_service(self):
        import direct_plink2_cost_model as plink
        self.assertGreater(plink.finite_main_calc_schedule(1,4,2,100,60),4)
        with self.assertRaises(ValueError):plink.finite_main_calc_schedule(-1,4,2,100)


@needs_fast
class FastInitializationTests(unittest.TestCase):
    def test_preloop_component_replaces_existing_fixed_terms(self):
        c=fast.Cohort(1000,100,2,500,10000,record_mix={0:1000},psam_bytes=800,
            pheno_bytes=1500,qcovar_bytes=1800,chromosome_variant_counts=(1000,))
        m=fast.Machine(1e9,1e10,1e11,1e9,8,threads=4,process_startup_seconds=9,
            initialization_seconds_excluding_pvar=1,pvar_read_parse_seconds_per_byte=1e-6,
            initialization_shape=(100,2,800,1500,1800),initialization_source='preloop-phases.json')
        result=fast.fastgwa_seconds(c,1,m)
        self.assertAlmostEqual(result['fixed_serial_seconds_per_pass'],1.01)
        with self.assertRaises(ValueError):fast.fastgwa_seconds(replace(c,covariates=3),1,m)
        with self.assertRaises(ValueError):fast.fastgwa_seconds(c,1,replace(m,pvar_read_parse_seconds_per_byte=-1))
