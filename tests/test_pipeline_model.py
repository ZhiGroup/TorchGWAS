from dataclasses import replace, asdict
from types import SimpleNamespace
import numpy as np
import unittest
from unittest.mock import patch
from torchgwas.pipeline_model import (Workload, InputProfile, Hardware, PipelinePlan,
                                      estimate, choose_plan, plan_parameters,
                                      apply_pipeline_profile, host_pinned_bytes,
                                      explain_memory, device_ring_bytes,
                                      explain_time, auto_chunk_variants,
                                      MAX_AUTO_CHUNK_VARIANTS,
                                      VALIDATED_MAX_CHUNK_VARIANTS,
                                      predicted_peak_bytes, auto_trait_block)


class PipelineModelTests(unittest.TestCase):
    def setUp(self):
        resolver = patch('torchgwas.scan_gpu.resolve_statistics_backend', return_value='torch')
        resolver.start()
        self.addCleanup(resolver.stop)
        self.w = Workload(10000, 1000, 8, 2)
        self.p = InputProfile('pgen-hardcall-gpu', 2500000, 250, 4,
                              gpu_decode_seconds_per_variant=1e-5,
                              max_stored_bytes_per_variant=250, decode_on_gpu=True)
        self.h = Hardware(1e9, 1e10, 1e11, 1e12, 1e11, 2**30, 2**30, 8, d2h_bytes_per_second=1e10)
        self.plan = PipelinePlan(1000, 1000, 1000, 2, 2)

    def test_host_staging_copy_accounts_read_and_write(self):
        for gpu in (False, True):
            direct = replace(self.p, decode_on_gpu=gpu, host_staging_copies=0)
            staged = replace(direct, host_staging_copies=1)
            a = estimate(self.w, direct, self.h, self.plan)
            b = estimate(self.w, staged, self.h, self.plan)
            expected = 2 * self.w.variants * self.p.transfer_bytes_per_variant / self.h.host_bytes_per_second
            self.assertAlmostEqual(b['resource_seconds']['host_memory'] - a['resource_seconds']['host_memory'], expected)
            self.assertEqual(a['resource_seconds']['h2d'], b['resource_seconds']['h2d'])
        for invalid in (-1, 0.5, True):
            with self.assertRaises(ValueError):
                estimate(self.w, replace(self.p, host_staging_copies=invalid), self.h, self.plan)

    def test_profile_rejects_unimplemented_direct_to_pinned_claim(self):
        source, profile = self.profile_fixture()
        profile['input']['host_staging_copies'] = 0
        with self.assertRaisesRegex(ValueError, 'host_staging_copies'):
            apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
        self.assertEqual(source.decode_batch_size, 8)

    def test_default_depth_search_can_use_many_cpu_workers(self):
        p = replace(self.p, decode_on_gpu=False, gpu_decode_seconds_per_variant=0,
                    cpu_decode_core_seconds_per_variant=.01)
        h = replace(self.h, cpu_workers=48, host_memory_bytes=2**32, device_memory_bytes=2**32)
        candidates = dict(read_variants=(100,), decode_variants=(100,), chunk_variants=(100,))
        automatic = choose_plan(self.w, p, h, **candidates)
        constrained = choose_plan(self.w, p, h, depths=(2, 3, 4), **candidates)
        self.assertEqual(automatic['plan']['workers'], 48)
        self.assertEqual(automatic['plan']['depth'], 48)
        self.assertTrue(automatic['memory_feasible'])
        self.assertLess(automatic['planning_seconds'], constrained['planning_seconds'] / 5)
        self.assertLessEqual(constrained['plan']['depth'], 4)

    def test_direct_fill_infers_staging_and_ties_geometry(self):
        source, profile = self.profile_fixture()
        source.allows_direct_native_fill = True
        source.decode_tile_matches_chunk = True
        source.native_host_staging_copies = 0
        del profile['input']['host_staging_copies']
        profile.pop('source_controls')
        profile['candidates']['read_variants'] = [2000]
        result = apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
        self.assertEqual(result['estimate']['effective_geometry']['decode_variants'], 1000)
        self.assertEqual(result['estimate']['effective_geometry']['read_variants'], 1000)
        self.assertEqual(source.decode_batch_size, 8)
        self.assertEqual(result['model_constraints']['host_staging_copies'], 0)
        self.assertTrue(result['all_source_controls_applied'])
        base = replace(self.p, decode_on_gpu=False, direct_native_fill=True,
                       decode_tile_matches_chunk=True, host_staging_copies=0)
        direct = estimate(self.w, base, self.h, self.plan)
        reordered = estimate(self.w, replace(base, host_staging_copies=1), self.h, self.plan)
        self.assertEqual(reordered['host_buffer_bytes'] - direct['host_buffer_bytes'],
                         2 * 1000 * 1000 * 4)
        self.assertGreater(reordered['resource_seconds']['host_memory'],
                           direct['resource_seconds']['host_memory'])

    def test_default_workers_include_memory_efficient_intermediate(self):
        p = replace(self.p, decode_on_gpu=False, gpu_decode_seconds_per_variant=0,
                    cpu_decode_core_seconds_per_variant=.01)
        h = replace(self.h, cpu_workers=48, host_memory_bytes=2**32, device_memory_bytes=2**32)
        candidates = dict(read_variants=(100,), decode_variants=(100,),
                          chunk_variants=(100,), depths=(16,))
        chosen = choose_plan(self.w, p, h, **candidates)
        oversubscribed = choose_plan(self.w, p, h, workers=(48,), **candidates)
        self.assertEqual(chosen['plan']['workers'], 16)
        self.assertEqual(chosen['plan']['depth'], 16)
        self.assertEqual(chosen['planning_seconds'], oversubscribed['planning_seconds'])
        self.assertLessEqual(chosen['host_buffer_bytes'], oversubscribed['host_buffer_bytes'])
        self.assertTrue(chosen['memory_feasible'])

    def test_cpu_decoded_input_conversion_can_dominate_gpu(self):
        p = replace(self.p, format='pgen-hardcall', decode_on_gpu=False,
                    gpu_decode_seconds_per_variant=0)
        base = estimate(self.w, p, self.h, self.plan)
        converted = estimate(self.w, replace(p, gpu_input_conversion_seconds_per_variant=.01),
                             self.h, self.plan)
        self.assertEqual(converted['gpu_input_conversion_seconds'], 100)
        self.assertEqual(converted['gpu_memory_bytes'],
                         self.w.variants * (32*self.w.samples + 16*(self.w.traits+self.w.covariates+1)))
        self.assertAlmostEqual(converted['resource_seconds']['gpu'], base['resource_seconds']['gpu']+100)
        self.assertEqual(converted['bottleneck'], 'gpu')
        for resource in ('storage', 'cpu_decode', 'h2d', 'host_memory'):
            self.assertEqual(converted['resource_seconds'][resource], base['resource_seconds'][resource])

    def test_native_fused_statistics_traffic_and_actual_selection(self):
        workload = replace(self.w, statistics_kernel='native_fused')
        for width in (1, 4):
            result = estimate(workload, replace(self.p, native_input_bytes_per_value=width), self.h, self.plan)
            self.assertEqual(result['gpu_memory_bytes'], self.w.variants *
                             ((2*width+8)*self.w.samples+16*(self.w.traits+self.w.covariates+1)))
            self.assertTrue(result['conversion_fused_into_statistics'])
        with self.assertRaises(ValueError):
            estimate(workload, self.p, self.h, self.plan)
        with self.assertRaises(ValueError):
            estimate(workload, replace(self.p, native_input_bytes_per_value=1,
                                      gpu_input_conversion_seconds_per_variant=.01), self.h, self.plan)
        source, profile = self.profile_fixture()
        with patch('torchgwas.scan_gpu.resolve_statistics_backend', return_value='native_fused'):
            result = apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
            self.assertEqual(result['workload']['statistics_kernel'], 'native_fused')
            self.assertEqual(result['estimate']['gpu_input_conversion_status'], 'fused_into_statistics')
            profile['workload'] = {'statistics_kernel': 'torch'}
            with self.assertRaisesRegex(ValueError, 'statistics_kernel'):
                apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)

    def test_packed_pgen_uses_padded_physical_rows(self):
        source, profile = self.profile_fixture()
        source.native_encoding = 'pgen_2bit'
        source.native_transfer_dtype = np.uint8
        source.native_row_width = 256  # ceil(1000/4)=250, round up64 ->256
        source.allows_direct_native_fill = True
        source.decode_tile_matches_chunk = True
        source.native_host_staging_copies = 0
        profile.pop('source_controls')
        profile['input'].update(native_encoding='pgen_2bit', transfer_bytes_per_variant=256,
                                host_staging_copies=0)
        with patch('torchgwas.scan_gpu.resolve_statistics_backend', return_value='native_fused'):
            result = apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
            self.assertEqual(result['workload']['samples'], 1000)
            self.assertEqual(result['estimate']['native_row_width'], 256)
            self.assertEqual(result['estimate']['input_transfer_bytes'], 2560000)
            self.assertEqual(result['estimate']['gpu_memory_bytes'], 10000*(2*256+8*1000+16*9))
            for field, value in [('native_row_width', 250), ('native_input_bytes_per_value', 4),
                                 ('transfer_bytes_per_variant', 1000)]:
                invalid = dict(profile, input=dict(profile['input'], **{field: value}))
                with self.assertRaises(ValueError):
                    apply_pipeline_profile(source, np.ones((1000, 8)), None, invalid)
        for n, padded in ((1,64), (256,64), (257,128), (22250,5568)):
            w = Workload(10, max(n,4), 1, 0, statistics_kernel='native_fused')
            if n == 1:
                w = replace(w, samples=4)
            p = replace(self.p, native_encoding='pgen_2bit', native_row_width=padded,
                        native_input_bytes_per_value=1, transfer_bytes_per_variant=padded, max_stored_bytes_per_variant=0,
                        decode_on_gpu=False, gpu_decode_seconds_per_variant=0)
            self.assertEqual(estimate(w,p,self.h,self.plan)['input_transfer_bytes'], 10*padded)

    def profile_fixture(self):
        source = SimpleNamespace(shape=(1000, 10000), input_bytes=2500000,
                    native_dtype=np.uint8, decode_batch_size=8, decode_workers=1,
                    reader_workers=1, prefetch_chunks=2, read_batch_bytes=1024,
                    read_workers=1, read_ahead_batches=2)
        profile = dict(hardware=asdict(self.h),
                    input=asdict(replace(self.p, decode_on_gpu=False,
                        transfer_bytes_per_variant=1000, decoded_bytes_per_value=1,
                        gpu_decode_seconds_per_variant=0,
                        cpu_decode_core_seconds_per_variant=1e-5)),
                    candidates=dict(read_variants=[1000], decode_variants=[500],
                        chunk_variants=[1000], workers=[2], depths=[2]),
                    source_controls=dict(read_workers=4, read_ahead_batches=4))
        return source, profile

    def test_profile_uses_actual_rank_and_applies_supported_controls(self):
        source, profile = self.profile_fixture()
        x = np.arange(1000, dtype=float)
        covariates = np.column_stack((np.ones(1000), x, 2 * x))
        result = apply_pipeline_profile(source, np.ones((1000, 8)), covariates, profile)
        self.assertEqual(result['workload'], dict(variants=10000, samples=1000,
                         traits=8, covariates=1, output_bytes_per_test=0.0, statistics_kernel='torch'))
        self.assertEqual(source.decode_batch_size, 500)
        self.assertEqual(source.decode_workers, 2)
        self.assertEqual(source.reader_workers, 2)
        self.assertEqual(source.read_batch_bytes, 250000)
        self.assertEqual(source.read_workers, 4)
        self.assertEqual(result['scan_kwargs']['chunk_size'], 1000)
        self.assertEqual(result['stored_bytes_evidence'], 'source.input_bytes')

    def test_profile_mismatches_rejected_before_source_mutation(self):
        for change in ('bytes', 'samples', 'workload', 'decode_cost', 'transfer'):
            source, profile = self.profile_fixture()
            phenotype = np.ones((1000, 8))
            if change == 'bytes': profile['input']['stored_bytes'] += 1
            elif change == 'samples': phenotype = phenotype[:-1]
            elif change == 'workload': profile['workload'] = {'variants': 99}
            elif change == 'transfer': profile['input']['transfer_bytes_per_variant'] = 4000
            else: del profile['input']['cpu_decode_core_seconds_per_variant']
            with self.assertRaises(ValueError):
                apply_pipeline_profile(source, phenotype, None, profile)
            self.assertEqual(source.decode_batch_size, 8)
            self.assertEqual(source.read_workers, 1)

    def test_gpu_profile_uses_actual_cpu_parse_workers(self):
        source, profile = self.profile_fixture()
        profile['input']['decode_on_gpu'] = True
        profile['candidates']['workers'] = [4]
        result = apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
        self.assertEqual(result['model_constraints']['cpu_parse_workers'], 4)
        self.assertEqual(source.read_workers, 4)
        self.assertEqual(result['estimate']['resource_seconds']['cpu_decode'], .05)
        source, profile = self.profile_fixture()
        profile['input']['decode_on_gpu'] = True
        with self.assertRaisesRegex(ValueError, 'read_workers'):
            apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)

    def test_declared_gpu_parse_mapping_and_rounded_tile_memory(self):
        source, profile = self.profile_fixture()
        del source.read_workers
        source.cpu_parse_worker_parameter = 'reader_workers'
        source.decode_tile_multiple_of_chunk = True
        source._n_bgen_samples = 1500
        profile['input']['decode_on_gpu'] = True
        profile['candidates']['workers'] = [4]
        profile['candidates']['decode_variants'] = [1500]
        result = apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
        self.assertEqual(source.reader_workers, 4)
        self.assertFalse(hasattr(source, 'read_workers'))
        self.assertEqual(source.decode_batch_size, 2000)
        self.assertEqual(result['estimate']['effective_geometry']['decode_variants'], 2000)
        self.assertEqual(result['estimate']['effective_geometry']['read_variants'], 2000)
        self.assertEqual(result['decode_input_samples'], 1500)
        self.assertGreaterEqual(result['estimate']['device_buffer_bytes'], 3 * 2000 * 1000 * 4)
        smaller = estimate(self.w, replace(self.p, decode_tile_multiple_of_chunk=True), self.h,
                           replace(self.plan, decode_variants=1000))
        bigger = estimate(self.w, replace(self.p, decode_tile_multiple_of_chunk=True), self.h,
                          replace(self.plan, decode_variants=1500))
        self.assertGreater(bigger['device_buffer_bytes'], smaller['device_buffer_bytes'])

    def test_profile_reports_unsupported_controls_as_advisory(self):
        source, profile = self.profile_fixture()
        del source.read_batch_bytes
        del source.decode_batch_size
        result = apply_pipeline_profile(source, np.ones((1000, 8)), None, profile)
        self.assertIn('read_batch_bytes', result['advisory_source_settings'])
        self.assertIn('decode_batch_size', result['advisory_source_settings'])
        self.assertNotIn('read_batch_bytes', result['applied_source_settings'])
        self.assertFalse(hasattr(source, 'read_batch_bytes'))
        self.assertEqual(result['estimate']['application_status'], 'conditional_on_advisory_controls')
        self.assertIsNone(result['estimate']['applied_io_bound_possible'])

    def test_scan_kwargs_match_selected_plan(self):
        result = plan_parameters(self.w, self.p, self.h, read_variants=(1000,),
                                 decode_variants=(1000,), chunk_variants=(1000,), workers=(2,), depths=(2,))
        self.assertEqual(result['scan_kwargs'], dict(chunk_size=1000, reader_workers=2, prefetch_chunks=2))
        self.assertEqual(result['source_recommendations'], dict(read_variants=1000, decode_variants=1000))

    def test_shared_gpu_service_is_added(self):
        r = estimate(self.w, self.p, self.h, self.plan)
        self.assertAlmostEqual(r['resource_seconds']['gpu'],
                               r['gpu_compute_seconds'] + .1)
        self.assertFalse(r['io_bound_possible'])

    def test_dimensions_and_linear_resource_scaling(self):
        a = estimate(self.w, self.p, self.h, self.plan)
        b = estimate(replace(self.w, variants=20000), replace(self.p, stored_bytes=5000000), self.h, self.plan)
        self.assertAlmostEqual(a['io_lower_bound_seconds'], .0025)
        for resource, value in a['resource_seconds'].items():
            self.assertAlmostEqual(b['resource_seconds'][resource], 2 * value)

    def test_cpu_and_dma_host_traffic_combined(self):
        p = replace(self.p, decode_on_gpu=False, transfer_bytes_per_variant=4000,
                    cpu_decode_core_seconds_per_variant=.0001, gpu_decode_seconds_per_variant=0)
        r = estimate(self.w, p, self.h, self.plan)
        self.assertAlmostEqual(r['resource_seconds']['host_memory'], (2500000 + 3 * 40000000 + 40000000 + 3 * 10000 * 69) / 1e11)
        self.assertAlmostEqual(r['resource_seconds']['cpu_decode'], .5)

    def test_format_does_not_change_statistical_compute(self):
        a = estimate(self.w, self.p, self.h, self.plan)
        b = estimate(self.w, replace(self.p, format='bgen-gpu', stored_bytes=5000000,
                     max_stored_bytes_per_variant=500, decoded_bytes_per_value=1), self.h, self.plan)
        self.assertEqual(a['gpu_compute_seconds'], b['gpu_compute_seconds'])
        self.assertEqual(2 * a['io_lower_bound_seconds'], b['io_lower_bound_seconds'])

    def test_results_transfer_without_durable_output(self):
        result = estimate(self.w, self.p, self.h, self.plan)
        self.assertEqual(result['result_d2h_bytes'], 10000 * (8 * 8 + 5))
        self.assertEqual(result['resource_seconds']['d2h'], 690000 / 1e10)
        self.assertFalse(result['memory_is_upper_bound'])

    def test_covariate_rank_excludes_intercept(self):
        with self.assertRaisesRegex(ValueError, 'degrees of freedom'):
            estimate(replace(self.w, samples=3, covariates=1), self.p, self.h, self.plan)

    def test_cpu_workers_cannot_exceed_pending_depth(self):
        p = replace(self.p, decode_on_gpu=False, cpu_decode_core_seconds_per_variant=.0001)
        two = estimate(self.w, p, self.h, self.plan)
        eight = estimate(self.w, p, self.h, replace(self.plan, workers=8))
        self.assertEqual(two['resource_seconds']['cpu_decode'], eight['resource_seconds']['cpu_decode'])

    def test_output_shares_storage(self):
        r = estimate(replace(self.w, output_bytes_per_test=4), self.p,
                     replace(self.h, output_bytes_per_second=1e6), self.plan)
        self.assertAlmostEqual(r['resource_seconds']['storage'], .0025 + .32)
        self.assertEqual(r['result_d2h_bytes'], 10000 * (16 * 8 + 5))

    def test_memory_constraint_and_shallow_depth(self):
        result = choose_plan(self.w, self.p, self.h, read_variants=(1000,),
                             decode_variants=(1000,), chunk_variants=(1000,), workers=(1,), depths=(2, 4))
        self.assertEqual(result['plan']['depth'], 2)
        with self.assertRaisesRegex(ValueError, 'no candidate'):
            choose_plan(self.w, self.p, replace(self.h, host_memory_bytes=1))

    def test_missing_output_bandwidth_and_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            estimate(replace(self.w, output_bytes_per_test=4), self.p, self.h, self.plan)
        with self.assertRaises(ValueError):
            estimate(self.w, self.p, replace(self.h, disk_bytes_per_second=float('nan')), self.plan)




class PeakMemoryTests(unittest.TestCase):
    """The host ring, which the device-only model could not see.

    Sizing trait blocks against `device_ring_bytes` alone chose a block that
    needed twice the host's memory in pinned pages. These pin the arithmetic
    that caught it.
    """

    # 33,417 subjects, chunk 4096, depth 32, packed two-bit rows.
    GEOMETRY = dict(chunk_variants=4096, depth=32, transfer_bytes_per_variant=8355.0)

    def test_result_ring_dominates_staging_at_high_trait_counts(self):
        staging = host_pinned_bytes(n_traits=1, reduction_width=None, **self.GEOMETRY)
        wide = host_pinned_bytes(n_traits=100000, reduction_width=None, **self.GEOMETRY)
        # beta and t are float32 per marker-trait cell: 8 bytes times chunk,
        # depth and traits, and nothing else grows with K.
        self.assertAlmostEqual(
            wide - staging,
            32 * 4096 * 8.0 * (100000 - 1), delta=1.0)
        # Compare the RINGS, not the totals: at K = 100,000 the result ring
        # is ~97x the staging ring, which is why sizing a trait block against
        # the device alone underestimates the host by an order of magnitude.
        staging_ring = 32 * 4096 * 8355.0
        self.assertGreater((wide - staging_ring) / staging_ring, 90.0)

    def test_logp_adds_float64_result_ring(self):
        without = host_pinned_bytes(
            n_traits=100, reduction_width=None, compute_log10_p=False,
            **self.GEOMETRY)
        with_logp = host_pinned_bytes(
            n_traits=100, reduction_width=None, compute_log10_p=True,
            **self.GEOMETRY)
        self.assertEqual(
            with_logp - without,
            32 * 4096 * 8 * 100)

    def test_reduction_replaces_the_trait_axis_in_the_result_ring(self):
        wide = host_pinned_bytes(n_traits=2085000, reduction_width=None, **self.GEOMETRY)
        narrow = host_pinned_bytes(n_traits=2085000, reduction_width=100, **self.GEOMETRY)
        # The reduced ring must not scale with K at all -- that is the whole
        # reason a voxel-scale scan is feasible reduced and not unreduced.
        self.assertEqual(narrow, host_pinned_bytes(
            n_traits=10, reduction_width=100, **self.GEOMETRY))
        self.assertGreater(wide / narrow, 1000)

    def test_device_native_source_stages_nothing_on_the_host(self):
        staged = host_pinned_bytes(n_traits=8, stage_on_host=True, **self.GEOMETRY)
        direct = host_pinned_bytes(n_traits=8, stage_on_host=False, **self.GEOMETRY)
        self.assertAlmostEqual(staged - direct, 32 * 4096 * 8355.0, delta=1.0)

    def test_trait_block_respects_the_host_ceiling_not_just_the_device(self):
        common = dict(n_samples=33417, n_traits=2085000, covariate_rank=27,
                      chunk_variants=4096, depth=32,
                      transfer_bytes_per_variant=8355.0,
                      device_memory_bytes=80 * 1024**3)
        device_only = auto_trait_block(**common)
        both = auto_trait_block(**common, host_memory_bytes=1000 * 1024**3,
                                trait_devices=4)
        self.assertLess(both, device_only)
        # And the block it chose must actually fit the ceiling it was given.
        self.assertLessEqual(
            host_pinned_bytes(chunk_variants=4096, depth=32, n_traits=both,
                              transfer_bytes_per_variant=8355.0) * 4,
            1000 * 1024**3 * 0.85)

    def test_reduced_scan_is_not_squeezed_by_the_host_ceiling(self):
        common = dict(n_samples=33417, n_traits=2085000, covariate_rank=27,
                      chunk_variants=4096, depth=32,
                      transfer_bytes_per_variant=8355.0,
                      device_memory_bytes=80 * 1024**3)
        device_only = auto_trait_block(**common)
        reduced = auto_trait_block(**common, host_memory_bytes=1000 * 1024**3,
                                   reduction_width=100, trait_devices=4)
        self.assertEqual(reduced, device_only)

    def test_peak_report_divides_device_and_multiplies_host_by_shards(self):
        args = dict(chunk_variants=4096, depth=32, n_samples=33417,
                    n_traits=200000, covariate_rank=27,
                    transfer_bytes_per_variant=8355.0, trait_block=50000)
        one = predicted_peak_bytes(**args, trait_devices=1)
        four = predicted_peak_bytes(**args, trait_devices=4)
        # Sharding does not change what one card holds; it changes how many
        # result rings are pinned in the one host address space.
        self.assertEqual(one['gpu_bytes_per_device'], four['gpu_bytes_per_device'])
        self.assertAlmostEqual(four['host_pinned_bytes'],
                               4 * one['host_pinned_bytes'], delta=1.0)
        self.assertEqual(four['trait_passes'], 4)

    def test_peak_report_never_blocks_wider_than_the_trait_axis(self):
        peak = predicted_peak_bytes(
            chunk_variants=4096, depth=32, n_samples=1000, n_traits=10,
            covariate_rank=2, transfer_bytes_per_variant=250.0,
            trait_block=50000)
        self.assertEqual(peak['trait_block'], 10)
        self.assertEqual(peak['trait_passes'], 1)

class MeasuredDevicePeakTests(unittest.TestCase):
    """The memory model against real measured peaks, per transport shape.

    Each case here is a bug the model actually had. They are pinned because
    every one of them was invisible until a format was measured at full scale:
    two were wrong by 2.5x and 3.2x while the model looked healthy on the
    formats that happened to share the default transport.
    """

    N = 35365
    C = 27

    def test_fused_packed_kernel_has_no_centred_genotype(self):
        """BED reads two-bit codes directly; there is no float32 chunk.

        Charging one put the model at 2.415 GB against a measured 0.803 GB --
        over by 3x, and in the dangerous direction: it would refuse chunk
        sizes the card can actually take.
        """
        common = dict(chunk_variants=4096, depth=18, n_samples=self.N,
                      n_traits=128, covariate_rank=self.C,
                      transfer_bytes_per_variant=8842.0, decode_on_gpu=True)
        charged = device_ring_bytes(**common)
        fused = device_ring_bytes(**common, fused_packed_statistics=True)
        self.assertAlmostEqual(charged - fused,
                               3.0 * 4096 * self.N * 4.0, delta=1.0)
        measured = 0.803e9
        self.assertLess(abs(measured / fused - 1.0), 0.25)

    def test_device_decoder_output_ring_is_the_genotype(self):
        """BGEN's device memory is its decoded float32 output ring.

        Measured exactly: 8 * 4096 * 35365 * 4 = 4,634,951,680 B against a
        measured peak of 4.635 G. Charging a centred copy as well put it at
        6.379 G -- over by precisely a three-slot centred term.
        """
        predicted = device_ring_bytes(
            chunk_variants=4096, depth=8, n_samples=self.N, n_traits=8,
            covariate_rank=self.C, transfer_bytes_per_variant=0.0,
            decode_on_gpu=True, gpu_decoded_bytes_per_variant=4.0 * self.N)
        self.assertLess(abs(4.635e9 / predicted - 1.0), 0.02)

    def test_host_staged_transports_scale_with_their_width(self):
        """The staging ring is what differs 16x between packed and float32."""
        packed = device_ring_bytes(
            chunk_variants=4096, depth=32, n_samples=self.N, n_traits=512,
            covariate_rank=self.C, transfer_bytes_per_variant=8896.0)
        dosage = device_ring_bytes(
            chunk_variants=4096, depth=32, n_samples=self.N, n_traits=512,
            covariate_rank=self.C, transfer_bytes_per_variant=141460.0)
        # Measured 2.47 G and 19.84 G respectively.
        self.assertLess(abs(2.47e9 / packed - 1.0), 0.15)
        self.assertLess(abs(19.84e9 / dosage - 1.0), 0.15)

    def test_explain_names_the_binding_term(self):
        """The decomposition must identify what to change, not just the total."""
        dosage = explain_memory(
            chunk_variants=4096, depth=32, n_samples=self.N, n_traits=512,
            covariate_rank=self.C, transfer_bytes_per_variant=141460.0)
        self.assertEqual(dosage['binding_term'], 'staging_ring')
        self.assertGreater(dosage['binding_fraction'], 0.9)
        # At K=512 the design is a rounding error, so blocking traits is not
        # the lever -- the model must say so rather than leave it to be found.
        self.assertLess(dosage['share']['design'], 0.01)

class TimeDecompositionTests(unittest.TestCase):
    """`explain_time` must name the binding resource and separate the terms.

    Each assertion here corresponds to a mistake the model actually made:
    folding the metadata parse into the scan term (it is reported separately as
    `open`), and assuming perfect overlap when the measurements sit nearer the
    serialised bound.
    """

    RATES = dict(disk_bytes_per_second=9.0e9, h2d_bytes_per_second=55.0e9,
                 d2h_bytes_per_second=55.0e9, gemm_flops_per_second=38.0e12,
                 write_bytes_per_second=1.46e9)

    def shape(self, **over):
        base = dict(variants=8_931_083, samples=35365, traits=128,
                    covariate_rank=27, chunk_variants=4096,
                    stored_bytes=35_932_909_454,
                    transfer_bytes_per_variant=8896.0,
                    output_bytes_per_test=8.0, **self.RATES)
        base.update(over)
        return explain_time(**base)

    def test_scan_excludes_the_metadata_parse(self):
        """The parse is reported as `open`, so folding it into scan double-counts."""
        without = self.shape()
        with_parse = self.shape(metadata_bytes=236e6,
                                text_parse_bytes_per_second=236e6 / 8.5)
        self.assertAlmostEqual(without['scan_seconds'],
                               with_parse['scan_seconds'], places=6)
        self.assertAlmostEqual(with_parse['open_seconds'], 8.5, places=1)
        # End to end DOES pay it.
        self.assertGreater(with_parse['end_to_end_seconds'],
                           without['end_to_end_seconds'])

    def test_overlap_interpolates_between_max_and_sum(self):
        """overlap=1 is the slowest stage; overlap=0 is the sum of all of them."""
        perfect = self.shape(overlap=1.0)
        serial = self.shape(overlap=0.0)
        terms = perfect['resource_seconds']
        self.assertAlmostEqual(perfect['steady_state_seconds'],
                               max(terms.values()), places=6)
        self.assertAlmostEqual(serial['steady_state_seconds'],
                               sum(terms.values()), places=6)
        self.assertGreater(serial['scan_seconds'], perfect['scan_seconds'])

    def test_binding_term_follows_the_rates_not_the_format(self):
        """Which resource binds must change when the machine changes.

        This is the property that makes the model portable: on a box with a
        slow bus the same workload should become transport-bound. Measured,
        the A100 has 6.66 GB/s H2D against the H100's 55.26 -- so a model that
        always names the same resource is not reading its inputs.
        """
        # The PACKED transport is where the flip actually happens. On the
        # float32 dosage transport H2D dominates at both rates (22.97 s vs
        # 189.7 s), so the binding term legitimately does not change -- picking
        # that case tests nothing, which is how this test first failed.
        fast_bus = self.shape()                       # 8,896 B/variant
        slow_bus = self.shape(h2d_bytes_per_second=6.66e9)
        self.assertEqual(fast_bus['binding_term'], 'disk')
        self.assertEqual(slow_bus['binding_term'], 'h2d')
        # And the term itself must scale with the rate, not merely reorder.
        self.assertAlmostEqual(
            slow_bus['resource_seconds']['h2d']
            / fast_bus['resource_seconds']['h2d'], 55.0 / 6.66, places=2)


class FinitePipelineTests(TimeDecompositionTests):
    def test_one_chunk_pays_each_stage_once(self):
        r = self.shape(variants=100, chunk_variants=4096)
        self.assertEqual(r['chunks'], 1)
        self.assertAlmostEqual(r['scan_seconds'], sum(r['resource_seconds'].values()))

    def test_serial_pipeline_has_no_extra_fill_charge(self):
        r = self.shape(variants=100, overlap=0)
        self.assertAlmostEqual(r['scan_seconds'], sum(r['resource_seconds'].values()))

    def test_partial_chunk_matches_explicit_schedule(self):
        from torchgwas.pipeline_model import finite_pipeline_seconds
        # Stages 2,3 seconds per full chunk; two full chunks and a half.
        # Completions: (2,5), (4,8), (5,9.5).
        self.assertAlmostEqual(finite_pipeline_seconds([5, 7.5], 5, 2), 9.5)

    def test_fully_overlapped_slow_writer_still_limits_total(self):
        r = self.shape(write_overlap=1, write_bytes_per_second=1)
        self.assertAlmostEqual(r['end_to_end_seconds'], r['write_seconds'])

    def test_launch_count_includes_the_partial_chunk(self):
        r = self.shape(variants=4097, chunk_variants=4096, per_chunk_seconds=1)
        self.assertEqual(r['chunks'], 2)
        self.assertEqual(r['per_chunk_seconds_total'], 2)


class GemmRateByWidthTests(unittest.TestCase):
    """Achieved FP32 depends strongly on design width, so one scalar is wrong.

    Measured with CUDA events around the statistics kernel, 400,000 variants:
    8.35 TFLOPS at width 36, 8.71 at 60, 20.57 at 156, 20.28 at 540 -- a 2.5x
    spread. The model carried a single 44.88 TFLOPS calibrated at width 540,
    which over-estimates the GPU about 2x at K=128 and 5x at K=8. It stayed
    invisible end to end because the GEMM rarely binds.
    """

    RATES = TimeDecompositionTests.RATES
    CURVE = {36: 8.35e12, 60: 8.71e12, 156: 20.57e12, 540: 20.28e12}

    def shape(self, **over):
        base = dict(variants=400_000, samples=35365, traits=128,
                    covariate_rank=27, chunk_variants=4096,
                    stored_bytes=400_000 * 8842.0,
                    transfer_bytes_per_variant=8842.0,
                    output_bytes_per_test=8.0, **self.RATES)
        base.update(over)
        return explain_time(**base)

    def test_a_scalar_still_works(self):
        """Every caller passed a float before this existed."""
        result = self.shape(gemm_flops_per_second=20.0e12)
        self.assertAlmostEqual(
            result['resource_seconds']['gemm'],
            400_000 * 2.0 * 35365 * 156 / 20.0e12, places=6)

    def test_a_curve_is_interpolated_at_the_design_width(self):
        """K=128 with 27 covariates is width 156, a measured point."""
        result = self.shape(gemm_flops_per_second=self.CURVE)
        self.assertAlmostEqual(
            result['resource_seconds']['gemm'],
            400_000 * 2.0 * 35365 * 156 / 20.57e12, places=6)

    def test_a_narrow_design_is_priced_much_slower_than_a_wide_one(self):
        """The whole point: 8.35 TFLOPS at width 36 against 20.28 at 540.

        With a scalar these two differ only by the FLOP count; with the curve
        the narrow design is additionally penalised for being memory-bound.
        """
        narrow = self.shape(traits=8, gemm_flops_per_second=self.CURVE)
        wide = self.shape(traits=512, gemm_flops_per_second=self.CURVE)
        flops_ratio = (512 + 28) / (8 + 28)
        seconds_ratio = (wide['resource_seconds']['gemm']
                         / narrow['resource_seconds']['gemm'])
        # Wide does far more work but at a better rate, so it takes less time
        # per FLOP -- the seconds ratio must fall short of the FLOP ratio.
        self.assertLess(seconds_ratio, flops_ratio)

    def test_widths_outside_the_measured_span_clamp_rather_than_extrapolate(self):
        """No evidence exists past the ends, so invent none."""
        from torchgwas.pipeline_model import gemm_rate_at_width
        self.assertEqual(gemm_rate_at_width(self.CURVE, 1), 8.35e12)
        self.assertEqual(gemm_rate_at_width(self.CURVE, 10_000), 20.28e12)

    def test_interpolation_between_measured_points(self):
        from torchgwas.pipeline_model import gemm_rate_at_width
        got = gemm_rate_at_width(self.CURVE, 108)   # midway 60 -> 156
        self.assertAlmostEqual(got, (8.71e12 + 20.57e12) / 2, delta=1e10)

    def test_an_empty_table_is_refused(self):
        from torchgwas.pipeline_model import gemm_rate_at_width
        with self.assertRaises(ValueError):
            gemm_rate_at_width({}, 156)


class DecodeTermTests(unittest.TestCase):
    """Compression trades disk bytes for CPU work; the model must see both.

    Without a decode term the model predicted the zstd hard-call store at
    7.63x where the measurement gives 1.84x -- a model with only a saving and
    no cost can reach no other answer. These numbers are the real ones:
    full genome, cold, two interleaved rounds.
    """

    RATES = TimeDecompositionTests.RATES

    # bed 78,968,635,889 B and the store 10,370,816,693 B over 8,931,083
    # variants; the store expands back to the bed's width.
    BED_PER_VARIANT = 8842.0
    STORE_PER_VARIANT = 1161.0
    # Achieved, not peak: measured bed K=1 scan of 17.18 s over 79.0 GB. The
    # 9.00 GB/s 16-reader figure is not what the scan gets.
    ACHIEVED_DISK = 4.7e9

    def shape(self, **over):
        base = dict(variants=8_931_083, samples=35365, traits=1,
                    covariate_rank=27, chunk_variants=4096,
                    stored_bytes=78_968_635_889,
                    transfer_bytes_per_variant=self.BED_PER_VARIANT,
                    output_bytes_per_test=8.0, **self.RATES)
        # Set the achieved rate BEFORE the override, so a caller can still
        # replace it. Doing it after with `.get` silently kept RATES' 9.00
        # GB/s peak, because that key is always present -- which made this
        # helper quietly measure the wrong machine.
        base['disk_bytes_per_second'] = self.ACHIEVED_DISK
        base.update(over)
        return explain_time(**base)

    def test_no_decode_term_unless_asked_for(self):
        """An uncompressed format must not grow a resource it does not use."""
        self.assertNotIn('decode', self.shape()['resource_seconds'])

    def test_a_store_pays_decode_on_the_bytes_it_produces(self):
        result = self.shape(
            stored_bytes=10_370_816_693,
            transfer_bytes_per_variant=self.STORE_PER_VARIANT,
            decoded_bytes_per_variant=self.BED_PER_VARIANT,
            decode_bytes_per_second=11.3e9)
        self.assertIn('decode', result['resource_seconds'])
        # 79.0 GB expanded at 11.3 GB/s.
        self.assertAlmostEqual(result['resource_seconds']['decode'], 7.0,
                               delta=0.5)

    def test_the_decode_term_can_bind(self):
        """It is the whole point: the store reads less and works more.

        Its disk term is 2.2 s against the bed's 16.8 s, so if decode could
        never bind, the model would keep predicting a speedup that the
        measurement says is not there.
        """
        result = self.shape(
            stored_bytes=10_370_816_693,
            transfer_bytes_per_variant=self.STORE_PER_VARIANT,
            decoded_bytes_per_variant=self.BED_PER_VARIANT,
            decode_bytes_per_second=11.3e9)
        self.assertEqual(result['binding_term'], 'decode')

    def test_the_store_no_longer_looks_like_a_7x_win_at_k1(self):
        """The measured ratio is 1.84x (17.18 s against 9.32 s).

        With the decode term the model lands in that neighbourhood instead of
        promising seven times. This is the regression that matters: the old
        model could not express the cost at all.
        """
        bed = self.shape()['scan_seconds']
        store = self.shape(
            stored_bytes=10_370_816_693,
            transfer_bytes_per_variant=self.STORE_PER_VARIANT,
            decoded_bytes_per_variant=self.BED_PER_VARIANT,
            decode_bytes_per_second=11.3e9)['scan_seconds']
        self.assertGreater(bed / store, 1.3)
        self.assertLess(bed / store, 3.0)


class SerialWriteTests(unittest.TestCase):
    """The sumstats write does not overlap the scan. Measured, not assumed.

    `end_to_end` was `max(scan, write)` on the assumption that the writer runs
    alongside. At K=512 on the full genome, enabling binary sumstats takes the
    scan from 14.90 s to 22.42 s -- **+7.52 s** -- where `max` predicts +0,
    because a 7.55 s write fits entirely inside a 14.90 s scan. And
    14.90 + 7.52 = 22.42 to the hundredth.
    """

    RATES = TimeDecompositionTests.RATES

    def shape(self, **over):
        base = dict(variants=8_931_083, samples=35365, traits=512,
                    covariate_rank=27, chunk_variants=4096,
                    stored_bytes=78_968_635_889,
                    transfer_bytes_per_variant=8842.0,
                    output_bytes_per_test=8.0, **self.RATES)
        base['disk_bytes_per_second'] = 78_968_635_889 / 14.21
        base['write_bytes_per_second'] = (8_931_083 * 512 * 8.0) / 7.52
        base['overlap'] = 0.93
        base.update(over)
        return explain_time(**base)

    def test_the_write_adds_to_the_scan_rather_than_hiding_inside_it(self):
        result = self.shape()
        self.assertGreater(result['scan_seconds'], result['write_seconds'],
                           'fixture must have a write that would fit inside '
                           'the scan, or this tests nothing')
        self.assertAlmostEqual(
            result['end_to_end_seconds'],
            result['scan_seconds'] + result['write_seconds'], places=6)

    def test_it_reproduces_the_measured_end_to_end(self):
        """22.42 s measured, cold, full genome, K=512, binary sumstats."""
        self.assertAlmostEqual(self.shape()['end_to_end_seconds'], 22.42,
                               delta=0.5)

    def test_full_write_overlap_recovers_the_old_behaviour(self):
        """Kept as a parameter for an implementation that does overlap it."""
        result = self.shape(write_overlap=1.0)
        self.assertAlmostEqual(result['end_to_end_seconds'],
                               result['scan_seconds'], places=6)

    def test_a_scan_with_no_output_pays_no_write(self):
        result = self.shape(output_bytes_per_test=0.0)
        self.assertEqual(result['write_seconds'], 0.0)
        self.assertAlmostEqual(result['end_to_end_seconds'],
                               result['scan_seconds'], places=6)


class ValidityRangeTests(unittest.TestCase):
    """The model must say where it stops knowing.

    A calculator meant to guide decisions is worse than useless when it returns
    a confident number outside the range it was checked in. Measured against
    `explain_time` plus the 2.10 s startup, the model is within ~10% between
    chunk 1,024 and 4,096 and then departs badly: +0.97 s at 8,192 on the quiet
    curve, +20.84 s at 65,536 on the cold one, with no term in the model for
    any of it. See `benchmarks/direct_chunk_curve_check.py`.
    """

    RATES = TimeDecompositionTests.RATES

    def shape(self, **over):
        return TimeDecompositionTests.shape(self, **over)

    def test_the_operating_point_is_inside_the_validated_range(self):
        self.assertTrue(
            self.shape(chunk_variants=MAX_AUTO_CHUNK_VARIANTS)
            ['chunk_within_validated_range'])

    def test_a_chunk_past_the_ceiling_is_flagged(self):
        result = self.shape(chunk_variants=65_536)
        self.assertFalse(result['chunk_within_validated_range'])
        self.assertTrue(any('validated ceiling' in c
                            for c in result['caveats']))

    def test_a_fitted_overlap_announces_itself_not_the_structural_one(self):
        """This assertion is INVERTED from what it used to be, deliberately.

        It used to require a caveat on the 1.0 default ("mildly optimistic
        against a measured 0.93") and none on a sub-1.0 value ("a caller that
        has measured it should not be nagged"). That had it backwards. 0.93 was
        not measured -- it was SOLVED FOR by inverting this model against wall
        clocks -- and re-checked once the decode rate was measured rather than
        inferred, it predicts worse than the structural 1.0 at every point on
        both regimes: 1.10-1.15x against 0.99-1.04x on the zstd store, and
        worse at all four rungs of the raw .bed ladder it was derived from.

        So the value that has to announce itself is the fitted one.
        """
        self.assertFalse(any('overlap' in c for c in self.shape()['caveats']),
                         "the structural 1.0 is the best available value and "
                         "should not be warned about")
        fitted = self.shape(overlap=0.5)['caveats']
        self.assertTrue(any('overlap' in c for c in fitted))
        self.assertIn('FITTED', " ".join(fitted))

    def test_a_plan_inside_the_range_carries_no_caveats(self):
        self.assertEqual(
            self.shape(chunk_variants=2048, overlap=1.0)['caveats'], [])

    def test_auto_chunk_never_returns_an_unvalidated_plan(self):
        """The cap must stay inside the validity ceiling.

        These used to be the same number, with the cap "placed at the edge of
        validity". They no longer are: the edge has been measured further out
        (8,192 and 16,384 now predict as well as 4,096 does), so the ceiling
        moved to 16,384 while the cap stayed at 4,096 -- not because larger is
        worse, it is flat, but because larger buys nothing measurable and costs
        device memory in proportion. The invariant that matters is unchanged
        and is what this asserts: if `auto_chunk_variants` could hand back a
        chunk the model has not been checked at, the planner would routinely
        predict outside its own range -- which is what happened when the cap
        was 65,536.
        """
        for device_bytes in (8 << 30, 40 << 30, 80 << 30, 640 << 30):
            chunk = auto_chunk_variants(
                n_samples=35365, n_traits=128, covariate_rank=27,
                transfer_bytes_per_variant=8896.0,
                device_memory_bytes=device_bytes, depth=32)
            self.assertLessEqual(chunk, VALIDATED_MAX_CHUNK_VARIANTS,
                                 f'{device_bytes} bytes gave chunk {chunk}')


if __name__ == '__main__':
    unittest.main()
