from __future__ import annotations
from dataclasses import replace
import unittest
from torchgwas.runtime import predict_runtime, predict_runtime_from_hardware
from torchgwas.pipeline_model import Workload, InputProfile, Hardware, PipelinePlan, estimate


class RuntimeModelTestCase(unittest.TestCase):
    def arguments(self):
        return dict(n_variants=100, n_samples=40, n_traits=2, covariate_rank=1,
                    chunk_size=20, decoded_bytes_per_value=4,
                    h2d_bytes_per_value=.25, compression_ratio=1, disk_gbps=1,
                    decode_gbps=1, h2d_gbps=1, measured_gemm_tflops=1)

    def test_delegates_to_canonical_resource_model(self):
        with self.assertWarns(DeprecationWarning):
            result = predict_runtime(**self.arguments())
        canonical = estimate(Workload(100, 40, 2, 1),
            InputProfile('legacy-cpu-decode', 16000, 10, 4,
                         cpu_decode_core_seconds_per_variant=160 / 1e9),
            Hardware(1e9, 1e9, 1e300, 1e12, 1e300, 2**62, 2**62, 1, d2h_bytes_per_second=1e9),
            PipelinePlan(20, 20, 20, 1, 2))
        self.assertEqual(result['resource_seconds'], canonical['resource_seconds'])
        self.assertEqual(result['scan_seconds'], canonical['resource_lower_bound_seconds'])
        self.assertEqual(result['decoded_gb'], 16000 / 1e9)
        self.assertEqual(result['h2d_gb'], 1000 / 1e9)
        self.assertEqual(result['h2d_seconds_isolated'], 1e-6)

    def test_empirical_joint_and_contention_overrides_retired(self):
        for fields in ({'joint_decode_chunks_per_second': 5},
                       {'joint_decode_chunks_per_second': 5,
                        'joint_h2d_chunks_per_second': 6,
                        'joint_compute_chunks_per_second': 7},
                       {'contention_factor': 1.2}):
            with self.assertRaisesRegex(ValueError, 'retired'):
                predict_runtime(**self.arguments(), **fields)

    def test_negative_overhead_rejected(self):
        with self.assertRaises(ValueError):
            predict_runtime(**self.arguments(), setup_seconds=-1)

    def test_hardware_wrapper_supplies_rates_to_same_model(self):
        with self.assertWarns(DeprecationWarning):
            result = predict_runtime_from_hardware(
                n_variants=100, n_samples=40, n_traits=2, covariate_rank=1,
                gpu_fp32_tflops=10, gpu_count=2, gemm_efficiency=.5,
                cpu_frequency_ghz=2, decode_threads=4,
                decode_gbps_per_core_at_reference=1,
                reference_cpu_frequency_ghz=2, decode_thread_efficiency=.5,
                disk_gbps=1, h2d_gbps=10)
        self.assertEqual(result['estimated_decode_gbps'], 2)
        self.assertEqual(result['estimated_sustained_gemm_tflops'], 10)
        self.assertEqual(result['model'], 'canonical_partial_resource_lower_bound')


if __name__ == '__main__':
    unittest.main()
