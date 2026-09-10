from __future__ import annotations

import unittest

from torchgwas.runtime import predict_runtime, predict_runtime_from_hardware


class RuntimeModelTestCase(unittest.TestCase):
    def test_joint_rates_reproduce_remote_contention_ceiling(self):
        result = predict_runtime(
            n_variants=8_086_101,
            n_samples=22_250,
            n_traits=128,
            covariate_rank=20,
            chunk_size=2_500,
            decoded_bytes_per_value=1,
            compression_ratio=3,
            disk_gbps=10,
            decode_gbps=50,
            h2d_gbps=25,
            measured_gemm_tflops=19.1,
            joint_decode_chunks_per_second=763.3,
            joint_h2d_chunks_per_second=482.8,
            joint_compute_chunks_per_second=482.8,
        )
        self.assertEqual(result["model"], "measured_joint_pipeline")
        self.assertEqual(result["n_chunks"], 3235)
        self.assertAlmostEqual(result["scan_seconds"], 3235 / 482.8)

    def test_hardware_model_reproduces_project_contended_ceiling(self):
        result = predict_runtime_from_hardware(
            n_variants=8_086_101,
            n_samples=22_250,
            n_traits=128,
            covariate_rank=27,
            chunk_size=2_500,
            gpu_fp32_tflops=19.1,
            gpu_count=2,
            gemm_efficiency=12.9 / 19.1,
            cpu_frequency_ghz=2.3,
            decode_threads=40,
            decode_gbps_per_core_at_reference=1.52,
            reference_cpu_frequency_ghz=2.3,
            decode_thread_efficiency=55.3 / (1.52 * 40),
            disk_gbps=100.0,
            h2d_gbps=42.5,
            non_gemm_ms_per_chunk_per_gpu=0.818,
            pipeline_contention_factor=765.1 / 482.8,
        )
        self.assertEqual(result["n_chunks"], 3235)
        self.assertAlmostEqual(result["estimated_decode_gbps"], 55.3, places=6)
        self.assertAlmostEqual(result["scan_seconds"], 6.71, delta=0.1)

    def test_partial_joint_rates_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "all three"):
            predict_runtime(
                n_variants=100,
                n_samples=10,
                n_traits=2,
                covariate_rank=1,
                chunk_size=20,
                decoded_bytes_per_value=4,
                compression_ratio=16,
                disk_gbps=1,
                decode_gbps=1,
                h2d_gbps=1,
                measured_gemm_tflops=1,
                joint_decode_chunks_per_second=5,
            )


if __name__ == "__main__":
    unittest.main()
