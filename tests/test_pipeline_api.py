"""Exercise profile application through the real public API and CUDA pipeline."""
from dataclasses import asdict
import unittest
import numpy as np
import torch

from torchgwas.api import run_linear_gwas
from torchgwas.pipeline_model import Hardware, InputProfile


class ArrayStream:
    supports_fused_qc = True
    native_dtype = np.float32
    ndim = 2

    def __init__(self, data):
        self.data = data
        self.sample_ids = np.asarray([f's{i}' for i in range(data.shape[0])])
        self.marker_ids = np.asarray([f'm{i}' for i in range(data.shape[1])])
        self.input_bytes = data.nbytes
        self.reader_workers = self.decode_workers = 1
        self.decode_batch_size = 4
        self.prefetch_chunks = 2
        self.calls = []

    @property
    def shape(self):
        return self.data.shape

    @property
    def genotype(self):
        return self

    def __getitem__(self, key):
        return self.data[key]

    def iter_chunks(self, chunk_size, dtype=np.float32, reader_workers=None, prefetch_chunks=None):
        self.calls.append((chunk_size, reader_workers, prefetch_chunks))
        for start in range(0, self.shape[1], chunk_size):
            end = min(self.shape[1], start + chunk_size)
            yield start, end, np.asarray(self.data[:, start:end], dtype=dtype)


class PipelineAPITests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(211)
        self.data = rng.integers(0, 3, size=(32, 13)).astype(np.float32)
        self.phenotype = rng.normal(size=(32, 2))
        self.covariates = rng.normal(size=(32, 2))
        self.profile = {
            'hardware': asdict(Hardware(1e9, 1e9, 1e10, 1e12, 1e11,
                               1 << 28, 1 << 28, 8, d2h_bytes_per_second=1e9)),
            'input': asdict(InputProfile('array-cpu', self.data.nbytes, 32 * 4, 4,
                        cpu_decode_core_seconds_per_variant=1e-9,
                        gpu_decode_seconds_per_variant=0,
                        max_stored_bytes_per_variant=32 * 4)),
            'candidates': dict(chunk_variants=[2], read_variants=[4],
                               decode_variants=[4], workers=[1], depths=[2]),
        }

    @unittest.skipUnless(torch.cuda.is_available(), 'CUDA required for real pipeline integration')
    def test_explicit_settings_override_profile_and_drive_actual_scan(self):
        source = ArrayStream(self.data)
        result = run_linear_gwas(source, self.phenotype, self.covariates,
                    device='cuda', compute_dtype='float32', pipeline_profile=self.profile,
                    chunk_size=5, reader_workers=3, prefetch_chunks=3)
        self.assertEqual(source.calls, [(5, 3, 3)])
        self.assertEqual(result.run_metadata['pipeline_plan']['scan_kwargs'],
                         dict(chunk_size=5, reader_workers=3, prefetch_chunks=3))
        self.assertEqual(result.run_metadata['reader_workers'], 3)
        self.assertEqual(result.run_metadata['prefetch_chunks'], 3)
        self.assertEqual(result.run_metadata['pipeline_plan']['workload']['covariates'], 2)
        reference = run_linear_gwas(ArrayStream(self.data), self.phenotype,
                    self.covariates, device='cuda', compute_dtype='float32', chunk_size=4)
        np.testing.assert_allclose([row['beta'] for row in result.table],
                                   [row['beta'] for row in reference.table], rtol=1e-4, atol=1e-5)
        # API must not mutate the caller's reusable candidate dictionary.
        self.assertEqual(self.profile['candidates']['chunk_variants'], [2])

    @unittest.skipUnless(torch.cuda.is_available(), 'CUDA required for backend/profile dispatch')
    def test_resolved_backend_mismatch_rejected_before_scan(self):
        class FallbackStream(ArrayStream):
            decode_backend = 'auto'
            backend_used = None
            def resolve_decode_backend(self, device):
                self.backend_used = 'cpu'
                return 'cpu'
        source = FallbackStream(self.data)
        self.profile['input']['decode_on_gpu'] = True
        with self.assertRaisesRegex(ValueError, 'decode placement'):
            run_linear_gwas(source, self.phenotype, self.covariates,
                            device='cuda', compute_dtype='float32', pipeline_profile=self.profile)
        self.assertEqual(source.backend_used, 'cpu')
        self.assertEqual(source.calls, [])

    def test_dense_profile_rejected_instead_of_ignored(self):
        with self.assertRaisesRegex(ValueError, 'streaming genotype source'):
            run_linear_gwas(self.data, self.phenotype, self.covariates,
                            device='cpu', pipeline_profile=self.profile)


if __name__ == '__main__':
    unittest.main()
