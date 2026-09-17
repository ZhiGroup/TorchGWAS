import unittest
import os
from unittest.mock import patch
from contextlib import contextmanager
import numpy as np
from scipy import stats
import torch
from torchgwas.linear import linear_scan, linear_scan_streaming_chunks
from torchgwas.api import run_linear_gwas


class Source:
    supports_fused_qc = True
    decode_workers = 1

    def __init__(self, values):
        self.values = values
        self.shape = values.shape
        self.sample_ids = np.asarray([str(i) for i in range(values.shape[0])])
        self.marker_ids = np.asarray([str(i) for i in range(values.shape[1])])
        self.passes = 0

    @property
    def genotype(self):
        return self

    def iter_chunks(self, chunk_size, dtype=np.float32, **kwargs):
        self.passes += 1
        for start in range(0, self.shape[1], chunk_size):
            end = min(self.shape[1], start + chunk_size)
            yield start, end, self.values[:, start:end].astype(dtype)


@unittest.skipUnless(torch.cuda.is_available(), "CUDA required")
class NativeScanTests(unittest.TestCase):
    @patch.dict(os.environ, {"TORCHGWAS_SCAN_PROFILE": "1"})
    def test_exact_ols_and_ring_reuse_across_chunk_sizes(self):
        rng = np.random.default_rng(12)
        g = rng.uniform(0, 2, (129, 83)).astype(np.float32)
        y = rng.normal(size=(129, 3))
        c = rng.normal(size=(129, 2))
        c = np.column_stack((c, c[:, 0]))  # Rank deficient covariates.
        reference = linear_scan(g, y, c, device="cpu", compute_dtype="float64")
        for chunk in (7, 19):
            source = Source(g)
            iterator, _ = linear_scan_streaming_chunks(source, y, c,
                chunk_size=chunk, device="cuda:0", prefetch_chunks=2)
            rows = list(iterator)  # Retention must survive result-ring reuse.
            for field in (2, 3, 4):
                actual = np.concatenate([row[field] for row in rows])
                np.testing.assert_allclose(actual, reference[field-2], rtol=2e-4, atol=2e-5)
            self.assertEqual(source.passes, 1)

    def test_api_single_pass_and_fused_missing_qc(self):
        rng = np.random.default_rng(7)
        g = rng.integers(0, 3, (97, 25)).astype(np.float32)
        g[:, 2] = 0
        g[:, 3] = 1
        g[:, 4] = 2
        g[:, 6] = np.float32(1.17)
        g[0, 5] = np.nan
        source = Source(g)
        result = run_linear_gwas(source, rng.normal(size=(97, 2)),
                                 chunk_size=7, device="cuda:0")
        self.assertEqual(source.passes, 1)
        # The missing call costs its own sample, not the variant, so nothing
        # is excluded for missingness; the four constant columns still are.
        self.assertEqual(result.qc_summary["genotype_exclusion_counts"],
                         {"missing": 0, "invariant": 4})

    def test_native_integer_transfer_matches_float64_ols(self):
        class Native(Source):
            def iter_native_chunks(self, chunk_size, **kwargs):
                self.passes += 1
                for start in range(0, self.shape[1], chunk_size):
                    end = min(self.shape[1], start + chunk_size)
                    yield start, end, self.values[:, start:end]

        rng = np.random.default_rng(31)
        y = rng.normal(size=(113, 3))
        c = rng.normal(size=(113, 2))
        for dtype, scale in ((np.int8, 1), (np.uint8, 127)):
            raw = rng.integers(0, 2 * scale + 1, (113, 43), dtype=dtype)
            source = Native(raw)
            source.native_dtype = dtype
            source.native_scale = scale
            expected = raw.astype(np.float64) / scale
            reference = linear_scan(expected, y, c, device="cpu", compute_dtype="float64")
            iterator, _ = linear_scan_streaming_chunks(source, y, c,
                chunk_size=7, device="cuda", prefetch_chunks=2)
            rows = list(iterator)
            for field in (2, 3, 4):
                np.testing.assert_allclose(np.concatenate([r[field] for r in rows]),
                    reference[field-2], rtol=2e-4, atol=2e-5)
            self.assertEqual(source.passes, 1)

    def test_nearly_fixed_dosage_is_stable_across_chunk_shapes(self):
        rng = np.random.default_rng(52)
        g = np.full((22250, 67), 2, dtype=np.float32)
        for j in range(g.shape[1]):
            indices = rng.choice(g.shape[0], 1 + j % 7, replace=False)
            g[indices, j] -= np.float32(0.01 + (j % 3) * 0.1)
        y = rng.normal(size=(g.shape[0], 3))
        c = rng.normal(size=(g.shape[0], 2))
        reference = linear_scan(g, y, c, device='cpu', compute_dtype='float64')
        for chunk in (7, 31):
            iterator, _ = linear_scan_streaming_chunks(Source(g), y, c,
                chunk_size=chunk, device='cuda', prefetch_chunks=2)
            rows = list(iterator)
            for field in (2, 3):
                np.testing.assert_allclose(np.concatenate([r[field] for r in rows]),
                    reference[field-2], rtol=2e-4, atol=2e-5)

    def test_direct_pinned_fill_preserves_values_missingness_and_close(self):
        class Direct(Source):
            allows_direct_native_fill = True
            validate_native_range = True
            native_missing_value = -9
            decode_workers = 3
            def iter_chunks(self, *args, **kwargs):
                raise AssertionError('direct fill should bypass allocating iterator')
            @contextmanager
            def native_reader_session(self):
                self.closed = False
                def fill(start, end, out):
                    np.copyto(out, self.values[:, start:end].T)
                try:
                    yield fill
                finally:
                    self.closed = True
        rng = np.random.default_rng(57)
        y = rng.normal(size=(129, 3))
        for dtype in (np.int8, np.float32):
            raw = (rng.integers(0, 3, (129, 43)).astype(dtype) if dtype == np.int8
                   else rng.uniform(0, 2, (129, 43)).astype(dtype))
            raw[1, 5] = -9
            source = Direct(raw)
            source.native_dtype = dtype
            # Column 5 carries one missing call. It is analysed on its
            # observed samples, which is the same estimand as replacing that
            # call with the column's observed mean -- so the reference covers
            # every column, including 5, rather than excluding it.
            imputed = raw.astype(np.float64)
            observed = imputed[:, 5] != -9
            imputed[~observed, 5] = imputed[observed, 5].mean()
            good = np.arange(raw.shape[1]) != 5
            reference = linear_scan(imputed, y, None,
                                    device='cpu', compute_dtype='float64')
            iterator, _ = linear_scan_streaming_chunks(source, y, None,
                chunk_size=7, device='cuda', prefetch_chunks=3)
            rows = list(iterator)
            # Column 5 has one missing call, so it spends 128 samples rather
            # than 129 and its t is scaled by sqrt(df_variant / df_full).
            # beta does not move: df enters only the standard error.
            n_samples = raw.shape[0]
            df_full = n_samples - 0 - 2
            t_scale = np.sqrt((df_full - 1) / df_full)
            for field in (2, 3, 4):
                actual = np.concatenate([r[field] for r in rows])
                self.assertTrue(np.isfinite(actual[5]).all())
                want = np.array(reference[field-2], copy=True)
                if field == 3:
                    want[5] = want[5] * t_scale
                np.testing.assert_allclose(actual[good], want[good],
                                           rtol=2e-4, atol=2e-5)
                if field == 2:
                    np.testing.assert_allclose(actual[5], want[5],
                                               rtol=2e-4, atol=2e-5)
            self.assertTrue(source.closed)
            self.assertEqual(source._last_scan_exclusion_counts,
                             {'missing': 0, 'invariant': 0})
            iterator, _ = linear_scan_streaming_chunks(source, y, None,
                chunk_size=7, device='cuda', prefetch_chunks=3)
            next(iterator)
            iterator.close()
            self.assertTrue(source.closed)
            raw[0, 3] = 3
            iterator, _ = linear_scan_streaming_chunks(source, y, None,
                chunk_size=7, device='cuda', prefetch_chunks=3)
            with self.assertRaisesRegex(ValueError, 'invalid native ALT1 dosage'):
                list(iterator)
            self.assertTrue(source.closed)

    def test_early_consumer_close_stops_producer(self):
        rng = np.random.default_rng(9)
        source = Source(rng.uniform(0, 2, (65, 100)).astype(np.float32))
        iterator, _ = linear_scan_streaming_chunks(source, rng.normal(size=(65, 1)),
            None, chunk_size=4, device="cuda:0", prefetch_chunks=2)
        next(iterator)
        iterator.close()

    @patch.dict(os.environ, {"TORCHGWAS_NATIVE_STATS": "1"})
    def test_packed_pgen_transport_matches_independent_ols(self):
        from torchgwas.scan_gpu import available
        if not available():
            self.skipTest('native statistics library required')
        class Packed(Source):
            allows_direct_native_fill = True
            validate_native_range = True
            native_dtype = np.dtype(np.int8)  # Public iterator contract.
            native_transfer_dtype = np.dtype(np.uint8)
            native_encoding = 'pgen_2bit'
            decode_workers = 3
            def iter_chunks(self, *args, **kwargs):
                raise AssertionError('packed transport must use direct fill')
            @contextmanager
            def native_reader_session(self):
                self.closed = False
                def fill(start, end, out):
                    assert out.dtype == np.uint8
                    assert out.shape == (end-start, self.native_row_width)
                    assert out.ctypes.data % 64 == 0
                    # Poison padding: it must never enter the statistical sample count.
                    out.fill(255)
                    for sample in range(self.shape[0]):
                        shift = 2 * (sample % 4)
                        calls = self.values[sample, start:end].astype(np.uint8)
                        out[:, sample // 4] &= np.uint8(255 ^ (3 << shift))
                        out[:, sample // 4] |= calls << shift
                try:
                    yield fill
                finally:
                    self.closed = True
        rng = np.random.default_rng(91)
        for samples in (97, 129):
            calls = rng.integers(0, 3, (samples, 43), dtype=np.uint8)
            calls[:, 2] = 0
            calls[:, 3] = 2
            calls[0, 5] = 3
            source = Packed(calls)
            source.native_row_width = ((samples + 3) // 4 + 63) // 64 * 64
            y = rng.normal(size=(samples, 3))
            c = rng.normal(size=(samples, 2))
            # Columns 2 and 3 are constant and stay excluded. Column 5 has
            # one missing call and is now analysed, so it leaves `good` and
            # joins the reference with that call set to the observed mean.
            good = ~np.isin(np.arange(43), [2, 3])
            imputed = calls.astype(np.float64)
            present = imputed[:, 5] != 3
            imputed[~present, 5] = imputed[present, 5].mean()
            reference = linear_scan(imputed[:, good], y, c,
                                    device='cpu', compute_dtype='float64')
            for chunk in (7, 19):
                iterator, _ = linear_scan_streaming_chunks(source, y, c,
                    chunk_size=chunk, device='cuda', prefetch_chunks=3)
                rows = list(iterator)
                # Column 5 carries a missing call and now uses its own df,
                # so its t is scaled; beta is untouched. The reference covers
                # `good`, which includes column 5, so scale that entry.
                df_full = samples - c.shape[1] - 2
                t_scale = np.sqrt((df_full - 1) / df_full)
                kept = np.flatnonzero(good)
                elsewhere = kept != 5
                for field in (2, 3, 4):
                    actual = np.concatenate([row[field] for row in rows])
                    self.assertTrue(np.isnan(actual[~good]).all())
                    self.assertTrue(np.isfinite(actual[5]).all())
                    np.testing.assert_allclose(
                        actual[good][elsewhere], reference[field-2][elsewhere],
                        rtol=2e-4, atol=2e-5)

                # Column 5 carries one missing call, so it spends one fewer
                # sample. beta is untouched -- df enters only the standard
                # error -- t scales by sqrt(df_variant / df_full), and the
                # p-value is the tail of a different distribution and has to
                # be recomputed rather than scaled.
                at = int(np.flatnonzero(kept == 5)[0])
                beta5 = np.concatenate([row[2] for row in rows])[5]
                t5 = np.concatenate([row[3] for row in rows])[5]
                p5 = np.concatenate([row[4] for row in rows])[5]
                np.testing.assert_allclose(beta5, reference[0][at],
                                           rtol=2e-4, atol=2e-5)
                want_t = reference[1][at] * t_scale
                np.testing.assert_allclose(t5, want_t, rtol=2e-4, atol=2e-5)
                np.testing.assert_allclose(
                    p5, 2.0 * stats.t.sf(np.abs(want_t), df_full - 1),
                    rtol=1e-3, atol=1e-6)
                self.assertTrue(source.closed)
                self.assertEqual(source._last_scan_exclusion_counts,
                                 {'missing': 0, 'invariant': 2})
            iterator, _ = linear_scan_streaming_chunks(source, y, c,
                chunk_size=7, device='cuda', prefetch_chunks=3)
            next(iterator)
            iterator.close()
            self.assertTrue(source.closed)

    def test_final_copy_release_error_is_not_lost(self):
        rng = np.random.default_rng(79)
        source = Source(rng.uniform(0, 2, (65, 7)).astype(np.float32))
        with patch('torchgwas.streaming.PinnedDosageLoader.release',
                   side_effect=RuntimeError('copy lease failed')):
            iterator, _ = linear_scan_streaming_chunks(source,
                rng.normal(size=(65, 2)), None, chunk_size=10,
                device='cuda:0', prefetch_chunks=2)
            with self.assertRaisesRegex(RuntimeError, 'copy lease failed'):
                list(iterator)

    def test_immediate_reader_error_is_propagated(self):
        class Broken(Source):
            def iter_chunks(self, *args, **kwargs):
                raise ValueError("reader initialization failed")
        source = Broken(np.ones((65, 20), dtype=np.float32))
        iterator, _ = linear_scan_streaming_chunks(source,
            np.arange(65, dtype=float)[:, None], None,
            chunk_size=4, device="cuda:0")
        with self.assertRaisesRegex(ValueError, "reader initialization failed"):
            next(iterator)


if __name__ == "__main__":
    unittest.main()
