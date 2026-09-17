"""Independent edge cases for optional native GPU reductions and finalization."""
import unittest
import numpy as np
import torch
from torchgwas.scan_gpu import available, prepare, finish
from torchgwas.linear import _dosage_statistics


@unittest.skipUnless(torch.cuda.is_available() and available(), 'native CUDA statistics required')
class NativeStatisticsTests(unittest.TestCase):
    def test_nan_inf_range_and_constant_status_match_torch(self):
        rng = np.random.default_rng(881)
        raw = rng.uniform(0, 2, (37, 129)).astype(np.float32)
        raw[1, 0] = np.nan
        raw[2, 0] = np.inf
        raw[3] = 1.17
        raw[4, 0] = 3
        raw[5, :2] = [np.nan, 3]
        raw[6, 0] = -np.inf
        device_raw = torch.as_tensor(raw, device='cuda')
        design = torch.as_tensor(rng.normal(size=(129, 6)).astype(np.float32) / np.float32(np.sqrt(129)), device='cuda')
        phenotype_ss = torch.ones(3, device='cuda')
        # 129 samples, 3 covariate columns: df = present - 3 - 2, so the
        # offset is -5 and a complete variant lands back on 124.
        centered, ss, lo, hi, present = prepare(device_raw)
        beta, statistic, status = finish(centered @ design, ss, lo, hi,
                                         phenotype_ss, present, -5.0, True)
        expected = _dosage_statistics(device_raw, design, phenotype_ss, 3, 124,
                                      True, covariate_rank=3)
        np.testing.assert_array_equal(status.cpu(), expected[2].cpu())
        # Rows 1 and 5 carry a missing call. Masking means row 1 is now a
        # normal variant analysed on its observed samples, and row 5 is
        # rejected on the merits of its remaining value (3, out of range)
        # rather than pre-empted by the missing call. Rows 2 and 6 hold an
        # infinity, which is not missing and is not masked.
        np.testing.assert_array_equal(status.cpu().numpy()[1:7], [0, 3, 2, 3, 3, 3])
        valid = status.cpu().numpy() == 0
        for actual, reference in zip((beta, statistic), expected[:2]):
            np.testing.assert_allclose(actual.cpu().numpy()[valid], reference.cpu().numpy()[valid], rtol=2e-4, atol=2e-5)

    def test_native_inputs_centered_ss_and_nondefault_stream(self):
        rng = np.random.default_rng(883)
        stream = torch.cuda.Stream()
        for dtype, scale in ((np.int8, 1.0), (np.uint8, 127.5), (np.float32, 1.0)):
            raw = (rng.uniform(0, 2, (19, 22250)).astype(dtype) if dtype == np.float32
                   else rng.integers(0, 3 if dtype == np.int8 else 256, (19, 22250), dtype=dtype))
            raw[0] = 2 if dtype != np.uint8 else 255
            raw[0, 3] -= 1 if dtype != np.float32 else np.float32(0.01)
            missing = -9 if dtype != np.uint8 else None
            if missing is not None:
                raw[1, 0] = -9
            expected = raw.astype(np.float64) / scale
            if missing is not None:
                expected[raw == missing] = np.nan
            # uint8 conversion rounds to float32 before native mean accumulation.
            rounded = (raw.astype(np.float32) / scale).astype(np.float64)
            if missing is not None:
                rounded[raw == missing] = np.nan
            # Missing calls are masked: they take no part in the mean, and
            # their centred value is exactly zero rather than NaN.
            observed = ~np.isnan(rounded)
            mean = np.nanmean(rounded, axis=1, keepdims=True).astype(np.float32)
            expected_centered = np.where(
                observed, rounded.astype(np.float32) - mean, np.float32(0.0))
            expected_ss = (expected_centered.astype(np.float64) ** 2).sum(axis=1).astype(np.float32)
            device_raw = torch.as_tensor(raw, device='cuda')
            stream.wait_stream(torch.cuda.current_stream())
            with torch.cuda.stream(stream):
                centered, ss, lo, hi, present = prepare(device_raw, scale, missing)
            torch.cuda.current_stream().wait_stream(stream)
            np.testing.assert_allclose(centered.cpu(), expected_centered, rtol=1e-6, atol=2e-7, equal_nan=True)
            np.testing.assert_allclose(ss.cpu(), expected_ss, rtol=1e-6, atol=1e-8, equal_nan=True)
            with np.errstate(invalid='ignore'):
                np.testing.assert_allclose(lo.cpu(), np.nanmin(expected, axis=1), rtol=1e-6, equal_nan=True)
                np.testing.assert_allclose(hi.cpu(), np.nanmax(expected, axis=1), rtol=1e-6, equal_nan=True)


if __name__ == '__main__':
    unittest.main()