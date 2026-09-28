"""triton_scan: scan_gpu's prepare/finish contract, and the scans that run on it by default."""
import numpy as np
import pytest
import torch

from torchgwas import scan_gpu

triton_scan = pytest.importorskip('torchgwas.triton_scan')
pytestmark = pytest.mark.skipif(not triton_scan.available(), reason='Triton kernels need a CUDA device')
DEVICE = 'cuda'


def _pack(codes, pgen=True):
    """(rows, n) two-bit codes to 64-byte padded packed rows, sample s at byte s//4, bits 2*(s%4)."""
    rows, n = codes.shape
    width = ((n + 3) // 4 + 63) // 64 * 64
    padded = np.zeros((rows, width * 4), dtype=np.uint8)
    padded[:, :n] = codes
    quads = padded.reshape(rows, width, 4).astype(np.uint16)
    return torch.as_tensor((quads[..., 0] | quads[..., 1] << 2 | quads[..., 2] << 4 | quads[..., 3] << 6)
                           .astype(np.uint8), device=DEVICE)


def _reference(values):
    """FP64 prepare from dosages (NaN missing): centred, ss, min, max, present."""
    observed = ~np.isnan(values)
    present = observed.sum(1)
    filled = np.where(observed, values, 0.0)
    mean = np.where(present > 0, filled.sum(1) / np.maximum(present, 1), 0.0).astype(np.float32).astype(np.float64)
    centred = np.where(observed, values - mean[:, None], 0.0)
    low = np.where(present > 0, np.where(observed, values, np.inf).min(1), np.nan)
    high = np.where(present > 0, np.where(observed, values, -np.inf).max(1), np.nan)
    return centred, (centred * centred).sum(1), low, high, present


def _dosages(rng, rows, n):
    values = rng.integers(0, 3, size=(rows, n)).astype(np.float64)
    values[rng.random(values.shape) < 0.05] = np.nan
    values[0] = np.nan          # nothing observed
    values[1] = 1.0             # invariant
    values[2, :5] = np.nan
    return values


def _check_prepare(got, values):
    centred, ss, low, high, present = _reference(values)
    np.testing.assert_allclose(got[0].double().cpu().numpy(), centred, rtol=0, atol=2e-6)
    np.testing.assert_allclose(got[1].double().cpu().numpy(), ss, rtol=2e-6)
    np.testing.assert_array_equal(got[2].cpu().numpy(), low.astype(np.float32))
    np.testing.assert_array_equal(got[3].cpu().numpy(), high.astype(np.float32))
    np.testing.assert_array_equal(got[4].cpu().numpy(), present)


@pytest.mark.parametrize('n', [37, 2500])
def test_prepare_decodes_every_encoding(n):
    rng = np.random.default_rng(n)
    values = _dosages(rng, 9, n)
    missing = np.isnan(values)
    as_int8 = torch.as_tensor(np.where(missing, -9, values).astype(np.int8), device=DEVICE)
    _check_prepare(triton_scan.prepare(as_int8), values)
    # uint8 codes with a scale, with and without a sentinel.
    scaled = torch.as_tensor(np.where(missing, 255, values * 127).astype(np.uint8), device=DEVICE)
    _check_prepare(triton_scan.prepare(scaled, 127.0, 255), values)
    complete = np.where(missing, 1.0, values)
    _check_prepare(triton_scan.prepare(torch.as_tensor((complete * 127).astype(np.uint8), device=DEVICE), 127.0),
                   complete)
    # float32: NaN, or an explicit sentinel.
    _check_prepare(triton_scan.prepare(torch.as_tensor(values.astype(np.float32), device=DEVICE)), values)
    _check_prepare(triton_scan.prepare(torch.as_tensor(np.where(missing, -9, values).astype(np.float32),
                                                       device=DEVICE), 1.0, -9.0), values)
    # Packed PGEN (3 missing) and PLINK1 (01 missing; codes 00, 10, 11 are 0, 1, 2).
    pgen = np.where(missing, 3, values).astype(np.uint8)
    _check_prepare(triton_scan.prepare(_pack(pgen), encoding='pgen_2bit', n_samples=n), values)
    plink = np.where(missing, 1, np.choose(np.nan_to_num(values).astype(int), [0, 2, 3])).astype(np.uint8)
    _check_prepare(triton_scan.prepare(_pack(plink), encoding='plink_2bit', n_samples=n), values)


@pytest.mark.parametrize('traits,covariates', [(5, 0), (1500, 3)])
def test_finish_is_scan_gpus_contract(traits, covariates):
    rng = np.random.default_rng(traits)
    n = 300
    values = _dosages(rng, 12, n)
    values[3, 7] = 2.5  # out of the dosage range, for validate_range
    raw = torch.as_tensor(values.astype(np.float32), device=DEVICE)
    centred, ss, low, high, present = triton_scan.prepare(raw)
    design = torch.as_tensor(rng.normal(size=(n, traits + covariates)) / n ** 0.5, dtype=torch.float32,
                             device=DEVICE)
    products = centred @ design
    phenotype_ss = (design[:, :traits] ** 2).sum(0) + 1
    offset = -float(covariates) - 2
    for validate in (False, True):
        beta, t, status = triton_scan.finish(products, ss, low, high, phenotype_ss, present, offset, validate)
        # The reference: finish_kernel's formula in FP64.
        p, s2 = products.double().cpu().numpy(), ss.double().cpu().numpy()
        residual = s2 - (p[:, traits:] ** 2).sum(1)
        df = present.cpu().numpy() + offset
        lo, hi = low.cpu().numpy(), high.cpu().numpy()
        valid = (residual > 1e-12) & (hi > lo) & (df > 0)
        want_status = np.where(np.isfinite(residual), np.where(valid, 0, 2), 1)
        if validate:
            want_status = np.where((lo < 0) | (hi > 2), 3, want_status)
        np.testing.assert_array_equal(status.cpu().numpy(), want_status)
        safe = np.maximum(residual, 1e-12)[:, None]
        gy = p[:, :traits]
        want_beta = gy / safe
        yss = np.maximum(phenotype_ss.double().cpu().numpy()[None, :] - gy * gy / safe, 1e-12)
        want_t = want_beta / np.sqrt(yss / np.maximum(df, 1)[:, None] / safe)
        want_beta[~valid], want_t[~valid] = 0, 0
        np.testing.assert_allclose(beta.cpu().numpy(), want_beta, rtol=2e-5, atol=1e-6)
        np.testing.assert_allclose(t.cpu().numpy(), want_t, rtol=2e-5, atol=1e-5)
    if scan_gpu.available():
        native = scan_gpu.finish(products, ss, low, high, phenotype_ss, present, offset, True)
        ours = triton_scan.finish(products, ss, low, high, phenotype_ss, present, offset, True)
        for a, b in zip(ours, native):
            np.testing.assert_allclose(a.float().cpu().numpy(), b.float().cpu().numpy(), rtol=1e-6, atol=1e-6)


def test_finish_min_p_ranks_as_variant_reduction():
    from torchgwas.reduce import VariantReduction
    rng = np.random.default_rng(5)
    n, traits = 200, 1500
    values = _dosages(rng, 16, n)
    raw = torch.as_tensor(values.astype(np.float32), device=DEVICE)
    centred, ss, low, high, present = triton_scan.prepare(raw)
    design = torch.as_tensor(rng.normal(size=(n, traits + 2)) / n ** 0.5, dtype=torch.float32, device=DEVICE)
    # A tie across trait blocks (columns 3 and 1203) and a NaN trait.
    design[:, 1203] = design[:, 3]
    products = centred @ design
    phenotype_ss = (design[:, :traits] ** 2).sum(0) + 1
    phenotype_ss[10] = float('nan')
    beta, t, status = triton_scan.finish(products, ss, low, high, phenotype_ss, present, -4.0)
    want = VariantReduction('min-p').reduce(beta, t, status, present.float() - 4, 1)
    got = triton_scan.finish_min_p(products, ss, low, high, phenotype_ss, present, -4.0)
    np.testing.assert_array_equal(got[2].cpu().numpy(), want[2].cpu().numpy())
    np.testing.assert_array_equal(got[3].cpu().numpy(), want[3].cpu().numpy())
    np.testing.assert_allclose(got[1].cpu().numpy(), want[1].cpu().numpy(), rtol=1e-6, equal_nan=True)
    np.testing.assert_allclose(got[0].cpu().numpy(), want[0].cpu().numpy(), rtol=1e-6, equal_nan=True)


def test_backend_resolution(monkeypatch):
    for key in ('TORCHGWAS_STATS_BACKEND', 'TORCHGWAS_NATIVE_STATS'):
        monkeypatch.delenv(key, raising=False)
    assert scan_gpu.resolve_statistics_backend() == 'triton'
    # An explicit TORCHGWAS_NATIVE_STATS=0 keeps the Torch statistics plans were calibrated on.
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS', '0')
    assert scan_gpu.resolve_statistics_backend() == 'torch'
    monkeypatch.setenv('TORCHGWAS_STATS_BACKEND', 'triton')
    assert scan_gpu.resolve_statistics_backend() == 'triton'
    monkeypatch.setenv('TORCHGWAS_STATS_BACKEND', 'torch')
    monkeypatch.delenv('TORCHGWAS_NATIVE_STATS')
    assert scan_gpu.resolve_statistics_backend() == 'torch'
    monkeypatch.setenv('TORCHGWAS_STATS_BACKEND', 'bogus')
    with pytest.raises(ValueError):
        scan_gpu.resolve_statistics_backend()


@pytest.mark.parametrize('fmt', ['pgen', 'bed'])
@pytest.mark.parametrize('reduce', [None, 'min-p'])
def test_the_default_scan_matches_torch_statistics(tmp_path, monkeypatch, fmt, reduce):
    from test_drop_subject import COMMON, _genotype, _panel
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats import open_binary_sumstats
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    calls, covariates, phenotype = _panel()
    phenotype = np.nan_to_num(phenotype)
    path, load = _genotype(tmp_path, 'g', fmt, calls, None)
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    monkeypatch.delenv('TORCHGWAS_PGEN_PACKED', raising=False)
    outputs = {}
    monkeypatch.delenv('TORCHGWAS_NATIVE_STATS', raising=False)
    monkeypatch.delenv('TORCHGWAS_BED_NATIVE', raising=False)
    for backend in ('triton', 'torch'):
        monkeypatch.setenv('TORCHGWAS_STATS_BACKEND', backend)
        run_linear_gwas(path, phenotype, covariates, device='cuda:0', compute_dtype='float32',
                        output_dir=tmp_path / backend, reduce=reduce, **load, **COMMON)
        if reduce is None:
            outputs[backend] = [np.asarray(a) for a in open_binary_sumstats(tmp_path / backend / 'sumstats')[:3]]
        else:
            _, parts = open_indexed_sumstats(tmp_path / backend / 'sumstats')
            parts = list(parts)
            outputs[backend] = [np.concatenate([p[key] for p in parts])
                                for key in ('variant_index', 'trait_index', 't_stat', 'neg_log10_p')]
    for got, want in zip(outputs['triton'], outputs['torch']):
        if got.dtype.kind in 'iu':
            np.testing.assert_array_equal(got, want)
        else:
            np.testing.assert_allclose(got, want, rtol=3e-5, atol=3e-6, equal_nan=True)


def test_packed_rows_on_a_device_without_triton_are_unpacked_by_torch(tmp_path, monkeypatch):
    # The source expects fused statistics and sends packed rows; the scan's
    # device cannot run the kernels, so Torch unpacks and runs the statistics.
    from test_drop_subject import COMMON, _genotype, _panel
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats import open_binary_sumstats
    calls, covariates, phenotype = _panel()
    phenotype = np.nan_to_num(phenotype)
    path, load = _genotype(tmp_path, 'g', 'pgen', calls, None)
    for key in ('TORCHGWAS_STATS_BACKEND', 'TORCHGWAS_NATIVE_STATS', 'TORCHGWAS_PGEN_PACKED'):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    monkeypatch.setattr(triton_scan, 'available', lambda device=None: False)
    options = dict(device='cuda:0', compute_dtype='float32', **load, **COMMON)
    run_linear_gwas(path, phenotype, covariates, output_dir=tmp_path / 'unpacked', **options)
    monkeypatch.setenv('TORCHGWAS_STATS_BACKEND', 'torch')
    run_linear_gwas(path, phenotype, covariates, output_dir=tmp_path / 'torch', **options)
    got, want = open_binary_sumstats(tmp_path / 'unpacked' / 'sumstats'), open_binary_sumstats(tmp_path / 'torch' / 'sumstats')
    for index in range(3):
        np.testing.assert_allclose(np.asarray(got[index]), np.asarray(want[index]), rtol=1e-6, atol=1e-7,
                                   equal_nan=True)
