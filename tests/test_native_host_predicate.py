"""Exact CPU predicate semantics, immutable builds and calculator price binding."""
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
import numpy as np
import pytest
from torchgwas import native_host_predicate as native
from torchgwas.host_significance import select_host_pairs, ceil_float32, host_selector


@pytest.mark.parametrize('layout', ['scalar', 'row', 'reversed_row', 'pair', 'trait'])
@pytest.mark.parametrize('beta', [True, False])
def test_native_and_fallback_paths_have_exact_owned_results(monkeypatch, layout, beta):
    rng = np.random.default_rng(92383)
    values = rng.standard_normal((17, 29), dtype=np.float32)
    values.flat[::11] = np.nan; values.flat[1::11] = np.inf
    critical = dict(scalar=np.array(1.), row=np.linspace(-2., 2., 17)[:, None],
        reversed_row=np.linspace(-2., 2., 17)[::-1, None], pair=rng.random(values.shape),
        trait=np.linspace(-1., 1., 29)[None, :])[layout]
    df = np.arange(17, dtype=np.float32)[:, None]+20
    b = values if beta else None
    monkeypatch.setenv('TORCHGWAS_HOST_PREDICATE', 'numpy')
    expected = select_host_pairs(b, values, df, critical, 23)
    monkeypatch.setenv('TORCHGWAS_HOST_PREDICATE', 'native')
    with patch.object(native, 'fill_mask', wraps=native.fill_mask) as fill:
        actual = select_host_pairs(b, values, df, critical, 23)
        assert fill.call_count == 1
    values.fill(-99); df.fill(-99)
    for got, want in zip(actual, expected):
        if want is None: assert got is None
        else: np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize('case', ['float64', 'strided', 'byteswapped', 'unaligned', 'pair'])
def test_unsupported_layout_is_declined_before_native_load(case):
    values = np.ones((3, 5), np.float32)
    limits = np.broadcast_to(np.ones((3, 1), np.float32), values.shape)
    if case == 'float64': values = values.astype(np.float64)
    if case == 'byteswapped': values = values.astype('>f4')
    if case == 'strided': values = np.ones((3, 10), np.float32)[:, ::2]
    if case == 'unaligned': values = np.ndarray((3, 5), np.float32, buffer=bytearray(61), offset=1)
    if case == 'pair': limits = np.ones((3, 5), np.float32)
    with patch.object(native, 'library', side_effect=AssertionError('Unexpected load')):
        assert native.fill_mask(values, limits, np.empty(values.shape, bool)) is False


def test_unsafe_mask_rejected():
    values = np.ones((3, 5), np.float32)
    limits = np.broadcast_to(np.ones((3, 1), np.float32), values.shape)
    for out in (np.empty((5, 3), bool), values.view(bool).reshape(3, 20)[:, :5],
                np.zeros((3, 5), np.int8), np.empty((5, 3), bool).T):
        with pytest.raises(ValueError, match='separate writable'):
            native.fill_mask(values, limits, out)


def test_scalar_and_reversed_threshold_strides_and_simultaneous_calls():
    values = np.arange(161*31, dtype=np.float32).reshape(161, 31)
    critical = np.arange(161, dtype=np.float32)[::-1, None]+40
    limits = np.broadcast_to(critical[::-1], values.shape)
    wanted = np.isfinite(values) & (np.abs(values) >= limits)
    def fill(_):
        out = np.empty(values.shape, bool)
        assert native.fill_mask(values, limits, out)
        np.testing.assert_array_equal(out, wanted)
    with ThreadPoolExecutor(4) as pool: list(pool.map(fill, range(32)))


def test_immutable_build_reuse_and_loaded_identity():
    native.library()
    path = native._directory(native._identity()[2])
    before = {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in path.iterdir()}
    record = native.build()
    after = {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in path.iterdir()}
    assert before == after
    assert native.context()['binary_sha256'] == record['binary_sha256']


def test_changed_loaded_source_or_binary_rejected(monkeypatch):
    native.library(); source, identity, key = native._identity()
    with patch.object(native, '_identity', return_value=(source, identity, 'changed')):
        with pytest.raises(ValueError, match='source changed'): native.context()
    class ChangedDigest:
        def digests(self, paths): return ['changed']
    with pytest.raises(ValueError, match='binary changed'): native.context(digest_cache=ChangedDigest())


def test_missing_build_is_explicit_error(tmp_path, monkeypatch):
    monkeypatch.setattr(native, '_LOADED', None)
    monkeypatch.setattr(native, '_directory', lambda key: tmp_path/'absent')
    with pytest.raises(RuntimeError, match='--build'): native.library()


def test_invalid_backend_rejected(monkeypatch):
    monkeypatch.setenv('TORCHGWAS_HOST_PREDICATE', 'automatic')
    with pytest.raises(ValueError, match='numpy or native'): host_selector()


def test_native_ledger_and_legacy_price_rejection(monkeypatch):
    from test_significant_host_model import bank
    from torchgwas.significant_host_work import host_significant_selection_work, host_selection_service
    from torchgwas.significant_host_model import significant_host_runtime
    monkeypatch.setenv('TORCHGWAS_HOST_PREDICATE', 'native')
    work = host_significant_selection_work(1024, 8193, 0, return_beta=False)
    assert work['predicate_calls'] == 1
    assert work['predicate_logical_bytes'] == 5*1024*8193+4*1024
    assert work['predicate_primitive'] == 'predicate_native'
    with pytest.raises(ValueError, match='predicate limit'):
        significant_host_runtime({}, bank(), occupancy='empty', host_serial_fraction=0.)
    with pytest.raises(ValueError, match='predicate_native'):
        host_selection_service(work, bank()['prices'], cpu_fraction=1., dram_bytes_per_second=1e30, host_serial_fraction=0.)
    prices = bank()['prices']; prices['predicate_native'] = dict(call_cpu_seconds=.01,
        unit_cpu_seconds=1e-9, dram_bytes_per_unit=5.)
    steps = host_selection_service(work, prices, cpu_fraction=1., dram_bytes_per_second=1e30, host_serial_fraction=0.)
    assert steps[3]['seconds'] == .01+1024*8193*1e-9
    assert steps[3]['seconds']*steps[3]['resources']['dram'] == work['predicate_logical_bytes']
