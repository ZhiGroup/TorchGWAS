"""Borrowed native buffers and acknowledgement across asynchronous shards."""
import threading
import time
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import torch

from torchgwas.linear import linear_scan_multigpu, linear_scan_streaming_chunks, linear_scan_streaming, linear_scan


@pytest.mark.parametrize('ordered', [False, True])
@pytest.mark.parametrize('end_early', [False, True])
def test_borrowed_shard_waits_for_consumer_before_overwriting(ordered, end_early):
    closed = []

    def scan(source, y, c, variant_range, borrow_results, **kwargs):
        assert borrow_results

        def chunks():
            buffer = np.empty(8)
            try:
                for i in range(*variant_range):
                    buffer.fill(i)
                    yield i, i + 1, buffer, None, None
            finally:
                closed.append(variant_range)
        return chunks(), c

    def preprocess(y, c, **kwargs):
        return y, c, np.full(y.shape[1], y.shape[0])

    with patch('torchgwas.linear.choose_device', return_value=torch.device('cpu')), \
         patch('torchgwas.linear.residualize_and_standardize', side_effect=preprocess), \
         patch('torchgwas.linear.linear_scan_streaming_chunks', side_effect=scan):
        chunks, _ = linear_scan_multigpu(SimpleNamespace(shape=(10, 12)),
            np.ones((10, 1)), devices=['cuda:1', 'cuda:2'], chunk_size=1,
            ordered=ordered, result_queue_depth=4, borrow_results=True)
        observed = []
        try:
            for start, end, beta, _, _ in chunks:
                # Give the producer time to overwrite if acknowledgement is absent.
                time.sleep(.01)
                np.testing.assert_array_equal(beta, np.full(8, start))
                observed.append(start)
                if end_early:
                    break
        finally:
            chunks.close()
    assert len(closed) == 2
    assert not any(t.name.startswith('torchgwas-shard-') for t in threading.enumerate())
    if not end_early:
        assert sorted(observed) == list(range(12))
        if ordered:
            assert observed == list(range(12))


def test_borrowed_worker_failure_cancels_peers_waiting_for_acknowledgement():
    ready = threading.Event()

    def scan(source, y, c, variant_range, **kwargs):
        def chunks():
            if variant_range[0]:
                assert ready.wait(5)
                raise RuntimeError('borrowed reader failed')
            ready.set()
            for i in range(*variant_range):
                yield i, i + 1, np.zeros(4), None, None
        return chunks(), c

    with patch('torchgwas.linear.choose_device', return_value=torch.device('cpu')), \
         patch('torchgwas.linear.residualize_and_standardize',
               side_effect=lambda y, c, **kw: (y, c, np.full(y.shape[1], y.shape[0]))), \
         patch('torchgwas.linear.linear_scan_streaming_chunks', side_effect=scan):
        chunks, _ = linear_scan_multigpu(SimpleNamespace(shape=(10, 12)),
            np.ones((10, 1)), devices=['cuda:1', 'cuda:2'], chunk_size=1,
            ordered=False, borrow_results=True)
        with pytest.raises(RuntimeError, match='borrowed reader failed'):
            list(chunks)
    assert not any(t.name.startswith('torchgwas-shard-') for t in threading.enumerate())


@pytest.mark.parametrize('count,ordered', [(1, False), (2, False), (2, True)])
@pytest.mark.parametrize('missing', [False, True])
def test_native_borrowed_matches_owned_with_ring_reuse_and_partial_shards(count, ordered, missing, monkeypatch):
    if not torch.cuda.is_available() or torch.cuda.device_count() < 3:
        pytest.skip('CUDA devices 1 and 2 required')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS', '0')
    rng = np.random.default_rng(531)
    n, m, k = 257, 73, 9
    values = rng.integers(0, 3, (n, m)).astype(np.float32)
    values[:, 3] = 1
    values[::7, 8] = np.nan
    y, c = rng.normal(size=(n, k)), rng.normal(size=(n, 2))
    if missing:
        y[:11, 0] = np.nan
        y[::9, 2] = np.nan

    class Source:
        supports_fused_qc = True
        decode_workers = 1
        shape = (n, m)

        def iter_chunks(self, chunk_size, dtype=np.float32, variant_range=None, **kwargs):
            lo, hi = variant_range or (0, m)
            for start in range(lo, hi, chunk_size):
                end = min(start + chunk_size, hi)
                yield start, end, values[:, start:end].astype(dtype)

    def scan(borrow):
        return linear_scan_multigpu(Source(), y, c, chunk_size=7,
            devices=['cuda:1', 'cuda:2'][:count], reader_workers=2,
            prefetch_chunks=2, ordered=ordered, borrow_results=borrow)[0]

    # Deliberately retain the default owned outputs beyond every ring reuse.
    owned = list(scan(False))
    expected = {row[0]: row for row in owned}
    coverage = np.zeros(m, dtype=int)
    addresses = set()
    for row in scan(True):
        start, end, beta, stat, p = row
        assert not beta.flags.owndata and not stat.flags.owndata
        time.sleep(.002)
        for actual, reference in zip(row[2:], expected[start][2:]):
            np.testing.assert_allclose(actual, reference, rtol=2e-5, atol=3e-6, equal_nan=True)
        coverage[start:end] += 1
        addresses.add(beta.__array_interface__['data'][0])
    assert (coverage == 1).all()
    assert len(addresses) <= count * 2


def test_public_borrow_option_is_forwarded_to_native_dispatch():
    source = SimpleNamespace(shape=(12, 10), supports_fused_qc=True)
    with patch('torchgwas.linear.choose_device', return_value=torch.device('cuda:1')), \
         patch('torchgwas.native_scan.dosage_cuda_iterator', return_value=iter(())) as scan:
        chunks, _ = linear_scan_streaming_chunks(source, np.ones((12, 2)), None,
            chunk_size=3, device='cuda:1', already_processed=True,
            borrow_results=True)
        list(chunks)
    assert scan.call_args.kwargs['borrow_results'] is True


def test_cpu_streaming_fallback_uses_requested_device():
    rng = np.random.default_rng(723)
    values = rng.integers(0, 3, (31, 17)).astype(float)
    y, c = rng.normal(size=(31, 2)), rng.normal(size=(31, 2))

    class Source:
        shape = values.shape

        def iter_chunks(self, chunk_size, dtype, **kwargs):
            for start in range(0, self.shape[1], chunk_size):
                end = min(start + chunk_size, self.shape[1])
                yield start, end, values[:, start:end].astype(dtype)

    expected = linear_scan(values, y, c, device='cpu', compute_dtype='float64')
    actual = linear_scan_streaming(Source(), y, c, chunk_size=5,
        device='cpu', compute_dtype='float64')
    for result, reference in zip(actual[:3], expected[:3]):
        np.testing.assert_allclose(result, reference, atol=1e-10, rtol=1e-10)
