"""Host genotype cache: later phenotype-tile rounds replay the first round's fills."""
import contextlib
import json

import numpy as np
import pytest
import torch

from torchgwas.genotype_cache import CachedFillSource, GenotypeFillCache


class FakeSource:
    allows_direct_native_fill = True
    native_dtype = np.int8

    def __init__(self, values):
        self.values, self.shape = values, (values.shape[1], values.shape[0])
        self.sessions = self.fills = 0

    @contextlib.contextmanager
    def native_reader_session(self):
        self.sessions += 1
        def fill(start, end, out):
            self.fills += 1
            np.copyto(out, self.values[start:end])
        yield fill


def test_second_pass_is_served_from_the_cache():
    values = np.random.default_rng(1).integers(0, 3, (50, 7), dtype=np.int8)  # variant-major rows
    source = FakeSource(values)
    cache = GenotypeFillCache(0, 50, 7, np.int8)
    view = CachedFillSource(source, cache)
    assert view.shape == source.shape  # everything else forwards
    for _ in range(3):
        out = np.empty((10, 7), np.int8)
        with view.native_reader_session() as fill:
            for start in range(0, 50, 10):
                fill(start, start + 10, out)
                np.testing.assert_array_equal(out, values[start:start + 10])
    assert source.fills == 5 and source.sessions == 1  # later passes never open the reader
    assert cache.audit()['hit_rows'] == 100 and cache.audit()['miss_rows'] == 50


@pytest.mark.parametrize('reduce', ['significant', None])
def test_tile_rounds_match_without_cache(tmp_path, reduce, monkeypatch):
    if not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    from test_empirical_autotune import _fixture, _pairs, _run
    from torchgwas.sumstats import open_binary_sumstats
    for key, value in dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0',
                           TORCHGWAS_NATIVE_STATS='0', TORCHGWAS_SCAN_PROFILE='0').items():
        monkeypatch.setenv(key, value)
    path = _fixture(tmp_path)  # 40 traits -> 5 tiles of 8 on one GPU: 5 rounds
    common = dict(device='cuda:0', chunk_size=64, trait_block=8, trait_devices=['cuda:0'])
    if reduce:
        common.update(reduce='significant', significance_threshold=1e-3)
    monkeypatch.setenv('TORCHGWAS_GENOTYPE_CACHE', '0')
    _run(tmp_path, 'plain', path, **common)
    monkeypatch.setenv('TORCHGWAS_GENOTYPE_CACHE', '1')
    _, _, run = _run(tmp_path, 'cached', path, **common)
    audit = run['genotype_cache']
    assert audit['miss_rows'] >= 8192 and audit['hit_rows'] >= 3 * 8192  # decoded once, replayed after
    if reduce:
        expected, actual = _pairs(tmp_path/'plain'), _pairs(tmp_path/'cached')
        assert expected and actual == expected
    else:
        a, b = open_binary_sumstats(tmp_path/'plain'/'sumstats'), open_binary_sumstats(tmp_path/'cached'/'sumstats')
        np.testing.assert_array_equal(np.asarray(a[1]), np.asarray(b[1]))
