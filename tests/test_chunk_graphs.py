"""TORCHGWAS_CHUNK_GRAPHS: full chunks replayed as CUDA graphs give the eager scan's output.

The graphs hold the same kernels the eager loop launches, so every stored
value must match bit for bit: dense output and fused min-p, with missing
calls (a df per variant), a sample selection (the packed missing mask), one
device and variant shards. The short last chunk stays eager.
"""
import numpy as np
import pytest
import torch

from torchgwas.api import run_linear_gwas
from torchgwas.sumstats import open_binary_sumstats
from test_min_p import _read
from test_pgen_native_reader import write_pgen

N, M, K = 203, 1100, 37
CHUNK = 128  # eight full chunks, then a short one of 76

pytestmark = pytest.mark.skipif(torch.cuda.device_count() < 1, reason='CUDA required')


def _inputs(tmp_path, drop_subjects):
    rng = np.random.default_rng(20261006)
    calls = rng.integers(0, 3, size=(M, N), dtype=np.uint8)
    calls[rng.random(calls.shape) < 0.03] = 3
    path = tmp_path / 'input.pgen'
    write_pgen(path, calls)
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(N)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'
                                         + ''.join(f'1\t{i + 1}\tv{i}\tA\tC\n' for i in range(M)))
    y = rng.normal(size=(N, K)).astype(np.float32)
    y[:, 4] += 0.7 * np.where(calls[9] == 3, 0, calls[9])
    if drop_subjects:
        y[[3, 50, 101], [0, 5, 6]] = np.nan
    return path, y, rng.normal(size=(N, 3)).astype(np.float32)


def _layouts():
    count = torch.cuda.device_count()
    first = 1 if count > 1 else 0
    layouts = [dict(device=f'cuda:{first}')]
    if count - first >= 2:
        layouts.append(dict(variant_devices=[f'cuda:{first}', f'cuda:{first + 1}']))
    return layouts


def _run(tmp_path, monkeypatch, name, graphs, path, y, covariates, layout, mode, drop_subjects):
    replays = []

    class CountingGraph(torch.cuda.CUDAGraph):
        def replay(self):
            replays.append(1)
            super().replay()

    monkeypatch.setattr(torch.cuda, 'CUDAGraph', CountingGraph)
    monkeypatch.setenv('TORCHGWAS_CHUNK_GRAPHS', '1' if graphs else '0')
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    run_linear_gwas(path, y, covariates, output_dir=tmp_path / name, genotype_format='pgen',
                    pgen_mode='hardcall', compute_dtype='float32', chunk_size=CHUNK, prefetch_chunks=3, reader_workers=2,
                    missing_phenotype='drop_subject' if drop_subjects else 'impute',
                    **(dict(reduce='min-p') if mode == 'min-p' else {}), **layout)
    return len(replays)


@pytest.mark.parametrize('drop_subjects', [False, True])
@pytest.mark.parametrize('mode', ['dense', 'min-p'])
@pytest.mark.parametrize('layout', _layouts(), ids=lambda layout: 'shards' if 'variant_devices' in layout else 'one')
def test_graph_replay_matches_eager(tmp_path, monkeypatch, layout, mode, drop_subjects):
    path, y, covariates = _inputs(tmp_path, drop_subjects)
    eager = _run(tmp_path, monkeypatch, 'eager', False, path, y, covariates, layout, mode, drop_subjects)
    graphed = _run(tmp_path, monkeypatch, 'graphs', True, path, y, covariates, layout, mode, drop_subjects)
    assert eager == 0
    # Per device: the first full chunk warms up eagerly and the last is short;
    # every other chunk replays one graph.
    shards = len(layout.get('variant_devices', [None]))
    assert graphed >= M // CHUNK - shards
    if mode == 'dense':
        want, got = open_binary_sumstats(tmp_path / 'eager/sumstats'), open_binary_sumstats(tmp_path / 'graphs/sumstats')
        for index in range(3):
            np.testing.assert_array_equal(np.asarray(got[index]), np.asarray(want[index]))
        assert got[3]['df'] == want[3]['df']
    else:
        _assert_same_min_p(tmp_path)


def _assert_same_min_p(tmp_path):
    want, got = _read(tmp_path / 'eager/sumstats')[1], _read(tmp_path / 'graphs/sumstats')[1]
    assert want.keys() == got.keys()
    for key in want:
        np.testing.assert_array_equal(got[key], want[key])


def test_graph_replay_under_the_chunk_tuner(tmp_path, monkeypatch):
    # Autotune observes every chunk. With one candidate size the tuner is
    # settled from the start, so graphs replay and both runs use one size.
    path, y, covariates = _inputs(tmp_path, False)
    tuned = dict(autotune=True, autotune_options=dict(devices=[_layouts()[0]['device']], cpus_per_device=1))
    eager = _run(tmp_path, monkeypatch, 'eager', False, path, y, covariates, tuned, 'min-p', False)
    graphed = _run(tmp_path, monkeypatch, 'graphs', True, path, y, covariates, tuned, 'min-p', False)
    assert eager == 0 and graphed > 0
    _assert_same_min_p(tmp_path)
