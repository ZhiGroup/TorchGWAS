"""Variant shards vs phenotype tiles, priced from work and transfer."""
import json

import numpy as np
import pytest
import torch

from torchgwas.layout_pricing import choose_split, split_costs

DEVICES = ['cuda:0', 'cuda:1', 'cuda:2', 'cuda:3']


def costs(**kwargs):
    base = dict(n_samples=35_000, n_traits=16_384, n_variants=10_000_000, genotype_bytes_per_variant=35_000,
                devices=DEVICES, h2d_bytes_per_second=50e9, peer_bytes_per_second={d: 300e9 for d in DEVICES[1:]},
                flops={d: 50e12 for d in DEVICES})
    base.update(kwargs)
    return split_costs(**base)


def test_ties_take_shards_and_slow_peer_copies_take_tiles():
    # H100-like: the genotype stream hides under the GEMM either way, a tie -> shards (measured faster).
    fast = costs()
    assert fast['phenotype_tiles']['stream'] < fast['phenotype_tiles']['gemm'] and choose_split(fast) == 'variants'
    # Shared uplink (A100 pair, 3.3 GB/s): every GPU streaming everything is the bottleneck -> shards.
    slow = costs(n_traits=4096, h2d_bytes_per_second=3.3e9, flops={d: 19.5e12 for d in DEVICES})
    assert slow['phenotype_tiles']['stream'] > slow['phenotype_tiles']['gemm'] and choose_split(slow) == 'variants'
    # A wide panel over a slow peer path makes the up-front copies decide -> tiles.
    wide = costs(n_traits=200_000, n_variants=200_000, h2d_bytes_per_second=12e9,
                 peer_bytes_per_second={d: 6e9 for d in DEVICES[1:]})
    assert wide['variant_shards']['upfront'] > wide['phenotype_tiles']['upfront'] and choose_split(wide) == 'traits'


@pytest.mark.parametrize('split', ['variants', 'auto'])
def test_autotuned_significant_split_matches_single_gpu(tmp_path, split, monkeypatch):
    if torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    from test_empirical_autotune import OPTIONS, RING, _fixture, _pairs, _run
    for key, value in dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0',
                           TORCHGWAS_NATIVE_STATS='0', TORCHGWAS_SCAN_PROFILE='0').items():
        monkeypatch.setenv(key, value)
    path = _fixture(tmp_path)
    common = dict(reduce='significant', significance_threshold=1e-3)
    _run(tmp_path, 'fixed', path, device='cuda:0', chunk_size=32, **common)
    _, _, run = _run(tmp_path, 'tuned', path, autotune=True,
                     autotune_options=dict(OPTIONS, devices=['cuda:0', 'cuda:1'], min_tile_traits=8,
                                           cpus_per_device=1, split=split), **RING, **common)
    layout = run['autotune']['layout']
    pricing = layout['split_pricing']
    assert pricing['requested'] == split and pricing['prices']['h2d_bytes_per_second'] > 0
    if split == 'variants':
        assert layout['variant_devices'] == ['cuda:0', 'cuda:1'] and layout['trait_block'] is None
    expected, actual = _pairs(tmp_path/'fixed'), _pairs(tmp_path/'tuned')
    assert expected and actual.keys() == expected.keys()
    np.testing.assert_allclose([actual[k] for k in expected], list(expected.values()), rtol=2e-4, atol=2e-4)
