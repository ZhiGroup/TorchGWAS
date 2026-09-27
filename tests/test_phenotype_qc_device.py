"""GPU phenotype column QC makes exactly the NumPy decisions."""
import numpy as np
import pytest
import torch

from torchgwas.preprocess import _phenotype_column_mask, _phenotype_qc, prepare_inputs_for_prep


def panel(seed=5, n=301, k=97, dtype=np.float32):
    rng = np.random.default_rng(seed)
    y = rng.normal(size=(n, k)).astype(dtype)
    y[rng.random((n, k)) < .1] = np.nan
    y[:, 3] = 2.5                      # constant
    y[:, 7] = np.nan                   # all missing
    y[1:, 11] = np.nan                 # one observation
    y[:, 13] = 1e6; y[0, 13] = 1e6+1   # near-constant with a large mean
    y[::2, 17] = np.nan; y[1::2, 17] = 4.0  # constant where observed
    return y


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('block_bytes', [512 << 20, 301*4*10])  # one block, and 10-column blocks
def test_device_qc_matches_numpy(dtype, block_bytes):
    if not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    y = panel(dtype=dtype)
    keep, counts = _phenotype_column_mask(y)
    got_keep, got_counts, missing = _phenotype_qc(y, dtype=dtype, device='cuda:0', block_bytes=block_bytes)
    np.testing.assert_array_equal(got_keep, keep)
    np.testing.assert_array_equal(got_counts, counts)
    assert missing == int(np.isnan(y).sum())
    assert not got_keep[[3, 7, 11, 17]].any() and got_keep[13]


def test_device_qc_refuses_infinities_and_matches_through_prepare():
    if not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    y = panel()
    genotype = np.zeros((y.shape[0], 5), dtype=np.float32)
    genotype[0] = 1  # nonconstant variants
    host = prepare_inputs_for_prep(genotype, y, dtype=np.float32, phenotype_block_size=32, validate_genotype=False)
    device = prepare_inputs_for_prep(genotype, y, dtype=np.float32, phenotype_block_size=32,
                                     validate_genotype=False, qc_device='cuda:0')
    assert host[2] == device[2]
    np.testing.assert_array_equal(np.asarray(host[0]), np.asarray(device[0]))
    y[5, 2] = np.inf
    with pytest.raises(ValueError, match='infinite'):
        prepare_inputs_for_prep(genotype, y, dtype=np.float32, validate_genotype=False, qc_device='cuda:0')
