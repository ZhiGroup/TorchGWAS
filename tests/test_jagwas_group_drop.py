"""drop_subject for JAGWAS groups: each group's joint test runs on its own subjects.

A subject missing (or outlier-masked in) one of a group's traits leaves that
group, not the others. Each group's statistic must be what a run of that
group alone gives on files without its dropped subjects. Groups that share
a trait share their drops.
"""
import numpy as np
import pytest
import torch

from torchgwas.api import run_linear_gwas
from torchgwas.sumstats_indexed import open_indexed_sumstats
from test_drop_subject import _genotype

N, M, K = 150, 17, 12
GROUPS = [('A', [0, 1, 2, 3, 4]), ('B', [5, 6, 7, 8, 9]), ('C', [8, 9, 10, 11])]
MISSING = [(7, 1), (50, 2), (20, 6), (33, 11)]  # (row, trait)
DROPPED = {'A': [7, 50], 'B': [20, 33], 'C': [20, 33]}  # B and C share traits 8 and 9


def _inputs(missing_calls):
    rng = np.random.default_rng(20260928)
    calls = rng.integers(0, 3, size=(N, M)).astype(np.uint8)
    calls[:3] = np.arange(3, dtype=np.uint8)[:, None]
    if missing_calls:
        calls[rng.random(calls.shape) < 0.03] = 3
    covariates = rng.normal(size=(N, 2))
    mixing = rng.normal(size=(K, K)) * 0.3 + np.eye(K)
    phenotype = (rng.normal(size=(N, K)) @ mixing + 0.3 * np.where(calls[:, :K] == 3, 1, calls[:, :K]))
    phenotype = phenotype.astype(np.float32)
    for row, trait in MISSING:
        phenotype[row, trait] = np.nan
    return calls, covariates, phenotype


def _chi2(directory):
    _manifest, parts = open_indexed_sumstats(directory / 'sumstats')
    parts = list(parts)
    index = np.concatenate([part['variant_index'] for part in parts])
    chi2 = np.concatenate([np.asarray(part['chi2']).reshape(len(part['variant_index']), -1) for part in parts])
    order = np.argsort(index)
    return index[order], chi2[order]


@pytest.mark.parametrize('device,missing_calls', [('cpu', False), ('cuda:0', False), ('cuda:0', True)])
def test_each_group_runs_on_its_own_subjects(tmp_path, monkeypatch, device, missing_calls):
    if device != 'cpu' and not torch.cuda.is_available():
        pytest.skip('CUDA required')
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    calls, covariates, phenotype = _inputs(missing_calls=missing_calls)
    options = dict(device=device, compute_dtype='float64' if device == 'cpu' else 'float32', chunk_size=4,
                   reader_workers=2, prefetch_chunks=2, reduce='jagwas')
    path, load = _genotype(tmp_path, 'full', 'pgen', calls, None)
    result = run_linear_gwas(path, phenotype, covariates, jagwas_groups=GROUPS, output_dir=tmp_path / 'grouped',
                             **load, **options)
    assert result.run_metadata['dropped_subjects_by_group'] == {name: len(rows) for name, rows in DROPPED.items()}
    index, chi2 = _chi2(tmp_path / 'grouped')
    tolerance = dict(rtol=1e-8, atol=1e-8) if device == 'cpu' else dict(rtol=2e-4, atol=2e-4)
    if missing_calls:
        # A call missing inside a group's subjects is imputed at the variant's
        # mean over every scanned subject, not over the group's own (the
        # genotype convention, complete_case.py): 0.3% here, 3% of calls
        # missing and 2 of 150 subjects dropped.
        tolerance = dict(rtol=1e-2, atol=1e-2)
    for position, (name, columns) in enumerate(GROUPS):
        kept = np.setdiff1d(np.arange(N), DROPPED[name])
        alone, _ = _genotype(tmp_path, f'alone_{name}', 'pgen', calls[kept], kept)
        run_linear_gwas(alone, phenotype[kept][:, columns], covariates[kept], output_dir=tmp_path / name,
                        **load, **options)
        want_index, want = _chi2(tmp_path / name)
        np.testing.assert_array_equal(index, want_index)
        np.testing.assert_allclose(chi2[:, position], want[:, 0], err_msg=name, **tolerance)


def test_outliers_leave_only_their_group(tmp_path, monkeypatch):
    from torchgwas.preprocess import mask_phenotype_outliers
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    calls, covariates, phenotype = _inputs(missing_calls=False)
    phenotype = np.nan_to_num(phenotype)
    phenotype[12, 3] = 60.0  # group A only
    masked, _ = mask_phenotype_outliers(phenotype, covariates, 5.0, whole_rows=False)
    flagged = {name: np.flatnonzero(np.isnan(masked[:, columns]).any(axis=1)).tolist() for name, columns in GROUPS}
    assert 12 in flagged['A'] and 12 not in flagged['B'] and 12 not in flagged['C']
    options = dict(device='cpu', compute_dtype='float64', chunk_size=4, reader_workers=2, prefetch_chunks=2,
                   reduce='jagwas')
    path, load = _genotype(tmp_path, 'full', 'pgen', calls, None)
    result = run_linear_gwas(path, phenotype, covariates, jagwas_groups=GROUPS, phenotype_outlier_sd=5.0,
                             output_dir=tmp_path / 'grouped', **load, **options)
    assert result.run_metadata['dropped_subjects_by_group']['A'] == len(flagged['A'])
    index, chi2 = _chi2(tmp_path / 'grouped')
    kept = np.setdiff1d(np.arange(N), flagged['A'])
    alone, _ = _genotype(tmp_path, 'alone', 'pgen', calls[kept], kept)
    run_linear_gwas(alone, phenotype[kept][:, GROUPS[0][1]], covariates[kept], output_dir=tmp_path / 'A',
                    **load, **options)
    np.testing.assert_allclose(chi2[:, 0], _chi2(tmp_path / 'A')[1][:, 0], rtol=1e-8, atol=1e-8)


def test_variant_shards_drop_by_group_as_one_device(tmp_path, monkeypatch):
    if not torch.cuda.is_available() or torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    calls, covariates, phenotype = _inputs(missing_calls=True)
    path, load = _genotype(tmp_path, 'full', 'pgen', calls, None)
    options = dict(compute_dtype='float32', chunk_size=4, reader_workers=2, prefetch_chunks=2, reduce='jagwas',
                   jagwas_groups=GROUPS, **load)
    run_linear_gwas(path, phenotype, covariates, device='cuda:0', output_dir=tmp_path / 'one', **options)
    run_linear_gwas(path, phenotype, covariates, variant_devices=['cuda:0', 'cuda:1'], output_dir=tmp_path / 'two',
                    **options)
    one, two = _chi2(tmp_path / 'one'), _chi2(tmp_path / 'two')
    np.testing.assert_array_equal(one[0], two[0])
    np.testing.assert_allclose(two[1], one[1], rtol=1e-5, atol=1e-5)
