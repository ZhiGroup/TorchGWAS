"""missing_phenotype='drop_subject': a run is the run on inputs without those samples.

A sample with any missing (or outlier-masked) phenotype value leaves the
analysis for every trait, so each backend must give what it gives when the
genotype, phenotype and covariate files simply do not contain that sample.
"""
import numpy as np
import pytest
import torch

from torchgwas.api import run_linear_gwas
from torchgwas.preprocess import mask_phenotype_outliers
from torchgwas.sumstats import open_binary_sumstats
from test_pgen_native_reader import write_pgen
from test_statistics import _write_bed

N, M, K = 83, 13, 5
# Three subjects miss one trait each and one misses its whole row.
MISSING = ([2, 17, 40, 61], [0, 3, 4, slice(None)])


def _panel(missing_calls=True, seed=20260928):
    rng = np.random.default_rng(seed)
    calls = rng.integers(0, 3, size=(N, M)).astype(np.uint8)
    calls[:3, :] = np.arange(3, dtype=np.uint8)[:, None]
    if missing_calls:
        # The CPU paths refuse missing calls, so only the CUDA paths carry them.
        calls[rng.random(calls.shape) < 0.04] = 3
    covariates = rng.normal(size=(N, 3))
    signal = np.where(calls[:, :K] == 3, 1, calls[:, :K])
    phenotype = (rng.normal(size=(N, K)) + 0.4 * signal).astype(np.float32)
    for row, column in zip(*MISSING):
        phenotype[row, column] = np.nan
    return calls, covariates, phenotype


def _write_pgen(path, calls, ids):
    ids = range(calls.shape[0]) if ids is None else ids
    write_pgen(path, calls.T)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
        f'1\t{i + 1}\tv{i}\tA\tC\n' for i in range(calls.shape[1])))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f'{i}\n' for i in ids))
    return path


def _dosage(calls):
    return np.where(calls == 3, np.nan, calls).astype(np.float64)


def _genotype(tmp_path, name, fmt, calls, ids):
    if fmt == 'pgen':
        return _write_pgen(tmp_path / f'{name}.pgen', calls, ids), dict(genotype_format='pgen', pgen_mode='hardcall')
    bed = _write_bed(tmp_path / name, _dosage(calls))
    if ids is not None:
        fam = bed.with_suffix('.fam')
        fam.write_text(''.join(f'F{i} I{i} 0 0 0 -9\n' for i in ids))
    return bed, dict(genotype_format='plink')


def _assert_same_output(actual_dir, expected_dir, rtol, atol):
    actual, expected = open_binary_sumstats(actual_dir / 'sumstats'), open_binary_sumstats(expected_dir / 'sumstats')
    for index in range(3):
        np.testing.assert_allclose(np.asarray(actual[index]), np.asarray(expected[index]),
                                   rtol=rtol, atol=atol, equal_nan=True)
    assert actual[3]['df'] == expected[3]['df']


BACKENDS = {
    'cpu': dict(device='cpu', compute_dtype='float64', env={}),
    # BED's Torch path gathers the selected samples; its native kernel (the
    # packed backend here) reads whole stored rows with the rest masked.
    'cuda-torch': dict(device='cuda:0', compute_dtype='float32',
                       env={'TORCHGWAS_NATIVE_STATS': '0', 'TORCHGWAS_BED_NATIVE': '0'}),
    'cuda-native-int8': dict(device='cuda:0', compute_dtype='float32',
                             env={'TORCHGWAS_NATIVE_STATS': '1', 'TORCHGWAS_PGEN_PACKED': '0'}),
    'cuda-native-packed': dict(device='cuda:0', compute_dtype='float32',
                               env={'TORCHGWAS_NATIVE_STATS': '1', 'TORCHGWAS_BED_NATIVE': '1'}),
}


def _backend(monkeypatch, name):
    backend = BACKENDS[name]
    if backend['device'].startswith('cuda'):
        if not torch.cuda.is_available():
            pytest.skip('CUDA required')
        if backend['env'].get('TORCHGWAS_NATIVE_STATS') == '1':
            from torchgwas import scan_gpu
            if not scan_gpu.available('cuda:0'):
                pytest.skip('native statistics kernels unavailable on this device')
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    for key, value in backend['env'].items():
        monkeypatch.setenv(key, value)
    if 'TORCHGWAS_PGEN_PACKED' not in backend['env']:
        monkeypatch.delenv('TORCHGWAS_PGEN_PACKED', raising=False)
    return dict(device=backend['device'], compute_dtype=backend['compute_dtype'])


def _tolerance(options):
    return (1e-10, 1e-10) if options['compute_dtype'] == 'float64' else (3e-5, 3e-6)


COMMON = dict(chunk_size=4, reader_workers=2, prefetch_chunks=2, sumstats_block_bytes=64)


@pytest.mark.parametrize('fmt', ['pgen', 'bed'])
@pytest.mark.parametrize('backend', list(BACKENDS))
def test_dropping_subjects_matches_inputs_without_them(tmp_path, monkeypatch, fmt, backend):
    if fmt == 'bed' and backend == 'cuda-native-int8':
        pytest.skip('int8 transport is PGEN-only')
    options = _backend(monkeypatch, backend)
    calls, covariates, phenotype = _panel(missing_calls=options['device'] != 'cpu')
    kept = np.setdiff1d(np.arange(N), MISSING[0])
    full, load = _genotype(tmp_path, 'full', fmt, calls, None)
    subset, _ = _genotype(tmp_path, 'subset', fmt, calls[kept], kept)
    result = run_linear_gwas(full, phenotype, covariates, output_dir=tmp_path / 'dropped', **load, **options, **COMMON)
    run_linear_gwas(subset, phenotype[kept], covariates[kept], output_dir=tmp_path / 'reference',
                    **load, **options, **COMMON)
    assert result.run_metadata['dropped_subjects'] == len(MISSING[0])
    assert result.run_metadata['missing_phenotype'] == 'drop_subject'
    _assert_same_output(tmp_path / 'dropped', tmp_path / 'reference', *_tolerance(options))


def test_in_memory_arrays_drop_subjects(tmp_path):
    calls, covariates, phenotype = _panel(missing_calls=False)
    kept = np.setdiff1d(np.arange(N), MISSING[0])
    genotype = _dosage(calls)
    options = dict(device='cpu', compute_dtype='float64', chunk_size=4)
    dropped = run_linear_gwas(genotype, phenotype, covariates, output_dir=tmp_path / 'dropped', **options)
    run_linear_gwas(genotype[kept], phenotype[kept], covariates[kept], output_dir=tmp_path / 'reference', **options)
    assert dropped.run_metadata['dropped_subjects'] == len(MISSING[0])
    _assert_same_output(tmp_path / 'dropped', tmp_path / 'reference', 1e-10, 1e-10)


def test_outlier_rows_leave_every_trait(tmp_path, monkeypatch):
    options = _backend(monkeypatch, 'cpu')
    calls, covariates, phenotype = _panel(missing_calls=False)
    phenotype = np.nan_to_num(phenotype)
    phenotype[[5, 30], [1, 2]] = 40.0
    masked, rows = mask_phenotype_outliers(phenotype, covariates, 4.0, whole_rows=True)
    assert rows[[5, 30]].all()
    kept = np.flatnonzero(~rows)
    full, load = _genotype(tmp_path, 'full', 'pgen', calls, None)
    subset, _ = _genotype(tmp_path, 'subset', 'pgen', calls[kept], kept)
    result = run_linear_gwas(full, phenotype, covariates, phenotype_outlier_sd=4.0,
                             output_dir=tmp_path / 'dropped', **load, **options, **COMMON)
    run_linear_gwas(subset, phenotype[kept], covariates[kept], output_dir=tmp_path / 'reference',
                    **load, **options, **COMMON)
    assert result.run_metadata['dropped_subjects'] == int(rows.sum())
    _assert_same_output(tmp_path / 'dropped', tmp_path / 'reference', 1e-10, 1e-10)


def test_complete_panel_is_untouched(tmp_path):
    calls, covariates, phenotype = _panel(missing_calls=False)
    phenotype = np.nan_to_num(phenotype)
    options = dict(device='cpu', compute_dtype='float64', chunk_size=4)
    dropped = run_linear_gwas(_dosage(calls), phenotype, covariates, output_dir=tmp_path / 'dropped', **options)
    imputed = run_linear_gwas(_dosage(calls), phenotype, covariates, missing_phenotype='impute',
                              output_dir=tmp_path / 'imputed', **options)
    assert dropped.run_metadata['dropped_subjects'] == 0
    _assert_same_output(tmp_path / 'dropped', tmp_path / 'imputed', 0, 0)


def test_impute_keeps_the_release_convention(tmp_path):
    calls, covariates, phenotype = _panel(missing_calls=False)
    options = dict(device='cpu', compute_dtype='float64', chunk_size=4)
    imputed = run_linear_gwas(_dosage(calls), phenotype, covariates, missing_phenotype='impute',
                              output_dir=tmp_path / 'imputed', **options)
    assert imputed.run_metadata['dropped_subjects'] == 0
    assert imputed.qc_summary['phenotype_missing_cells'] == 3 + K


def test_unknown_policy_is_refused():
    with pytest.raises(ValueError, match='missing_phenotype'):
        run_linear_gwas(np.zeros((4, 2)), np.zeros((4, 1)), missing_phenotype='drop')


@pytest.mark.parametrize('packed', [True, False])
def test_a_majority_selection_travels_as_whole_rows(tmp_path, monkeypatch, packed):
    from torchgwas.pgen import PgenGenotype
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS', '1' if packed else '0')
    monkeypatch.delenv('TORCHGWAS_PGEN_PACKED', raising=False)
    calls, _, _ = _panel()
    path = _write_pgen(tmp_path / 'g.pgen', calls, None)
    order = [i for i in range(N) if i not in (0, 2, 5, 9)]
    order[:3] = order[2::-1]
    source = PgenGenotype(path, mode='hardcall', selected_sample_ids=[str(i) for i in order])
    assert source.native_physical_samples == N
    assert source.native_sample_positions.tolist() == order
    assert source.native_unselected_samples.tolist() == [0, 2, 5, 9]
    if packed:
        assert source.native_encoding == 'pgen_2bit'
        assert source.native_row_width == ((N + 3) // 4 + 63) // 64 * 64
        pairs = np.unpackbits(source.native_packed_missing_mask, bitorder='little').reshape(-1, 2)[:N]
        assert np.flatnonzero(pairs.all(axis=1)).tolist() == [0, 2, 5, 9]
        assert (pairs.all(axis=1) | ~pairs.any(axis=1)).all()
    else:
        assert source.native_encoding == 'dosage' and source.native_row_width == N
    # A minority is decoded on the host; every sample again is the ordinary transport.
    source.select_samples([str(i) for i in range(30)])
    assert source.native_physical_samples is None and source.native_row_width == 30
    source.select_samples(None)
    assert source.native_physical_samples is None and source.native_packed_missing_mask is None


def test_reader_gathers_a_host_selection(tmp_path, monkeypatch):
    from torchgwas.pgen_native_reader import NativePgenReader
    calls, _, _ = _panel()
    path = _write_pgen(tmp_path / 'g.pgen', calls, None)
    subset = np.asarray([1, 4, 8, 30, 31, 60, 82])
    full = np.empty((M, N), np.int8)
    part = np.empty((M, subset.size), np.int8)
    with NativePgenReader(path, raw_sample_ct=N, variant_ct=M) as reader:
        reader.read_range(0, M, full)
    with NativePgenReader(path, raw_sample_ct=N, variant_ct=M, sample_subset=subset) as reader:
        reader.read_range(0, M, part)
        np.testing.assert_array_equal(part, full[:, subset])
        # A shorter read reuses the scratch the first one sized.
        reader.read_range(2, 5, part[:3])
        np.testing.assert_array_equal(part[:3], full[2:5, subset])


def test_packed_selection_in_another_order_matches_int8(tmp_path, monkeypatch):
    _backend(monkeypatch, 'cuda-native-packed')
    calls, covariates, phenotype = _panel()
    phenotype = np.nan_to_num(phenotype)
    path, load = _genotype(tmp_path, 'full', 'pgen', calls, None)
    order = np.random.default_rng(7).permutation(N)[:60]  # a majority, reordered
    options = dict(device='cuda:0', compute_dtype='float32', **load, **COMMON)
    run_linear_gwas(path, phenotype[order], covariates[order], sample_ids=[str(i) for i in order],
                    output_dir=tmp_path / 'packed', **options)
    monkeypatch.setenv('TORCHGWAS_PGEN_PACKED', '0')
    run_linear_gwas(path, phenotype[order], covariates[order], sample_ids=[str(i) for i in order],
                    output_dir=tmp_path / 'int8', **options)
    _assert_same_output(tmp_path / 'packed', tmp_path / 'int8', 3e-5, 3e-6)


@pytest.mark.parametrize('backend', ['cuda-torch', 'cuda-native-int8'])
def test_a_minority_selection_matches_its_own_file(tmp_path, monkeypatch, backend):
    options = _backend(monkeypatch, backend)
    calls, covariates, phenotype = _panel()
    phenotype = np.nan_to_num(phenotype)
    chosen = np.sort(np.random.default_rng(11).choice(N, 30, replace=False))
    full, load = _genotype(tmp_path, 'full', 'pgen', calls, None)
    subset, _ = _genotype(tmp_path, 'subset', 'pgen', calls[chosen], chosen)
    run_linear_gwas(full, phenotype[chosen], covariates[chosen], sample_ids=[str(i) for i in chosen],
                    output_dir=tmp_path / 'selected', **load, **options, **COMMON)
    run_linear_gwas(subset, phenotype[chosen], covariates[chosen], output_dir=tmp_path / 'reference',
                    **load, **options, **COMMON)
    _assert_same_output(tmp_path / 'selected', tmp_path / 'reference', 3e-5, 3e-6)
