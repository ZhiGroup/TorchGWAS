"""trait_workers='process': one child process per trait device writes the store one process's threads write."""
import json

import numpy as np
import pytest
import torch

from torchgwas import trait_processes
from torchgwas.api import run_linear_gwas


def _inputs(tmp_path, n=203, m=37, k=13):
    from test_pgen_native_reader import write_pgen
    rng = np.random.default_rng(4)
    calls = rng.integers(0, 3, size=(n, m))
    y = rng.normal(size=(n, k))
    # Strong effects in both halves of the traits, so the default Bonferroni
    # threshold (5e-8 over every kept trait) selects pairs on each device.
    for variant, trait in ((3, 1), (10, 2), (20, 8), (30, 11), (5, 12)):
        y[:, trait] += 1.2 * (calls[:, variant] - calls[:, variant].mean())
    y[rng.random(y.shape) < 0.05] = np.nan
    y[:, 6] = 1.0                                   # constant: QC drops it, so kept != input columns
    covariates = rng.normal(size=(n, 3))
    path = tmp_path / 'input.pgen'
    write_pgen(path, calls.T.astype(np.uint8))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    return path, y.astype(np.float32), covariates.astype(np.float32)


def _run(path, y, covariates, out, workers, **extra):
    return run_linear_gwas(path, y, covariates, output_dir=out, genotype_format='pgen', pgen_mode='hardcall',
                           compute_dtype='float32', chunk_size=8, reader_workers=4, prefetch_chunks=2,
                           trait_devices=['cuda:0', 'cuda:1'], trait_block=4, trait_workers=workers, **extra)


def _pairs(store):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    manifest, parts = open_indexed_sumstats(store)
    parts = list(parts)
    return manifest, {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}


def _same_run(thread_out, process_out):
    qc = [json.loads((out / 'qc.json').read_text()) for out in (thread_out, process_out)]
    for key in ('phenotype_columns_input', 'phenotype_columns_kept', 'phenotype_kept_column_indices',
                'phenotype_observed_counts', 'phenotype_missing_cells', 'n_samples'):
        assert qc[0][key] == qc[1][key], key
    run = [json.loads((out / 'run.json').read_text()) for out in (thread_out, process_out)]
    assert run[0]['trait_columns'] == run[1]['trait_columns']
    assert run[0]['significance_threshold'] == run[1]['significance_threshold']
    assert run[1]['trait_workers'] == 'process' and len(run[1]['trait_processes']) == 2
    assert not (process_out / trait_processes.WORK_DIRECTORY).exists()


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason='two CUDA devices required')
@pytest.mark.parametrize('missing', ['exact', 'impute'])
def test_significant_pairs_match_the_thread_run(tmp_path, missing):
    path, y, covariates = _inputs(tmp_path)
    for workers in ('thread', 'process'):
        _run(path, y, covariates, tmp_path / workers, workers, reduce='significant', missing_phenotype=missing)
    _same_run(tmp_path / 'thread', tmp_path / 'process')
    want_manifest, want = _pairs(tmp_path / 'thread' / 'sumstats')
    got_manifest, got = _pairs(tmp_path / 'process' / 'sumstats')
    assert got_manifest['traits'] == want_manifest['traits'] and got_manifest['shape'] == want_manifest['shape']
    assert got_manifest.get('row_order') == want_manifest.get('row_order')
    # Pairs on both devices, chosen at one threshold over all 12 kept traits.
    assert len(want['variant_index']) >= 5 and len(set(want['trait_index'] < 6)) == 2
    for key in ('variant_index', 'trait_index', 'df'):
        np.testing.assert_array_equal(got[key], want[key], err_msg=key)
    for key in ('beta', 't_stat', 'neg_log10_p'):
        np.testing.assert_allclose(got[key], want[key], rtol=1e-5, err_msg=key)


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason='two CUDA devices required')
def test_min_p_winners_match_the_thread_run(tmp_path):
    path, y, covariates = _inputs(tmp_path)
    for workers in ('thread', 'process'):
        _run(path, y, covariates, tmp_path / workers, workers, reduce='min-p')
    _same_run(tmp_path / 'thread', tmp_path / 'process')
    _, want = _pairs(tmp_path / 'thread' / 'sumstats')
    _, got = _pairs(tmp_path / 'process' / 'sumstats')
    order = np.argsort(want['variant_index'])
    np.testing.assert_array_equal(got['variant_index'], want['variant_index'][order])
    np.testing.assert_array_equal(got['trait_index'], want['trait_index'][order])
    assert len(set(got['trait_index'] < 6)) == 2                 # winners from both partitions
    for key in ('t_stat', 'neg_log10_p', 'df'):
        np.testing.assert_allclose(got[key], want[key][order], rtol=1e-5, err_msg=key)


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason='two CUDA devices required')
@pytest.mark.parametrize('missing', [True, False])
def test_full_output_tiles_match_the_thread_run(tmp_path, missing):
    from torchgwas.sumstats import open_binary_df, open_binary_sumstats
    path, y, covariates = _inputs(tmp_path)
    if not missing:
        y = np.nan_to_num(y)
    for workers in ('thread', 'process'):
        _run(path, y, covariates, tmp_path / workers, workers)
    _same_run(tmp_path / 'thread', tmp_path / 'process')
    want, got = (open_binary_sumstats(tmp_path / workers / 'sumstats') for workers in ('thread', 'process'))
    for a, b in zip(got[:3], want[:3]):
        np.testing.assert_allclose(np.asarray(a), np.asarray(b), rtol=1e-5)
    manifest = json.loads((tmp_path / 'process' / 'sumstats' / 'manifest.json').read_text())
    assert manifest['layout'] == 'trait_tiles' and manifest['traits'] == want[3]['traits']
    if missing:
        assert manifest['df'] == want[3]['df']
    else:
        np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path / 'process' / 'sumstats')),
                                      np.asarray(open_binary_df(tmp_path / 'thread' / 'sumstats')))


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason='two CUDA devices required')
def test_a_phenotype_table_is_split_by_its_columns(tmp_path):
    path, y, covariates = _inputs(tmp_path)
    table = tmp_path / 'phenotype.tsv'
    table.write_text('IID\t' + '\t'.join(f'p{j}' for j in range(y.shape[1])) + '\n' + ''.join(
        f's{i}\t' + '\t'.join('nan' if np.isnan(v) else repr(float(v)) for v in row) + '\n' for i, row in enumerate(y)))
    for workers in ('thread', 'process'):
        run_linear_gwas(path, None, covariates, phenotype_table=table, output_dir=tmp_path / workers,
                        genotype_format='pgen', pgen_mode='hardcall', compute_dtype='float32', chunk_size=8,
                        reader_workers=4, trait_devices=['cuda:0', 'cuda:1'], trait_block=4, trait_workers=workers,
                        reduce='significant')
    _same_run(tmp_path / 'thread', tmp_path / 'process')
    want_manifest, want = _pairs(tmp_path / 'thread' / 'sumstats')
    got_manifest, got = _pairs(tmp_path / 'process' / 'sumstats')
    assert got_manifest['traits'] == want_manifest['traits'] and 'p6' not in got_manifest['traits']
    for key in ('variant_index', 'trait_index', 'df'):
        np.testing.assert_array_equal(got[key], want[key], err_msg=key)


def test_a_failed_child_stops_the_run_with_its_log(tmp_path):
    y = np.random.default_rng(0).normal(size=(20, 4)).astype(np.float32)
    with pytest.raises(RuntimeError, match='trait process .* exited with'):
        run_linear_gwas(tmp_path / 'missing.pgen', y, None, output_dir=tmp_path / 'out', genotype_format='pgen',
                        trait_devices=['cuda:0', 'cuda:1'], trait_workers='process', reduce='significant')
    logs = sorted((tmp_path / 'out' / trait_processes.WORK_DIRECTORY).glob('part_*.log'))
    assert len(logs) == 2 and any('Error' in log.read_text() for log in logs)    # kept for diagnosis


@pytest.mark.parametrize('extra, message', [
    (dict(missing_phenotype='drop_subject'), 'drop_subject'),
    (dict(reduce='jagwas'), 'jagwas'),
    (dict(output_dir=None), 'output_dir'),
    (dict(variant_devices=['cuda:2']), 'variant_devices'),
])
def test_unsupported_runs_are_refused(tmp_path, extra, message):
    arguments = dict(output_dir=tmp_path / 'out', trait_devices=['cuda:0', 'cuda:1'], trait_workers='process')
    arguments.update(extra)
    with pytest.raises(ValueError, match=message):
        run_linear_gwas(tmp_path / 'input.pgen', np.zeros((4, 2)), None, **arguments)


def test_each_child_sees_its_own_card_as_cuda0(monkeypatch):
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '2,6,1')
    assert [trait_processes._physical(f'cuda:{i}') for i in range(3)] == ['2', '6', '1']
    monkeypatch.delenv('CUDA_VISIBLE_DEVICES')
    assert trait_processes._physical('cuda:5') == '5'
    assert trait_processes._partitions(13, 2) == [(0, 7), (7, 13)]
    assert trait_processes._partitions(2, 4) == [(0, 1), (1, 2)]
