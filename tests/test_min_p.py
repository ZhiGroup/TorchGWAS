"""reduce='min-p': one row per variant, the trait with the smallest exact p.

Every comparison is against the full scan it replaces: the dense store's t,
the df the scan used (per variant from genotype missingness, per pair once
phenotypes are missing) and the exact tail recomputed on the host.
"""
import json

import numpy as np
import pytest
import torch

from torchgwas.api import run_linear_gwas
from torchgwas.min_p import MinPReduction
from torchgwas.sumstats import open_binary_sumstats
from torchgwas.sumstats_indexed import open_indexed_sumstats
from torchgwas.tails import upper_tail_log10_from_t
from test_pgen_native_reader import write_pgen


def _inputs(tmp_path, *, missing_calls, missing_pheno, n=97, m=80, k=11, seed=5150):
    rng = np.random.default_rng(seed)
    calls = rng.integers(0, 3, size=(m, n), dtype=np.uint8)
    if missing_calls:
        # A different observed count, hence df, at every variant.
        for i in range(m):
            calls[i, :i] = 3
    path = tmp_path/'input.pgen'
    write_pgen(path, calls)
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'
                                         + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    y = rng.normal(size=(n, k)).astype(np.float32)
    # A planted effect, so some winners are real and some are noise.
    y[:, 3] += 0.8 * (calls[5].astype(np.float32) * (calls[5] != 3))
    if missing_pheno:
        y[:5, 0] = np.nan
        y[:9, 2] = np.nan
        # Heavily missing traits, whose pairs have far fewer df.
        y[:60, 7] = np.nan
        y[:75, 9] = np.nan
    covariates = rng.normal(size=(n, 2)).astype(np.float32)
    return path, y, covariates, calls


def _options(device, missing_phenotype='impute'):
    # A missing-phenotype convention that keeps the samples, which _expected recomputes.
    return dict(genotype_format='pgen', pgen_mode='hardcall', device=device, compute_dtype='float32',
                chunk_size=4, reader_workers=2, prefetch_chunks=2, sumstats_queue_depth=1,
                missing_phenotype=missing_phenotype)


def _expected(tmp_path, path, y, covariates, calls, options, missing_pheno):
    """The full scan's per-variant minimum p, recomputed exactly on the host."""
    run_linear_gwas(path, y, covariates, output_dir=tmp_path/'full', **options)
    beta, t, _logp, _ = open_binary_sumstats(tmp_path/'full/sumstats')
    beta, t = np.asarray(beta), np.asarray(t)
    n = y.shape[0]
    # Two covariates, the intercept and the genotype.
    if options.get('missing_phenotype') == 'exact':
        # Each pair's own samples (call and phenotype observed) less the two
        # covariates, the intercept and the genotype: complete-case OLS.
        df = ((calls != 3).astype(np.int64) @ np.isfinite(y).astype(np.int64) - 4).astype(np.float64)
    else:
        # The release: variant df times trait_df / df.
        variant_df = np.count_nonzero(calls != 3, axis=1)[:, None].astype(np.float64) - 4
        df = (variant_df * (np.isfinite(y).sum(0)[None, :] - 4) / float(n - 4) if missing_pheno
              else np.broadcast_to(variant_df, t.shape))
    exact = upper_tail_log10_from_t(t.astype(np.float64), df)
    exact = np.where(np.isfinite(exact), exact, -np.inf)
    winner = exact.argmax(axis=1)
    rows = np.arange(t.shape[0])
    return dict(variant_index=rows, trait_index=winner, beta=beta[rows, winner], t_stat=t[rows, winner],
                df=df[rows, winner], neg_log10_p=exact[rows, winner]), exact


def _read(directory):
    manifest, parts = open_indexed_sumstats(directory)
    parts = list(parts)
    values = {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}
    order = np.argsort(values['variant_index'], kind='stable')
    return manifest, {key: value[order] for key, value in values.items()}


def _assert_min_p(got, want):
    np.testing.assert_array_equal(got['variant_index'], want['variant_index'])
    np.testing.assert_array_equal(got['trait_index'], want['trait_index'])
    np.testing.assert_allclose(got['beta'], want['beta'], rtol=3e-5, atol=3e-6)
    np.testing.assert_allclose(got['t_stat'], want['t_stat'], rtol=3e-5, atol=3e-6)
    np.testing.assert_allclose(got['df'], want['df'], rtol=1e-6)
    np.testing.assert_allclose(got['neg_log10_p'], want['neg_log10_p'], rtol=3e-5, atol=1e-6)


@pytest.mark.parametrize('device', ['cpu', 'cuda:1'])
@pytest.mark.parametrize('block', [None, 4])
@pytest.mark.parametrize('missing_pheno,convention', [(False, 'impute'), (True, 'impute'), (True, 'exact')])
def test_min_p_is_the_full_scans_smallest_p_per_variant(tmp_path, monkeypatch, device, block, missing_pheno,
                                                         convention):
    if device != 'cpu' and torch.cuda.device_count() < 2:
        pytest.skip('second CUDA device required')
    path, y, covariates, calls = _inputs(tmp_path, missing_calls=device != 'cpu', missing_pheno=missing_pheno)
    for key, value in [('TORCHGWAS_PGEN_BACKEND', 'native'), ('TORCHGWAS_PGEN_PACKED', '0'),
                       ('TORCHGWAS_NATIVE_STATS', '0')]:
        monkeypatch.setenv(key, value)
    options = _options(device, convention)
    want, exact = _expected(tmp_path, path, y, covariates, calls, options, missing_pheno)
    if missing_pheno:
        # Some winners are partly missing traits, so their df is per pair.
        # (Where the largest |t| is not the smallest p is pinned below, in
        # test_reduce_ranks_by_the_exact_tail_when_pair_df_differ: with the
        # mean-imputed t, a partly missing trait rarely wins by a hair.)
        assert np.any(np.isin(want['trait_index'], [0, 2, 7, 9]))
    run_linear_gwas(path, y, covariates, output_dir=tmp_path/'minp', reduce='min-p', trait_block=block, **options)
    manifest, got = _read(tmp_path/'minp/sumstats')
    assert manifest['kind'] == 'linear' and manifest['reduction'] == 'min-p'
    assert manifest['df'] == dict(layout='per_part', axis='pair', field='df')
    _assert_min_p(got, want)
    run = json.loads((tmp_path/'minp/run.json').read_text())
    assert run['reduce'] == 'min-p'


def test_min_p_on_the_packed_native_statistics(tmp_path):
    if torch.cuda.device_count() < 2:
        pytest.skip('second CUDA device required')
    path, y, covariates, calls = _inputs(tmp_path, missing_calls=True, missing_pheno=True, m=41)
    options = _options('cuda:1')
    want, _ = _expected(tmp_path, path, y, covariates, calls, options, True)
    run_linear_gwas(path, y, covariates, output_dir=tmp_path/'minp', reduce='min_p', **options)
    _assert_min_p(_read(tmp_path/'minp/sumstats')[1], want)


@pytest.mark.parametrize('missing_pheno', [False, True])
def test_min_p_variant_shards_match_one_device(tmp_path, missing_pheno):
    if torch.cuda.device_count() < 3:
        pytest.skip('CUDA devices 1 and 2 required')
    path, y, covariates, _calls = _inputs(tmp_path, missing_calls=True, missing_pheno=missing_pheno, m=61)
    options = _options('cuda:1')
    run_linear_gwas(path, y, covariates, output_dir=tmp_path/'one', reduce='min-p', **options)
    options.pop('device')
    run_linear_gwas(path, y, covariates, output_dir=tmp_path/'two', reduce='min-p',
                    variant_devices=['cuda:1', 'cuda:2'], **options)
    one, two = _read(tmp_path/'one/sumstats')[1], _read(tmp_path/'two/sumstats')[1]
    assert one.keys() == two.keys()
    for key in one:
        np.testing.assert_array_equal(two[key], one[key])
    run = json.loads((tmp_path/'two/run.json').read_text())
    assert run['variant_devices'] == ['cuda:1', 'cuda:2']


def test_reduce_ranks_by_the_exact_tail_when_pair_df_differ():
    # Trait 0 has the larger |t| but a quarter of trait 1's df: its p is larger.
    t = torch.tensor([[4.0, 3.9]])
    beta = torch.tensor([[1.0, 2.0]])
    status = torch.zeros(1, dtype=torch.uint8)
    variant_df = torch.tensor([40.0])
    # Complete-case pair df: trait 0 observed on a quarter of trait 1's samples.
    pair_df = torch.tensor([[10.0, 40.0]])
    reduction = MinPReduction()
    kept = reduction.reduce(beta, t, status, variant_df, 1, log10_p=(pair_df, torch.float64))
    assert int(kept[2][0, 0]) == 1
    np.testing.assert_allclose(kept[6].numpy(), [[40.0]])
    np.testing.assert_allclose(kept[5].numpy(), upper_tail_log10_from_t(np.array([[3.9]], np.float32).astype(np.float64), 40.0),
                               rtol=1e-10)
    # A complete panel ranks by |t| and prices the winner alone.
    kept = reduction.reduce(beta, t, status, variant_df, 1, log10_p=(None, torch.float32))
    assert int(kept[2][0, 0]) == 0 and kept[5].dtype == torch.float32
    np.testing.assert_allclose(kept[6].numpy(), [[40.0]])
    # Without log10_p it is VariantReduction('min-p'), the internal |t| ranking.
    assert len(reduction.reduce(beta, t, status, variant_df, 1)) == 5


def test_merge_ranks_min_p_blocks_by_their_log10_p():
    reduction = MinPReduction()
    status, variant_df = torch.zeros(1, dtype=torch.uint8), torch.tensor([40.0])

    def block(t, df):
        logp = torch.as_tensor(upper_tail_log10_from_t(np.array([[t]]), df))
        return (torch.tensor([[1.0]]), torch.tensor([[t]]), None, torch.tensor([[0]], dtype=torch.int32),
                status, variant_df, logp, torch.tensor([[df]]))

    running = reduction.merge(None, block(4.0, 10.0), 0)
    running = reduction.merge(running, block(3.9, 40.0), 1)
    assert int(running[3][0, 0]) == 1 and float(running[7][0, 0]) == 40.0
    with pytest.raises(ValueError, match='all carry'):
        reduction.merge(running, block(3.0, 40.0)[:6], 2)


def test_min_p_argument_guards(tmp_path):
    genotype = np.random.default_rng(0).integers(0, 3, size=(40, 12)).astype(np.float32)
    phenotype = np.random.default_rng(1).normal(size=(40, 3))
    with pytest.raises(ValueError, match='requires output_dir'):
        run_linear_gwas(genotype, phenotype, reduce='min-p')
    with pytest.raises(ValueError, match='reduce_top_k'):
        run_linear_gwas(genotype, phenotype, reduce='min-p', reduce_top_k=2, output_dir=tmp_path/'a')
    with pytest.raises(ValueError, match='significance_threshold'):
        run_linear_gwas(genotype, phenotype, reduce='min-p', significance_threshold=1e-3, output_dir=tmp_path/'b')
    with pytest.raises(ValueError, match='p_value_threshold'):
        run_linear_gwas(genotype, phenotype, reduce='min-p', p_value_threshold=1e-3, output_dir=tmp_path/'c')


@pytest.mark.skipif(not torch.cuda.is_available(), reason='CUDA required')
def test_compiled_tail_takes_a_column_or_a_row():
    from torchgwas import tails
    device = torch.device('cuda', torch.cuda.device_count() - 1)
    tails.prepare_device_tail(device)
    # A column of winners (odd: padded; even; too small for the build), and a
    # scan's one-row last chunk, whose per-variant df broadcasts along it.
    for shape, df_shape in (((4097, 1), (4097, 1)), ((4096, 1), (4096, 1)), ((3, 1), (3, 1)),
                            ((1, 6), (1, 1)), ((1, 7), (1, 1)), ((1, 6), (1, 6))):
        cells = shape[0] * shape[1]
        t = torch.linspace(-40, 40, cells, dtype=torch.float32, device=device).reshape(shape)
        df = torch.linspace(5, 5000, df_shape[0] * df_shape[1], dtype=torch.float32, device=device).reshape(df_shape)
        got = tails.neg_log10_p_device(t, df, out=torch.empty(shape, dtype=torch.float64, device=device))
        want = upper_tail_log10_from_t(t.cpu().double().numpy(), np.broadcast_to(df.cpu().double().numpy(), shape))
        np.testing.assert_allclose(got.cpu().numpy(), want, rtol=1e-9, atol=1e-300, err_msg=str(shape))


def test_planner_shards_min_p_like_jagwas():
    from torchgwas.empirical_autotune import plan_layout
    seven = [f'cuda:{i}' for i in range(1, 8)]
    base = dict(n_samples=35_365, covariate_rank=27, n_variants=1_048_576, capacity=2048, depth=4,
                transfer_bytes_per_variant=35_365., device_free_bytes=80 * 2**30, host_free_bytes=500 * 2**30)
    wide = plan_layout(mode='min-p', n_traits=8192, devices=seven, cpus=40, reduction_width=1,
                       gpu_seconds_per_variant=8e-6, shard_setup_seconds=0.5, **base)
    assert wide['variant_devices'] and wide['trait_block'] is None
    assert wide['shard_model']['best'] == len(wide['variant_devices'])
    # A narrow panel: at most two shards, as for JAGWAS (tiny output either way).
    narrow = plan_layout(mode='min-p', n_traits=512, devices=seven, cpus=40, reduction_width=1,
                         gpu_seconds_per_variant=5e-7, shard_setup_seconds=0.1, **base)
    assert narrow['variant_devices'] == seven[:2]
    # A panel that does not fit one GPU is tiled and merged per variant.
    huge = plan_layout(mode='min-p', n_traits=2_000_000, devices=seven[:2], cpus=40, reduction_width=1,
                       **dict(base, device_free_bytes=20 * 2**30))
    assert huge['trait_block'] and huge['trait_devices'] == seven[:2]


@pytest.mark.parametrize('missing_pheno', [False, True])
def test_min_p_on_packed_bed_input_matches_pgen(tmp_path, monkeypatch, missing_pheno):
    # The packed BED path once reduced without staging the winners' -log10 P
    # and df, so its store fell back to the writer's scalar-df branch.
    if torch.cuda.device_count() < 2:
        pytest.skip('second CUDA device required')
    from test_statistics import _write_bed
    path, y, covariates, calls = _inputs(tmp_path, missing_calls=True, missing_pheno=missing_pheno, m=41)
    dosage = np.where(calls.T == 3, np.nan, calls.T).astype(float)
    bed = _write_bed(tmp_path/'input', dosage)
    for key, value in [('TORCHGWAS_PGEN_BACKEND', 'native'), ('TORCHGWAS_PGEN_PACKED', '0'),
                       ('TORCHGWAS_NATIVE_STATS', '0')]:
        monkeypatch.setenv(key, value)
    options = _options('cuda:1')
    run_linear_gwas(path, y, covariates, output_dir=tmp_path/'pgen', reduce='min-p', **options)
    options['genotype_format'] = 'plink'
    run_linear_gwas(bed, y, covariates, output_dir=tmp_path/'bed', reduce='min-p', **options)
    manifest, got = _read(tmp_path/'bed/sumstats')
    assert manifest['df'] == dict(layout='per_part', axis='pair', field='df')
    _, want = _read(tmp_path/'pgen/sumstats')
    np.testing.assert_array_equal(got['variant_index'], want['variant_index'])
    np.testing.assert_array_equal(got['trait_index'], want['trait_index'])
    np.testing.assert_allclose(got['df'], want['df'], rtol=1e-6)
    np.testing.assert_allclose(got['neg_log10_p'], want['neg_log10_p'], rtol=3e-5, atol=1e-6)
