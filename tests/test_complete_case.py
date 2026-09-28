"""Missing phenotypes: each trait is tested on its own observed samples (complete-case OLS).

The reference is FP64 least squares on the trait's observed rows (intercept,
covariates, genotype), with df = those rows less rank([1, covariates] there)
less one, as plink2's --glm does per phenotype.
"""
import numpy as np
import pytest
import torch
from scipy import special

from torchgwas import complete_case
from torchgwas.kernels import linear_chunk_statistics
from torchgwas.linear import linear_scan
from torchgwas.preprocess import residualize_and_standardize


def _panel(seed=11, n=240, m=30):
    rng = np.random.default_rng(seed)
    calls = rng.integers(0, 3, size=(n, m)).astype(np.float64)
    covariates = np.column_stack([rng.normal(size=n), rng.normal(size=n), (rng.random(n) < 0.2).astype(float)])
    effects = rng.normal(scale=0.3, size=(m, 8)) * (rng.random((m, 8)) < 0.3)
    y = rng.normal(size=(n, 8)) + (calls - calls.mean(0)) @ effects + 0.5 * covariates[:, :1]
    y[rng.choice(n, 40, replace=False), 1] = np.nan
    y[:, 2] = np.where(np.isnan(y[:, 1]), np.nan, y[:, 2])            # the same pattern as trait 1
    y[covariates[:, 0] > np.quantile(covariates[:, 0], 0.6), 3] = np.nan  # missing by a covariate
    y[7, 4] = np.nan                                                   # one missing value
    y[covariates[:, 2] == 1, 5] = np.nan                               # the indicator is constant on S
    y[rng.choice(n, 150, replace=False), 6] = np.nan                   # mostly missing
    return calls, y, covariates


def _reference(calls, y, covariates):
    n, m = calls.shape
    t = np.full((m, y.shape[1]), np.nan)
    df = np.full((m, y.shape[1]), np.nan)
    for j in range(y.shape[1]):
        keep = ~np.isnan(y[:, j])
        base = np.column_stack([np.ones(keep.sum()), covariates[keep]])
        rank = np.linalg.matrix_rank(base)
        for v in range(m):
            x = np.column_stack([base, calls[keep, v]])
            coef = np.linalg.pinv(x) @ y[keep, j]
            resid = y[keep, j] - x @ coef
            df[v, j] = keep.sum() - rank - 1
            cov = np.linalg.pinv(x.T @ x) * (resid @ resid / df[v, j])
            t[v, j] = coef[-1] / np.sqrt(cov[-1, -1])
    return t, df


@pytest.mark.parametrize('dtype', [np.float64, np.float32])
def test_the_scan_is_complete_case_ols_for_every_pattern(dtype):
    calls, y, covariates = _panel()
    want_t, want_df = _reference(calls, y, covariates)
    # float32 inputs give a float32 covariate basis: a covariate constant on a
    # subset then leaves a noise direction both projections must drop.
    tolerance = 1e-8 if dtype == np.float64 else 2e-5
    _, t, logp, _ = linear_scan(calls, y.astype(dtype), covariates.astype(dtype), device='cpu', missing_phenotype='exact',
                                compute_dtype='float64', return_log10_p=True)
    np.testing.assert_allclose(t, want_t, rtol=tolerance, atol=tolerance)
    want_logp = -np.log10(2 * special.stdtr(want_df, -np.abs(want_t)))
    np.testing.assert_allclose(logp, want_logp, rtol=tolerance, atol=tolerance)
    if dtype == np.float32:
        return
    _, _, p, _ = linear_scan(calls, y, covariates, device='cpu', missing_phenotype='exact', compute_dtype='float64')
    np.testing.assert_allclose(p, 2 * special.stdtr(want_df, -np.abs(want_t)), rtol=1e-7, atol=1e-300)


def test_the_plan_groups_patterns_and_counts_the_subset_rank():
    calls, y, covariates = _panel()
    _, q, counts, plan = residualize_and_standardize(y, covariates, return_observed_counts=True,
                                                     missing='complete_case', return_plan=True)
    np.testing.assert_array_equal(plan.traits, [1, 2, 3, 4, 5, 6])
    assert plan.pattern_of_trait[0] == plan.pattern_of_trait[1] and len(plan.rows) == 5
    # Trait 5 drops every sample of the indicator: one rank fewer on its subset.
    assert plan.subset_rank[plan.pattern_of_trait[4]] == q.shape[1]
    assert plan.subset_rank[plan.pattern_of_trait[0]] == q.shape[1] + 1
    np.testing.assert_array_equal(plan.observed_counts[plan.pattern_of_trait], counts[plan.traits])
    # The imputed panel is unchanged; its plan is the release's rescaling.
    imputed = residualize_and_standardize(y, covariates, return_plan=True)
    assert isinstance(imputed[-1], complete_case.ImputedPlan)
    np.testing.assert_array_equal(imputed[-1].traits, plan.traits)


def test_small_blocks_and_split_patterns_give_the_same_answer(monkeypatch):
    calls, y, covariates = _panel()
    _, want, _, _ = linear_scan(calls, y, covariates, device='cpu', missing_phenotype='exact', compute_dtype='float64')
    monkeypatch.setattr(complete_case, 'BLOCK_CELLS', 16)
    monkeypatch.setattr(complete_case, 'PATTERN_BLOCK', 2)
    _, got, _, _ = linear_scan(calls, y, covariates, device='cpu', missing_phenotype='exact', compute_dtype='float64')
    np.testing.assert_allclose(got, want, rtol=1e-10, atol=1e-12)


def test_missing_calls_leave_the_pair_with_its_own_count():
    calls, y, covariates = _panel()
    calls[np.random.default_rng(3).random(calls.shape) < 0.05] = np.nan
    pheno, q, _, plan = residualize_and_standardize(y, covariates, return_observed_counts=True,
                                                    missing='complete_case', return_plan=True)
    rank = q.shape[1]
    geno = torch.as_tensor(calls)
    _, t, variant_df, pair_df = linear_chunk_statistics(
        geno, torch.as_tensor(pheno), torch.as_tensor(q), calls.shape[0] - rank - 2,
        covariate_rank=rank, complete_case=plan)
    both = (~np.isnan(calls)).astype(int).T @ (~np.isnan(y)).astype(int)
    subset_rank = np.full(y.shape[1], rank + 1)
    subset_rank[plan.traits] = plan.subset_rank[plan.pattern_of_trait]
    np.testing.assert_array_equal(pair_df.numpy(), both - subset_rank[None, :] - 1)
    np.testing.assert_array_equal(variant_df.numpy(), (~np.isnan(calls)).sum(0) - rank - 2)
    assert np.isfinite(t.numpy()).all()


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason='second CUDA device required')
@pytest.mark.parametrize('backend', ['pgen-torch', 'pgen-native-dosage', 'pgen-native-packed', 'bed'])
def test_every_gpu_backend_writes_complete_case_statistics(tmp_path, monkeypatch, backend):
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats import open_binary_sumstats
    from test_pgen_native_reader import write_pgen
    from test_statistics import _write_bed
    calls, y, covariates = _panel(n=203, m=37)
    want_t, want_df = _reference(calls, y, covariates)
    n, m = calls.shape
    if backend == 'bed':
        path, fmt = _write_bed(tmp_path/'input', calls), 'plink'
    else:
        path, fmt = tmp_path/'input.pgen', 'pgen'
        write_pgen(path, calls.T.astype(np.uint8))
        path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'
                                             + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
        path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
        monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
        monkeypatch.setenv('TORCHGWAS_NATIVE_STATS', '0' if backend == 'pgen-torch' else '1')
        monkeypatch.setenv('TORCHGWAS_PGEN_PACKED', '1' if backend == 'pgen-native-packed' else '0')
    run_linear_gwas(path, y.astype(np.float32), covariates.astype(np.float32), output_dir=tmp_path/'out', missing_phenotype='exact',
                    genotype_format=fmt, pgen_mode='hardcall', device='cuda:1', compute_dtype='float32',
                    chunk_size=8, reader_workers=2, prefetch_chunks=2, sumstats_queue_depth=1)
    _, t, logp, _ = open_binary_sumstats(tmp_path/'out/sumstats')
    np.testing.assert_allclose(np.asarray(t), want_t, rtol=2e-4, atol=2e-4)
    want_logp = -np.log10(2 * special.stdtr(want_df, -np.abs(want_t)))
    np.testing.assert_allclose(np.asarray(logp), want_logp, rtol=2e-4, atol=2e-4)


@pytest.mark.skipif(torch.cuda.device_count() < 2, reason='second CUDA device required')
def test_min_p_ranks_complete_case_p_values(tmp_path, monkeypatch):
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    from test_pgen_native_reader import write_pgen
    calls, y, covariates = _panel(n=203, m=37)
    want_t, want_df = _reference(calls, y, covariates)
    n, m = calls.shape
    path = tmp_path/'input.pgen'
    write_pgen(path, calls.T.astype(np.uint8))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    run_linear_gwas(path, y.astype(np.float32), covariates.astype(np.float32), output_dir=tmp_path/'out', missing_phenotype='exact',
                    reduce='min-p', genotype_format='pgen', pgen_mode='hardcall', device='cuda:1',
                    compute_dtype='float32', chunk_size=8, reader_workers=2, prefetch_chunks=2)
    _, parts = open_indexed_sumstats(tmp_path/'out/sumstats')
    parts = list(parts)
    got = {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}
    order = np.argsort(got['variant_index'])
    want_logp = -np.log10(2 * special.stdtr(want_df, -np.abs(want_t)))
    np.testing.assert_array_equal(got['trait_index'][order], want_logp.argmax(axis=1))
    np.testing.assert_allclose(got['neg_log10_p'][order], want_logp.max(axis=1), rtol=2e-4)
    np.testing.assert_array_equal(got['df'][order], want_df[np.arange(m), want_logp.argmax(axis=1)])


@pytest.mark.skipif(torch.cuda.device_count() < 3, reason='CUDA devices 1 and 2 required')
def test_variant_shards_store_device_log10_p_for_missing_phenotypes(tmp_path, monkeypatch):
    # The shard writer reads -log10 P's presence from the chunk's length; a
    # chunk without its df once made it recompute every tail on the host
    # (19 s against 1.9 s for 200k variants at K = 512).
    from torchgwas import tails
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats import open_binary_sumstats
    from test_pgen_native_reader import write_pgen
    calls, y, covariates = _panel(n=203, m=37)
    want_t, want_df = _reference(calls, y, covariates)
    n, m = calls.shape
    path = tmp_path/'input.pgen'
    write_pgen(path, calls.T.astype(np.uint8))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')

    def refuse(*args, **kwargs):
        raise AssertionError('the store recomputed -log10 P on the host')
    monkeypatch.setattr(tails, 'upper_tail_log10_from_t', refuse)
    run_linear_gwas(path, y.astype(np.float32), covariates.astype(np.float32), output_dir=tmp_path/'out', missing_phenotype='exact',
                    genotype_format='pgen', pgen_mode='hardcall', variant_devices=['cuda:1', 'cuda:2'],
                    compute_dtype='float32', chunk_size=8, reader_workers=4, prefetch_chunks=2)
    _, t, logp, _ = open_binary_sumstats(tmp_path/'out/sumstats')
    np.testing.assert_allclose(np.asarray(t), want_t, rtol=2e-4, atol=2e-4)
    want_logp = -np.log10(2 * special.stdtr(want_df, -np.abs(want_t)))
    np.testing.assert_allclose(np.asarray(logp), want_logp, rtol=2e-4, atol=2e-4)
