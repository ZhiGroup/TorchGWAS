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


def _pattern_reference(observed, q):
    """The plan pattern by pattern: each trait's missing rows, each pattern's own eigendecomposition."""
    basis = complete_case.complete_case_basis(observed.shape[0], q)
    gram = basis.T @ basis
    traits = np.flatnonzero(~observed.all(axis=0))
    index, rows, of_trait, inverse, rank = {}, [], [], [], []
    for trait in traits:
        missing = np.flatnonzero(~observed[:, trait])
        if missing.tobytes() not in index:
            index[missing.tobytes()] = len(rows)
            rows.append(missing)
            values, vectors = np.linalg.eigh(gram - basis[missing].T @ basis[missing])
            keep = values > 1e-9 * max(values.max(), 1e-300)
            inverse.append((vectors[:, keep] / values[keep]) @ vectors[:, keep].T)
            rank.append(int(keep.sum()))
        of_trait.append(index[missing.tobytes()])
    return traits, rows, np.array(of_trait), np.array(inverse), np.array(rank)


@pytest.mark.parametrize('cells', [complete_case.DOWNDATE_CELLS, 7])
def test_the_batched_plan_matches_a_pattern_by_pattern_reference(monkeypatch, cells):
    # 7 cells a batch: most patterns share batches, the 150-row one has its own.
    monkeypatch.setattr(complete_case, 'DOWNDATE_CELLS', cells)
    calls, y, covariates = _panel()
    y = np.column_stack([y, y[:, [6, 3]]])                 # more traits sharing patterns
    from torchgwas.preprocess import _covariate_basis
    q = _covariate_basis(covariates)
    observed = ~np.isnan(y)
    plan = complete_case.CompleteCasePlan(observed, q)
    traits, rows, of_trait, inverse, rank = _pattern_reference(observed, q)
    np.testing.assert_array_equal(plan.traits, traits)
    np.testing.assert_array_equal(plan.pattern_of_trait, of_trait)
    assert len(plan.rows) == len(rows) and all(np.array_equal(a, b) for a, b in zip(plan.rows, rows))
    np.testing.assert_allclose(plan.inverse, inverse, rtol=1e-9, atol=1e-12)
    np.testing.assert_array_equal(plan.subset_rank, rank)
    np.testing.assert_array_equal(plan.observed_counts, [observed.shape[0] - len(r) for r in rows])
    # The panel read a few columns at a time gives the same plan; a complete one none.
    panel = complete_case.CompleteCasePlan.from_panel(y, q, block=3)
    np.testing.assert_array_equal(panel.pattern_of_trait, plan.pattern_of_trait)
    np.testing.assert_array_equal(panel.inverse, plan.inverse)
    assert complete_case.CompleteCasePlan.from_panel(np.nan_to_num(y), q) is None
    # A block of traits residualized at once is each trait residualized alone.
    block = plan.residualize_block(y[:, plan.traits], plan.traits)
    for column, trait in enumerate(plan.traits):
        np.testing.assert_allclose(block[:, column], plan.residualize(y[:, trait], trait), rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize('trait_block', [1, 3, 4096])
def test_residualizing_in_blocks_is_trait_by_trait(trait_block):
    calls, y, covariates = _panel()
    got, q, plan = residualize_and_standardize(y, covariates, missing='complete_case', return_plan=True,
                                               trait_block=trait_block)
    for trait in plan.traits:
        np.testing.assert_allclose(got[:, trait], plan.residualize(y[:, trait], trait), rtol=1e-10, atol=1e-12)
    complete = np.setdiff1d(np.arange(y.shape[1]), plan.traits)
    want = residualize_and_standardize(y[:, complete], covariates)[0]
    np.testing.assert_allclose(got[:, complete], want, rtol=1e-10, atol=1e-12)


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


@pytest.mark.skipif(not torch.cuda.is_available(), reason='CUDA required')
@pytest.mark.parametrize('missing', ['exact', 'impute'])
@pytest.mark.parametrize('stats', ['native', 'torch'])
def test_significant_pairs_select_on_the_device_at_each_pairs_df(tmp_path, monkeypatch, stats, missing):
    # Missing-phenotype pairs leave the device already selected, each at its
    # own df (whole numbers under 'exact', fractional under 'impute'): the
    # same pairs, t and df as the host selector, and under 'exact' the reference.
    from torchgwas.api import run_linear_gwas
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    from test_pgen_native_reader import write_pgen
    calls, y, covariates = _panel(n=203, m=37)
    want_t, want_df = _reference(calls, y, covariates)
    n, m = calls.shape
    path = tmp_path / 'input.pgen'
    written = calls.T.astype(np.uint8)
    if missing == 'impute':
        # Missing calls give each variant its own df, so 'impute' pair df turn fractional.
        written[np.random.default_rng(8).random(written.shape) < 0.03] = 3
    write_pgen(path, written)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND', 'native')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS', '1' if stats == 'native' else '0')
    threshold = 0.2
    found, calls = {}, []
    selector = complete_case.device_significant_pairs_by_pair_df

    def spy(*args, **kwargs):
        calls.append(backend)
        yield from selector(*args, **kwargs)
    monkeypatch.setattr(complete_case, 'device_significant_pairs_by_pair_df', spy)
    for backend in ('device', 'host'):
        monkeypatch.setenv('TORCHGWAS_SIGNIFICANCE_BACKEND', backend)
        result = run_linear_gwas(path, y.astype(np.float32), covariates.astype(np.float32),
                                 output_dir=tmp_path / backend, missing_phenotype=missing, reduce='significant',
                                 significance_threshold=threshold, genotype_format='pgen', pgen_mode='hardcall',
                                 device='cuda:0', compute_dtype='float32', chunk_size=8, reader_workers=2,
                                 prefetch_chunks=2, sumstats_queue_depth=1)
        _, parts = open_indexed_sumstats(tmp_path / backend / 'sumstats')
        parts = list(parts)
        found[backend] = {key: np.concatenate([part[key] for part in parts]) for key in parts[0]}
    assert calls and set(calls) == {'device'}       # the device selected; the host run did not
    device, host = found['device'], found['host']
    for key in ('variant_index', 'trait_index', 'df'):
        np.testing.assert_array_equal(device[key], host[key], err_msg=key)
    np.testing.assert_allclose(device['t_stat'], host['t_stat'], rtol=1e-6, atol=1e-6)
    if missing == 'impute':
        assert (device['df'] != np.floor(device['df'])).any()   # fractional df were selected
        return
    # The reference's passing pairs, away from the threshold's float32 edge.
    p = 2 * special.stdtr(want_df, -np.abs(want_t))
    v, k = device['variant_index'], device['trait_index']
    np.testing.assert_array_equal(device['df'], want_df[v, k])
    np.testing.assert_allclose(device['t_stat'], want_t[v, k], rtol=2e-4, atol=2e-4)
    clear = (p < threshold * 0.99) & np.isfinite(p)
    assert clear.sum() > 10 and set(zip(*np.nonzero(clear))) <= set(zip(v.tolist(), k.tolist()))
    assert not (p[v, k] > threshold * 1.01).any()


@pytest.mark.parametrize('device', ['cpu'] + (['cuda:0'] if torch.cuda.is_available() else []))
def test_fractional_pair_df_select_exactly_at_the_threshold(device):
    # |t| just either side of each pair's own critical value, at fractional
    # and whole-number df: the selection is the FP64 two-sided p <= threshold.
    from scipy import special as sp
    from torchgwas.reduce import SignificantPairs, device_significance_critical
    threshold, n = 1e-6, 400
    rng = np.random.default_rng(3)
    df = np.round(rng.uniform(50, 390, size=(64, 32)), 0)
    df[:, ::2] = rng.uniform(50, 390, size=(64, 16))              # half the columns fractional
    critical = np.abs(sp.stdtrit(df, threshold / 2))
    t = critical * (1 + rng.choice([-1e-5, 1e-5, -0.2, 0.2], size=df.shape)) * rng.choice([-1, 1], size=df.shape)
    t = t.astype(np.float32)
    want = 2 * sp.stdtr(df, -np.abs(t.astype(np.float64))) <= threshold
    if device == 'cpu':
        got_critical = SignificantPairs(threshold).critical_abs_t(df, 32)
        np.testing.assert_allclose(got_critical, critical, rtol=1e-12)
        return
    table = device_significance_critical(SignificantPairs(threshold), n, 32, device)
    tt = torch.as_tensor(t, device=device)
    chunks = list(complete_case.device_significant_pairs_by_pair_df(
        tt, tt, torch.zeros(64, dtype=torch.int8, device=device), torch.as_tensor(df, dtype=torch.float32, device=device),
        table, threshold=threshold))
    got = np.zeros(df.shape, dtype=bool)
    for chunk in chunks:
        got[chunk[2], chunk[3]] = True                            # rows are absolute (start=0)
    df32 = df.astype(np.float32).astype(np.float64)               # the df the device saw
    want = 2 * sp.stdtr(df32, -np.abs(t.astype(np.float64))) <= threshold
    np.testing.assert_array_equal(got, want)
    assert want.any() and not want.all()


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
