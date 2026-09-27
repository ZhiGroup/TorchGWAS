import math

import numpy as np
import pytest
import torch

from torchgwas.jagwas_projection import JagwasGroups, JagwasRankSelection, JagwasReduction, triangular_blocks

U32 = 2.0 ** -24


def _panel(samples, traits, seed=7):
    rng = np.random.default_rng(seed)
    base = rng.standard_normal((samples, traits))
    panel = base + 0.4 * base.mean(axis=1, keepdims=True)
    return torch.as_tensor((panel - panel.mean(0)) / panel.std(0), dtype=torch.float32)


def _standardize(panel):
    return torch.as_tensor((panel - panel.mean(0)) / panel.std(0), dtype=torch.float32)


def _collinear_panel(samples=3000, independent=30, seed=11):
    """30 independent traits, then an exact duplicate (of 5) and an exact combination (of 0 and 1)."""
    rng = np.random.default_rng(seed)
    base = rng.standard_normal((samples, independent))
    return _standardize(np.column_stack([base, base[:, 5], base[:, 0] + base[:, 1]]))


def _score(t, df):
    """The statistic's input, z = t / sqrt(1 + t^2 / df), with the reduction's FP32 operations."""
    return t * torch.rsqrt(t.square() / df.unsqueeze(1) + 1.0)


def _df(rows, value, device='cpu'):
    return torch.full((rows,), float(value), dtype=torch.float32, device=device)


def _reduce(reduction, t, df):
    beta = torch.zeros_like(t)
    status = torch.zeros(t.shape[0], dtype=torch.uint8, device=t.device)
    return reduction.reduce(beta, t, status, _df(t.shape[0], df, t.device), 1)[1][:, 0]


def _exact_quadratic_form(panel, z):
    """z' R^-1 z in FP64 with R the exact Gram of the FP32 panel (independent of the reduction)."""
    p = panel.double()
    r = p.T @ p / p.shape[0]
    return (z.double() @ torch.linalg.inv(r) * z.double()).sum(1)


def test_blocks_cover_every_row_once():
    for traits in (1, 7, 511, 512, 1500, 8192, 40_000):
        blocks = triangular_blocks(traits)
        assert blocks[0][0] == 0 and blocks[-1][1] == traits
        assert all(a[1] == b[0] for a, b in zip(blocks, blocks[1:]))
        assert len(blocks) <= 16


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
@pytest.mark.parametrize('traits', [3, 40, 1100])
def test_projection_matches_the_dense_quadratic_form(traits, device):
    if device == 'cuda' and not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    panel = _panel(3000, traits)
    reduction = JagwasReduction(rcond=0).prepare(panel, device=device)
    rng = np.random.default_rng(traits)
    t = torch.as_tensor(rng.standard_normal((257, traits)) * 3, dtype=torch.float32)
    t[5, traits // 2] = float('nan')
    t[9, 0] = float('inf')
    t = t.to(device)
    beta = torch.zeros_like(t)
    status = torch.zeros(257, dtype=torch.uint8, device=device)
    status[11] = 2
    df = _df(257, 2997, device)
    got = reduction.reduce(beta, t, status, df, 1)
    # Dense ||L^-1 z||^2 with the same factor: the blocks change only the summation order.
    z = _score(torch.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), df).double()
    expected = (z @ reduction._inverse_cholesky.T).square().sum(1)
    finite = torch.isfinite(got[1][:, 0])
    assert finite.tolist() == [i not in (5, 9, 11) for i in range(257)]
    assert got[1].device == t.device and got[1].dtype == torch.float32
    torch.testing.assert_close(got[1][finite, 0].double(), expected[finite], rtol=1e-6, atol=0)
    assert torch.equal(got[2], torch.zeros_like(got[2])) and torch.equal(got[3], status)
    assert torch.isnan(got[0]).all()
    assert reduction.degrees_of_freedom == traits


def test_triangular_blocks_are_views_of_one_factor():
    reduction = JagwasReduction(rcond=0).prepare(_panel(2000, 1100))
    base = reduction._inverse_cholesky.untyped_storage().data_ptr()
    assert all(rows.untyped_storage().data_ptr() == base for *_, rows in reduction._blocks)


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
def test_full_rank_panel_keeps_the_unpivoted_factor(device):
    if device == 'cuda' and not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    panel = _panel(9000, 40)  # three FP64 Gram sample blocks
    reduction = JagwasReduction(rcond=0).prepare(panel, device=device)
    exact = panel.double().T @ panel.double() / 9000.
    report = reduction.rank_report()
    assert report['method'] == 'cholesky' and report['rank'] == 40 and report['dropped'] == []
    assert report['target'] == 0.01 and report['precision'] == pytest.approx(U32 * math.sqrt(9000), rel=1e-12)
    # 2 eps_z sqrt(tr R^-1): the first-order null rms rounding error of T.
    assert report['rounding_error'] == pytest.approx(
        2 * report['precision'] * math.sqrt(float(torch.linalg.inv(exact).trace())), rel=1e-9)
    assert reduction._kept is None and reduction.degrees_of_freedom == 40
    # R is the exact (FP64) Gram of the FP32 panel.
    expected = torch.linalg.solve_triangular(torch.linalg.cholesky(exact), torch.eye(40, dtype=torch.float64), upper=False)
    torch.testing.assert_close(reduction._inverse_cholesky.cpu(), expected, rtol=1e-12, atol=1e-12)


def test_fast_path_and_pivoted_selection_agree():
    # tr R_S^-1 does not depend on the order: a set that meets the target in
    # input order keeps every trait in the greedy order too, at the same error.
    panel = _panel(4000, 60)
    reduction = JagwasReduction(rcond=0).prepare(panel)
    correlation = (panel.double().T @ panel.double()) / 4000.
    precision, limit = reduction._selection.trace_limit(4000, 60, torch.float32)
    pivoted = reduction._pivoted_selection(correlation, limit, precision)
    assert pivoted['rank'] == 60 and pivoted['dropped'] == []
    assert pivoted['rounding_error'] == pytest.approx(reduction.rank_report()['rounding_error'], rel=1e-9)
    # A limit just below the full trace keeps one trait fewer.
    trace = float(torch.linalg.inv(correlation).trace())
    assert reduction._pivoted_selection(correlation, trace * (1 - 1e-9), precision)['rank'] == 59


@pytest.mark.parametrize('device', ['cpu', 'cuda'])
def test_exactly_collinear_traits_are_dropped(device):
    if device == 'cuda' and not torch.cuda.is_available():
        pytest.skip('CUDA device required')
    panel = _collinear_panel()
    with pytest.warns(UserWarning, match='2 of 32 traits are collinear'):
        reduction = JagwasReduction(rcond=0).prepare(panel, device=device)
    report = reduction.rank_report([f'y{i}' for i in range(32)])
    assert report['method'] == 'pivoted_cholesky' and report['rank'] == 30 == reduction.degrees_of_freedom
    assert report['rounding_error'] <= report['target']
    dropped = [item['index'] for item in report['dropped']]
    assert len(dropped) == 2 and len({5, 30} & set(dropped)) == 1 and len({0, 1, 31} & set(dropped)) == 1
    assert all(abs(item['residual_variance']) < 1e-10 for item in report['dropped'])
    assert [item['trait'] for item in report['dropped']] == [f'y{i}' for i in dropped]
    kept = [i for i in range(32) if i not in dropped]
    assert reduction._kept.tolist() == kept
    # T is the quadratic form over the kept traits; a dropped trait does not enter it.
    t = torch.as_tensor(np.random.default_rng(3).standard_normal((300, 32)) * 2, dtype=torch.float32)
    t[7, dropped[0]] = float('nan')
    t[8, kept[3]] = float('nan')
    got = _reduce(reduction, t.to(device), 2997).cpu()
    expected = _exact_quadratic_form(panel[:, kept], _score(torch.nan_to_num(t[:, kept]), _df(300, 2997)))
    assert torch.isnan(got[8]) and torch.isfinite(got[7])
    finite = torch.isfinite(got)
    torch.testing.assert_close(got[finite].double(), expected[finite], rtol=2.5e-7, atol=0)


def test_fp64_statistics_still_drop_exact_collinearity():
    # With FP64 statistics the rounding target alone allows tr R^-1 ~ 1e25; an
    # exact duplicate (pivot ~ eps64) must still go: R's FP64 numerical rank.
    rng = np.random.default_rng(12)
    base = rng.standard_normal((3000, 30))
    panel = np.column_stack([base, base[:, 5], base[:, 0] + base[:, 1]])
    panel = torch.as_tensor((panel - panel.mean(0)) / panel.std(0), dtype=torch.float64)
    with pytest.warns(UserWarning, match='2 of 32 traits are collinear'):
        reduction = JagwasReduction(rcond=0).prepare(panel)
    assert reduction.degrees_of_freedom == 30
    precision, limit = reduction._selection.trace_limit(3000, 32, torch.float64)
    assert limit == pytest.approx(1 / (32 * np.finfo(np.float64).eps))
    # Full rank in FP64 keeps everything with the unpivoted factor.
    assert JagwasReduction(rcond=0).prepare(panel[:, :30]).rank_report()['method'] == 'cholesky'


def test_residual_variance_matches_explicit_schur_complement():
    rng = np.random.default_rng(5)
    base = rng.standard_normal((4000, 12))
    near = base[:, 3] + 1e-4 * rng.standard_normal(4000)  # 1 - R^2 ~ 1e-8: VIF ~ 1e8
    panel = _standardize(np.column_stack([base, near, base[:, :4] @ [1., -2., .5, 3.]]))
    reduction = JagwasReduction(rcond=0).prepare(panel)
    report = reduction.rank_report()
    assert report['rank'] == 12 and len(report['dropped']) == 2
    # R is the exact (FP64) Gram of the FP32 panel; LAPACK and the factor read its lower triangle.
    lower = ((panel.double().T @ panel.double()) / 4000.).numpy()
    correlation = np.tril(lower) + np.tril(lower, -1).T
    kept = reduction._kept.numpy()
    near = next(item for item in report['dropped'] if item['index'] == 12)
    assert 5e-9 < near['residual_variance'] < 2e-8
    for item in report['dropped']:
        j = item['index']
        residual = correlation[j, j] - correlation[j, kept] @ np.linalg.solve(correlation[np.ix_(kept, kept)],
                                                                               correlation[kept, j])
        assert item['residual_variance'] == pytest.approx(residual / correlation[j, j], abs=1e-12)
        if item['residual_variance'] > 0:
            assert item['vif'] == pytest.approx(1 / item['residual_variance'], rel=1e-12)
        else:
            assert item['vif'] is None


def test_target_moves_the_cutoff(monkeypatch):
    rng = np.random.default_rng(8)
    base = rng.standard_normal((20000, 8))
    noise = rng.standard_normal((20000, 2))
    # Residual variances given trait 0: ~1e-2 and ~1e-4, so tr R^-1 ~ 2e4.
    panel = _standardize(np.column_stack([base, base[:, 0] + .1 * noise[:, 0], base[:, 0] + .01 * noise[:, 1]]))
    monkeypatch.delenv('TORCHGWAS_JAGWAS_T_ROUNDING', raising=False)
    assert JagwasReduction(rcond=0).prepare(panel).degrees_of_freedom == 10  # 2 eps_z sqrt(2e4) ~ 2.4e-3 <= 0.01
    with pytest.warns(UserWarning):
        tight = JagwasReduction(JagwasRankSelection(target=1e-3), rcond=0).prepare(panel)
    assert tight.degrees_of_freedom == 9 and tight.rank_report()['rounding_error'] <= 1e-3
    monkeypatch.setenv('TORCHGWAS_JAGWAS_T_ROUNDING', '1e-4')
    with pytest.warns(UserWarning):
        assert JagwasReduction(rcond=0).prepare(panel).degrees_of_freedom == 8
    for bad in ('0', '-1', 'nan'):
        monkeypatch.setenv('TORCHGWAS_JAGWAS_T_ROUNDING', bad)
        with pytest.raises(ValueError, match='T_ROUNDING'):
            JagwasReduction(rcond=0)


def test_null_statistic_is_chi_square_on_the_kept_rank():
    panel = _collinear_panel(samples=2000)
    with pytest.warns(UserWarning):
        reduction = JagwasReduction(rcond=0).prepare(panel)
    # Null z-scores with exactly the panel correlation: z = Y'g / sqrt(N) for g ~ N(0, I).
    genotypes = torch.as_tensor(np.random.default_rng(21).standard_normal((2000, 20000)), dtype=torch.float32)
    z = (panel.T @ genotypes / np.sqrt(2000.)).T.contiguous()
    statistic = _reduce(reduction, z, 1997).double().numpy()
    from scipy import stats
    assert abs(statistic.mean() - 30) < 4 * np.sqrt(2 * 30 / 20000)
    assert stats.kstest(statistic, 'chi2', args=(30,)).pvalue > 1e-3
    assert stats.kstest(statistic, 'chi2', args=(32,)).pvalue < 1e-6


def test_score_form_is_invariant_to_redundant_traits_under_a_strong_effect():
    rng = np.random.default_rng(31)
    samples, df = 5000, 4998
    base = rng.standard_normal((samples, 6))
    independent = _standardize(base)
    # Trait 0 replaced by a combination of traits 0 and 1: the same span, a different basis.
    swapped = _standardize(np.column_stack([base[:, 0] + 2 * base[:, 1], base[:, 1:]]))
    redundant = _standardize(np.column_stack([base, base[:, 0] + 2 * base[:, 1]]))
    # Strong effects on trait 0 (z ~ 20 and up), t exactly as the scan forms it from r.
    genotype = 0.3 * independent[:, :1].numpy() + rng.standard_normal((samples, 32)) * np.linspace(1, 3, 32)
    genotype = (genotype - genotype.mean(0)) / genotype.std(0)

    def t_of(panel):
        r = panel.double().T @ torch.as_tensor(genotype) / samples
        return (np.sqrt(df) * r / torch.sqrt(1 - r * r)).T.float().contiguous()

    assert t_of(independent).abs().max() > 18
    with pytest.warns(UserWarning, match='1 of 7 traits are collinear'):
        pivoted = JagwasReduction(rcond=0).prepare(redundant)
    values = {name: _reduce(JagwasReduction(rcond=0).prepare(panel), t_of(panel), df).double()
              for name, panel in (('independent', independent), ('swapped', swapped))}
    values['redundant'] = _reduce(pivoted, t_of(redundant), df).double()
    torch.testing.assert_close(values['swapped'], values['independent'], rtol=2e-5, atol=0)
    torch.testing.assert_close(values['redundant'], values['independent'], rtol=2e-5, atol=0)
    # A quadratic form of t itself changes with the basis for the same effects.
    legacy = {name: _exact_quadratic_form(panel, t_of(panel)) for name, panel in
              (('independent', independent), ('swapped', swapped))}
    assert ((legacy['swapped'] - legacy['independent']).abs() / legacy['independent']).max() > 1e-3


def test_devices_share_one_kept_set():
    panel = _collinear_panel()
    first = JagwasReduction(rcond=0)
    with pytest.warns(UserWarning):
        first.prepare(panel)
    second = first.spawn()
    assert second is not first and second._selection is first._selection
    # A device whose own panel looks full rank still adopts the run's decision.
    second.prepare(_panel(3000, 32))
    assert second._kept.tolist() == first._kept.tolist() and second.degrees_of_freedom == 30
    unprepared = first.spawn()
    assert unprepared.degrees_of_freedom == 30
    assert JagwasReduction(rcond=0).rank_report() is None
    with pytest.raises(ValueError, match='not prepared'):
        JagwasReduction(rcond=0).degrees_of_freedom


def test_groups_match_separate_reductions():
    independent, collinear = _panel(3000, 12, seed=3), _collinear_panel()
    groups = JagwasGroups([('a', range(12)), ('b', range(12, 44))], rcond=0)
    with pytest.warns(UserWarning, match='jagwas group b: 2 of 32 traits'):
        groups.prepare(torch.cat([independent, collinear], dim=1))
    alone_a, alone_b = JagwasReduction(rcond=0), JagwasReduction(rcond=0)
    alone_a.prepare(independent)
    with pytest.warns(UserWarning, match='^jagwas: 2 of 32'):
        alone_b.prepare(collinear)
    t = torch.as_tensor(np.random.default_rng(5).standard_normal((64, 44)) * 3, dtype=torch.float32)
    t[7, 20] = float('nan')  # invalidates group b only
    beta, status = torch.zeros_like(t), torch.zeros(64, dtype=torch.uint8)
    statistic = groups.reduce(beta, t, status, _df(64, 2998), groups.resolved_width(44))[1]
    assert statistic.shape == (64, 2) and groups.degrees_of_freedom == [12, 30]
    torch.testing.assert_close(statistic[:, 0], _reduce(alone_a, t[:, :12], 2998), rtol=1e-6, atol=0)
    torch.testing.assert_close(statistic[:, 1], _reduce(alone_b, t[:, 12:], 2998), rtol=1e-6, atol=0,
                               equal_nan=True)
    assert torch.isfinite(statistic[7, 0]) and torch.isnan(statistic[7, 1])
    names = [f'a{i}' for i in range(12)] + [f'b{i}' for i in range(32)]
    report_a, report_b = groups.rank_report(names)
    assert report_a['group'] == 'a' and report_a['rank'] == 12 and report_a['dropped'] == []
    assert report_b == dict(group='b', **alone_b.rank_report(names[12:]))
    spawned = groups.spawn()
    assert all(mine._selection is theirs._selection for mine, theirs in zip(spawned.reductions, groups.reductions))
    assert spawned.degrees_of_freedom == [12, 30]


def test_default_cutoff_is_eigen_truncation(monkeypatch):
    monkeypatch.delenv('TORCHGWAS_JAGWAS_RCOND', raising=False)
    monkeypatch.delenv('TORCHGWAS_JAGWAS_MIN_RESIDUAL', raising=False)
    assert (JagwasReduction().rcond, JagwasReduction().min_residual) == (1e-3, None)
    assert JagwasReduction(rcond=0).rcond is None  # the rounding cutoff over traits
    monkeypatch.setenv('TORCHGWAS_JAGWAS_RCOND', '0')
    assert JagwasReduction().rcond is None
    monkeypatch.delenv('TORCHGWAS_JAGWAS_RCOND')
    monkeypatch.setenv('TORCHGWAS_JAGWAS_MIN_RESIDUAL', '0.01')
    assert (JagwasReduction().rcond, JagwasReduction().min_residual) == (None, 0.01)
    monkeypatch.setenv('TORCHGWAS_JAGWAS_RCOND', '1e-3')
    with pytest.raises(ValueError, match='exclusive'):
        JagwasReduction()


def test_eigen_factor_is_upper_trapezoidal_and_the_truncated_pseudo_inverse():
    panel = _collinear_panel()
    with pytest.warns(UserWarning, match='eigen-directions'):
        reduction = JagwasReduction(rcond=1e-3).prepare(panel)
    factor = reduction._inverse_cholesky
    assert factor.shape == (reduction.degrees_of_freedom, 32)
    assert torch.equal(factor, torch.triu(factor))
    p = panel.double()
    pinv = torch.linalg.pinv(p.T @ p / p.shape[0], rtol=1e-3, hermitian=True)
    torch.testing.assert_close(factor.T @ factor, pinv, rtol=1e-9, atol=1e-9)
    # Row block [s, e) of an upper-trapezoidal factor needs only columns [s, K).
    assert all(first == start and last == 32 for start, _end, first, last, _rows in reduction._blocks)


def test_eigen_truncation_matches_numpy_pinv():
    panel = _collinear_panel()
    p = panel.double().numpy()
    correlation = p.T @ p / p.shape[0]
    values = np.linalg.eigvalsh(correlation)
    reduction = JagwasReduction(rcond=1e-3)
    with pytest.warns(UserWarning, match='eigen-directions'):
        reduction.prepare(panel)
    rank = int((values > 1e-3 * values.max()).sum())
    report = reduction.rank_report()
    assert reduction.degrees_of_freedom == rank == report['rank'] < 32
    assert report['method'] == 'eigen' and report['bound'] == 'rcond' and report['dropped'] == []
    assert report['rcond'] == 1e-3 and report['largest_eigenvalue'] == pytest.approx(values.max())
    t = torch.as_tensor(np.random.default_rng(8).standard_normal((64, 32)) * 3, dtype=torch.float32)
    z = _score(t, _df(64, 2998)).double().numpy()
    expected = np.einsum('ij,jk,ik->i', z, np.linalg.pinv(correlation, rcond=1e-3, hermitian=True), z)
    np.testing.assert_allclose(_reduce(reduction, t, 2998).double().numpy(), expected, rtol=1e-5)


def test_eigen_truncation_nests_and_stays_within_the_rounding_target():
    rng = np.random.default_rng(12)
    base = rng.standard_normal((3000, 30))
    # One trait within 1e-4 of a combination of two others: an eigenvalue near 1e-8.
    panel = _standardize(np.column_stack([base, base[:, 0] + base[:, 1] + 1e-4 * rng.standard_normal(3000)]))
    t = torch.as_tensor(rng.standard_normal((64, 31)) * 3, dtype=torch.float32)
    statistics = []
    for rcond in (1e-1, 1e-2, 1e-3):
        reduction = JagwasReduction(rcond=rcond)
        reduction.prepare(panel)
        statistics.append(_reduce(reduction, t, 2998).double())
    # The same eigenvectors at every level, so a lower cutoff only adds directions.
    assert all((later >= earlier * (1 - 1e-5)).all() for earlier, later in zip(statistics, statistics[1:]))
    tiny = JagwasReduction(rcond=1e-15)
    with pytest.warns(UserWarning, match='keeping 30 of 31'):
        tiny.prepare(panel)
    assert tiny.rank_report()['bound'] == 'rounding' and tiny.rank_report()['rounding_error'] <= 0.01


def test_min_residual_drops_traits_mostly_explained_by_the_others():
    rng = np.random.default_rng(21)
    base = rng.standard_normal((3000, 30))
    # About 0.0025 of trait 30's variance is its own given trait 0: VIF about 400.
    panel = _standardize(np.column_stack([base, base[:, 0] + 0.05 * rng.standard_normal(3000)]))
    loose = JagwasReduction(min_residual=1e-3)
    loose.prepare(panel)
    assert loose.degrees_of_freedom == 31 and loose.rank_report()['method'] == 'cholesky'
    strict = JagwasReduction(min_residual=1e-2)
    with pytest.warns(UserWarning, match='less than 0.01 of each'):
        strict.prepare(panel)
    report = strict.rank_report()
    assert report['rank'] == 30 and report['bound'] == 'min_residual' and report['min_residual'] == 1e-2
    (dropped,) = report['dropped']
    assert dropped['index'] in (0, 30) and dropped['vif'] == pytest.approx(400, rel=0.15)
    t = torch.as_tensor(rng.standard_normal((64, 31)) * 3, dtype=torch.float32)
    kept = [index for index in range(31) if index != dropped['index']]
    expected = _exact_quadratic_form(panel[:, kept], _score(t, _df(64, 2998))[:, kept])
    torch.testing.assert_close(_reduce(strict, t, 2998).double(), expected, rtol=1e-5, atol=1e-4)
    with pytest.raises(ValueError, match='not both'):
        JagwasReduction(rcond=1e-3, min_residual=1e-2)
    groups = JagwasGroups([('a', range(31), {'min_residual': 1e-2}), ('b', range(31), 1e-3)])
    assert [(r.min_residual, r.rcond) for r in groups.reductions] == [(1e-2, None), (None, 1e-3)]
    with pytest.raises(ValueError, match='unknown cutoff'):
        JagwasGroups([('a', range(3), {'vif': 10})])


def test_groups_take_their_own_cutoff():
    panel = _collinear_panel()
    groups = JagwasGroups([('rounding', range(32), {'rcond': 0}), ('eigen', range(32), 1e-3)])
    with pytest.warns(UserWarning):
        groups.prepare(panel)
    rounding, eigen = JagwasReduction(rcond=0), JagwasReduction(rcond=1e-3)
    with pytest.warns(UserWarning):
        rounding.prepare(panel)
        eigen.prepare(panel)
    t = torch.as_tensor(np.random.default_rng(9).standard_normal((32, 32)) * 3, dtype=torch.float32)
    statistic = groups.reduce(torch.zeros_like(t), t, torch.zeros(32, dtype=torch.uint8), _df(32, 2998), 2)[1]
    torch.testing.assert_close(statistic[:, 0], _reduce(rounding, t, 2998), rtol=1e-6, atol=0)
    torch.testing.assert_close(statistic[:, 1], _reduce(eigen, t, 2998), rtol=1e-6, atol=0)
    assert [report['method'] for report in groups.rank_report()] == ['pivoted_cholesky', 'eigen']
    assert groups.spawn().reductions[1].rcond == 1e-3
    with pytest.raises(ValueError, match='rcond'):
        JagwasReduction(rcond=1.5)


def test_groups_refuse_bad_definitions():
    with pytest.raises(ValueError, match='unique'):
        JagwasGroups([('a', [0]), ('a', [1])])
    with pytest.raises(ValueError, match='distinct'):
        JagwasGroups({'a': [0, 0]})
    with pytest.raises(ValueError, match='outside the 3 scanned traits'):
        JagwasGroups({'a': [0, 3]}).prepare(_panel(100, 3))




# --- Public API (single device) --------------------------------------------


def _api_inputs(samples=240, variants=30, traits=5, seed=4):
    rng = np.random.default_rng(seed)
    genotype = rng.integers(0, 3, size=(samples, variants)).astype(np.float64)
    phenotype = rng.normal(size=(samples, traits))
    phenotype[:, 0] += 0.4 * genotype[:, 3]  # one real association
    return genotype, phenotype


def _score_reference(genotype, phenotype):
    """Independent FP64 OLS (intercept only) and the score-form quadratic form."""
    g = genotype - genotype.mean(0)
    y = phenotype - phenotype.mean(0)
    df = len(y) - 2
    correlation = np.corrcoef(y.T)
    values = []
    for column in g.T:
        ss = column @ column
        beta = column @ y / ss
        residual = y - column[:, None] * beta
        t = beta / np.sqrt((residual * residual).sum(0) / df / ss)
        z = t / np.sqrt(1 + t * t / df)
        values.append(float(z @ np.linalg.solve(correlation, z)))
    return np.asarray(values)


def _api_run(genotype, phenotype, output, **extra):
    from test_reduce_modes import _ArrayStream
    from torchgwas.api import run_linear_gwas
    return run_linear_gwas(_ArrayStream(genotype), phenotype, reduce='jagwas', output_dir=output,
                           chunk_size=8, device='cpu', compute_dtype='float64', **extra)


def _api_statistic(output):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    manifest, parts = open_indexed_sumstats(output / 'sumstats')
    values = np.full(manifest['shape'][0], np.nan)
    for part in parts:
        values[np.asarray(part['variant_index'])] = np.asarray(part['chi2'])
    return manifest, values


def test_api_matches_an_independent_fp64_score_reference(tmp_path):
    genotype, phenotype = _api_inputs()
    result = _api_run(genotype, phenotype, tmp_path / 'out')
    manifest, values = _api_statistic(tmp_path / 'out')
    assert manifest['df'] == 5 and manifest['jagwas_rank']['method'] == 'eigen'  # the default cutoff
    assert result.run_metadata['sumstats_write']['jagwas_rank'] == manifest['jagwas_rank']
    np.testing.assert_allclose(values, _score_reference(genotype, phenotype), rtol=1e-6, atol=1e-9)


def test_api_drops_a_duplicated_trait_and_reports_the_kept_rank(tmp_path):
    genotype, phenotype = _api_inputs()
    _api_run(genotype, phenotype, tmp_path / 'independent', jagwas_rcond=0)
    with pytest.warns(UserWarning, match='1 of 6 traits are collinear'):
        _api_run(genotype, np.column_stack([phenotype, phenotype[:, 2]]), tmp_path / 'duplicated', jagwas_rcond=0)
    reference, independent = _api_statistic(tmp_path / 'independent')
    manifest, duplicated = _api_statistic(tmp_path / 'duplicated')
    assert reference['df'] == 5 and manifest['df'] == 5 and manifest['shape'][1] == 6
    rank = manifest['jagwas_rank']
    assert rank['rank'] == 5 and rank['traits'] == 6 and rank['method'] == 'pivoted_cholesky'
    assert rank['rounding_error'] <= rank['target'] == 0.01
    assert [item['trait'] for item in rank['dropped']] in (['trait_2'], ['trait_5'])
    np.testing.assert_allclose(duplicated, independent, rtol=1e-9, atol=0)


def test_rank_callback_reports_df_before_the_first_result_chunk(tmp_path):
    genotype, phenotype = _api_inputs()
    events = []
    _api_run(genotype, np.column_stack([phenotype, phenotype[:, 1] + phenotype[:, 3]]), tmp_path / 'callback',
             sumstats_format='none', result_chunk_callback=lambda chunk: events.append(('chunk', int(chunk[0]))),
             jagwas_rank_callback=lambda report: events.append(('rank', report['rank'])))
    assert events[0] == ('rank', 5) and [kind for kind, _ in events].count('rank') == 1
    assert [value for kind, value in events if kind == 'chunk'] == [0, 8, 16, 24]
    with pytest.raises(TypeError, match='jagwas_rank_callback'):
        from test_reduce_modes import _ArrayStream
        from torchgwas.api import run_linear_gwas
        run_linear_gwas(_ArrayStream(genotype), phenotype, reduce='significant', output_dir=tmp_path / 'bad',
                        device='cpu', jagwas_rank_callback=print)


def test_api_groups_match_separate_scans(tmp_path):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    genotype, phenotype = _api_inputs(traits=6)
    duplicated = np.column_stack([phenotype[:, 3:], phenotype[:, 4]])
    reports = []
    with pytest.warns(UserWarning, match='jagwas group b: 1 of 4 traits'):
        result = _api_run(genotype, np.column_stack([phenotype[:, :3], duplicated]), tmp_path / 'grouped',
                          jagwas_groups=[('a', [0, 1, 2]), ('b', [3, 4, 5, 6])], jagwas_rank_callback=reports.append,
                          jagwas_rcond=0)
    _api_run(genotype, phenotype[:, :3], tmp_path / 'a', jagwas_rcond=0)
    with pytest.warns(UserWarning):
        _api_run(genotype, duplicated, tmp_path / 'b', jagwas_rcond=0)
    manifest, parts = open_indexed_sumstats(tmp_path / 'grouped' / 'sumstats')
    values = np.full((manifest['shape'][0], 2), np.nan)
    for part in parts:
        values[np.asarray(part['variant_index'])] = np.asarray(part['chi2'])
    assert manifest['groups'] == ['a', 'b'] and manifest['df'] == [3, 3]
    assert reports == [manifest['jagwas_rank']] and result.run_metadata['jagwas_groups'] == ['a', 'b']
    rank_a, rank_b = manifest['jagwas_rank']
    assert rank_a['group'] == 'a' and rank_a['rank'] == 3 and rank_b['rank'] == 3 and rank_b['traits'] == 4
    assert [item['trait'] for item in rank_b['dropped']] in (['trait_4'], ['trait_6'])
    for column, name in enumerate('ab'):
        separate, alone = _api_statistic(tmp_path / name)
        assert separate['df'] == 3
        np.testing.assert_allclose(values[:, column], alone, rtol=1e-12, atol=0)


def test_groups_are_residualised_as_if_alone():
    from torchgwas.preprocess import residualize_and_standardize
    rng = np.random.default_rng(3)
    phenotype = rng.normal(size=(300, 9)).astype(np.float32)
    covariates = rng.normal(size=(300, 4))
    grouped, _ = residualize_and_standardize(phenotype, covariates, column_groups=[[0, 1, 2], [5, 6]])
    for columns in ([0, 1, 2], [5, 6], [3, 4, 7, 8]):  # the last: columns in no group
        alone, _ = residualize_and_standardize(phenotype[:, columns], covariates)
        np.testing.assert_array_equal(grouped[:, columns], alone)


def test_api_eigen_truncation_ignores_a_duplicated_trait(tmp_path):
    genotype, phenotype = _api_inputs()
    _api_run(genotype, phenotype, tmp_path / 'independent')
    with pytest.warns(UserWarning, match='keeping 5 of 6 eigen-directions'):
        result = _api_run(genotype, np.column_stack([phenotype, phenotype[:, 2]]), tmp_path / 'eigen',
                          jagwas_rcond=1e-3)
    _, independent = _api_statistic(tmp_path / 'independent')
    manifest, eigen = _api_statistic(tmp_path / 'eigen')
    assert manifest['df'] == 5 and manifest['jagwas_rank']['method'] == 'eigen'
    assert result.run_metadata['jagwas_rcond'] == 1e-3
    np.testing.assert_allclose(eigen, independent, rtol=1e-8, atol=0)


def test_api_jagwas_takes_missing_phenotype_rows(tmp_path):
    from torchgwas.preprocess import residualize_and_standardize
    genotype, phenotype = _api_inputs()
    masked = phenotype.copy()
    masked[[3, 50, 200]] = np.nan  # whole rows, as the pipeline drops an outlier's phenotype
    _api_run(genotype, masked, tmp_path / 'masked')
    manifest, values = _api_statistic(tmp_path / 'masked')
    # The scan regresses the mean-imputed panel; its R is that panel's Gram, so
    # the reference is the ordinary score statistic of the imputed panel.
    imputed, _ = residualize_and_standardize(masked, None)
    assert manifest['df'] == 5
    np.testing.assert_allclose(values, _score_reference(genotype, imputed), rtol=1e-6, atol=1e-9)


def test_api_refuses_groups_without_jagwas(tmp_path):
    from test_reduce_modes import _ArrayStream
    from torchgwas.api import run_linear_gwas
    genotype, phenotype = _api_inputs()
    with pytest.raises(TypeError, match='jagwas_groups'):
        run_linear_gwas(_ArrayStream(genotype), phenotype, reduce='significant', output_dir=tmp_path / 'bad',
                        device='cpu', jagwas_groups={'a': [0, 1]})


def test_api_groups_follow_phenotype_qc(tmp_path):
    # A constant trait in group a is removed by phenotype QC; the groups index
    # the input panel, so b must stay b (before: IndexError past the last column).
    genotype, phenotype = _api_inputs(traits=6)
    panel = phenotype.copy()
    panel[:, 1] = 1.0
    with pytest.warns(UserWarning, match='jagwas group a: 1 of 3 traits removed by phenotype QC'):
        _api_run(genotype, panel, tmp_path / 'grouped', jagwas_groups=[('a', [0, 1, 2]), ('b', [3, 4, 5])])
    _api_run(genotype, phenotype[:, 3:], tmp_path / 'b')
    import json
    manifest = json.loads((tmp_path / 'grouped' / 'sumstats' / 'manifest.json').read_text())
    assert manifest['df'] == [2, 3] and manifest['groups'] == ['a', 'b']
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    _, parts = open_indexed_sumstats(tmp_path / 'grouped' / 'sumstats')
    grouped = np.full(manifest['shape'][0], np.nan)
    for part in parts:
        grouped[np.asarray(part['variant_index'])] = np.asarray(part['chi2'])[:, 1]
    _, alone = _api_statistic(tmp_path / 'b')
    np.testing.assert_allclose(grouped, alone, rtol=1e-6, atol=0)
