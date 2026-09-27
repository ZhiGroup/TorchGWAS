"""Candidate selector parity at representable threshold boundaries and layouts."""
import numpy as np
import pytest
from torchgwas.host_significance import select_host_pairs, predicate_block_shape


def baseline(beta, values, row_df, critical, start=0):
    keep = np.isfinite(values) & (np.abs(values) >= critical)
    rows, columns = np.nonzero(keep)
    selected_df = np.broadcast_to(row_df, values.shape)[rows, columns]
    return (rows.astype(np.int64) + start, columns.astype(np.int64),
            beta[rows, columns], values[rows, columns], selected_df)


ALGORITHMS = dict(production=select_host_pairs)


@pytest.fixture(params=['numpy', 'native'], autouse=True)
def predicate_backend(request, monkeypatch):
    monkeypatch.setenv('TORCHGWAS_HOST_PREDICATE', request.param)
    return request.param


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('layout', ['c', 'f', 'reversed', 'strided'])
@pytest.mark.parametrize('df_layout', ['scalar', 'row', 'trait', 'pair'])
def test_full_selection_preserves_values_df_indices_and_ownership(dtype, layout, df_layout):
    rng = np.random.default_rng(9219171)
    values = rng.standard_normal((17, 38)).astype(dtype)
    values.flat[::19] = np.nan
    values.flat[1::19] = np.inf
    values.flat[2::19] = -np.inf
    if layout == 'f': values = np.asfortranarray(values)
    elif layout == 'reversed': values = values[::-1, ::-1]
    elif layout == 'strided': values = values[:, ::2]
    beta = np.arange(values.size, dtype=dtype).reshape(values.shape)
    df = {'scalar': np.asarray(31.), 'row': np.arange(1, 18.)[:, None],
          'trait': np.arange(1, values.shape[1]+1.)[None, :],
          'pair': np.arange(1, values.size+1.).reshape(values.shape)}[df_layout]
    critical = np.broadcast_to(1.234567890123, df.shape)
    expected = baseline(beta, values, df, critical, 123)
    for function in ALGORITHMS.values():
        actual = function(beta, values, df, critical, 123)
        for got, want in zip(actual, expected):
            np.testing.assert_array_equal(got, want)
            assert not np.shares_memory(got, values)
            assert not np.shares_memory(got, beta)
            assert not np.shares_memory(got, df)
        assert actual[0].dtype == actual[1].dtype == np.int64
        assert not np.shares_memory(actual[0], actual[1])


@pytest.mark.parametrize('dtype', [np.float32, np.float64])
def test_exact_boundary_extremes(dtype):
    critical = np.array([0., 1e-50, 1e-40, 1., np.nextafter(1., 2.),
        7.123456789, np.finfo(np.float32).max, 1e40, np.inf, np.nan])[:, None]
    with np.errstate(over='ignore', invalid='ignore'):
        center = critical.astype(dtype)
        values = np.concatenate([center, np.nextafter(center, dtype(-np.inf)),
            np.nextafter(center, dtype(np.inf)), -center,
            -np.nextafter(center, dtype(-np.inf)), -np.nextafter(center, dtype(np.inf))], axis=1)
    beta = np.arange(values.size, dtype=dtype).reshape(values.shape)
    df = np.arange(1, len(critical)+1, dtype=dtype)[:, None]
    expected = baseline(beta, values, df, critical)
    for function in ALGORITHMS.values():
        for got, want in zip(function(beta, values, df, critical), expected):
            np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize('shape', [(0, 7), (3, 0), (1, 1), (1, (1 << 20)+3), (513, 2049)])
@pytest.mark.parametrize('retained', [False, True])
def test_empty_dense_and_block_tails(shape, retained):
    values = np.full(shape, 3. if retained else .5, np.float32)
    beta = np.ones_like(values)
    df = np.arange(shape[0], dtype=np.float32)[:, None] + 1
    critical = np.ones((shape[0], 1), np.float64)
    expected = baseline(beta, values, df, critical, 23)
    for function in ALGORITHMS.values():
        for got, want in zip(function(beta, values, df, critical, 23), expected):
            np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize('shape', [(0, 7), (3, 0), (513, 2049), (3, (1 << 20)+13)])
def test_predicate_temporaries_obey_the_declared_bound(monkeypatch, shape, predicate_backend):
    from torchgwas import host_significance as module
    original = np.abs
    sizes = []
    def observed(values):
        sizes.append(values.size)
        return original(values)
    monkeypatch.setattr(module.np, 'abs', observed)
    values = np.zeros(shape, np.float32)
    select_host_pairs(values, values, np.asarray(31.), np.asarray([1.]))
    if predicate_backend == 'native':
        assert sizes == []
        return
    height, width, calls = predicate_block_shape(*shape)
    assert len(sizes) == calls
    assert sum(sizes) == values.size
    assert all(0 < size <= module.PREDICATE_MAX_CELLS for size in sizes)


def test_owned_selected_arrays_survive_reused_dense_ring():
    values = np.ones((3, 5), np.float32)
    beta = np.arange(15, dtype=np.float32).reshape(3, 5)
    df = np.full((3, 1), 31., np.float32)
    result = select_host_pairs(beta, values, df, np.asarray([0.]), 19)
    expected = [a.copy() for a in result]
    beta.fill(-99); values.fill(-99); df.fill(-99)
    for got, want in zip(result, expected):
        np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize('value_layout',['c','f','reverse'])
@pytest.mark.parametrize('beta_layout',['c','f','reverse','wider','none'])
def test_independent_payload_layouts_keep_flat_coordinates_and_owned_results(value_layout,beta_layout):
    values=np.arange(35,dtype=np.float32).reshape(5,7)-17
    if value_layout=='f':values=np.asfortranarray(values)
    if value_layout=='reverse':values=values[:,::-1]
    beta=np.arange(35 if beta_layout!='wider' else 55,dtype=np.float32).reshape(5,-1)+100
    if beta_layout=='f':beta=np.asfortranarray(beta)
    if beta_layout=='reverse':beta=beta[::-1,::-1]
    actual_beta=None if beta_layout=='none' else beta
    df=np.arange(5,dtype=np.float32)[::-1,None]+20
    critical=np.asarray(6.123456789,np.float64)
    expected=list(baseline(beta,values,df,critical,31))
    if beta_layout=='none':expected[2]=None
    result=select_host_pairs(actual_beta,values,df,critical,31)
    values.fill(-99);beta.fill(-99);df.fill(-99)
    for got,want in zip(result,expected):
        if want is None:assert got is None
        else:np.testing.assert_array_equal(got,want)


def test_noncontiguous_dense_inputs_are_never_flattened_by_copy():
    class NoFlatten(np.ndarray):
        def reshape(self,*a,**kw):raise AssertionError('Noncontiguous input flattened')
        def ravel(self,*a,**kw):raise AssertionError('Noncontiguous input flattened')
        def flatten(self,*a,**kw):raise AssertionError('Noncontiguous input flattened')
    values=np.arange(70,dtype=np.float32).reshape(5,14)[:,::2].view(NoFlatten)
    beta=(np.arange(70,dtype=np.float32).reshape(5,14)+100)[:,::-2].view(NoFlatten)
    df=np.arange(5,dtype=np.float32)[:,None]+20
    critical=np.asarray(20.,np.float64)
    expected=baseline(np.asarray(beta),np.asarray(values),df,critical)
    for got,want in zip(select_host_pairs(beta,values,df,critical),expected):
        # The input guard must not also intercept the assertion helper's
        # legitimate flattening of its small, already selected result.
        np.testing.assert_array_equal(np.asarray(got),want)
