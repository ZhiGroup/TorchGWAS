"""Reconcile source factor semantics with counted, executable scalar algebra."""
import numpy as np
import pytest
import torch
from unittest.mock import patch
from torchgwas.reduce import JagwasReduction
from torchgwas.reduction_tensor_work import jagwas_factor_arithmetic_work


def reference(matrix):
    k=len(matrix);lower=np.zeros((k,k),dtype=np.float64)
    counts=dict(chol_flops=0,chol_div=0,chol_sqrt=0,solve_flops=0,solve_div=0)
    for column in range(k):
        for row in range(column,k):
            value=float(matrix[row,column])
            for previous in range(column):
                value-=lower[row,previous]*lower[column,previous]
                counts['chol_flops']+=2
            if row==column:
                lower[row,column]=np.sqrt(value);counts['chol_sqrt']+=1
            else:
                lower[row,column]=value/lower[column,column];counts['chol_div']+=1
    inverse=np.zeros_like(lower)
    for rhs in range(k):
        for row in range(k):
            value=float(row==rhs)
            for previous in range(row):
                value-=lower[row,previous]*inverse[previous,rhs];counts['solve_flops']+=2
            inverse[row,rhs]=value/lower[row,row];counts['solve_div']+=1
    return inverse,counts


@pytest.mark.parametrize('k',[1,2,7,17])
@pytest.mark.parametrize('dtype',['float32','float64'])
def test_counted_dense_factor_matches_source_and_ledger(k,dtype):
    n=2*k+9
    y=np.random.default_rng(912).normal(size=(n,k)).astype(dtype)
    # Use the source's correlation precision/order, then scalar FP64 algebra.
    # The source's correlation: the FP64 Gram of the scan-precision values.
    matrix=torch.as_tensor(y).double();correlation=((matrix.T@matrix)/float(n)).numpy()
    inverse,counts=reference(correlation)
    source=JagwasReduction(rcond=0).prepare(y,device='cpu')._inverse_cholesky.numpy()
    np.testing.assert_allclose(inverse,source,rtol=1e-12,atol=1e-12)
    # Cholesky reads the lower triangle. Independent FP32 Gram entries can
    # differ slightly across the two triangles; do not change the tolerance.
    symmetric=np.tril(correlation)+np.tril(correlation,-1).T
    np.testing.assert_allclose(inverse@symmetric@inverse.T,np.eye(k),rtol=1e-12,atol=1e-12)
    work=jagwas_factor_arithmetic_work(n,k,compute_dtype=dtype,method='rounding')
    assert work['cholesky']['multiply_add_flops']==counts['chol_flops']
    assert work['cholesky']['divisions']==counts['chol_div']
    assert work['cholesky']['square_roots']==counts['chol_sqrt']
    assert work['triangular_solve']['multiply_add_flops']==counts['solve_flops']
    assert work['triangular_solve']['divisions']==counts['solve_div']
    assert work['h2d_bytes']==y.nbytes
    assert work['persistent_factor_bytes']==source.nbytes
    assert work['correlation']['useful_flops']==2*n*k*k
    assert not work['prediction_complete']


def test_source_change_requires_reconciliation():
    from torchgwas.reduction_tensor_work import jagwas_tensor_work
    for method,dropped in [('rounding','aten.linalg_solve_triangular.default'),('eigen','aten.linalg_qr.default')]:
        changed=jagwas_tensor_work(257,1,7,phase='prepare',method=method)
        changed['steps']=[s for s in changed['steps'] if s['op']!=dropped]
        with patch('torchgwas.reduction_tensor_work.jagwas_tensor_work',return_value=changed):
            with pytest.raises(ValueError,match='source operations changed'):
                jagwas_factor_arithmetic_work(257,7,method=method)


@pytest.mark.parametrize('k',[1,2,7,17])
@pytest.mark.parametrize('dtype',['float32','float64'])
def test_eigen_factor_ledger_matches_source(k,dtype):
    n=2*k+9
    y=np.random.default_rng(913).normal(size=(n,k)).astype(dtype)
    matrix=torch.as_tensor(y).double();correlation=(matrix.T@matrix)/float(n)
    source=JagwasReduction().prepare(y,device='cpu')._inverse_cholesky
    # R'R is the truncated pseudo-inverse; R is upper trapezoidal (k x K, here k = K).
    assert torch.equal(source,torch.triu(source))
    torch.testing.assert_close(source.T@source,torch.linalg.pinv(correlation,rtol=1e-3,hermitian=True),
                               rtol=1e-9,atol=1e-9)
    work=jagwas_factor_arithmetic_work(n,k,compute_dtype=dtype)
    assert work['method']=='eigen' and work['correlation']['useful_flops']==2*n*k*k
    assert work['eigh']['multiply_add_flops']==(2*k**3)//3+k**3
    assert work['qr']['multiply_add_flops']==k**3-(k**3)//3
    assert work['scale']['square_roots']==k and work['spectrum']['d2h_bytes']==8*k
    assert work['h2d_bytes']==y.nbytes and work['persistent_factor_bytes']==8*k*k
    assert set(work['phase_logical_bytes'])=={'correlation','eigh','scale','qr'}|({'cast'} if dtype=='float32' else set())
    assert not work['prediction_complete']
