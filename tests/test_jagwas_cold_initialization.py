"""Cold multi-GPU linalg initialization and correct factor failure reporting."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import subprocess
import sys
import threading

import numpy as np
import pytest
import torch

import torchgwas.jagwas_projection as reduction


def test_concurrent_first_use_loads_once_then_factors_need_no_lock(monkeypatch):
    monkeypatch.setattr(reduction,'_CUDA_LINALG_READY',False)
    monkeypatch.setattr(torch.cuda,'device',lambda device:nullcontext())
    monkeypatch.setattr(torch,'eye',lambda *a,**k:object())
    start=threading.Barrier(4);entered=threading.Event();release=threading.Event();calls=[]
    def cold_call(matrix):
        calls.append(matrix);entered.set()
        assert release.wait(5), 'Initialization test did not release the loader'
    monkeypatch.setattr(torch.linalg,'cholesky',cold_call)
    def worker(index):
        start.wait(timeout=5)
        reduction._ensure_cuda_linalg_initialized(torch.device('cuda:'+str(index%2)))
    with ThreadPoolExecutor(4) as executor:
        futures=[executor.submit(worker,index) for index in range(4)]
        try:
            assert entered.wait(5)
            assert not any(future.done() for future in futures)
        finally:
            release.set()
        for future in futures:future.result(timeout=5)
    assert len(calls)==1 and reduction._CUDA_LINALG_READY
    # A warmed call must not acquire the initialization lock again.
    with reduction._CUDA_LINALG_LOCK:
        reduction._ensure_cuda_linalg_initialized(torch.device('cuda:0'))
    assert len(calls)==1


def test_failed_initialization_remains_retryable_and_preserves_runtime_error(monkeypatch):
    monkeypatch.setattr(reduction,'_CUDA_LINALG_READY',False)
    monkeypatch.setattr(torch.cuda,'device',lambda device:nullcontext())
    monkeypatch.setattr(torch,'eye',lambda *a,**k:object())
    error=RuntimeError('CUDA library cannot be loaded');calls=[]
    def failing(matrix):
        calls.append(matrix)
        if len(calls)==1:raise error
    monkeypatch.setattr(torch.linalg,'cholesky',failing)
    with pytest.raises(RuntimeError) as raised:
        reduction._ensure_cuda_linalg_initialized(torch.device('cuda:0'))
    assert raised.value is error and not reduction._CUDA_LINALG_READY
    reduction._ensure_cuda_linalg_initialized(torch.device('cuda:1'))
    assert reduction._CUDA_LINALG_READY and len(calls)==2


def test_cpu_and_meta_do_not_initialize_cuda(monkeypatch):
    monkeypatch.setattr(reduction,'_CUDA_LINALG_READY',False)
    def forbidden(*args,**kwargs):raise AssertionError('Unexpected CUDA initialization')
    monkeypatch.setattr(torch.cuda,'device',forbidden)
    for device in ('cpu','meta'):
        reduction._ensure_cuda_linalg_initialized(torch.device(device))
    assert not reduction._CUDA_LINALG_READY


@pytest.mark.parametrize('error',[RuntimeError('lazy wrapper should be called at most once'),
                                 torch.OutOfMemoryError('CUDA allocation failed')])
def test_non_numerical_cholesky_error_is_not_labeled_collinearity(monkeypatch,error):
    def fail(matrix):raise error
    monkeypatch.setattr(torch.linalg,'cholesky_ex',fail)
    with pytest.raises(type(error)) as raised:
        reduction.JagwasReduction(rcond=0).prepare(np.ones((32,2),np.float32),device='cpu')
    assert raised.value is error


@pytest.mark.parametrize('error',[RuntimeError('lazy wrapper should be called at most once'),
                                 torch.OutOfMemoryError('CUDA allocation failed')])
def test_non_numerical_eigh_error_propagates(monkeypatch,error):
    def fail(matrix):raise error
    monkeypatch.setattr(torch.linalg,'eigh',fail)
    with pytest.raises(type(error)) as raised:
        reduction.JagwasReduction().prepare(np.ones((32,2),np.float32),device='cpu')
    assert raised.value is error


def test_singular_phenotype_still_has_specific_joint_error():
    with pytest.raises(ValueError,match='no jagwas trait has nonzero variance'):
        reduction.JagwasReduction().prepare(np.zeros((32,2),np.float32),device='cpu')


@pytest.mark.parametrize('repeat',range(2))
def test_two_gpu_first_factors_succeed_in_a_fresh_interpreter(repeat):
    if not torch.cuda.is_available() or torch.cuda.device_count()<2:
        pytest.skip('Two CUDA devices required')
    program='''
from concurrent.futures import ThreadPoolExecutor
import threading
import numpy as np
import torch
from torchgwas.reduce import JagwasReduction

torch.set_num_threads(2)
torch.set_num_interop_threads(1)
torch.backends.cuda.matmul.allow_tf32=False
rng=np.random.default_rng(922451)
y=rng.normal(size=(257,64)).astype(np.float32)
correlation=y.astype(np.float64).T@y.astype(np.float64)/len(y)
truth=np.linalg.pinv(correlation,rcond=1e-3,hermitian=True)
ready=threading.Barrier(2)
def prepare(index):
    device='cuda:'+str(index)
    torch.cuda.set_device(device)
    # CUDA context creation is deliberately complete, but no CUDA linalg has
    # run before both workers enter the production factor preparation.
    torch.empty(1,device=device)
    ready.wait(timeout=30)
    factor=JagwasReduction().prepare(y,device=device)._inverse_cholesky.cpu().numpy()
    np.testing.assert_allclose(factor.T@factor,truth,rtol=1e-9,atol=1e-9)
    return factor
with ThreadPoolExecutor(2) as pool:
    factors=list(pool.map(prepare,range(2)))
np.testing.assert_array_equal(*factors)
print('fresh concurrent factors verified')
'''
    result=subprocess.run([sys.executable,'-c',program],text=True,capture_output=True,timeout=90)
    assert result.returncode==0,result.stdout+'\n'+result.stderr
    assert 'fresh concurrent factors verified' in result.stdout
