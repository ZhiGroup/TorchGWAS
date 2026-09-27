"""Output lifetime and memory regression for one-block GPU residualization."""
import gc
import numpy as np
import pytest
import torch
import torchgwas.preprocess as prep


@pytest.mark.skipif(torch.cuda.device_count()<2, reason='CUDA required')
@pytest.mark.parametrize('dtype',[np.float32,np.float64])
@pytest.mark.parametrize('covariates',[False,True])
@pytest.mark.parametrize('width',[1,17])
def test_download_owns_cpu_storage_without_extra_numpy_result(monkeypatch,dtype,covariates,width):
    rng=np.random.default_rng(13829)
    full=rng.normal(size=(128,2*width)).astype(dtype)
    y=full[:,::2];before=full.copy();y.flags.writeable=False
    cov=rng.normal(size=(128,3)).astype(dtype) if covariates else None
    q=prep._covariate_basis(cov) if covariates else None
    original=np.empty
    def refuse(shape,*args,**kwargs):
        if tuple(np.atleast_1d(shape))==y.shape:
            raise AssertionError('Redundant whole-phenotype output allocation')
        return original(shape,*args,**kwargs)
    with monkeypatch.context() as patch:
        patch.setattr(np,'empty',refuse)
        actual=prep._residualize_on_device(y,q,'cuda:1')
    expected=prep.residualize_and_standardize(y,cov,device='cpu')[0]
    np.testing.assert_allclose(actual,expected,rtol=1e-5 if dtype==np.float32 else 1e-12,atol=2e-6 if dtype==np.float32 else 1e-12)
    np.testing.assert_array_equal(full,before)
    assert actual.dtype==dtype and actual.flags.writeable and actual.flags.c_contiguous
    assert not np.shares_memory(actual,full)
    assert isinstance(actual.base,torch.Tensor) and actual.base.device.type=='cpu'
    saved=actual.copy();del y,full,q,cov;gc.collect()
    np.testing.assert_array_equal(actual,saved)
    actual[0,0]=123.
    assert actual[0,0]==123.


@pytest.mark.skipif(torch.cuda.device_count()<2, reason='two CUDA devices required')
def test_multiple_blocks_and_tail_still_assemble(monkeypatch):
    rng=np.random.default_rng(385)
    y=rng.normal(size=(128,17)).astype(np.float32)
    q=prep._covariate_basis(rng.normal(size=(128,3)).astype(np.float32))
    one=prep._residualize_on_device(y,q,'cuda:1')
    monkeypatch.setattr(prep,'_device_trait_block',lambda *args,**kwargs:7)
    many=prep._residualize_on_device(y,q,'cuda:1')
    np.testing.assert_allclose(many,one,rtol=1e-5,atol=2e-6)
    assert many.flags.owndata and many.flags.c_contiguous
    empty=prep._residualize_on_device(y[:,:0],q,'cuda:1')
    assert empty.shape==(128,0) and empty.dtype==np.float32
