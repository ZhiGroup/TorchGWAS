"""Only counts established by input QC may bypass phenotype rescanning."""
from types import SimpleNamespace
import numpy as np
import pytest
from torchgwas.preprocess import residualize_and_standardize


def test_prevalidated_complete_gpu_path_does_not_rescan_input(monkeypatch):
    import torchgwas.preprocess as module
    x=np.arange(45,dtype=np.float32).reshape(15,3)
    counts=np.full(3,15,dtype=np.int64);counts.flags.writeable=False
    monkeypatch.setattr(module,'_residualize_on_device',lambda p,q,d:p.copy())
    def forbidden(*args,**kwargs):raise AssertionError('Input was scanned again')
    monkeypatch.setattr(module.np,'isnan',forbidden)
    actual,q,used=module.residualize_and_standardize(x,None,device=SimpleNamespace(type='cuda'),
        return_observed_counts=True,_prevalidated_observed_counts=counts)
    np.testing.assert_array_equal(actual,x)
    np.testing.assert_array_equal(used,counts)
    assert q is None


@pytest.mark.parametrize('missing',[False,True])
@pytest.mark.parametrize('dtype',[np.float32,np.float64])
def test_cached_counts_preserve_cpu_results_and_inputs(missing,dtype):
    rng=np.random.default_rng(925361);x=rng.normal(size=(37,7)).astype(dtype)
    if missing:x[::3,::2]=np.nan
    cov=rng.normal(size=(37,2)).astype(dtype)
    counts=np.sum(~np.isnan(x),axis=0);before=x.copy()
    expected=residualize_and_standardize(x,cov,return_observed_counts=True,trait_block=3)
    actual=residualize_and_standardize(x,cov,return_observed_counts=True,trait_block=3,
        _prevalidated_observed_counts=counts)
    for a,b in zip(actual,expected):np.testing.assert_array_equal(a,b)
    np.testing.assert_array_equal(x,before)
    np.testing.assert_array_equal(counts,expected[2])


@pytest.mark.parametrize('counts',[[3,3],[3,3,3,3],[[3,3,3]],[3.,3.,3.],[True,True,True],[-1,3,3],[4,3,3]])
def test_invalid_prevalidated_counts_are_rejected(counts):
    with pytest.raises(ValueError,match='Prevalidated observed counts'):
        residualize_and_standardize(np.ones((3,3)),None,_prevalidated_observed_counts=counts)


def test_entirely_missing_prevalidated_trait_still_fails():
    with pytest.raises(ValueError,match='entirely missing'):
        residualize_and_standardize(np.array([[np.nan,1],[np.nan,2]]),None,
            _prevalidated_observed_counts=np.array([0,2]))


def test_failed_gpu_preparation_keeps_cpu_fallback(monkeypatch):
    import torchgwas.preprocess as module
    rng=np.random.default_rng(925369);x=rng.normal(size=(31,5))
    expected=residualize_and_standardize(x,None,return_observed_counts=True)
    def failed(*args):raise RuntimeError('simulated CUDA allocation failure')
    monkeypatch.setattr(module,'_residualize_on_device',failed)
    actual=residualize_and_standardize(x,None,device=SimpleNamespace(type='cuda'),return_observed_counts=True,
        _prevalidated_observed_counts=np.full(5,31))
    np.testing.assert_array_equal(actual[0],expected[0])
    np.testing.assert_array_equal(actual[2],expected[2])


def test_full_tile_api_passes_retained_count_slices_and_tail(tmp_path,monkeypatch):
    import torchgwas.api as api
    from test_statistics import _write_bed
    rng=np.random.default_rng(925376);n=41
    bed=_write_bed(tmp_path/'input',rng.integers(0,3,size=(n,9)).astype(float))
    y=rng.normal(size=(n,6));y[:,2]=1.
    original=api.linear_scan_streaming_chunks;seen=[]
    def record(*args,**kwargs):
        seen.append((args[1].shape[1],kwargs['_prevalidated_observed_counts'].copy()))
        return original(*args,**kwargs)
    monkeypatch.setattr(api,'linear_scan_streaming_chunks',record)
    api.run_linear_gwas(str(bed),y,genotype_format='plink',device='cpu',chunk_size=4,
        reader_workers=2,prefetch_chunks=2,trait_block=2,trait_devices=['cpu'],output_dir=tmp_path/'out')
    assert [width for width,counts in seen]==[2,2,1]
    for width,counts in seen:np.testing.assert_array_equal(counts,np.full(width,n))
    from torchgwas.sumstats import read_manifest
    assert read_manifest(tmp_path/'out'/'sumstats')['traits']==['trait_0','trait_1','trait_3','trait_4','trait_5']