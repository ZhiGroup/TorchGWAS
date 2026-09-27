"""Reuse rank-validated covariates across full-output tiles and shards."""
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.preprocess import _covariate_basis,residualize_and_standardize


@pytest.mark.parametrize('kind',['absent','zero_rank','full_rank','dependent'])
@pytest.mark.parametrize('dtype',[np.float32,np.float64])
def test_prevalidated_basis_matches_recomputed_reference_and_stays_readonly(kind,dtype):
    rng=np.random.default_rng(951)
    y=rng.normal(size=(67,7)).astype(dtype);y[3,2]=np.nan
    c=None if kind=='absent' else rng.normal(size=(67,3)).astype(dtype)
    if kind=='zero_rank':c.fill(1)
    if kind=='dependent':c[:,2]=2*c[:,0]
    basis=None if c is None else _covariate_basis(c)
    expected=residualize_and_standardize(y,c,return_observed_counts=True)
    saved=None if basis is None else basis.copy()
    if basis is not None:basis.setflags(write=False)
    with patch('torchgwas.preprocess._covariate_basis',side_effect=AssertionError('recomputed basis')):
        actual=residualize_and_standardize(y,c,return_observed_counts=True,
            _prevalidated_covariate_basis=basis)
    np.testing.assert_array_equal(actual[0],expected[0])
    np.testing.assert_array_equal(actual[2],expected[2])
    assert actual[1] is basis
    if basis is not None:np.testing.assert_array_equal(basis,saved)


@pytest.mark.parametrize('shape',[(66,2),(67,4),(67,)])
def test_prevalidated_basis_dimensions_are_checked(shape):
    with pytest.raises(ValueError,match='Prevalidated covariate basis'):
        residualize_and_standardize(np.ones((67,7)),np.ones((67,3)),
            _prevalidated_covariate_basis=np.ones(shape))


@pytest.mark.parametrize('axis,devices',[
    ('trait',['cpu']),('trait',['cuda:1','cuda:2']),('variant',['cuda:1','cuda:2'])])
def test_full_output_computes_basis_once_and_all_workers_share_it(tmp_path,axis,devices):
    import torch
    from test_statistics import _write_bed
    from torchgwas.api import run_linear_gwas
    from torchgwas.bed import PlinkBedGenotype
    from torchgwas.sumstats import open_binary_sumstats,open_binary_df
    from torchgwas import linear
    if devices[0].startswith('cuda') and (not torch.cuda.is_available() or torch.cuda.device_count()<3):
        pytest.skip('three visible CUDA devices required')
    rng=np.random.default_rng(1952);n,m,k=97,19,7
    bed=_write_bed(tmp_path/'input',rng.integers(0,3,size=(n,m)).astype(float))
    y=rng.normal(size=(n,k)).astype(np.float32);c=rng.normal(size=(n,3)).astype(np.float32)
    kwargs=dict(device=devices[0],compute_dtype='float32',chunk_size=4,reader_workers=4,
        prefetch_chunks=4,sumstats_block_bytes=64,sumstats_queue_depth=1,sumstats_fsync=True)
    run_linear_gwas(PlinkBedGenotype(bed),y,c,output_dir=tmp_path/'reference',**kwargs)
    passed=[];original=linear.residualize_and_standardize
    def observe(*args,**kwargs):
        passed.append(kwargs['_prevalidated_covariate_basis'])
        return original(*args,**kwargs)
    partition=dict(trait_block=3,trait_devices=devices) if axis=='trait' else dict(variant_devices=devices)
    with patch('torchgwas.preprocess._covariate_basis',wraps=_covariate_basis) as basis:
        with patch('torchgwas.linear.residualize_and_standardize',side_effect=observe):
            run_linear_gwas(PlinkBedGenotype(bed),y,c,output_dir=tmp_path/'partition',**partition,**kwargs)
    assert basis.call_count==1
    assert len(passed)==(3 if axis=='trait' else len(devices))
    assert all(value is passed[0] for value in passed)
    ref_b,ref_t,_=open_binary_sumstats(tmp_path/'reference'/'sumstats')
    actual_b,actual_t,_=open_binary_sumstats(tmp_path/'partition'/'sumstats')
    for actual,expected in [(actual_b,ref_b),(actual_t,ref_t)]:
        np.testing.assert_allclose(np.asarray(actual),expected,rtol=3e-5,atol=3e-6)
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path/'partition'/'sumstats')),
        np.broadcast_to(open_binary_df(tmp_path/'reference'/'sumstats'),
            (m,k) if axis=='trait' else (m,1)))


def test_model_excludes_shared_api_basis_from_each_tile(input_path):
    from test_trait_tiling_model import candidate
    from torchgwas.trait_tiling_model import _prepare_graph
    from torchgwas.setup_work import setup_work
    p=candidate(input_path)['tiles'][0]['profile']
    p['process_units']['covariate_basis_work']=1e9
    graph,_=_prepare_graph(setup_work(32,2,8),p,1.,0,0)
    assert graph.nodes['covariate_basis'][0]==0.


from test_trait_tiling_model import input_path
