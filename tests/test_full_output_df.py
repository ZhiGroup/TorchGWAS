"""Binary p-value reconstruction must use the df actually used by the scan."""
import numpy as np
import pytest
from scipy import special

from torchgwas.api import run_linear_gwas
from torchgwas.bed import PlinkBedGenotype
from torchgwas.linear import linear_scan_streaming_chunks
from torchgwas.sumstats import open_binary_sumstats
from test_statistics import _write_bed


@pytest.mark.parametrize('device', ['cpu', 'cuda:1'])
def test_binary_df_reconstructs_scan_probabilities_with_missing_calls(tmp_path, device):
    import torch
    if device.startswith('cuda') and (not torch.cuda.is_available() or torch.cuda.device_count()<2):
        pytest.skip('two CUDA devices required')
    rng=np.random.default_rng(9317)
    n,m,k=67,11,3
    genotype=rng.integers(0,3,size=(n,m)).astype(float)
    if device!='cpu':
        for column in range(m):
            genotype[:3*column,column]=np.nan
    phenotype=rng.normal(size=(n,k))
    phenotype[:,0]+=np.nan_to_num(genotype[:,6],nan=1.)*.8
    covariates=rng.normal(size=(n,2))
    bed=_write_bed(tmp_path/'missing',genotype)
    source=PlinkBedGenotype(bed,reader_workers=2,prefetch_chunks=2)
    chunks,_=linear_scan_streaming_chunks(source,phenotype,covariates,chunk_size=4,
                                         device=device,compute_dtype='float32')
    direct=list(chunks)
    expected_p=np.concatenate([row[4] for row in direct])
    out=tmp_path/'output'
    run_linear_gwas(source,phenotype,covariates,output_dir=out,device=device,
                    compute_dtype='float32',chunk_size=4,reader_workers=2,prefetch_chunks=2,
                    sumstats_fsync=True)
    _,tstat,manifest=open_binary_sumstats(out/'sumstats')
    assert isinstance(manifest['df'],dict), 'binary output lost variant-specific degrees of freedom'
    from torchgwas.sumstats import open_binary_df
    df=open_binary_df(out/'sumstats')
    np.testing.assert_array_equal(np.asarray(df).ravel(), np.sum(np.isfinite(genotype),axis=0)-4)
    restored=2*special.stdtr(df,-np.abs(tstat))
    np.testing.assert_allclose(restored,expected_p,rtol=8e-5,atol=2e-7)
