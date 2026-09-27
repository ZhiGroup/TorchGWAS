"""QC temporary reuse must preserve missingness, column selection and inputs."""
import tracemalloc
import numpy as np
import pytest
from torchgwas.preprocess import _phenotype_column_mask


@pytest.mark.parametrize('dtype',[np.float32,np.float64])
@pytest.mark.parametrize('layout',['c','fortran','strided'])
def test_qc_keeps_counts_masks_and_readonly_inputs(dtype,layout):
    x=np.zeros((12,6),dtype=dtype)
    x[:,0]=np.arange(12)/4
    x[:,1]=1
    x[:,2]=np.nan
    x[:,3]=np.nan;x[:2,3]=[1,2]
    x[:,4]=np.nan;x[0,4]=3
    x[:,5]=np.arange(12);x[::3,5]=np.nan
    if layout=='fortran':x=np.asfortranarray(x)
    if layout=='strided':
        storage=np.full((12,12),-123,dtype=dtype);storage[:,::2]=x;x=storage[:,::2]
    else:storage=x
    before=storage.copy();x.flags.writeable=False
    keep,counts=_phenotype_column_mask(x)
    np.testing.assert_array_equal(keep,[True,False,False,True,False,True])
    np.testing.assert_array_equal(counts,[12,12,0,2,1,8])
    np.testing.assert_array_equal(storage,before)


@pytest.mark.parametrize('value',[np.inf,-np.inf])
def test_qc_still_rejects_infinities_among_missing_values(value):
    x=np.array([[np.nan,0],[value,1],[1,2]])
    with pytest.raises(ValueError,match='infinite'):_phenotype_column_mask(x)


def test_qc_avoids_two_full_centered_temporaries():
    # Input is allocated before tracing. 13 MiB accommodates mask/reduction
    # buffers but not two simultaneous 8 MiB centered/squared temporaries.
    x=np.random.default_rng(923791).standard_normal((1024,1024),dtype=np.float32)
    tracemalloc.start()
    try:
        keep,counts=_phenotype_column_mask(x)
        _,peak=tracemalloc.get_traced_memory()
    finally:tracemalloc.stop()
    assert keep.all() and np.all(counts==1024)
    assert peak<13*1024**2