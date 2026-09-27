"""The prospective real-input fixture retains raw records across index blocks."""
import os
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks'))
from direct_real_pgen_longrun_20260921 import copy_record_prefix
from torchgwas.pgen_reader import read_header
from torchgwas.pgen_native_reader import NativePgenReader
from torchgwas import pgen_native


@pytest.fixture(scope='module')
def original(tmp_path_factory):
    pgenlib=pytest.importorskip('pgenlib')
    path=tmp_path_factory.mktemp('raw_prefix')/'source.pgen'
    n,m=257,65540
    rng=np.random.default_rng(518)
    base=rng.binomial(2,.3,size=n).astype(np.int8)
    with pgenlib.PgenWriter(os.fsencode(path),n,variant_ct=m,nonref_flags=False) as writer:
        for at in range(0,m,1024):
            count=min(1024,m-at);block=np.tile(base,(count,1))
            for i in range(count):block[i,(at+i)%n]=(at+i)%3
            writer.append_biallelic_batch(block)
    return path,n,m


@pytest.mark.parametrize('count',[1,7,65536,65539,65540])
@pytest.mark.skipif(not pgen_native.available(),reason='native PGEN decoder is not built')
def test_raw_prefix_standard_index_and_both_readers(tmp_path,original,count):
    import pgenlib
    source,n,m=original;destination=tmp_path/'prefix.pgen'
    provenance=copy_record_prefix(source,destination,count)
    h=read_header(destination)
    assert h.sample_ct==n and h.variant_ct==count
    assert provenance['records_byte_identical']
    with pgenlib.PgenReader(os.fsencode(source)) as expected,pgenlib.PgenReader(os.fsencode(destination)) as reference,NativePgenReader(destination) as native:
        for at in sorted({0,count-1,max(0,count-3),min(65535,count-1),min(65536,count-1)}):
            hi=min(count,at+3)
            want=np.empty((hi-at,n),np.int8);a=np.empty_like(want);b=np.empty_like(want)
            expected.read_range(at,hi,want);reference.read_range(at,hi,a);native.read_range(at,hi,b)
            np.testing.assert_array_equal(a,want);np.testing.assert_array_equal(b,want)
