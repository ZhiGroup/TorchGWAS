"""Packed PGEN native preparation against independent unpacked dosage."""
import unittest
import warnings

import numpy as np
import torch
from torchgwas.scan_gpu import prepare, available

@unittest.skipUnless(torch.cuda.is_available() and available(), 'native CUDA required')
class PackedPreparation(unittest.TestCase):
    def test_tail_padding_missing_and_stream(self):
        for n in (1,3,4,5,255,256,257,22250):
            rng=np.random.default_rng(n)
            codes=rng.integers(0,3,size=(5,n),dtype=np.uint8)
            codes[0,:]=0
            codes[1,:]=2
            codes[2,-1]=3
            width=((n+3)//4+63)//64*64
            packed=np.full((5,width),255,np.uint8)
            for i in range(n):
                shift=2*(i%4)
                packed[:,i//4]=(packed[:,i//4]&np.uint8(255^(3<<shift)))|(codes[:,i]<<shift)
            ref=codes.astype(np.float64);ref[codes==3]=np.nan
            # Missing calls are masked: excluded from the mean and the range,
            # and centred to exactly zero. A row with nothing observed (n==1
            # here, where the only call is missing) keeps an empty range.
            observed=~np.isnan(ref)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                mean=np.nanmean(ref,axis=1,keepdims=True)
                row_min=np.nanmin(ref,axis=1)
                row_max=np.nanmax(ref,axis=1)
            centered=np.where(observed,ref-mean,0.0).astype(np.float32)
            ss=(centered.astype(np.float64)**2).sum(axis=1).astype(np.float32)
            stream=torch.cuda.Stream()
            with torch.cuda.stream(stream):
                outputs=prepare(torch.tensor(packed,device='cuda'),encoding='pgen_2bit',n_samples=n)
            stream.synchronize()
            actual=[x.cpu().numpy() for x in outputs]
            np.testing.assert_allclose(actual[0],centered,atol=2e-7,rtol=2e-6,equal_nan=True)
            np.testing.assert_allclose(actual[1],ss,atol=2e-6,rtol=2e-6,equal_nan=True)
            np.testing.assert_allclose(actual[2],row_min,equal_nan=True)
            np.testing.assert_allclose(actual[3],row_max,equal_nan=True)
    def test_the_two_bit_tables_are_not_interchangeable(self):
        """PGEN and PLINK1 two-bit codes collide, so both tables need pinning.

        PGEN reads 0/1/2 as the ALT1 count and reserves 3 for missing; PLINK1
        reserves **1** for missing and maps 0/2/3 to 0/1/2. The same byte
        therefore decodes differently under each, with no error and entirely
        plausible numbers -- which is why the encoding is a compile-time
        template parameter reached through two distinct exported symbols
        (`tg_scan_prepare_pgen2`, `tg_scan_prepare_bed2`) and never a runtime
        flag. Only `pgen_2bit` had coverage; `plink_2bit` is what the BED scan
        actually runs.

        The assertions are on the observed count and the value range, because
        those separate the two tables where checking decoded values would not:
        codes 0,1,2,3 give the same min and max under both.
        """
        n=4
        cases={
            # codes        pgen (present,min,max)  plink1 (present,min,max)
            (1,1,2,2):     ((4,1.0,2.0),           (2,1.0,1.0)),
            (3,3,0,0):     ((2,0.0,0.0),           (4,0.0,2.0)),
            (0,1,2,3):     ((3,0.0,2.0),           (3,0.0,2.0)),
        }
        width=((n+3)//4+63)//64*64
        for codes,(want_pgen,want_plink) in cases.items():
            packed=np.zeros((1,width),np.uint8)
            for i,code in enumerate(codes):
                packed[0,i//4]|=np.uint8(code<<(2*(i%4)))
            for encoding,(present,low,high) in (('pgen_2bit',want_pgen),
                                                ('plink_2bit',want_plink)):
                with self.subTest(codes=codes,encoding=encoding):
                    out=prepare(torch.tensor(packed,device='cuda'),
                                encoding=encoding,n_samples=n)
                    self.assertEqual(int(out[4].cpu().numpy()[0]),present)
                    self.assertAlmostEqual(float(out[2].cpu().numpy()[0]),low,places=6)
                    self.assertAlmostEqual(float(out[3].cpu().numpy()[0]),high,places=6)

    def test_an_unknown_encoding_is_refused_rather_than_guessed(self):
        raw=torch.zeros((2,64),dtype=torch.uint8,device='cuda')
        for encoding in ('bed_2bit','plink1','PGEN_2BIT',''):
            with self.subTest(encoding=encoding),self.assertRaises(ValueError):
                prepare(raw,encoding=encoding,n_samples=4)

    def test_invalid_contract(self):
        raw=torch.zeros((2,1),dtype=torch.uint8,device='cuda')
        for kwargs in ({'n_samples':5},{'n_samples':0},{'n_samples':True},{'n_samples':4,'scale':2},{'n_samples':4,'missing_value':3}):
            with self.assertRaises(ValueError):prepare(raw,encoding='pgen_2bit',**kwargs)
        with self.assertRaises(ValueError):prepare(raw,n_samples=4)
