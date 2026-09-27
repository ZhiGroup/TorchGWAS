import json
import threading
from unittest.mock import patch

import numpy as np
import pytest

from torchgwas.api import run_linear_gwas
from torchgwas.bed import PlinkBedGenotype
from torchgwas.sumstats import open_binary_sumstats,open_binary_df,read_manifest
from torchgwas.sumstats_tiled import write_trait_tiled_sumstats,TiledSumstatsArray
from test_statistics import _write_bed
from torchgwas.variant_source import store_variant_ids


def test_cli_full_output_tiling_runs_and_preserves_trait_tail(tmp_path):
    from torchgwas.cli import _build_parser,_run_linear
    rng=np.random.default_rng(297)
    bed=_write_bed(tmp_path/'input',rng.integers(0,3,size=(41,9)).astype(float))
    y=rng.normal(size=(41,5)).astype(np.float32)
    np.save(tmp_path/'y.npy',y)
    args=_build_parser().parse_args(['linear','--genotype',str(bed),
        '--genotype-format','plink','--phenotype',str(tmp_path/'y.npy'),
        '--device','cpu','--trait-devices','cpu','--trait-block','2',
        '--chunk-size','4','--reader-workers','2','--prefetch-chunks','2',
        '--output-dir',str(tmp_path/'cli'),'--sumstats-format','binary'])
    assert _run_linear(args)==0
    _,actual,manifest=open_binary_sumstats(tmp_path/'cli'/'sumstats')
    run_linear_gwas(str(bed),y,genotype_format='plink',device='cpu',chunk_size=4,
        reader_workers=2,prefetch_chunks=2,output_dir=tmp_path/'reference')
    _,expected,_=open_binary_sumstats(tmp_path/'reference'/'sumstats')
    np.testing.assert_allclose(np.asarray(actual),expected,rtol=3e-5,atol=3e-6)
    assert [t['trait_range'] for t in manifest['tiles']]==[[0,2],[2,4],[4,5]]


@pytest.mark.parametrize('store_beta',[True,False])
def test_tiled_store_slices_preserve_coordinates_and_borrowed_lifetimes(tmp_path,store_beta):
    m,k=11,7
    beta=np.arange(m*k,dtype=np.float32).reshape(m,k)
    t=beta/3
    df=np.arange(m,dtype=np.float32)[:,None]+23
    assignments=[]
    def scan(offset,width,device,workers):
        assignments.append((offset,width,device,workers))
        bbuf=np.empty((4,width),dtype=np.float32)
        tbuf=bbuf.copy();dbuf=np.empty((4,1),dtype=np.float32)
        for start in range(0,m,4):
            end=min(start+4,m);size=end-start
            bbuf[:size]=beta[start:end,offset:offset+width]
            tbuf[:size]=t[start:end,offset:offset+width]
            dbuf[:size]=df[start:end]
            yield start,end,bbuf[:size],tbuf[:size],None,dbuf[:size]
            bbuf.fill(-999);tbuf.fill(-999);dbuf.fill(-999)
    summary=write_trait_tiled_sumstats(tmp_path,n_variants=m,trait_names=list('abcdefg'),n_samples=40,
        df=37,trait_block=3,devices=['cuda:1','cuda:2'],reader_workers=5,scan_factory=scan,
        block_bytes=32,queue_depth=1,fsync=True,store_beta=store_beta)
    stored_b,stored_t,manifest=open_binary_sumstats(tmp_path)
    assert isinstance(stored_t,TiledSumstatsArray)
    if store_beta:
        np.testing.assert_array_equal(np.asarray(stored_b),beta)
    else:
        assert stored_b is None
    for key in [(slice(None),slice(None)),(slice(2,9,2),slice(1,7,2)),(-1,-2),
                (slice(None,None,-1),slice(None,None,-2)),(3,slice(None)),
                (slice(0,0),slice(None)),(slice(None),2)]:
        np.testing.assert_array_equal(stored_t[key],t[key])
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path)),np.broadcast_to(df,t.shape))
    assert sorted(assignments)==[(0,3,'cuda:1',3),(3,3,'cuda:2',2),(6,1,'cuda:1',3)]
    assert summary['cells']==m*k
    assert summary['payload_bytes']==m*k*(8 if store_beta else 4)+m*4*3
    assert manifest['genotype_passes']==3


def test_failed_tile_joins_workers_and_does_not_publish_completion(tmp_path):
    closed=[]
    peer_started=threading.Event()
    def scan(offset,width,device,workers):
        try:
            if offset:
                peer_started.set()
                raise RuntimeError('intentional tile failure')
            assert peer_started.wait(3)
            for row in range(100):
                yield row,row+1,np.zeros((1,width)),np.zeros((1,width)),None,np.full((1,1),98)
        finally:
            closed.append(offset)
    with pytest.raises(RuntimeError):
        write_trait_tiled_sumstats(tmp_path,n_variants=100,trait_names=list('abcd'),n_samples=100,
            df=98,trait_block=2,devices=['cuda:1','cuda:2'],reader_workers=4,scan_factory=scan,
            block_bytes=32,queue_depth=1)
    assert sorted(closed)==[0,2]
    assert not (tmp_path/'manifest.json').exists()
    assert not any(t.name.startswith(('torchgwas-fulltile','torchgwas-write-')) for t in threading.enumerate())


@pytest.mark.parametrize('devices',[['cpu'],['cuda:1'],['cuda:1','cuda:2']])
@pytest.mark.parametrize('store_beta',[True,False])
def test_full_api_tiling_matches_untiled_with_genotype_missingness_and_tails(tmp_path,devices,store_beta):
    import torch
    if devices[0].startswith('cuda') and (not torch.cuda.is_available() or torch.cuda.device_count()<3):
        pytest.skip('three visible CUDA devices required')
    rng=np.random.default_rng(9351)
    n,m,k=97,19,8
    genotype=rng.integers(0,3,size=(n,m)).astype(float)
    if devices[0]!='cpu':
        for column in range(m):
            genotype[:column,column]=np.nan
    y=rng.normal(size=(n,k));y[:,2]=1.
    cov=rng.normal(size=(n,2));cov[:,1]=cov[:,0]*2
    bed=_write_bed(tmp_path/'input',genotype)
    pheno_path=tmp_path/'phenotype.npy';np.save(pheno_path,y)
    kwargs=dict(device=devices[0],compute_dtype='float32',chunk_size=4,reader_workers=4,
                prefetch_chunks=2,variant_range=(1,18),sumstats_fields='beta+t' if store_beta else 't',
                sumstats_block_bytes=64,sumstats_queue_depth=1,sumstats_fsync=True)
    with pytest.warns(RuntimeWarning,match='linearly dependent'):
        run_linear_gwas(PlinkBedGenotype(bed),pheno_path,cov,output_dir=tmp_path/'whole',**kwargs)
    with pytest.warns(RuntimeWarning,match='linearly dependent'):
        run_linear_gwas(PlinkBedGenotype(bed),pheno_path,cov,output_dir=tmp_path/'tiled',
                        trait_block=3,trait_devices=devices,**kwargs)
    whole_b,whole_t,whole_meta=open_binary_sumstats(tmp_path/'whole'/'sumstats')
    tile_b,tile_t,tile_meta=open_binary_sumstats(tmp_path/'tiled'/'sumstats')
    assert tile_meta['traits']==whole_meta['traits']==[f'trait_{i}' for i in range(k) if i!=2]
    assert tile_meta['shape']==[17,7]
    for actual,expected in [(tile_t,whole_t),(tile_b,whole_b)] if store_beta else [(tile_t,whole_t)]:
        np.testing.assert_allclose(np.asarray(actual),expected,rtol=3e-5,atol=3e-6,equal_nan=True)
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path/'tiled'/'sumstats')),
                                 np.broadcast_to(open_binary_df(tmp_path/'whole'/'sumstats'),tile_t.shape))
    np.testing.assert_array_equal(store_variant_ids(tmp_path/'tiled'/'sumstats'),store_variant_ids(tmp_path/'whole'/'sumstats'))


def test_blocked_input_qc_never_materializes_whole_filtered_panel(tmp_path):
    from torchgwas.preprocess import prepare_inputs_for_prep,PhenotypeColumnView,_phenotype_column_mask
    n,k=31,13
    y=np.arange(n*k,dtype=np.float64).reshape(n,k);y[:,5]=1.
    path=tmp_path/'large.npy';np.save(path,y)
    mapped=np.load(path,mmap_mode='r')
    with patch('torchgwas.preprocess._phenotype_column_mask',wraps=_phenotype_column_mask) as check:
        selected,_,qc=prepare_inputs_for_prep(np.ones((n,2)),mapped,validate_genotype=False,
                                             dtype=np.float32,phenotype_block_size=4)
    assert isinstance(selected,PhenotypeColumnView)
    assert np.shares_memory(selected.values,mapped)
    assert max(call.args[0].shape[1] for call in check.call_args_list)<=4
    np.testing.assert_array_equal(selected[:,3:7],y[:,[3,4,6,7]])
    assert qc['phenotype_kept_column_indices']==[i for i in range(k) if i!=5]


def test_reader_budget_rejected_before_any_output_mutation(tmp_path):
    sentinel=tmp_path/'manifest.json';sentinel.write_text('unchanged')
    with pytest.raises(ValueError,match='at least one reader'):
        write_trait_tiled_sumstats(tmp_path,n_variants=2,trait_names=['a','b'],n_samples=40,df=38,
            trait_block=1,devices=['cuda:1','cuda:2'],reader_workers=1,scan_factory=None)
    assert sentinel.read_text()=='unchanged'


def test_source_view_does_not_own_or_close_shared_file_handles(tmp_path):
    import gc,os
    from torchgwas.sumstats_tiled import ScanSourceView
    bed=_write_bed(tmp_path/'source',np.arange(31*7).reshape(31,7)%3)
    source=PlinkBedGenotype(bed)
    source._last_scan_exclusion_counts={'missing':0}
    view=ScanSourceView(source)
    view._last_scan_exclusion_counts={'missing':2}
    assert source._last_scan_exclusion_counts=={'missing':0}
    del view;gc.collect()
    assert os.fstat(source._fd).st_size>0
    del source
    gc.collect()


def test_native_pgen_full_tiles_preserve_variant_df_and_source_lifetime(tmp_path,monkeypatch):
    import torch
    from test_pgen_native_reader import write_pgen
    from torchgwas.io import load_genotype
    if not torch.cuda.is_available() or torch.cuda.device_count()<3:
        pytest.skip('three visible CUDA devices required')
    for key,value in dict(TORCHGWAS_NATIVE_STATS='0',TORCHGWAS_PGEN_PACKED='0',TORCHGWAS_PGEN_BACKEND='native').items():
        monkeypatch.setenv(key,value)
    n,m,k=129,23,13
    rng=np.random.default_rng(9379)
    calls=rng.integers(0,3,size=(m,n),dtype=np.uint8)
    for row in range(m):calls[row,:row]=3
    prefix=tmp_path/'native'
    write_pgen(prefix.with_suffix('.pgen'),calls)
    prefix.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    prefix.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    source,*_=load_genotype(prefix,genotype_format='pgen',pgen_mode='hardcall',reader_workers=4)
    y=rng.normal(size=(n,k));c=rng.normal(size=(n,2))
    common=dict(device='cuda:1',chunk_size=4,reader_workers=4,prefetch_chunks=2,
                compute_dtype='float32',sumstats_block_bytes=64,sumstats_queue_depth=1)
    run_linear_gwas(source,y,c,output_dir=tmp_path/'whole',**common)
    run_linear_gwas(source,y,c,output_dir=tmp_path/'tiles',trait_block=5,trait_devices=['cuda:1','cuda:2'],**common)
    # Reuse the very same source after every tile has closed its scan.
    run_linear_gwas(source,y,c,output_dir=tmp_path/'after',**common)
    ref_b,ref_t,_=open_binary_sumstats(tmp_path/'whole'/'sumstats')
    for name in ('tiles','after'):
        b,t,_=open_binary_sumstats(tmp_path/name/'sumstats')
        np.testing.assert_allclose(np.asarray(b),ref_b,rtol=3e-5,atol=3e-6)
        np.testing.assert_allclose(np.asarray(t),ref_t,rtol=3e-5,atol=3e-6)
    expected=np.sum(calls!=3,axis=1)[:,None]-4
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path/'tiles'/'sumstats')),
                                 np.broadcast_to(expected,(m,k)))
