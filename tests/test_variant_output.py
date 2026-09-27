import copy
import json
import threading
from pathlib import Path

import numpy as np
import pytest

from torchgwas.api import run_linear_gwas
from torchgwas.sumstats import open_binary_sumstats,open_binary_df,read_manifest
from torchgwas.sumstats_sharded import write_variant_sharded_sumstats,VariantShardedArray
from torchgwas.variant_source import store_variant_ids


def synthetic(directory, *, store_beta=True,before_publish=None):
    m,k=23,7
    b=np.arange(m*k,dtype=np.float32).reshape(m,k);t=b/3;df=(np.arange(m)+31).astype(np.float32)[:,None]
    calls=[]
    def scan(first,last,device,workers):
        calls.append((first,last,device,workers))
        bb=np.empty((4,k),np.float32);tt=bb.copy();dd=np.empty((4,1),np.float32)
        for start in range(first,last,4):
            stop=min(last,start+4);size=stop-start
            bb[:size]=b[start:stop];tt[:size]=t[start:stop];dd[:size]=df[start:stop]
            yield start-first,stop-first,bb[:size],tt[:size],None,dd[:size]
            bb.fill(-999);tt.fill(-999);dd.fill(-999)
    summary=write_variant_sharded_sumstats(directory,n_variants=m,trait_names=list('abcdefg'),
        n_samples=100,df=90,chunk_size=4,devices=['cuda:1','cuda:2'],reader_workers=5,
        scan_factory=scan,queue_depth=1,block_bytes=48,fsync=True,store_beta=store_beta,before_publish=before_publish)
    assert sorted(calls)==[(0,12,'cuda:1',3),(12,23,'cuda:2',2)]
    assert summary['payload_bytes']==m*k*(8 if store_beta else 4)+m*4
    assert summary['genotype_passes']==1 and summary['df_payload_bytes']==m*4
    return b,t,df,summary


@pytest.mark.parametrize('store_beta',[True,False])
def test_variant_store_preserves_borrowed_buffers_coordinates_and_slices(tmp_path,store_beta):
    b,t,df,_=synthetic(tmp_path,store_beta=store_beta)
    actual_b,actual_t,manifest=open_binary_sumstats(tmp_path);actual_df=open_binary_df(tmp_path)
    assert isinstance(actual_t,VariantShardedArray) and actual_df.shape==(23,1)
    if store_beta:np.testing.assert_array_equal(np.asarray(actual_b),b)
    else:assert actual_b is None
    np.testing.assert_array_equal(np.asarray(actual_t),t)
    np.testing.assert_array_equal(np.asarray(actual_df),df)
    indices=[(slice(None),slice(None)),(slice(None,None,-1),slice(None,None,-1)),
        (slice(21,1,-3),slice(6,0,-2)),(slice(2,22,3),slice(1,7,2)),
        (slice(11,15),slice(None)),(slice(0,0),slice(None)),(3,slice(None)),
        (slice(None),-2),(-1,-1),(slice(None),slice(0,0)),Ellipsis]
    rng=np.random.default_rng(492)
    for _ in range(80):
        limits=rng.integers(-35,36,size=4);steps=rng.choice([-5,-2,-1,1,2,5],size=2)
        indices.append((slice(int(limits[0]),int(limits[1]),int(steps[0])),
                        slice(int(limits[2]),int(limits[3]),int(steps[1]))))
    for key in indices:np.testing.assert_array_equal(actual_t[key],t[key])
    np.testing.assert_array_equal(actual_df[::-3,0],df[::-3,0])
    for key in [True,(slice(None),True),([1,2],slice(None)),(30,0),(0,9)]:
        with pytest.raises(IndexError):actual_t[key]
    with pytest.raises(ValueError):actual_t.__array__(copy=False)


def test_shard_failure_joins_peers_and_never_publishes_root(tmp_path):
    closed=[];peer=threading.Event()
    def scan(first,last,device,workers):
        try:
            if first:
                peer.set();raise RuntimeError('intentional shard failure')
            assert peer.wait(5)
            for i in range(last-first):
                yield i,i+1,np.zeros((1,2)),np.zeros((1,2)),None,np.ones((1,1))*30
        finally:closed.append(first)
    with pytest.raises(RuntimeError):
        write_variant_sharded_sumstats(tmp_path,n_variants=16,trait_names=['a','b'],n_samples=40,
            df=38,chunk_size=4,devices=['cuda:1','cuda:2'],reader_workers=2,scan_factory=scan,queue_depth=1)
    assert sorted(closed)==[0,8]
    assert not (tmp_path/'manifest.json').exists()
    assert not any(t.name.startswith(('torchgwas-variant','torchgwas-write-')) for t in threading.enumerate())


def test_publication_failure_does_not_advertise_completed_root(tmp_path):
    def fail():raise RuntimeError('variant ID publication failed')
    with pytest.raises(RuntimeError,match='variant ID'):synthetic(tmp_path,before_publish=fail)
    assert not (tmp_path/'manifest.json').exists()


@pytest.mark.parametrize('change,match',[
    (lambda p:p['shards'][1].update(variant_range=[11,23]),'overlap'),
    (lambda p:p['shards'][1].update(variant_range=[13,23]),'gap'),
    (lambda p:p['shards'].pop(),'Incomplete'),
    (lambda p:p['shards'][0].update(directory='../outside'),'leaves'),
    (lambda p:p.update(traits=['a']),'Invalid'),
    (lambda p:p.update(byte_order='big'),'Invalid')])
def test_bad_root_coverage_or_metadata_is_refused(tmp_path,change,match):
    synthetic(tmp_path);manifest=read_manifest(tmp_path);change(manifest)
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match=match):open_binary_sumstats(tmp_path)


def test_truncated_or_foreign_child_is_refused(tmp_path):
    synthetic(tmp_path);manifest=read_manifest(tmp_path)
    child=tmp_path/manifest['shards'][0]['directory'];record=read_manifest(child)
    record['df']['array']='../outside.f32'
    (child/'manifest.json').write_text(json.dumps(record))
    with pytest.raises(ValueError,match='leaves'):open_binary_df(tmp_path)
    record['df']['array']='df.f32';(child/'manifest.json').write_text(json.dumps(record))
    (child/'tstat.f32').write_bytes(b'bad')
    with pytest.raises(ValueError,match='payload length'):open_binary_sumstats(tmp_path)


def test_reader_budget_refusal_precedes_any_output_change(tmp_path):
    old=tmp_path/'manifest.json';old.write_text('preserved')
    with pytest.raises(ValueError,match='at least one reader'):
        write_variant_sharded_sumstats(tmp_path,n_variants=16,trait_names=['a'],n_samples=40,
            df=38,chunk_size=4,devices=['cuda:1','cuda:2'],reader_workers=1,scan_factory=None)
    assert old.read_text()=='preserved'


@pytest.mark.parametrize('fmt',['plink','pgen'])
@pytest.mark.parametrize('devices',[['cpu'],['cuda:1'],['cuda:1','cuda:2']])
@pytest.mark.parametrize('fields',['beta+t','t'])
def test_full_api_variant_shards_match_serial_with_missing_genotypes_and_ranges(tmp_path,fmt,devices,fields,monkeypatch):
    import torch
    if devices[0].startswith('cuda') and (not torch.cuda.is_available() or torch.cuda.device_count()<3):
        pytest.skip('three visible CUDA devices required')
    from test_statistics import _write_bed
    from test_pgen_native_reader import write_pgen
    rng=np.random.default_rng(2901);n,m,k=97,19,7
    calls=rng.integers(0,3,size=(n,m)).astype(float)
    if devices[0]!='cpu':
        for column in range(m):calls[:column,column]=np.nan
        calls[:,7]=1.;calls[:,13]=np.nan
    if fmt=='plink':path=_write_bed(tmp_path/'input',calls)
    else:
        path=tmp_path/'input.pgen';write_pgen(path,np.where(np.isnan(calls.T),3,calls.T).astype(np.uint8))
        path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
        path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
        monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native');monkeypatch.setenv('TORCHGWAS_PGEN_PACKED','0')
    y=rng.normal(size=(n,k)).astype(np.float32);y[:,2]=1
    cov=rng.normal(size=(n,2)).astype(np.float32);np.save(tmp_path/'y.npy',y)
    options=dict(genotype_format=fmt,pgen_mode='hardcall',compute_dtype='float32',chunk_size=4,
        reader_workers=4,prefetch_chunks=2,variant_range=(1,18),sumstats_fields=fields,
        sumstats_queue_depth=1,sumstats_block_bytes=64)
    serial=run_linear_gwas(path,tmp_path/'y.npy',cov,device=devices[0],output_dir=tmp_path/'serial',**options)
    parallel=run_linear_gwas(path,tmp_path/'y.npy',cov,variant_devices=devices,output_dir=tmp_path/'sharded',**options)
    expected=open_binary_sumstats(tmp_path/'serial'/'sumstats');actual=open_binary_sumstats(tmp_path/'sharded'/'sumstats')
    for index in ([0,1] if fields=='beta+t' else [1]):
        np.testing.assert_allclose(np.asarray(actual[index]),np.asarray(expected[index]),rtol=3e-5,atol=3e-6,equal_nan=True)
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path/'sharded'/'sumstats')),
                                  np.asarray(open_binary_df(tmp_path/'serial'/'sumstats')))
    assert actual[2]['shape']==[17,6]
    assert actual[2]['traits']==[f'trait_{i}' for i in range(k) if i!=2]
    np.testing.assert_array_equal(store_variant_ids(tmp_path/'sharded'/'sumstats'),store_variant_ids(tmp_path/'serial'/'sumstats'))
    assert parallel.run_metadata['variant_devices']==devices
    for key in ['dropped_genotype_columns','genotype_columns_kept']:
        assert parallel.qc_summary[key]==serial.qc_summary[key]


def _indexed_pairs(directory):
    from torchgwas.sumstats_indexed import open_indexed_sumstats
    manifest,parts=open_indexed_sumstats(directory)
    rows={}
    for part in parts:
        for i,(v,j) in enumerate(zip(part['variant_index'],part['trait_index'])):
            assert (int(v),int(j)) not in rows  # no pair written twice across shards
            rows[int(v),int(j)]=(float(part['beta'][i]),float(part['t_stat'][i]),float(part['df'][i]))
    assert list(rows)==sorted(rows)  # published in (variant, trait) order
    return manifest,rows


@pytest.mark.parametrize('devices',[['cuda:0','cuda:1'],['cuda:0','cuda:1','cuda:2']])
@pytest.mark.parametrize('span',[None,(1,18)])
def test_significant_variant_shards_match_single_gpu(tmp_path,devices,span,monkeypatch):
    # Pairs are selected cell by cell, so splitting variants over GPUs must
    # write exactly the single-GPU pair set (missing calls, constant and
    # all-missing variants, a variant range, and an uneven last shard).
    import torch
    if not torch.cuda.is_available() or torch.cuda.device_count()<len(devices):
        pytest.skip(f'{len(devices)} visible CUDA devices required')
    from test_pgen_native_reader import write_pgen
    rng=np.random.default_rng(2902);n,m,k=97,23,9
    calls=rng.integers(0,3,size=(n,m)).astype(float)
    for column in range(m):calls[:column%5,column]=np.nan
    calls[:,7]=1.;calls[:,13]=np.nan
    y=rng.normal(size=(n,k)).astype(np.float32);y[:,:3]+=.8*np.nan_to_num(calls[:,:3],nan=1.)
    path=tmp_path/'input.pgen';write_pgen(path,np.where(np.isnan(calls.T),3,calls.T).astype(np.uint8))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native');monkeypatch.setenv('TORCHGWAS_PGEN_PACKED','0')
    cov=rng.normal(size=(n,2)).astype(np.float32);np.save(tmp_path/'y.npy',y)
    options=dict(pgen_mode='hardcall',compute_dtype='float32',chunk_size=4,reader_workers=6,prefetch_chunks=2,
                 variant_range=span,reduce='significant',significance_threshold=.3)
    run_linear_gwas(path,tmp_path/'y.npy',cov,device='cuda:0',output_dir=tmp_path/'serial',**options)
    sharded=run_linear_gwas(path,tmp_path/'y.npy',cov,variant_devices=devices,output_dir=tmp_path/'sharded',**options)
    expected_manifest,expected=_indexed_pairs(tmp_path/'serial'/'sumstats')
    manifest,actual=_indexed_pairs(tmp_path/'sharded'/'sumstats')
    assert expected and actual.keys()==expected.keys()
    for key,(beta,t,df) in expected.items():
        np.testing.assert_allclose(actual[key][:2],(beta,t),rtol=3e-5,atol=3e-6)
        assert actual[key][2]==df
    assert manifest['shape']==expected_manifest['shape'] and manifest['traits']==expected_manifest['traits']
    assert sharded.run_metadata['variant_devices']==devices
    assert sharded.run_metadata['sumstats_write']['execution_layout']['partition_axis']=='variant'


def test_cli_variant_devices_forwarding_and_early_conflicts(tmp_path,monkeypatch):
    import torchgwas.cli as cli
    seen=[];monkeypatch.setattr(cli,'run_linear_gwas',lambda **k:seen.append(k))
    args=cli._build_parser().parse_args(['linear','--genotype','x.pgen','--phenotype','y.npy',
        '--variant-devices','cuda:1','cuda:2','--output-dir','out'])
    assert cli._run_linear(args)==0 and seen[0]['variant_devices']==['cuda:1','cuda:2']
    # reduce='significant' is allowed: pairs are selected per cell on each shard.
    for change in [dict(trait_block=2),dict(trait_devices=['cuda:1']),
                   dict(sumstats_format='none'),dict(p_value_threshold=.01),dict(pipeline_profile={})]:
        with pytest.raises(ValueError,match='variant_devices'):
            run_linear_gwas('missing.pgen','missing.npy',variant_devices=['cuda:1'],output_dir=tmp_path/'out',**change)
    assert not (tmp_path/'out').exists()

@pytest.mark.parametrize('dtype',[np.int8,np.float32,np.float64])
def test_pgen_cpu_ranges_respect_decode_boundaries_and_global_indices(tmp_path,dtype):
    from test_pgen_native_reader import write_pgen
    from torchgwas.pgen import PgenGenotype
    rng=np.random.default_rng(499);calls=rng.integers(0,3,size=(19,32),dtype=np.uint8)
    path=tmp_path/'input.pgen';write_pgen(path,calls)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(19)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(32)))
    source=PgenGenotype(path,mode='hardcall',decode_batch_size=3,reader_workers=2,prefetch_chunks=2)
    try:
        chunks=list(source.iter_chunks(4,dtype=dtype,variant_range=(5,16)))
        assert [(a,b) for a,b,_ in chunks]==[(5,9),(9,13),(13,16)]
        np.testing.assert_array_equal(np.concatenate([x for _,_,x in chunks],axis=1),calls[5:16].T)
        assert list(source.iter_chunks(4,dtype=dtype,variant_range=(7,7)))==[]
    finally:source.close()

@pytest.mark.parametrize('devices',[['cpu'],['cuda:1'],['cuda:1','cuda:2']])
@pytest.mark.parametrize('span',[None,(1,2)])
def test_bgen_dosage_variant_shards_match_serial_and_trim_idle_devices(tmp_path,devices,span):
    import torch
    if devices[0].startswith('cuda') and (not torch.cuda.is_available() or torch.cuda.device_count()<3):
        pytest.skip('three visible CUDA devices required')
    from test_bgen_direct import write_bgen
    path=tmp_path/'input.bgen';write_bgen(path,bits=8,missing=devices[0]!='cpu')
    y=np.random.default_rng(734).normal(size=(4,3)).astype(np.float32)
    options=dict(genotype_format='bgen',compute_dtype='float32',chunk_size=1,
                 reader_workers=2,prefetch_chunks=2,variant_range=span,sumstats_queue_depth=1)
    run_linear_gwas(path,y,device=devices[0],output_dir=tmp_path/'serial',**options)
    result=run_linear_gwas(path,y,variant_devices=devices,output_dir=tmp_path/'shards',**options)
    expected=open_binary_sumstats(tmp_path/'serial'/'sumstats')
    actual=open_binary_sumstats(tmp_path/'shards'/'sumstats')
    for index in (0,1):
        np.testing.assert_allclose(np.asarray(actual[index]),np.asarray(expected[index]),rtol=3e-5,atol=3e-6,equal_nan=True)
    np.testing.assert_array_equal(np.asarray(open_binary_df(tmp_path/'shards'/'sumstats')),
                                  np.asarray(open_binary_df(tmp_path/'serial'/'sumstats')))
    active=devices[:1] if span else devices
    assert result.run_metadata['variant_devices']==active
    np.testing.assert_array_equal(store_variant_ids(tmp_path/'shards'/'sumstats'),store_variant_ids(tmp_path/'serial'/'sumstats'))


@pytest.mark.parametrize('layout',['shards','tiles'])
@pytest.mark.parametrize('devices',[['cpu'],['cuda:1','cuda:2']])
def test_full_api_partitions_keep_the_missing_phenotype_contract(tmp_path,layout,devices,monkeypatch):
    # Missing phenotype cells: tiles and shards must write what one device
    # writes -- t of the mean-imputed panel, per-trait df, no variant df.
    import torch
    if devices[0].startswith('cuda') and (not torch.cuda.is_available() or torch.cuda.device_count()<3):
        pytest.skip('three visible CUDA devices required')
    from test_pgen_native_reader import write_pgen
    rng=np.random.default_rng(3107);n,m,k=97,19,7
    calls=rng.integers(0,3,size=(n,m)).astype(np.uint8)
    path=tmp_path/'input.pgen';write_pgen(path,calls.T)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native');monkeypatch.setenv('TORCHGWAS_PGEN_PACKED','0')
    y=rng.normal(size=(n,k)).astype(np.float32);y[:9,1]=np.nan;y[40:47,4]=np.nan;y[3,6]=np.nan
    np.save(tmp_path/'y.npy',y)
    options=dict(genotype_format='pgen',pgen_mode='hardcall',compute_dtype='float32',chunk_size=4,
        reader_workers=4,prefetch_chunks=2,sumstats_queue_depth=1,sumstats_block_bytes=64)
    serial=run_linear_gwas(path,tmp_path/'y.npy',device=devices[0],output_dir=tmp_path/'serial',**options)
    layout_options=(dict(variant_devices=devices) if layout=='shards' else dict(trait_block=3,trait_devices=devices))
    run_linear_gwas(path,tmp_path/'y.npy',output_dir=tmp_path/'parts',**layout_options,**options)
    expected=open_binary_sumstats(tmp_path/'serial'/'sumstats');actual=open_binary_sumstats(tmp_path/'parts'/'sumstats')
    for index in [0,1]:
        np.testing.assert_allclose(np.asarray(actual[index]),np.asarray(expected[index]),rtol=3e-5,atol=3e-6,equal_nan=True)
    assert isinstance(expected[2]['df'],list) and actual[2]['df']==expected[2]['df']
    np.testing.assert_array_equal(open_binary_df(tmp_path/'parts'/'sumstats'),open_binary_df(tmp_path/'serial'/'sumstats'))
    assert serial.qc_summary['phenotype_missing_cells']==17
