"""All significant pairs, correct df, bounded device selection and API parity."""
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import torch
from scipy import special
from torchgwas.reduce import SignificantPairs,device_significance_critical,device_significant_pairs
from torchgwas.linear import _significant_pairs_iterator
from torchgwas.api import run_linear_gwas
from torchgwas.sumstats import open_binary_sumstats,open_binary_df
from torchgwas.sumstats_indexed import open_indexed_sumstats
from test_pgen_native_reader import write_pgen
from torchgwas.variant_source import store_variant_ids


def selected_arrays(chunks):
    chunks=list(chunks)
    result=[np.concatenate([np.asarray(c[i]) for c in chunks]) for i in range(2,7)]
    order=np.lexsort((result[1],result[0]))
    return [a[order] for a in result]


@pytest.mark.parametrize('device',['cpu','cuda:1'])
@pytest.mark.parametrize('threshold',[1.,.05,1e-14,None])
@pytest.mark.parametrize('limit',[3,17,1000])
def test_device_filter_preserves_threshold_boundary_all_pairs_and_df(device,threshold,limit):
    if device!='cpu' and torch.cuda.device_count()<2:pytest.skip('second CUDA device required')
    sig=SignificantPairs(threshold);df=np.array([1,2,3,20,77,98],np.float32)
    critical=sig.critical_abs_t(df[:,None],11)
    center=critical.astype(np.float32)
    t=np.concatenate([np.nextafter(center,-np.inf),center,np.nextafter(center,np.inf),
        -np.nextafter(center,-np.inf),-center,-np.nextafter(center,np.inf),
        np.zeros_like(center),np.full_like(center,np.nan),np.full_like(center,np.inf)],axis=1)
    beta=np.arange(t.size,dtype=np.float32).reshape(t.shape)
    status=np.array([0,0,1,0,0,0],np.uint8)
    expected=np.nonzero(np.isfinite(t)&(np.abs(t).astype(np.float64)>=critical)&(status[:,None]==0))
    table=device_significance_critical(sig,100,11,device)
    actual=selected_arrays(device_significant_pairs(*[torch.as_tensor(v,device=device) for v in [beta,t,status,df]],
        table,start=13,max_cells=limit))
    for got,want in zip(actual,[expected[0]+13,expected[1],beta[expected],t[expected],df[expected[0]]]):
        np.testing.assert_array_equal(got,want)


def test_host_filter_uses_produced_df_instead_of_nominal_and_keeps_empty_chunks():
    t=np.array([[3.,3.],[3.,np.nan]],np.float32);beta=np.ones_like(t)
    chunks=[(4,6,beta,t,None,np.array([[2.],[98.]],np.float32)),(6,7,beta[:1],t[:1],None,np.array([[2.]],np.float32))]
    values=list(_significant_pairs_iterator(chunks,SignificantPairs(.05),2,98))
    assert len(values)==2 and values[1][2].size==0
    np.testing.assert_array_equal(values[0][2],[5]);np.testing.assert_array_equal(values[0][-1],[98.])


def read_pairs(directory):
    manifest,parts=open_indexed_sumstats(directory)
    chunks=[(0,0,p['variant_index'],p['trait_index'],p.get('beta',np.full(len(p['t_stat']),np.nan)),p['t_stat'],p['df']) for p in parts]
    assert manifest['df']==dict(layout='per_part',axis='pair',field='df')
    return selected_arrays(chunks) if chunks else [np.empty(0) for _ in range(5)]


@pytest.mark.parametrize('device',['cpu','cuda:1'])
@pytest.mark.parametrize('block',[None,4])
@pytest.mark.parametrize('missing_pheno,convention',[(False,'impute'),(True,'impute'),(True,'exact')])
def test_api_selected_pairs_match_full_pgen_statistics_with_missing_calls(tmp_path,monkeypatch,device,block,missing_pheno,convention):
    if device!='cpu' and torch.cuda.device_count()<2:pytest.skip('second CUDA device required')
    rng=np.random.default_rng(813);n,m,k=97,19,11
    calls=rng.integers(0,3,size=(m,n),dtype=np.uint8)
    if device!='cpu':
        for i in range(m):calls[i,:i]=3
        calls[7]=1;calls[13]=3
    path=tmp_path/'input.pgen';write_pgen(path,calls)
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    y=rng.normal(size=(n,k)).astype(np.float32);cov=rng.normal(size=(n,2)).astype(np.float32)
    if missing_pheno:y[:5,0]=np.nan;y[:9,2]=np.nan
    for key,value in [('TORCHGWAS_PGEN_BACKEND','native'),('TORCHGWAS_PGEN_PACKED','0'),('TORCHGWAS_NATIVE_STATS','0')]:monkeypatch.setenv(key,value)
    # A convention that keeps the samples (pair df), which the reference below recomputes.
    options=dict(genotype_format='pgen',pgen_mode='hardcall',device=device,compute_dtype='float32',chunk_size=4,missing_phenotype=convention,
        reader_workers=2,prefetch_chunks=2,variant_range=(1,18),sumstats_queue_depth=1)
    run_linear_gwas(path,y,cov,output_dir=tmp_path/'full',**options)
    beta,t,_logp,_=open_binary_sumstats(tmp_path/'full/sumstats');df=np.broadcast_to(open_binary_df(tmp_path/'full/sumstats'),t.shape)
    if missing_pheno and convention=='exact':
        # Dense stores keep only each trait's df here. The selected pairs carry
        # their own complete-case df: call and phenotype observed, less 4.
        df=((calls[1:18]!=3).astype(np.int64)@np.isfinite(y).astype(np.int64)-4).astype(np.float64)
    elif missing_pheno:
        # The release's pair df: the variant's df times trait_df / df.
        df=(np.count_nonzero(calls[1:18]!=3,axis=1)[:,None]-4)*(np.isfinite(y).sum(0)[None,:]-4)/float(n-4)
    keep=np.nonzero(np.isfinite(t)&(2*special.stdtr(df,-np.abs(t))<=.2))
    expected=[keep[0],keep[1],beta[keep],t[keep],df[keep]]
    assert len(keep[0])>0
    for backend in ['host','device']:
        monkeypatch.setenv('TORCHGWAS_SIGNIFICANCE_BACKEND',backend)
        with patch('torchgwas.reduce.device_significant_pairs',wraps=device_significant_pairs) as selector:
            run_linear_gwas(path,y,cov,output_dir=tmp_path/backend,reduce='significant',significance_threshold=.2,trait_block=block,**options)
        if backend=='device' and device!='cpu' and not missing_pheno:assert selector.call_count>0
        if backend=='host' or device=='cpu':assert selector.call_count==0
        actual=read_pairs(tmp_path/backend/'sumstats')
        for index,(got,want) in enumerate(zip(actual,expected)):
            # The pair df is staged as float32: exact for 'exact' (integers), 6e-8 for 'impute'.
            if index==4 and missing_pheno:np.testing.assert_allclose(got,want,rtol=1e-14 if convention=='exact' else 1e-7,atol=0)
            elif index in (0,1,4):np.testing.assert_array_equal(got,want)
            else:np.testing.assert_allclose(got,want,rtol=3e-5,atol=3e-6)


def test_multigpu_selected_pair_queue_owns_outputs_and_keeps_global_trait_indices(tmp_path,monkeypatch):
    if torch.cuda.device_count()<3:pytest.skip('CUDA devices 1 and 2 required')
    from torchgwas.api import _trait_blocked_significant_chunks
    n,m,k=40,7,11;sig=SignificantPairs(1.)
    def scan(offset,width,device):
        dev=device or 'cuda:1';table=device_significance_critical(sig,n,k,dev)
        beta=torch.arange(m*width,dtype=torch.float32,device=dev).reshape(m,width)+offset
        t=torch.ones_like(beta)*3;status=torch.zeros(m,dtype=torch.uint8,device=dev);df=torch.full((m,),38.,device=dev)
        yield from device_significant_pairs(beta,t,status,df,table,max_cells=5)
        beta.fill_(-999);t.fill_(-999)
    actual=selected_arrays(_trait_blocked_significant_chunks(scan,sig,k,4,38,devices=['cuda:1','cuda:2']))
    assert len(actual[0])==m*k
    np.testing.assert_array_equal(actual[0],np.repeat(np.arange(m),k));np.testing.assert_array_equal(actual[1],np.tile(np.arange(k),m))
    np.testing.assert_array_equal(actual[3],np.full(m*k,3.));np.testing.assert_array_equal(actual[4],np.full(m*k,38.))

def test_device_selection_rejects_precision_outside_its_exact_threshold_contract():
    sig=SignificantPairs(.05);table=device_significance_critical(sig,40,3,'cpu')
    with pytest.raises(ValueError,match='FP32'):
        list(device_significant_pairs(torch.ones((2,3),dtype=torch.float64),torch.ones((2,3),dtype=torch.float64),
            torch.zeros(2,dtype=torch.uint8),torch.ones(2)*38,table))

def test_threshold_one_includes_exact_zero_and_rejects_invalid_df():
    sig=SignificantPairs(1.);df=np.array([1.,20.,98.,0.],np.float32)
    t=np.zeros((4,3),np.float32);beta=np.ones_like(t)
    host=selected_arrays(_significant_pairs_iterator([(0,4,beta,t,None,df[:,None])],sig,3,98))
    assert len(host[0])==9
    table=device_significance_critical(sig,100,3,'cpu')
    selected=selected_arrays(device_significant_pairs(torch.from_numpy(beta),torch.from_numpy(t),
        torch.zeros(4,dtype=torch.uint8),torch.from_numpy(df),table,max_cells=2))
    for a,b in zip(host,selected):np.testing.assert_array_equal(a,b)
    assert np.all(2*special.stdtr(host[-1],-np.abs(host[-2]))==1.)


@pytest.mark.parametrize('backend', ['host', 'device'])
@pytest.mark.parametrize('threshold', [.2, 1.])
def test_real_native_multigpu_significance_tiles_isolate_qc_and_share_setup(tmp_path, monkeypatch, backend, threshold):
    if torch.cuda.device_count() < 3:
        pytest.skip('CUDA devices 1 and 2 required')
    from torchgwas.pgen import PgenGenotype
    from torchgwas import preprocess, linear
    rng = np.random.default_rng(92261)
    n, m, k = 129, 67, 17
    calls = rng.integers(0, 3, size=(m, n), dtype=np.uint8)
    for index in range(m):
        calls[index, :index % 19] = 3
    calls[7] = 1
    calls[13] = 3
    path = tmp_path / 'input.pgen'
    write_pgen(path, calls)
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    for key, value in [('TORCHGWAS_PGEN_BACKEND', 'native'), ('TORCHGWAS_PGEN_PACKED', '0'),
                       ('TORCHGWAS_NATIVE_STATS', '0'), ('TORCHGWAS_SIGNIFICANCE_BACKEND', backend)]:
        monkeypatch.setenv(key, value)
    source = PgenGenotype(path, mode='hardcall', reader_workers=3)
    y = rng.normal(size=(n, k)).astype(np.float32)
    cov = rng.normal(size=(n, 2)).astype(np.float32)
    opts = dict(device='cuda:1', compute_dtype='float32', chunk_size=7,
                reader_workers=3, prefetch_chunks=2, variant_range=(3, 62), sumstats_queue_depth=1)
    full = run_linear_gwas(source, y, cov, output_dir=tmp_path/'full', **opts)
    beta, t, _logp, _ = open_binary_sumstats(tmp_path/'full/sumstats')
    df = np.broadcast_to(open_binary_df(tmp_path/'full/sumstats'), t.shape)
    keep = np.nonzero(np.isfinite(t) & (2 * special.stdtr(df, -np.abs(t)) <= threshold))
    expected = [keep[0], keep[1], beta[keep], t[keep], df[keep]]
    with patch('torchgwas.preprocess._covariate_basis', wraps=preprocess._covariate_basis) as basis, \
         patch('scipy.special.stdtrit', wraps=special.stdtrit) as inverse, \
         patch('torchgwas.api.linear_scan_streaming_chunks', wraps=linear.linear_scan_streaming_chunks) as scan:
        result = run_linear_gwas(source, y, cov, output_dir=tmp_path/'selected', reduce='significant',
            significance_threshold=threshold, trait_block=6, trait_devices=['cuda:1', 'cuda:2'], **opts)
    assert basis.call_count == 1
    assert scan.call_count == 3
    assert inverse.call_count == (0 if threshold == 1. else 1)
    assert all(call.kwargs['borrow_results'] for call in scan.call_args_list)
    views = [call.args[0] for call in scan.call_args_list]
    assert len({id(view) for view in views}) == 3 and all(view is not source for view in views)
    assert sorted(call.kwargs['_reader_worker_limit'] for call in scan.call_args_list) == [1, 2, 2]
    assert result.qc_summary['genotype_exclusion_counts'] == full.qc_summary['genotype_exclusion_counts'] == dict(missing=0, invariant=2)
    assert result.qc_summary['n_variants_excluded'] == 2
    actual = read_pairs(tmp_path/'selected/sumstats')
    assert len(actual[0]) > 0
    for index, (got, want) in enumerate(zip(actual, expected)):
        if index in (0, 1, 4):
            np.testing.assert_array_equal(got, want)
        else:
            np.testing.assert_allclose(got, want, rtol=4e-5, atol=3e-6)
    manifest, _ = open_indexed_sumstats(tmp_path/'selected/sumstats')
    np.testing.assert_array_equal(store_variant_ids(tmp_path/'selected/sumstats'), [f'v{i}' for i in range(3, 62)])
    profile = source._last_scan_profile
    assert profile['layout']['queue_depth'] == 1
    assert profile['layout']['readers_per_device'] == [2, 1]
    assert [record['trait_range'] for record in profile['tiles']] == [[0, 6], [6, 12], [12, 17]]
    if backend == 'device':
        assert profile['result_payload_bytes'] == 3 * (62 - 3) + 28 * len(actual[0])

@pytest.mark.parametrize('dtype', [np.float32, np.float64])
def test_significant_npy_input_and_qc_stay_trait_bounded(tmp_path, monkeypatch, dtype):
    from torchgwas import preprocess
    rng = np.random.default_rng(7126)
    n, m, k = 97, 19, 11
    calls = rng.integers(0, 3, size=(m, n), dtype=np.uint8)
    path = tmp_path/'input.pgen'
    write_pgen(path, calls)
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(n)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    y = rng.normal(size=(n, k)).astype(dtype)
    y[:, 2] = 1.  # Keep a lazy filtered column view after QC.
    y[:5, 4] = np.nan
    panel = tmp_path/'phenotype.npy'
    np.save(panel, y)
    loaded = []
    original_load = np.load
    def load(filename, *args, **kwargs):
        result = original_load(filename, *args, **kwargs)
        if Path(filename) == panel:
            assert kwargs['mmap_mode'] == 'r'
            loaded.append(result)
        return result
    prepared = []
    original_prepare = preprocess.prepare_inputs_for_prep
    def prepare(*args, **kwargs):
        result = original_prepare(*args, **kwargs)
        prepared.append(result[0])
        return result
    with patch('numpy.load', side_effect=load), \
         patch('torchgwas.preprocess._phenotype_column_mask', wraps=preprocess._phenotype_column_mask) as qc, \
         patch('torchgwas.api.prepare_inputs_for_prep', side_effect=prepare):
        # 'impute' keeps the memory-mapped panel as loaded; dropping subjects copies the kept rows.
        result = run_linear_gwas(path, panel, device='cpu', genotype_format='pgen', pgen_mode='hardcall',
            compute_dtype='float32', trait_block=4, chunk_size=4, reader_workers=2, missing_phenotype='impute',
            reduce='significant', significance_threshold=1., output_dir=tmp_path/'selected')
    assert len(loaded) == 1 and isinstance(loaded[0], np.memmap)
    assert max(call.args[0].shape[1] for call in qc.call_args_list) <= 4
    assert np.shares_memory(prepared[0].values, loaded[0])
    assert result.qc_summary['dropped_phenotype_columns'] == 1
    pairs = read_pairs(tmp_path/'selected/sumstats')
    assert len(pairs[0]) == m * (k - 1)
    assert pairs[2].dtype == np.float32 and pairs[3].dtype == np.float32

def test_integer_df_cache_is_exact_readonly_and_threshold_scoped():
    from scipy import special
    sig = SignificantPairs(.05)
    dfs = np.array([-4., 0., 1., 2., 3., 20., 97.])
    expected = np.abs(special.stdtrit(dfs, .025))
    with patch('scipy.special.stdtrit', wraps=special.stdtrit) as inverse:
        table = sig.prepare_integer_df(100, 11)
        assert inverse.call_count == 1 and not table.flags.writeable
        assert sig.prepare_integer_df(100, 11) is table
        for _ in range(3):
            np.testing.assert_array_equal(sig.critical_abs_t(dfs, 11), expected)
        assert inverse.call_count == 1
        fractional = np.array([2.5, 13.7])
        actual = sig.critical_abs_t(fractional, 11)
        assert inverse.call_count == 2
    np.testing.assert_array_equal(actual, np.abs(special.stdtrit(fractional, .025)))
    sig.threshold = .2
    np.testing.assert_array_equal(sig.critical_abs_t(dfs, 11), np.abs(special.stdtrit(dfs, .1)))
    assert sig.prepare_integer_df(100, 11) is not table
    default = SignificantPairs()
    default.prepare_integer_df(100, 11)
    np.testing.assert_array_equal(default.critical_abs_t(dfs, 17), np.abs(special.stdtrit(dfs, 5e-8/17/2)))