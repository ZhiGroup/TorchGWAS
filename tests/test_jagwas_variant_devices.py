"""Joint statistics must survive variant partitioning and independent GPU state."""
import threading
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pytest
import torch
from torchgwas.api import run_linear_gwas
from torchgwas.linear import linear_scan_multigpu, linear_scan_streaming_chunks
from torchgwas.reduce import JagwasReduction
from torchgwas.sumstats_indexed import open_indexed_sumstats


def test_factory_instances_and_scan_counters_are_private():
    original = SimpleNamespace(shape=(32, 24), decode_workers=17,
        _last_scan_exclusion_counts={'invariant': 999})
    seen = []
    registered = []
    def scan(source, y, c, *, variant_range, reduction, reader_workers, **kwargs):
        assert source is not original
        assert source._last_scan_exclusion_counts == {}
        seen.append((source, reduction, variant_range, reader_workers))
        def generate():
            source._last_scan_exclusion_counts = {'invariant': 1}
            source._last_scan_profile = {'tag': variant_range[0]}
            yield variant_range[0], variant_range[1], None, np.zeros((variant_range[1]-variant_range[0],1)), None, None
        return generate(), c
    with patch('torchgwas.linear.residualize_and_standardize',
            side_effect=lambda y,c,**kw: (y,c,np.full(y.shape[1],len(y)))) as prep, \
         patch('torchgwas.linear.linear_scan_streaming_chunks', side_effect=scan):
        chunks, _ = linear_scan_multigpu(original, np.ones((32,3)), devices=['cpu:0','cpu:1'],
            chunk_size=4, variant_range=(3,19), reader_workers=5, ordered=False,
            reduction_factory=JagwasReduction, shared_queue_depth=1,
            result_queue_registration=lambda q,bounds,devices,finished:
                registered.append((q.maxsize,list(bounds),list(devices),finished)))
        result = list(chunks)
    assert prep.call_count == 1
    assert len(registered)==1 and registered[0][:3]==(1,[(3,11),(11,19)],['cpu:0','cpu:1'])
    assert len({id(row[0]) for row in seen}) == len({id(row[1]) for row in seen}) == 2
    assert sorted((row[2],row[3]) for row in seen) == [((3,11),3),((11,19),2)]
    assert sorted((row[0],row[1]) for row in result) == [(3,11),(11,19)]
    assert original.decode_workers == 17
    assert original._last_scan_exclusion_counts == {'invariant':2}
    assert original._last_scan_profile['result_queue_capacity'] == 1
    assert [row['profile']['tag'] for row in original._last_scan_profile['shards']] == [3,11]


@pytest.mark.parametrize('device_count,chunk',[(1,5),(2,32)])
def test_single_active_device_enforces_reader_budget_on_real_native_scan(tmp_path,monkeypatch,device_count,chunk):
    if not torch.cuda.is_available():pytest.skip('CUDA device required')
    from torchgwas.pgen import PgenGenotype
    import torchgwas.native_scan as native
    for key,value in dict(TORCHGWAS_PGEN_BACKEND='native',TORCHGWAS_PGEN_PACKED='0',TORCHGWAS_NATIVE_STATS='0').items():
        monkeypatch.setenv(key,value)
    path,calls,y,c=fixture(tmp_path,'pgen',missing=False)
    source=PgenGenotype(path,mode='hardcall',reader_workers=24)
    observed=[];loader=native.PinnedDosageLoader
    def load(*args,**kwargs):
        value=loader(*args,**kwargs)
        observed.append((value.decode_workers_requested,value.decode_workers_effective))
        return value
    monkeypatch.setattr(native,'PinnedDosageLoader',load)
    chunks,_=linear_scan_multigpu(source,y[:,[0,1,3,4,5,6]],c,
        devices=['cuda:0','cuda:1'][:device_count],chunk_size=chunk,
        reader_workers=1,prefetch_chunks=4,compute_dtype='float32',
        compute_p_values=False,reduction_factory=JagwasReduction)
    actual={}
    for first,last,_,values,*_ in chunks:
        for index,value in zip(range(first,last),np.asarray(values).reshape(-1)):
            assert index not in actual
            actual[index]=float(value)
    expected=reference(calls,y,c,(0,23))
    assert actual.keys()==expected.keys()
    np.testing.assert_allclose(list(actual.values()),[expected[i] for i in actual],rtol=3e-4,atol=3e-4)
    assert observed==[(1,1)] and source.decode_workers==24


@pytest.mark.parametrize('kind',['shared','duplicate','not_callable'])
def test_unsafe_factory_rejected_before_preprocessing(kind):
    instance = JagwasReduction()
    options = {'reduction':instance} if kind == 'shared' else {
        'reduction_factory': (lambda:instance) if kind == 'duplicate' else 3}
    with patch('torchgwas.linear.residualize_and_standardize',side_effect=AssertionError('must refuse first')):
        with pytest.raises(ValueError,match='reduction'):
            linear_scan_multigpu(SimpleNamespace(shape=(32,12)),np.ones((32,3)),
                devices=['cpu:0','cpu:1'],chunk_size=4,**options)


def test_invalid_joint_shape_is_refused_before_factor_allocation():
    n, k = 8, 8
    y = np.zeros((n,k),np.float32)
    reduction = JagwasReduction()
    with patch.object(reduction,'prepare',side_effect=AssertionError('factor allocated too early')):
        with pytest.raises(ValueError,match='residual phenotype rank'):
            linear_scan_streaming_chunks(SimpleNamespace(shape=(n,4)),y,None,device='cpu',
                already_processed=True,observed_counts=np.full(k,n),reduction=reduction)


def test_missing_phenotype_values_reach_the_joint_factor():
    # The mean-imputed panel's t with the common df matches R, its Gram.
    n, k = 12, 3
    y = np.zeros((n,k),np.float32)
    observed = np.full(k,n); observed[0] -= 1
    reduction = JagwasReduction()
    with patch.object(reduction,'prepare',side_effect=AssertionError('prepared')):
        with pytest.raises(AssertionError,match='prepared'):
            linear_scan_streaming_chunks(SimpleNamespace(shape=(n,4)),y,None,device='cpu',
                already_processed=True,observed_counts=observed,reduction=reduction)


def fixture(tmp_path, fmt, *, missing=True):
    from test_statistics import _write_bed
    from test_pgen_native_reader import write_pgen
    rng=np.random.default_rng(9219121); n,m,k=129,23,7
    calls=rng.integers(0,3,(n,m)).astype(np.float64)
    if missing:
        calls[:,7]=1.
        for column in range(m): calls[:column//3,column]=np.nan
        calls[:,13]=np.nan
    y=rng.normal(size=(n,k)).astype(np.float32); y[:,2]=1.
    y[:,0]+=.25*np.nan_to_num(calls[:,2],nan=1.)
    c=rng.normal(size=(n,2)).astype(np.float32)
    if fmt=='plink': path=_write_bed(tmp_path/'input',calls)
    else:
        path=tmp_path/'input.pgen'
        write_pgen(path,np.where(np.isnan(calls.T),3,calls.T).astype(np.uint8))
        path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
        path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    return path,calls,y,c


def rows(path):
    manifest,parts=open_indexed_sumstats(path/'sumstats')
    values={}
    for part in parts:
        assert set(part)=={'variant_index','chi2'}
        for index,value in zip(part['variant_index'],part['chi2']):
            assert int(index) not in values
            values[int(index)]=float(value)
    return manifest,values


def reference(calls,y,c,span):
    # Independent FP64 OLS and correlation quadratic form of the score
    # z = t / sqrt(1 + t^2 / df) (= sqrt(df) r); no production kernels.
    y=np.asarray(y[:,[i for i in range(y.shape[1]) if i!=2]],np.float64)
    x=np.column_stack([np.ones(len(y)),c.astype(np.float64)])
    yr=y-x@np.linalg.lstsq(x,y,rcond=None)[0]
    correlation=np.corrcoef(yr.T)
    rank=np.linalg.matrix_rank(x); result={}
    for variant in range(*span):
        g=calls[:,variant].copy();valid=np.isfinite(g)
        if not valid.any() or np.ptp(g[valid])==0:continue
        g[~valid]=g[valid].mean(); gr=g-x@np.linalg.lstsq(x,g,rcond=None)[0]
        ss=gr@gr; beta=(gr@yr)/ss
        residual=yr-gr[:,None]*beta
        se=np.sqrt(np.sum(residual*residual,axis=0)/(valid.sum()-rank-1)/ss)
        t=beta/se
        t=t/np.sqrt(1+t*t/(valid.sum()-rank-1))
        result[variant-span[0]]=float(t@np.linalg.solve(correlation,t))
    return result


@pytest.mark.parametrize('fmt',['plink','pgen'])
@pytest.mark.parametrize('cuda',[False,True])
@pytest.mark.parametrize('chunk',[4,9])
def test_joint_api_matches_serial_and_fp64_ols(tmp_path,monkeypatch,fmt,cuda,chunk):
    if cuda and (not torch.cuda.is_available() or torch.cuda.device_count()<2):pytest.skip('two CUDA devices required')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS','0'); monkeypatch.setenv('TORCHGWAS_PGEN_PACKED','0')
    if fmt=='pgen':monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    devices=['cuda:0','cuda:1'] if cuda else ['cpu:0','cpu:1']
    path,calls,y,c=fixture(tmp_path,fmt,missing=cuda); span=(2,22)
    options=dict(genotype_format=fmt,compute_dtype='float32',chunk_size=chunk,reader_workers=3,
        prefetch_chunks=2,variant_range=span,reduce='jagwas',sumstats_queue_depth=1)
    serial=run_linear_gwas(path,y,c,device=devices[0],output_dir=tmp_path/'serial',**options)
    joint=JagwasReduction; factors=[]; original_prepare=joint.prepare
    def prepare(self,phenotype,device=None):
        result=original_prepare(self,phenotype,device=device)
        assert phenotype.shape[1] == 6
        assert tuple(self._inverse_cholesky.shape) == (6, 6)
        factors.append((id(self),str(self._inverse_cholesky.device)))
        return result
    with patch.object(joint,'prepare',prepare):
        parallel=run_linear_gwas(path,y,c,variant_devices=devices,output_dir=tmp_path/'parallel',**options)
    expected_manifest,expected=rows(tmp_path/'serial'); manifest,actual=rows(tmp_path/'parallel')
    assert set(actual)==set(expected)==set(reference(calls,y,c,span))
    ordered=sorted(actual)
    np.testing.assert_allclose([actual[i] for i in ordered],[expected[i] for i in ordered],rtol=6e-5,atol=2e-5)
    truth=reference(calls,y,c,span)
    np.testing.assert_allclose([actual[i] for i in ordered],[truth[i] for i in ordered],rtol=3e-4,atol=3e-5)
    assert manifest['df']==expected_manifest['df']==6
    assert manifest['shape']==[20,6]
    assert len(factors)==len({r[0] for r in factors})==2
    if cuda:assert {r[1] for r in factors}==set(devices)
    assert parallel.run_metadata['variant_devices']==devices
    assert parallel.run_metadata['sumstats_write']['execution_layout']['shared_result_queue_depth']==1
    for key in ['dropped_genotype_columns','genotype_columns_kept']:
        assert parallel.qc_summary[key]==serial.qc_summary[key]
    metadata=tmp_path/'parallel'/'sumstats'/'variant_metadata.npz'
    if metadata.exists():
        with np.load(metadata) as archive:assert all(len(archive[key])==20 for key in archive.files)
    np.testing.assert_array_equal(np.load(tmp_path/'parallel'/'sumstats'/'variant_ids.npy'),
                                  np.load(tmp_path/'serial'/'sumstats'/'variant_ids.npy'))


def test_small_range_trims_idle_devices_and_keeps_global_trait_df(tmp_path,monkeypatch):
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    path,_,y,c=fixture(tmp_path,'pgen',missing=False)
    result=run_linear_gwas(path,y,c,variant_devices=['cpu:0','cpu:1'],output_dir=tmp_path/'out',
        genotype_format='pgen',compute_dtype='float32',chunk_size=8,variant_range=(2,3),
        reduce='jagwas',reader_workers=1,sumstats_queue_depth=1)
    manifest,actual=rows(tmp_path/'out')
    assert set(actual)=={0} and manifest['df']==6
    assert result.run_metadata['variant_devices']==['cpu:0']
    assert result.run_metadata['sumstats_write']['execution_layout']['shared_result_queue_depth']==0


def test_writer_failure_closes_all_joint_shards_without_manifest(tmp_path,monkeypatch):
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    path,_,y,c=fixture(tmp_path,'pgen',missing=False)
    with patch('torchgwas.sumstats_indexed.np.savez',side_effect=OSError('injected writer failure')):
        with pytest.raises(OSError,match='injected'):
            run_linear_gwas(path,y,c,variant_devices=['cpu:0','cpu:1'],output_dir=tmp_path/'out',
                genotype_format='pgen',compute_dtype='float32',chunk_size=4,
                reduce='jagwas',reader_workers=2,sumstats_queue_depth=1)
    assert not (tmp_path/'out'/'sumstats'/'manifest.json').exists()
    assert not any(t.name.startswith('torchgwas-shard-') for t in threading.enumerate())


def test_cpu_missing_genotype_contract_is_not_changed(tmp_path,monkeypatch):
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    path,_,y,c=fixture(tmp_path,'pgen',missing=True)
    with pytest.raises(ValueError,match='missing/non-finite'):
        run_linear_gwas(path,y,c,variant_devices=['cpu:0','cpu:1'],output_dir=tmp_path/'out',
            genotype_format='pgen',compute_dtype='float32',reduce='jagwas',chunk_size=4)


def test_joint_factor_failure_closes_peer_and_prevents_publication(tmp_path,monkeypatch):
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    path,_,y,c=fixture(tmp_path,'pgen',missing=False)
    joint=JagwasReduction; original=joint.prepare
    def fail_second(self,phenotype,device=None):
        if torch.device(device).index==1:raise RuntimeError('injected factor failure')
        return original(self,phenotype,device=device)
    with patch.object(joint,'prepare',fail_second):
        with pytest.raises(RuntimeError,match='injected factor'):
            run_linear_gwas(path,y,c,variant_devices=['cpu:0','cpu:1'],output_dir=tmp_path/'out',
                genotype_format='pgen',compute_dtype='float32',chunk_size=4,
                reduce='jagwas',reader_workers=2,sumstats_queue_depth=1)
    assert not (tmp_path/'out'/'sumstats'/'manifest.json').exists()
    assert not any(t.name.startswith('torchgwas-shard-') for t in threading.enumerate())


def test_joint_rank_impossible_after_qc_is_refused_before_preprocessing(tmp_path,monkeypatch):
    from test_statistics import _write_bed
    rng=np.random.default_rng(18)
    path=_write_bed(tmp_path/'rank',rng.integers(0,3,(16,8)).astype(float))
    y=rng.normal(size=(16,17)).astype(np.float32)
    with patch('torchgwas.linear.residualize_and_standardize',side_effect=AssertionError('refuse before joint preprocessing')):
        with pytest.raises(ValueError,match='residual phenotype rank'):
            run_linear_gwas(path,y,variant_devices=['cpu:0','cpu:1'],output_dir=tmp_path/'out',
                genotype_format='plink',reduce='jagwas',chunk_size=2,reader_workers=2)
    assert not (tmp_path/'out').exists()


@pytest.mark.parametrize('compute_dtype',['float32','float64'])
@pytest.mark.parametrize('from_path',[False,True])
@pytest.mark.parametrize('cuda',[False,True])
def test_joint_dtype_normalization_survives_blocked_qc(tmp_path,monkeypatch,compute_dtype,from_path,cuda):
    if cuda and (not torch.cuda.is_available() or torch.cuda.device_count()<2):pytest.skip('two CUDA devices required')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS','0');monkeypatch.setenv('TORCHGWAS_PGEN_PACKED','0')
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    path,calls,y,c=fixture(tmp_path,'pgen',missing=False)
    y=y.astype(np.float64);y[:,0]+=np.arange(len(y))*1e-9
    if from_path:
        input_path=tmp_path/'phenotype.npy';np.save(input_path,y);phenotype=input_path
    else:phenotype=y
    devices=['cuda:0','cuda:1'] if cuda else ['cpu:0','cpu:1']
    joint=JagwasReduction;seen=[];original=joint.prepare
    def prepare(self,phenotype,device=None):
        seen.append(np.dtype(str(phenotype.dtype).removeprefix('torch.')) if isinstance(phenotype,torch.Tensor)
                    else np.asarray(phenotype).dtype)
        return original(self,phenotype,device=device)
    options=dict(genotype_format='pgen',compute_dtype=compute_dtype,chunk_size=9,
        reader_workers=2,prefetch_chunks=2,reduce='jagwas',sumstats_queue_depth=1)
    with patch.object(joint,'prepare',prepare):
        run_linear_gwas(path,phenotype,c,device=devices[0],output_dir=tmp_path/'serial',**options)
        run_linear_gwas(path,phenotype,c,variant_devices=devices,output_dir=tmp_path/'parallel',**options)
    assert seen==[np.dtype(compute_dtype)]*3
    _,serial=rows(tmp_path/'serial');_,parallel=rows(tmp_path/'parallel')
    assert set(serial)==set(parallel)==set(range(calls.shape[1]))
    np.testing.assert_allclose([parallel[i] for i in sorted(parallel)],
        [serial[i] for i in sorted(serial)],rtol=6e-5,atol=2e-5)
    truth=reference(calls,y,c,(0,calls.shape[1]))
    np.testing.assert_allclose([parallel[i] for i in sorted(parallel)],
        [truth[i] for i in sorted(truth)],rtol=3e-4,atol=3e-5)


def test_factor_capacity_uses_retained_traits_and_only_active_devices(tmp_path,monkeypatch):
    monkeypatch.setenv('TORCHGWAS_PGEN_BACKEND','native')
    path,_,y,c=fixture(tmp_path,'pgen',missing=False)
    with patch('torchgwas.reduction_tensor_work.require_jagwas_factor_capacity',
               side_effect=ValueError('capacity gate')) as gate, \
         patch('torchgwas.linear.residualize_and_standardize',side_effect=AssertionError('gate too late')):
        with pytest.raises(ValueError,match='capacity gate'):
            run_linear_gwas(path,y,c,variant_devices=['cpu:0','cpu:1'],output_dir=tmp_path/'out',
                genotype_format='pgen',compute_dtype='float64',chunk_size=8,variant_range=(2,3),
                reduce='jagwas',reader_workers=1)
    gate.assert_called_once_with(129,6,['cpu:0'],compute_dtype='float64',method='eigen')
    assert not (tmp_path/'out').exists()


@pytest.mark.parametrize('partition', [
    {'trait_block': 1},
    {'trait_block': 1_000_000},  # Even a nominal no-op is not a tuning axis.
    {'trait_devices': ['cuda:0']},
    {'trait_devices': ['cuda:0', 'cuda:1']},
    {'trait_block': 1, 'trait_devices': ['cuda:0', 'cuda:1']},
])
@pytest.mark.parametrize('variant_devices', [None, ['cuda:0', 'cuda:1']])
def test_joint_phenotype_partition_refused_before_inputs_or_cuda(
        tmp_path, partition, variant_devices):
    with patch('torchgwas.api.load_genotype', side_effect=AssertionError('loaded input')), \
         patch('torchgwas.api._coerce_array_or_path', side_effect=AssertionError('loaded phenotype')), \
         patch('torchgwas.api.choose_device', side_effect=AssertionError('initialized device')):
        with pytest.raises(ValueError, match='full phenotype panel and joint-test state'):
            run_linear_gwas(tmp_path/'unopened.pgen', tmp_path/'unopened.npy',
                reduce='jagwas', output_dir=tmp_path/'out',
                variant_devices=variant_devices, **partition)
    assert not (tmp_path/'out').exists()
