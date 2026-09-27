"""T-only transport, ownership and independently bound component accounting."""
import copy
from unittest.mock import patch

import numpy as np
import pytest
import torch

from torchgwas.linear import linear_scan_streaming_chunks, linear_scan_multigpu, _significant_pairs_iterator
from torchgwas.native_control_work import native_control_work
from torchgwas.owned_result_work import owned_result_work
from torchgwas.pinned_work import pinned_scan_work
from torchgwas.tensor_memory import eager_scan_memory
from torchgwas.mechanistic_torch import torch_scan_work
from torchgwas.result_service import result_finish_prices
from torchgwas.significant_host_work import host_significant_selection_work
from torchgwas.reduce import SignificantPairs
from torchgwas.sumstats_indexed import write_indexed_sumstats, open_indexed_sumstats


class Source:
    supports_fused_qc=True
    decode_workers=1

    def __init__(self, values):
        self.values=values
        self.shape=values.shape
        self.sample_ids=np.arange(self.shape[0]).astype(str)
        self.marker_ids=np.arange(self.shape[1]).astype(str)

    @property
    def genotype(self):return self

    def iter_chunks(self,chunk_size,dtype=np.float32,variant_range=None,**kwargs):
        lo,hi=variant_range or (0,self.shape[1])
        for start in range(lo,hi,chunk_size):
            end=min(start+chunk_size,hi)
            yield start,end,self.values[:,start:end].astype(dtype)


def data():
    rng=np.random.default_rng(9922)
    g=rng.integers(0,3,(129,47)).astype(np.float32)
    g[:,4]=1;g[::7,5]=np.nan;g[:,6]=np.nan
    return g,rng.normal(size=(129,7)),rng.normal(size=(129,2))


@pytest.mark.parametrize('device',['cpu','cuda:1'])
@pytest.mark.parametrize('borrow',[False,True])
@pytest.mark.parametrize('missing',[False,True])
def test_transport_retains_statistics_df_and_ownership(device,borrow,missing,monkeypatch):
    if device.startswith('cuda') and torch.cuda.device_count()<2:pytest.skip('CUDA 1 required')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS','0')
    g,y,c=data()
    if missing:y[::13,0]=np.nan
    def scan(flag,owned=False):
        source=Source(g)
        rows,_=linear_scan_streaming_chunks(source,y,c,chunk_size=7,device=device,
            prefetch_chunks=2,return_beta=flag,return_df=True,borrow_results=False if owned else borrow)
        return source,rows
    _,rows=scan(True,owned=True);reference=list(rows)
    source,rows=scan(False);retained=[]
    for actual,expected in zip(rows,reference,strict=True):
        assert actual[:2]==expected[:2] and actual[2] is None
        for a,e in zip(actual[3:],expected[3:]):np.testing.assert_allclose(a,e,rtol=1e-6,atol=1e-7,equal_nan=True)
        if device.startswith('cuda'):assert actual[3].flags.owndata != borrow
        retained.append(actual)
    if not borrow:
        for actual,expected in zip(retained,reference):np.testing.assert_array_equal(actual[3],expected[3])
    if device.startswith('cuda'):
        assert source._last_scan_profile['result_payload_bytes']==47*(4*7+5)
        assert source._last_scan_profile['return_beta'] is False


@pytest.mark.parametrize('borrow',[False,True])
def test_two_gpu_t_only_covers_every_variant_once(borrow,monkeypatch):
    if torch.cuda.device_count()<3:pytest.skip('CUDA 1 and 2 required')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS','0')
    g,y,c=data()
    baseline,_=linear_scan_multigpu(Source(g),y,c,chunk_size=7,devices=['cuda:1','cuda:2'],
        reader_workers=2,prefetch_chunks=2,return_df=True,ordered=False)
    reference={row[0]:row for row in baseline};coverage=np.zeros(47,int)
    rows,_=linear_scan_multigpu(Source(g),y,c,chunk_size=7,devices=['cuda:1','cuda:2'],
        reader_workers=2,prefetch_chunks=2,return_df=True,ordered=False,return_beta=False,borrow_results=borrow)
    for row in rows:
        assert row[2] is None
        coverage[row[0]:row[1]]+=1
        for a,e in zip(row[3:],reference[row[0]][3:]):np.testing.assert_allclose(a,e,rtol=1e-6,atol=1e-7,equal_nan=True)
    assert (coverage==1).all()


@pytest.mark.parametrize('threshold',[1.,.05,1e-50])
def test_host_selected_t_only_preserves_coordinates_df_and_archives(tmp_path,threshold):
    t=np.array([[0.,3.,np.nan],[-4.,1.,20.]],np.float32)
    beta=np.ones_like(t);df=np.array([[50.],[30.]],np.float32)
    significance=SignificantPairs(threshold=threshold)
    def select(value):return list(_significant_pairs_iterator([(3,5,value,t,None,df)],significance,3,50))[0]
    a,b=select(None),select(beta)
    assert a[4] is None
    for i in (2,3,5,6):np.testing.assert_array_equal(a[i],b[i])
    count,_=write_indexed_sumstats(tmp_path,[str(i) for i in range(5)],['a','b','c'],54,
        [a],kind='significant',df=50,store_beta=False)
    manifest,parts=open_indexed_sumstats(tmp_path)
    assert count==a[2].size
    for part in parts:
        assert 'beta' not in part
        np.testing.assert_array_equal(part['t_stat'],b[5])


def test_t_only_memory_and_operation_ledgers():
    full=pinned_scan_work(129,7,512,4);small=pinned_scan_work(129,7,512,4,return_beta=False)
    assert full['requested_bytes']-small['requested_bytes']==4*7*512*4
    assert full['allocation_count']-small['allocation_count']==4
    work=owned_result_work(7,512,return_beta=False)
    assert work['copy_bytes']==7*(4*512+5) and work['allocation_calls']==3
    assert 'beta' not in work['array_bytes']
    assert native_control_work(return_beta=False)['result_submit']['copy_d2h']==3
    full=eager_scan_memory(129,7,512,2,4);small=eager_scan_memory(129,7,512,2,4,return_beta=False)
    beta=((4*7*512+511)//512)*512
    assert full['retained_output_bytes']-small['retained_output_bytes']==3*beta
    assert small['previous_beta_bytes']==beta
    assert full['tensor_storage_budget']-small['tensor_storage_budget']==2*beta
    selected=host_significant_selection_work(7,512,19,return_beta=False)
    assert selected['matrix_gather_calls']==1 and selected['matrix_gather_elements']==19
    assert selected['selected_array_bytes']==24*19


def test_finish_evidence_cannot_silently_cross_output_layouts():
    from test_result_service import fixture
    probe,kw=fixture()
    with pytest.raises(ValueError,match='layout'):result_finish_prices(probe,**kw,return_beta=False)
    probe.update(return_beta=False,return_df=False,result_layout=owned_result_work(32,1,return_beta=False)['array_bytes'])
    prices=result_finish_prices(probe,**kw,return_beta=False)['prices']
    assert prices['owned']['baseline_copy_bytes']==288
    assert prices['borrowed']['baseline_copy_bytes']==0
    assert prices['owned']['result_arrays']==3
    with pytest.raises(ValueError,match='layout'):result_finish_prices(probe,**kw)
    probe['result_layout']['beta']=128
    with pytest.raises(ValueError,match='three-array'):result_finish_prices(probe,**kw,return_beta=False)


def test_candidate_preparation_uses_output_choice_without_relabeling_evidence(tmp_path):
    from test_trait_candidate_space import spec,prepare
    from test_mechanistic_shapes import component
    value=spec(tmp_path);value['output']['store_beta']=False
    space=prepare(value)
    for candidate in space['candidates']:
        for tile in candidate['tiles']:
            assert tile['profile']['return_beta'] is False
            assert 'return_beta' not in tile['profile']['result_finish_service']
            with pytest.raises(ValueError,match='three-array'):
                torch_scan_work(tile['data'],tile['profile'])


def test_owned_t_only_copy_service_requires_an_explicit_finish_baseline():
    from torchgwas.owned_result_work import owned_result_copy_service
    work=owned_result_work(1024,9,return_beta=False)
    scenario=dict(resident_cpu_seconds_per_byte=1e-9,fresh_cpu_seconds_per_byte=2e-9,fresh_fraction=0.)
    with pytest.raises(ValueError,match='explicit'):owned_result_copy_service(work,scenario)
    result=owned_result_copy_service(work,scenario,baseline_copy_bytes=288)
    assert result['additional_copy_cpu_seconds']==pytest.approx((work['copy_bytes']-288)*1e-9)


@pytest.mark.parametrize('borrow',[False,True])
def test_scan_model_demands_correct_finish_then_counts_three_arrays(borrow):
    from test_mechanistic_shapes import fixture,component
    data,profile=fixture(k=512)
    profile.update(return_beta=False,result_ownership='borrowed' if borrow else 'owned')
    with pytest.raises(ValueError,match='three-array'):torch_scan_work(data,profile)
    profile['control_primitives']={key:1e-6 for row in native_control_work().values() for key in row}
    profile['result_finish_service']=dict(cpu_seconds=1e-6,serial_cpu_seconds=1e-6,
        baseline_copy_bytes=0 if borrow else 288,return_beta=False,result_arrays=3,baseline_rows=32,return_df=False,
        replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        result=torch_scan_work(data,profile)
    assert sum(row['d2h_bytes'] for row in result['blocks'])==10*(4*512+5)
    for row in result['owned_result_work'].values():
        assert 'beta' not in row['work'].get('borrowed_array_bytes',row['work']['array_bytes'])
    profile['result_finish_service']['baseline_copy_bytes']=416
    with pytest.raises(ValueError,match='coverage'):torch_scan_work(data,profile)
