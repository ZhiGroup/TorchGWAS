"""Header windows feed the existing scan equations without exact-count claims."""
from copy import deepcopy
from unittest.mock import patch
import pytest

from torchgwas.decoder_work import decoder_work
from torchgwas.mechanistic_torch import torch_scan_work,torch_scan_header_work
from torchgwas.pgen_work_bounds import PgenHeaderWork,WINDOW_KIND
from torchgwas.pgen_work_census import census
from test_pgen_work_bounds import fixture as write_fixture
from test_mechanistic_shapes import fixture as model_fixture,component


def fixture(tmp_path,chunk=4,span=(0,15)):
    path=tmp_path/'input.pgen';write_fixture(path,129)
    index=PgenHeaderWork(path);window=index.window(*span,chunk)
    data,profile=model_fixture(k=2,c=8)
    data.update(samples=129,markers=span[1]-span[0],encoded=window)
    profile.update(chunk_markers=chunk,kernel_geometry=[dict(N=129,B=b,K=2,C=8,kernels=[])
        for b in {hi-lo for lo,hi in window['chunk_ranges']}])
    names=sorted({name for row in window['chunks'] for name in row['source_units']})
    profile['decode_units']={name:(i+1)*1e-8 for i,name in enumerate(names)}
    return path,index,data,profile


@pytest.mark.parametrize('span',[(0,15),(3,15),(2,9)])
@pytest.mark.parametrize('buffered',[False,True])
def test_scenarios_enclose_exact_cpu_service_and_keep_other_work(tmp_path,span,buffered):
    path,index,data,profile=fixture(tmp_path,span=span)
    if buffered:profile['input_read_cpu_prices']=dict(cpu_seconds_per_byte=1e-9,cpu_seconds_per_call=1e-6)
    actual=deepcopy(data);actual['encoded']=census(path,4,variant_range=span,include_chunks=True)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        exact=torch_scan_work(actual,profile)
        low=torch_scan_header_work(data,profile,endpoint='lower')
        high=torch_scan_header_work(data,profile,endpoint='upper')
    assert low['cpu_work_seconds']<=exact['cpu_work_seconds']+1e-12<=high['cpu_work_seconds']+2e-12
    for lo,real,hi in zip(low['blocks'],exact['blocks'],high['blocks']):
        assert lo['decode_seconds']<=real['decode_seconds']<=hi['decode_seconds']
        for resource in ('cpu','dram'):
            amounts=[b['decode_seconds']*b['decode_resources'][resource] for b in (lo,real,hi)]
            assert amounts[0]-1e-9<=amounts[1]<=amounts[2]+1e-9
        # This includes read, transfers, tensor work, result ownership and all
        # finish/control work; only the decoder interval varies.
        for key in set(real)-{'decode_seconds','decode_resources'}:
            assert lo[key]==real[key]==hi[key],key
    assert low['host_workspace']==exact['host_workspace']==high['host_workspace']
    assert low['components']==exact['components']==high['components']
    assert low['source_scenario']['variant_range']==list(span)
    assert 'not elapsed-time bounds' in low['source_scenario']['scope']
    for call in (lambda:torch_scan_work(data,profile),lambda:decoder_work(data['encoded'],'torch_native_int8')):
        with pytest.raises(ValueError,match='exact'):call()


def test_resumed_window_does_not_charge_reader_initialization_again(tmp_path):
    _,_,data,profile=fixture(tmp_path,span=(3,15))
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        cold=torch_scan_header_work(data,profile,endpoint='upper')
        warm=torch_scan_header_work(data,profile,endpoint='upper',issued_chunks=4)
    assert all('reader_init_seconds' not in b for b in warm['blocks'])
    assert cold['cpu_work_seconds']-warm['cpu_work_seconds']==pytest.approx(2*15*profile['process_units']['pgen_index_records'])


def test_window_budgets_precede_per_chunk_work_and_input_changes_fail(tmp_path):
    path,index,data,profile=fixture(tmp_path)
    with patch.object(index,'bounds',side_effect=AssertionError('budget must reject first')):
        with pytest.raises(ValueError,match='chunk budget'):index.window(0,15,1,max_chunks=8)
        with pytest.raises(ValueError,match='record budget'):index.window(0,15,4,max_records=14)
    with patch('torchgwas.mechanistic_torch._torch_scan_work',side_effect=AssertionError('must not price')):
        with pytest.raises(ValueError,match='bounded|Bounded'):torch_scan_header_work(data,profile,endpoint='upper',max_chunks=2)
        with pytest.raises(ValueError,match='budget'):torch_scan_header_work(data,profile,endpoint='upper',max_records=14)
    with path.open('ab') as stream:stream.write(b'changed')
    with pytest.raises(ValueError,match='changed'):index.window(0,15,4)


@pytest.mark.parametrize('change',[
    lambda d:d['encoded']['chunks'][0].update(variant_range=[1,4]),
    lambda d:d['encoded']['chunks'][0].update(input_identity={}),
    lambda d:d['encoded'].update(path='/unbound'),
    lambda d:d['encoded'].update(file_markers=1),
    lambda d:d['encoded']['chunks'][0].update(read_bytes=-1),
    lambda d:d['encoded']['chunks'][0].update(decode_input_bytes=float('inf')),
    lambda d:d['encoded'].update(kind='exact'),
    lambda d:d['encoded'].update(chunk_markers=8)])
def test_window_binding_mismatch_is_rejected(tmp_path,change):
    _,_,data,profile=fixture(tmp_path);change(data)
    with pytest.raises(ValueError):torch_scan_header_work(data,profile,endpoint='upper')


def test_uncounted_possible_unit_refuses_pricing(tmp_path):
    _,_,data,profile=fixture(tmp_path)
    del profile['decode_units'][next(iter(profile['decode_units']))]
    with pytest.raises(ValueError,match='Unpriced'):torch_scan_header_work(data,profile,endpoint='lower')


def test_real_joint_geometry_and_tail_match_exact_scan_work(tmp_path):
    import numpy as np
    from test_pgen_native_reader import write_pgen
    from test_jagwas_actual_candidate import actual_candidate
    path=tmp_path/'joint.pgen'
    write_pgen(path,(np.arange(1025*2049,dtype=np.uint32).reshape(1025,2049)%3).astype(np.uint8))
    tile=actual_candidate(path,128,1)['tiles'][0]
    data=deepcopy(tile['data']);data['markers']=897
    data['encoded']=census(path,128,variant_range=(128,1025),include_chunks=True)
    exact=torch_scan_work(data,tile['profile'])
    data['encoded']=PgenHeaderWork(path).window(128,1025,128)
    for endpoint in ('lower','upper'):
        scenario=torch_scan_header_work(data,tile['profile'],endpoint=endpoint)
        assert scenario['blocks']==exact['blocks']
        assert scenario['components']==exact['components']
        assert sum(b['d2h_bytes'] for b in scenario['blocks'])==17*897
        assert scenario['components'][1]['reduction_component']['gemm']['useful_flops']==2*512**2
