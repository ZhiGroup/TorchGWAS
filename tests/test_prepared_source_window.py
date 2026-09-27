"""Bounded windows retain admission identity and independent mutable ownership."""
from copy import deepcopy
import pytest

from test_pgen_work_bounds import fixture
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.window_model import prepared_source_window


@pytest.fixture
def source(tmp_path):
    path=tmp_path/'input.pgen';fixture(path,129);header=PgenHeaderWork(path)
    geometry=[dict(shape=[129,3,4],kernels=[dict(name='example')])]
    tile=dict(device='cuda:0',trait_range=[0,4],
        data=dict(samples=129,markers=15,traits_analyzed=4,covariates=2,
            encoded=dict(kind='torchgwas.pgen_memory_layout.v1',path=str(path),samples=129,
                variant_range=[0,15],discarded_full_index=[dict(record=i) for i in range(15)]),
            fixed=dict(columns=['x','y'])),
        profile=dict(chunk_markers=4,kernel_geometry=geometry,other_geometry=geometry,
            process_units=dict(bytes=1e-9)))
    return tile,header


def window(tile,header,**kw):
    return prepared_source_window(tile,header,**dict(dict(start=2,stop=13,chunk_markers=3,
        issued_chunks=7,expected_input_identity=header.input_identity),**kw))


@pytest.mark.parametrize('embedded',[False,True])
def test_exact_old_construction_with_independent_mutable_ownership(source,embedded):
    tile,header=source
    if embedded:tile['data']['encoded']['input_identity']=deepcopy(header.input_identity)
    before=deepcopy(tile)
    expected={key:deepcopy(tile[key]) for key in ('device','trait_range','data','profile')}
    expected['issued_chunks']=7;expected['profile']['chunk_markers']=3
    expected['data'].update(markers=11,encoded=header.window(2,13,3))
    result=window(tile,header,expected_input_identity=None if embedded else header.input_identity)
    assert result==expected and tile==before
    assert result['profile']['kernel_geometry'] is result['profile']['other_geometry']
    result['profile']['kernel_geometry'][0]['shape'][0]=-1
    result['data']['fixed']['columns'].clear();result['trait_range'][0]=1
    result['data']['encoded']['input_identity'].clear()
    assert tile==before and window(tile,header)==expected


def test_discarded_full_source_description_is_never_copied(source):
    class NoCopy:
        def __deepcopy__(self,memo):raise AssertionError('Discarded encoded source was copied')
    tile,header=source;tile['data']['encoded']['discarded_full_index']=NoCopy()
    result=window(tile,header)
    assert result['data']['encoded']['variant_range']==[2,13]


@pytest.mark.parametrize('change',[
    lambda t,h:t['data']['encoded'].update(input_identity={}),
    lambda t,h:t['data'].update(samples=130),
    lambda t,h:t['data']['encoded'].update(samples=130),
    lambda t,h:t['data']['encoded'].update(variant_range=[3,15]),
    lambda t,h:t['data']['encoded'].update(variant_range=[0,16]),
    lambda t,h:t['data']['encoded'].update(variant_range=[False,15]),
    lambda t,h:t['data']['encoded'].update(variant_range=[0]),
])
def test_mismatched_admission_refused(source,change):
    tile,header=source;change(tile,header)
    with pytest.raises(ValueError):window(tile,header)


@pytest.mark.parametrize('options',[
    dict(expected_input_identity=None),dict(expected_input_identity={}),
    dict(start=-1),dict(stop=16),dict(issued_chunks=-1),
    dict(max_chunks=3),dict(max_records=10),dict(max_signatures=1),
])
def test_identity_and_work_limits_refused(source,options):
    with pytest.raises(ValueError):window(*source,**options)


def test_other_path_and_replaced_admitted_input_refused(source,tmp_path):
    tile,header=source
    other=tmp_path/'other.pgen';fixture(other,129)
    tile['data']['encoded']['path']=str(other)
    with pytest.raises(ValueError,match='admitted source'):window(tile,header)
    tile['data']['encoded']['path']=str(header.path)
    old=deepcopy(header.input_identity)
    with header.path.open('ab') as stream:stream.write(b'\0')
    with pytest.raises(ValueError,match='changed'):window(tile,header)
    replacement=PgenHeaderWork(header.path)
    with pytest.raises(ValueError,match='admitted source'):
        window(tile,replacement,expected_input_identity=old)
