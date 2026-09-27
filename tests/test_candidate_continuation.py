"""Counterfactual future chunks retain the started multi-GPU model prefix."""
import copy
import numpy as np
import pytest

from torchgwas.jagwas_candidate import jagwas_candidate_runtime
from torchgwas.adaptive_candidate import future_chunk_candidate
from torchgwas.pgen_work_census import census,scheduled_census
from test_jagwas_actual_candidate import actual_candidate,writer_prices
from test_jagwas_candidate import preparation
from test_pgen_native_reader import write_pgen
from test_adaptive_candidate import joint_fixture


def candidate_graph(choice):
    return jagwas_candidate_runtime(choice,writer_prices(),preparation=preparation(choice),
        occupancy='dense',host_serial_fraction=.5,return_graph=True)


def future_choice(choice,source,state,size):
    prefixes=[]
    for i,tile in enumerate(choice['tiles']):
        prefix=f'tile:{i}:submit_decode:'
        issued=max((int(name[len(prefix):]) for name in state['start'] if name.startswith(prefix)),default=-1)+1
        prefixes.append(issued)
    return future_chunk_candidate(choice,source_census=source,chunk_sizes=[128,256,512],
        next_size=size,issued_chunks=prefixes,reduction='jagwas'),prefixes


def test_actual_shape_counterfactual_preserves_issued_ranges_and_paid_factors(tmp_path):
    n,m=2049,4097;path=tmp_path/'input.pgen'
    write_pgen(path,(np.arange(n*m,dtype=np.uint32).reshape(m,n)%3).astype(np.uint8))
    choice=actual_candidate(path,512,2);source=census(path,128,include_chunks=True)
    original=candidate_graph(choice);full=original.solve()
    at=min(value for name,value in full['end'].items() if name.endswith('jagwas:0:0:write:done'))
    state=original.checkpoint(at);before=copy.deepcopy(state)
    assert state['remaining_service_seconds'] or state['fifo']
    original_remaining=original.resume(state)
    assert original_remaining['seconds']==pytest.approx(full['seconds'])
    for size in [128,256,512]:
        candidate,issued=future_choice(choice,source,state,size)
        assert sum(issued)<sum(len(t['data']['encoded']['chunks']) for t in choice['tiles'])
        changed=candidate_graph(candidate);result=changed.resume(state)
        assert result['remaining_seconds']>0
        for name in state['start']:
            assert result['start'][name]==state['start'][name]
        for name in state['end']:
            assert result['end'][name]==state['end'][name]
        for i,count in enumerate(issued):
            assert candidate['tiles'][i]['data']['encoded']['chunks'][:count]==choice['tiles'][i]['data']['encoded']['chunks'][:count]
        if size==512:assert result['seconds']==pytest.approx(full['seconds'])
    assert state==before


@pytest.mark.parametrize('damage',['issued_count','boolean','length','capacity','source','prefix','budget'])
def test_future_candidate_rejects_inconsistent_issued_prefix(tmp_path,damage):
    _,choice,source=joint_fixture(tmp_path)
    args=dict(source_census=source,chunk_sizes=[128,256,512],next_size=128,
              issued_chunks=[1],reduction='jagwas')
    if damage=='issued_count':args['issued_chunks']=[999]
    elif damage=='boolean':args['issued_chunks']=[True]
    elif damage=='length':args['issued_chunks']=[]
    elif damage=='capacity':args['chunk_sizes']=[128,256]
    elif damage=='source':source['file_bytes']+=1
    elif damage=='prefix':choice['tiles'][0]['data']['encoded']['chunks'][0]['record_payload_bytes']+=1
    else:args['max_source_chunks']=1
    with pytest.raises(ValueError):future_chunk_candidate(choice,**args)


def test_successive_future_revisions_preserve_issued_chunks_and_do_not_alias(tmp_path):
    _,choice,source=joint_fixture(tmp_path);old=copy.deepcopy((choice,source))
    args=dict(source_census=source,chunk_sizes=[128,256,512],reduction='jagwas')
    first=future_chunk_candidate(choice,next_size=128,issued_chunks=[1],**args)
    before=copy.deepcopy(first)
    second=future_chunk_candidate(first,next_size=256,issued_chunks=[3],**args)
    assert second['tiles'][0]['data']['encoded']['chunks'][:3]==first['tiles'][0]['data']['encoded']['chunks'][:3]
    second['tiles'][0]['profile']['chunk_markers']=123
    second['tiles'][0]['data']['encoded']['chunks'][0]['record_payload_bytes']=0
    assert first==before and (choice,source)==old
