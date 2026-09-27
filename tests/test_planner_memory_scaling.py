"""Memory admission must not materialize the scan/writer execution schedule."""
import itertools
from unittest.mock import patch

import numpy as np
import pytest

from test_trait_tiling_model import candidate, input_path, plan
from torchgwas.binary_output_work import binary_output_memory, binary_output_work
from torchgwas.trait_tiling_model import trait_tiled_memory, _tile_memory


def test_direct_memory_matches_enumerated_writer_for_empty_tail_and_borrowed_chunks():
    for m,k,b,block,beta,df,borrow in itertools.product(
            [0,1,7,11],[1,3,127],[1,4,20],[None,1,17,4096],
            [False,True],[False,True],[False,True]):
        memory=binary_output_memory(m,k,b,block,2,beta,df)
        ledger=binary_output_work(m,k,b,block,2,borrow,beta,store_variant_df=df)
        for key in ['block_bytes','queue_depth','arrays','pooled_buffers_per_array',
                    'allocated_staging_bytes','zero_initialization_bytes']:
            assert memory[key]==ledger[key], (m,k,b,block,beta,df,borrow,key)


@pytest.mark.parametrize('m,k,b,expected',[
    (0,600000,512,16<<20),
    (1,1,512,1<<20),
    (1,600000,512,2400000),
    (8931083,4096,512,8<<20),
    (8931083,600000,512,16<<20),
    (2**60+1,600000,512,16<<20),
])
def test_auto_memory_depends_only_on_first_chunk(m,k,b,expected):
    with patch('torchgwas.binary_output_work.binary_output_work',side_effect=AssertionError('expanded schedule')):
        result=binary_output_memory(m,k,b,queue_depth=3,store_variant_df=True)
    assert result['block_bytes_by_stream']==dict(beta=expected,t=expected,logp=expected,df=min(expected,1<<20))
    assert result['allocated_staging_bytes']==4*(3*expected+min(expected,1<<20))


@pytest.mark.parametrize('beta,df,borrow',itertools.product([False,True],repeat=3))
def test_memory_agrees_with_production_writer_allocations(tmp_path,monkeypatch,beta,df,borrow):
    from torchgwas import sumstats
    allocations=[]
    def allocate(size):
        allocations.append(size)
        return bytearray(size)
    monkeypatch.setattr(sumstats,'bytearray',allocate,raising=False)
    m,k,b,block,depth=7,3,4,17,2
    writer=sumstats.BinarySumstatsWriter(tmp_path,m,list(range(k)),32,22,
        block_bytes=block,queue_depth=depth,borrow_chunks=borrow,store_beta=beta,
        store_variant_df=df,fsync=False,writeback_bytes=0)
    for start in range(0,m,b):
        stop=min(start+b,m)
        writer.write_chunk(start,stop,np.ones((stop-start,k)),np.ones((stop-start,k)),
            variant_df=np.full((stop-start,1),22) if df else None)
    writer.close()
    memory=binary_output_memory(m,k,b,block,depth,beta,df)
    assert sum(allocations)==memory['allocated_staging_bytes']
    assert len(allocations)==memory['arrays']*(depth+1)


def test_memory_reuses_extents_and_shared_census_without_expanding_writer(input_path):
    c=candidate(input_path,width=2,count=2,traits=9)
    class Counted(list):
        iterations=0
        def __iter__(self):
            self.iterations+=1
            return super().__iter__()
    encoded=c['tiles'][0]['data']['encoded']
    encoded['chunks']=Counted(encoded['chunks'])
    for tile in c['tiles']:tile['data']['encoded']=encoded
    with patch('torchgwas.trait_tiling_model.binary_output_work',side_effect=AssertionError('expanded schedule')):
        with patch('torchgwas.trait_tiling_model._tile_memory',wraps=_tile_memory) as estimate:
            before=trait_tiled_memory(c)
            assert estimate.call_count==3  # full width on two devices and one tail
        assert encoded['chunks'].iterations==1
        for row in encoded['chunks']:row['record_payload_bytes']+=1000
        after=trait_tiled_memory(c)
    assert after['host_bytes']>before['host_bytes']  # no identity cache between calls
    assert after['device_bytes']==before['device_bytes']
    assert after['pinned_cache_bytes_by_device']==before['pinned_cache_bytes_by_device']


@pytest.mark.parametrize('reason',['device_memory','host_memory','shared_reader_budget'])
def test_graph_budget_counts_only_resource_feasible_candidates(input_path,reason):
    good=candidate(input_path,width=5,count=1)
    bad=candidate(input_path,width=1,count=2)
    kwargs=dict(max_scenario_evaluations=2,max_tiles=2,max_chunk_evaluations=6)
    if reason=='device_memory':kwargs['device_memory_bytes']={'cuda:0':2**30,'cuda:1':1}
    elif reason=='host_memory':
        bad['output']['block_bytes']=1<<30
        kwargs['host_memory_bytes']=trait_tiled_memory(good)['host_bytes']
    else:
        for tile in bad['tiles']:tile['profile']['decode_workers']=2
    expected=plan([good],**kwargs)
    actual=plan([bad,good],**kwargs)
    assert actual['selected']['candidate_index']==1
    assert actual['selected']['worst_supplied_scenario_seconds']==expected['selected']['worst_supplied_scenario_seconds']
    assert actual['rejected'][0]['reason']==reason
    assert actual['candidates_evaluated']==2 and actual['candidates_feasible']==1


@pytest.mark.parametrize('limit',['max_scenario_evaluations','max_tiles','max_chunk_evaluations'])
def test_complete_feasible_budget_is_checked_before_any_graph(input_path,limit):
    proposals=[candidate(input_path,count=1),candidate(input_path)]
    with patch('torchgwas.trait_tiling_plan.torch_trait_tiled_runtime',side_effect=AssertionError('started graph')):
        with pytest.raises(ValueError,match=limit):plan(proposals,**{limit:1})


def test_rejected_candidate_order_remains_input_order(input_path):
    unsupported=candidate(input_path,count=1)
    unsupported['tiles'][0]['profile'].pop('setup_primitives')
    memory=candidate(input_path)
    good=candidate(input_path,count=1)
    result=plan([unsupported,memory,good],device_memory_bytes={'cuda:0':2**30,'cuda:1':1})
    assert [(row['candidate_index'],row['reason']) for row in result['rejected']]==[
        (0,'unsupported_model_context'),(1,'device_memory')]
