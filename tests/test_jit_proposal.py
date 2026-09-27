import copy

import pytest

from torchgwas.adaptive_candidate import future_chunk_candidate
from torchgwas.jit_proposal import analytical_chunk_proposal
from torchgwas.remaining_work import remaining_work_bounds,source_wait_chains
from test_adaptive_candidate import joint_fixture
from test_candidate_continuation import candidate_graph


def snapshot(candidate):
    rows=[]
    for i,tile in enumerate(candidate['tiles']):
        lo,hi=tile['data']['encoded']['variant_range']
        ranges=[[lo,lo+128],[lo+128,lo+256],[lo+256,lo+512]]
        rows.append(dict(id=str(i),device=tile['device'],variant_range=[lo,hi],
            trait_range=tile['trait_range'],ranges=ranges,cursor=lo+512,issued_chunks=3))
    return dict(partitions=rows,prefix_complete=True,finished=None,first_written=1.,
                issued_revision=3*len(rows),current_chunk_size=256)


def test_actual_prefix_replaces_regular_model_prefix_without_reading_input(tmp_path):
    _,candidate,source=joint_fixture(tmp_path)
    state=snapshot(candidate);before=copy.deepcopy((candidate,source,state))
    arguments=dict(source_census=source,chunk_sizes=[128,256,512],issued_chunks=[3],
        issued_ranges=[state['partitions'][0]['ranges']],reduction='jagwas')
    choice=future_chunk_candidate(candidate,next_size=512,**arguments)
    ranges=[row['variant_range'] for row in choice['tiles'][0]['data']['encoded']['chunks']]
    assert ranges[:3]==state['partitions'][0]['ranges']
    assert ranges[3][0]==512 and ranges[-1][1]==source['markers']
    assert (candidate,source,state)==before
    # The number of actual small reserved chunks can exceed a coarse plan.
    ranges=[];cursor=0
    while cursor<source['markers']:
        hi=min(cursor+128,source['markers']);ranges.append([cursor,hi]);cursor=hi
    future_chunk_candidate(candidate,next_size=512,**dict(arguments,issued_chunks=[len(ranges)],issued_ranges=[ranges]))


def test_source_calculator_generates_one_conditional_proposal(tmp_path):
    _,candidate,source=joint_fixture(tmp_path)
    # This sparse LD fixture exercises additional source units. Supply named
    # synthetic component prices explicitly, never infer them from GWAS time.
    for tile in candidate['tiles']:
        tile['profile']['decode_units'].update({name:1e-9 for name in
            ['uleb1','uleb2','set_category','difflist_group_absolute_id',
             'difflist_category_extract','difflist_record_header']})
    state=snapshot(candidate);before=copy.deepcopy((candidate,source,state));graphs=[]
    def build(choice):
        graph=candidate_graph(choice);graphs.append((choice,graph));return graph
    proposal,audit=analytical_chunk_proposal(candidate,source_census=source,snapshot=state,
        chunk_sizes=[128,256,512],next_size=512,reduction='jagwas',graph_factory=build)
    assert len(graphs)==2
    assert proposal['baseline_seconds']==audit['baseline']['lower_seconds']
    assert proposal['candidate_seconds']==audit['candidate']['upper_seconds']
    for choice,graph in graphs:
        assert [p['variant_range'] for p in choice['tiles'][0]['data']['encoded']['chunks'][:3]]==state['partitions'][0]['ranges']
        # For a real model checkpoint immediately before unissued submission,
        # the live-prefix bound contains the solved remaining model time.
        solved=graph.solve();at=solved['start']['tile:0:submit_decode:3']
        checkpoint=graph.checkpoint(max(0.,at-1e-10))
        remainder=graph.resume(checkpoint)['remaining_seconds']
        bound=remaining_work_bounds(graph,unissued_nodes=['tile:0:submit_decode:3'],
            wait_chains=source_wait_chains(graph,choice,indexed=True))
        assert bound['lower_seconds']<=remainder+1e-8<=bound['upper_seconds']+1e-8
    assert (candidate,source,state)==before


def test_missing_component_prices_are_not_filled_by_the_proposer(tmp_path):
    _,candidate,source=joint_fixture(tmp_path)
    with pytest.raises(ValueError,match='Unpriced decoder'):
        analytical_chunk_proposal(candidate,source_census=source,snapshot=snapshot(candidate),
            chunk_sizes=[128,256,512],next_size=512,reduction='jagwas',graph_factory=candidate_graph)


def test_productive_bridge_accepts_source_generated_proposal_only_after_write(tmp_path):
    import time
    from torchgwas.planning_session import IncrementalPlanningBudget
    from torchgwas.productive_run import ProductiveTuningRun
    from torchgwas.sumstats_indexed import IndexedChunkWrite
    _,candidate,source=joint_fixture(tmp_path)
    for tile in candidate['tiles']:
        tile['profile']['decode_units'].update({name:1e-9 for name in
            ['uleb1','uleb2','set_category','difflist_group_absolute_id',
             'difflist_category_extract','difflist_record_header']})
    tile=candidate['tiles'][0]
    run=ProductiveTuningRun([dict(id='tile',device=tile['device'],
        variant_range=tile['variant_range'],trait_range=tile['trait_range'])],
        chunk_sizes=[128,256,512],initial=128,
        budget=IncrementalPlanningBudget(max_cpu_seconds=5.,max_window_seconds=10.))
    audits=[]
    def build(state):
        proposal,audit=analytical_chunk_proposal(candidate,source_census=source,snapshot=state,
            chunk_sizes=[128,256,512],next_size=512,reduction='jagwas',graph_factory=candidate_graph)
        audits.append(audit);return proposal
    # A deliberately generous horizon tests wiring, not real planner payoff.
    forecasts=dict(remaining_seconds=1000.,expected_cpu_seconds=.01,expected_wall_seconds=.01)
    assert not run.planning_step(build,**forecasts)['evaluated'] and not audits
    control=run.for_partition('tile')
    for start in [0,128,256]:assert control(start,source['markers'],512)==128
    now=time.perf_counter()
    run.output_written(IndexedChunkWrite(0,128,'jagwas',128,123,'part.npz',now,now,True))
    result=run.planning_step(build,**forecasts)
    assert result['evaluated'] and 'error' in result and result['error'] is None
    assert audits[0]['issued_chunks']==[3]
    assert result['value']['baseline_seconds']==audits[0]['baseline']['lower_seconds']
    assert result['value']['candidate_seconds']==audits[0]['candidate']['upper_seconds']
    assert result['applied']==(result['usable_for_decision'] and audits[0]['conditional_gain_floor_seconds']>result['wall_seconds'])
    run.finish(successful=False)


@pytest.mark.parametrize('reduction,occupancy',[(None,None),('significant','empty'),('significant','dense')])
@pytest.mark.parametrize('devices',[1,2])
def test_dense_and_significant_trait_tiles_use_exact_prefix_and_serial_worker_chains(tmp_path,reduction,occupancy,devices):
    from torchgwas.pgen_work_census import census
    from test_scheduled_census import make_dense,dense_runtime,significant_runtime
    candidate=make_dense(tmp_path,count=devices)
    source=census(candidate['tiles'][0]['data']['encoded']['path'],2,include_chunks=True)
    state=dict(prefix_complete=True,finished=None,first_written=1.,issued_revision=1,current_chunk_size=2,
        partitions=[dict(id=str(i),device=tile['device'],trait_range=tile['trait_range'],variant_range=[0,10],
            ranges=[[0,2]] if i==0 else [],issued_chunks=1 if i==0 else 0,cursor=2 if i==0 else 0)
            for i,tile in enumerate(candidate['tiles'])])
    graphs=[]
    def build(choice):
        graph=(dense_runtime(choice,return_graph=True) if reduction is None else
               significant_runtime(choice,occupancy,return_graph=True))
        graphs.append((choice,graph));return graph
    proposal,audit=analytical_chunk_proposal(candidate,source_census=source,snapshot=state,
        chunk_sizes=[2,4],next_size=4,reduction=reduction,graph_factory=build)
    assert audit['issued_chunks']==[1,0,0] and len(graphs)==2
    assert all(row['wait_groups']<=devices*3 for row in [audit['baseline'],audit['candidate']])
    for choice,graph in graphs:
        bound=remaining_work_bounds(graph,unissued_nodes=[],wait_chains=source_wait_chains(graph,choice,indexed=reduction is not None))
        assert graph.solve()['seconds']<=bound['upper_seconds']+1e-8
    assert candidate['output']['store_beta'] is True and candidate['output']['fsync'] is True


@pytest.mark.parametrize('fault',['closed','truncated','before_output','cursor','count','range','device','duplicate','size','budget'])
def test_bad_prefix_or_outside_action_fails_before_graph_build(tmp_path,fault):
    _,candidate,source=joint_fixture(tmp_path);state=snapshot(candidate)
    arguments=dict(source_census=source,snapshot=state,chunk_sizes=[128,256,512],next_size=512,
                   reduction='jagwas',graph_factory=lambda choice:pytest.fail('unvalidated graph'))
    if fault=='closed':state['finished']={}
    elif fault=='truncated':state['prefix_complete']=False
    elif fault=='before_output':state['first_written']=None
    elif fault=='cursor':state['partitions'][0]['cursor']+=1
    elif fault=='count':state['partitions'][0]['issued_chunks']+=1
    elif fault=='range':state['partitions'][0]['ranges'][1][0]+=1
    elif fault=='device':state['partitions'][0]['device']='cuda:2'
    elif fault=='duplicate':state['partitions'].append(copy.deepcopy(state['partitions'][0]))
    elif fault=='size':arguments['next_size']=100
    else:arguments['max_source_chunks']=1
    with pytest.raises(ValueError):analytical_chunk_proposal(candidate,**arguments)
