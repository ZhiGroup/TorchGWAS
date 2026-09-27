"""A decision must retain paid setup, pending FIFOs, tokens and running work."""
import copy
import json
import random

import pytest

from torchgwas.execution_graph import ExecutionGraph
from test_streamed_graphs import graph
from test_jagwas_schedule import run,shard


def equivalent(g,when):
    before=copy.deepcopy(g.__dict__);full=g._solve_shared_python()
    checkpoint=g.checkpoint(when)
    saved=copy.deepcopy(checkpoint)
    # Round trips use the same types as persisted calibration/provenance JSON.
    result=g.resume(json.loads(json.dumps(checkpoint)))
    assert result['seconds']==pytest.approx(full['seconds'],rel=1e-11,abs=1e-12)
    assert result['remaining_seconds']==pytest.approx(full['seconds']-checkpoint['at_seconds'],abs=1e-12)
    for key in ('start','end','wait_resource_seconds'):
        assert result[key]==pytest.approx(full[key],rel=1e-11,abs=1e-12)
    for name,row in result['conditional_delays'].items():
        assert row['blocked']==full['conditional_delays'][name]['blocked']
        assert row['extra_elapsed_service_seconds']==full['conditional_delays'][name]['extra_elapsed_service_seconds']
    assert g.__dict__==before and checkpoint==saved
    return checkpoint,result


@pytest.mark.parametrize('policy',['fluid','held-first','held-last'])
@pytest.mark.parametrize('seed',range(12))
def test_resource_token_fifo_and_wakeup_state_resume_exactly(seed,policy):
    g=graph(seed)
    if policy!='fluid':g=g.with_serial_sections(policy)
    full=g._solve_shared_python()
    for fraction in [0.,.01,.2,.5,.85,1.,1.2]:equivalent(g,fraction*full['seconds'])
    for at in sorted(set(full['end'].values())):equivalent(g,at)


@pytest.mark.parametrize('seed',range(16))
def test_general_shared_dags_and_repeated_continuation(seed):
    rng=random.Random(seed);g=ExecutionGraph();g.capacities={'cpu':1.3,'dram':3.}
    g.token_capacities={'critical':2}
    for index in range(25):
        deps=rng.sample(list(g.nodes),min(len(g.nodes),rng.randrange(4)))
        resources={name:rng.choice([0.,.01,.3,1.,7.]) for name in g.capacities}
        name=str(index);g.add(name,rng.choice([0.,.0001,.002,.2]),deps,resources)
        if rng.random()<.5:g.token_actions[name]={'acquire':{'critical':1},'release_finish':{'critical':1}}
    full=g._solve_shared_python();state=None
    for fraction in [.1,.3,.5,.7,.9]:
        state=g.checkpoint(full['seconds']*fraction,initial=state)
        solved=g.resume(state)
        assert solved['end']==pytest.approx(full['end'],rel=1e-11,abs=1e-12)


@pytest.mark.parametrize('retained',[0,1])
def test_joint_multi_gpu_pending_parts_and_shared_preparation_survive(retained):
    g=run([shard(chunks=4,kernel=10.,writer=40.,retained=retained),
           shard('cuda:2',chunks=4,kernel=1.,writer=40.,retained=retained)],return_graph=True)
    result=g._solve_shared_python();saw_queue=saw_inflight=False
    for at in sorted(set(result['start'].values())|set(result['end'].values())):
        state,continued=equivalent(g,at)
        saw_queue|=bool(state['fifo']);saw_inflight|=bool(state['remaining_service_seconds'])
        if 'shared_prepare:complete' in state['end']:
            assert continued['end']['shared_prepare:complete']==7.
    assert saw_inflight and (saw_queue or not retained)


def test_running_work_is_nominal_service_and_capacity_changes_apply_only_from_now():
    g=ExecutionGraph();g.capacities={'cpu':1.}
    g.add('a',10.,resources={'cpu':1.});g.add('b',10.,resources={'cpu':1.})
    state=g.checkpoint(6.)
    assert state['remaining_service_seconds']=={'a':7.,'b':7.}
    faster=copy.deepcopy(g);faster.capacities['cpu']=2.
    result=faster.resume(state)
    assert result['remaining_seconds']==7. and result['seconds']==13.
    assert result['start']=={'a':0.,'b':0.}
    # A new scan would be ten seconds; neither replay nor wall subtraction is right.
    assert faster.solve()['seconds']==10.


def test_future_action_can_change_without_repaying_setup_or_restarting_running_work():
    g=ExecutionGraph();g.add('setup',5.);g.add('issued',4.,['setup']);g.add('future',10.,['issued'])
    state=g.checkpoint(6.)
    assert state['end']=={'setup':5.} and state['remaining_service_seconds']=={'issued':3.}
    baseline=g.resume(state)
    changed=copy.deepcopy(g);changed.nodes['future']=(2.,('issued',))
    result=changed.resume(state)
    assert baseline['remaining_seconds']==13. and result['remaining_seconds']==5.
    assert result['end']['issued']==9. and result['start']['setup']==0.
    # Replacing an unstarted future suffix may add nodes, while its paid prefix stays fixed.
    changed.nodes.pop('future');changed.add('new1',1.,['issued']);changed.add('new2',3.,['new1'])
    assert changed.resume(state)['remaining_seconds']==7.


def test_completed_checkpoint_does_no_work_and_empty_graph_is_valid():
    g=ExecutionGraph();g.add('done',3.)
    state=g.checkpoint(100.)
    assert state['at_seconds']==3. and g.resume(state)['remaining_seconds']==0.
    assert ExecutionGraph().resume(ExecutionGraph().checkpoint(0.))['remaining_seconds']==0.


@pytest.mark.parametrize('damage',['duration','demand','dependency','token_capacity','past_wait'])
def test_started_contract_cannot_be_rewritten(damage):
    g=graph(4);state=g.checkpoint(.5*g._solve_shared_python()['seconds'])
    started=next(iter(state['start']))
    if damage=='duration':g.nodes[started]=(g.nodes[started][0]+1,g.nodes[started][1])
    elif damage=='demand':g.demands[started]={'cpu':.123}
    elif damage=='dependency':g.nodes[started]=(g.nodes[started][0],('missing',))
    elif damage=='token_capacity':g.token_capacities[next(iter(g.token_capacities))]+=1
    else:
        name=next(iter(state['past_resource_waits']))
        g.resource_waits[name]['resources']['cpu']+=1.
    with pytest.raises(ValueError):g.resume(state)


@pytest.mark.parametrize('damage',['version','fields','time','remaining','missing_running','unknown','token','ready','conditional'])
def test_invalid_checkpoint_is_refused(damage):
    g=graph(7);state=g.checkpoint(.5*g._solve_shared_python()['seconds'])
    if damage=='version':state['schema']='other'
    elif damage=='fields':state['extra']=True
    elif damage=='time':state['at_seconds']=-1.
    elif damage=='remaining':state['remaining_service_seconds'][next(iter(state['remaining_service_seconds']))]=1e9
    elif damage=='missing_running':state['remaining_service_seconds'].clear()
    elif damage=='unknown':state['start']['unknown']=0.
    elif damage=='token':state['token_free'][next(iter(state['token_free']))]+=1
    elif damage=='ready':state['ready']=[next(iter(state['start']))]
    else:state['conditional_delays']['unknown']={}
    with pytest.raises(ValueError):g.resume(state)


def test_checkpoint_time_is_finite_forward_and_boolean_is_not_a_time():
    g=ExecutionGraph();g.add('a',5.)
    for at in [-1.,float('inf'),float('nan'),True]:
        with pytest.raises(ValueError):g.checkpoint(at)
    with pytest.raises(ValueError):g.checkpoint(1.,initial=g.checkpoint(2.))
