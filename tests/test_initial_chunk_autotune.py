"""Deferred public lifecycle: useful output, exact partitions, cost and failure."""
from copy import deepcopy
import json
import threading
import time
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import torch

from torchgwas.api import run_linear_gwas
from torchgwas.adaptive_chunks import ChunkObservation
from torchgwas.initial_chunk_autotune import (PublicInitialChunkTuning,validate_initial_chunk_config,
    productive_api_lifecycle,_ACTIVE)
from torchgwas.linear import multigpu_variant_ranges
from torchgwas.sumstats import open_binary_sumstats,open_binary_df
from torchgwas.sumstats_indexed import open_indexed_sumstats


def config(axis='trait',width=3):
    return dict(context='fixture',chunk_size=4,partition_axis=axis,trait_block=width,
        window_markers=[8,16,24],budget=dict(max_steps=2,max_cpu_seconds=2.,max_window_seconds=20.),
        cost_forecasts=dict(remaining_seconds=100.,expected_cpu_seconds=.01,expected_wall_seconds=.01,
            switching_seconds=0.,publication_seconds=0.,reserve_seconds=.01),
        forecast_options=dict(boundary_adjustments=dict(baseline=[0.,0.],candidate=[0.,0.]),
            relative_model_error=.05,max_slope_change=.1,max_extrapolation=20.,
            assumptions={key:'Declared synthetic control.' for key in
                ['source_work','output_occupancy','resource_capacity','partition_balance']}))


@pytest.mark.parametrize('fault',['initial','context','axis','width','horizons','budget','cost','nan','forecast','assumptions','scenarios','occupancy','cache','stage'])
def test_invalid_deferred_controls_fail_before_execution(fault):
    cfg=config();bounds=dict(chunks=[4,8],trait_blocks=[3]);contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'host':dict(host_serial_fraction=.5,host_serial_policy='fluid')},occupancy_scenarios={'all':'dense'})
    if fault=='initial':cfg['chunk_size']=5
    if fault=='context':cfg['context']='unbound'
    if fault=='axis':cfg['partition_axis']='variant'
    if fault=='width':cfg['trait_block']=1
    if fault=='horizons':cfg['window_markers']=[8,16,23]
    if fault=='budget':cfg['budget']['max_steps']=0
    if fault=='cost':cfg['cost_forecasts'].pop('switching_seconds')
    if fault=='nan':cfg['cost_forecasts']['reserve_seconds']=float('nan')
    if fault=='forecast':cfg['forecast_options']['relative_model_error']=1.
    if fault=='assumptions':cfg['forecast_options']['assumptions']={}
    if fault=='scenarios':joint['host_scenarios']={}
    if fault=='occupancy':joint['occupancy_scenarios']={'sample':'infer_from_p_value'}
    if fault=='cache':cfg['structural_cache_dir']=' '
    if fault=='stage':cfg['stage_observations']=dict(max_chunks_per_device=0,warmup_chunks=1,
        stride=1,max_window_seconds=10.,cuda_events=True,measurement_reserve_seconds=.01)
    with pytest.raises(ValueError):validate_initial_chunk_config(cfg,bounds,contexts,'significant',joint)


def test_jagwas_deferred_configuration_never_tiles_phenotypes():
    bounds=dict(chunks=[4,8]);contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')},occupancy_scenarios={'all':'dense'})
    validate_initial_chunk_config(config('variant',None),bounds,contexts,'jagwas',joint)
    with pytest.raises(ValueError,match='JAGWAS'):validate_initial_chunk_config(config(),bounds,contexts,'jagwas',joint)


@pytest.mark.parametrize('reduction,axis,width',[
    ('significant','trait',3),('jagwas','variant',None)])
def test_deferred_configuration_accepts_explicit_partial_output(reduction,axis,width):
    bounds=dict(chunks=[4,8],trait_blocks=[3])
    contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')},
        occupancy_scenarios={'partial':dict(retained_fraction=[1,100],placement='spread')})
    validate_initial_chunk_config(config(axis,width),bounds,contexts,reduction,joint)


def test_conditional_capacity_scenarios_are_bounded_and_include_nominal():
    bounds=dict(chunks=[4,8],trait_blocks=[3]);contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')})
    cfg=config();nominal=dict(cpu=1.,dram=1.,input=1.,output=1.)
    cfg['capacity_scenarios']={'nominal':nominal,'cpu_half':dict(nominal,cpu=.5)}
    validate_initial_chunk_config(cfg,bounds,contexts,None,joint)
    for scenarios in ({'half':dict(nominal,cpu=.5)},
            {'nominal':nominal,'invalid':dict(nominal,cpu=0.)},
            {'nominal':nominal,'invalid':dict(nominal,input=1.01)}):
        cfg['capacity_scenarios']=scenarios
        with pytest.raises(ValueError,match='capacity scenarios'):
            validate_initial_chunk_config(cfg,bounds,contexts,None,joint)
    cfg['capacity_scenarios']={'nominal':nominal,**{str(i):dict(nominal,cpu=.5) for i in range(4)}}
    with pytest.raises(ValueError,match='One to four'):
        validate_initial_chunk_config(cfg,bounds,contexts,None,joint)


def test_background_planning_configuration_is_explicit_boolean():
    bounds=dict(chunks=[4,8],trait_blocks=[3]);contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')})
    cfg=config();cfg['background_planning']=True
    validate_initial_chunk_config(cfg,bounds,contexts,None,joint)
    cfg['background_planning']='yes'
    with pytest.raises(ValueError,match='background_planning'):
        validate_initial_chunk_config(cfg,bounds,contexts,None,joint)


def test_native_staged_screen_config_requires_bounded_admitted_chunks():
    bounds=dict(chunks=[4,8],trait_blocks=[3]);contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,
        host_serial_policy='fluid')})
    cfg=config()
    cfg['source_staging']=dict(records_per_step=8,max_steps=32,
        max_cpu_seconds=10.,max_window_seconds=10.,
        max_retained_bytes=1<<20,extra_host_reserve_bytes=2<<20)
    cfg['staged_screen']=dict(chunk_sizes=[4,8],occupancy_scenario=None,
        max_partitions=2,max_unique_records=100,
        max_chunks_per_partition=100,max_cpu_seconds=5.,
        max_wall_seconds=10.,max_rebases=1)
    validate_initial_chunk_config(cfg,bounds,contexts,None,joint)
    invalid=deepcopy(cfg);invalid['staged_screen']['chunk_sizes']=[8]
    with pytest.raises(ValueError,match='bounded admitted chunks'):
        validate_initial_chunk_config(invalid,bounds,contexts,None,joint)
    invalid=deepcopy(cfg);invalid['staged_screen']['max_rebases']=3
    with pytest.raises(ValueError,match='rebases'):
        validate_initial_chunk_config(invalid,bounds,contexts,None,joint)
    invalid=deepcopy(cfg);del invalid['source_staging']
    with pytest.raises(ValueError,match='source staging'):
        validate_initial_chunk_config(invalid,bounds,contexts,None,joint)


def test_first_output_stages_source_without_uncharged_layout_switch(monkeypatch):
    from torchgwas.sumstats import DenseWriteProgress
    cfg=config();cfg['forecast_options']['max_extrapolation']=2.
    cfg['source_staging']=dict(records_per_step=8,max_steps=32,
        max_cpu_seconds=10.,max_window_seconds=10.,
        max_retained_bytes=1<<20,extra_host_reserve_bytes=2<<20)
    bounds=dict(chunks=[4,8],trait_blocks=[3]);contexts=[dict(name='fixture')]
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')})
    validate_initial_chunk_config(cfg,bounds,contexts,None,joint)
    calls=[]
    class Stage:
        def __init__(self,*args,**kwargs):
            assert callable(kwargs['on_step_observation'])
            calls.append('created')
        def output_written(self):calls.append('written')
        def finish(self):calls.append('finished');return dict(complete=False)
    monkeypatch.setattr('torchgwas.initial_chunk_autotune.ProductiveSourceStage',Stage)
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',reduction=None)
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});bound=tuner.for_partition('cuda:0',[0,257],[0,3])
    assert bound(0,257,8)==4
    tuner.output_written(DenseWriteProgress(0,4,(0,3),48,4,time.perf_counter(),'test','cuda:0'))
    assert bound(4,257,8)==4
    tuner.output_written(DenseWriteProgress(4,8,(0,3),48,8,time.perf_counter(),'test','cuda:0'))
    assert calls==['created','written','written']
    observed=tuner.run.snapshot()
    assert observed['stop_reason'] is None and observed['prefix_complete']
    assert observed['partitions'][0]['ranges']==[[0,4],[4,8]]
    assert observed['current_chunk_size']==4
    tuner.finish(successful=False)
    assert tuner.audit['productive']['source_staging']==dict(complete=False)
    assert calls==['created','written','written','finished']


def test_staged_first_chunks_keep_exact_frontier_beyond_old_tile_limit(monkeypatch):
    from torchgwas.sumstats import DenseWriteProgress
    cfg=config()
    cfg['source_staging']=dict(records_per_step=8,max_steps=32,
        max_cpu_seconds=10.,max_window_seconds=10.,
        max_retained_bytes=1<<20,extra_host_reserve_bytes=2<<20)
    class Stage:
        def __init__(self,*args,**kwargs):
            assert callable(kwargs['on_step_observation'])
        def output_written(self):pass
        def finish(self):return dict(complete=False)
    monkeypatch.setattr('torchgwas.initial_chunk_autotune.ProductiveSourceStage',Stage)
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',reduction=None)
    parts=[dict(id=str(i),device='cuda:0',variant_range=[0,257],trait_range=[i,i+1])
           for i in range(9)]
    start=dict(partitions=parts,chunk_sizes=[4,8],initial_size=4,
               input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{}, {})
    bound=tuner.for_partition('cuda:0',[0,257],[0,1])
    assert bound(0,257,8)==4
    tuner.output_written(DenseWriteProgress(0,4,(0,1),32,4,
        time.perf_counter(),'test','cuda:0'))
    assert bound(4,257,8)==4
    state=tuner.run.snapshot()
    assert state['stop_reason'] is None and state['prefix_complete']
    assert state['partitions'][0]['ranges']==[[0,4],[4,8]]
    assert state['current_chunk_size']==4
    tuner.finish(successful=False)


def test_background_proposal_releases_writer_and_discards_stale_source(monkeypatch):
    from torchgwas.sumstats import DenseWriteProgress
    cfg=config();cfg['background_planning']=True
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',
        reduction=None)
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});bound=tuner.for_partition('cuda:0',[0,257],[0,3])
    assert bound(0,257,8)==4
    entered=threading.Event();resume=threading.Event()
    def propose(snapshot,size):
        entered.set()
        assert resume.wait(5.)
        return dict(chunk_size=size,baseline_seconds=20.,candidate_seconds=10.)
    monkeypatch.setattr(tuner,'_propose',propose)
    try:
        tuner.output_written(DenseWriteProgress(0,4,(0,3),48,4,time.perf_counter(),'test'))
        assert entered.wait(5.)
        assert bound(4,257,8)==4
    finally:
        resume.set()
    tuner._planner_worker.join(5.)
    decision=tuner.run.snapshot()['decisions'][-1]
    assert decision['stale_frontier'] and not decision['applied']
    assert decision['planned_issued_revision']==1 and decision['issued_revision']==2
    assert tuner.run.snapshot()['current_chunk_size']==4
    assert tuner._retry_size==8
    # The same admitted candidate is retried against a new exact frontier.
    tuner.output_written(DenseWriteProgress(4,8,(0,3),48,4,time.perf_counter(),'test'))
    tuner.finish(successful=False)
    decisions=tuner.audit['productive']['decisions']
    assert len(decisions)==2 and decisions[-1]['applied']
    assert tuner.audit['productive']['current_chunk_size']==8


def test_productive_capacity_scenarios_keep_profiles_fixed_and_pair_gains(monkeypatch):
    cfg=config();nominal=dict(cpu=4.,dram=100.,input=100.,output=100.)
    cfg['capacity_scenarios']={'nominal':dict.fromkeys(nominal,1.),
        'cpu_half':dict(cpu=.5,dram=1.,input=1.,output=1.)}
    joint=dict(host_scenarios={'h':dict(host_serial_fraction=0.,host_serial_policy='fluid')})
    owner=SimpleNamespace(config=dict(initial_chunks=cfg,joint=joint),input_path='test.pgen',
        reduction=None,profile=dict(source_sha256={}))
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,input_file_identity={},
        candidate=dict(tiles=[dict(profile=dict(cpu_available_cores=4.))],output={},shared_capacities=nominal))
    original=deepcopy(start)
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    monkeypatch.setattr(tuner,'_check',lambda:None)
    monkeypatch.setattr('torchgwas.pgen_work_bounds.PgenHeaderWork',lambda path:SimpleNamespace(input_identity={}))
    monkeypatch.setattr('torchgwas.window_model.prepared_source_window',lambda *a,**kw:kw)
    seen=[]
    def compare(*args,**kwargs):
        seen.append(kwargs['shared_capacities']['cpu'])
        return dict(cpu=kwargs['shared_capacities']['cpu'])
    monkeypatch.setattr('torchgwas.window_model.compare_prepared_windows',compare)
    monkeypatch.setattr('torchgwas.price_binding.validate_comparison_prices',lambda *a:{'checked':True})
    def proposal(snapshot,comparisons,**kwargs):
        gain=10. if comparisons[0]['cpu']==4. else 2.
        return dict(chunk_size=8,baseline_seconds=50.,candidate_seconds=50.-gain),dict(forecast_status='stable_scenario')
    monkeypatch.setattr('torchgwas.productive_forecast.productive_window_proposal',proposal)
    result=tuner._propose(tuner.run.snapshot(),8)
    assert seen==[4.]*3+[2.]*3
    assert start==original
    assert result==dict(chunk_size=8,baseline_seconds=50.,candidate_seconds=48.)
    assert [row['capacity_scenario'] for row in tuner._attempts[0]['scenarios']]==['nominal','cpu_half']
    tuner.finish(successful=False)


def test_productive_stage_observer_uses_first_partition_per_gpu():
    cfg=config();cfg['stage_observations']=dict(max_chunks_per_device=2,warmup_chunks=0,
        stride=1,max_window_seconds=10.,cuda_events=True,measurement_reserve_seconds=.01)
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction=None)
    partitions=[dict(id='a',device='cuda:0',variant_range=[0,64],trait_range=[0,3]),
                dict(id='b',device='cuda:0',variant_range=[0,64],trait_range=[3,6]),
                dict(id='c',device='cuda:1',variant_range=[0,64],trait_range=[0,3])]
    start=dict(partitions=partitions,chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});sample=tuner.stage_sample
    assert tuner.stage_observer('cuda:0',(0,64),(0,3)) is sample
    assert tuner.stage_observer('cuda:0',(0,64),(3,6)) is None
    assert tuner.for_scan('cuda:1',(0,64),3,reader_workers=1,capacity=8,depth=2) is sample
    assert set(sample.snapshot()['devices'])=={'cuda:0','cuda:1'}
    tuner.finish(successful=False)
    assert tuner.audit['productive']['stage_sample']['new_measurements_stopped']


@pytest.mark.parametrize('declared,expected',[ (.01,.07), (.10,.11) ])
def test_productive_cost_gate_charges_known_probe_wall_if_reserve_is_too_small(
        declared,expected):
    cfg=config();cfg['stage_observations']=dict(max_chunks_per_device=2,warmup_chunks=0,
        stride=1,max_window_seconds=10.,cuda_events=False,
        measurement_reserve_seconds=declared)
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction=None)
    start=dict(partitions=[dict(id='a',device='cuda:0',variant_range=[0,64],
        trait_range=[0,3])],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});sample=tuner.stage_sample
    assert sample.reserve_read(0,4,'cuda:0')
    sample(ChunkObservation(start=0,end=4,capacity=4,device='cuda:0',
        read_started=1.,read_finished=2.,submitted=3.,first_result=4.,
        completed=5.,consumer_seconds=.1,result_blocks=1,result_bytes=16,
        read_probe_wall_seconds=.02,consumer_probe_wall_seconds=.04))
    assert sample.known_probe_wall_seconds()==pytest.approx(.06)
    captured={}
    def planning_step(build,**costs):
        captured.update(costs)
        return dict(evaluated=False,usable_for_decision=False)
    with patch.object(tuner.run,'planning_step',planning_step):
        tuner._proposal_step(8,release_issue_frontier=False)
    assert captured['reserve_seconds']==pytest.approx(expected)
    tuner.finish(successful=False)
    assert tuner.audit['productive']['stage_sample']['known_probe_wall_seconds']==pytest.approx(.06)


def test_planning_history_directory_can_reuse_file_digests(tmp_path):
    cfg=config()
    cfg['planning_cost_history']=dict(cache_dir=str(tmp_path/'history'),
        max_age_seconds=3600.,publication_seconds=.01)
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction=None)
    start=dict(partitions=[dict(id='a',device='cuda:0',variant_range=[0,64],
        trait_range=[0,3])],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{}, {})
    cache=tuner._digest_cache()
    assert cache.cache.directory==tmp_path/'history'/'binding-digests-v1'
    tuner.finish(successful=False)

    cfg['binding_digest_cache_dir']=None
    tuner=PublicInitialChunkTuning(owner,start,{}, {})
    assert tuner._digest_cache() is None
    tuner.finish(successful=False)


def test_lifecycle_closes_all_registered_runs_on_exception_and_nested_api():
    closed=[]
    @productive_api_lifecycle
    def inner():
        _ACTIVE.get().append(SimpleNamespace(finish=lambda **kw:closed.append(('inner',kw))))
        raise RuntimeError('writer failed')
    @productive_api_lifecycle
    def outer():
        _ACTIVE.get().append(SimpleNamespace(finish=lambda **kw:closed.append(('outer',kw))))
        inner()
    with pytest.raises(RuntimeError,match='writer failed'):outer()
    assert closed==[('inner',dict(successful=False)),('outer',dict(successful=False))]
    assert _ACTIVE.get() is None


def pgen(tmp_path):
    from test_pgen_native_reader import write_pgen
    rng=np.random.default_rng(31122);n,m,k=97,257,7
    values=rng.integers(0,3,(m,n),dtype=np.uint8);values[4]=1;values[5,::9]=3
    path=tmp_path/'input.pgen';write_pgen(path,values)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    y=rng.normal(size=(n,k)).astype(np.float32);c=rng.normal(size=(n,2)).astype(np.float32)
    return path,y,c


def output(path,reduction):
    if reduction:
        _,parts=open_indexed_sumstats(path/'sumstats');parts=list(parts)
        values={key:np.concatenate([p[key] for p in parts]) for key in parts[0]}
        order=(np.lexsort((values['trait_index'],values['variant_index'])) if reduction=='significant'
               else np.argsort(values['variant_index']))
        return {key:value[order] for key,value in values.items()}
    beta,t,_=open_binary_sumstats(path/'sumstats')
    return dict(t_stat=np.asarray(t),df=np.asarray(open_binary_df(path/'sumstats')),
                **({} if beta is None else dict(beta=np.asarray(beta))))


@pytest.mark.parametrize('mode',['traits','variants','significant','jagwas'])
@pytest.mark.parametrize('decision',['profitable','stale','unprofitable','writer_failure'])
def test_public_deferred_switch_preserves_written_results_and_full_source(tmp_path,monkeypatch,mode,decision):
    if torch.cuda.device_count()<3:pytest.skip('CUDA 1 and 2 required')
    monkeypatch.setenv('TORCHGWAS_NATIVE_STATS','0');monkeypatch.setenv('TORCHGWAS_SIGNIFICANCE_BACKEND','host')
    path,y,c=pgen(tmp_path);n,k=y.shape;m=257;devices=['cuda:1','cuda:2']
    reduction=mode if mode in ('jagwas','significant') else None
    axis='variant' if mode in ('variants','jagwas') else 'trait'
    width=None if axis=='variant' else 3
    cfg=config(axis,width);audit={};tuned=[]
    if decision=='profitable':
        cfg['stage_observations']=dict(max_chunks_per_device=1,warmup_chunks=0,
            stride=1,max_window_seconds=20.,cuda_events=True,measurement_reserve_seconds=.01)
    settings=dict(chunk_size=8,reader_workers=2,prefetch_chunks=2,device='cuda:1')
    if axis=='trait':settings.update(trait_block=3,trait_devices=devices)
    else:settings['variant_devices']=devices
    partitions=[]
    if axis=='variant':
        for i,span in enumerate(multigpu_variant_ranges(m,8,2)):
            partitions.append(dict(id=str(i),device=devices[i],variant_range=list(span),trait_range=[0,k]))
    else:
        for i,lo in enumerate(range(0,k,3)):
            partitions.append(dict(id=str(i),device=devices[i%2],variant_range=[0,m],trait_range=[lo,min(k,lo+3)]))
    class Tuner:
        qc_trait_block=3
        devices=['cuda:1','cuda:2']
        productive=None
        def __init__(self,*a,**kw):pass
        def validate_inputs(self,*a,**kw):pass
        def select(self,genotype,phenotype,covariates,qc,*,output):
            owner=SimpleNamespace(config=dict(initial_chunks=cfg),reduction=reduction)
            start=dict(partitions=partitions,chunk_sizes=[4,8],initial_size=4)
            self.productive=PublicInitialChunkTuning(owner,start,audit,{})
            _ACTIVE.get().append(self.productive);tuned.append(self.productive)
            return settings,audit,None
    def propose(self,snapshot,size):
        assert snapshot['first_written'] is not None and snapshot['written_events']>0
        assert snapshot['current_chunk_size']==4
        self.test_prefix=deepcopy(snapshot['partitions'])
        if decision=='stale':raise ValueError('immutable price expired')
        return dict(chunk_size=size,baseline_seconds=90. if decision=='profitable' else 0.,candidate_seconds=0.)
    common=dict(genotype_format='pgen',pgen_mode='hardcall',compute_dtype='float32',sumstats_fields='t',
        reduce=reduction,significance_threshold=1. if reduction=='significant' else None,
        sumstats_block_bytes=None if reduction else 128,sumstats_queue_depth=1)
    run_linear_gwas(path,y,c,output_dir=tmp_path/'reference',**common,**settings)
    monkeypatch.setattr('torchgwas.detailed_autotune.DetailedAutotune',Tuner)
    monkeypatch.setattr(PublicInitialChunkTuning,'_propose',propose)
    if decision=='writer_failure':
        # Raise inside the real writer after output exists; native producers
        # and all registered productive controllers must still be closed.
        def failed_progress(self,event):
            self.run.output_written(event)
            raise OSError('injected completed-write callback failure')
        monkeypatch.setattr(PublicInitialChunkTuning,'output_written',failed_progress)
        with pytest.raises((OSError,RuntimeError)) as raised:
            run_linear_gwas(path,y,c,output_dir=tmp_path/'deferred',autotune_profile={},autotune_config={},**common)
        error=raised.value
        while error.__cause__ is not None:error=error.__cause__
        assert isinstance(error,OSError) and 'completed-write callback failure' in str(error)
        state=tuned[0].audit['productive']
        assert state['finished']['successful'] is False and state['written_events']>0
        assert _ACTIVE.get() is None
        import threading
        assert not any(t.name.startswith('torchgwas-') for t in threading.enumerate())
        return
    result=run_linear_gwas(path,y,c,output_dir=tmp_path/'deferred',autotune_profile={},autotune_config={},**common)
    state=result.run_metadata['autotune']['productive']
    assert state['finished']['successful'] is True
    assert state['output_boundary']['written_events']==state['written_events']
    assert state['output_boundary']['invalid_events']==0
    if decision=='profitable':
        stage=state['stage_sample']
        assert not stage['pending']
        assert {row['device'] for row in stage['observations']}==set(devices)
        assert all(row['cuda'] is not None for row in stage['observations'])
        assert state['decisions'][0]['total_tuning_cost_seconds']>=.02
    assert state['current_chunk_size']==(8 if decision=='profitable' else 4)
    assert state['planning']['steps'] and state['written_events']>0
    for p in state['partitions']:
        assert p['cursor']==p['variant_range'][1]
        prefix=next(row for row in tuned[0].test_prefix if row['id']==p['id'])['ranges']
        assert p['ranges'][:len(prefix)]==prefix
        assert all(hi-lo==4 for lo,hi in prefix)
        cursor=p['variant_range'][0]
        for lo,hi in p['ranges']:
            assert lo==cursor;cursor=hi
        if decision=='profitable':assert any(hi-lo==8 for lo,hi in p['ranges'])
    before,after=output(tmp_path/'reference',reduction),output(tmp_path/'deferred',reduction)
    assert before.keys()==after.keys()
    for key in before:
        np.testing.assert_allclose(after[key],before[key],rtol=3e-4,atol=3e-4,equal_nan=True)
    recorded=json.loads((tmp_path/'deferred'/'run.json').read_text())
    assert recorded['autotune']['productive']['finished']['successful'] is True
    import threading
    assert not any(t.name.startswith('torchgwas-') for t in threading.enumerate())


@pytest.mark.parametrize('fault',['expired','artifact','source','context','input','gpu','host','host_credit'])
def test_public_live_check_rejects_changed_evidence_without_changing_chunk(tmp_path,monkeypatch,fault):
    from torchgwas.calibration_cache import CalibrationParameterCache
    from torchgwas.detailed_calibration import bind_detailed_profile,sha256_file
    from torchgwas.sumstats import DenseWriteProgress
    now=[100.];monkeypatch.setattr('time.time',lambda:now[0])
    source={'model.py':'fixed-source'};execution=dict(devices={'cuda:0':dict(uuid='card-zero')})
    deps=dict(source_sha256=source,execution_context=execution,measurement_protocol=dict(operation='copy'))
    record=CalibrationParameterCache(tmp_path/'prices').store('cpu_capacity','copy',dict(rate=3.),
        dependencies=deps,provenance=dict(job='first'),observed_unix_seconds=90.,max_age_seconds=30.)
    contexts=[dict(name='fixture',devices=['cuda:0'],profiles={'cuda:0':dict(rate=3.)})]
    profile=bind_detailed_profile(contexts,execution,sources=source,limitations=['Synthetic test coefficient'],
        component_artifacts={record['path']:sha256_file(record['path'])},price_bindings=[dict(
            artifact=record['path'],kind='cpu_capacity',name='copy',dependencies=deps,max_age_seconds=None,
            targets=[dict(context_path=[0,'profiles','cuda:0','rate'],value_path=['rate'])])])
    owner=SimpleNamespace(config=dict(initial_chunks=config()),devices=['cuda:0'],input_path='test.pgen',
        output_path=tmp_path,profile=profile,reduction=None)
    start=dict(partitions=[dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])],
        chunk_sizes=[4,8],initial_size=4,input_file_identity={'unchanged':True},
        memory=dict(device_bytes={'cuda:0':100},host_bytes=100,
            retained_index_bases_bytes=20 if fault=='host_credit' else 0))
    tuner=PublicInitialChunkTuning(owner,start,{},dict.fromkeys(owner.devices,0))
    monkeypatch.setattr('torchgwas.detailed_calibration.source_identity',lambda:source)
    monkeypatch.setattr('torchgwas.detailed_calibration.execution_context',lambda *a,**kw:execution)
    monkeypatch.setattr('torchgwas.analytical_plan_cache.input_identity',lambda path:{'unchanged':True})
    monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(1000,2000))
    monkeypatch.setattr(torch.cuda,'memory_reserved',lambda d:0)
    monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:1000)
    tuner._check()
    original=Path(record['path']).read_bytes()
    if fault=='expired':now[0]=120.
    if fault=='artifact':Path(record['path']).write_text('changed')
    if fault=='source':source={'model.py':'new-source'}
    if fault=='context':execution=dict(devices={'cuda:0':dict(uuid='different-card')})
    if fault=='input':monkeypatch.setattr('torchgwas.analytical_plan_cache.input_identity',lambda path:{'unchanged':False})
    if fault=='gpu':monkeypatch.setattr(torch.cuda,'mem_get_info',lambda d:(50,2000))
    if fault=='host':monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:50)
    if fault=='host_credit':
        monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:80)
        tuner._check()  # The exact retained allocation was already paid for.
        monkeypatch.setattr('torchgwas.api._available_host_bytes',lambda:79)
    bound=tuner.for_partition('cuda:0',[0,257],[0,3]);assert bound(0,257,8)==4
    tuner.output_written(DenseWriteProgress(0,4,(0,3),48,4,time.perf_counter(),'test'))
    state=tuner.run.snapshot()
    assert state['current_chunk_size']==4 and state['stop_reason']=='planning_error'
    assert state['decisions'][-1]['error']=='ValueError'
    # Rejection is optional planning state; remaining source still advances.
    assert bound(4,257,8)==4
    if fault!='artifact':assert Path(record['path']).read_bytes()==original
    tuner.finish(successful=False)


@pytest.mark.parametrize('unstable',[False,True])
def test_public_forecast_requires_gain_across_every_declared_scenario(monkeypatch,unstable):
    cfg=config('variant',None)
    joint=dict(host_scenarios={'fast':dict(host_serial_fraction=0.),'slow':dict(host_serial_fraction=1.)},
        occupancy_scenarios={'none':'empty','all':'dense'})
    owner=SimpleNamespace(config=dict(initial_chunks=cfg,joint=joint),input_path='test.pgen',
        reduction='jagwas',significance_threshold=None,profile=dict(source_sha256={}),
        price_evidence=dict(record=dict(value=dict(writer_prices={}))))
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,input_file_identity={},
        candidate=dict(tiles=[{}],output={},shared_capacities=dict(cpu=1.,dram=1.,input=1.,output=1.)))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    checks=[];seen=[]
    monkeypatch.setattr(tuner,'_check',lambda:checks.append(True))
    monkeypatch.setattr('torchgwas.pgen_work_bounds.PgenHeaderWork',lambda path:SimpleNamespace(input_identity={}))
    monkeypatch.setattr('torchgwas.window_model.prepared_source_window',lambda *a,**kw:kw)
    def compare(*args,**kwargs):
        seen.append(kwargs)
        return dict(slow=kwargs['host_serial_fraction'],retained=kwargs['survivor_evidence']['bins'][0]['retained'])
    monkeypatch.setattr('torchgwas.window_model.compare_prepared_windows',compare)
    monkeypatch.setattr('torchgwas.price_binding.validate_comparison_prices',lambda *a:{'checked':True})
    def proposal(snapshot,comparisons,**kwargs):
        slow=comparisons[0]['slow'];dense=bool(comparisons[0]['retained'])
        return (dict(chunk_size=8,baseline_seconds=90.-20.*slow,candidate_seconds=10.+20.*dense),
            dict(forecast_status='unstable' if unstable and slow and dense else 'stable_scenario'))
    monkeypatch.setattr('torchgwas.productive_forecast.productive_window_proposal',proposal)
    result=tuner._propose(tuner.run.snapshot(),8)
    assert len(seen)==12 and len(checks)==2
    assert len(tuner._attempts[0]['scenarios'])==4
    assert result==(dict(chunk_size=4,baseline_seconds=0.,candidate_seconds=0.) if unstable else
                   dict(chunk_size=8,baseline_seconds=70.,candidate_seconds=30.))
    tuner.finish(successful=False)


@pytest.mark.parametrize('reduction',[None,'significant','jagwas'])
@pytest.mark.parametrize('case',['gain','loss','cost','prior_cost','unstable','horizon'])
@pytest.mark.parametrize('reverse_order',[False,True])
def test_productive_payback_uses_matched_conditions_and_all_costs(monkeypatch,reduction,case,reverse_order):
    # Conditional model controls, not empirical capacity or runtime evidence.
    clock=[100.]
    monkeypatch.setattr(time,'perf_counter',lambda:clock[0])
    monkeypatch.setattr(time,'thread_time',lambda:clock[0])
    cfg=config('variant' if reduction=='jagwas' else 'trait',None if reduction=='jagwas' else 3)
    cfg['cost_forecasts'].update(remaining_seconds=250.,switching_seconds=.125,
        publication_seconds=.125,reserve_seconds=.25)
    if case=='cost':cfg['cost_forecasts']['reserve_seconds']=1.5
    if case=='horizon':cfg['cost_forecasts']['remaining_seconds']=50.
    hosts=[('fast',dict(host_serial_fraction=0.)),('slow',dict(host_serial_fraction=1.))]
    occupancies=[('none','empty'),('all','dense')]
    joint=dict(host_scenarios=dict(hosts[::-1] if reverse_order else hosts))
    if reduction is not None:joint['occupancy_scenarios']=dict(occupancies[::-1] if reverse_order else occupancies)
    owner=SimpleNamespace(config=dict(initial_chunks=cfg,joint=joint),input_path='test.pgen',
        reduction=reduction,significance_threshold=.01 if reduction=='significant' else None,
        profile=dict(source_sha256={}),price_evidence=None if reduction is None else
        dict(record=dict(value=dict(writer_prices={}) if reduction=='jagwas' else {})))
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,input_file_identity={},
        candidate=dict(tiles=[{}],output={},shared_capacities=dict(cpu=1.,dram=1.,input=1.,output=1.)))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    def check():clock[0]+=.125
    monkeypatch.setattr(tuner,'_check',check)
    monkeypatch.setattr('torchgwas.pgen_work_bounds.PgenHeaderWork',lambda path:SimpleNamespace(input_identity={}))
    monkeypatch.setattr('torchgwas.window_model.prepared_source_window',lambda *a,**kw:kw)
    def compare(*args,**kwargs):
        evidence=kwargs['survivor_evidence']
        return dict(slow=bool(kwargs['host_serial_fraction']),
            dense=False if evidence is None else bool(evidence['bins'][0]['retained']))
    monkeypatch.setattr('torchgwas.window_model.compare_prepared_windows',compare)
    monkeypatch.setattr('torchgwas.price_binding.validate_comparison_prices',lambda *a:{'checked':True})
    def proposal(snapshot,comparisons,**kwargs):
        slow,dense=comparisons[0]['slow'],comparisons[0]['dense']
        before,after={(False,False):(10.,8.),(False,True):(20.,17.),
            (True,False):(100.,80.),(True,True):(200.,150.)}[slow,dense]
        if case=='loss' and not slow and not dense:after=11.
        return (dict(chunk_size=8,baseline_seconds=before,candidate_seconds=after),
            dict(forecast_status='unstable_marginal_cost' if case=='unstable' and slow else 'stable_scenario'))
    monkeypatch.setattr('torchgwas.productive_forecast.productive_window_proposal',proposal)
    bound=tuner.for_partition('cuda:0',[0,257],[0,3])
    assert bound(0,257,8)==4
    from torchgwas.sumstats import DenseWriteProgress
    from torchgwas.sumstats_indexed import IndexedChunkWrite,IndexedOutputPartition
    def written_event(lo,hi):
        if reduction is None:
            return DenseWriteProgress(lo,hi,(0,3),(hi-lo)*3*4,hi,clock[0],
                                      'test','cuda:0')
        return IndexedChunkWrite(lo,hi,reduction,0,0,None,clock[0],clock[0],
            False,IndexedOutputPartition('cuda:0',(0,257),(0,3)),(lo,hi))
    event=written_event(0,4)
    next_start=4;expected_prefix=[[0,4]]
    if case=='prior_cost':
        tuner.run.output_written(event)
        def previous(snapshot):
            clock[0]+=1.5
            return dict(chunk_size=4,baseline_seconds=0.,candidate_seconds=0.)
        prior=tuner.run.planning_step(previous,**cfg['cost_forecasts'])
        assert not prior['applied'] and prior['usable_for_decision']
        assert bound(4,257,8)==4
        next_start=8;expected_prefix.append([4,8])
        event=written_event(4,8)
    tuner.output_written(event)
    state=tuner.run.snapshot();decision=state['decisions'][-1]
    assert decision['applied'] is (case=='gain')
    assert state['current_chunk_size']==(8 if case=='gain' else 4)
    attempt=tuner._attempts[-1]
    if case=='horizon':
        assert decision['error']=='ValueError' and 'remaining horizon' in decision['message']
        # The limiting-gain baseline is only 10 seconds; another scenario is
        # what violates this horizon and must not disappear in aggregation.
        assert attempt['matched_scenarios']['maximum_baseline_seconds']>50.
    elif case=='unstable':
        assert decision['forecast_gain_seconds']==0. and 'matched_scenarios' not in attempt
    else:
        gain=-1. if case=='loss' else 2.
        assert decision['forecast_gain_seconds']==gain
        paired=attempt['matched_scenarios'];worst=attempt['scenarios'][paired['limiting_scenario_index']]
        assert paired['gain_floor_seconds']==gain and min(paired['scenario_gains_seconds'])==gain
        assert worst['host_scenario']=='fast'
        if reduction is not None:assert worst['occupancy_scenario']=='none'
        assert decision['total_tuning_cost_seconds']==(2. if case=='cost' else 2.25 if case=='prior_cost' else .75)
    # No issued work is redone, and a rejected decision does not stop science.
    assert state['partitions'][0]['ranges']==expected_prefix
    width=8 if case=='gain' else 4
    assert bound(next_start,257,8)==width
    assert bound(next_start+width,257,8)==width
    tuner.finish(successful=False)
    json.dumps(tuner.audit)  # The complete decision audit remains serializable.


@pytest.mark.parametrize('reduction,occupancy',[(None,None),('significant','dense'),
    ('significant',dict(retained_fraction=[1,100],placement='spread'))])
def test_public_controller_builds_real_dense_and_significant_forecasts(tmp_path,monkeypatch,reduction,occupancy):
    from test_adaptive_candidate import joint_fixture,memory_profile
    from test_significant_host_model import bank
    from torchgwas.adaptive_start import prepare_adaptive_start
    from torchgwas.calibration_cache import CalibrationParameterCache
    from torchgwas.detailed_calibration import bind_detailed_profile,sha256_file,validate_detailed_profile
    path,candidate,_=joint_fixture(tmp_path)
    template=deepcopy(candidate['tiles'][0]['profile'])
    template.update(reduction=None,result_ownership='borrowed',borrow_results=True,return_beta=True)
    # Joint owned-result finish measurements do not describe this dense ring.
    from torchgwas.native_control_work import native_control_work
    template['result_finish_service']=dict(cpu_seconds=1e-6,serial_cpu_seconds=0.,baseline_copy_bytes=0,
        replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False)
    template['control_primitives']={key:1e-6 for counts in native_control_work().values() for key in counts}
    template['writeback_service'].update(submit_seconds=1e-6,wait_seconds=1e-6,fadvise_seconds=1e-6)
    template['process_units']['bytearray_zero_bytes']=1e-10
    template['decode_units'].update({name:1e-8 for name in ['uleb1','uleb2','uleb3','uleb4','uleb5',
        'set_category','difflist_group_absolute_id','difflist_category_extract','difflist_record_header']})
    context=dict(name='fixture',devices=['cuda:0'],profiles={'cuda:0':template},
        shared_capacities=candidate['shared_capacities'])
    workload=dict(genotype=str(path),samples=2049,markers=1025,traits=512,covariates=2,
        matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
    output=dict(block_bytes=None if reduction else 1<<20,queue_depth=2,store_beta=True,fsync=True)
    start=prepare_adaptive_start(workload,context,chunk_sizes=[128,256],initial_size=128,
        partition_axis='trait',trait_block=512,reduction=reduction,output=output,cpu_workers=8,
        host_memory_bytes=1<<40,device_memory_bytes={'cuda:0':1<<40},host_reserve_bytes=1<<20,
        device_reserve_bytes=1<<20,device_memory_profiles={'cuda:0':memory_profile()},
        significance_threshold=.01 if reduction else None)
    execution=dict(devices={'cuda:0':dict(uuid='synthetic',sm_count=template['gpu_resources']['sm_count'],
        **template.get('reduction_gpu_properties',{}))});source={'test':'fixed'}
    deps=dict(source_sha256=source,execution_context=execution,measurement_protocol={'operation':'synthetic-control'})
    record=CalibrationParameterCache(tmp_path/'prices').store('cpu_capacity','copy',
        dict(rate=template['process_units']['numpy_copy_bytes']),dependencies=deps,
        provenance={'scope':'Accounting test, not empirical calibration'},
        observed_unix_seconds=time.time(),max_age_seconds=3600.)
    profile=bind_detailed_profile([context],execution,sources=source,limitations=['Synthetic prices'],
        component_artifacts={record['path']:sha256_file(record['path'])},price_bindings=[dict(
            artifact=record['path'],kind='cpu_capacity',name='copy',dependencies=deps,max_age_seconds=None,
            targets=[dict(context_path=[0,'profiles','cuda:0','process_units','numpy_copy_bytes'],value_path=['rate'])])])
    cfg=config('trait',512);cfg.update(chunk_size=128,window_markers=[256,512,768])
    cfg['capacity_scenarios']={'nominal':dict(cpu=1.,dram=1.,input=1.,output=1.),
        'cpu_half':dict(cpu=.5,dram=1.,input=1.,output=1.)}
    joint=dict(host_scenarios={'one':dict(host_serial_fraction=.5,host_serial_policy='fluid')})
    if reduction:joint['occupancy_scenarios']={'mode':occupancy}
    owner=SimpleNamespace(config=dict(initial_chunks=cfg,joint=joint),input_path=path,reduction=reduction,
        significance_threshold=.01 if reduction else None,profile=profile,
        price_evidence=None if reduction is None else dict(record=dict(value=bank())))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=workload),{})
    # Hardware is simulated in this accounting test. Source windows, tensor
    # launch geometry, writer graph, forecasts and immutable binding are real.
    monkeypatch.setattr(tuner,'_check',lambda:validate_detailed_profile(profile,execution,sources=source))
    from torchgwas.sumstats_indexed import IndexedChunkWrite
    from torchgwas.sumstats import DenseWriteProgress
    assert tuner.run.for_partition('0')(0,1025,256)==128
    now=time.perf_counter()
    event=(IndexedChunkWrite(0,128,'significant',128*512,1,'part.npz',now,now,True) if reduction else
           DenseWriteProgress(0,128,(0,512),128*512*8,128,now,'test'))
    tuner.run.output_written(event)
    proposal=tuner._propose(tuner.run.snapshot(),256)
    audit=tuner._attempts[0]['scenarios'][0]
    assert [row['capacity_scenario'] for row in tuner._attempts[0]['scenarios']]==['nominal','cpu_half']
    assert audit['price_evidence']['status']=='declared_targets_verified'
    assert audit['remaining_pairs']==897*512
    assert set(audit['forecasts'])=={'baseline','candidate'}
    assert proposal['chunk_size'] in (128,256)
    tuner.finish(successful=False)


def test_short_remaining_work_skips_source_window_construction(monkeypatch):
    cfg=config()
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',
        reduction=None)
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,dict(workload=dict(traits=3)),{})
    checked=[]
    monkeypatch.setattr(tuner,'_check',lambda:checked.append(True))
    monkeypatch.setattr('torchgwas.pgen_work_bounds.PgenHeaderWork',
                        lambda *a,**kw:(_ for _ in ()).throw(AssertionError('source window built')))
    snapshot=dict(written_events=0,partitions=[dict(cursor=249,variant_range=[0,257])])
    with pytest.raises(ValueError,match='Remaining source'):
        tuner._propose(snapshot,8)
    assert checked==[] and tuner._header is None
    tuner.finish(successful=False)


def test_first_output_skips_prices_when_issue_frontier_is_short(monkeypatch):
    from torchgwas.sumstats import DenseWriteProgress
    cfg=config()
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',
        reduction=None)
    partition=dict(id='0',device='cuda:0',variant_range=[0,20],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});bound=tuner.for_partition('cuda:0',[0,20],[0,3])
    assert bound(0,20,8)==4
    monkeypatch.setattr(tuner,'_check',lambda:pytest.fail('short job must not validate prices'))
    tuner.output_written(DenseWriteProgress(0,4,(0,3),48,4,time.perf_counter(),'test'))
    state=tuner.run.snapshot()
    assert state['stop_reason']=='insufficient_unissued_work'
    assert state['planning']['steps']==[]
    assert state['current_chunk_size']==4
    tuner.finish(successful=False)


def test_first_output_skips_model_when_large_job_exceeds_extrapolation(monkeypatch):
    from torchgwas.sumstats import DenseWriteProgress
    cfg=config();cfg['forecast_options']['max_extrapolation']=2.
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',
        reduction=None)
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});bound=tuner.for_partition('cuda:0',[0,257],[0,3])
    assert bound(0,257,8)==4
    monkeypatch.setattr(tuner,'_check',lambda:pytest.fail('ineligible extrapolation must not validate prices'))
    tuner.output_written(DenseWriteProgress(0,4,(0,3),48,4,time.perf_counter(),'test'))
    state=tuner.run.snapshot()
    assert state['stop_reason']=='extrapolation_limit' and state['planning']['steps']==[]
    assert tuner.audit['productive_extrapolation_gate']['largest_remaining_to_horizon']==pytest.approx(253/24)
    assert bound(4,257,8)==4
    tuner.finish(successful=False)


@pytest.mark.parametrize('scenarios,message',[
    ({'empty':'empty'},'empty-only'),
    ({'empty':'empty','dense':'dense'},'partial output')])
def test_productive_output_rejects_refuted_empty_scenario_before_model(monkeypatch,scenarios,message):
    from torchgwas.sumstats_indexed import IndexedChunkWrite,IndexedOutputPartition
    cfg=config()
    owner=SimpleNamespace(config=dict(initial_chunks=cfg,joint=dict(
        occupancy_scenarios=scenarios)),input_path='never-opened.pgen',
        reduction='significant')
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{});bound=tuner.for_partition('cuda:0',[0,257],[0,3])
    assert bound(0,257,8)==4
    monkeypatch.setattr(tuner,'_check',lambda:pytest.fail('price validation should not start'))
    now=time.perf_counter()
    event=IndexedChunkWrite(0,4,'significant',1,100,'part.npz',now-.001,now,True,
        IndexedOutputPartition('cuda:0',(0,257),(0,3)),(0,4))
    tuner.output_written(event)
    state=tuner.run.snapshot()
    assert state['stop_reason']=='planning_error'
    assert state['current_chunk_size']==4
    assert message in state['decisions'][-1]['message']
    assert tuner.output_sample.snapshot()['counts']['partial']==1
    tuner.finish(successful=False)


def test_unbound_dense_output_stops_optional_jit_before_model(monkeypatch):
    from torchgwas.sumstats import DenseWriteProgress
    cfg=config()
    owner=SimpleNamespace(config=dict(initial_chunks=cfg),input_path='never-opened.pgen',
        reduction=None)
    partition=dict(id='0',device='cuda:0',variant_range=[0,257],trait_range=[0,3])
    start=dict(partitions=[partition],chunk_sizes=[4,8],initial_size=4,
        input_file_identity={},memory=dict(retained_index_bases_bytes=0))
    tuner=PublicInitialChunkTuning(owner,start,{},{})
    assert tuner.for_partition('cuda:0',[0,257],[0,3])(0,257,8)==4
    monkeypatch.setattr(tuner,'_check',
        lambda:pytest.fail('Unbound output must stop before model validation'))
    tuner.output_written(DenseWriteProgress(
        0,4,(0,3),48,4,time.perf_counter(),'test','cuda:1'))
    state=tuner.run.snapshot()
    assert state['stop_reason']=='planning_error'
    assert state['current_chunk_size']==4
    assert 'output boundary' in state['decisions'][-1]['message']
    assert tuner.boundary.snapshot()['invalid_events']==1
    tuner.finish(successful=False)
