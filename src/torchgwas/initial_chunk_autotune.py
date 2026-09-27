"""Deferred calculator decisions inside the public, productive scan lifecycle.

The starting layout is explicit and memory-admitted. Only future genotype
chunks change here; phenotype/device reassignment remains a separate move.
"""
from contextvars import ContextVar
from copy import deepcopy
from functools import wraps
import itertools
import math
import threading
import time
import weakref

from .adaptive_chunks import aligned_chunk_sizes,InitialChunkMeasurements
from .calibration_cache import _digest
from .planning_session import IncrementalPlanningBudget
from .productive_output_sample import ProductiveOutputSample
from .productive_boundary import ProductiveBoundaryProgress
from .productive_planning_cost import ProductivePlanningCostHistory,validate_planning_cost_history
from .productive_run import ProductiveTuningRun
from .productive_occupancy import survivor_bins,validate_occupancy_scenarios
from .productive_source_stage import ProductiveSourceStage,validate_source_stage_config


_ACTIVE=ContextVar('torchgwas_productive_api_runs',default=None)


def productive_api_lifecycle(function):
    """Close optional planning state on every public API failure path."""
    @wraps(function)
    def run(*args,**kwargs):
        active=[];token=_ACTIVE.set(active)
        try:return function(*args,**kwargs)
        finally:
            try:
                for owner in active:owner.finish(successful=False)
            finally:_ACTIVE.reset(token)
    return run


def validate_initial_chunk_config(config, bounds, contexts, reduction, joint):
    required={'context','chunk_size','partition_axis','trait_block','window_markers',
              'budget','cost_forecasts','forecast_options'}
    optional={'structural_cache_dir','resident_copy_refresh','binding_digest_cache_dir','reuse_binding_digests',
              'stage_observations','planning_cost_history','capacity_scenarios','background_planning',
              'source_staging','staged_screen'}
    if not isinstance(config,dict) or set(config)-optional!=required:
        raise ValueError('initial_chunks requires an explicit starting layout, horizons, budget and forecast assumptions')
    if type(config.get('reuse_binding_digests',True)) is not bool:
        raise ValueError('reuse_binding_digests must be boolean')
    if type(config.get('background_planning',False)) is not bool:
        raise ValueError('background_planning must be boolean')
    if 'source_staging' in config:
        validate_source_stage_config(config['source_staging'],config['chunk_size'])
    if 'structural_cache_dir' in config and (not isinstance(config['structural_cache_dir'],str)
            or not config['structural_cache_dir'].strip()):
        raise ValueError('structural_cache_dir must be a nonempty string path')
    if 'binding_digest_cache_dir' in config and (not isinstance(config['binding_digest_cache_dir'],str)
            or not config['binding_digest_cache_dir'].strip()):
        raise ValueError('binding_digest_cache_dir must be a nonempty string path')
    if 'resident_copy_refresh' in config:
        from .resident_copy_refresh import validate_refresh_config
        validate_refresh_config(config['resident_copy_refresh'])
    if 'planning_cost_history' in config:
        validate_planning_cost_history(config['planning_cost_history'])
    if 'stage_observations' in config:
        stage=config['stage_observations']
        if not isinstance(stage,dict) or set(stage)!={'max_chunks_per_device','warmup_chunks','stride',
                                                     'max_window_seconds','cuda_events',
                                                     'measurement_reserve_seconds'}:
            raise ValueError('Explicit bounded productive stage observation settings required')
        if type(stage['max_chunks_per_device']) is not int or not 1<=stage['max_chunks_per_device']<=32:
            raise ValueError('At most 32 productive stage observations per GPU required')
        if (isinstance(stage['measurement_reserve_seconds'],bool) or
                not isinstance(stage['measurement_reserve_seconds'],(int,float)) or
                not math.isfinite(stage['measurement_reserve_seconds']) or
                stage['measurement_reserve_seconds']<=0):
            raise ValueError('Positive declared stage measurement cost reserve required')
        InitialChunkMeasurements(['cuda:0'],**{key:value for key,value in stage.items()
            if key!='measurement_reserve_seconds'})
    sizes=aligned_chunk_sizes(bounds['chunks'])
    if type(config['chunk_size']) is not int or config['chunk_size'] not in sizes:
        raise ValueError('Initial chunk must belong to the admitted chunk bounds')
    if 'staged_screen' in config:
        screen=config['staged_screen']
        required_screen={'chunk_sizes','occupancy_scenario','max_partitions',
            'max_unique_records','max_chunks_per_partition','max_cpu_seconds',
            'max_wall_seconds','max_rebases'}
        if ('source_staging' not in config or not isinstance(screen,dict) or
                set(screen)!=required_screen or
                not isinstance(screen['chunk_sizes'],list) or
                not 0<len(screen['chunk_sizes'])<=4 or
                len(set(screen['chunk_sizes']))!=len(screen['chunk_sizes']) or
                any(type(size) is not int or size not in sizes
                    for size in screen['chunk_sizes']) or
                config['chunk_size'] not in screen['chunk_sizes']):
            raise ValueError('Staged screen requires source staging and bounded admitted chunks')
        for key in ('max_partitions','max_unique_records','max_chunks_per_partition'):
            if type(screen[key]) is not int or screen[key]<1:
                raise ValueError('Positive staged screen '+key+' required')
        if screen['max_partitions']>16:
            raise ValueError('At most 16 staged screen partitions supported')
        for key in ('max_cpu_seconds','max_wall_seconds'):
            value=screen[key]
            if (isinstance(value,bool) or not isinstance(value,(int,float)) or
                    not math.isfinite(value) or value<=0):
                raise ValueError('Positive finite staged screen '+key+' required')
        if type(screen['max_rebases']) is not int or not 0<=screen['max_rebases']<=2:
            raise ValueError('At most two staged screen rebases supported')
        if ((reduction is None and screen['occupancy_scenario'] is not None) or
                (reduction is not None and screen['occupancy_scenario'] is None)):
            raise ValueError('Staged screen occupancy must match output mode')
    if len([c for c in contexts if c.get('name')==config['context']])!=1:
        raise ValueError('Initial context must identify one calibrated context')
    axis=config['partition_axis'];width=config['trait_block']
    if axis not in ('trait','variant') or (reduction=='jagwas' and axis!='variant') or (reduction=='significant' and axis!='trait'):
        raise ValueError('Initial partition axis differs from output mode; JAGWAS never partitions phenotypes')
    if axis=='variant' and width is not None:
        raise ValueError('Variant partitioning cannot set a phenotype block')
    if axis=='trait' and (type(width) is not int or width not in bounds['trait_blocks']):
        raise ValueError('Initial phenotype block must belong to the admitted tile bounds')
    horizons=config['window_markers']
    if (not isinstance(horizons,list) or len(horizons)!=3 or
            any(type(h) is not int or h<1 or h%sizes[-1] for h in horizons) or
            horizons!=sorted(set(horizons)) or horizons[-1]//sizes[0]>32):
        raise ValueError('Three increasing horizons aligned to chunk capacity, at most 32 smallest chunks, required')
    budget=config['budget']
    if not isinstance(budget,dict) or set(budget)!={'max_steps','max_cpu_seconds','max_window_seconds'}:
        raise ValueError('Explicit productive step, CPU and wall-window budgets required')
    IncrementalPlanningBudget(**budget)
    if 'planning_cost_history' in config and budget['max_steps']>32:
        raise ValueError('Planning cost history requires at most 32 productive steps')
    costs=config['cost_forecasts']
    if not isinstance(costs,dict) or set(costs)!={'remaining_seconds','expected_cpu_seconds','expected_wall_seconds',
            'switching_seconds','publication_seconds','reserve_seconds'}:
        raise ValueError('Explicit remaining horizon and all productive cost forecasts required')
    for key,value in costs.items():
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (
                key in ('remaining_seconds','expected_cpu_seconds','expected_wall_seconds') and value==0):
            raise ValueError('Invalid productive cost forecast: '+key)
    forecast=config['forecast_options']
    if not isinstance(forecast,dict) or set(forecast)!={'boundary_adjustments','relative_model_error',
            'max_slope_change','max_extrapolation','assumptions'}:
        raise ValueError('Explicit bounded forecast options required')
    from .window_forecast import _interval,_number
    if not isinstance(forecast['boundary_adjustments'],dict) or set(forecast['boundary_adjustments'])!={'baseline','candidate'}:
        raise ValueError('Both continuation boundary scenarios required')
    for key,value in forecast['boundary_adjustments'].items():_interval(value,key)
    for key in ('relative_model_error','max_slope_change'):
        if _number(forecast[key],key,True)>=1:raise ValueError(key+' must be below one')
    if _number(forecast['max_extrapolation'],'max_extrapolation')<1:raise ValueError('Invalid extrapolation limit')
    if (not isinstance(forecast['assumptions'],dict) or set(forecast['assumptions'])!=
            {'source_work','output_occupancy','resource_capacity','partition_balance'} or
            any(not isinstance(v,str) or not v.strip() for v in forecast['assumptions'].values())):
        raise ValueError('Explicit source, output, resource and balance assumptions required')
    hosts=joint['host_scenarios']
    occupancy=joint.get('occupancy_scenarios',{'full':None})
    capacities=config.get('capacity_scenarios',{'nominal':dict.fromkeys(('cpu','dram','input','output'),1.)})
    if (not isinstance(capacities,dict) or not capacities or
            any(not isinstance(name,str) or not name.strip() or not isinstance(row,dict) or
                set(row)!={'cpu','dram','input','output'} or
                any(isinstance(value,bool) or not isinstance(value,(int,float)) or
                    not math.isfinite(value) or not 0<value<=1 for value in row.values())
                for name,row in capacities.items()) or
            not any(all(value==1 for value in row.values()) for row in capacities.values())):
        raise ValueError('Explicit bounded capacity scenarios require one nominal case and positive fractions at most one')
    if (not isinstance(hosts,dict) or not hosts or not isinstance(occupancy,dict) or not occupancy
            or len(hosts)*len(occupancy)*len(capacities)>4):
        raise ValueError('One to four explicit host/output/capacity scenarios required for deferred tuning')
    if reduction is not None:
        validate_occupancy_scenarios(occupancy)
    for host in hosts.values():
        if not isinstance(host,dict) or set(host)!={'host_serial_fraction','host_serial_policy'}:
            raise ValueError('Explicit host sharing and serial policy required')
        fraction=host['host_serial_fraction']
        if _number(fraction,'host_serial_fraction',True)>1 or host['host_serial_policy'] not in ('fluid','held-first','held-last'):
            raise ValueError('Invalid host sharing scenario')


def prepare_public_initial_chunks(owner,workload,output,basis,started,*,_prepared_header=None):
    """Admit one starting layout; perform no runtime search or payload census."""
    import os
    import torch
    from .adaptive_start import prepare_adaptive_start
    from .api import _available_host_bytes
    from .detailed_calibration import execution_context,validate_detailed_profile
    from .price_binding import validate_price_bindings
    config=owner.config['initial_chunks'];joint=owner.config['joint']
    refresh=getattr(owner,'refresh',None)
    validation={} if refresh is None else dict(deferred_bindings=refresh.deferred_bindings)
    probe_bytes=0 if refresh is None else refresh.host_bytes
    context=next(c for c in owner.profile['contexts'] if c['name']==config['context'])
    if owner.reduction is None and output['block_bytes'] is None:
        raise ValueError('Deferred dense output requires an explicit writer block size')
    # Layout identity is part of the profile, not an implied relabeling of its
    # measured finish service. Preparation must preserve that fixed identity.
    for profile in context['profiles'].values():
        expected=True if owner.reduction=='jagwas' else output['store_beta']
        if profile.get('return_beta') is not expected:
            raise ValueError('Deferred profile must explicitly bind the requested beta layout')
    status=validate_price_bindings(owner.profile,**validation)['status']
    if status!=('declared_targets_verified' if refresh is None else 'pending_refresh'):
        raise ValueError('Deferred tuning requires declared immutable component price bindings')
    devices=context['devices'];initial_allocated={d:torch.cuda.memory_allocated(d) for d in devices}
    available={d:int(torch.cuda.mem_get_info(d)[0]+torch.cuda.memory_reserved(d)-initial_allocated[d]) for d in devices}
    host_available=_available_host_bytes()
    prepared_headers=[];prepared_indexes=[];prepared_caches=[]
    start=prepare_adaptive_start(workload,context,chunk_sizes=owner.config['bounds']['chunks'],
        initial_size=config['chunk_size'],partition_axis=config['partition_axis'],trait_block=config['trait_block'],
        reduction=owner.reduction,output=output,cpu_workers=joint['cpu_workers'],
        host_memory_bytes=min(joint['host_memory_bytes'],host_available)-probe_bytes,
        device_memory_bytes={d:min(joint['device_memory_bytes'][d],available[d]) for d in devices},
        host_reserve_bytes=joint['host_reserve_bytes'],device_reserve_bytes=joint['device_reserve_bytes'],
        device_memory_profiles={d:joint['device_memory_profiles'][d] for d in devices},
        significance_threshold=owner.significance_threshold,
        max_tiles=owner.config['bounds'].get('max_candidate_tiles',10000),
        max_census_chunks=owner.config['bounds'].get('max_census_chunks',1000000),
        _header_receiver=prepared_headers.append,_prepared_header=_prepared_header,
        _index_receiver=prepared_indexes.append,_compact_memory=True,
        _compact_cache_dir=config.get('structural_cache_dir'),
        _compact_cache_receiver=prepared_caches.append,
        _defer_compact_bases=True)
    if len(prepared_headers)!=1 or prepared_headers[0][0]!=start['input_file_identity']:
        raise ValueError('Adaptive admission did not retain one bound source header')
    if (len(prepared_indexes)>1 or (prepared_indexes and
            (prepared_indexes[0][0]!=start['input_file_identity'] or
             prepared_indexes[0][1] is not prepared_headers[0][1]))):
        raise ValueError('Adaptive admission returned an incompatible source index')
    retained_bases_bytes=int(prepared_indexes[0][2].nbytes) if prepared_indexes else 0
    admission_cache=prepared_caches[0] if prepared_caches else None
    retained_cache_bytes=(0 if admission_cache is None else
        admission_cache.owned_extra_bytes(retained_bases_bytes))
    start['memory']['host_bytes']+=retained_bases_bytes+retained_cache_bytes
    start['memory']['retained_index_bases_bytes']=retained_bases_bytes
    if admission_cache is not None:
        start['memory']['retained_admission_cache_bytes']=retained_cache_bytes
        start['structural_work']['admission_cache']=admission_cache.audit()
    if probe_bytes:
        start['memory']['host_bytes']+=probe_bytes
        start['memory']['resident_copy_refresh_bytes']=probe_bytes
        start['resource_budgets']['host_memory_bytes']+=probe_bytes
    if 'source_staging' in config:
        reserve=config['source_staging']['extra_host_reserve_bytes']
        start['memory']['host_bytes']+=reserve
        start['memory']['source_staging_reserve_bytes']=reserve
    if start['memory']['host_bytes']>min(joint['host_memory_bytes'],host_available):
        raise ValueError('Retained PGEN index exceeds starting host memory budget')
    if start['timing_geometry_missing']:
        raise ValueError('Deferred timing requires geometry for every admitted chunk and source tail')
    for key,value in start['required_environment'].items():
        if os.getenv(key)!=value:raise ValueError('Deferred execution environment differs: '+key)
    if torch.backends.cuda.matmul.allow_tf32 is not False:raise ValueError('Deferred tuning requires TF32 disabled')
    current=execution_context(owner.devices,input_path=owner.input_path,output_path=owner.output_path)
    validate_detailed_profile(owner.profile,current,**validation)
    if owner.reduction is not None:owner.price_evidence=owner._load_reduction_prices()
    audit=dict(status='deferred_initial_chunks',profile_sha256=_digest(owner.profile),
        config=deepcopy(owner.config),workload=deepcopy(workload),output=deepcopy(output),
        selected=dict(api_kwargs=deepcopy(start['api_kwargs']),memory=deepcopy(start['memory']),
            devices=list(devices),initial_chunk_size=config['chunk_size']),
        structural_work=deepcopy(start['structural_work']),candidates_evaluated=0,
        binding_seconds=owner.binding_seconds,planning_and_validation_seconds=time.perf_counter()-started,
        context_matches=True,selection_validated=False,runtime_prediction_validated=False,
        productive=None,scope='Explicit admitted starting layout; bounded calculator work begins after actual written output. Only future genotype chunk sizes may change. Source/model uncertainty and unbound prices remain explicit limitations.')
    # Memory admission consumed the complete fine-grid layout. Pricing later
    # needs only fixed template identity and the retained index, so release
    # tens of thousands of per-chunk Python records before the scan starts.
    for tile in start['candidate']['tiles']:
        encoded=tile['data']['encoded']
        tile['data']['encoded']={key:encoded[key] for key in ('path','samples','variant_range')}
        tile['data']['encoded']['input_identity']=start['input_file_identity']
    start['source_layout']=None
    # The metadata pass can itself take seconds on a large PGEN. Recheck live
    # capacity after releasing its temporary layout, before any scan starts.
    current_host=_available_host_bytes()
    if start['memory']['host_bytes']>min(joint['host_memory_bytes'],current_host+retained_bases_bytes+retained_cache_bytes):
        raise ValueError('Current host capacity changed during adaptive admission')
    current_devices={d:int(torch.cuda.mem_get_info(d)[0]+torch.cuda.memory_reserved(d)-initial_allocated[d])
        for d in devices}
    if any(start['memory']['device_bytes'][d]>min(joint['device_memory_bytes'][d],current_devices[d])
           for d in devices):
        raise ValueError('Current GPU capacity changed during adaptive admission')
    audit['live_capacity_after_admission']=dict(host_available_bytes=current_host,
        owned_retained_index_bases_bytes=retained_bases_bytes,
        owned_retained_admission_cache_bytes=retained_cache_bytes,
        device_available_bytes=current_devices)
    owner.productive=PublicInitialChunkTuning(owner,start,audit,initial_allocated,
        _prepared_header=prepared_headers[0],
        _prepared_index=prepared_indexes[0] if prepared_indexes else None,
        _admission_cache=admission_cache)
    active=_ACTIVE.get()
    if active is not None:active.append(owner.productive)
    return start['api_kwargs'],audit,basis


class PublicInitialChunkTuning:
    def __init__(self,owner,start,audit,initial_allocated,_prepared_header=None,
                 _prepared_index=None,_admission_cache=None):
        self.owner=owner;self.start=start;self.audit=audit;self.initial_allocated=initial_allocated
        self.config=deepcopy(owner.config['initial_chunks'])
        self.run=ProductiveTuningRun(start['partitions'],chunk_sizes=start['chunk_sizes'],initial=start['initial_size'],
            budget=IncrementalPlanningBudget(**self.config['budget']),structural_cache_dir=self.config.get('structural_cache_dir'))
        self.output_sample=ProductiveOutputSample(start['partitions'],getattr(owner,'reduction',None))
        self.boundary=ProductiveBoundaryProgress(start['partitions'],getattr(owner,'reduction',None))
        stage_config=self.config.get('stage_observations')
        self.stage_sample=(None if stage_config is None else InitialChunkMeasurements(
            list(dict.fromkeys(p['device'] for p in start['partitions'])),
            **{key:value for key,value in stage_config.items() if key!='measurement_reserve_seconds'}))
        self._stage_scans={}
        self._lock=threading.RLock();self._finished=False;self._finishing=False
        self._header=None;self._attempts=[]
        self._prepared_header=_prepared_header;self._prepared_index=_prepared_index
        self._admission_cache=_admission_cache
        self._retained_bases_bytes=int(start.get('memory',{}).get('retained_index_bases_bytes',0))
        self._retained_cache_bytes=int(start.get('memory',{}).get('retained_admission_cache_bytes',0))
        self.refresh=getattr(owner,'refresh',None)
        self._planning_done=False
        self._planning_cost_history=None
        self._binding_digests=None;self._binding_digest_initialized=False;self._binding_digest_error=None
        self._validation=dict(checks=0,wall_seconds=0.,cpu_seconds=0.)
        self._alternatives=iter(size for size in start['chunk_sizes'] if size!=start['initial_size'])
        self._planner_worker=None
        self._staged_screen_request=None
        self._staged_screen_worker=None
        self._staged_screen_evidence=None
        self._source_stage=None
        if 'source_staging' in self.config:
            def stage_header():
                from .pgen_work_bounds import PgenHeaderWork
                options=dict(max_cached_signatures=self.config['source_staging'].get(
                    'max_cached_signatures',1024),max_cached_bounds=self.config[
                    'source_staging'].get('max_cached_bounds',64))
                if _prepared_header is not None:options['_prepared_header']=_prepared_header
                if _prepared_index is not None:options['_prepared_index']=_prepared_index
                return PgenHeaderWork(owner.input_path,**options)
            self._source_stage=ProductiveSourceStage(stage_header,
                start['input_file_identity'],start['initial_size'],
                self.config['source_staging'],
                on_step_observation=self._source_writer_observation,
                on_complete=self._staged_source_complete)
        self._retry_size=None
        self._live_writers=weakref.WeakValueDictionary()
        self._writer_devices={}
        self._indexed_result_queue_ref=None
        self._indexed_result_queue_spec=None
        self._indexed_writer_ref=None
        if 'staged_screen' in self.config:
            self.register_native_staged_chunk_screen(self.config['staged_screen'])

    def register_native_staged_chunk_screen(self,config):
        """Register measured fixed-partition chunks without cold PGEN work."""
        from .productive_staged_candidates import fixed_partition_chunk_candidates
        context=next(row for row in self.owner.profile['contexts']
                     if row['name']==self.start['context'])
        shared=context['shared_capacities']
        options=dict(shared_source_capacities={name:shared[name]
                     for name in ('cpu','dram','input')},
            occupancy_scenario=config['occupancy_scenario'],
            max_candidates=len(config['chunk_sizes']),
            max_partitions=config['max_partitions'],
            max_unique_records=config['max_unique_records'],
            max_chunks_per_partition=config['max_chunks_per_partition'],
            max_cpu_seconds=config['max_cpu_seconds'],
            max_wall_seconds=config['max_wall_seconds'],
            max_rebases=config['max_rebases'])
        rank=self.start['candidate']['tiles'][0]['data']['covariates']
        output=self.start['candidate']['output']
        def factory(frontier):
            active=next(row for row in self.owner.profile['contexts']
                        if row['name']==self.start['context'])
            writer_options=None
            if self.owner.reduction is None:
                from .sumstats import _SYNC_FILE_RANGE
                with self._lock:
                    writers=list(self._live_writers.values())
                if not writers:
                    raise ValueError('Observed native dense writer required for staged candidates')
                settings=[]
                for writer in writers:
                    settings.append(dict(block_bytes=writer.block_bytes,
                        queue_depth=writer.queue_depth,
                        borrow_chunks=writer.borrow_chunks,fsync=writer.fsync,
                        writeback_bytes=writer.writeback_bytes,
                        sync_file_range=_SYNC_FILE_RANGE is not None,
                        store_variant_df=writer.store_variant_df))
                if any(row!=settings[0] for row in settings[1:]):
                    raise ValueError('Dense writers have different staged settings')
                writer_options=settings[0]
            prices=(None if self.owner.reduction is None else
                    self.owner.price_evidence['record']['value'])
            return fixed_partition_chunk_candidates(frontier,active,
                chunk_sizes=config['chunk_sizes'],
                partition_axis=self.config['partition_axis'],
                covariate_rank=rank,output=output,
                writer_options=writer_options,reduction_prices=prices,
                significance_threshold=(self.owner.significance_threshold
                    if self.owner.reduction=='significant' else None))
        self.register_staged_screen(factory,options)

    def register_staged_screen(self, candidate_factory, options):
        """Opt in to an evidence-only screen after the first-chunk source stage.

        The factory receives the exact unissued frontier and returns already
        priced candidates. Registration does no PGEN or calculator work. This
        experimental hook never applies a layout, even for a current screen.
        """
        required={'shared_source_capacities','occupancy_scenario',
                  'max_candidates','max_partitions','max_unique_records',
                  'max_chunks_per_partition','max_cpu_seconds','max_wall_seconds'}
        if (self._source_stage is None or not callable(candidate_factory) or
                not isinstance(options,dict) or
                set(options) not in (required,required|{'max_rebases'})):
            raise ValueError('Explicit bounded staged screen registration required')
        if (type(options.get('max_rebases',0)) is not int or
                not 0<=options.get('max_rebases',0)<=2):
            raise ValueError('At most two staged screen rebases required')
        for key in ('max_candidates','max_partitions','max_unique_records',
                    'max_chunks_per_partition'):
            if type(options[key]) is not int or options[key]<1:
                raise ValueError('Positive staged screen '+key+' required')
        for key in ('max_cpu_seconds','max_wall_seconds'):
            value=options[key]
            if (isinstance(value,bool) or not isinstance(value,(int,float)) or
                    not math.isfinite(value) or value<=0):
                raise ValueError('Positive finite staged screen '+key+' required')
        capacities=options['shared_source_capacities']
        if (not isinstance(capacities,dict) or set(capacities)!={'cpu','dram','input'} or
                any(isinstance(value,bool) or not isinstance(value,(int,float)) or
                    not math.isfinite(value) or value<=0 for value in capacities.values())):
            raise ValueError('Positive staged screen shared capacities required')
        with self._lock:
            if (self._finished or self._finishing or
                    self._staged_screen_request is not None or
                    self.run.revision_token()['written_events']):
                raise ValueError('Staged screen must register once before first output')
            self._staged_screen_request=(candidate_factory,deepcopy(options))

    def _staged_source_complete(self):
        with self._lock:
            if (self._finished or self._finishing or
                    self._staged_screen_request is None or
                    self._staged_screen_worker is not None):
                return
            try:
                self._staged_screen_worker=threading.Thread(
                    target=self._screen_staged_evidence,
                    name='torchgwas-staged-screen',daemon=True)
                self._staged_screen_worker.start()
            except (OSError,RuntimeError) as error:
                self._staged_screen_evidence=dict(status='worker_start_error',
                    error=type(error).__name__+': '+str(error),
                    prediction_complete=False,selection_validated=False)

    def _screen_staged_evidence(self):
        started=time.perf_counter();cpu_started=time.thread_time()
        before=None;attempts=[]
        try:
            from .layout_frontier import unissued_frontier
            from .productive_staged_screen import productive_staged_partial_screen
            from .productive_staged_source_binding import audit_staged_source_price_binding
            from .productive_staged_work_binding import audit_staged_work_price_binding
            with self._lock:factory,options=self._staged_screen_request
            screen_options={key:value for key,value in options.items()
                            if key!='max_rebases'}
            for attempt in range(options.get('max_rebases',0)+1):
                with self._lock:
                    snapshot=self.run.snapshot()
                    before={key:snapshot[key] for key in (
                        'issued_revision','written_events','current_chunk_size',
                        'stop_reason','finished')}
                    bound=(self.boundary.bind(snapshot)
                           if self.boundary.events==snapshot['written_events']
                           and snapshot['written_events'] else None)
                if attempt and not any(
                        p['cursor']<p['variant_range'][1]
                        for p in snapshot['partitions']):
                    evidence['rebase_stop_reason']='source_fully_issued'
                    break
                if bound is None or bound.get('valid') is not True:
                    raise ValueError('A bound written-output checkpoint is required')
                if self.refresh is not None and self.refresh.pending:
                    raise ValueError('Pending refresh cannot bind staged screen prices')
                self._check()
                profile_before=_digest(self.owner.profile)
                partitions=self.start['partitions']
                frontier=unissued_frontier(snapshot,
                    source_identity=self.start['input_file_identity'],
                    reduction=self.owner.reduction,
                    total_traits=self.audit['workload']['traits'],
                    job_variant_range=[min(p['variant_range'][0] for p in partitions),
                                       max(p['variant_range'][1] for p in partitions)])
                candidates=factory(frontier)
                source_binding=audit_staged_source_price_binding(
                    self.owner.profile,self.start['context'],candidates,
                    options['shared_source_capacities'],profile_sha256=profile_before)
                work_binding=audit_staged_work_price_binding(
                    self.owner.profile,self.start['context'],candidates,
                    reduction_price_evidence=(None if self.owner.reduction is None
                        else self.owner.price_evidence),profile_sha256=profile_before)
                remaining=dict(screen_options)
                remaining['max_cpu_seconds']=max(1e-9,
                    options['max_cpu_seconds']-(time.thread_time()-cpu_started))
                remaining['max_wall_seconds']=max(1e-9,
                    options['max_wall_seconds']-(time.perf_counter()-started))
                screen=productive_staged_partial_screen(frontier,self._source_stage,
                    candidates,output_boundary=bound,**remaining)
                self._check()
                if self.refresh is not None and self.refresh.pending:
                    raise ValueError('Refresh became pending during staged screening')
                profile_after=_digest(self.owner.profile)
                after=self.run.revision_token()
                status=('stale_binding' if profile_before!=profile_after else
                        'stale_frontier' if before!=after else
                        'screen_budget' if screen['stop_reason']!='complete' else
                        'source_prices_unbound' if
                        source_binding['status']!='declared_source_prices_verified' else
                        'work_prices_unbound' if
                        work_binding['status']!='declared_work_prices_verified' else
                        'current_evidence')
                attempts.append(dict(issue_output_token=before,
                    current_issue_output_token=after,status=status,
                    screen_stop_reason=screen['stop_reason'],
                    screen_cpu_seconds=screen.get('screen_cpu_seconds'),
                    screen_wall_seconds=screen.get('screen_wall_seconds')))
                evidence=dict(status=status,
                    issue_output_token=before,current_issue_output_token=after,
                    profile_sha256_before=profile_before,
                    profile_sha256_after=profile_after,
                    source_price_binding=source_binding,
                    work_price_binding=work_binding,
                    screen=screen,prediction_complete=False,selection_validated=False,
                    scope='First-chunk background calculator evidence. Candidate prices are supplied by the caller; a current exact frontier is not a completion forecast or layout-switch authorization.')
                if status!='stale_frontier' or screen['stop_reason']!='complete':
                    break
                with self._lock:finishing=self._finishing or self._finished
                if finishing:
                    evidence['rebase_stop_reason']='job_finishing'
                    break
                if attempt==options.get('max_rebases',0):
                    evidence['rebase_stop_reason']='attempt_limit'
                    break
                if (time.thread_time()-cpu_started>=options['max_cpu_seconds'] or
                        time.perf_counter()-started>=options['max_wall_seconds']):
                    evidence['rebase_stop_reason']='budget'
                    break
                if not any(p['cursor']<p['variant_range'][1]
                           for p in self.run.snapshot()['partitions']):
                    evidence['rebase_stop_reason']='source_fully_issued'
                    break
        except Exception as error:
            evidence=dict(status='screen_error',error=type(error).__name__+': '+str(error),
                issue_output_token=before,prediction_complete=False,
                selection_validated=False)
        evidence['rebase_attempts']=attempts
        evidence['wall_seconds']=time.perf_counter()-started
        evidence['cpu_seconds']=time.thread_time()-cpu_started
        with self._lock:self._staged_screen_evidence=evidence

    def register_writer(self,writer,device):
        """Keep only weak, active native writer references for live queue reads."""
        from .sumstats import BinarySumstatsWriter
        if not isinstance(writer,BinarySumstatsWriter) or writer.on_write_progress is None:
            raise ValueError('An observed native dense writer is required')
        device=str(device)
        if device not in {row['device'] for row in self.start['partitions']}:
            raise ValueError('Dense writer device differs from admitted partitions')
        directory=str(writer.directory.resolve())
        with self._lock:
            if self._finished:
                raise ValueError('Cannot register a writer after productive finish')
            self._writer_devices={key:value for key,value in self._writer_devices.items()
                                  if key in self._live_writers}
            prior=self._live_writers.get(directory)
            if prior is not None and prior is not writer:
                raise ValueError('Two active dense writers share one directory')
            self._live_writers[directory]=writer
            self._writer_devices[directory]=device

    def capture_live_writer_queues(self,*,max_writers=16):
        """Bracket current dense writer queues at one common instant.

        Two passes make cross-writer byte bounds possible while callbacks and
        GPU work continue. Issue/output revisions must remain unchanged, but
        this still does not observe progress inside already-issued GPU work.
        """
        from .productive_dense_writer_queue_bracket import bracket_dense_writer_queues
        if type(max_writers) is not int or max_writers<1:
            raise ValueError('Positive live writer observation budget required')
        started=time.perf_counter()
        with self._lock:
            before=self.run.revision_token()
            registered=[(path,self._writer_devices[path],writer)
                        for path,writer in self._live_writers.items()]
        if len(registered)>max_writers:
            return dict(kind='torchgwas.live_dense_writer_queues.v1',
                observation_valid=False,reason='writer_budget',
                registered_writers=len(registered),max_writers=max_writers,
                issued_revision=before['issued_revision'],
                written_events=before['written_events'])
        observations=[]
        errors={}
        for pass_index in range(2):
            observed={}
            for path,device,writer in registered:
                try:observed[path]=dict(device=device,observation=writer.queue_snapshot())
                except Exception as error:errors[path]=type(error).__name__+': '+str(error)
            observations.append(observed)
            if pass_index==0:
                with self._lock:
                    middle=self.run.revision_token()
                    anchor=time.perf_counter()
        with self._lock:
            after=self.run.revision_token()
            same_writers=(len(self._live_writers)==len(registered) and
                all(self._live_writers.get(path) is writer and
                    self._writer_devices.get(path)==device
                    for path,device,writer in registered))
        finished=time.perf_counter()
        stable=before==middle==after and same_writers
        bracket=None
        if not errors and registered:
            try:bracket=bracket_dense_writer_queues(*observations,anchor)
            except ValueError as error:errors['bracket']=str(error)
        valid=stable and bracket is not None and not errors
        return dict(kind='torchgwas.live_dense_writer_queues.v1',
            observation_valid=valid,stable_revision=stable,
            issue_output_token=before,
            registered_writers=len(registered),max_writers=max_writers,
            capture_started_seconds=started,capture_anchor_seconds=anchor,
            capture_finished_seconds=finished,
            writers={path:dict(device=device,first=observations[0].get(path),
                               second=observations[1].get(path))
                     for path,device,_ in registered},
            bracket=bracket,errors=errors,
            prediction_complete=False,selection_validated=False,
            scope='Two live passes bracket accepted-minus-written bytes across all registered writers at one common time if issue/output revisions and writer identities stay stable. Producer/GPU progress inside issued work, unregistered or future writers, dirty writeback and final fsync are not synchronized. No completion bound or JIT switch authorization.')

    def register_indexed_result_queue(self,result_queue,bounds,devices,finished):
        # The shared queue is created before JAGWAS producer threads start.
        import queue
        if self.owner.reduction!='jagwas' or not isinstance(result_queue,queue.Queue):
            raise ValueError('JAGWAS shared result queue required')
        claimed={(str(device),tuple(span)) for span,device in zip(bounds,devices)}
        expected={(row['device'],tuple(row['variant_range']))
                  for row in self.start['partitions']}
        if (len(bounds)!=len(devices) or claimed!=expected or
                len(claimed)<2 or result_queue.maxsize<1):
            raise ValueError('JAGWAS result queue differs from admitted shards')
        with self._lock:
            if self._finished or self._indexed_result_queue_ref is not None:
                raise ValueError('JAGWAS queue registration must be unique and active')
            self._indexed_result_queue_ref=weakref.ref(result_queue)
            self._indexed_result_queue_spec=(tuple(map(tuple,bounds)),tuple(devices),finished)

    def register_indexed_writer(self,writer):
        from .indexed_writer_progress import IndexedWriterProgress
        if (self.owner.reduction not in ('jagwas','significant') or
                not isinstance(writer,IndexedWriterProgress) or
                writer.kind!=self.owner.reduction):
            raise ValueError('Observed indexed writer must match reduction mode')
        with self._lock:
            if self._finished or self._indexed_writer_ref is not None:
                raise ValueError('JAGWAS indexed writer registration must be unique and active')
            self._indexed_writer_ref=weakref.ref(writer)

    def capture_indexed_result_queue(self):
        from .productive_indexed_result_queue import snapshot_indexed_result_queue
        from .productive_indexed_queue_join import bind_jagwas_result_queue
        with self._lock:
            before=self.run.revision_token()
            reference=self._indexed_result_queue_ref
            spec=self._indexed_result_queue_spec
            current=None if reference is None else reference()
            writer_reference=self._indexed_writer_ref
            writer=None if writer_reference is None else writer_reference()
            bound=None
            if before['written_events'] and self.boundary.events==before['written_events']:
                issued=self.run.snapshot()
                if issued['prefix_complete']:
                    try:bound=self.boundary.bind(issued)
                    except ValueError:pass
        if current is None or spec is None:
            return dict(kind='torchgwas.live_indexed_result_queue.v1',
                observation_valid=False,reason='no_live_shared_queue',
                issue_output_token=before)
        error=None;observation=None
        writer_first=None if writer is None else writer.snapshot()
        try:observation=snapshot_indexed_result_queue(current,*spec)
        except Exception as caught:error=type(caught).__name__+': '+str(caught)
        writer_second=None if writer is None else writer.snapshot()
        with self._lock:
            after=self.run.revision_token()
            same=(self._indexed_result_queue_ref is reference and reference() is current
                  and self._indexed_writer_ref is writer_reference and
                  (writer_reference is None or writer_reference() is writer))
        stable=before==after and same
        active=None
        if (writer_first is not None and writer_second is not None and
                writer_first['active'] is not None and
                writer_first['active']==writer_second['active'] and
                writer_first['phase'] in ('selecting','emitting_part') and
                writer_second['phase'] in ('selecting','emitting_part')):
            active=writer_first['active']
        joined=None
        if stable and observation is not None and bound is not None and bound['valid']:
            try:joined=bind_jagwas_result_queue(bound,observation,active_writer=active)
            except ValueError as caught:error=type(caught).__name__+': '+str(caught)
        return dict(kind='torchgwas.live_indexed_result_queue.v1',
            observation_valid=stable and observation is not None and error is None,
            checkpoint_valid=stable and joined is not None and error is None,
            stable_revision=stable,issue_output_token=before,
            observation=observation,issued_output_boundary=bound,
            indexed_writer_bracket=dict(first=writer_first,second=writer_second,
                active_at_queue_anchor=active),
            joined_issued_queue=joined,error=error,
            prediction_complete=False,selection_validated=False,
            scope='Atomic bounded JAGWAS shared-result queue snapshot joined to an unchanged source-issue/indexed-completion frontier when available. A currently consumed part and upstream producer/GPU work remain conservative; no completion bound or JIT switch authorization.')

    def _price_live_writer_capture(self,capture):
        if not capture.get('observation_valid'):
            return dict(status='unavailable',reason='unstable_or_missing_writer_capture')
        try:
            from .productive_dense_writer_queue_service import price_bracketed_dense_writer_queues
            matched=[row for row in self.owner.profile['contexts']
                     if row['name']==self.start['context']]
            if len(matched)!=1:raise ValueError('One current writer profile context required')
            devices={row['device'] for row in capture['bracket']['writers'].values()}
            profiles={device:matched[0]['profiles'][device] for device in devices}
            return price_bracketed_dense_writer_queues(capture['bracket'],profiles)
        except (AttributeError,KeyError,TypeError,ValueError) as error:
            return dict(status='unavailable',error=type(error).__name__+': '+str(error))

    def capture_significant_indexed_writer(self):
        if self.owner.reduction!='significant':
            raise ValueError('Significant-pair writer observation required')
        with self._lock:
            before=self.run.revision_token()
            output=self.boundary.snapshot()
            reference=self._indexed_writer_ref
            writer=None if reference is None else reference()
        observed=None if writer is None else writer.snapshot()
        with self._lock:
            after=self.run.revision_token()
            same=self._indexed_writer_ref is reference and (
                reference is None or reference() is writer)
        return dict(kind='torchgwas.staged_significant_output_observation.v1',
            issue_output_token=before,output_observations=output,
            writer=observed,observation_valid=observed is not None and
                before==after and same,
            status=('indexed_writer_observed_host_selector_inflight_unobserved'
                if observed is not None and before==after and same else
                'indexed_writer_unavailable_host_selector_inflight_unobserved'),
            prediction_complete=False,selection_validated=False,
            scope='Atomic active indexed-writer range and publication phase at a stable source/output revision. Host selector, shared result queue and final durability service remain unpriced; no switch authorization.')

    def _source_writer_observation(self):
        # Compact output evidence after one background source-metadata step.
        if self.owner.reduction=='jagwas':
            return self.capture_indexed_result_queue()
        if self.owner.reduction=='significant':
            return self.capture_significant_indexed_writer()
        capture=self.capture_live_writer_queues()
        return dict(kind='torchgwas.staged_dense_writer_observation.v1',
            issue_output_token=capture.get('issue_output_token'),
            observation_valid=capture['observation_valid'],
            stable_revision=capture.get('stable_revision'),
            registered_writers=capture['registered_writers'],
            bracket=capture.get('bracket'),
            priced_service=self._price_live_writer_capture(capture),
            scope='Bounded evidence after one useful-output source step. No source/GPU or durable-output completion bound and no switch authorization.')

    def for_partition(self,device,variant_range,trait_range):
        matches=[p['id'] for p in self.start['partitions'] if p['device']==str(device)
            and p['variant_range']==list(variant_range) and p['trait_range']==list(trait_range)]
        if len(matches)!=1:raise ValueError('Public scan does not identify one admitted productive partition')
        return self.run.for_partition(matches[0])

    def stage_observer(self,device,variant_range,trait_range):
        """Measure the first admitted scan partition per GPU, at most once."""
        if self.stage_sample is None:return None
        device=str(device)
        matches=[p['id'] for p in self.start['partitions'] if p['device']==device
            and p['variant_range']==list(variant_range) and p['trait_range']==list(trait_range)]
        if len(matches)!=1:raise ValueError('Stage observation differs from admitted scan partition')
        with self._lock:
            prior=self._stage_scans.setdefault(device,matches[0])
            return self.stage_sample if prior==matches[0] else None

    def __call__(self,observation):
        raise ValueError('Bind the productive stage observer to a scan first')

    def for_scan(self,device,variant_range,n_traits,**_scan_settings):
        """Native variant-shard route for the same first-partition observer."""
        matches=[p for p in self.start['partitions'] if p['device']==str(device)
            and p['variant_range']==list(variant_range)
            and p['trait_range'][1]-p['trait_range'][0]==n_traits]
        if len(matches)!=1:raise ValueError('Stage scan does not identify one admitted partition')
        return self.stage_observer(device,variant_range,matches[0]['trait_range'])

    def _check(self):
        began=time.perf_counter();cpu=time.thread_time()
        try:self._check_current()
        finally:
            self._validation['checks']+=1
            self._validation['wall_seconds']+=time.perf_counter()-began
            self._validation['cpu_seconds']+=time.thread_time()-cpu

    def _cost_history(self):
        options=self.config.get('planning_cost_history')
        if options is None:return None
        if self._planning_cost_history is None:
            joint=self.owner.config['joint']
            dependencies=dict(protocol='torchgwas.productive_planning_cost.v1',
                input_identity=deepcopy(self.start['input_file_identity']),
                profile_sha256=_digest(self.owner.profile),
                workload=deepcopy(self.audit['workload']),output=deepcopy(self.audit['output']),
                reduction=self.owner.reduction,partitions=deepcopy(self.start['partitions']),
                chunk_sizes=list(self.start['chunk_sizes']),initial_chunk_size=self.start['initial_size'],
                window_markers=list(self.config['window_markers']),
                host_scenarios=deepcopy(joint['host_scenarios']),
                occupancy_scenarios=deepcopy(joint.get('occupancy_scenarios')),
                capacity_scenarios=deepcopy(self.config.get('capacity_scenarios')),
                forecast_options=deepcopy(self.config['forecast_options']),
                planning_budget=deepcopy(self.config['budget']),
                structural_cache_dir=self.config.get('structural_cache_dir'),
                binding_digest_cache_dir=self.config.get('binding_digest_cache_dir'),
                reuse_binding_digests=self.config.get('reuse_binding_digests',True),
                resident_copy_refresh=deepcopy(self.config.get('resident_copy_refresh')),
                stage_observations=deepcopy(self.config.get('stage_observations')))
            dependencies['background_planning']=self.config.get('background_planning',False)
            dependencies['source_staging']=deepcopy(self.config.get('source_staging'))
            self._planning_cost_history=ProductivePlanningCostHistory(options,dependencies,
                max_steps=self.config['budget']['max_steps'])
        return self._planning_cost_history

    def _digest_cache(self):
        # Initialization/lookup happens only inside a charged productive step.
        # Startup retains full-byte validation. Only file digests are reused.
        if not self._binding_digest_initialized:
            self._binding_digest_initialized=True
            if 'binding_digest_cache_dir' in self.config:
                directory=self.config['binding_digest_cache_dir']
            else:
                directory=(self.config.get('resident_copy_refresh',{}).get('cache_dir')
                    or self.config.get('structural_cache_dir')
                    or self.config.get('planning_cost_history',{}).get('cache_dir'))
            if directory is not None and self.config.get('reuse_binding_digests',True):
                from .binding_digests import BindingDigestCache
                try:self._binding_digests=BindingDigestCache(directory)
                except (OSError,ValueError) as error:self._binding_digest_error=type(error).__name__+': '+str(error)
        return self._binding_digests

    def _check_current(self):
        import torch
        from .api import _available_host_bytes
        from .analytical_plan_cache import input_identity
        from .detailed_calibration import execution_context,source_identity,validate_detailed_profile
        phases={};started=time.perf_counter()
        cache=self._digest_cache()
        context=execution_context(self.owner.devices,input_path=self.owner.input_path,output_path=self.owner.output_path,
            **({} if cache is None else dict(digest_cache=cache)))
        phases['execution_context']=time.perf_counter()-started;started=time.perf_counter()
        options={} if self.refresh is None else dict(deferred_bindings=self.refresh.deferred_bindings)
        if cache is not None:options['sources']=source_identity(digest_cache=cache)
        validate_detailed_profile(self.owner.profile,context,**options)
        phases['profile_and_prices']=time.perf_counter()-started;started=time.perf_counter()
        if cache is not None and not cache.unchanged():
            raise ValueError('Source or library identity changed during productive validation')
        phases['digest_recheck']=time.perf_counter()-started;started=time.perf_counter()
        if input_identity(self.owner.input_path)!=self.start['input_file_identity']:
            raise ValueError('Deferred source input changed')
        if self.owner.reduction is not None:self.owner.price_evidence=self.owner._load_reduction_prices()
        phases['input_and_reduction_prices']=time.perf_counter()-started;started=time.perf_counter()
        # Fixed rings were admitted at capacity. Reclaimable allocator cache and
        # this run's own allocations are available to its admitted envelope;
        # unrelated allocations present before this run remain excluded.
        for device,needed in self.start['memory']['device_bytes'].items():
            available=torch.cuda.mem_get_info(device)[0]+torch.cuda.memory_reserved(device)-self.initial_allocated[device]
            if needed>available:raise ValueError('Current GPU capacity no longer admits deferred shapes')
        if self.start['memory']['host_bytes']>_available_host_bytes()+self._retained_bases_bytes+self._retained_cache_bytes:
            raise ValueError('Current host capacity no longer admits deferred shapes')
        phases['live_memory']=time.perf_counter()-started
        rows=self._validation.setdefault('phase_rows',[])
        if len(rows)<8:rows.append(phases)

    def _propose(self,snapshot,size):
        from .pgen_work_bounds import PgenHeaderWork
        from .window_model import prepared_source_window,compare_prepared_windows
        from .productive_forecast import productive_window_proposal
        from .price_binding import validate_comparison_prices
        # Background planning can run while writer callbacks append progress.
        # Take one atomic observation of the held issue/output checkpoint.
        with self._lock:
            bound_output=None
            if self.boundary.events==snapshot['written_events'] and self.boundary.events:
                bound_output=self.boundary.bind(snapshot)
                if not bound_output['valid']:
                    raise ValueError('Productive output boundary cannot bind to issued source')
        live_writer_queues=(self.capture_live_writer_queues()
                            if self.owner.reduction is None else None)
        if self.refresh is not None and self.refresh.pending:
            raise ValueError('Pending refresh cannot authorize a runtime comparison')
        if self.output_sample.reduction is not None:
            with self._lock:
                self.output_sample.check_scenarios(self.owner.config['joint'].get(
                    'occupancy_scenarios',{'full':None}))
        pending=[p for p in snapshot['partitions'] if p['cursor']<p['variant_range'][1]]
        horizons=self.config['window_markers']
        if len(snapshot['partitions'])>8 or not pending or min(p['variant_range'][1]-p['cursor'] for p in pending)<horizons[-1]:
            raise ValueError('Remaining source does not support the bounded forecast horizons')
        phase_started=time.perf_counter();phases={}
        self._check()
        phases['initial_validation']=time.perf_counter()-phase_started
        phase_started=time.perf_counter()
        if self._header is None:
            if self._prepared_header is None:
                self._header=PgenHeaderWork(self.owner.input_path)
            elif self._prepared_index is None:
                self._header=PgenHeaderWork(self.owner.input_path,_prepared_header=self._prepared_header)
            else:
                self._header=PgenHeaderWork(self.owner.input_path,_prepared_header=self._prepared_header,
                    _prepared_index=self._prepared_index)
            self._prepared_header=None;self._prepared_index=None
        phases['header_setup']=time.perf_counter()-phase_started
        phase_started=time.perf_counter()
        k=self.audit['workload']['traits'];reduction=self.owner.reduction;candidate=self.start['candidate']
        live_writer_service=(None if live_writer_queues is None else
                             self._price_live_writer_capture(live_writer_queues))
        templates={p['id']:tile for p,tile in zip(self.start['partitions'],candidate['tiles'])}
        layouts=[]
        for count in horizons:
            pair=[]
            for chunk in (snapshot['current_chunk_size'],size):
                windows=[prepared_source_window(templates[p['id']],self._header,start=p['cursor'],stop=p['cursor']+count,
                    chunk_markers=chunk,issued_chunks=p['issued_chunks'],expected_input_identity=self.start['input_file_identity'],
                    max_chunks=32) for p in pending]
                pair.append(dict(windows=windows,partition_axis=self.config['partition_axis']))
            layouts.append(pair)
        phases['source_windows']=time.perf_counter()-phase_started
        phase_started=time.perf_counter()
        prices=None if reduction is None else self.owner.price_evidence['record']['value']
        if reduction=='jagwas':prices=prices['writer_prices']
        common=dict(total_traits=k,reduction=reduction,output=candidate['output'],shared_capacities=candidate['shared_capacities'],
            endpoint='upper',max_source_chunks=256,max_survivor_bins=256,prices=prices,
            model_identity=dict(source_sha256=self.owner.profile['source_sha256'],price_profile_sha256=_digest(self.owner.profile)))
        for key in ('shared_links','shared_storage_bytes_per_second'):
            if key in candidate:common[key]=candidate[key]
        if reduction=='significant':common['significance_threshold']=self.owner.significance_threshold
        proposals=[];audits=[];scenario_phases=[]
        occupancies=self.owner.config['joint'].get('occupancy_scenarios',{'full':None})
        # The source/occupancy ledger is independent of host sharing and
        # capacity fractions. Reuse it across those matched scenarios so a
        # clustered global pattern is walked at most once per horizon.
        occupancy_bins=({(name,count):survivor_bins(
            pending,count,fine_chunk=self.start['chunk_sizes'][0],
            reduction=reduction,total_traits=k,scenario=value)
            for name,value in occupancies.items() for count in horizons}
            if reduction is not None else {})
        phases['scenario_setup']=time.perf_counter()-phase_started
        capacities=self.config.get('capacity_scenarios',{'nominal':dict.fromkeys(('cpu','dram','input','output'),1.)})
        for (host_name,host),(occupancy_name,occupancy),(capacity_name,fractions) in itertools.product(
                self.owner.config['joint']['host_scenarios'].items(),occupancies.items(),capacities.items()):
            scenario_started=time.perf_counter()
            comparisons=[]
            scenario_common=dict(common,shared_capacities={resource:common['shared_capacities'][resource]*fraction
                for resource,fraction in fractions.items()})
            for pair,count in zip(layouts,horizons):
                evidence=None
                if reduction is not None:
                    evidence=dict(input_identity=self._header.input_identity,reduction=reduction,total_traits=k,
                        significance_threshold=self.owner.significance_threshold,
                        bins=occupancy_bins[(occupancy_name,count)])
                comparisons.append(compare_prepared_windows(*pair,survivor_evidence=evidence,**scenario_common,**host))
            comparison_seconds=time.perf_counter()-scenario_started
            scenario_started=time.perf_counter()
            evidence=validate_comparison_prices(self.owner.profile,comparisons)
            price_validation_seconds=time.perf_counter()-scenario_started
            scenario_started=time.perf_counter()
            proposal,audit=productive_window_proposal(snapshot,comparisons,chunk_sizes=self.start['chunk_sizes'],**self.config['forecast_options'])
            forecast_seconds=time.perf_counter()-scenario_started
            scenario_phases.append(dict(host_scenario=host_name,occupancy_scenario=occupancy_name,
                capacity_scenario=capacity_name,
                comparison_seconds=comparison_seconds,price_validation_seconds=price_validation_seconds,
                forecast_seconds=forecast_seconds))
            audit.update(host_scenario=host_name,occupancy_scenario=occupancy_name,
                capacity_scenario=capacity_name,shared_capacity_fractions=deepcopy(fractions),price_evidence=evidence)
            proposals.append(proposal);audits.append(audit)
        phase_started=time.perf_counter()
        self._check()  # Final original-age and current-resource check before apply.
        phases['final_validation']=time.perf_counter()-phase_started
        attempt=dict(issued_revision=snapshot['issued_revision'],candidate_size=size,scenarios=audits,
            phase_wall_seconds=phases,scenario_phase_wall_seconds=scenario_phases,
            header_bounds_cache=getattr(self._header,'bounds_cache_info',lambda:None)())
        attempt['output_boundary']=(bound_output if bound_output is not None else
            dict(valid=False,reason='output_progress_not_observed'))
        if live_writer_queues is not None:attempt['live_writer_queues']=live_writer_queues
        if live_writer_service is not None:attempt['live_writer_queue_service']=live_writer_service
        self._attempts.append(attempt)
        if any(a['forecast_status']!='stable_scenario' for a in audits):
            return dict(chunk_size=snapshot['current_chunk_size'],baseline_seconds=0.,candidate_seconds=0.)
        # Both layouts face the same declared host/output scenario. Crossing
        # the fastest baseline with the slowest candidate compares different
        # conditions and can discard a gain present in every matched pair.
        gains=[p['baseline_seconds']-p['candidate_seconds'] for p in proposals]
        worst=min(range(len(proposals)),key=gains.__getitem__)
        maximum_baseline=max(p['baseline_seconds'] for p in proposals)
        attempt['matched_scenarios']=dict(method='minimum_paired_gain',
            gain_floor_seconds=gains[worst],scenario_gains_seconds=gains,
            limiting_scenario_index=worst,maximum_baseline_seconds=maximum_baseline,
            scope='Minimum baseline lower minus candidate upper within each declared host/output scenario. '
                  'Conditional model savings, not a hardware-time guarantee.')
        # The budget checks the returned pair too, but the limiting-gain pair
        # need not have the longest baseline. Keep every scenario in horizon.
        if maximum_baseline>self.config['cost_forecasts']['remaining_seconds']:
            raise ValueError('A matched scenario exceeds the declared remaining horizon')
        return dict(chunk_size=size,baseline_seconds=proposals[worst]['baseline_seconds'],
            candidate_seconds=proposals[worst]['candidate_seconds'])

    def _refresh_step(self,snapshot):
        self._check()
        profile=self.refresh.advance(written_events=snapshot['written_events'],issued_revision=snapshot['issued_revision'],validate=self._check)
        if profile is not None:
            self.owner.profile=profile
            self.refresh.patch_tiles(self.start['candidate']['tiles'],self.start['context'])
            self.audit['active_profile_sha256']=_digest(profile)
        self._check()
        # A probe is charged planning work, never itself a chunk-size proposal.
        return dict(chunk_size=snapshot['current_chunk_size'],baseline_seconds=0.,candidate_seconds=0.)

    def _stage_measurement_reserve(self):
        if self.stage_sample is None:return 0.
        return max(self.config['stage_observations']['measurement_reserve_seconds'],
                   self.stage_sample.known_probe_wall_seconds())

    def _proposal_step(self,size,*,release_issue_frontier):
        costs=dict(self.config['cost_forecasts'])
        costs['reserve_seconds']+=self._stage_measurement_reserve()
        history_started=time.perf_counter()
        history=self._cost_history()
        if history is not None:
            costs=history.costs(costs)
            preparation=time.perf_counter()-history_started
            history.charge_preparation(preparation)
            costs['reserve_seconds']+=preparation
        observed=time.time()
        result=self.run.planning_step(lambda snapshot:self._propose(snapshot,size),
            release_issue_frontier=release_issue_frontier,**costs)
        if history is not None:history.observe(result,observed_unix_seconds=observed)
        with self._lock:
            if result.get('stale_frontier',False):
                # An asynchronous result consumed budget but did not test
                # this size at the current source/output revision. Retain it
                # for the next written event, ahead of later alternatives.
                self._retry_size=size
            elif not result.get('usable_for_decision',False):
                self._planning_done=True
        return result

    def _background_proposal(self,size):
        try:
            self._proposal_step(size,release_issue_frontier=True)
        except Exception as error:
            # Optional planning must not turn a completed scientific output
            # into a writer failure. Retain the error and stop future tuning.
            with self._lock:
                self.audit['background_planning_error']=dict(type=type(error).__name__,message=str(error))
                self._planning_done=True
                self.run.stop_planning('background_error')

    def output_written(self,event):
        # Writer threads never wait on another writer or GPU. This lock merely
        # serializes optional calculations and ends before any scan is joined.
        with self._lock:
            self.run.output_written(event)
            self.output_sample.observe(event)
            self.boundary.observe(event)
            if self._source_stage is not None:self._source_stage.output_written()
            if self._planning_done:return
            if self._planner_worker is not None and self._planner_worker.is_alive():return
            if self._source_stage is not None:
                # Staged first-chunk evidence uses the finite continuation,
                # not the old three-window/at-most-eight-partition gate. A
                # pending independent refresh may still spend one charged step.
                if self.refresh is not None and self.refresh.pending:
                    costs=dict(self.config['cost_forecasts'],
                        expected_cpu_seconds=self.refresh.config['expected_cpu_seconds'],
                        expected_wall_seconds=self.refresh.config['expected_wall_seconds'])
                    costs['reserve_seconds']+=self._stage_measurement_reserve()
                    result=self.run.planning_step(self._refresh_step,**costs)
                    if not result.get('usable_for_decision',False):
                        self._planning_done=True;self.refresh.close()
                # The old short-window switch never runs in staged mode. Keep
                # the exact issue prefix until the normal range/window budget
                # ends it, so later chunks still have a bindable frontier.
                return
            # Admission is already complete; if the real issue frontier has
            # passed the largest horizon, no price validation, cache lookup or
            # source-window graph can repay a switch in this job.
            frontier=self.run.snapshot()['partitions']
            pending=[p for p in frontier if p['cursor']<p['variant_range'][1]]
            if len(frontier)>8 or not pending or min(p['variant_range'][1]-p['cursor'] for p in pending)<self.config['window_markers'][-1]:
                reason=('partition_budget' if len(frontier)>8 else
                        'source_fully_issued' if not pending else 'insufficient_unissued_work')
                self.run.stop_planning(reason)
                self._planning_done=True
                if self.refresh is not None:self.refresh.close()
                return
            if self.refresh is not None and self.refresh.pending:
                costs=dict(self.config['cost_forecasts'],
                    expected_cpu_seconds=self.refresh.config['expected_cpu_seconds'],
                    expected_wall_seconds=self.refresh.config['expected_wall_seconds'])
                costs['reserve_seconds']+=self._stage_measurement_reserve()
                result=self.run.planning_step(self._refresh_step,**costs)
                if not result.get('usable_for_decision',False):
                    self._planning_done=True;self.refresh.close()
                return
            # Every modeled partition spans the same largest horizon. The
            # forecast's upper balance envelope is at least the longest
            # remaining partition divided by that horizon. If it already
            # exceeds the declared cap, no graph/price check can make this
            # callback eligible for a switch. Complete scheduled component
            # refresh first so this job can still publish reusable evidence.
            horizon=self.config['window_markers'][-1]
            ratio=max((p['variant_range'][1]-p['cursor'])/horizon for p in pending)
            if ratio>self.config['forecast_options']['max_extrapolation']:
                self.audit['productive_extrapolation_gate']=dict(
                    largest_remaining_to_horizon=ratio,
                    max_extrapolation=self.config['forecast_options']['max_extrapolation'],
                    largest_horizon_markers=horizon,
                    scope='Exact issue-frontier extent exceeds the declared finite-window extrapolation cap before model construction.')
                self.run.stop_planning('extrapolation_limit')
                self._planning_done=True
                return
            size=self._retry_size
            if size is None:size=next(self._alternatives,None)
            self._retry_size=None
            if size is None:
                self._planning_done=True
                return
            if self.config.get('background_planning',False):
                self._planner_worker=threading.Thread(target=self._background_proposal,args=(size,),
                    name='torchgwas-productive-planner',daemon=True)
                self._planner_worker.start()
            else:self._proposal_step(size,release_issue_frontier=False)

    def finish(self,*,successful):
        with self._lock:self._finishing=True
        source_stage=None if self._source_stage is None else self._source_stage.finish()
        # The optional worker owns calculator/cache state until it returns.
        # Join outside the controller lock so it can perform its final audit.
        worker=self._planner_worker
        if worker is not None and worker is not threading.current_thread():worker.join()
        screen_worker=self._staged_screen_worker
        if screen_worker is not None and screen_worker is not threading.current_thread():
            screen_worker.join()
        with self._lock:
            if not self._finished:
                if self.refresh is not None:self.refresh.close()
                try:state=self.run.finish(successful=successful)
                except BaseException:
                    if self._admission_cache is not None:self._admission_cache.publish(successful=False)
                    self._finish_binding_digests(successful=False)
                    raise
                admission=None
                if self._admission_cache is not None:
                    self._admission_cache.publish(successful=successful)
                    admission=self._admission_cache.audit()
                bindings=self._finish_binding_digests(successful=successful)
                history=(None if self._planning_cost_history is None else
                    self._planning_cost_history.finish(successful=successful))
                self.audit['productive']=dict(state,forecast_attempts=deepcopy(self._attempts),
                    validation=deepcopy(self._validation),binding_digests=bindings,
                    admission_cache=admission,output_sample=self.output_sample.snapshot(),
                    output_boundary=self.boundary.snapshot())
                if history is not None:self.audit['productive']['planning_cost_history']=history
                if source_stage is not None:self.audit['productive']['source_staging']=source_stage
                if self._staged_screen_request is not None:
                    self.audit['productive']['staged_screen_evidence']=(
                        deepcopy(self._staged_screen_evidence) if self._staged_screen_evidence is not None
                        else dict(status='stage_not_complete',prediction_complete=False,
                                  selection_validated=False))
                if self.stage_sample is not None:
                    self.stage_sample.stop()
                    self.audit['productive']['stage_sample']=self.stage_sample.snapshot()
                    self.audit['productive']['stage_sample']['measurement_reserve_seconds']=self.config[
                        'stage_observations']['measurement_reserve_seconds']
                    self.audit['productive']['stage_sample']['known_probe_wall_seconds']=self.stage_sample.known_probe_wall_seconds()
                if self.refresh is not None:self.audit['productive']['resident_copy_refresh']=self.refresh.snapshot()
                self._header=None;self._prepared_header=None;self._prepared_index=None
                self._admission_cache=None;self._finished=True
            return self.audit['productive']

    def _finish_binding_digests(self,*,successful):
        cache=self._binding_digests
        if cache is None:return dict(enabled=False,initialization_error=self._binding_digest_error)
        began=time.perf_counter();cpu=time.thread_time()
        try:publication=cache.publish(successful=successful)
        except Exception as error:
            publication=dict(status='cache_error',error=type(error).__name__,stored=[])
        finally:cache.close()
        return dict(enabled=True,state=cache.snapshot(),publication=publication,
            publication_wall_seconds=time.perf_counter()-began,publication_cpu_seconds=time.thread_time()-cpu)
