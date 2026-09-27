"""Guarded execution bridge for the independently priced detailed calculator.

The search is explicit and finite. Unsupported inputs or stale calibration fail
before a writer exists; this module never substitutes the coarse calculator.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np

from .detailed_calibration import (execution_context, read_detailed_profile,
    validate_detailed_profile)
from .mechanistic_plan import _integer
from .trait_candidate_space import bounded_trait_plan, _axis
from .significant_candidate_space import bounded_significant_host_plan
from .jagwas_candidate_space import bounded_jagwas_plan


def _digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),
                                    allow_nan=False).encode()).hexdigest()


def validate_autotune_request(*, genotype, phenotype, output_dir, options):
    """Check mode conflicts before loading a potentially huge phenotype."""
    if not isinstance(genotype,(str,Path)) or Path(genotype).suffix.lower()!='.pgen':
        raise ValueError('Detailed autotune requires an explicit .pgen input path')
    if isinstance(phenotype,(str,Path)) and Path(phenotype).suffix.lower()!='.npy':
        raise ValueError('Detailed autotune requires mmap-capable .npy phenotype input')
    if output_dir is None:
        raise ValueError('Detailed autotune requires output_dir')
    required=dict(genotype_format=('auto','pgen'),pgen_mode=('hardcall',),
        compute_dtype=('auto','float32'),device=('auto',),sumstats_format=('binary',),
        sumstats_fsync=(True,),reduce=(None,'significant','jagwas'))
    for key,allowed in required.items():
        if options[key] not in allowed or (key=='sumstats_fsync' and options[key] is not True):
            raise ValueError('Detailed autotune does not support '+key+'='+repr(options[key]))
    for key in ['pipeline_profile','chunk_size','trait_block','trait_devices',
                'reader_workers','prefetch_chunks','pgen_decode_workers','variant_devices',
                '_internal_reduction','variant_range','p_value_threshold',
                'phenotype_table','covariates_table','hardcall_store','zstd_read_workers']:
        if options[key] is not None:
            raise ValueError('Detailed autotune controls or does not support '+key)
    if options['reduce'] in ('significant','jagwas') and options['sumstats_block_bytes'] is not None:
        raise ValueError('Detailed autotune reductions require indexed output without dense block coalescing')


class DetailedAutotune:
    """Validate one live context, then choose and audit one executable plan."""
    def __init__(self, profile, config, *, input_path, output_path, reduction=None, significance_threshold=None):
        started=time.perf_counter()
        if reduction not in (None,'significant','jagwas'):
            raise ValueError('Detailed autotune does not support this reduction mode')
        if reduction=='significant' and significance_threshold is not None and (
                isinstance(significance_threshold,bool) or not isinstance(significance_threshold,(int,float))
                or not math.isfinite(significance_threshold) or not 0<significance_threshold<=1):
            raise ValueError('Detailed autotune significance threshold must be in (0,1] or None')
        if reduction=='jagwas' and significance_threshold is not None:
            raise ValueError('JAGWAS emits every valid joint statistic; significance_threshold applies only to significant pairs')
        self.reduction=reduction
        self.significance_threshold=significance_threshold if reduction=='significant' else None
        self.profile=(read_detailed_profile(profile) if isinstance(profile,(str,Path))
                      else copy.deepcopy(profile))
        if isinstance(config,(str,Path)):
            config=json.loads(Path(config).read_text())
        self.config=copy.deepcopy(config)
        self.productive=None
        self.refresh=None
        price_fields={'significant':'significant_host_prices','jagwas':'jagwas_services'}
        if not isinstance(config,dict) or set(config)-({'plan_cache_dir','initial_chunks'}|set(price_fields.values()))!={'bounds','joint','qc_trait_block'}:
            raise ValueError('autotune_config requires exactly bounds, joint and qc_trait_block, with optional plan_cache_dir and mode-specific component artifact')
        price_field=price_fields.get(reduction)
        if any(field in config and field!=price_field for field in price_fields.values()):
            raise ValueError('Component artifact applies only to its matching reduction mode')
        self.price_artifact=config.get(price_field)
        if price_field is not None:
            if not isinstance(self.price_artifact,(str,Path)) or not str(self.price_artifact).strip():
                raise ValueError('Detailed '+reduction+' autotune requires a bound '+price_field+' artifact')
            self.price_artifact=str(Path(self.price_artifact).expanduser().resolve(strict=True))
            self.config[price_field]=self.price_artifact
        self.cache_directory=config.get('plan_cache_dir')
        if self.cache_directory is not None and (not isinstance(self.cache_directory,(str,Path)) or not str(self.cache_directory).strip()):
            raise ValueError('plan_cache_dir must be a nonempty path or None')
        if isinstance(self.cache_directory,Path):self.config['plan_cache_dir']=str(self.cache_directory)
        bounds=config['bounds'];joint=config['joint']
        required_axes={'chunks'} if reduction=='jagwas' else {'chunks','trait_blocks'}
        allowed_bounds=required_axes|{'max_candidates','max_candidate_tiles','max_census_chunks'}
        if reduction is None:allowed_bounds.add('partition_axes')
        if not isinstance(bounds,dict) or not required_axes<=set(bounds):
            raise ValueError('Explicit chunk and, outside JAGWAS, trait axes required')
        if set(bounds)-allowed_bounds:
            raise ValueError('Unknown detailed autotune bound; JAGWAS excludes phenotype partitioning')
        _axis('chunks',bounds['chunks'])
        _integer('qc_trait_block',config['qc_trait_block'])
        if reduction!='jagwas':
            _axis('trait_blocks',bounds['trait_blocks'])
            if config['qc_trait_block']>min(bounds['trait_blocks']):
                raise ValueError('QC trait block must not exceed the narrowest requested tile')
        required={'host_scenarios','cpu_workers','host_memory_bytes','device_memory_bytes',
                  'host_reserve_bytes','device_reserve_bytes','device_memory_profiles'}
        allowed=required|{'max_candidates','max_scenario_evaluations','max_tiles',
                         'max_chunk_evaluations','shortlist_size','max_slowdown_fraction'}
        if reduction in ('significant','jagwas'):
            required=required|{'occupancy_scenarios'}
            allowed=required|{'max_candidates','max_scenario_evaluations','max_source_chunks'}
            if reduction=='significant':allowed.add('max_selection_blocks')
        if not isinstance(joint,dict) or not required<=set(joint) or set(joint)-allowed:
            raise ValueError('Explicit search scenarios, capacities, reserves and device profiles required')
        if 'initial_chunks' in config:
            from .initial_chunk_autotune import validate_initial_chunk_config
            validate_initial_chunk_config(config['initial_chunks'],bounds,self.profile['contexts'],reduction,joint)
            if self.cache_directory is not None:
                raise ValueError('Deferred tuning does not reuse an upfront ranking; omit plan_cache_dir')
        for key in ['cpu_workers','host_memory_bytes','host_reserve_bytes','device_reserve_bytes']:
            _integer(key,joint[key],0 if key.endswith('reserve_bytes') else 1)
        self.devices=list(dict.fromkeys(d for c in self.profile['contexts'] for d in c['devices']))
        if set(joint['device_memory_bytes'])!=set(self.devices) or set(joint['device_memory_profiles'])!=set(self.devices):
            raise ValueError('Exactly one memory budget and workspace profile per calibrated device required')
        for value in joint['device_memory_bytes'].values():_integer('device_memory_bytes',value)
        self.input_path=str(Path(input_path).resolve(strict=True))
        self.output_path=str(Path(output_path).absolute())
        self.context=execution_context(self.devices,input_path=self.input_path,output_path=self.output_path)
        refresh_config=config.get('initial_chunks',{}).get('resident_copy_refresh')
        if refresh_config is not None:
            from .resident_copy_refresh import ResidentCopyRefresh
            self.refresh=ResidentCopyRefresh(self.profile,refresh_config,
                preserve_artifacts=[] if self.price_artifact is None else [self.price_artifact])
        options={} if self.refresh is None else dict(deferred_bindings=self.refresh.deferred_bindings)
        self.checked=validate_detailed_profile(self.profile,self.context,**options)
        if self.context['torch_default_dtype']!='torch.float32':
            raise ValueError('Detailed autotune requires the float32 default tensor dtype for its priced pinned buffers')
        for device in self.devices:
            memory=joint['device_memory_profiles'][device];gpu=self.context['devices'][device]
            for field in ['sm_count','max_threads_per_sm','compute_capability']:
                if memory.get(field)!=gpu[field]:raise ValueError('Device memory properties differ: '+device+'.'+field)
            if memory.get('torch_version')!=self.context['torch_version']:
                raise ValueError('Device memory PyTorch version differs')
            if memory.get('cublas_workspace_config')!=self.context['environment']['CUBLAS_WORKSPACE_CONFIG']:
                raise ValueError('Device memory cuBLAS workspace configuration differs')
            if memory.get('cublas_handle_stream_pairs',0)<2:
                raise ValueError('Detailed executor requires at least two cuBLAS handle/stream pairs')
        self.price_evidence=self._load_reduction_prices() if reduction is not None else None
        self.binding_seconds=time.perf_counter()-started

    def _load_reduction_prices(self):
        from .calibration_cache import read_calibration_record
        from .detailed_calibration import sha256_file
        digest=self.profile.get('component_artifacts',{}).get(self.price_artifact)
        if digest is None or sha256_file(self.price_artifact)!=digest:
            raise ValueError('Reduction component artifact is absent from or differs from the bound profile')
        name='significant_host_components' if self.reduction=='significant' else 'jagwas_components'
        evidence=read_calibration_record(self.price_artifact,kind='cpu_capacity',name=name,
            dependencies=dict(source_sha256=self.profile['source_sha256'],execution_context=self.context))
        prices=evidence['record']['value']
        if self.reduction=='significant':
            from .numpy_nonzero_work import validate_host_price_protocol
            if not isinstance(prices,dict) or set(prices)!={'host_selector','predicate_max_cells','nonzero_protocol','prices','archive','queue_cpu_seconds'}:
                raise ValueError('Bound significant prices must describe the current host selector and archive/queue services')
            validate_host_price_protocol(prices)
        elif (not isinstance(prices,dict) or set(prices)!={'writer_prices','preparation_services'}
              or not isinstance(prices['writer_prices'],dict)
              or set(prices['writer_prices'])!={'prices','archive','queue_cpu_seconds'}
              or not isinstance(prices['preparation_services'],dict) or not prices['preparation_services']):
            raise ValueError('Bound JAGWAS services must describe factor preparation and indexed writer services')
        evidence['artifact_sha256']=digest
        return evidence

    @property
    def qc_trait_block(self):
        return self.config['qc_trait_block']

    def validate_inputs(self, genotype, phenotype, covariates):
        from .pgen import PgenGenotype
        if not isinstance(genotype,PgenGenotype) or genotype.mode!='hardcall' or genotype.pgen_backend!='native':
            raise ValueError('Detailed autotune requires native hardcall PGEN')
        if (genotype._reader_sample_subset is not None or not genotype._reorder_is_identity
            or genotype.native_dtype!=np.dtype('int8') or genotype.native_encoding=='pgen_2bit'):
            raise ValueError('Detailed autotune requires native int8 and every sample in file order')
        if str(Path(genotype.genotype_path).resolve(strict=True))!=self.input_path:
            raise ValueError('Detailed autotune input path changed')
        if (not isinstance(phenotype,np.ndarray) or phenotype.ndim!=2
            or phenotype.dtype!=np.dtype('float32') or not phenotype.flags.c_contiguous):
            raise ValueError('Detailed autotune currently requires a C-contiguous float32 phenotype panel')
        if phenotype.shape[0]!=genotype.shape[0]:raise ValueError('Phenotype rows differ from PGEN samples')
        if covariates is not None and (np.ndim(covariates)!=2 or np.shape(covariates)[0]!=genotype.shape[0]):
            raise ValueError('Detailed autotune requires aligned two-dimensional covariates')
        if self.reduction!='jagwas':
            _axis('trait_blocks',self.config['bounds']['trait_blocks'],phenotype.shape[1])

    def select(self, genotype, phenotype, covariates, qc, *, output):
        from .preprocess import _covariate_basis
        import torch
        started=time.perf_counter()
        if qc['phenotype_missing_cells'] or qc['dropped_phenotype_columns']:
            raise ValueError('Detailed autotune currently requires complete retained phenotype/covariate columns')
        self.validate_inputs(genotype,phenotype,covariates)
        columns=0 if covariates is None else int(covariates.shape[1])
        basis=_covariate_basis(covariates) if columns else None
        rank=0 if basis is None else int(basis.shape[1])
        if columns>=genotype.shape[0]-2:raise ValueError('Detailed autotune requires fewer covariate columns than N-2')
        if self.reduction=='jagwas' and phenotype.shape[1]>genotype.shape[0]-rank-1:
            raise ValueError('JAGWAS trait count exceeds the residual phenotype rank')
        from .setup_work import setup_primitive_bank
        for context in self.profile['contexts']:
            for profile in context['profiles'].values():setup_primitive_bank(profile,rank)
        workload=dict(genotype=self.input_path,samples=int(genotype.shape[0]),
            markers=int(genotype.shape[1]),traits=int(phenotype.shape[1]),covariates=rank,
            matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
        if columns!=rank:workload['covariate_columns']=columns
        # Bounded input QC may take longer than a short empirical lifetime.
        # Recheck the original immutable record, never substitute a newer one.
        if self.profile.get('price_bindings') is not None:
            from .price_binding import validate_price_bindings
            options={} if self.refresh is None else dict(deferred_bindings=self.refresh.deferred_bindings)
            validate_price_bindings(self.profile,**options)
        if self.reduction is not None:self.price_evidence=self._load_reduction_prices()
        if 'initial_chunks' in self.config:
            from .initial_chunk_autotune import prepare_public_initial_chunks
            return prepare_public_initial_chunks(self,workload,output,basis,started,
                _prepared_header=getattr(genotype,'_prepared_native_header',None))
        cache=None;plan=None
        if self.cache_directory is not None:
            from .analytical_plan_cache import AnalyticalPlanCache,input_identity,input_is_stable
            cache_input=input_identity(self.input_path)
            cache=AnalyticalPlanCache(self.cache_directory,dict(profile_sha256=_digest(self.profile),
                config={key:value for key,value in self.config.items() if key!='plan_cache_dir'},
                workload=workload,output=output,input=cache_input,reduction=self.reduction,
                significance_threshold=self.significance_threshold))
            if input_is_stable(cache_input):plan=cache.load()
            else:cache.state='recent_input'
            # Older dense cache entries retained only a short display list.
            # Rebuild their bounded search once so live admission can consider
            # every previously feasible layout, without repricing at runtime.
            if self.reduction is None and plan is not None:
                admission=plan.get('admission_candidates')
                if (not isinstance(admission,list) or not admission
                        or len(admission)!=plan.get('candidates_feasible')
                        or not all(isinstance(row,dict) for row in admission)
                        or sum(row.get('candidate_index')==plan['selected'].get('candidate_index')
                               for row in admission)!=1):
                    plan=None;cache.state='schema_upgrade'
        if plan is None:
            if self.reduction=='significant':
                plan=bounded_significant_host_plan(workload,self.profile['contexts'],bounds=self.config['bounds'],
                    joint=self.config['joint'],output=output,prices=self.price_evidence['record']['value'],
                    significance_threshold=self.significance_threshold)
            elif self.reduction=='jagwas':
                services=self.price_evidence['record']['value']
                plan=bounded_jagwas_plan(workload,self.profile['contexts'],bounds=self.config['bounds'],
                    joint=self.config['joint'],output=output,prices=services['writer_prices'],
                    preparation_services=services['preparation_services'])
            else:
                plan=bounded_trait_plan(workload,self.profile['contexts'],bounds=self.config['bounds'],
                                        joint=self.config['joint'],output=output)
        # Keep cached analytical rankings immutable. Fresh resource admission
        # chooses the first still-feasible reduced-output candidate on every run.
        choices=plan['candidates'] if self.reduction is not None else plan['admission_candidates']
        devices=list(dict.fromkeys(d for candidate in choices for d in candidate['devices']))
        available={device:int(torch.cuda.mem_get_info(device)[0]+torch.cuda.memory_reserved(device)
                              -torch.cuda.memory_allocated(device)) for device in devices}
        from .api import _available_host_bytes
        host_available=_available_host_bytes();selected=None;live_rejected=[];live_feasible=[]
        for candidate in choices:
            if host_available<=0 or candidate['memory']['host_bytes']>host_available:
                reason='current host capacity'
            elif any(candidate['memory']['device_bytes'][d]>available[d] for d in candidate['devices']):
                reason='current device capacity'
            else:
                live_feasible.append(candidate)
                if self.reduction is not None:
                    selected=candidate;break
                continue
            live_rejected.append(dict(candidate_index=candidate.get('candidate_index'),reason=reason))
        if self.reduction is None and live_feasible:
            best=min(row['worst_supplied_scenario_seconds'] for row in live_feasible)
            maximum=best*(1+plan['max_slowdown_fraction'])
            selected=min((row for row in live_feasible if row['worst_supplied_scenario_seconds']<=maximum),
                key=lambda row:(len(row['devices']),row['worst_supplied_scenario_seconds'],row['candidate_index']))
        if selected is None:
            raise ValueError('No detailed autotune candidate fits '+', '.join(sorted({r['reason'] for r in live_rejected})))
        for key,value in selected['required_environment'].items():
            if os.getenv(key)!=value:raise ValueError('Detailed autotune execution environment differs: '+key)
        if torch.backends.cuda.matmul.allow_tf32 is not False:
            raise ValueError('Detailed autotune requires TF32 disabled')
        # Repeat context/source checks after planning, before starting any writer.
        current=execution_context(self.devices,input_path=self.input_path,output_path=self.output_path)
        final_checked=validate_detailed_profile(self.profile,current)
        if self.reduction is not None:self.price_evidence=self._load_reduction_prices()
        identity=plan['search_space']['input_file_identity'];stat=Path(self.input_path).stat()
        if any(identity[key]!=value for key,value in dict(bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,
                                                       device=stat.st_dev,inode=stat.st_ino).items()):
            raise ValueError('PGEN changed after detailed autotune census')
        if cache is not None:
            if input_identity(self.input_path)!=cache_input:
                raise ValueError('PGEN changed during analytical plan reuse')
            if cache.state!='hit':
                if input_is_stable(cache_input):cache.store(plan)
                else:cache.write_status='recent_input'
        # Keep the audit compact: per-operation execution graphs can be much
        # larger than the association metadata and are reproducible from inputs.
        def row(value):
            return {k:v for k,v in value.items() if k!='scenarios'}
        audit=dict(status='guarded_detailed_autotune',profile_sha256=_digest(self.profile),
            config=copy.deepcopy(self.config),config_sha256=_digest(self.config),
            output=copy.deepcopy(output),workload=workload,selected=row(selected),
            shortlist=[row(r) for r in (plan['candidates'] if self.reduction is not None else plan['shortlist'])],rejected=plan['rejected'],
            candidates_evaluated=plan['candidates_evaluated'],candidates_feasible=plan['candidates_feasible'],
            input_file_identity=identity,source_sha256=self.profile['source_sha256'],
            objective=plan['objective'],selection_validated=False,runtime_prediction_validated=False,
            unresolved_timing_terms=plan['unpriced_terms'] if self.reduction is not None else plan['unresolved_timing_terms'],
            context_matches=True,host_available_bytes=host_available,device_available_bytes=available,
            binding_seconds=self.binding_seconds,planning_and_validation_seconds=time.perf_counter()-started,
            scope=plan['scope'])
        audit.update(reduction=self.reduction,significance_threshold=self.significance_threshold,
            live_admission_rejections=live_rejected,planned_candidate_index=plan['selected'].get('candidate_index'))
        if 'price_evidence' in final_checked:audit['price_evidence']=copy.deepcopy(final_checked['price_evidence'])
        if self.price_evidence is not None:
            record=self.price_evidence['record']
            audit['reduction_calibration']={key:self.price_evidence[key]
                for key in ['record_sha256','artifact_sha256','age_seconds']}
            audit['reduction_calibration'].update({key:copy.deepcopy(record[key])
                for key in ['observed_unix_seconds','max_age_seconds','provenance']})
        if cache is not None:audit['analytical_cache']=cache.audit()
        return selected['api_kwargs'],audit,basis
