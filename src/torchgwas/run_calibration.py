"""Public-run binding for bounded initial production-chunk observations.

This collects/validates interval evidence; it neither ranks configurations nor
turns overlapping intervals into independent hardware capacities. Only the
first scan partition on each GPU is measured, within one shared early window.
"""
from copy import deepcopy
from pathlib import Path
import threading
import time
import uuid

from .calibration_cache import CalibrationParameterCache, _age
from .initial_calibration import InitialCalibrationController


OPTIONS = {'validation_chunks_per_device', 'max_chunks_per_device',
    'max_age_seconds', 'ratio_limit', 'absolute_floor_seconds', 'warmup_chunks',
    'stride', 'max_window_seconds', 'cuda_events', 'max_decision_cpu_seconds','reuse_binding_digests'}


def calibration_options(config):
    if not isinstance(config,dict) or set(config)-OPTIONS != {'cache_dir'}:
        raise ValueError('initial_calibration requires cache_dir and supported bounded-window options')
    directory=config['cache_dir']
    if not isinstance(directory,(str,Path)) or not str(directory).strip():
        raise ValueError('initial_calibration cache_dir must be a nonempty path')
    if type(config.get('reuse_binding_digests',True)) is not bool:
        raise ValueError('reuse_binding_digests must be boolean')
    options={key:deepcopy(value) for key,value in config.items() if key not in ('cache_dir','reuse_binding_digests')}
    _age(options.get('max_age_seconds',300.))
    # Reuse the controller's complete numerical validation without cache I/O,
    # CUDA initialization, source work, or a test association.
    class EmptyCache:
        def lookup(self,*args,**kwargs):return dict(hit=False,reason='validation')
    InitialCalibrationController(['cuda:0'],cache=EmptyCache(),dependencies={'validation':True},
        provenance={'validation':True},**options)
    return str(Path(directory).expanduser().absolute()),options


def validate_initial_calibration(config, *, genotype, options):
    calibration_options(config)
    if not isinstance(genotype,(str,Path)) or Path(genotype).suffix.lower()!='.pgen':
        raise ValueError('initial_calibration requires an explicit native hardcall .pgen path')
    if options['pgen_mode']!='hardcall' or options['compute_dtype'] not in ('auto','float32'):
        raise ValueError('initial_calibration requires hardcall PGEN and float32 computation')
    if options['output_dir'] is None or options['sumstats_format']!='binary':
        raise ValueError('initial_calibration requires productive binary output')
    if options['phenotype_table'] is not None or options['covariates_table'] is not None:
        raise ValueError('initial_calibration currently accepts aligned arrays or .npy inputs')
    from .adaptive_chunks import _positive_size
    _positive_size(options['chunk_size'],'explicit initial_calibration chunk_size')
    if options['pipeline_profile'] is not None or options['autotune_profile'] is not None:
        raise ValueError('initial_calibration starts from an explicit configuration without an upfront profile search')
    if options['_internal_reduction'] is not None:
        raise ValueError('initial_calibration supports public full, significant and jagwas output')


class _PartitionWindow(InitialCalibrationController):
    def __init__(self, owner, *args, **kwargs):
        self.owner=owner
        super().__init__(*args,**kwargs)

    def reserve_read(self, start, end, device):
        if not self.owner.reserve_window():return False
        return super().reserve_read(start,end,device)


class RunCalibration:
    def __init__(self, config, *, inputs, output_path, started_perf_counter=None):
        self._created=time.perf_counter() if started_perf_counter is None else started_perf_counter
        directory,self.options=calibration_options(config)
        self.cache=CalibrationParameterCache(directory)
        self._reuse_binding_digests=config.get('reuse_binding_digests',True)
        self._binding_digests=None
        self._inputs=dict(inputs);self.output_path=str(Path(output_path).absolute())
        self.run_id=str(uuid.uuid4());self._lock=threading.RLock()
        self._identities=self._input_identities()
        self._started=None;self._finished=None;self._windows={};self._scans={}
        self._errors={};self._context=None;self._source=None;self._binding_seconds=0.
        self._request=None;self._devices=[]
        self._lookup_seconds=0.
        self._output_bins=[];self._output_bin_counts={};self._output_observed=None
        self._output_binding=None;self._output_lookup=None;self._output_cache_seconds=0.
        self._output_sample_status='not_started';self._output_unbound=0
        self._output=dict(events=0,indexed_chunks=0,dense_ranges=0,indexed_rows=0,dense_statistic_cells=0,
            indexed_part_bytes=0,dense_statistic_bytes=0,first_written=None,first_material=None,
            first_fsynced_part=None,last_written=None)

    def __call__(self, observation):
        raise RuntimeError('Bind the calibration observer to its scan before delivery')

    def for_scan(self, device, variant_range, n_traits, *, reader_workers, capacity, depth):
        """Full-panel binding used by the native variant-sharding driver."""
        return self.observer(device,variant_range=variant_range,trait_range=(0,n_traits),
            reader_workers=reader_workers,capacity=capacity,depth=depth)

    def _input_identities(self):
        result={}
        for name,value in self._inputs.items():
            if value is None:result[name]=None
            elif isinstance(value,(str,Path)):
                path=Path(value).expanduser().resolve(strict=True);stat=path.stat()
                result[name]=dict(path=str(path),device=stat.st_dev,inode=stat.st_ino,
                    bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,ctime_ns=stat.st_ctime_ns)
            else:
                # Avoid an extra full read/hash of a potentially huge panel.
                # Mutable in-memory inputs deliberately cannot hit across jobs.
                result[name]=dict(in_memory=True,run_id=self.run_id)
        return result

    def prepare(self, genotype, *, devices, request):
        from .adaptive_chunks import validate_chunk_control
        from .detailed_calibration import execution_context,source_identity
        import torch
        started=time.perf_counter()
        devices=list(dict.fromkeys(map(str,devices)))
        if not devices or any(not d.startswith('cuda:') for d in devices):
            raise ValueError('initial_calibration requires explicit CUDA scan devices')
        for device in devices:
            validate_chunk_control(genotype,torch.device(device),'float32',request['capacity'],None,lambda row:None)
        if self._input_identities()!=self._identities:
            raise ValueError('Calibration input identity changed before the scan')
        if self._reuse_binding_digests:
            from .binding_digests import BindingDigestCache
            try:self._binding_digests=BindingDigestCache(self.cache.directory)
            except (OSError,ValueError) as error:
                self._errors['binding_digests']=dict(stage='cache_initialization',error=str(error))
        self._source=source_identity(digest_cache=self._binding_digests)
        self._context=execution_context(devices,input_path=self._inputs['genotype'],output_path=self.output_path,
            digest_cache=self._binding_digests)
        self._devices=devices;self._request=deepcopy(request)
        if request.get('reduction') in ('significant','jagwas'):
            self._output_binding=dict(protocol='torchgwas.initial_output_survivors.v1',
                source_sha256=self._source,execution_context=self._context,inputs=self._identities,request=self._request)
        self._binding_seconds=time.perf_counter()-started

    def reserve_window(self):
        with self._lock:
            if self._finished is not None:return False
            now=time.perf_counter()
            if self._started is None:self._started=now
            return now-self._started<self.options.get('max_window_seconds',10.)

    def observer(self, device, *, variant_range, trait_range, reader_workers, capacity, depth):
        """Bind once per GPU; later phenotype partitions add no event samples."""
        device=str(device)
        scan=dict(device=device,variant_range=list(variant_range),trait_range=list(trait_range),
                  reader_workers=int(reader_workers),capacity=int(capacity),depth=int(depth))
        with self._lock:
            if self._context is None or device not in self._devices:
                raise ValueError('Calibration scan differs from the prepared CUDA context')
            if self._finished is not None:raise ValueError('Calibration job is finished')
            if device in self._scans:
                return self._windows.get(device) if self._scans[device]==scan else None
            self._scans[device]=scan
            dependencies=dict(protocol='torchgwas.public_initial_chunks.v1',source_sha256=self._source,
                execution_context=self._context,inputs=self._identities,request=self._request,scan=scan)
            started=time.perf_counter()
            try:
                window=_PartitionWindow(self,[device],cache=self.cache,dependencies=dependencies,
                    provenance=dict(run_id=self.run_id,output_path=self.output_path,
                        scope='First production partition per GPU; loaded intervals, not hardware capacities'),**self.options)
            except OSError as error:
                self._errors[device]=dict(stage='cache_lookup',error=str(error));return None
            finally:
                self._lookup_seconds+=time.perf_counter()-started
            self._windows[device]=window
            return window

    def output_written(self, event):
        """Constant-space output progress audit, never a capacity measurement."""
        from .sumstats import DenseWriteProgress
        from .sumstats_indexed import IndexedChunkWrite
        if not isinstance(event,(DenseWriteProgress,IndexedChunkWrite)):
            raise ValueError('Typed writer completion event required')
        with self._lock:
            if self._finished is not None:raise ValueError('Calibration job is finished')
            row=self._output;row['events']+=1
            dense=isinstance(event,DenseWriteProgress)
            if dense:
                row['dense_ranges']+=1;row['dense_statistic_cells']+=event.rows
                row['dense_statistic_bytes']+=event.statistic_bytes
            else:
                row['indexed_chunks']+=1;row['indexed_rows']+=event.rows
                row['indexed_part_bytes']+=event.part_bytes
                self._observe_output_counts(event)
            for key,present in [('first_written',True),('first_material',event.rows>0),
                ('first_fsynced_part',not dense and event.part_file_fsynced)]:
                if present:row[key]=event.completed if row[key] is None else min(row[key],event.completed)
            row['last_written']=event.completed if row['last_written'] is None else max(row['last_written'],event.completed)

    @staticmethod
    def _output_bin_key(row):
        return (row['device'],*row['partition_variant_range'],*row['trait_range'],*row['variant_range'])

    def _validated_output_bins(self,value,*,previous=()):
        """Validate a bounded immutable sample before exposing it as evidence."""
        if not isinstance(value,dict) or set(value)!={'bins'} or not isinstance(value['bins'],list) or not 0<len(value['bins'])<=128:
            raise ValueError('Bounded output survivor sample required')
        rows=value['bins'];seen=list(previous)
        for row in rows:
            if not isinstance(row,dict) or set(row)!={'device','partition_variant_range','variant_range','trait_range','retained','part_bytes','part_file_fsynced'}:
                raise ValueError('Explicit output survivor coordinates required')
            if row['device'] not in self._devices:raise ValueError('Unknown output producer device')
            for name in ('partition_variant_range','variant_range','trait_range'):
                span=row[name]
                if not isinstance(span,list) or len(span)!=2 or any(type(x) is not int for x in span) or not 0<=span[0]<span[1]:
                    raise ValueError('Invalid output survivor range')
            lo,hi=row['variant_range'];a,b=row['partition_variant_range'];c,d=row['trait_range']
            outer=self._request['variant_range'];k=self._request['phenotype_shape'][1]
            if not outer[0]<=a<=lo<hi<=b<=outer[1] or d>k:raise ValueError('Output survivor range differs from the request')
            if self._request['reduction']=='jagwas' and [c,d]!=[0,k]:raise ValueError('JAGWAS requires a complete phenotype panel')
            cells=(hi-lo)*(d-c if self._request['reduction']=='significant' else 1)
            if (type(row['retained']) is not int or not 0<=row['retained']<=cells or
                    type(row['part_bytes']) is not int or row['part_bytes']<0 or type(row['part_file_fsynced']) is not bool or
                    (row['retained']==0 and (row['part_bytes']!=0 or row['part_file_fsynced']))):
                raise ValueError('Invalid output survivor count or durability')
            for other in seen:
                x,y=other['variant_range'];u,v=other['trait_range']
                if max(lo,x)<min(hi,y) and max(c,u)<min(d,v):raise ValueError('Overlapping output survivor evidence')
            seen.append(row)
        return sorted(deepcopy(rows),key=self._output_bin_key)

    def _observe_output_counts(self,event):
        if self._output_binding is None:return
        if self._output_sample_status=='invalid':return
        partition=event.partition
        if partition is None or event.source_variant_range is None:
            self._output_unbound+=1;return
        maximum=self.options.get('max_chunks_per_device',8)
        if len(self._output_bins)>=128 or self._output_bin_counts.get(partition.device,0)>=maximum:return
        start=self._started if self._started is not None else self._created
        if time.perf_counter()-start>=self.options.get('max_window_seconds',10.):return
        row=dict(device=partition.device,partition_variant_range=list(partition.variant_range),
            variant_range=list(event.source_variant_range),trait_range=list(partition.trait_range),
            retained=event.rows,part_bytes=event.part_bytes,part_file_fsynced=event.part_file_fsynced)
        try:
            if event.kind!=self._request['reduction']:raise ValueError('Output mode differs from request')
            self._validated_output_bins(dict(bins=[row]),previous=self._output_bins)
        except (ValueError,KeyError,TypeError):
            self._output_sample_status='invalid';self._output_bins.clear();return
        if self._output_lookup is None:
            began=time.perf_counter()
            try:
                self._output_lookup=self.cache.lookup('stage_observations','initial_output_survivors.v1',
                    dependencies=self._output_binding,max_age_seconds=self.options.get('max_age_seconds',300.))
                if self._output_lookup['hit']:
                    self._validated_output_bins(self._output_lookup['record']['value'])
            except (OSError,ValueError,KeyError,TypeError) as error:
                self._output_lookup=dict(hit=False,reason='invalid_or_unavailable_output_cache',error=str(error))
            finally:self._output_cache_seconds+=time.perf_counter()-began
        self._output_bins.append(row);self._output_bin_counts[partition.device]=self._output_bin_counts.get(partition.device,0)+1
        # Age the aggregate from its oldest included completion, not job finish
        # or a later cache-validation/publication time.
        observed=time.time()-max(0.,time.perf_counter()-event.completed)
        self._output_observed=observed if self._output_observed is None else min(self._output_observed,observed)
        self._output_sample_status='sampled'

    def _finish_output_counts(self,successful):
        rows=sorted(deepcopy(self._output_bins),key=self._output_bin_key)
        previous=self._output_lookup or dict(hit=False,reason='not_looked_up')
        status='not_published';publication=None;error=None
        if successful and rows and self._output_sample_status!='invalid':
            value=dict(bins=rows)
            old=None if not previous['hit'] else self._validated_output_bins(previous['record']['value'])
            if rows==old:status='reused_without_renewal'
            else:
                try:
                    publication=self.cache.store('stage_observations','initial_output_survivors.v1',value,
                        dependencies=self._output_binding,provenance=dict(run_id=self.run_id,output_path=self.output_path,
                            scope='Completed output survivor counts on explicit producer partitions; observed ranges only'),
                        max_age_seconds=self.options.get('max_age_seconds',300.),observed_unix_seconds=self._output_observed)
                    status='published'
                except (OSError,ValueError) as caught:status='cache_error';error=str(caught)
        return dict(status=status,sample_status=self._output_sample_status,bins=rows,
            reduction=None if self._request is None else self._request.get('reduction'),
            total_traits=None if self._request is None else self._request.get('phenotype_shape',[None,None])[1],
            significance_threshold=None if self._request is None else self._request.get('significance_threshold'),
            max_bins=128,max_bins_per_device=self.options.get('max_chunks_per_device',8),
            unbound_events=self._output_unbound,observed_unix_seconds=self._output_observed,
            cache_hit=previous['hit'],cache_reason=previous.get('reason'),
            previous_record_sha256=previous.get('record_sha256'),cache_lookup_seconds=self._output_cache_seconds,
            publication=publication,cache_error=error,
            scope='Bounded completed output counts in source-file and post-QC phenotype coordinates. Dependency/age-bound observations, not immutable constants, independent writer capacity or future selectivity. Cache reuse never renews original age.')

    def _output_report(self):
        row=deepcopy(self._output)
        for key in ('first_written','first_material','first_fsynced_part','last_written'):
            value=row.pop(key);row[key+'_seconds']=None if value is None else value-self._created
        row['elapsed_to_report_seconds']=time.perf_counter()-self._created
        row['scope']='Times from API entry (or controller construction for direct callers). Dense events follow beta/t writes and may precede df flushing and fsync. Indexed part fsync is distinct from final manifest durability. Final report rewriting and function return follow this snapshot. Output progress is not independent writer capacity.'
        return row

    def finish(self, *, successful):
        from .detailed_calibration import source_identity
        if type(successful) is not bool:raise ValueError('successful must be boolean')
        with self._lock:
            if self._finished is not None:return deepcopy(self._finished)
            started=time.perf_counter()
            unchanged=False
            if successful:
                try:
                    unchanged=(self._input_identities()==self._identities and
                        source_identity(digest_cache=self._binding_digests)==self._source and
                        (self._binding_digests is None or self._binding_digests.unchanged()))
                except (OSError,ValueError):unchanged=False
            windows={}
            for device,window in self._windows.items():
                window.finish(successful=successful and unchanged)
                windows[device]=window.snapshot()
            output_counts=self._finish_output_counts(successful and unchanged)
            binding_digests=None
            if self._binding_digests is not None:
                began=time.perf_counter();cpu=time.thread_time()
                try:
                    publication=self._binding_digests.publish(successful=successful and unchanged)
                except Exception as error:
                    publication=dict(status='cache_error',error=type(error).__name__,stored=[])
                finally:self._binding_digests.close()
                binding_digests=dict(state=self._binding_digests.snapshot(),publication=publication,
                    publication_wall_seconds=time.perf_counter()-began,publication_cpu_seconds=time.thread_time()-cpu)
            self._finished=dict(successful=successful,publication_inputs_unchanged=unchanged,
                run_id=self.run_id,binding_seconds=self._binding_seconds,
                cache_lookup_seconds=self._lookup_seconds,finalization_seconds=time.perf_counter()-started,
                shared_window_started=self._started,first_partitions=deepcopy(self._scans),windows=windows,
                cache_errors=deepcopy(self._errors),unobserved_devices=sorted(set(self._devices)-set(windows)),
                binding_digests=binding_digests,
                output_progress=self._output_report(),
                output_occupancy=output_counts,
                reuse_policy='File metadata and execution dependencies must match; in-memory inputs are job-specific. Reuse never renews original observation age.',
                scope='Bounded useful production-chunk calibration only. No upfront benchmark or candidate search, automatic configuration change, or hardware-capacity inference. Decision CPU budget is per active GPU; reservation wall window is shared.')
            return deepcopy(self._finished)
