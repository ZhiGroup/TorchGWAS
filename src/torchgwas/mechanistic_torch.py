"""Full-major-stage Torch runtime candidate from source work and resources.

This is not a calibrated predictor or an accuracy guarantee. Every approximation
is returned; no observed association time is accepted by the interface.
"""
from __future__ import annotations
import math
from .first_principles import positive,service
from .environment_service import environment_service
from .input_read_work import buffered_read_service
from .output_write_work import writer_copy_cost
from .decoder_work import decoder_work,decoder_chunk_work,price_identified_units,census_chunk_ranges,_native_decoder_memory_bytes
from .tensor_work import eager_statistics_work
from .tensor_service import DeviceService,tensor_stage_service
from .binary_output_work import binary_output_work
from .binary_schedule import BinaryWriterSchedule
from .execution_graph import torch_scan_schedule
from .native_control_work import native_control_service
from .setup_work import setup_work,setup_service
from .pinned_work import pinned_scan_work
from .owned_result_work import owned_result_work,owned_result_copy_service,owned_result_allocation_policy,owned_result_copy_threading

def torch_runtime(data,profile):
 if profile.get('reduction') is not None:raise ValueError('Reduced output requires its indexed writer and preparation schedule')
 n,m,c,k=(data[t] for t in ['samples','markers','covariates','traits_analyzed'])
 if c!=8 or not data['matching_sample_order'] or k<1:raise ValueError('Candidate requires C8, complete phenotypes and matching sample order')
 if k!=1 and 'setup_primitives' not in profile:raise ValueError('Wide phenotypes require independent blocked setup primitives')
 if n<32 or m<1:raise ValueError('Unsupported dimensions')
 q=profile['cpu_fraction'];cores=profile['cpu_available_cores'];gpu_fraction=profile['gpu_resources']['gpu_fraction']
 if not 0<=q<=1:raise ValueError('CPU share must be in [0,1]')
 required=['cpu_fraction','cpu_available_cores','read_bytes_per_second','write_bytes_per_second','shared_dram_bytes_per_second','h2d_bytes_per_second','d2h_bytes_per_second']
 for name in required:
  if positive(name,profile[name],True)==0:return dict(estimated_seconds=None,status='zero_available_capacity',resource=name)
 if gpu_fraction==0:return dict(estimated_seconds=None,status='zero_available_capacity',resource='gpu_fraction')
 units=profile['process_units'];b=profile['chunk_markers'];depth=profile['depth'];workers=min(depth,profile['decode_workers'],math.ceil(m/b))
 # CPU scheduling share caps each active serial task; total process CPU
 # capacity is enforced by the finite graph's active-demand resource sharing.
 effective_q=q
 output=binary_output_work(m,k,b,store_variant_df=True,store_beta=profile.get('return_beta',True))
 data_bytes=n*(data['traits_in_file']+c)*profile['npy_itemsize']+2*128
 metadata_bytes=data['tables']['pvar']['bytes']+data['tables']['psam']['bytes']+data_bytes
 # Source parses index for scope and initial reader, then creates one reader
 # per decoder worker. Worker index initialization is charged in its first job.
 index_bytes=data['encoded'].get('index_bytes',data['encoded']['file_bytes']-data['encoded']['record_payload_bytes'])
 file_markers=data['encoded'].get('file_markers',m)
 setup_cpu=units['pvar_rows']*file_markers+units['psam_rows']*n
 setup_cpu+=units['phenotype_qc_cells']*n*k+units['covariate_qc_cells']*n*c
 setup_cpu+=units['covariate_basis_work']*(n*c*c+c**3)
 setup_cpu+=units['pgen_index_records']*file_markers*profile['initial_index_parses']
 # Cast loaded FP64 phenotype/covariates to FP32 and promote C for the rank SVD.
 setup_cpu+=units['numpy_copy_bytes']*(12*n*(k+c)+12*n*c)
 ready_boundary=profile.get('timing_boundary','process')=='environment-ready'
 context=0. if ready_boundary else profile['cuda_context_seconds']
 startup=dict(seconds=0.) if ready_boundary else environment_service(profile.get('environment_events'))
 import_service=startup['seconds']
 setup_io=metadata_bytes/profile['read_bytes_per_second']+index_bytes/profile['read_bytes_per_second']
 pin_work=pinned_scan_work(n,b,k,depth,return_beta=profile.get('return_beta',True))
 pinned_bytes=pin_work['requested_bytes']
 pinned=pin_work['allocation_pages']*profile['pin_cpu_seconds_per_page']/q+profile['pin_driver_seconds_per_page']*pin_work['allocation_pages']
 # Tiny source setup is fixed API/library overhead; only additional bulk work
 # beyond N32 is added. Source setup runs serially before the scan pipeline.
 dn=n-32
 setup_gpu=profile['tiny_setup_cpu_seconds']/q+profile['tiny_setup_non_cpu_seconds']/gpu_fraction
 setup_gpu+=max(4*dn*c*k/profile['gpu_resources']['fp32_flops_per_second'],
                4*dn*(6*c+18*k+4)/profile['gpu_resources']['hbm_bytes_per_second'])/gpu_fraction
 setup_gpu+=4*dn*(2*c+3*k)/profile['h2d_bytes_per_second']+4*dn*k/profile['d2h_bytes_per_second']
 setup_ledger=None;setup_estimate=None
 if 'setup_primitives' in profile:
  setup_ledger=setup_work(n,k,c);setup_estimate=setup_service(setup_ledger,profile);setup_gpu=setup_estimate['seconds']
 stages=dict(environment=import_service,cuda_context=context,metadata_io=setup_io,
             metadata_qc_basis=setup_cpu/q,pinned_allocation=pinned,gpu_design=setup_gpu,
             first_use=profile['first_use_seconds'])
 if 'design_first_use_seconds' in profile:
  stages['design_first_use']=positive('design_first_use_seconds',profile['design_first_use_seconds'],True)
 scan_work=torch_scan_work(data,profile)
 if scan_work.get('status')=='zero_available_capacity':return scan_work
 blocks=scan_work['blocks'];components=scan_work['components'];decode=scan_work['decoder']
 # Active write calls share one explicit output resource. A solitary writer
 # receives full capacity, while concurrent writes contend within the graph.
 writeback_service=profile.get('writeback_service')
 if writeback_service is not None:
  writeback_service={key:value if key=='storage_seconds_per_byte' else value/q for key,value in writeback_service.items()}
 copy_service=writer_copy_cost(profile)
 writer=BinaryWriterSchedule(output,copy_seconds_per_byte=copy_service['cpu_seconds_per_byte']/effective_q,
   copy_seconds_per_call=copy_service['cpu_seconds_per_call']/effective_q,
   zero_seconds_per_byte=units['bytearray_zero_bytes']/effective_q,
   write_seconds_per_byte=1/profile['write_bytes_per_second'],fsync_seconds_per_array=profile['fsync_seconds'],
   append_seconds=profile['executor_cpu_seconds']/effective_q,
   handoff_seconds=profile['executor_cpu_seconds']/effective_q,cpu_fraction=q,write_capacity=profile['write_bytes_per_second'],
   writeback_service=writeback_service)
 schedule=torch_scan_schedule(blocks,depth=depth,decode_workers=workers,consumer=writer,shared_capacities={'cpu':cores,'dram':profile['shared_dram_bytes_per_second'],'input':profile['read_bytes_per_second'],'output':profile['write_bytes_per_second']})
 stages['scan_binary_close']=schedule['seconds']
 sidecar_bytes=data['tables']['pvar']['field_characters']+4096
 stages['sidecars_exit']=m*units['variant_id_rows']/q+sidecar_bytes/profile['write_bytes_per_second']
 return dict(status='development_candidate',estimated_seconds=sum(stages.values()),stage_seconds=stages,
  prediction_complete=False,unpriced_terms=scan_work['unpriced_terms']+scan_work['allocator_unpriced_terms']+['pinned allocator cache/driver backing and fixed-call service','variant-df output validation/metadata and atomic manifest/directory fsync service']+(['kernel early writeback, dirty-page throttling and page-cache DRAM contention'] if profile.get('writeback_service') is not None else [])+(setup_estimate['unpriced_terms'] if setup_estimate else [])+([] if 'design_first_use_seconds' in profile else ['cold design-library initialization']),
  source_work=dict(setup=setup_ledger,setup_service=setup_estimate,input_array_bytes=data_bytes,pinned_bytes=pinned_bytes,pinned_allocations=pin_work,decoder=decode,host_workspace=scan_work['host_workspace'],encoded_work_distribution=scan_work['encoded_work_distribution'],binary_output={a:v for a,v in output.items() if a not in ['chunks','events_per_array','stream_work']}),
  schedule=dict(blocks=len(blocks),depth=depth,decode_workers=workers,seconds=schedule['seconds'],writer_payload_bytes=writer.payload,writer_copy_bytes=writer.copied,writer_write_calls=writer.write_calls),
  component_summaries={str(size):{key:component[key] for key in ['estimated_span_seconds','host_dispatch_cpu_seconds','kernel_service_seconds','gemm']} for size,component in components.items()},
  capacity=dict(cpu_fraction=q,effective_role_cpu_fraction=effective_q,cpu_available_cores=cores),
  assumptions=[
   'Source-major-stage candidate, not exact instruction timing or an approved crossover calculator.',
   'Independent fixed parser and dense-SVD source services scale by rows/cells or N*C^2+C^3; call overhead and instruction shape effects remain approximate.',
   'Encoded input work uses '+scan_work['encoded_work_distribution']+'; native LD restarts require exact base records and contiguous read prefixes.',
   'CPU and I/O service use proportional active-demand fluid sharing, an explicit scheduling scenario. Idle roles consume no capacity; OS priority and burstiness are not resolved.',
   'Active writer streams share aggregate write capacity. A private 1MiB commit control supplies fsync latency; actual dirty-page writeback is load dependent.',
   'Fresh pinned allocation scales with page count. Repeated reads of parsed indexes use resident page-cache memory after the first cold index read.',
   'Import and thread-pool startup uses explicit independent elapsed event services including waiting; scan resource overrides do not automatically change its environment state.',
   'Fixed tiny setup and generic lazy-library service are environment primitives; warm large-shape allocator/driver behavior remains approximate.',
   'GPU launch geometry is untimed compiled work, never an N-indexed service rate. Arbitrary missing geometry is refused, not interpolated.',
   'CPU decoder traffic and ideal GPU cache traces are analytical traffic approximations; periodic binary writeback requires explicit page-cache, storage and syscall service. Its schedule assumes ranges begin storage at submission; kernel early flushing and dirty throttling remain unresolved.',
   'Metadata lacks long-ID and missing-phenotype branches. All inputs remain explicit and no GWAS timing coefficient is used.'])


def _shape_component(n,size,k,c,profile,gpu,statistics_geometry,joint_geometry):
 """One exact tensor shape; retain source and rate bindings across JIT steps."""
 import torch
 from . import tensor_work,tensor_service,reduction_tensor_work,host_serial_work
 from .linear import _dosage_statistics
 from .jagwas_projection import JagwasReduction
 from .planning_session import cached_planning_work
 reduction=profile.get('reduction');q=profile['cpu_fraction']
 def build():
  work=eager_statistics_work(n,size,k,c,validate_range=profile.get('validate_range',False))
  component=tensor_stage_service(work,gpu,statistics_geometry['kernels'],host_primitives=profile['host_primitives'])
  if component.get('status')=='zero_available_capacity':return component
  if reduction=='jagwas':
   projection=reduction_tensor_work.jagwas_projection_service(n,size,k,gpu,joint_geometry['kernels'],
       host_primitives=profile['joint_host_primitives'])
   if projection.get('status')=='zero_available_capacity':return projection
   component=reduction_tensor_work.append_jagwas_projection(component,projection,q)
  if profile.get('host_serial_primitives') is not None:
   component=host_serial_work.attach_host_serial_work(component,profile['host_serial_primitives'],q)
  return component
 bindings=dict(shape=[n,size,k,c],validate_range=profile.get('validate_range',False),
     torch_default_dtype=str(torch.get_default_dtype()),gpu_resources=profile['gpu_resources'],
     cpu_fraction=q,reduction=reduction,statistics_geometry=statistics_geometry,
     joint_geometry=joint_geometry,host_primitives=profile['host_primitives'],
     joint_host_primitives=profile.get('joint_host_primitives'),
     host_serial_primitives=profile.get('host_serial_primitives'))
 return cached_planning_work('native_tensor_component.v1',bindings,build,
     implementation=(_shape_component,eager_statistics_work,tensor_work._trace_statistics_work,
         _dosage_statistics,tensor_stage_service,tensor_service.tensor_stage_service,
         reduction_tensor_work.jagwas_projection_service,reduction_tensor_work.jagwas_tensor_work,
         reduction_tensor_work._trace_jagwas_tensor_work,
         reduction_tensor_work.append_jagwas_projection,JagwasReduction.prepare,JagwasReduction.reduce,
         JagwasReduction._set_factor,
         host_serial_work.attach_host_serial_work))


def torch_scan_work(data,profile):
 """Exact-census scan work; header intervals require the explicit scenario API."""
 if data.get('encoded',{}).get('kind') in ('torchgwas.pgen_header_work_bounds.v1','torchgwas.pgen_header_window.v1'):
  raise ValueError('Header work bounds cannot substitute for an exact scan census')
 return _torch_scan_work(data,profile)


def torch_scan_header_work(data,profile,*,endpoint,issued_chunks=0,max_chunks=8,max_records=65536):
 """One bounded source-window scenario using the same scan component equations.

 Lower/upper select decoder CPU service endpoints, not makespan bounds. The
 resulting blocks preserve tensor, transfer, finish and reduction semantics;
 output consumption must still be composed by the mode-specific writer model.
 A future-work forecast additionally needs an explicit extrapolation policy.
 """
 from .pgen_work_bounds import WINDOW_KIND,native_decoder_service_bounds
 from .planning_session import planning_work_scope
 if endpoint not in ('lower','upper'):raise ValueError('Explicit lower or upper decoder endpoint required')
 for name,value,minimum in [('issued_chunks',issued_chunks,0),('max_chunks',max_chunks,1),('max_records',max_records,1)]:
  if type(value) is not int or value<minimum:raise ValueError('Invalid '+name)
 encoded=data['encoded'];rows=encoded.get('chunks')
 if encoded.get('kind')!=WINDOW_KIND or not isinstance(rows,list) or not rows or len(rows)>max_chunks:
  raise ValueError('Bounded typed header window required')
 if (type(encoded.get('markers')) is not int or not 0<encoded['markers']<=max_records
     or encoded.get('chunk_markers')!=profile['chunk_markers']):
  raise ValueError('Header window record budget or chunk size differs')
 if (type(encoded.get('file_markers')) is not int or encoded['file_markers']<encoded['variant_range'][1]
     or encoded.get('path')!=encoded.get('input_identity',{}).get('path')):
  raise ValueError('Header window file binding differs')
 ranges=census_chunk_ranges(encoded,profile['chunk_markers'])
 if len(ranges)!=len(rows) or encoded['markers']!=sum(hi-lo for lo,hi in ranges):
  raise ValueError('Header window coverage differs')
 intervals=[]
 for span,row in zip(ranges,rows):
  if (row.get('variant_range')!=span or row.get('markers')!=span[1]-span[0]
      or row.get('samples')!=encoded['samples'] or row.get('input_identity')!=encoded['input_identity']):
   raise ValueError('Header window source binding differs')
  priced=native_decoder_service_bounds(row,profile)
  intervals.append(priced['decode_resource_work']['cpu'])
 with planning_work_scope():
  result=_torch_scan_work(data,profile,header_scenario=dict(endpoint=endpoint,intervals=intervals,issued_chunks=issued_chunks))
 result['source_scenario']=dict(kind=WINDOW_KIND,endpoint=endpoint,variant_range=list(encoded['variant_range']),
     input_identity=encoded['input_identity'],issued_chunks=issued_chunks,decoder_cpu_intervals=intervals,
     scope='Conditional component service scenario on this source window only. Endpoint schedules are not elapsed-time bounds; output and future-source extrapolation are not included.')
 return result


def _torch_scan_work(data,profile,*,header_scenario=None):
 """General N/B/K/C scan service from source work and exact compiled geometry.

 Input is the native int8 eager path with complete phenotypes and matching
 sample order. Source loading/setup and output consumption are composed by
 callers; no K1 setup extrapolation is performed here.
 """
 if data.get('encoded',{}).get('kind')=='torchgwas.pgen_memory_layout.v1':
  raise ValueError('Memory-only PGEN layout cannot price scan work; exact decoder work evidence is required')
 reduction=profile.get('reduction')
 if reduction not in (None,'jagwas','device_significant'):raise ValueError('Unsupported scan reduction')
 return_beta=profile.get('return_beta',True)
 if type(return_beta) is not bool or (not return_beta and reduction is not None):
  raise ValueError('Beta omission requires an unreduced result layout')
 device_selection=reduction=='device_significant'
 ownership=profile.get('result_ownership','borrowed' if profile.get('borrow_results',False) else 'owned')
 if ownership not in ('owned','borrowed') or (profile.get('borrow_results',False) and ownership!='borrowed'):
  raise ValueError('Invalid result ownership')
 borrowed=ownership=='borrowed';finish_service=profile.get('result_finish_service')
 baseline_bytes=544 if reduction=='jagwas' else 32*(9+4*int(return_beta))
 if not return_beta and (finish_service is None or profile.get('control_primitives') is None
     or finish_service.get('return_beta') is not False or finish_service.get('result_arrays')!=3
     or finish_service.get('baseline_rows')!=32 or finish_service.get('return_df') is not False):
  raise ValueError('T-only results require an independent matching three-array finish service')
 if finish_service is not None and finish_service.get('return_beta',True)!=return_beta:
  raise ValueError('Finish service beta layout mismatch')
 if device_selection:
  if borrowed or profile.get('compute_dtype','float32')!='float32':
   raise ValueError('Device significance requires owned float32 selected results')
  if finish_service is not None:
   raise ValueError('Device significance has no dense result finish service')
  if profile.get('control_primitives') is None or profile.get('device_status_service') is None:
   raise ValueError('Independent native control and device status services required')
  if any(profile.get(key) is not None for key in ['owned_result_copy_scenario','owned_result_allocation_policy','owned_result_allocator']):
   raise ValueError('Device selection owns payloads in its source graph, not the dense result copier')
 elif reduction=='jagwas':
  if borrowed:raise ValueError('JAGWAS indexed scheduling requires owned results')
  if profile.get('compute_dtype','float32')!='float32':raise ValueError('Joint native calculator requires float32 statistics geometry')
  if (finish_service is None or profile.get('control_primitives') is None
      or finish_service.get('reduction')!='jagwas' or finish_service.get('result_arrays')!=5
      or finish_service.get('baseline_rows')!=32 or finish_service.get('return_df') is not False):
   raise ValueError('JAGWAS requires an independent matching five-array finish service')
 elif finish_service is not None and finish_service.get('reduction') is not None:
  raise ValueError('Finish service reduction mismatch')
 if borrowed and (finish_service is None or profile.get('control_primitives') is None):
  raise ValueError('Borrowed results require independent finish and ready-event service plus acknowledgement scheduling')
 if finish_service is not None:
  total=positive('finish CPU',finish_service['cpu_seconds'],True)
  held=positive('finish serial CPU',finish_service['serial_cpu_seconds'],True)
  if held>total or finish_service.get('replaces_fixed_finish_and_tensor_conversion') is not True or finish_service.get('includes_ready_cuda_event') is not False or finish_service.get('baseline_copy_bytes')!=(0 if borrowed else baseline_bytes):
   raise ValueError('Finish service coverage or ownership mismatch')
 n,m,c,k=(data[t] for t in ['samples','markers','covariates','traits_analyzed'])
 if not data['matching_sample_order']:raise ValueError('Matching sample order required')
 if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in (n,m,k)) or not isinstance(c,int) or c<0 or n<=c+2:
  raise ValueError('Invalid scan dimensions or residual degrees of freedom')
 if reduction=='jagwas' and k>n-c-1:raise ValueError('Joint phenotype rank exceeds residual degrees of freedom')
 b=profile['chunk_markers'];depth=profile['depth']
 if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in (b,depth,profile['decode_workers'])) or depth<2:
  raise ValueError('Invalid chunk, depth or worker count')
 effective_q=profile['cpu_fraction']
 if not 0<effective_q<=1:raise ValueError('Positive CPU fraction no greater than one required')
 for name in ['cpu_available_cores','read_bytes_per_second','shared_dram_bytes_per_second','h2d_bytes_per_second','d2h_bytes_per_second']:
  if positive(name,profile[name],True)==0:raise ValueError('Zero available resource '+name)
 storage_terms=['independence of supplied input capacity from per-reader buffered CPU service']
 if profile.get('input_storage_service') is not None:
  storage=profile['input_storage_service']
  if storage.get('service_kind')!='independent_direct_read_aggregate' or storage['bytes_per_second']!=profile['read_bytes_per_second']:
   raise ValueError('Input storage service provenance differs from supplied capacity')
  storage_terms=list(storage['unpriced_terms'])
 units=profile['process_units']
 census=data['encoded']
 ranges=census_chunk_ranges(census,b)
 workers=min(depth,profile['decode_workers'],len(ranges)+(header_scenario['issued_chunks'] if header_scenario else 0))
 if device_selection and 'chunk_ranges' in census:
  raise ValueError('Explicit schedules are not yet supported by the device-significance writer model')
 if (census['samples'],census['markers'])!=(n,m):raise ValueError('Encoded census dimensions differ from scan')
 if header_scenario is None and any(int(form) in (2,3) and count for form,count in census['record_form_counts'].items()) and 'chunks' not in census:
  raise ValueError('Exact per-chunk LD replay census required for parallel chunks')
 if header_scenario is None:
  decode=decoder_work(census,'torch_native_int8',restart_ld_bases=True)
  priced=price_identified_units(decode,profile['decode_units'])
  if priced['unpriced_source_units']:raise ValueError('Unpriced decoder units: '+str(priced['unpriced_source_units']))
  decoded_cpu=priced['identified_cpu_seconds']
  chunk_decoders=decoder_chunk_work(census,b)
  distribution='exact_encoded_chunks' if chunk_decoders is not None else 'uniform_aggregate_scenario'
 else:
  decode=dict(kind=census['kind'],native_ld_replay_records=sum(row['ld_replay'] is not None for row in census['chunks']))
  decoded_cpu=None;chunk_decoders=None;distribution='header_interval_'+header_scenario['endpoint']
 gpu=DeviceService(**dict(profile['gpu_resources'],host_cpu_fraction=effective_q))
 geometry={}
 for r in profile['kernel_geometry']:
  key=(r['N'],r['B'],r.get('K',1),r.get('C',8),r.get('validate_range',False))
  if key in geometry:raise ValueError('Duplicate compiled geometry '+str(key))
  geometry[key]=r
 components={};blocks=[];cpu_total=0.;result_work={}
 joint_geometry={}
 if reduction=='jagwas':
  for row in profile.get('joint_kernel_geometry',[]):
   key=(row['N'],row['B'],row['K'],row['compute_dtype'])
   if key in joint_geometry:raise ValueError('Duplicate joint compiled geometry')
   joint_geometry[key]=row
  if profile.get('joint_host_primitives') is None:raise ValueError('Independent joint host primitives required')
 # Use per-chunk source counts when supplied. Legacy aggregate-only inputs
 # retain their explicitly labeled uniform-work scenario.
 for i,(start,end) in enumerate(ranges):
  size=end-start
  if size not in components:
   if (n,size,k,c,profile.get('validate_range',False)) not in geometry:raise ValueError('Supply duration-free compiled geometry for N,B,K,C='+str((n,size,k,c)))
   key=(n,size,k,'float32')
   if reduction=='jagwas' and key not in joint_geometry:
    raise ValueError('Supply exact duration-free joint geometry for '+str(key))
   component=_shape_component(n,size,k,c,profile,gpu,
       geometry[n,size,k,c,profile.get('validate_range',False)],
       joint_geometry[key] if reduction=='jagwas' else None)
   if component.get('status')=='zero_available_capacity':return component
   components[size]=component
  component=components[size];portion=size/m
  if header_scenario is not None:
   row=census['chunks'][i];cpu=header_scenario['intervals'][i][int(header_scenario['endpoint']=='upper')]
   input_bytes=row['read_bytes'];decode_input_bytes=row['decode_input_bytes']
   base_update_bytes=row['native_ld_base_update_bytes'];replay_packed_bytes=row['native_ld_replay_packed_bytes']
  elif chunk_decoders is None:
   cpu=decoded_cpu*portion
   input_bytes=portion*census['record_payload_bytes']
   base_update_bytes=portion*decode['native_ld_base_update_bytes']
   decode_input_bytes=input_bytes;replay_packed_bytes=0
  else:
   exact=chunk_decoders[i]
   chunk_price=price_identified_units(exact['decoder'],profile['decode_units'])
   if chunk_price['unpriced_source_units']:raise ValueError('Unpriced chunk decoder units: '+str(chunk_price['unpriced_source_units']))
   cpu=chunk_price['identified_cpu_seconds']
   input_bytes=exact['census']['record_payload_bytes']+exact['decoder'].get('native_ld_replay_read_bytes',0)
   decode_input_bytes=exact['census']['record_payload_bytes']+exact['decoder'].get('native_ld_replay_decoded_input_bytes',0)
   replay_packed_bytes=exact['decoder'].get('native_ld_replay_packed_bytes',0)
   base_update_bytes=exact['decoder']['native_ld_base_update_bytes']
  # Each worker opens/parses its reader before issuing its first payload read.
  # Keeping this in post-read decode incorrectly changes initial overlap.
  ordinal=i+(header_scenario['issued_chunks'] if header_scenario else 0)
  reader_init_cpu=census.get('file_markers',m)*units['pgen_index_records'] if ordinal<workers else 0.
  # Packed staging write/read, int8 output with CPU write allocation, and
  # encoded input buffering. No claimed exact hardware transaction count.
  decoder_memory=_native_decoder_memory_bytes(n,size,input_bytes,decode_input_bytes,base_update_bytes,replay_packed_bytes,
      separate_buffered_read=profile.get('input_read_cpu_prices') is not None)
  read_cpu=0.;read_resources={'input':profile['read_bytes_per_second']}
  read_seconds=input_bytes/profile['read_bytes_per_second']
  if profile.get('input_read_cpu_prices') is not None:
   read=buffered_read_service(input_bytes,1,profile['input_read_cpu_prices'],
       storage_bytes_per_second=profile['read_bytes_per_second'],
       dram_bytes_per_second=profile['shared_dram_bytes_per_second'],cpu_fraction=effective_q)
   read_seconds=read['seconds'];read_cpu=read['cpu_seconds'];read_resources=read['resources']
   # The old two-access input-buffer approximation included writing the blob.
   # Its write is now in the read node. Decode retains one payload read.
  decode_seconds=max(cpu/effective_q,decoder_memory/profile['shared_dram_bytes_per_second'])
  executor=profile['executor_cpu_seconds']
  host_submit=component['host_dispatch_cpu_seconds']/effective_q
  if device_selection:
   from .native_status_work import native_status_service
   status=native_status_service(size,profile['device_status_service'],cpu_fraction=effective_q,
       dram_bytes_per_second=profile['shared_dram_bytes_per_second'])
   control=native_control_service(profile['control_primitives'],reused=ordinal>=depth,reduction=reduction)['cpu_seconds']
   block=dict(markers=size,variant_range=[start,end],decode_seconds=decode_seconds,decode_read_seconds=read_seconds,
       decode_resources={'cpu':cpu/decode_seconds,'dram':decoder_memory/decode_seconds},read_resources=read_resources,
       host_resources={'cpu':effective_q},h2d_bytes=n*size,d2h_bytes=size,
       h2d_seconds=n*size/profile['h2d_bytes_per_second'],operations=component['operations'],
       host_submit_seconds=host_submit,consumer_seconds=0.,result_contract='device_significant_status')
   for phase in ['decode_submit','publish','transfer_submit','result_submit','release','resolve']:
    block[phase+'_seconds']=control[phase]/effective_q
   for field in ['status_submit_seconds','status_submit_resources','d2h_seconds','finish_seconds','finish_operations','result_wait_resources']:
    block[field]=status[field]
   if ordinal<workers:block.update(reader_init_seconds=reader_init_cpu/effective_q,reader_init_resources={'cpu':effective_q})
   if 'host_dispatch_serial_cpu_seconds' in component:
    block['host_submit_serial_cpu_seconds']=component['host_dispatch_serial_cpu_seconds']
   if profile.get('handoff_wakeup_seconds') is not None:block['handoff_wakeup_seconds']=dict(profile['handoff_wakeup_seconds'])
   wait_fraction=profile.get('event_wait_cpu_fraction')
   if wait_fraction is not None:
    if isinstance(wait_fraction,bool) or not math.isfinite(wait_fraction) or not 0<=wait_fraction<=1:
     raise ValueError('Event-wait CPU fraction must be in [0,1]')
    block['event_wait_resources']={'cpu':effective_q*wait_fraction}
   blocks.append(block)
   cpu_total+=read_cpu+cpu+reader_init_cpu+component['host_dispatch_cpu_seconds']+sum(control.values())+status['cpu_seconds']
   result_work[size]=dict(ownership='owned',status=status,work=dict(array_bytes={},array_elements={},allocation_calls=0,copy_bytes=0),
       service=None,copy_threading=dict(unpriced_terms=[]),allocation_policy=dict(unpriced_terms=[]))
   continue
  result_bytes=17*size if reduction=='jagwas' else (4*k*(1+int(return_beta))+5)*size
  owned=owned_result_work(size,k,reduction=reduction,return_beta=return_beta)
  if borrowed:
   owned=dict(array_bytes={},array_elements={},allocation_calls=0,copy_bytes=0,borrowed_array_bytes=owned['array_bytes'])
  owned_service=None
  if not borrowed and profile.get('owned_result_copy_scenario') is not None:
   owned_service=owned_result_copy_service(owned,profile['owned_result_copy_scenario'],baseline_copy_bytes=baseline_bytes)
  threading=owned_result_copy_threading(owned,profile.get('owned_result_allocation_policy')) if not borrowed else dict(gil_released_by_array={},unpriced_terms=[])
  policy=owned_result_allocation_policy(owned,profile.get('owned_result_allocation_policy')) if not borrowed else dict(advice_compaction_exposure=False,unpriced_terms=[])
  result_work[size]=dict(ownership=ownership,work=owned,service=owned_service,copy_threading=threading,allocation_policy=policy)
  allocator=None
  if not borrowed and profile.get('owned_result_allocator') is not None:
   from .allocator_service import owned_allocator_service
   allocator_input=profile['owned_result_allocator']
   allocator=owned_allocator_service(owned,allocator_input['prices'],allocator_input['geometry'])
   result_work[size]['allocator']=allocator
  bulk_copy_bytes=0 if borrowed else max(0,result_bytes-baseline_bytes)
  bulk_copy_cpu=owned_service['additional_copy_cpu_seconds'] if owned_service else bulk_copy_bytes*units['numpy_copy_bytes']
  finish_baseline=finish_service['cpu_seconds'] if finish_service is not None else units['finish_fixed_calls']
  finish_cpu=finish_baseline+bulk_copy_cpu;finish_control_cpu=0.
  finish=finish_cpu/effective_q
  cpu_total+=read_cpu+cpu+reader_init_cpu+component['host_dispatch_cpu_seconds']+finish_cpu+3*executor+profile['stream_submission_cpu_seconds']
  blocks.append(dict(markers=size,variant_range=[start,end],decode_seconds=decode_seconds,decode_read_seconds=read_seconds,decode_resources={'cpu':cpu/decode_seconds,'dram':decoder_memory/decode_seconds},read_resources=read_resources,host_resources={'cpu':effective_q},decode_submit_seconds=executor/effective_q,
    h2d_bytes=n*size,d2h_bytes=result_bytes,h2d_seconds=n*size/profile['h2d_bytes_per_second'],transfer_submit_seconds=profile['stream_submission_cpu_seconds']/effective_q,
    operations=component['operations'],host_submit_seconds=host_submit,
    result_submit_seconds=2*executor/effective_q,d2h_seconds=result_bytes/profile['d2h_bytes_per_second'],finish_seconds=finish,consumer_seconds=0.))
  if ordinal<workers:
   blocks[-1].update(reader_init_seconds=reader_init_cpu/effective_q,
                     reader_init_resources={'cpu':effective_q})
  if 'host_dispatch_serial_cpu_seconds' in component:
   blocks[-1]['host_submit_serial_cpu_seconds']=component['host_dispatch_serial_cpu_seconds']
  if profile.get('control_primitives') is not None:
   control=native_control_service(profile['control_primitives'],reused=ordinal>=depth,reduction=reduction,return_beta=return_beta)['cpu_seconds']
   if finish_service is not None:
    # The exact finish probe already includes this mode's slices/numpy views.
    # Its dummy CUDA event was removed, so retain only the real ready event.
    control['finish']=profile['control_primitives']['event_synchronize_ready']
   finish_control_cpu=control['finish']
   cpu_total+=sum(control.values())-(3*executor+profile['stream_submission_cpu_seconds'])
   for phase in ['decode_submit','publish','transfer_submit','result_submit','release','resolve']:
    blocks[-1][phase+'_seconds']=control[phase]/effective_q
   blocks[-1]['finish_seconds']+=control['finish']/effective_q
  # Keep the fixed tiny/control service separate from numeric copy loops.
  # The tiny byte baseline is distributed proportionally across arrays; this
  # preserves the independently priced total without inventing allocation cost.
  phases=([dict(seconds=finish_control_cpu/effective_q),
           dict(seconds=finish_baseline/effective_q,host_serial_fraction=finish_service['serial_cpu_seconds']/finish_baseline if finish_baseline else 0.)]
          if finish_service is not None else [dict(seconds=max(0.,blocks[-1]['finish_seconds']-bulk_copy_cpu/effective_q))])
  for name,amount in owned['array_bytes'].items():
   if allocator and allocator['arrays'][name]['allocate_cpu_seconds']:
    phases.append(dict(seconds=allocator['arrays'][name]['allocate_cpu_seconds']/effective_q,host_serial_fraction=1.))
   portion=amount/owned['copy_bytes'];seconds=bulk_copy_cpu*portion/effective_q
   released=threading['gil_released_by_array'][name]
   phases.append(dict(seconds=seconds,host_serial_fraction=None if released is None else (0. if released else 1.),
                      resources={'cpu':effective_q,'dram':2*bulk_copy_bytes*portion/seconds if seconds else 0.}))
  if allocator:
   extra=allocator['allocate_cpu_seconds']+allocator['worker_release_cpu_seconds']
   if allocator['worker_release_cpu_seconds']:
    phases.append(dict(seconds=allocator['worker_release_cpu_seconds']/effective_q,host_serial_fraction=1.))
   blocks[-1]['finish_seconds']+=extra/effective_q
   blocks[-1]['discard_seconds']=allocator['consumer_release_cpu_seconds']/effective_q
   cpu_total+=extra+allocator['consumer_release_cpu_seconds']
  blocks[-1]['finish_operations']=phases
  if profile.get('handoff_wakeup_seconds') is not None:
   blocks[-1]['handoff_wakeup_seconds']=dict(profile['handoff_wakeup_seconds'])
  wait_fraction=profile.get('event_wait_cpu_fraction')
  if wait_fraction is not None:
   if isinstance(wait_fraction,bool) or not math.isfinite(wait_fraction) or not 0<=wait_fraction<=1:
    raise ValueError('Event-wait CPU fraction must be in [0,1]')
   blocks[-1]['event_wait_resources']={'cpu':effective_q*wait_fraction}
 allocator_terms={term for result in result_work.values() for term in result.get('allocator',{}).get('unpriced_terms',[])}
 if device_selection:
  allocator_terms.update(term for result in result_work.values() for term in result['status']['unpriced_terms'])
  allocator_terms.add('Device-selected tensor operations, count barriers and owned payloads require the indexed selection graph')
 if finish_service is not None:allocator_terms.add('fixed-32-row status/QC and returned-view metadata destruction timing beyond the finish probe')
 workspace=dict(packed_decoder_bytes_upper=min(workers,len(blocks))*min(b,m)*((n+3)//4),
     packed_decoder_allocation_calls_upper=min(workers,len(blocks)),
     scope='One grow-only packed workspace per active native reader, retained until reader close. Excludes record buffers, metadata and LD replay; initial page faults and close release remain unpriced.')
 if decode.get('native_ld_replay_records'):
  from .decoder_work import native_reader_workspace
  if header_scenario is None:_,extra=native_reader_workspace(census['chunks'])
  else:extra=max(((n+3)//4)+8*(row['variant_range'][0]-row['ld_replay']['base_variant'])
      for row in census['chunks'] if row['ld_replay'] is not None)
  workspace['ld_replay_workspace_bytes_upper']=min(workers,len(blocks))*extra
  workspace['scope']='Grow-only packed workspace plus one LD-base scratch row and prefix offsets per active reader. Record buffers and ordinary metadata require separate memory accounting; initial page faults and close release remain unpriced.'
  allocator_terms.add('LD-base restart Python/ctypes dispatch, scratch allocation and backward index search')
 return dict(blocks=blocks,components=components,decoder=decode,encoded_work_distribution=distribution,cpu_work_seconds=cpu_total,result_ownership=ownership,host_workspace=workspace,allocator_unpriced_terms=sorted(allocator_terms),
             workers=workers,depth=depth,owned_result_work=result_work,unpriced_terms=sorted({term for component in components.values() for term in component.get('unpriced_terms',[])}|{term for result in result_work.values() for part in ['allocation_policy','copy_threading'] for term in result[part]['unpriced_terms']}|({'CUDA completion-event waiting CPU demand'} if profile.get('event_wait_cpu_fraction') is None else set())|set(profile.get('decode_pricing_unpriced_terms',['Native decoder primitive/compiler equivalence']))|({'input-read destination allocation/faults, short-read retries and page-cache fill traffic'} if profile.get('input_read_cpu_prices') is not None else {'input-read syscall/copy CPU service'})|({'handoff parking/resume CPU work and loaded-context latency'} if profile.get('handoff_wakeup_seconds') is not None else {'blocked queue/future handoff wakeup latency'})|({'owned-result allocation/free service and actual page-reuse state'} if not (borrowed or device_selection) else set())|{'completion-event wakeup and worker-start latency'}|set(storage_terms)))


def torch_scan_runtime(data,profile):
 """Finite scan-to-discard candidate. No preprocessing or durable output."""
 work=torch_scan_work(data,profile)
 if work.get('status')=='zero_available_capacity':return work
 schedule=torch_scan_schedule(work['blocks'],depth=work['depth'],decode_workers=work['workers'],
     shared_capacities={'cpu':profile['cpu_available_cores'],'dram':profile['shared_dram_bytes_per_second'],
                        'input':profile['read_bytes_per_second']})
 scope=('Native int8 eager scan with prepared JAGWAS factor to owned narrow-result discard; excludes shared preprocessing, per-device factor preparation, p-values and indexed output.' if profile.get('reduction')=='jagwas' else
        'Native int8 eager scan to '+work['result_ownership']+' result discard; excludes setup, reduction, p-values and durable output. Independently supplied resources and exact N/B/K/C launch geometry; no shape interpolation.')
 return dict(status='development_scan_candidate',estimated_scan_seconds=schedule['seconds'],
     prediction_complete=False,unpriced_terms=work['unpriced_terms']+work['allocator_unpriced_terms'],
     components={str(size):{key:component[key] for key in ['estimated_span_seconds','host_dispatch_cpu_seconds','kernel_service_seconds','gemm']} for size,component in work['components'].items()},
     source_work=dict(decoder=work['decoder'],host_workspace=work['host_workspace'],encoded_work_distribution=work['encoded_work_distribution'],cpu_work_seconds=work['cpu_work_seconds'],wait_resource_seconds=schedule.get('wait_resource_seconds',{}),handoffs=handoff_summary(schedule),owned_result_work=work['owned_result_work']),
     schedule=dict(blocks=len(work['blocks']),depth=work['depth'],decode_workers=work['workers']),
     result_ownership=work['result_ownership'],
     scope=scope,
     assumptions=['Encoded input work uses '+work['encoded_work_distribution']+'; native LD restarts require exact base records and contiguous read prefixes.',
                  'Proportional fluid resource sharing is a scheduling scenario.',
                  'Source-derived traffic/cache/issue approximations are not validated runtime guarantees.'])



def torch_multigpu_scan_runtime(shards, shared_capacities, shared_links=(), ordered=False, host_serial_fraction=None, queue_service=None,host_serial_policy='fluid', result_queue_depth=4):
 """Price caller-supplied per-device workloads; never invent shard LD censuses.

 Each shard is {device, data, profile}. Setup/trait-block transitions are outside
 this simultaneous-scan boundary. Shared capacities are aggregate, not per GPU.
 """
 from .execution_graph import torch_multigpu_schedule
 graphs=[];unpriced=set();reports=[];ownerships=set();ack_services=[]
 for shard in shards:
  work=torch_scan_work(shard['data'],shard['profile'])
  if work.get('status')=='zero_available_capacity':return work
  ownerships.add(work['result_ownership'])
  ack_services.append(shard['profile'].get('acknowledgement_service'))
  graphs.append(dict(device=shard['device'],blocks=work['blocks'],depth=work['depth'],decode_workers=work['workers']))
  unpriced.update(work['unpriced_terms'])
  unpriced.update(work['allocator_unpriced_terms'])
  reports.append(dict(device=shard['device'],blocks=len(work['blocks']),host_workspace=work['host_workspace'],encoded_work_distribution=work['encoded_work_distribution'],cpu_work_seconds=work['cpu_work_seconds'],owned_result_work=work['owned_result_work']))
 if len(ownerships)!=1:raise ValueError('All device shards must use the same result ownership')
 borrowed=ownerships=={'borrowed'}
 if borrowed and len(shards)>1 and (ack_services[0] is None or any(s!=ack_services[0] for s in ack_services)):
  raise ValueError('Consistent independent acknowledgement service required')
 schedule=torch_multigpu_schedule(graphs,shared_capacities,shared_links,ordered=ordered,result_queue_depth=result_queue_depth,host_serial_fraction=host_serial_fraction,queue_service=queue_service,host_serial_policy=host_serial_policy,
     borrow_results=borrowed,acknowledgement_service=ack_services[0] if borrowed else None)
 if borrowed and len(shards)>1:unpriced.add('acknowledgement parking/resume CPU placement, timeout retries and loaded-context transfer')
 unpriced.add('GIL transition order, unobserved host paths and driver contention remain unresolved')
 if len(shards)>1:unpriced.add('result queue wakeup, timeout, sentinel and thread lifecycle service')
 if len(shards)>1 and queue_service is None:unpriced.add('multiGPU result queue put/get CPU service')
 return dict(status='development_multigpu_scan_candidate',estimated_scan_seconds=schedule['seconds'],
             prediction_complete=False,result_ownership=next(iter(ownerships)),unpriced_terms=sorted(unpriced),shards=reports,wait_resource_seconds=schedule.get('wait_resource_seconds',{}),handoffs=handoff_summary(schedule),
             scope=schedule['scope'],resource_policy=schedule.get('resource_policy'))







def handoff_summary(schedule):
 """Summed conditional service across pipelines, not added elapsed wall time."""
 rows=list(schedule.get('conditional_delays',{}).values())
 return dict(possible_waits=len(rows),blocked_waits=sum(row['blocked'] for row in rows),
             extra_elapsed_service_seconds=sum(row['extra_elapsed_service_seconds'] for row in rows),
             scope='Sum across handoffs; overlapping delays must not be added to scan wall time.')
