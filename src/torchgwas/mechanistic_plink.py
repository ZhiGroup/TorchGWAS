"""PLINK2 K=1 source-work/runtime candidate; no GWAS timing coefficients."""
from __future__ import annotations
import math
from .first_principles import positive,service,plink_missing_branches
from .tensor_service import IntervalLRU
from .decoder_work import decoder_work
from .execution_graph import plink_block_schedule

def plink_variant_trace(n,c,observed,branch,profile):
 """Ideal source-array reuse with explicit arithmetic and scalar-loop work.

 The packed reader has a separate work census. BLAS cache traffic uses a
 unique-operand scan approximation; internal packing/tile reloads are not
 decoded from proprietary MKL machine code. No N-specific throughput is used.
 """
 if branch not in ['gram','opener','sparse']:raise ValueError('Unknown PLINK branch')
 p=c+2;u=profile['plink_units'];nobs=observed if branch=='gram' else n
 # Each source vector occupies cache lines. Expected missingness can produce
 # a fractional nobs; line rounding applies to its expected workspace size.
 sizes={'source':8*n*(c+1),'design':8*nobs*p,'phenotype':8*nobs,
        'packed':math.ceil(n/4),'mask':math.ceil(n/8),'small':8*(5*p*p+4*p)}
 layout={};addr=0
 for key,size in sizes.items():layout[key]=(addr,addr+math.ceil(size/64));addr+=math.ceil(size/64)+1
 caches={k:IntervalLRU(int(profile['cache_bytes'][k])//64) for k in ['L1','L2','L3']}
 # L3 is an explicit per-worker capacity scenario set by the profile builder.
 stages=[]
 def add(name,reads,writes,instruction=0.,flops=0.):stages.append((name,reads,writes,instruction,flops))
 if branch in ['gram','opener']:
  if branch=='gram':add('nonmissing_mask',['packed'],['mask'],math.ceil(n/128)*u['mask128'])
  add('expand',['packed'],['design'],math.ceil(n/128)*u['expand128'])
  add('gather',['source','mask'],['design','phenotype'],(c+1)*nobs*u['gather126']/126)
  if branch=='gram':
   add('gram',['design'],['small'],flops=nobs*p*(p+1))
   add('vif',['small'],['small'],u['vif'])
   add('xty',['design','phenotype'],['small'],flops=2*nobs*p)
   add('solve',['small'],['small'],u['solve'])
  else:
   add('crossproducts',['design','phenotype'],['small'],flops=2*n*(c+1))
 else:
  # A minor-allele carrier supplies one genotype*c covariate/phenotype FMA.
  carriers=n*profile['carrier_fraction']
  add('sparse_carriers',['packed','source'],['small'],flops=2*carriers*(c+1))
 if branch!='gram':
  add('vif_update',['small'],['small'],u['vif_nm'])
  add('inverse_update',['small'],['small'],u['rank1'])
  add('coefficient_product',['small'],['small'],flops=2*p*p)
 # RSS dot, triangle scaling, coefficient/SE validity and sqrt. Arithmetic
 # equivalent work is explicit; branch/div/sqrt dependency costs approximate.
 add('post_regression',['small'],['small'],flops=2*p+p*(p+1)+6*p)
 results=[]
 for repetition in range(4):
  rows=[]
  for name,reads,writes,instruction,flops in stages:
   # Expansion writes only genotype column. Gather writes C design columns
   # and one phenotype, not the already-expanded genotype and intercept.
   def extent(key,write):
    lo,hi=layout[key]
    if write and key=='design' and name=='expand':hi=lo+math.ceil(8*nobs/64)
    elif write and key=='design' and name=='gather':lo+=math.ceil(16*nobs/64)
    return lo,hi
   logical=0.;below={k:0 for k in caches}
   for writing,keys in [(False,reads),(True,writes)]:
    for key in keys:
     lo,hi=extent(key,writing);logical+=(hi-lo)*64
     for level,cache in caches.items():
      row=cache.access(lo,hi,write=writing,write_allocate=writing);below[level]+=64*(row['read_fill']+row['writeback'])
   traffic=dict(L1=logical,L2=below['L1'],L3=below['L2'],DRAM=below['L3'])
   mem=max(service(value,profile['bandwidth_per_cpu_second'][level]) for level,value in traffic.items())
   arithmetic_rate=profile['fp64_flops_per_cpu_second']
   library={'gram':'syrk','xty':'vector_gemm','crossproducts':'gemv'}.get(name)
   if library and 'blas_flops_per_cpu_second' in profile:arithmetic_rate=profile['blas_flops_per_cpu_second'][library]
   arithmetic=service(flops,arithmetic_rate)
   rows.append(dict(stage=name,cpu_seconds=max(instruction,mem,arithmetic),instruction_cpu_seconds=instruction,flops=flops,traffic_bytes=traffic))
  results.append(dict(cpu_seconds=sum(r['cpu_seconds'] for r in rows),dram_bytes=sum(r['traffic_bytes']['DRAM'] for r in rows),stages=rows))
 return dict(first=results[0],steady=results[-1],branch=branch,observed_samples=observed)

def plink_runtime(data,profile):
 n,m,c=(data[k] for k in ['samples','markers','covariates']);p=c+2;u=profile['plink_units'];q=profile['cpu_fraction']
 if data['traits_analyzed']!=1 or not data['matching_sample_order']:raise ValueError('Unpriced cohort branch')
 if p!=profile['small_matrix_dimension']:raise ValueError('Require independent small-matrix primitive at this covariate dimension')
 for key in ['cpu_fraction','cpu_available_cores','read_bytes_per_second','write_bytes_per_second','shared_dram_bytes_per_second']:
  if positive(key,profile[key],True)==0:return dict(estimated_seconds=None,status='zero_available_capacity',resource=key)
 if q>1:raise ValueError('Scheduling share exceeds one')
 threads=min(profile['compute_workers'],m);workers=min(threads*q,profile['cpu_available_cores']);block_size=profile['block_markers']
 rate=profile['missing_rate'];segments=sum(min(threads,min(block_size,m-lo)) for lo in range(0,m,block_size))
 branches=plink_missing_branches(n,m,rate,segments);gram=branches['gram'];opener=branches['opener'];sparse=branches['sparse']
 observed_sum=n*(gram-m*rate);observed=observed_sum/gram if gram else n
 if observed<=p:raise ValueError('Scenario includes invalid regression degrees of freedom')
 placements=profile.get('worker_cache_bytes',[profile['cache_bytes']]*threads)
 if len(placements)!=threads:raise ValueError('Cache placement must specify every compute worker')
 worker_traces=[];trace_cache={}
 for cache in placements:
  key=tuple((level,int(cache[level])) for level in ['L1','L2','L3'])
  if key not in trace_cache:
   local=dict(profile,cache_bytes=cache)
   trace_cache[key]={branch:plink_variant_trace(n,c,observed,branch,local) for branch in ['gram','opener','sparse']}
  worker_traces.append(trace_cache[key])
 # Retain one full trace per worker: sources are shared within each L3 domain,
 # but the equal capacity partition below is a declared conservative-capacity
 # scenario and still counts shared operand residency separately per worker.
 traces=worker_traces[0]
 tables=data['tables'];metadata_bytes=sum(t['bytes'] for t in tables.values())
 # Source scans full input rows but numeric-converts only selected phenotype
 # columns. Positional integer fields and identifier retention are separate.
 token_characters=sum(t['field_characters']+t['fields'] for t in tables.values())
 parse=token_characters*u['token16']/16+n*(c+1)*u['scan_double']+m*u['scan_uint']
 retained=64*sum(t['fields'] for t in tables.values())
 metadata_memory=service(retained,profile['bandwidth_per_cpu_second']['DRAM'])
 # Ordered IDs still require dictionary/key comparison work; model the byte
 # scanning/copy passes explicitly, not as a measured whole-process residual.
 matching_chars=data['tables']['psam']['field_characters']*3
 matching=matching_chars*u['token16']/16
 setup_flops=2*n*(c+1)*(c+2)+2*n*(c+1)+2*(c+1)**3
 setup=profile['process_launch_seconds']+service(metadata_bytes,profile['read_bytes_per_second'])+(parse+metadata_memory+matching+service(setup_flops,profile['fp64_flops_per_cpu_second']))/q
 encoded=data['encoded'];decoder=decoder_work(encoded,'pgenlib_avx2_bmi2',restart_ld_bases=True)
 if decoder['uncounted_mechanisms']:raise ValueError('Uncounted PGEN replay/record work: '+str(decoder['uncounted_mechanisms']))
 counts=decoder['source_units'];decoder_terms={}
 for name,amount in counts.items():
  if name in ['copy_packed_byte','fill_packed_byte','invert_packed_byte']:continue
  if name not in profile['decode_units']:raise ValueError('Missing decoder primitive '+name)
  decoder_terms[name]=amount*profile['decode_units'][name]
 packed=math.ceil(n/4)*m
 decode_bytes=encoded['record_payload_bytes']+packed+2*counts.get('copy_packed_byte',0)+counts.get('fill_packed_byte',0)+2*counts.get('invert_packed_byte',0)
 decode_cpu=max(sum(decoder_terms.values()),service(decode_bytes,profile['bandwidth_per_cpu_second']['L2']))
 worker_cpu=[]
 for wt in worker_traces:
  demand=sum(branches[b]*wt[b]['steady']['cpu_seconds'] for b in wt)/threads+decode_cpu/threads
  demand+=max(0.,wt['gram']['first']['cpu_seconds']-wt['gram']['steady']['cpu_seconds'])*segments/threads
  worker_cpu.append(demand)
 compute_cpu=sum(worker_cpu)
 compute_wall=max(max(worker_cpu)/q,compute_cpu/profile['cpu_available_cores'])
 # Current --glm hide-covar default: A1_FREQ,BETA,SE,T_STAT and P, plus POS/N.
 # Source output header is verified by the profile builder.
 format_cpu=m*(4*u['dtoa']+u['ln_format']+u['t_tail']+2*u['u32toa'])
 identifier_bytes=data['tables']['pvar']['field_characters']/m
 output_bytes=m*(identifier_bytes+4*12+14+len(str(n))+15)
 output=service(output_bytes,profile['write_bytes_per_second'])
 blocks=[]
 for lo in range(0,m,block_size):
  fraction=min(block_size,m-lo)/m
  blocks.append(dict(read_seconds=encoded['file_bytes']*fraction/profile['read_bytes_per_second'],
   compute_seconds=compute_wall*fraction,format_write_seconds=format_cpu*fraction/q+output*fraction))
 finite=plink_block_schedule(blocks)
 dram=decode_bytes+sum(sum(branches[b]*wt[b]['steady']['dram_bytes'] for b in wt) for wt in worker_traces)/threads
 cpu_floor=(compute_cpu+format_cpu)/profile['cpu_available_cores'];dram_floor=dram/profile['shared_dram_bytes_per_second']
 scan=max(finite['seconds'],cpu_floor,dram_floor)
 return dict(estimated_seconds=setup+scan,status='full_major_stage_candidate_with_approximations',setup_seconds=setup,scan_seconds=scan,
  stage_seconds=dict(process_launch=profile['process_launch_seconds'],metadata_read=metadata_bytes/profile['read_bytes_per_second'],metadata_parse=parse/q,
   metadata_memory=metadata_memory/q,sample_match=matching/q,setup_arithmetic=setup_flops/profile['fp64_flops_per_cpu_second']/q,
   genotype_read=encoded['file_bytes']/profile['read_bytes_per_second'],decoder_pool=decode_cpu/workers,
   regression_pool=(compute_cpu-decode_cpu)/workers,format_main=format_cpu/q,buffered_output=output),
  branch_counts=branches,gram_observed_sample_sum=observed_sum,traces=traces,worker_cpu_seconds=worker_cpu,worker_cache_bytes=placements,decoder_cpu_terms=decoder_terms,
  worker_segments=segments,workers=threads,queue=dict(seconds=finite['seconds'],blocks=len(blocks)),
  shared_resource_seconds=dict(cpu=cpu_floor,dram=dram_floor),output_bytes_scenario=output_bytes,
  scope='K1, additive autosomal biallelic hardcalls, C8, same-order complete phenotype/covariates and all variants valid. Source work plus generic capacities, no GWAS runtime fitting.',
  assumptions=[
   'Missingness and carrier fractions are independent data-statistic parameters. Expected iid branch counts include every worker/block restart; actual clustered missingness needs an exact branch census.',
   'Fixed generic 64x512 library throughput, when supplied, prices the matching stored-source BLAS wrapper; this is an independent attained-capacity approximation, not an N-sized timing table or proof of rank-independent library efficiency.',
   'Generic scalar and 10x10 library primitives are fixed independent work. T-tail uses a generic t/df distribution; data-dependent continued-fraction iterations can differ.',
   'Fully associative cache and unique BLAS operand scans approximate actual MKL packing/tile traffic and instruction issue. Cache sharing uses explicit worker placement and an equal capacity partition within each physical domain; common input residency is conservatively duplicated by this partition approximation.',
   'Text scanning uses fixed-length token character service and actual input byte/field counts; hash-table and metadata bookkeeping use a scanning/copy approximation, not exact compiled instruction simulation.',
   'Decoder record header/group/category instructions use generic ULEB1 service proxies, recorded in the profile. No unknown LD base replay is silently omitted.',
   'Within-block worker work is balanced in expectation. Scalar branch/error handling and OS/thread creation/teardown are not fully instruction-resolved.',
   'Native output close, without a new common fsync. Read/format share the main thread with finite compute overlap.',
   'Independent resource measurements describe their observed load, not a guarantee of future capacity. This candidate is not a validated crossover or bound.'
  ])
