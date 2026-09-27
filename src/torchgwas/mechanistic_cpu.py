"""Source-stage CPU estimates and finite queues. No GWAS timing coefficients."""
from __future__ import annotations
import math
from .first_principles import positive,service
from .tensor_service import IntervalLRU


def reader_analysis_schedule(reader,analysis,depth=3):
 """One sequential producer and consumer; a slot is held through analysis.

 Durations must already include explicit resource availability. A block's slot
 cannot be reused until that block's analysis/output finishes. Tail blocks are
 real blocks, not a fractional startup or free pipeline drain.
 """
 if len(reader)!=len(analysis) or not isinstance(depth,int) or depth<1:raise ValueError('Invalid queue')
 end_reader=0.;end_main=0.;rows=[]
 for i,(read,main) in enumerate(zip(reader,analysis)):
  positive('reader duration',read,True);positive('analysis duration',main,True)
  start=max(end_reader,rows[i-depth]['analysis_end'] if i>=depth else 0.)
  end_reader=start+read;begin=max(end_main,end_reader);end_main=begin+main
  rows.append(dict(reader_start=start,reader_end=end_reader,analysis_start=begin,analysis_end=end_main))
 return dict(seconds=end_main,blocks=len(rows),first_block=rows[:1],last_block=rows[-1:])


def fastgwa_variant_trace(n,p,units,profile,iterations=4):
 """Source arrays plus Eigen 3.4 alias-safe x - X*(H*x) temporary.

 Each cache is an ideal fully-associative LRU stack, tracking the same logical
 references. CPU stores write-allocate. Operand-major order approximates actual
 tiled/interleaved instructions. This is an explicit cache model, not an N grid.
 """
 line=64;names=['L1','L2','L3'];caches=[IntervalLRU(int(profile['cache_bytes'][k])//line) for k in names]
 sizes={'H':8*n*p,'X':8*n*p,'y':8*n,'temporary_y':8*n,'Hy':8*p,'phenotype':8*n}
 layout={};addr=0
 for name,size in sizes.items():layout[name]=(addr,addr+math.ceil(size/line));addr+=math.ceil(size/line)+1
 packed=(n+3)//4
 # Allocation/free bookkeeping only. First touch and zero fill are separate.
 def malloc_cost(size):
  bins=sorted((int(k),v) for k,v in profile['malloc_cpu_seconds'].items())
  found=next((cost for ceiling,cost in bins if ceiling>=size),None)
  if found is None:raise ValueError('Allocator primitive size domain exceeded')
  return found
 count_words=(n+31)//32;full60,remain=divmod(count_words,120);full6,tail=divmod(remain,12)
 count_cpu=full60*units['count60']+full6*units['count6']+tail*units['packed_count_word']
 stages=[('count',['packed'],[],count_cpu,0),
   ('zero_genotype',[],['y'],malloc_cost(8*n),0),
   ('expand',['packed'],['y'],((n+1)//2)*units['expand2_double'],0),
   ('alias_safe_copy',['y'],['temporary_y'],malloc_cost(8*n),0),
   ('zero_Hy',[],['Hy'],malloc_cost(8*p),0),
   ('gemv_H',['H','y','Hy'],['Hy'],0,2*n*p),
   ('gemv_X_subtract',['X','Hy','temporary_y'],['temporary_y'],0,2*n*p+2*n),
   ('copy_result',['temporary_y'],['y'],0,0),
   ('dot_xx',['y'],[],0,2*n),
   ('dot_xy',['y','phenotype'],[],0,2*n)]
 output=[]
 for variant in range(iterations):
  layout['packed']=(addr+variant*(math.ceil(packed/line)+1),addr+variant*(math.ceil(packed/line)+1)+math.ceil(packed/line))
  costs=[]
  for name,reads,writes,instruction_cpu,flops in stages:
   logical=sum((packed if k=='packed' else sizes[k]) for k in reads+writes)
   below=[0.,0.,0.]
   for ci,cache in enumerate(caches):
    for key in reads:
     v=cache.access(*layout[key]);below[ci]+=line*(v['read_fill']+v['writeback'])
    for key in writes:
     # Explicit read-for-ownership for write misses in CPU caches.
     v=cache.access(*layout[key],write=True,write_allocate=True);below[ci]+=line*(v['read_fill']+v['writeback'])
   traffic=dict(L1=logical,L2=below[0],L3=below[1],DRAM=below[2])
   mem=max(service(traffic[k],profile['bandwidth_per_cpu_second'][k]) for k in traffic)
   arithmetic=service(flops,profile['fp64_flops_per_cpu_second'])
   if name=='gemv_H' and 'vector_mac_cpu_seconds' in profile:
    # Eigen 3.4 GeneralMatrixVector.h: 8/4/3/2/1 packets, then scalar tail.
    # SSE2 has two FP64 lanes. The multiply operands are independent of the
    # accumulator; each output packet is one addition recurrence.
    remaining=p;arithmetic=0.
    for width in [16,8,6,4,2,1]:
     groups=remaining//width;remaining%=width
     if not groups:continue
     chain_rate=(2/profile['scalar_mac_cpu_seconds'] if width==1 else (width//2)*4/profile['vector_mac_cpu_seconds'])
     arithmetic+=service(2*n*width*groups,min(profile['fp64_flops_per_cpu_second'],chain_rate))
   cpu=max(instruction_cpu,arithmetic,mem)
   costs.append(dict(stage=name,cpu_seconds=cpu,traffic_bytes=traffic,logical_bytes=logical,flops=flops,arithmetic_cpu_seconds=arithmetic,memory_cpu_seconds=mem))
  output.append(dict(cpu_seconds=sum(x['cpu_seconds'] for x in costs),dram_bytes=sum(x['traffic_bytes']['DRAM'] for x in costs),stages=costs))
 if output[-1]['cpu_seconds']!=output[-2]['cpu_seconds'] or output[-1]['dram_bytes']!=output[-2]['dram_bytes']:
  raise ValueError('Cache trace did not reach a stationary repeating cost; increase trace length')
 return dict(first=output[0],steady=output[-1],source='GCTA 1.95.3 FastFAM/Geno; Eigen 3.4 ProductEvaluators alias rules',
  scope='Operand-major LRU and generic capacities; not exact compiled CPU instruction or set-associative cache simulation.')


def fastgwa_runtime(data,profile):
 """Full-stage candidate: exact input work, independent primitives, finite queue.

 Covers the current K=1, same-order, autosomal hardcall, all-variants-valid
 benchmark branch. This is an approximation with explicit scope, not a bound.
 Profile capacities and load are inputs; no elapsed GWAS data are accepted.
 """
 from .decoder_work import decoder_work
 n,m,c=data['samples'],data['markers'],data['covariates'];p=c+1;u=profile['units']
 if data['traits_analyzed']!=1 or not data['matching_sample_order']:raise ValueError('Unpriced analysis/matching branch')
 for key in ['cpu_fraction','cpu_available_cores','read_bytes_per_second','write_bytes_per_second','shared_dram_bytes_per_second']:
  if positive(key,profile[key],True)==0:return dict(estimated_seconds=None,status='zero_available_capacity',resource=key)
 q=profile['cpu_fraction']
 if q>1:raise ValueError('CPU scheduling fraction exceeds one')
 # A profile defines the available capacity after background contention.
 # Main and reader are each serial even when --thread-num is greater than one.
 trace=fastgwa_variant_trace(n,p,u,profile)
 tables=data['tables'];text_cpu=0.;copy_cpu=0.
 for key,t in tables.items():
  mode='split'+str(t['columns'])+('long' if t['long_fields']>t['fields']/4 else '')
  if mode not in u:raise ValueError('Missing independent parser primitive: '+mode)
  text_cpu+=t['rows']*u[mode]
  # readTxtList's persistent column copies; pheno/covar parsing also retains ids.
  short=t['fields']-t['long_fields']
  copy_cpu+=short*u['string_short']+t['long_fields']*u['string_construct']
 # Marker copies ID/ref/alt; numeric input parsing reads every phenotype column,
 # even when --mpheno chooses just one trait.
 copy_cpu+=3*m*u['string_short']
 numeric=n*(data['traits_in_file']+c)*u['strtod']+2*m*u['stoi']
 sort=data['sort_counts']
 matching=2*sort['string_sort_compares']*u['string_compare']
 # Three vector_commonIndex_sorted1 calls still sort the integer index array;
 # the ordered-id branch avoids the two string-index sorts inside each call.
 matching+=3*sort['sorted_integer_sort_compares']*u['string_compare']
 matching+=(2*(n-1)+3*n)*u['string_compare']
 setup_flops=4*n*p*p+(2/3)*p**3+4*n*p
 table_memory=64*sum(t['fields'] for t in tables.values())+16*n*p
 setup_cpu=text_cpu+copy_cpu+numeric+matching+max(setup_flops/profile['fp64_flops_per_cpu_second'],table_memory/profile['bandwidth_per_cpu_second']['DRAM'])
 setup_read=sum(t['bytes'] for t in tables.values())/profile['read_bytes_per_second']
 setup=profile['process_launch_seconds']+setup_read+setup_cpu/q
 encoded=data['encoded'];work=decoder_work(encoded,'pgenlib_sse2');counts=work['source_units']
 decode_terms={}
 mapping={'onebit128_sse2':'onebit128_sse2','onebit32_pgenlib_sse2':'onebit32_pgenlib','set_category':'set_category'}
 for length in range(1,6):mapping['uleb'+str(length)]='uleb'+str(length)
 for name,unit in mapping.items():
  if counts.get(name,0):
   if unit not in u:raise ValueError('Missing decoder unit '+unit)
   decode_terms[name]=counts[name]*u[unit]
 # Header/group/category extraction are explicit instruction proxies, never
 # a decoder residual fitted from full records. Their uncertainty is reported.
 decode_terms['group_absolute_load']=counts.get('difflist_group_absolute_id',0)*u['uleb1']
 decode_terms['category_extract']=counts.get('difflist_category_extract',0)*u['uleb1']
 decode_terms['record_dispatch']=m*u['uleb1']
 decode_terms['invert']=math.ceil(counts.get('invert_packed_byte',0)/16)*u['invert16']
 decode_memory=2*counts.get('copy_packed_byte',0)+counts.get('fill_packed_byte',0)+2*counts.get('invert_packed_byte',0)+encoded['record_payload_bytes']+m*((n+3)//4)
 decode_cpu=max(sum(decode_terms.values()),decode_memory/profile['bandwidth_per_cpu_second']['L2'])
 read_wall=encoded['file_bytes']/profile['read_bytes_per_second']+decode_cpu/q
 # Four floating-point fields, one count, and two marker integer conversions.
 format_cpu=m*(4*u['ostream_double']+u['ostream_uint']+2*u['to_string_uint']+8*u['string_short'])
 tail_cpu=m*u['chi1']
 # General-format text uses at most precision+6 chars per floating field in
 # the currently observed exponent range. This is an explicit format scenario.
 marker_chars=tables['pvar']['field_characters']/m
 output_bytes=m*(marker_chars+len(str(n))+4*(6+6)+10)
 output_write=output_bytes/profile['write_bytes_per_second']
 math_cpu=trace['first']['cpu_seconds']+(m-1)*trace['steady']['cpu_seconds']
 analysis_wall=(math_cpu+format_cpu+tail_cpu)/q+output_write
 reader=[];analysis=[]
 for low in range(0,m,1024):
  length=min(1024,m-low);fraction=length/m
  reader.append(read_wall*fraction)
  # Preserve the first-variant cache cost in the first block only.
  block_math=length*trace['steady']['cpu_seconds']+(trace['first']['cpu_seconds']-trace['steady']['cpu_seconds'] if low==0 else 0.)
  analysis.append((block_math+(format_cpu+tail_cpu)*fraction)/q+output_write*fraction)
 schedule=reader_analysis_schedule(reader,analysis,3)
 shared_cpu=(decode_cpu+math_cpu+format_cpu+tail_cpu)/profile['cpu_available_cores']
 dram=decode_memory+trace['first']['dram_bytes']+(m-1)*trace['steady']['dram_bytes']
 shared_dram=dram/profile['shared_dram_bytes_per_second']
 scan=max(schedule['seconds'],shared_cpu,shared_dram)
 return dict(estimated_seconds=setup+scan,status='full_stage_candidate_with_explicit_approximations',
  setup_seconds=setup,scan_seconds=scan,reader_seconds=read_wall,analysis_seconds=analysis_wall,
  stage_seconds=dict(process_launch=profile['process_launch_seconds'],metadata_read=setup_read,
    text_parse=text_cpu/q,retained_strings=copy_cpu/q,numeric_parse=numeric/q,sample_matching=matching/q,
    other_setup=(setup_cpu-text_cpu-copy_cpu-numeric-matching)/q,
    genotype_read=encoded['file_bytes']/profile['read_bytes_per_second'],decoder=decode_cpu/q,
    analysis_math=math_cpu/q,chi_square=tail_cpu/q,formatting=format_cpu/q,buffered_write=output_write),
  decoder_instruction_cpu_seconds=decode_terms,variant_trace=trace,queue=schedule,
  shared_resource_floor_seconds=dict(cpu=shared_cpu,dram=shared_dram),output_bytes_scenario=output_bytes,
  assumptions=[
   'Exact encoded/text counts for this input; all samples remain in source order and all autosomal variants pass QC.',
   'One reader plus one analysis worker in the official GCTA 1.95.3 build; no four-thread division.',
   'Generic fixed-row Boost parsing represents short and long field classes; allocation/cache effects of large metadata tables are approximated.',
   'Integer-sort comparisons use the generic string-comparison cost as an explicit service proxy, not an exact instruction match.',
   'Group header loads, category extraction and record dispatch use ULEB1 primitive service proxies; compiler and branch costs are approximate.',
   'Packed scalar tails, rare error paths and lower-order scalar bookkeeping are not fully instruction-resolved.',
   'Fully associative operand-major CPU cache trace, write allocation and reuse of freed vector storage; actual associativity and tile order may differ.',
   'Startup control has accepted options but exits at no-analysis selection; data-specific setup is priced separately. Loader cache state can change.',
   'Three-slot finite FIFO; within-block decoder/analysis work is spread uniformly. Shared CPU/DRAM capacity enforced as additional whole-scan constraints.',
   'Native buffered output/process-exit boundary, without added fsync. Formatting scenario is not a universal maximum byte count.',
   'Generic resource profiles are observations under their recorded load, not guarantees of future throughput.'
  ],scope='Full major-stage runtime estimate for validation; no GWAS time, residual coefficient or N-indexed timing interpolation used. Not yet a validated crossover or guarantee.')
