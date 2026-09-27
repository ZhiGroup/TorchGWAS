"""Mechanistic CUDA stage estimate from tensor work and independent resources.

No N-indexed service timings, association residuals, interpolation or host
correction factors. Cache behavior is an explicit ideal LRU trace assumption.
This estimates a component, not process runtime or a crossover guarantee.
"""
from __future__ import annotations
from dataclasses import dataclass
import math,re
from .first_principles import positive,service

class IntervalLRU:
 """Exact fully-associative LRU for contiguous line-address range accesses.

 Resident intervals are oldest-to-newest. Access order is ascending addresses
 within each interval. By default stores overwrite full lines without read-for-ownership;
 setting write_allocate=True counts CPU write-miss ownership reads in the
 same reference stream as stores (not as a separate synthetic read pass);
 evicted dirty lines write back. Complexity depends on intervals, not bytes.
 """
 def __init__(self,capacity_lines):
  if not isinstance(capacity_lines,int) or capacity_lines<0:raise ValueError('Invalid cache capacity')
  self.capacity=capacity_lines;self.ranges=[]
 def access(self,lo,hi,write=False,write_allocate=False):
  if not 0<=lo<=hi:raise ValueError('Invalid line range')
  count=hi-lo
  if not count:return dict(read_fill=0,writeback=0,hits=0)
  if not self.capacity:return dict(read_fill=count if (not write or write_allocate) else 0,writeback=count if write else 0,hits=0)
  # Walk the access range at resident interval boundaries. Processing an
  # entire range at once would incorrectly preserve old tail hits after the
  # early misses have already evicted them (a classic cyclic-scan error).
  cursor=lo;hits=0;wb=0;fill=0
  while cursor<hi:
   resident=next(((a,b,d) for a,b,d in self.ranges if a<=cursor<b),None)
   if resident:
    stop=min(hi,resident[1]);dirty=resident[2] or write;hits+=stop-cursor
   else:
    next_start=min([a for a,b,d in self.ranges if a>cursor]+[hi])
    stop=min(hi,next_start);dirty=write
    if not write or write_allocate:fill+=stop-cursor
   updated=[]
   for a,b,d in self.ranges:
    if b<=cursor or a>=stop:updated.append((a,b,d))
    else:
     if a<cursor:updated.append((a,cursor,d))
     if b>stop:updated.append((stop,b,d))
   updated.append((cursor,stop,dirty))
   excess=sum(b-a for a,b,d in updated)-self.capacity
   while excess>0:
    a,b,d=updated[0];remove=min(excess,b-a);wb+=remove if d else 0
    if remove==b-a:updated.pop(0)
    else:updated[0]=(a+remove,b,d)
    excess-=remove
   # Coalesce adjacent intervals only when they have adjacent recency too.
   self.ranges=[]
   for a,b,d in updated:
    if self.ranges and self.ranges[-1][1]==a and self.ranges[-1][2]==d:
     old=self.ranges.pop();self.ranges.append((old[0],b,d))
    else:self.ranges.append((a,b,d))
   cursor=stop
  return dict(read_fill=fill,writeback=wb,hits=hits)

@dataclass(frozen=True)
class DeviceService:
 hbm_bytes_per_second: float
 l2_bytes_per_second: float
 fp32_flops_per_second: float
 kernel_launch_seconds: float
 host_dispatch_cpu_seconds: float
 available_l2_bytes: int
 sm_count: int
 gpu_fraction: float=1.0
 host_cpu_fraction: float=1.0
 cache_line_bytes: int=128
 fp64_flops_per_second: float | None=None
 fp64_tensor_flops_per_second: float | None=None
 def __post_init__(self):
  for k in ['hbm_bytes_per_second','l2_bytes_per_second','fp32_flops_per_second','kernel_launch_seconds','host_dispatch_cpu_seconds']:
   positive(k,getattr(self,k),True)
  for k in ['fp64_flops_per_second','fp64_tensor_flops_per_second']:
   if getattr(self,k) is not None:positive(k,getattr(self,k),True)
  for k in ['gpu_fraction','host_cpu_fraction']:
   if not 0<=getattr(self,k)<=1:raise ValueError(k+' must be in [0,1]')
  for k in ['sm_count','cache_line_bytes']:
   if not isinstance(getattr(self,k),int) or getattr(self,k)<1:raise ValueError(k+' must be positive integer')
  if not isinstance(self.available_l2_bytes,int) or self.available_l2_bytes<0:raise ValueError('Invalid L2 capacity')

def _gemv_work(samples,markers,width,main,reductions,*,dtype="float32"):
 """Logical service floor for the observed singleton-output cuBLAS path.

 Grid/block geometry identifies launches, not issued FLOPs or workspace
 padding inside this closed-source specialization. Keep those gaps explicit.
 """
 kernel=main[0];name=re.sub(r'\s+','',kernel['name'])
 scalar='double' if dtype=='float64' else 'float'
 element_bytes=8 if dtype=='float64' else 4
 signature='internal::gemvx::kernel<int,int,'+','.join([scalar]*4)+',false,true,false,false,'
 grid=kernel['geometry']['grid'];block=kernel['geometry'].get('block')
 # Duration-free captures expose four lane arrangements of this singleton
 # path. The x lanes partition output columns; z counts observed splits.
 # Accept only captured signature/block pairs, not arbitrary cuBLAS kernels.
 layouts=({(8,32,1):8} if dtype=="float64" else {(8,32,1):8,(32,16,1):9,(64,8,1):9,(128,4,1):9})
 specialization=layouts.get(tuple(block or ()))
 if (markers!=1 or specialization is None or signature+str(specialization)+',false,' not in name
     or 'cublasGemvParamsEx<int,' not in name
     or len(grid)!=3 or any(type(x) is not int or x<1 for x in grid)
     or grid[:2]!=[math.ceil(width/block[0]),1]):
  raise ValueError('Unsupported singleton GEMV shape, specialization or launch geometry')
 split=grid[2]
 if dtype=='float64' and (split!=1 or reductions):
  raise ValueError('FP64 singleton GEMV split reduction requires its own captured census')
 if bool(reductions)!=(split>1):raise ValueError('GEMV split grid requires an explicit reduction launch')
 if reductions:
  reduction=reductions[0]
  if ('cublasLt::splitKreduce_kernel<32,16,int,float,float,float,float,true,false,false>' not in re.sub(r'\s+','',reduction['name'])
      or reduction['geometry']['grid']!=[1,math.ceil(width/16),1]
      or reduction['geometry'].get('block')!=[32,16,1]):
   raise ValueError('Unsupported GEMV reduction launch geometry')
 workspace=4*width*split if reductions else 0
 return dict(useful_flops=2*samples*width,issued_flops=2*samples*width,
     arithmetic_accounting='logical_floor_not_verified_issued_work',
     grid_ctas=math.prod(grid),launched_ctas=math.prod(grid),noop_ctas=0,
     split_k=split,split_k_mode='separate_reduction' if reductions else 'none',
     swizzle=1,tile=None,input_l2_bytes=element_bytes*samples*(width+1),
     workspace_write_bytes=workspace,workspace_reduce_bytes=workspace,
     accumulation_logical_bytes=0,reduce_adds=width*(split-1),kernel_count=1+bool(reductions),
     unpriced_terms=['GEMV internal arithmetic padding, intra-CTA reduction and workspace layout',
       'GEMV physical operand rereads and reduction synchronization',
       'GEMV host dispatch equivalence to the independent tiny GEMM primitive'],
     inner_k_tile_verified=False,
     policy='Observed singleton GEMV specialization only. Logical dot-product, operand and split-workspace floors; actual launch count and grid. Issued FLOPs and physical traffic are not recovered from kernel names. No runtime fields accepted.')


def _dot_gemv_work(samples,markers,width,main,reductions):
 """Logical floors for the captured batched-dot singleton GEMV pair."""
 kernel=main[0];name=re.sub(r'\s+','',kernel['name']);geometry=kernel['geometry']
 grid=geometry['grid']
 if (markers!=1 or 'dot_kernel<float,128,0,cublasDotParams<' not in name
     or 'cublasGemvTensorStridedBatched<floatconst>' not in name
     or geometry.get('block')!=[128,1,1] or len(grid)!=3
     or any(type(v) is not int or v<1 for v in grid) or grid[1:]!=[1,width]):
  raise ValueError('Unsupported singleton dot GEMV launch geometry')
 if len(reductions)!=1:raise ValueError('Dot GEMV requires one explicit reduction launch')
 reduction=reductions[0]
 if ('reduce_1Block_kernel<float,128,7,cublasGemvTensorStridedBatched<float>,' not in re.sub(r'\s+','',reduction['name'])
     or reduction['geometry']['grid']!=[1,1,width] or reduction['geometry'].get('block')!=[128,1,1]):
  raise ValueError('Unsupported dot GEMV reduction launch geometry')
 partials=grid[0];workspace=4*width*partials
 return dict(useful_flops=2*samples*width,issued_flops=2*samples*width,
     arithmetic_accounting='logical_floor_not_verified_issued_work',
     grid_ctas=math.prod(grid),launched_ctas=math.prod(grid),noop_ctas=0,
     split_k=partials,split_k_mode='separate_reduction',swizzle=1,tile=None,
     input_l2_bytes=4*samples*(width+1),workspace_write_bytes=workspace,
     workspace_reduce_bytes=workspace,accumulation_logical_bytes=0,
     reduce_adds=width*(partials-1),kernel_count=2,inner_k_tile_verified=False,
     unpriced_terms=['Dot GEMV internal padding, per-CTA accumulation and workspace layout',
       'Dot GEMV physical operand rereads and reduction synchronization',
       'Dot GEMV host dispatch equivalence to the independent tiny GEMM primitive'],
     policy='Captured batched-dot singleton path only. Logical dot-product and partial workspace floors; '
            'actual grid and launch count. Closed-source issued work and physical traffic remain unresolved.')


def _nsp_gemv_work(samples,markers,width,main,reductions):
 """Logical floors for the two captured A100 FP32 singleton NSP layouts."""
 kernel=main[0];name=re.sub(r'\s+','',kernel['name']);geometry=kernel['geometry']
 # These are duration-free observed specialization/layout pairs, not a timing
 # grid. Do not infer support for unseen sample counts, splits or lane layouts.
 layout={(2049,1,515):(16,[32,24,1]),(2049,1,4):(32,[8,32,1])}.get((samples,markers,width))
 if layout is None:raise ValueError('Unsupported singleton NSP GEMV launch census')
 lanes,block=layout
 signature=f'gemvNSP_kernel<float,float,float,float,1,{lanes},4,1024,false,cublasGemvParamsEx<int,'
 if (signature not in name or geometry.get('block')!=block
     or geometry.get('grid')!=[math.ceil(width/block[0]),1,8] or len(reductions)!=1):
  raise ValueError('Unsupported singleton NSP GEMV launch census')
 reduction=reductions[0]
 if ('cublasLt::splitKreduce_kernel<32,16,int,float,float,float,float,true,false,false>' not in re.sub(r'\s+','',reduction['name'])
     or reduction['geometry'].get('grid')!=[1,math.ceil(width/16),1]
     or reduction['geometry'].get('block')!=[32,16,1]):
  raise ValueError('Unsupported singleton NSP GEMV reduction census')
 split=geometry['grid'][2];workspace=4*width*split
 return dict(useful_flops=2*samples*width,issued_flops=2*samples*width,
     arithmetic_accounting='logical_floor_not_verified_issued_work',
     grid_ctas=math.prod(geometry['grid']),launched_ctas=math.prod(geometry['grid']),noop_ctas=0,
     split_k=split,split_k_mode='separate_reduction',swizzle=1,tile=None,
     input_l2_bytes=4*samples*(width+1),workspace_write_bytes=workspace,
     workspace_reduce_bytes=workspace,accumulation_logical_bytes=0,
     reduce_adds=width*(split-1),kernel_count=2,inner_k_tile_verified=False,
     unpriced_terms=['NSP GEMV internal padding, partial workspace layout and per-CTA reduction',
       'NSP GEMV physical operand rereads and synchronization',
       'NSP GEMV host dispatch equivalence to the independent tiny GEMM primitive'],
     policy='Captured N2049/B1/W4 and W515 layouts only. Logical arithmetic and partial workspace floors; no runtime fields.')

def gemm_work(samples,markers,width,kernels,*,dtype="float32"):
 """Issued tile work from duration-free compiled launch geometry.

 CUTLASS identity swizzle: grid=(logical_m*s, ceil(logical_n/s), split).
 See NVIDIA CUTLASS v3.5.1 threadblock/threadblock_swizzle.h. Swizzle padding
 launches no-op CTAs; it must not be counted as full matrix multiplication.
 A z dimension >1 does not imply a separate reduction launch. Its existence
 comes from the census, not from the split count.
 """
 for value in (samples,markers,width):
  if type(value) is not int or value<1:raise ValueError('Matrix dimensions must be positive integers')
 if dtype not in ('float32','float64'):raise ValueError('Unsupported GEMM precision')
 element_bytes=4 if dtype=='float32' else 8
 if dtype=='float64':
  # CUTLASS v3.5.1 generator names encode the FP64 accumulator (d),
  # 8x8x4 instruction, threadblock shape and stages. Admit only the captured
  # Ampere family; other cuBLAS FP64 algorithms need their own census.
  main=[k for k in kernels if 'cutlass_80_tensorop_d884gemm_' in k['name'] or 'internal::gemvx::kernel<int,int,double,double,double,double,' in re.sub(r'\s+','',k['name'])]
 else:
  main=[k for k in kernels if ('sgemm_' in k['name'] or 'execute_kernel' in k['name'] or 'internal::gemvx::kernel<' in k['name'] or 'dot_kernel<' in k['name'] or 'gemvNSP_kernel<' in k['name']) and 'split' not in k['name'].lower()]
 reductions=[k for k in kernels if ('split_k' in k['name'].lower() or 'splitkreduce' in k['name'].lower() or 'reduce_1Block_kernel<' in k['name'])]
 if len(main)!=1 or len(reductions)>1:raise ValueError('Unsupported GEMM launch census; supply its counted work explicitly')
 if 'gemvNSP_kernel<' in main[0]['name']:return _nsp_gemv_work(samples,markers,width,main,reductions)
 if 'dot_kernel<' in main[0]['name']:return _dot_gemv_work(samples,markers,width,main,reductions)
 if 'internal::gemvx::kernel<' in main[0]['name']:
  # A singleton product is the same GEMV whichever operand is the vector:
  # the dense JAGWAS projection has one output row, a triangular block one
  # output column. The GEMV's output length is the non-singleton dimension.
  if width==1 and markers!=1:markers,width=1,markers
  result=_gemv_work(samples,markers,width,main,reductions,dtype=dtype)
  if dtype=='float64':result.update(dtype='float64',arithmetic_kind='fp64_scalar')
  return result
 kernel=main[0];name=kernel['name'];grid=kernel['geometry']['grid']
 if len(grid)!=3 or any(not isinstance(x,int) or x<1 for x in grid):raise ValueError('Invalid GEMM grid')
 swizzle=1
 if dtype=='float64' and reductions:raise ValueError('Separate FP64 split-K reduction requires a verified census')
 m=re.search(r'tilesize(\d+)x(\d+)x(\d+)',name)
 if m:
  bm,bn,bk=map(int,m.groups());gm,gn,split=grid
 elif 'cutlass' in name:
  m=re.search(r'(?:sgemm|d884gemm)_(\d+)x(\d+)_(\d+)x',name)
  if not m:raise ValueError('Unknown CUTLASS tile spelling')
  bn,bm,bk=map(int,m.groups());gm,gn=math.ceil(markers/bm),math.ceil(width/bn);split=grid[2]
  matches=[s for s in (1,2,4,8) if grid[0]==gm*s and grid[1]==math.ceil(gn/s)]
  if len(matches)!=1:raise ValueError('CUTLASS grid does not match identity-swizzled problem dimensions')
  swizzle=matches[0]
 else:
  m=re.search(r'sgemm_(\d+)x(\d+)',name)
  if not m:raise ValueError('Unknown GEMM tile spelling')
  bn,bm=map(int,m.groups());gn,gm,split=grid;bk=1
 if gm!=math.ceil(markers/bm) or gn!=math.ceil(width/bn):raise ValueError('GEMM grid does not match problem dimensions')
 if reductions and split==1:raise ValueError('Separate split-K reduction without split-K work')
 separate=bool(reductions)
 mode='separate_reduction' if separate else ('in_kernel_reduction_unresolved' if split>1 else 'none')
 issued_k=math.ceil(samples/(split*bk))*split*bk
 padded_m,padded_n=gm*bm,gn*bn
 # Without a separate reduction, no separate reduction workspace is asserted.
 # Partial sums must still be combined; count logical output read/write work,
 # but expose unknown atomics/semaphore/synchronization rather than price it zero.
 accumulation=2*element_bytes*markers*width*(split-1) if split>1 and not separate else 0
 result=dict(useful_flops=2*samples*markers*width,issued_flops=2*issued_k*padded_m*padded_n,
     arithmetic_accounting='compiled_tile_work',
     grid_ctas=gm*gn*split,launched_ctas=math.prod(grid),noop_ctas=math.prod(grid)-gm*gn*split,
     split_k=split,split_k_mode=mode,swizzle=swizzle,tile=[bm,bn,bk],
     input_l2_bytes=element_bytes*samples*(markers*gn+width*gm),
     workspace_write_bytes=element_bytes*padded_m*padded_n*split if separate else 0,
     workspace_reduce_bytes=element_bytes*padded_m*padded_n*split if separate else 0,
     accumulation_logical_bytes=accumulation,
     reduce_adds=markers*width*(split-1),kernel_count=1+int(separate),
     unpriced_terms=['in-kernel split-K synchronization/atomic service and physical accumulation traffic'] if mode=='in_kernel_reduction_unresolved' else [],
     inner_k_tile_verified=bool(re.search(r'tilesize|cutlass',name)),
     policy='Logical tiles reconstructed from compiled grid and identity swizzle; reduction launches counted explicitly. In-kernel reduction protocol remains unresolved. No runtime fields accepted.')
 if dtype=='float64':result.update(dtype='float64',arithmetic_kind='fp64_tensor')
 return result

GEMM_OPS=('aten.mm.default','aten.mm.out')


def _is_reduction_kernel(name):
 return 'split_k' in name.lower() or 'splitkreduce' in name.lower() or 'reduce_1Block_kernel<' in name


def _combined_gemm(products):
 """Totals over several products; `gemms` keeps each one's own census."""
 kinds={p.get('arithmetic_kind') for p in products}
 combined={key:sum(p[key] for p in products) for key in
     ('useful_flops','issued_flops','reduce_adds','kernel_count','input_l2_bytes','workspace_write_bytes',
      'workspace_reduce_bytes','accumulation_logical_bytes','grid_ctas') if all(key in p for p in products)}
 combined.update(gemms=products,products=len(products),arithmetic_kind=kinds.pop() if len(kinds)==1 else 'mixed',
     split_k=[p.get('split_k') for p in products],
     unpriced_terms=sorted({term for p in products for term in p.get('unpriced_terms',[])}))
 return combined


def tensor_stage_service(work,resources:DeviceService,kernels,initial_cache=None,host_primitives=None,*,gemm_dimensions=None,gemm_dtype="float32",host_primitive_resolver=None):
 """Evaluate ordered kernel service, with separate host dispatch capacity.

 Logical operands are scanned in order; actual GPU tiles interleave operands.
 Hence ideal fully-associative operand-major LRU is a declared approximation,
 not a claim to recover physical hardware traffic exactly. A zero-cache run
 supplies a useful no-reuse sensitivity, not a guaranteed runtime upper bound.
 """
 if gemm_dimensions is None:
  n,b,k,c=(work[t] for t in ['samples','markers','traits','covariates'])
  gemm_dimensions=(n,b,k+c+1)
 # One (inner, rows, columns) triple, or a list of them: one per matrix
 # product in execution order (the block-triangular JAGWAS projection).
 many=isinstance(gemm_dimensions,list) and all(isinstance(d,(tuple,list)) for d in gemm_dimensions)
 dimension_list=[tuple(d) for d in gemm_dimensions] if many else [gemm_dimensions]
 if not dimension_list or any(not isinstance(d,(tuple,list)) or len(d)!=3 for d in dimension_list):
  raise ValueError('Explicit inner/rows/columns GEMM dimensions required')
 active=[s for s in work['steps'] if not(s['alias_only'] or s['allocation_only'])]
 gemm_steps=[i for i,s in enumerate(active) if s['op'] in GEMM_OPS]
 if many and len(gemm_steps)!=len(dimension_list):
  raise ValueError(f'{len(gemm_steps)} traced matrix products but {len(dimension_list)} GEMM dimensions')
 # Each product's own launches: its main kernel, plus a separate split-K
 # reduction kernel when the census has one right after it.
 gemms={};position=0
 for i,step in enumerate(active):
  if i in gemm_steps:
   count=2 if position+1<len(kernels) and _is_reduction_kernel(kernels[position+1]['name']) else 1
   dims=dimension_list[gemm_steps.index(i)] if many else dimension_list[0]
   gemms[i]=gemm_work(*dims,kernels[position:position+count] if many else kernels,dtype=gemm_dtype)
   position+=gemms[i]['kernel_count'] if many else count
  else:position+=1
 if not gemm_steps:gemms[None]=gemm_work(*dimension_list[0],kernels,dtype=gemm_dtype)
 gemm=next(iter(gemms.values())) if len(gemms)==1 else _combined_gemm(list(gemms.values()))
 double_work=any(t['dtype']=='torch.float64' and t.get('device')=='meta'
                 for step in active for t in step['inputs']+step['outputs'])
 tensor_fp64=any(g.get('arithmetic_kind')=='fp64_tensor' for g in gemms.values())
 double_rates=(['fp64_flops_per_second'] if double_work else [])+(['fp64_tensor_flops_per_second'] if tensor_fp64 else [])
 if any(getattr(resources,name) is None for name in double_rates):
  raise ValueError('Missing independent FP64 arithmetic resources')
 expected=len(active)-len(gemm_steps)+sum(gemms[i]['kernel_count'] for i in gemm_steps)
 if expected!=len(kernels):raise ValueError(f'Unaccounted compiled kernels: graph expects {expected}, census contains {len(kernels)}')
 blocked=[name for name in ['gpu_fraction','hbm_bytes_per_second','l2_bytes_per_second','fp32_flops_per_second']+double_rates if getattr(resources,name)==0]
 if resources.host_cpu_fraction==0 and (resources.host_dispatch_cpu_seconds>0 or (host_primitives is not None and any(host_primitives.values()))):blocked.append('host_cpu_fraction')
 if blocked:return dict(estimated_span_seconds=None,status='zero_available_capacity',blocked_resources=blocked,
     kernel_count=expected,gemm=gemm,source_sha256=work['source_sha256'],scope='No finite service can be predicted for positive work with zero available capacity.')
 line=resources.cache_line_bytes;cache=initial_cache or IntervalLRU(resources.available_l2_bytes//line)
 layout={};address=0
 for step in work['steps']:
  for t in step['inputs']+step['outputs']:
   if t['storage'] not in layout:
    layout[t['storage']]=address;address+=math.ceil(t['storage_bytes']/line)+1
 def extent(t):
  # The compact K=1 views have strides below a cache line. For wide trait
  # views this bounding span is conservative traffic and is identified.
  item=t['bytes']//max(1,math.prod(t['shape']))
  span=item*(1+sum((d-1)*s for d,s in zip(t['shape'],t['stride']))) if t['bytes'] else 0
  start=layout[t['storage']]+t['offset_bytes']//line
  return start,layout[t['storage']]+math.ceil((t['offset_bytes']+span)/line)
 rows=[];index=0;gpu_end=0.;host_end=0.
 host_cpu=0.;call_ready={};host_rows=[]
 if host_primitives is not None:
  if 'host_calls' not in work:raise ValueError('Typed dispatch requires source host-call census')
  for call in work['host_calls']:
   primitive=(host_primitive_resolver or host_primitive_name)(call)
   if primitive not in host_primitives:raise ValueError('Missing independent host primitive: '+primitive)
   cpu=positive(primitive,host_primitives[primitive],True)
   host_cpu+=cpu;ready=service(host_cpu,resources.host_cpu_fraction)
   call_ready[call['id']]=ready
   host_rows.append(dict(call_id=call['id'],primitive=primitive,cpu_seconds=cpu,submit_finish=ready))
 for position,step in enumerate(active):
  is_gemm=step['op'] in GEMM_OPS;product=gemms.get(position) if is_gemm else None
  count=product['kernel_count'] if is_gemm else 1
  geometry=kernels[index:index+count];index+=count
  hbm=0
  read={} if step['shape_only_inputs'] else {t['storage']:t for t in step['inputs'] if t['device']=='meta'}
  write={t['storage']:t for t in step['outputs'] if t['device']=='meta'}
  for t in read.values():
   v=cache.access(*extent(t));hbm+=line*(v['read_fill']+v['writeback'])
  for t in write.values():
   v=cache.access(*extent(t),write=True);hbm+=line*(v['read_fill']+v['writeback'])
  l2=step['logical_bytes'];flops=0
  if is_gemm:
   l2=product['input_l2_bytes']+step['write_bytes']+product['workspace_write_bytes']+product['workspace_reduce_bytes']+product['accumulation_logical_bytes']
   # Split-K workspace is an additional explicitly accounted stream. Its
   # assumed non-reuse is conservative traffic, not a fitted GEMM penalty.
   hbm+=product['workspace_write_bytes']+product['workspace_reduce_bytes']+product['accumulation_logical_bytes']
   flops=product['issued_flops']+product['reduce_adds']
  elif step['op'].startswith(('aten.sum','aten.amin','aten.amax')):
   flops=max(0,math.prod(step['inputs'][0]['shape'])-math.prod(step['outputs'][0]['shape']))
  else:flops=max([math.prod(t['shape']) for t in step['outputs']]+[0])
  blocks=product['grid_ctas'] if is_gemm else math.prod(geometry[0]['geometry']['grid'])
  active_sm_fraction=min(1.,blocks/resources.sm_count)
  double_step=any(t['dtype']=='torch.float64' and t.get('device')=='meta' for t in step['inputs']+step['outputs'])
  arithmetic_rate=(resources.fp64_tensor_flops_per_second if is_gemm and product.get('arithmetic_kind')=='fp64_tensor' else
                   resources.fp64_flops_per_second if double_step else resources.fp32_flops_per_second)
  math_seconds=service(flops,arithmetic_rate*active_sm_fraction*resources.gpu_fraction)
  memory_seconds=max(service(hbm,resources.hbm_bytes_per_second*resources.gpu_fraction),
                     service(l2,resources.l2_bytes_per_second*resources.gpu_fraction))
  kernel_seconds=count*service(resources.kernel_launch_seconds,resources.gpu_fraction)+max(math_seconds,memory_seconds)
  submit=service(resources.host_dispatch_cpu_seconds,resources.host_cpu_fraction)
  # One tensor dispatch may enqueue two cuBLAS kernels. CPU serialization and
  # GPU serialization overlap in the same order as the eager implementation.
  if host_primitives is None:host_end+=submit
  else:host_end=call_ready[step['host_call_id']]
  gpu_start=max(gpu_end,host_end);gpu_end=gpu_start+kernel_seconds
  rows.append(dict(op=step['op'],phase=step['phase'],kernel_count=count,logical_bytes=step['logical_bytes'],
      modeled_l2_bytes=l2,modeled_hbm_bytes=hbm,issued_flops=flops,
      memory_seconds=memory_seconds,math_seconds=math_seconds,kernel_service_seconds=kernel_seconds,
      host_submit_finish=host_end,gpu_finish=gpu_end,grid_ctas=blocks))
 return dict(estimated_span_seconds=gpu_end,kernel_service_seconds=sum(r['kernel_service_seconds'] for r in rows),
    host_dispatch_cpu_seconds=host_cpu if host_primitives is not None else len(active)*resources.host_dispatch_cpu_seconds,
    unpriced_terms=gemm['unpriced_terms'],host_calls=host_rows,host_dispatch_mode='source_grouped_typed' if host_primitives is not None else 'generic_per_active_operation',
    unpriced_host_property_accesses=work.get('property_access_count',0),
    conversion_service_seconds=sum(r['kernel_service_seconds'] for r in rows if r['phase']=='conversion'),
    statistics_service_seconds=sum(r['kernel_service_seconds'] for r in rows if r['phase']=='statistics'),
    modeled_hbm_bytes=sum(r['modeled_hbm_bytes'] for r in rows),logical_bytes=work['logical_bytes'],
    source_sha256=work['source_sha256'],kernel_count=expected,gemm=gemm,operations=rows,
    scope='Analytical CUDA-component estimate under explicit cache, bandwidth, issue and load assumptions; not a process-runtime bound or crossover. No association timing used.',
    assumptions=['Fully associative operand-major LRU; ideal full-line write stores without read-for-ownership.',
      'Resource ceilings characterized by independent generic operations; throughput attainment is a model approximation.',
      'Single CUDA stream for conversion/statistics; host dispatch is a separate serial resource.',
      'Cache-line bounding spans for strided views; no assumed cache hit fraction from observed GWAS timings.',
      'GEMM tile work uses compiled launch geometry; singleton GEMV uses explicitly labelled logical floors with unresolved internal work. No duration data.',
      ('FP64 tensor operations use separate scalar and Tensor Core issue capacities; integer work uses FP32-equivalent capacity. Operation-specific instruction latency is unresolved.' if double_work else
       'Scalar/integer/reduction math uses FP32-equivalent issue capacity; operation-specific instruction latency is not resolved.')])


def host_primitive_name(call):
 """Match API semantics, not N or elapsed time, to a fixed tiny primitive."""
 name=call['name'];dtypes=call['input_dtypes'];shapes=call['input_shapes'];kw=call['kwargs']
 if name=='__getitem__':return 'view_index'
 if name=='squeeze':return 'view_squeeze'
 if name=='__eq__' and dtypes==['torch.int8']:return 'int8_compare'
 if name=='to':return {'torch.int8':'int8_to_fp32','torch.int64':'long_to_float'}[dtypes[0]]
 if name=='where':
  if len(dtypes)==2:return 'where_scalar'
  if dtypes[1]=='torch.uint8':return 'where_u8'
  return 'where_tensor'
 if name in ('zeros_like','ones_like','full_like'):
  return name+('_u8' if kw.get('dtype')=='torch.uint8' else '')
 if name=='sum':return 'sum_bool' if dtypes[0]=='torch.bool' else 'sum_float'
 if name=='clamp':return 'clamp_int64' if dtypes[0]=='torch.int64' else 'clamp_float'
 if name=='div':return 'div_mixed' if 'torch.int64' in dtypes else 'div_float'
 if name=='sub':return 'sub_scalar' if len(dtypes)==1 else ('sub_broadcast' if shapes[0]!=shapes[1] else 'sub_float')
 if name=='gt':return 'gt_scalar' if len(dtypes)==1 else 'gt_tensor'
 if name=='lt':return 'lt_scalar' if len(dtypes)==1 else 'lt_tensor'
 direct={'isnan':'isnan','__invert__':'invert_bool','amin':'amin','amax':'amax',
  'matmul':'gemm_fp32','mul':'mul_float','__and__':'and_bool','__or__':'or_bool','isfinite':'isfinite','sqrt':'sqrt_float'}
 if name not in direct:raise ValueError('Unpriced host API '+name)
 return direct[name]


