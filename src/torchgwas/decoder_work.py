"""Reader-specific source units, independent of GWAS runtimes or N grids.

Logical bytes are not automatically physical memory traffic. Units identify the
actual source branch so an SSE2 vector cannot be priced as four scalar words.
"""
from __future__ import annotations
from .first_principles import positive

IMPLEMENTATIONS=('torch_native_int8','pgenlib_sse2','pgenlib_avx2_bmi2')

def census_chunk_ranges(census, chunk_markers):
 """Exact file-global ranges; chunk_markers remains the allocated capacity.

 Explicit layouts are contiguous, nonempty and bounded by that capacity.
 They describe a complete scan extent, never an observation of elapsed time.
 """
 m=census['markers'];span=census.get('variant_range',[0,m])
 if (type(chunk_markers) is not int or chunk_markers<1 or type(m) is not int or m<1
     or not isinstance(span,(list,tuple)) or len(span)!=2
     or any(type(v) is not int or v<0 for v in span) or span[1]-span[0]!=m
     or census.get('chunk_markers',chunk_markers)!=chunk_markers):
  raise ValueError('Chunk census parent range or geometry mismatch')
 if 'chunk_ranges' not in census:
  return [[lo,min(lo+chunk_markers,span[1])] for lo in range(span[0],span[1],chunk_markers)]
 ranges=census['chunk_ranges']
 if not isinstance(ranges,list) or not ranges:raise ValueError('Nonempty explicit chunk ranges required')
 cursor=span[0]
 for row in ranges:
  if (not isinstance(row,(list,tuple)) or len(row)!=2 or any(type(v) is not int for v in row)
      or row[0]!=cursor or not 0<row[1]-row[0]<=chunk_markers or row[1]>span[1]):
   raise ValueError('Explicit chunk ranges must cover the scan once within capacity')
  cursor=row[1]
 if cursor!=span[1]:raise ValueError('Explicit chunk ranges do not cover the complete scan')
 return [list(row) for row in ranges]


def scan_chunk_count(data, profile):
 """Cheap expansion budget before constructing any per-chunk work graph."""
 if 'chunk_ranges' in data['encoded']:
  ranges=data['encoded']['chunk_ranges']
  if not isinstance(ranges,list) or not ranges:raise ValueError('Nonempty explicit chunk ranges required')
  return len(ranges)
 return (data['markers']+profile['chunk_markers']-1)//profile['chunk_markers']


def require_regular_memory(candidate):
 """Fixed-run peaks do not bound mixed-shape adaptive allocation lifetimes."""
 if any('chunk_ranges' in tile['data']['encoded'] for tile in candidate['tiles']):
  raise ValueError('Explicit chunk schedules require adaptive_candidate_memory on the fixed-capacity candidate')


def native_ld_replays(census):
 """Validate base-only decode work separately from the contiguous read prefix."""
 count=census.get('ld_records_at_chunk_starts',0)
 if type(count) is not int or count<0:raise ValueError('Invalid native LD replay count')
 rows=census.get('native_ld_replays')
 if rows is None:return None if count else []
 if not isinstance(rows,list) or len(rows)!=count:raise ValueError('Native LD replay count mismatch')
 if not rows:return []
 span=census['variant_range'];chunk=census['chunk_markers'];previous=-1;payload=0
 if (type(chunk) is not int or chunk<1 or len(span)!=2
     or any(type(v) is not int or v<0 for v in span) or span[1]<=span[0]):
  raise ValueError('Invalid native LD replay geometry')
 starts={row[0] for row in census_chunk_ranges(census,chunk)} if 'chunk_ranges' in census else None
 for row in rows:
  if not isinstance(row,dict) or set(row)!={'chunk_start','base_variant','read_prefix_bytes','skipped_prefix_ld_records','base_record'}:
   raise ValueError('Unknown native LD replay fields')
  at,base,amount,skipped=[row[k] for k in ('chunk_start','base_variant','read_prefix_bytes','skipped_prefix_ld_records')]
  if any(type(v) is not int or v<0 for v in (at,base,amount,skipped)):
   raise ValueError('Invalid native LD replay extent')
  aligned=at in starts if starts is not None else (at-span[0])%chunk==0
  if not span[0]<=at<span[1] or not aligned or at<=previous or base>=at or skipped!=at-base-1:
   raise ValueError('Native LD replay range mismatch')
  record=row['base_record'];forms={int(k):v for k,v in record['record_form_counts'].items()}
  if (record['samples']!=census['samples'] or record['markers']!=1 or record['variant_range']!=[base,base+1]
      or any(type(v) is not int or v<0 for v in forms.values())
      or sum(forms.values())!=1 or set(forms)-{0,1,4,6,7} or record.get('ld_records_at_chunk_starts',0)
      or record.get('native_ld_replays')):
   raise ValueError('Native LD replay requires one non-LD base record')
  for key in ('path','file_markers','file_bytes','index_bytes'):
   if record.get(key)!=census.get(key):raise ValueError('Native LD replay file mismatch: '+key)
  size=record['record_payload_bytes']
  if type(size) is not int or size<1 or amount<size:raise ValueError('Native LD prefix must include its base payload')
  payload+=size;previous=at
 if payload!=census['additional_base_record_bytes_if_every_chunk_restarts']:
  raise ValueError('Native LD base payload does not conserve census')
 return rows


def native_read_layout(census):
 """Reader payload extent plus additional LD scratch/relative-offset storage."""
 if census.get('kind')=='torchgwas.pgen_memory_layout.v1':
  from .pgen_memory_layout import reader_layout
  return reader_layout(census)
 rows=native_ld_replays(census)
 if rows is None:raise ValueError('Native LD restart census is incomplete')
 return dict(read_bytes=census['record_payload_bytes']+sum(r['read_prefix_bytes'] for r in rows),
             decode_input_bytes=census['record_payload_bytes']+sum(r['base_record']['record_payload_bytes'] for r in rows),
             extra_workspace_bytes=max((((census['samples']+3)//4)+8*(r['chunk_start']-r['base_variant']) for r in rows),default=0))


def native_reader_workspace(chunks):
 """Conservative per-reader maxima, including grow-only input and LD scratch."""
 read=extra=0
 for chunk in chunks:
  layout=native_read_layout(chunk)
  read=max(read,layout['read_bytes'])
  extra=max(extra,layout['extra_workspace_bytes'])
 return read,extra


def _native_header_units(n, forms, *, expand=True):
 """Native loops known from sample count and record forms, without payloads."""
 m=sum(forms.values());packed=(n+3)//4;one=forms.get(1,0)
 base_bytes=packed*(m-forms.get(2,0)-forms.get(3,0))
 units=dict(copy_packed_byte=packed*sum(forms.get(f,0) for f in (0,2,3))+base_bytes,
  fill_packed_byte=packed*sum(forms.get(f,0) for f in (4,6,7))+one*(packed-2*(n//8)),
  invert_packed_byte=packed*forms.get(3,0),onebit8_native=one*(n//8),
  native_onebit_tail_bit_test=one*(n%8),difflist_record_header=m-forms.get(0,0))
 if expand:units.update(expand4_int8=m*(n//4),expand_int8_tail_sample=m*(n%4))
 return {key:value for key,value in units.items() if value},base_bytes


def _native_decoder_memory_bytes(n, markers, input_bytes, decode_input_bytes,
                                 base_update_bytes, replay_packed_bytes, *, separate_buffered_read=False):
 """Existing logical traffic approximation, shared by exact and interval work.

 Packed staging write/read, int8 CPU write allocation, input buffering and LD
 base refresh. A separate buffered-read node already charges the blob write.
 These are model bytes, not a guarantee of physical memory transactions.
 """
 return ((0 if separate_buffered_read else input_bytes)+decode_input_bytes
         +2*((n+3)//4)*markers+2*n*markers+2*base_update_bytes+replay_packed_bytes)


def decoder_work(census, implementation, *, restart_ld_bases=False):
 if census.get('kind')=='torchgwas.pgen_memory_layout.v1':
  raise ValueError('Memory-only PGEN layout cannot price decoder work; collect exact counts productively or load valid cached counts')
 if census.get('kind') in ('torchgwas.pgen_header_work_bounds.v1','torchgwas.pgen_header_window.v1',
                          'torchgwas.pgen_header_schedule_bounds.v1'):
  raise ValueError('Header work bounds are intervals, not an exact decoder census')
 if implementation not in IMPLEMENTATIONS:raise ValueError('Unknown decoder implementation')
 n,m=census['samples'],census['markers']
 if not isinstance(n,int) or n<1 or not isinstance(m,int) or m<1:raise ValueError('Invalid dimensions')
 forms={int(k):int(v) for k,v in census['record_form_counts'].items()}
 if set(forms)-{0,1,2,3,4,6,7} or any(v<0 for v in forms.values()) or sum(forms.values())!=m:
  raise ValueError('Unsupported or inconsistent record forms')
 units={};unknown=[];packed=(n+3)//4;one=forms.get(1,0)
 def add(name,amount):
  positive(name,amount,True)
  if amount:units[name]=amount
 add('copy_packed_byte',packed*(forms.get(0,0)+forms.get(2,0)+forms.get(3,0)))
 add('fill_packed_byte',packed*sum(forms.get(f,0) for f in [4,6,7]))
 add('invert_packed_byte',packed*forms.get(3,0))
 if implementation=='torch_native_int8':
  # decode_range refreshes its packed LD base after every non-LD record,
  # including files with no LD references. This is additional to form0 input
  # copies and form2/3 copies from the existing base.
  units,base_update_bytes=_native_header_units(n,forms)
  if n%8 and one:
   if 'native_onebit_tail_high_count' not in census:unknown.append('native one-bit tail high-bit count')
   else:add('set_category',census['native_onebit_tail_high_count'])
  onebit_span=one*n
 elif implementation=='pgenlib_sse2':
  vectors=n//128;words=(n+31)//32-4*vectors
  add('onebit128_sse2',one*vectors)
  add('onebit32_pgenlib_sse2',one*words)
  onebit_span=one*(128*vectors+32*words)
 else:
  add('onebit32_pgenlib_avx2_bmi2',one*((n+31)//32))
  onebit_span=one*32*((n+31)//32)
 if implementation!='torch_native_int8' and n%32:
  add('pgenlib_partial_word_load',one)
 diff_records=m-forms.get(0,0)
 if diff_records:
  hist=census.get('total_varint_lengths')
  if hist is None:unknown.append('exact variable-integer length histogram')
  else:
   for length,count in hist.items():
    if not 1<=int(length)<=5:raise ValueError('Invalid 32-bit variable integer length')
    add('uleb'+str(length),count)
   if sum(int(k)*v for k,v in hist.items())!=census['source_work']['varint_bytes']:
    raise ValueError('Variable-integer byte counts do not conserve work')
  add('set_category',units.get('set_category',0)+census['source_work']['difflist_entries'])
  add('difflist_group_absolute_id',census['total_difflist_groups'])
  add('difflist_category_extract',census['source_work']['difflist_entries'])
  add('difflist_record_header',diff_records)
 replay_rows=[]
 if restart_ld_bases and census['ld_records_at_chunk_starts']:
  replay_rows=native_ld_replays(census) if implementation=='torch_native_int8' else None
  if replay_rows is None:
   unknown.append('full source census of replayed LD-base records at reader restarts')
  else:
   for row in replay_rows:
    replay=decoder_work(row['base_record'],implementation)
    unknown.extend(replay['uncounted_mechanisms'])
    for name,amount in replay['source_units'].items():
     # The base is decoded into packed scratch, never expanded to int8.
     if name not in ('expand4_int8','expand_int8_tail_sample'):add(name,units.get(name,0)+amount)
    base_update_bytes+=replay['native_ld_base_update_bytes']
 result=dict(implementation=implementation,samples=n,markers=m,source_units=units,
  onebit_covered_or_padded_samples=onebit_span,
  logical_final_packed_bytes=packed*m,
  native_ld_base_update_bytes=base_update_bytes if implementation=='torch_native_int8' else 0,
  logical_final_int8_bytes=n*m if implementation=='torch_native_int8' else 0,
  uncounted_mechanisms=unknown,
  source='native/pgen_decode.c and bundled pgenlib ParseOnebitUnsafe/ParseAndApplyDifflist',
  scope='Instruction-loop counts and logical storage only; record dispatch, cache-line traffic, synchronization and initialization still require service accounting. No runtime or crossover is asserted.')
 if replay_rows:
  result.update(native_ld_replay_records=len(replay_rows),native_ld_replay_packed_bytes=packed*len(replay_rows),
                native_ld_replay_read_bytes=sum(r['read_prefix_bytes'] for r in replay_rows),
                native_ld_replay_decoded_input_bytes=sum(r['base_record']['record_payload_bytes'] for r in replay_rows))
 return result


def decoder_chunk_work(census, chunk_markers, implementation='torch_native_int8'):
 """Validate exact source counts for the caller's contiguous chunk geometry.

 Legacy aggregate censuses return None and require an explicit uniform-work
 scenario. Supplied chunk counts must conserve payload, record forms and every
 priced source unit; malformed or mismatched geometry is never reweighted.
 """
 if 'chunks' not in census:
  if 'chunk_ranges' in census:raise ValueError('Explicit chunk ranges require exact per-chunk census')
  return None
 from collections import Counter
 if isinstance(chunk_markers,bool) or not isinstance(chunk_markers,int) or chunk_markers<1:
  raise ValueError('Positive integer chunk size required')
 n,m=census['samples'],census['markers']
 span=census.get('variant_range',[0,m])
 if (len(span)!=2 or any(isinstance(v,bool) or not isinstance(v,int) or v<0 for v in span)
     or span[1]-span[0]!=m or census.get('chunk_markers')!=chunk_markers):
  raise ValueError('Chunk census parent range or geometry mismatch')
 ranges=census_chunk_ranges(census,chunk_markers)
 chunks=census['chunks']
 if not isinstance(chunks,list) or len(chunks)!=len(ranges):
  raise ValueError('Chunk census does not cover the complete scan')
 aggregate=decoder_work(census,implementation,restart_ld_bases=True)
 units=Counter();forms=Counter();payload=0;parts=[];replays=[]
 for index,chunk in enumerate(chunks):
  lo,hi=ranges[index]
  if (chunk.get('variant_range')!=[lo,hi] or chunk.get('samples')!=n
      or chunk.get('markers')!=hi-lo or chunk.get('chunk_markers')!=chunk_markers):
   raise ValueError('Chunk census ranges, dimensions or geometry mismatch')
  for key in ('path','file_markers','file_bytes','index_bytes'):
   if key in census and chunk.get(key)!=census[key]:
    raise ValueError('Chunk census file context mismatch: '+key)
  amount=chunk.get('record_payload_bytes')
  if isinstance(amount,bool) or not isinstance(amount,int) or amount<0:
   raise ValueError('Chunk census payload must be a nonnegative integer')
  work=decoder_work(chunk,implementation,restart_ld_bases=True)
  units.update(work['source_units']);forms.update({int(k):v for k,v in chunk['record_form_counts'].items()})
  payload+=amount;parts.append(dict(census=chunk,decoder=work))
  replays.extend(chunk.get('native_ld_replays',[]))
 if (dict(units)!=aggregate['source_units'] or
     dict(forms)!={int(k):v for k,v in census['record_form_counts'].items()} or
     payload!=census['record_payload_bytes']):
  raise ValueError('Chunk census does not conserve aggregate encoded work')
 if replays!=census.get('native_ld_replays',[]):raise ValueError('Chunk census does not conserve native LD replays')
 return parts

def price_identified_units(work, cpu_seconds_per_unit):
 """Partial CPU service with missing costs explicit, never silently zeroed."""
 terms={};unpriced=list(work['uncounted_mechanisms'])
 for name,amount in work['source_units'].items():
  if name not in cpu_seconds_per_unit:unpriced.append(name);continue
  cost=positive(name,cpu_seconds_per_unit[name],True)
  terms[name]=amount*cost
 return dict(cpu_seconds_by_unit=terms,identified_cpu_seconds=sum(terms.values()),
  unpriced_source_units=unpriced,
  scope='Independent primitive-based partial CPU service; additive primitive timings approximate instruction scheduling. This omits physical traffic and is not a runtime prediction or a guaranteed lower bound.')


def verified_native_unit_prices(onebit, expansion, verification, *, library_sha256,
                                source_sha256, cpu_affinity):
 """Accept independently measured units only for the verified native binary.

 This validates provenance, not universal service bounds. Untested decoder
 units keep their previous explicit approximation; no scan observations enter.
 """
 import math
 import statistics
 libraries=verification['libraries']
 if (verification.get('executable_text_identical') is not True or
     libraries['production']['text_sha256']!=libraries['rebuilt']['text_sha256']):
  raise ValueError('Production decoder executable text was not verified')
 for digest in (libraries['production']['sha256'],onebit['production_library_sha256'],expansion['library_sha256']):
  if digest!=library_sha256:raise ValueError('Native decoder binary differs from primitive provenance')
 for report in (onebit,expansion,verification):
  if report['source_sha256']['native/pgen_decode.c']!=source_sha256:
   raise ValueError('Native decoder source differs from primitive provenance')
 if onebit['compiler']!=verification['compiler']:
  raise ValueError('Isolated loop compiler differs from verified rebuild')
 optimization=[flag for flag in verification['compile_command'] if flag.startswith('-O')]
 if len(optimization)!=1 or optimization[0] not in ('-O2','-O3'):
  raise ValueError('Unsupported or conflicting verified native optimization flags')
 flag=optimization[0];level=flag[1:]
 required={'-std=c99',flag,'-fPIC','-shared','-fno-strict-aliasing','-march=native'}
 commands=[command for command in onebit['commands'] if flag in command]
 if len(commands)!=1 or not required.issubset(commands[0]) or not required.issubset(verification['compile_command']):
  raise ValueError('Native primitive compiler flags do not match the verified configuration')
 for command in (commands[0],verification['compile_command']):
  if [option for option in command if option.startswith('-O')]!=[flag] or '-fstrict-aliasing' in command:
   raise ValueError('Conflicting native primitive compiler flags')
  allowed=required|{'-o','-Wall','-Wextra','-Wpedantic'}
  if any(option.startswith('-') and option not in allowed for option in command):
   raise ValueError('Unverified native primitive compiler flags')
 if level not in onebit['libraries']:
  raise ValueError('Missing primitive context for verified optimization level')
 for report in (onebit['libraries'][level],expansion):
  for when in ('cpu_before','cpu_after'):
   if report[when]['affinity']!=cpu_affinity:
    raise ValueError('Native primitive CPU affinity differs from requested context')
 if onebit.get('loop_units_per_call')!=1024:
  raise ValueError('Unexpected isolated onebit primitive extent')
 loop=[row for row in onebit['rows'] if row['compiler_level']==level and row['primitive']=='onebit8_native']
 expand=[row for row in expansion['rows'] if row['primitive']=='expand4_int8']
 prices={}
 for name,rows,extent_key in [('onebit8_native',loop,'loop_units_per_call'),('expand4_int8',expand,'packed_bytes_per_call')]:
  if len(rows)<3:raise ValueError('Insufficient independent primitive repeats')
  if any(row[extent_key]!=1024 for row in rows):raise ValueError('Mixed primitive extents')
  values=[row['cpu_seconds_per_unit'] for row in rows]
  if any(isinstance(value,bool) or not math.isfinite(value) or value<=0 for value in values):
   raise ValueError('Invalid independent native primitive CPU price')
  prices[name]=statistics.median(values)
 return dict(prices=prices,production_library_sha256=library_sha256,compiler_level=level,
             unpriced_terms=['Compiler and live-memory equivalence of remaining decoder primitive units'],
             scope='Independent fixed-unit native pricing with binary/source/compiler/context provenance; external call overhead remains amortized over1024 units. Not a workload timing fit or guaranteed bound.')
