"""Binary writer work from sumstats.py block/borrow/close semantics."""
from __future__ import annotations
import math
from .binary_writeback_work import binary_writeback_work


def _compact_stream(markers, traits, chunk_markers, block_bytes, *, borrow_chunks,
                    writeback_bytes, sync_file_range):
 """Count a fixed-chunk stream without materializing chunks or write events."""
 regular, tail = divmod(markers, chunk_markers)
 payload = chunk_markers*traits*4
 copied = calls = events = filled = 0
 if regular:
  if borrow_chunks and payload >= block_bytes:
   events = regular
  else:
   copied = regular*payload
   full, filled = divmod(copied, block_bytes)
   # Every chunk makes one slice assignment, and each block boundary inside
   # a chunk makes another. Boundaries coinciding with chunk ends add none.
   period = block_bytes//math.gcd(payload, block_bytes)
   calls = regular + full - regular//period
   events = full
 if tail:
  amount = tail*traits*4
  if borrow_chunks and filled == 0 and amount >= block_bytes:
   events += 1
  else:
   copied += amount
   calls += 1 + (filled+amount-1)//block_bytes
   full, filled = divmod(filled+amount, block_bytes)
   events += full
 events += bool(filled)
 total = markers*traits*4
 active = bool(sync_file_range and writeback_bytes)
 submits = total//writeback_bytes if active else 0
 waits = max(0, submits-1)
 return dict(payload_bytes=total, staging_copy_bytes=copied,
             staging_copy_calls=calls, write_calls_minimum=events,
             payload_queued_only_at_close_bytes=filled,
             writeback_submit_calls=submits, writeback_wait_calls=waits,
             writeback_fadvise_calls=waits,
             writeback_submitted_bytes=submits*writeback_bytes,
             writeback_waited_bytes=waits*writeback_bytes,
             writeback_unsubmitted_tail_bytes=total-submits*writeback_bytes)


def compact_binary_output_work(markers, traits=1, chunk_markers=2048,
                               block_bytes=None, queue_depth=3,
                               borrow_chunks=True, store_beta=True, fsync=True,
                               writeback_bytes=64<<20, sync_file_range=True,
                               store_variant_df=False):
 """Exact aggregate dense-writer counts in constant work per stream.

 This covers fixed chunks plus one tail. It deliberately omits the event
 order needed by a finite writer graph, short writes and durability service.
 """
 for name, value in [('markers', markers), ('traits', traits),
                     ('chunk_markers', chunk_markers), ('queue_depth', queue_depth)]:
  if type(value) is not int or value < (0 if name == 'markers' else 1):
   raise ValueError('Invalid '+name)
 if block_bytes is not None and (type(block_bytes) is not int or block_bytes < 1):
  raise ValueError('Invalid block size')
 if type(writeback_bytes) is not int or writeback_bytes < 0:
  raise ValueError('Invalid writeback interval')
 for name, value in [('borrow_chunks',borrow_chunks),('store_beta',store_beta),
                     ('fsync',fsync),('sync_file_range',sync_file_range),
                     ('store_variant_df',store_variant_df)]:
  if type(value) is not bool:raise ValueError('Boolean '+name+' required')
 memory = binary_output_memory(markers, traits, chunk_markers, block_bytes,
                               queue_depth, store_beta, store_variant_df)
 size = memory['block_bytes']
 streams = {}
 for name in (['beta','t'] if store_beta else ['t']):
  streams[name] = _compact_stream(markers, traits, chunk_markers, size,
                                 borrow_chunks=borrow_chunks,
                                 writeback_bytes=writeback_bytes,
                                 sync_file_range=sync_file_range)
 if store_variant_df:
  streams['df'] = _compact_stream(markers, 1, chunk_markers, min(size,1<<20),
                                  borrow_chunks=borrow_chunks,
                                  writeback_bytes=writeback_bytes,
                                  sync_file_range=sync_file_range)
 result = dict(kind='torchgwas.compact_binary_output_work.v1',
               markers=markers, traits=traits, chunk_markers=chunk_markers,
               block_bytes=size, queue_depth=queue_depth,
               arrays=len(streams), streams=streams,
               allocated_staging_bytes=memory['allocated_staging_bytes'],
               zero_initialization_bytes=memory['zero_initialization_bytes'],
               fsync_calls=len(streams) if fsync else 0,
               scope='Exact fixed-chunk dense writer payload, staging, minimum write and writeback request counts without event expansion. Short writes, writer CPU/storage service, queue contention, final metadata and durable completion remain unpriced.')
 for key in ('payload_bytes','staging_copy_bytes','staging_copy_calls',
             'write_calls_minimum','payload_queued_only_at_close_bytes',
             'writeback_submit_calls','writeback_wait_calls','writeback_fadvise_calls',
             'writeback_submitted_bytes','writeback_waited_bytes',
             'writeback_unsubmitted_tail_bytes'):
  result[key] = sum(row[key] for row in streams.values())
 return result

def binary_output_memory(markers,traits=1,chunk_markers=2048,block_bytes=None,
                         queue_depth=3,store_beta=True,store_variant_df=False):
 """Exact staging allocation without constructing chunks or write events.

 The first payload selects beta/t block size; df uses its own 1 MiB cap.
 Every stream allocates queue_depth+1 blocks even when chunks are borrowed.
 Empty writers open default-sized streams on close. Python/queue objects,
 caller-owned results and the filesystem page cache are outside this ledger.
 """
 for name,value in [('markers',markers),('traits',traits),('chunk_markers',chunk_markers),('queue_depth',queue_depth)]:
  if not isinstance(value,int) or value<(0 if name=='markers' else 1):raise ValueError('Invalid '+name)
 if block_bytes is not None and (not isinstance(block_bytes,int) or block_bytes<1):raise ValueError('Invalid block size')
 # Empty writer opens default-sized streams at close. Otherwise the initial
 # chunk selects a minimum 1 MiB and maximum 16 MiB auto block per array.
 first=min(markers,chunk_markers)*traits*4
 size=block_bytes if block_bytes is not None else (min(16<<20,max(1<<20,first)) if markers else 16<<20)
 streams={name:size for name in (['beta','t'] if store_beta else ['t'])}
 if store_variant_df:streams['df']=min(size,1<<20)
 allocated=(queue_depth+1)*sum(streams.values())
 return dict(block_bytes=size,queue_depth=queue_depth,arrays=len(streams),
  block_bytes_by_stream=streams,pooled_buffers_per_array=queue_depth+1,
  allocated_staging_bytes=allocated,zero_initialization_bytes=allocated)


def binary_output_work(markers,traits=1,chunk_markers=2048,block_bytes=None,
                       queue_depth=3,borrow_chunks=True,store_beta=True,fsync=True,
                       writeback_bytes=64<<20,sync_file_range=True,store_variant_df=False,*,chunk_rows=None):
 if chunk_rows is not None:
  if (type(chunk_markers) is not int or chunk_markers<1 or not isinstance(chunk_rows,(list,tuple))
      or any(type(v) is not int or not 0<v<=chunk_markers for v in chunk_rows)
      or sum(chunk_rows)!=markers):
   raise ValueError('Explicit writer chunks must cover all markers within capacity')
 first=chunk_rows[0] if chunk_rows else chunk_markers
 memory=binary_output_memory(markers,traits,first,block_bytes,queue_depth,store_beta)
 size=memory['block_bytes']
 filled=0;events=[];copies=0;copy_calls=0;chunks=[]
 counts=(chunk_rows if chunk_rows is not None else
         (min(chunk_markers,markers-start) for start in range(0,markers,chunk_markers)))
 for index,count in enumerate(counts):
  payload=count*traits*4;copy=0;calls=0
  if borrow_chunks and filled==0 and payload>=size:
   events.append(dict(chunk=index,bytes=payload,pooled=False,at_close=False))
  else:
   remaining=payload
   while remaining:
    take=min(size-filled,remaining);filled+=take;remaining-=take;copy+=take;calls+=1
    if filled==size:
     events.append(dict(chunk=index,bytes=size,pooled=True,at_close=False));filled=0
  copies+=copy;copy_calls+=calls
  chunks.append(dict(markers=count,payload_bytes_per_array=payload,staging_copy_bytes_per_array=copy,staging_copy_calls_per_array=calls))
 if filled:events.append(dict(chunk=len(chunks),bytes=filled,pooled=True,at_close=True))
 arrays=memory['arrays']
 result=dict(markers=markers,traits=traits,chunk_markers=chunk_markers,arrays=arrays,
  block_bytes=size,queue_depth=queue_depth,pooled_buffers_per_array=queue_depth+1,
  allocated_staging_bytes=memory['allocated_staging_bytes'],zero_initialization_bytes=memory['zero_initialization_bytes'],
  binary_payload_bytes=arrays*markers*traits*4,staging_copy_bytes=arrays*copies,
  staging_copy_calls=arrays*copy_calls,staging_view_creations=2*arrays*copy_calls,
  staging_copy_mechanism='memoryview_slice',
  staging_logical_memory_bytes=2*arrays*copies,binary_write_calls_minimum=arrays*len(events),
  payload_queued_only_at_close_bytes=arrays*sum(e['bytes'] for e in events if e['at_close']),
  fsync_calls=arrays if fsync else 0,chunks=chunks,events_per_array=events,
  writeback_per_array=binary_writeback_work(events,writeback_bytes,sync_file_range),
  scope='Exact source-level work for owned contiguous FP32 chunks. Short writes can add syscalls. Sidecars, manifest serialization, allocator service, page-cache writeback and storage commit latency need separate explicit service inputs.')
 if store_variant_df:
  primary=dict(result)
  df=binary_output_work(markers,1,chunk_markers,min(size,1<<20),queue_depth,
      borrow_chunks,False,fsync,writeback_bytes,sync_file_range,chunk_rows=chunk_rows)
  result['stream_work']={name:primary for name in (['beta','t'] if store_beta else ['t'])}
  result['stream_work']['df']=df
  result['arrays']+=1
  for name in ['allocated_staging_bytes','zero_initialization_bytes','binary_payload_bytes',
               'staging_copy_bytes','staging_copy_calls','staging_view_creations','staging_logical_memory_bytes','binary_write_calls_minimum',
               'payload_queued_only_at_close_bytes','fsync_calls']:
   result[name]+=df[name]
  result['df_payload_bytes']=df['binary_payload_bytes']
  result['scope']='Exact per-stream beta/t/variant-df work, including the independent 1 MiB df block cap. Manifest serialization, validation, atomic replacement/directory fsync and allocator/storage latency require separate service.'
 return result
