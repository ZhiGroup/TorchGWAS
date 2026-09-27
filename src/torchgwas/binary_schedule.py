"""Finite production binary-writer dependencies, with externally priced service.

Counts stage-copy/queue/pool-credit/borrow-credit/ordered close mechanisms.
Writer service must account for shared storage capacity; the DAG does not
invent independent full-device bandwidth for two simultaneous streams.
"""
from __future__ import annotations
from .binary_output_work import binary_output_work
from .binary_writeback_work import binary_writeback_work
from .binary_writeback_schedule import RangeWritebackSchedule

class BinaryWriterSchedule:
 def __init__(self,work,*,copy_seconds_per_byte,zero_seconds_per_byte,
              write_seconds_per_byte,fsync_seconds_per_array,
              append_seconds=0.,handoff_seconds=0.,open_seconds=0.,
              inflight_bytes=512<<20,writeback_bytes=64<<20,cpu_fraction=None,write_capacity=None,writeback_service=None,
              copy_seconds_per_call=0.):
  stream_work=work.get('stream_work') or {a:work for a in (['beta','t'] if work['arrays']==2 else ['t'])}
  if writeback_service is None and writeback_bytes and any(w['binary_payload_bytes']//w['arrays']>=writeback_bytes for w in stream_work.values()):
   raise ValueError('Periodic sync_file_range needs an explicit storage-writeback schedule')
  self.cpu_resources={'cpu':cpu_fraction} if cpu_fraction else None
  self.write_resources={'output':write_capacity} if write_capacity else None
  self.work=work;self.copy=copy_seconds_per_byte;self.zero=zero_seconds_per_byte
  self.copy_call=copy_seconds_per_call
  self.write=write_seconds_per_byte;self.fsync=fsync_seconds_per_array
  self.append_cost=append_seconds;self.handoff=handoff_seconds;self.open=open_seconds
  for value in [self.copy,self.copy_call,self.zero,self.write,self.fsync,self.append_cost,self.handoff,self.open]:
   if value<0:raise ValueError('Negative writer service')
  self.stream_work=stream_work;self.arrays=list(stream_work)
  self.range_schedules={}
  if writeback_service is not None:
   if not work['fsync_calls']:raise ValueError('Explicit writeback schedule currently requires closing fsync')
   self.range_schedules={a:RangeWritebackSchedule(binary_writeback_work(w['events_per_array'],writeback_bytes),writeback_service,cpu_resources=self.cpu_resources) for a,w in stream_work.items()}
   self.write=writeback_service['pagecache_seconds_per_byte']
   self.write_resources=self.cpu_resources
  self.states={a:dict(filled=0,writes=[],pooled=[],borrowed=[],serial=0) for a in self.arrays}
  self.credits={a:max(inflight_bytes,2*w['block_bytes']) for a,w in stream_work.items()};self.opened=False
  self.payload=0;self.copied=0;self.copy_calls=0;self.write_calls=0;self.close_payload=0
 def _open(self,g,after):
  name=g.add('writer:open',self.open+self.work['zero_initialization_bytes']*self.zero,after,resources=self.cpu_resources)
  self.opened=True;return name
 def _queue(self,g,a,size,pooled,after,at_close=False):
  st=self.states[a];j=st['serial'];st['serial']+=1;deps=[after]
  # Borrowed blocks have a byte credit budget, independent of pooled buffers.
  # The source permits an oversized single block only when no previous block
  # is in flight. Completion order is FIFO within each writer thread.
  if not pooled:
   outstanding=size
   for prior_size,prior_write in reversed(st['borrowed']):
    if outstanding+prior_size>self.credits[a]:
     deps.append(prior_write);break
    outstanding+=prior_size
  queued=g.add(f'writer:{a}:queue:{j}',self.handoff,deps,resources=self.cpu_resources)
  written=g.add(f'writer:{a}:write:{j}',size*self.write,[queued]+st['writes'][-1:],resources=self.write_resources)
  if a in self.range_schedules:written=self.range_schedules[a].after_write(g,a,j,written)
  st['writes'].append(written);self.payload+=size;self.write_calls+=1
  if at_close:self.close_payload+=size
  if pooled:
   count=len(st['pooled']);st['pooled'].append(written)
   # depth+1 total buffers: one staging, depth available before first handoff.
   wait=[queued]
   if count>=self.work['queue_depth']:wait.append(st['pooled'][count-self.work['queue_depth']])
   return g.add(f'writer:{a}:reclaim:{j}',0.,wait)
  st['borrowed'].append((size,written));return queued
 def append(self,g,index,after):
  if index>=len(self.work['chunks']):raise ValueError('Writer/scan chunk mismatch')
  last=g.add(f'writer:append:{index}',self.append_cost,after,resources=self.cpu_resources)
  if not self.opened:last=self._open(g,[last])
  for a in self.arrays:
   work=self.stream_work[a];chunk=work['chunks'][index];size=chunk['payload_bytes_per_array'];block=work['block_bytes']
   st=self.states[a]
   # Work ledger records a borrowed full chunk as zero staging-copy bytes.
   if size and chunk['staging_copy_bytes_per_array']==0:
    if st['filled']:raise ValueError('Borrow cannot bypass partially staged block')
    last=self._queue(g,a,size,False,last);continue
   remaining=size;part=0
   while remaining:
    take=min(block-st['filled'],remaining)
    last=g.add(f'writer:{a}:copy:{index}:{part}',take*self.copy+self.copy_call,[last],resources=self.cpu_resources);part+=1
    self.copy_calls+=1
    self.copied+=take;st['filled']+=take;remaining-=take
    if st['filled']==block:
     last=self._queue(g,a,block,True,last);st['filled']=0
  return g.add(f'consume:{index}',0.,[last])
 def close(self,g,after):
  last=g.add('writer:close:start',0.,after)
  if not self.opened:last=self._open(g,[last])
  # Source closes beta completely (join, fsync) before t.close() is called.
  # Thus t's partial staging block is not enqueued until beta's fsync returns.
  for a in self.arrays:
   st=self.states[a]
   if st['filled']:
    last=self._queue(g,a,st['filled'],True,last,at_close=True);st['filled']=0
   last=g.add(f'writer:{a}:drain',0.,[last]+st['writes'][-1:])
   if a in self.range_schedules:last=self.range_schedules[a].before_fsync(g,a,last)
   last=g.add(f'writer:{a}:fsync',self.fsync if self.work['fsync_calls'] else 0.,[last])
  if self.payload!=self.work['binary_payload_bytes'] or self.copied!=self.work['staging_copy_bytes']:
   raise ValueError('Writer payload/copy work not conserved')
  if self.copy_calls!=self.work['staging_copy_calls']:
   raise ValueError('Writer copy-call work not conserved')
  if self.close_payload!=self.work['payload_queued_only_at_close_bytes']:
   raise ValueError('Writer close-only work not conserved')
  return last
