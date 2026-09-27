"""Binary output for selected pairs and joint statistics.

Full scans retain the dense beta/t store. Reductions use indexed NumPy parts,
so variable-length results need neither text nor a dense M-by-K allocation.
"""
from pathlib import Path
from dataclasses import dataclass
import json,os,time
import numpy as np
from scipy import special


@dataclass(frozen=True)
class IndexedOutputPartition:
 """Producer identity in source-file variant and post-QC phenotype coordinates."""
 device: str
 variant_range: tuple[int,int]
 trait_range: tuple[int,int]

 def __post_init__(self):
  if not isinstance(self.device,str) or not self.device:
   raise ValueError('Explicit output producer device required')
  for name in ('variant_range','trait_range'):
   span=getattr(self,name)
   if not isinstance(span,tuple) or len(span)!=2 or any(type(v) is not int for v in span) or not 0<=span[0]<span[1]:
    raise ValueError('Immutable nonempty output partition ranges required')


@dataclass(frozen=True)
class PartitionedIndexedChunk:
 """Carry producer identity through a shared queue without copying arrays."""
 payload: tuple
 partition: IndexedOutputPartition

 def __post_init__(self):
  if not isinstance(self.payload,tuple) or not isinstance(self.partition,IndexedOutputPartition):
   raise ValueError('Tuple payload and typed output partition required')

 def __len__(self):return len(self.payload)
 def __getitem__(self,key):return self.payload[key]
 def __iter__(self):return iter(self.payload)
 def with_payload(self,payload):return PartitionedIndexedChunk(payload,self.partition)


@dataclass(frozen=True)
class IndexedChunkWrite:
 """One consumed indexed chunk, after optional part-file fsync.

 Empty chunks create no part. File fsync is distinct from final manifest and
 directory durability. This event is not a scan-queue acceptance timestamp.
 """
 start: int
 end: int
 kind: str
 rows: int
 part_bytes: int
 part_file: str | None
 started: float
 completed: float
 part_file_fsynced: bool
 partition: IndexedOutputPartition | None = None
 source_variant_range: tuple[int,int] | None = None


COALESCE_ROWS=1<<18  # rows per part when no per-chunk durability is observed
MERGE_WINDOW_VARIANTS=1<<16  # variant span merged at a time when parts overlap


def write_indexed_sumstats(directory,marker_names,trait_names,n_samples,chunks,
                           *,kind,df,chi2_df=None,extra_manifest=None,
                           p_value_threshold=None,variant_metadata=None,fsync=True,store_beta=True,before_publish=None,
                           on_chunk_written=None,partition_for_range=None,variant_offset=0,live_progress=None,
                           coalesce_rows=None,variant_source=None,embed_variant_ids=None):
 """variant_source (variant_source.variant_source_record): rows keep their
 variant_index and the manifest records the input they index; the IDs and
 variant metadata are then the genotype's own and are not rewritten.
 embed_variant_ids forces variant_ids.npy (and variant_metadata.npz) into the
 store; the default embeds only when no source is recorded."""
 if on_chunk_written is not None and (not callable(on_chunk_written) or kind not in ('jagwas','significant')):
  raise ValueError('Indexed chunk completion requires a callback and streaming jagwas/significant output')
 if live_progress is not None:
  from .indexed_writer_progress import IndexedWriterProgress
  if (not isinstance(live_progress,IndexedWriterProgress) or
      live_progress.kind != kind or on_chunk_written is None or
      (partition_for_range is None and kind!='significant')):
   raise ValueError('Live indexed writer requires a bound observed reduction')
 if type(variant_offset) is not int or variant_offset<0:
  raise ValueError('Nonnegative source variant offset required')
 if partition_for_range is not None and (not callable(partition_for_range) or on_chunk_written is None):
  raise ValueError('Partition resolver requires writer completion observation')
 directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
 (directory/'manifest.json').unlink(missing_ok=True)
 parts=[];total=0;started=time.perf_counter()
 # Pair stores publish in (variant, trait) order. The internal per-variant
 # reduction keeps its own row order (descending |t| within a variant).
 ordered_rows=kind in ('significant','filtered')
 def emit(values,start,end):
  nonlocal total
  if not store_beta: values.pop('beta', None)
  count=len(values['variant_index'])
  if not count:return None
  if ordered_rows and count>1:
   # Device selection walks blocks of variants x trait strips; store each
   # part in (variant, trait) order so the published store can be.
   order=np.lexsort((np.asarray(values['trait_index']),np.asarray(values['variant_index'])))
   values={key:np.asarray(value)[order] for key,value in values.items()}
  name=f'part_{len(parts):06d}.npz'
  with (directory/name).open('wb') as f:
   np.savez(f,**values);f.flush()
   if fsync:os.fsync(f.fileno())
   size=f.tell() if on_chunk_written is not None else 0
  parts.append({'file':name,'rows':count,'variant_range':[int(start),int(end)]});total+=count
  return (name,count,size) if on_chunk_written is not None else None
 # coalesce_rows (callers without a per-chunk durability observer): write one
 # part and one fsync per coalesce_rows rows instead of per chunk (196 fsyncs
 # cost ~6 s on lab-2080ti's output disk). The default stays one part per
 # chunk, which is what the JIT calculator's write model describes. A lane is
 # one producer's contiguous run: a chunk extends the lane that ends where it
 # starts, so interleaved variant shards or phenotype tiles keep their own.
 if coalesce_rows is not None and (on_chunk_written is not None or type(coalesce_rows) is not int or coalesce_rows<1):
  raise ValueError('coalesce_rows needs a positive row count and no per-chunk observer')
 lanes={}
 def flush(lane):
  if lane['values']:
   emit({key:np.concatenate([v[key] for v in lane['values']]) for key in lane['values'][0]},
        lane['start'],lane['end'])
 def add(values,start,end):
  if coalesce_rows is None:return emit(values,start,end)
  waiting=lanes.get(start)
  lane=waiting.pop() if waiting else dict(start=start,end=start,values=[],rows=0)
  if waiting is not None and not waiting:del lanes[start]
  count=len(values['variant_index'])
  if count:lane['values'].append(values);lane['rows']+=count
  lane['end']=end
  if lane['rows']>=coalesce_rows:
   flush(lane);lane=dict(start=end,end=end,values=[],rows=0)
  lanes.setdefault(end,[]).append(lane)
  return None
 def completed(start,end,began,part,partition):
  if on_chunk_written is not None:
   name,count,size=part if part is not None else (None,0,0)
   if live_progress is not None:
    live_progress.complete_chunk(count,size,bool(fsync and name is not None))
   on_chunk_written(IndexedChunkWrite(int(start),int(end),kind,count,size,name,
                    began,time.perf_counter(),bool(fsync and name is not None),partition,
                    (int(start)+variant_offset,int(end)+variant_offset)))
 critical=(abs(float(special.stdtrit(df,p_value_threshold/2))) if p_value_threshold is not None else None)
 iterator=iter(chunks)
 try:
  for chunk in iterator:
   start,end=chunk[:2]
   began=time.perf_counter() if on_chunk_written is not None else None
   partition=None
   if on_chunk_written is not None:
    source_span=(int(start)+variant_offset,int(end)+variant_offset)
    partition=chunk.partition if isinstance(chunk,PartitionedIndexedChunk) else None
    if partition_for_range is not None:
     resolved=partition_for_range(*source_span)
     if not isinstance(resolved,IndexedOutputPartition):raise ValueError('Typed indexed output partition required')
     if partition is not None and partition!=resolved:raise ValueError('Conflicting indexed output partition')
     partition=resolved
    if partition is not None:
     if (not isinstance(partition,IndexedOutputPartition) or not
         partition.variant_range[0]<=source_span[0]<source_span[1]<=partition.variant_range[1] or
         partition.trait_range[1]>len(trait_names)):
      raise ValueError('Indexed output lies outside its producer partition')
     if kind=='jagwas' and partition.trait_range!=(0,len(trait_names)):
      raise ValueError('JAGWAS output partition must retain the complete phenotype panel')
    if live_progress is not None:
     live_progress.begin_chunk(source_span,partition)
   if kind=='jagwas':
    values=np.asarray(chunk[3],dtype=np.float64).reshape(-1);keep=np.flatnonzero(np.isfinite(values))
    if live_progress is not None:live_progress.begin_emit()
    part=add({'variant_index':start+keep,'chi2':values[keep]},start,end)
    completed(start,end,began,part,partition);continue
   if kind=='significant':
    _,_,vi,ti,beta,t,_chunk_df=chunk
    if store_beta and beta is None:raise ValueError('Stored beta output requires supplied effect sizes')
    values={'variant_index':np.asarray(vi,dtype=np.int64),'trait_index':np.asarray(ti,dtype=np.int64)}
    if store_beta:values['beta']=np.asarray(beta)
    values.update(t_stat=np.asarray(t),df=np.broadcast_to(np.asarray(_chunk_df),np.asarray(t).shape).copy())
    if live_progress is not None:live_progress.begin_emit()
    part=add(values,start,end)
    completed(start,end,began,part,partition)
    continue
   beta=np.asarray(chunk[2]);t=np.asarray(chunk[3]);keep=np.isfinite(t)
   if critical is not None:keep &= np.abs(t)>=critical
   vi,ti=np.nonzero(keep)
   trait_index=np.asarray(chunk[5])[vi,ti] if kind=='reduced' else ti
   add({'variant_index':(start+vi).astype(np.int64),'trait_index':np.asarray(trait_index,dtype=np.int64),'beta':beta[vi,ti],'t_stat':t[vi,ti]},start,end)
  for waiting in list(lanes.values()):
   for lane in waiting:flush(lane)
 except BaseException:
  if live_progress is not None:live_progress.failed()
  raise
 finally:
  if hasattr(iterator,'close'):iterator.close()
 if live_progress is not None:live_progress.begin_finalization()
 # Per-chunk mode (the JIT calculator's write model; observed runs have also
 # reported their files) keeps parts as written and only orders the list.
 # Publication steps are timed: they follow the scan, so their fsyncs are
 # exposed to other writers' data on the same filesystem (ext4 data=ordered).
 publication={};step=time.perf_counter()
 parts,ordering=_order_parts(directory,parts,fsync,merge=ordered_rows and coalesce_rows is not None)
 publication['order_parts_seconds']=time.perf_counter()-step
 embed=variant_source is None if embed_variant_ids is None else bool(embed_variant_ids)
 if embed:
  step=time.perf_counter()
  with (directory/'variant_ids.npy').open('wb') as handle:
   np.save(handle,np.asarray(marker_names,dtype=str),allow_pickle=False);handle.flush()
   if fsync:os.fsync(handle.fileno())
  if variant_metadata is not None:
   with (directory/'variant_metadata.npz').open('wb') as handle:
    np.savez(handle,**{key:np.asarray(v,dtype=np.int64 if key=='position' else str) for key,v in variant_metadata.items()});handle.flush()
    if fsync:os.fsync(handle.fileno())
  publication['variant_ids_seconds']=time.perf_counter()-step
 manifest={'format':'torchgwas-indexed-sumstats','version':1,'kind':'jagwas' if kind=='jagwas' else 'linear','shape':[len(marker_names),len(trait_names)],'n_samples':int(n_samples),'df':int((chi2_df() if callable(chi2_df) else chi2_df) if kind=='jagwas' else df),'traits':list(trait_names),'parts':parts,'rows':total,'p_value':'derived from chi2 and df' if kind=='jagwas' else 'two-sided Student t on t_stat with df'}
 if embed:manifest['variant_ids']='variant_ids.npy'
 if variant_source is not None:manifest['variant_source']=dict(variant_source)
 if extra_manifest is not None:
  manifest.update(extra_manifest() if callable(extra_manifest) else extra_manifest)
 if kind=='significant':
  manifest.update(version=2,df={'layout':'per_part','axis':'pair','field':'df'},nominal_df=int(df),p_value='two-sided Student t on t_stat with per-pair df in each part')
 if (ordered_rows or kind=='jagwas') and ordering['globally_ordered']:
  manifest['row_order']='variant_index then trait_index, across parts in manifest order'
 if before_publish is not None:before_publish()
 if live_progress is not None:live_progress.begin_publish()
 from .sumstats import write_manifest
 step=time.perf_counter()
 write_manifest(directory,manifest,fsync=fsync)
 publication['manifest_seconds']=time.perf_counter()-step
 if live_progress is not None:live_progress.published()
 return total,{'cells':total,'directory':str(directory),'indexed':True,'write_seconds':time.perf_counter()-started,
               'ordering':ordering,'publication':publication,'variant_ids_embedded':embed}


def _order_parts(directory,parts,fsync,*,merge=True):
 """Publish parts so rows run in (variant, trait) order across the store.

 Every part is sorted internally and covers one producer's contiguous
 variant run. Disjoint ranges (one GPU, variant shards, JAGWAS) only need the
 part list sorted, which moves no data. Phenotype tiles write overlapping
 ranges; those are merged in fixed variant windows, so memory is bounded by
 the parts that intersect one window, and rewritten once as large parts,
 before the manifest exists.
 """
 began=time.perf_counter()
 parts=sorted(parts,key=lambda part:tuple(part['variant_range']))
 overlapping=any(b['variant_range'][0]<a['variant_range'][1] for a,b in zip(parts,parts[1:]))
 if not merge or not overlapping:
  return parts,{'merged':False,'globally_ordered':not overlapping,'seconds':time.perf_counter()-began}
 def load(part):
  with np.load(directory/part['file'],allow_pickle=False) as values:
   return {key:values[key] for key in values.files}
 merged,pending,index,out=[],[],0,[]
 def write(rows,first,last):
  name=f'ordered_{len(merged):06d}.npz'
  with (directory/name).open('wb') as f:
   np.savez(f,**rows);f.flush()
   if fsync:os.fsync(f.fileno())
  merged.append({'file':name,'rows':int(len(rows['variant_index'])),'variant_range':[int(first),int(last)]})
 def drain(last):
  if out:
   write({key:np.concatenate([rows[key] for _,rows in out]) for key in out[0][1]},out[0][0],last)
   out.clear()
 buffered=0
 cursor=parts[0]['variant_range'][0]
 while index<len(parts) or pending:
  if not pending:cursor=max(cursor,parts[index]['variant_range'][0])
  end=cursor+MERGE_WINDOW_VARIANTS
  while index<len(parts) and parts[index]['variant_range'][0]<end:
   pending.append(load(parts[index]));index+=1
  taken,rest=[],[]
  for values in pending:
   cut=int(np.searchsorted(values['variant_index'],end,side='left'))
   if cut:taken.append({key:value[:cut] for key,value in values.items()})
   if cut<len(values['variant_index']):rest.append({key:value[cut:] for key,value in values.items()})
  pending=rest
  if taken:
   rows={key:np.concatenate([values[key] for values in taken]) for key in taken[0]}
   order=np.lexsort((rows['trait_index'],rows['variant_index']))
   out.append((cursor,{key:value[order] for key,value in rows.items()}));buffered+=len(order)
   if buffered>=COALESCE_ROWS:
    drain(end);buffered=0
  cursor=end
 drain(cursor)
 for part in parts:(directory/part['file']).unlink()
 return merged,{'merged':True,'globally_ordered':True,'parts_written':len(parts),'parts_published':len(merged),
                'seconds':time.perf_counter()-began}


def open_indexed_sumstats(directory):
 """Return the manifest and an iterator of bounded binary result parts."""
 directory=Path(directory)
 manifest=json.loads((directory/'manifest.json').read_text())
 if manifest.get('format')!='torchgwas-indexed-sumstats':
  raise ValueError('not an indexed binary sumstats store')
 def parts():
  for part in manifest['parts']:
   with np.load(directory/part['file'],allow_pickle=False) as values:
    data={key:values[key] for key in values.files}
   if len(data['variant_index'])!=part['rows']:raise ValueError('indexed part row count mismatch')
   yield data
 return manifest,parts()
