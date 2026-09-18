"""Binary output for selected pairs and joint statistics.

Full scans retain the dense beta/t store. Reductions use indexed NumPy parts,
so variable-length results need neither text nor a dense M-by-K allocation.
"""
from pathlib import Path
import heapq,json,os,time
import numpy as np
from scipy import special
from .tails import upper_tail_log10_from_t

def write_indexed_sumstats(directory,marker_names,trait_names,n_samples,chunks,
                           *,kind,df,chi2_df=None,topk_per_trait=None,
                           p_value_threshold=None,variant_metadata=None,fsync=True,store_beta=True):
 directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
 parts=[];total=0;started=time.perf_counter()
 def emit(values):
  nonlocal total
  if not store_beta: values.pop('beta', None)
  count=len(values['variant_index'])
  if not count:return
  name=f'part_{len(parts):06d}.npz'
  with (directory/name).open('wb') as f:
   np.savez(f,**values);f.flush()
   if fsync:os.fsync(f.fileno())
  parts.append({'file':name,'rows':count});total+=count
 heaps=[[] for _ in trait_names] if topk_per_trait is not None else None
 critical=(abs(float(special.stdtrit(df,p_value_threshold/2))) if p_value_threshold is not None else None)
 for chunk in chunks:
  start,end=chunk[:2]
  if kind=='jagwas':
   values=np.asarray(chunk[3],dtype=np.float64).reshape(-1);keep=np.flatnonzero(np.isfinite(values))
   emit({'variant_index':start+keep,'chi2':values[keep]});continue
  if kind=='significant':
   _,_,vi,ti,beta,t,_chunk_df=chunk
   emit({'variant_index':np.asarray(vi,dtype=np.int64),'trait_index':np.asarray(ti,dtype=np.int64),'beta':np.asarray(beta),'t_stat':np.asarray(t),'neg_log10_p':upper_tail_log10_from_t(t,_chunk_df).astype(np.float32)})
   continue
  beta=np.asarray(chunk[2]);t=np.asarray(chunk[3]);keep=np.isfinite(t)
  if critical is not None:keep &= np.abs(t)>=critical
  vi,ti=np.nonzero(keep)
  trait_index=np.asarray(chunk[5])[vi,ti] if kind=='reduced' else ti
  if heaps is None:
   selected_t=t[vi,ti]
   emit({'variant_index':(start+vi).astype(np.int64),'trait_index':np.asarray(trait_index,dtype=np.int64),'beta':beta[vi,ti],'t_stat':selected_t,'neg_log10_p':upper_tail_log10_from_t(selected_t,df).astype(np.float32)})
  else:
   for v,j,trait in zip(vi,ti,trait_index):
    item=(abs(float(t[v,j])),int(start+v),float(beta[v,j]),float(t[v,j]))
    heap=heaps[int(trait)]
    if len(heap)<topk_per_trait:heapq.heappush(heap,item)
    elif item[0]>heap[0][0]:heapq.heapreplace(heap,item)
 if heaps is not None:
  values={'variant_index':[],'trait_index':[],'beta':[],'t_stat':[]}
  for trait,heap in enumerate(heaps):
   for _,variant,beta,t in sorted(heap,key=lambda x:(x[0],x[1]),reverse=True):
    values['variant_index'].append(variant);values['trait_index'].append(trait);values['beta'].append(beta);values['t_stat'].append(t)
  packed={key:np.asarray(v,dtype=np.int64 if key.endswith('_index') else np.float64) for key,v in values.items()}
  packed['neg_log10_p']=upper_tail_log10_from_t(packed['t_stat'],df).astype(np.float32)
  emit(packed)
 np.save(directory/'variant_ids.npy',np.asarray(marker_names,dtype=str),allow_pickle=False)
 if variant_metadata is not None:
  np.savez(directory/'variant_metadata.npz',**{key:np.asarray(v,dtype=np.int64 if key=='position' else str) for key,v in variant_metadata.items()})
 manifest={'format':'torchgwas-indexed-sumstats','version':2,'kind':'jagwas' if kind=='jagwas' else 'linear','shape':[len(marker_names),len(trait_names)],'n_samples':int(n_samples),'df':int(chi2_df if kind=='jagwas' else df),'traits':list(trait_names),'parts':parts,'rows':total,'variant_ids':'variant_ids.npy','significance':'chi-square tail derived from chi2 and df' if kind=='jagwas' else 'neg_log10_p is -log10 of the exact two-sided Student-t tail'}
 (directory/'manifest.json').write_text(json.dumps(manifest,indent=2))
 return total,{'cells':total,'directory':str(directory),'indexed':True,'write_seconds':time.perf_counter()-started}


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
