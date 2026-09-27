"""Durable contiguous variant shards, with one independent writer per GPU.

Borrowed scan buffers are copied into each shard's bounded writer before its
iterator advances. A root completion manifest joins the files without a full
matrix gather or a cross-device result queue.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
import threading
import time

import numpy as np

from .linear import multigpu_variant_ranges
from .sumstats import BinarySumstatsWriter,FORMAT_VERSION,MANIFEST_NAME,read_manifest,write_manifest


def write_variant_sharded_sumstats(directory, *, n_variants, trait_names, n_samples,
        df, chunk_size, devices, reader_workers, scan_factory, block_bytes=None,
        queue_depth=3, fsync=True, store_beta=True, before_publish=None, on_write_progress=None, on_writer_open=None,
        trait_df=None, extra_manifest=None):
    """scan_factory(start,end,device,workers) must yield shard-local row indices.

    trait_df: per-trait df for a panel with missing phenotypes, the
    single-device contract (manifest df list, no per-variant sidecar). The
    scan's sixth element is then the pair df of the imputed panel and is not
    stored. Without it every shard keeps a per-variant df sidecar.
    """
    if on_write_progress is not None and not callable(on_write_progress):
        raise ValueError('Dense write progress requires a callback')
    if on_writer_open is not None and not callable(on_writer_open):
        raise ValueError('Dense writer registration requires a callback')
    for name,value in [('n_variants',n_variants),('chunk_size',chunk_size),
                       ('reader_workers',reader_workers),('queue_depth',queue_depth)]:
        if isinstance(value,bool) or not isinstance(value,int) or value<1:
            raise ValueError(name+' must be a positive integer')
    devices=[str(device) for device in devices]
    if not devices or len(set(devices))!=len(devices):
        raise ValueError('Unique explicit variant devices required')
    if not trait_names:raise ValueError('At least one phenotype required')
    if trait_df is not None:
        trait_df=[int(value) for value in trait_df]
        if len(trait_df)!=len(trait_names):raise ValueError('One df per phenotype required')
    ranges=multigpu_variant_ranges(n_variants,chunk_size,len(devices))
    devices=devices[:len(ranges)]
    if reader_workers<len(devices):
        raise ValueError('reader_workers must provide at least one reader per active device')
    if block_bytes is not None and (isinstance(block_bytes,bool) or not isinstance(block_bytes,int) or block_bytes<1):
        raise ValueError('block_bytes must be a positive integer')
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    (directory/MANIFEST_NAME).unlink(missing_ok=True)
    stop=threading.Event();per_device,remainder=divmod(reader_workers,len(devices))
    failure=[];failure_lock=threading.Lock()
    def failed(error):
        with failure_lock:
            if not failure:failure.append(error)
        stop.set()

    def worker(index):
        start,end=ranges[index];device=devices[index];workers=per_device+(index<remainder)
        relative=f'variants_{start:012d}_{end:012d}'
        writer=None;iterator=None;began=time.perf_counter()
        try:
            if stop.is_set():raise RuntimeError('variant shard cancelled after peer failure')
            writer=BinarySumstatsWriter(directory/relative,end-start,list(trait_names),n_samples,
                df if trait_df is None else trait_df,
                block_bytes=block_bytes,queue_depth=queue_depth,fsync=fsync,borrow_chunks=False,
                store_beta=store_beta,store_variant_df=trait_df is None,
                on_write_progress=(None if on_write_progress is None else
                    lambda event: on_write_progress(replace(event,start=start+event.start,end=start+event.end,
                        variant_df_complete_to=(None if event.variant_df_complete_to is None else start+event.variant_df_complete_to),
                        device=device))))
            if on_writer_open is not None:on_writer_open(writer,device)
            iterator=iter(scan_factory(start,end,device,workers))
            for chunk in iterator:
                if stop.is_set():raise RuntimeError('variant shard cancelled after peer failure')
                # (start, end, beta, t, p, -log10 P, df); without -log10 P
                # (six long) the writer computes it.
                if len(chunk) not in (6,7):raise ValueError('Full variant shards require scan df metadata')
                first,last,beta,t_stat=chunk[:4];variant_df=chunk[-1];logp=chunk[5] if len(chunk)==7 else None
                writer.write_chunk(first,last,beta,t_stat,logp,variant_df=variant_df if trait_df is None else None)
            summary=writer.close()
            return dict(variant_range=[start,end],directory=relative,device=device,
                reader_workers=workers,seconds=time.perf_counter()-began,write=summary)
        except BaseException as error:
            failed(error)
            if writer is not None:writer.abort()
            raise
        finally:
            if iterator is not None and hasattr(iterator,'close'):iterator.close()

    completed=[]
    with ThreadPoolExecutor(max_workers=len(devices),thread_name_prefix='torchgwas-variant') as pool:
        futures=[pool.submit(worker,index) for index in range(len(devices))]
        try:
            for future in as_completed(futures):completed.append(future.result())
        except BaseException as error:
            failed(error)
            if not isinstance(error,Exception):raise
            raise failure[0]
    completed.sort(key=lambda row:row['variant_range'][0])
    if [row['variant_range'] for row in completed]!=[list(span) for span in ranges]:
        raise RuntimeError('Incomplete variant shard coverage')
    if before_publish is not None:before_publish()
    manifest=dict(format='torchgwas-binary-sumstats',version=2,layout='variant_shards',
        shape=[n_variants,len(trait_names)],dtype='float32',byte_order='little',
        arrays=['beta','t_stat','neg_log10_p'] if store_beta else ['t_stat','neg_log10_p'],traits=list(trait_names),
        n_samples=int(n_samples),shards=completed,
        **(dict(df=dict(layout='per_shard',axis='variant'),
                p_value='not stored; two-sided Student t using the matching variant df sidecar')
           if trait_df is None else
           dict(df=trait_df,p_value='not stored; two-sided Student t on t_stat with per-trait df')),
        excluded_convention='NaN marks an excluded variant in stored beta/t arrays',
        genotype_passes=1,reader_workers=reader_workers,
        scope='Disjoint variant ranges, whole phenotype panel per device. No full-matrix gather or cross-device result queue.',
        **(extra_manifest or {}))
    write_manifest(directory,manifest,fsync=fsync)
    return dict(directory=str(directory),cells=n_variants*len(trait_names),
        payload_bytes=sum(row['write']['payload_bytes'] for row in completed),
        df_payload_bytes=sum(row['write'].get('df_payload_bytes',0) for row in completed),
        genotype_passes=1,chunk_size=chunk_size,devices=devices,reader_workers=reader_workers,
        shards=completed,layout='variant_shards')


def _basic_indices(key,shape):
    if key is Ellipsis:key=(slice(None),slice(None))
    elif not isinstance(key,tuple):key=(key,slice(None))
    if len(key)!=2:raise IndexError('Sharded arrays require one or two basic indices')
    ranges=[];scalars=[]
    for value,size in zip(key,shape):
        if value is Ellipsis:value=slice(None)
        scalar=isinstance(value,(int,np.integer)) and not isinstance(value,(bool,np.bool_))
        if scalar:
            value=int(value);value=value+size if value<0 else value
            if not 0<=value<size:raise IndexError('Sharded array index out of range')
            value=slice(value,value+1)
        if not isinstance(value,slice):raise IndexError('Sharded arrays support integer and slice indexing')
        ranges.append(range(*value.indices(size)));scalars.append(scalar)
    return ranges,scalars


def _intersect_rows(rows,start,end):
    """Output-position interval for one shard, using O(1) indexing storage."""
    if rows.step>0:
        first=max(0,-((rows.start-start)//rows.step))
        last=min(len(rows),-((rows.start-end)//rows.step))
    else:
        step=-rows.step
        first=max(0,(rows.start-end)//step+1)
        last=min(len(rows),(rows.start-start)//step+1)
    if first>=last:return None
    local=rows.start+first*rows.step-start
    stop=local+(last-first)*rows.step
    return first,last,slice(local,None if rows.step<0 and stop<0 else stop,rows.step)


class VariantShardedArray:
    """Lazy basic slicing; only the requested result is materialized."""
    def __init__(self,directory,manifest,field):
        self.directory=Path(directory);self.manifest=manifest;self.field=field
        self.shape=(manifest['shape'][0],1 if field=='df' else manifest['shape'][1])
        self.dtype=np.dtype('<f4');self.ndim=2;self.size=self.shape[0]*self.shape[1]

    def __len__(self):return self.shape[0]

    def __array__(self,dtype=None,copy=None):
        if copy is False:raise ValueError('A sharded store cannot expose a contiguous array without a copy')
        value=self[:,:]
        return value.astype(dtype,copy=False) if dtype is not None else value

    def __getitem__(self,key):
        (rows,columns),scalars=_basic_indices(key,self.shape)
        result=np.empty((len(rows),len(columns)),dtype=self.dtype)
        if len(rows) and len(columns):
            column_slice=slice(columns.start,None if columns.step<0 and columns.stop<0 else columns.stop,columns.step)
            from .sumstats import open_binary_sumstats,open_binary_df
            for shard in self.manifest['shards']:
                match=_intersect_rows(rows,*shard['variant_range'])
                if match is None:continue
                first,last,local_rows=match
                path=self.directory/shard['directory']
                if self.field=='df':values=open_binary_df(path)
                else:
                    b,t,logp,_=open_binary_sumstats(path);values={'beta':b,'t_stat':t,'neg_log10_p':logp}[self.field]
                result[first:last]=values[local_rows,column_slice]
        if scalars[0]:result=result[0]
        if scalars[1]:result=result[...,0]
        return result


def _child_path(root,relative):
    path=root/relative
    if Path(relative).is_absolute() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Shard path leaves store directory')
    return path


def open_variant_sharded_sumstats(directory,manifest):
    directory=Path(directory);shape=manifest['shape']
    if (manifest.get('version')!=2 or len(shape)!=2 or min(shape)<1
        or len(manifest['traits'])!=shape[1] or manifest.get('dtype')!='float32'
        or manifest.get('byte_order')!='little'
        or manifest.get('arrays') not in [['t_stat'],['beta','t_stat'],['t_stat','neg_log10_p'],['beta','t_stat','neg_log10_p']]
        or not (manifest.get('df')==dict(layout='per_shard',axis='variant')
                or (isinstance(manifest.get('df'),list) and len(manifest['df'])==shape[1]))):
        raise ValueError('Invalid sharded sumstats manifest')
    per_trait=isinstance(manifest['df'],list)
    offset=0
    from .sumstats import open_binary_df
    for shard in manifest['shards']:
        start,end=shard['variant_range']
        if start!=offset or not start<end<=shape[0]:raise ValueError('Variant shards overlap or leave a gap')
        path=_child_path(directory,shard['directory']);child=read_manifest(path)
        # Children are version 2; per-trait df children written before
        # neglog10p.f32 was stored are version 1.
        if (child.get('version') not in ((1, FORMAT_VERSION) if per_trait else (FORMAT_VERSION,)) or child.get('layout') is not None
                or child.get('format')!='torchgwas-binary-sumstats'):
            raise ValueError('Invalid shard child format')
        if (child['shape']!=[end-start,shape[1]] or child['traits']!=manifest['traits']
            or child['n_samples']!=manifest['n_samples'] or child.get('dtype')!='float32'
            or child.get('byte_order')!='little' or set(child['arrays'])!=set(manifest['arrays'])):
            raise ValueError('Shard metadata differs from root manifest')
        for filename in child['arrays'].values():
            if _child_path(path,filename).stat().st_size!=(end-start)*shape[1]*4:
                raise ValueError('Shard payload length differs from manifest')
        if per_trait:
            if child['df']!=manifest['df']:raise ValueError('Shard per-trait df differs from root manifest')
        else:
            if not isinstance(child['df'],dict):raise ValueError('Variant df sidecar required in every shard')
            _child_path(path,child['df']['array']);open_binary_df(path)
        offset=end
    if offset!=shape[0]:raise ValueError('Incomplete variant coverage')
    return (VariantShardedArray(directory,manifest,'beta') if 'beta' in manifest['arrays'] else None,
            VariantShardedArray(directory,manifest,'t_stat'),
            VariantShardedArray(directory,manifest,'neg_log10_p') if 'neg_log10_p' in manifest['arrays'] else None,
            manifest)
