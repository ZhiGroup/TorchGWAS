"""Full summary statistics split into contiguous trait tiles.

Each active device writes one bounded, sequential store. A completion manifest
joins the stores without copying or assembling the full variant-by-trait array.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
import threading
import time

import numpy as np

from .sumstats import BinarySumstatsWriter, read_manifest, write_manifest, MANIFEST_NAME


class ScanSourceView:
    """Share read methods without acquiring ownership of the source handles.

    A shallow copy is unsafe: source destructors can close the shared file
    descriptor when a tile finishes. This view owns only per-scan attributes.
    """
    def __init__(self,source):
        self._source=source

    def __getattr__(self,name):
        return getattr(self._source,name)


def write_trait_tiled_sumstats(directory, *, n_variants, trait_names, n_samples,
                               df, trait_block, devices, reader_workers, scan_factory,
                               block_bytes=None, queue_depth=3, fsync=True, store_beta=True,
                               before_publish=None, on_write_progress=None, on_writer_open=None,
                               trait_df=None, extra_manifest=None):
    """Run independent full scans per trait tile, with a shared reader budget.

    trait_df: per-trait df for a panel with missing phenotypes (the
    single-device contract; see sumstats_sharded). Tiles then store their
    slice of it instead of a per-variant df sidecar.

    scan_factory(offset, width, device, workers) returns unreduced chunks with
    a sixth, per-variant df column. Chunks may be borrowed: write_chunk copies
    them into bounded staging before requesting another one. A failed tile
    cancels its peers and prevents publication of the root completion manifest.
    """
    if on_write_progress is not None and not callable(on_write_progress):
        raise ValueError('Dense write progress requires a callback')
    if on_writer_open is not None and not callable(on_writer_open):
        raise ValueError('Dense writer registration requires a callback')
    for name,value in [('trait_block',trait_block),('reader_workers',reader_workers),('queue_depth',queue_depth)]:
        if isinstance(value,bool) or not isinstance(value,int) or value<1:
            raise ValueError(name+' must be a positive integer')
    devices=[str(device) for device in devices]
    if not devices or len(set(devices))!=len(devices):
        raise ValueError('Unique explicit tile devices required')
    if not trait_names:
        raise ValueError('At least one phenotype required')
    if trait_df is not None:
        trait_df=[int(value) for value in trait_df]
        if len(trait_df)!=len(trait_names):raise ValueError('One df per phenotype required')
    tiles=[(start,min(trait_block,len(trait_names)-start)) for start in range(0,len(trait_names),trait_block)]
    devices=devices[:len(tiles)]
    if reader_workers<len(devices):
        raise ValueError('reader_workers must provide at least one reader per active device')
    if block_bytes is not None and block_bytes<=0:
        raise ValueError('block_bytes must be positive')
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    (directory/MANIFEST_NAME).unlink(missing_ok=True)
    stop=threading.Event()
    failure=[];failure_lock=threading.Lock()
    def failed(error):
        # Record before publishing cancellation: a peer may finish aborting
        # before the originating worker, and must not replace its diagnosis.
        with failure_lock:
            if not failure:failure.append(error)
        stop.set()
    per_device,remainder=divmod(reader_workers,len(devices))

    def worker(index,device):
        completed=[]
        workers=per_device+(index<remainder)
        for offset,width in tiles[index::len(devices)]:
            if stop.is_set():
                break
            relative=f'traits_{offset:08d}_{offset+width:08d}'
            writer=BinarySumstatsWriter(directory/relative,n_variants,
                list(trait_names[offset:offset+width]),n_samples,
                df if trait_df is None else trait_df[offset:offset+width],
                block_bytes=block_bytes,queue_depth=queue_depth,fsync=fsync,
                borrow_chunks=False,store_beta=store_beta,store_variant_df=trait_df is None,
                on_write_progress=(None if on_write_progress is None else
                    lambda event,offset=offset,width=width: on_write_progress(replace(event,
                        trait_range=(offset,offset+width),device=device))))
            iterator=None
            started=time.perf_counter()
            try:
                if on_writer_open is not None:on_writer_open(writer,device)
                iterator=iter(scan_factory(offset,width,device,workers))
                for chunk in iterator:
                    if stop.is_set():
                        raise RuntimeError('trait tile cancelled after peer failure')
                    if len(chunk)!=6:
                        raise ValueError('Full trait tiles require scan df metadata')
                    start,end,beta,t_stat,_,variant_df=chunk
                    writer.write_chunk(start,end,beta,t_stat,variant_df if trait_df is None else None)
                summary=writer.close()
                completed.append(dict(trait_range=[offset,offset+width],directory=relative,
                    device=device,reader_workers=workers,seconds=time.perf_counter()-started,write=summary))
            except BaseException as error:
                failed(error)
                writer.abort()
                raise
            finally:
                if iterator is not None and hasattr(iterator,'close'):
                    iterator.close()
        return completed

    completed=[]
    # Waiting for every worker on failure also joins its scan and writer pools.
    with ThreadPoolExecutor(max_workers=len(devices),thread_name_prefix='torchgwas-fulltile') as pool:
        futures=[pool.submit(worker,index,device) for index,device in enumerate(devices)]
        try:
            for future in as_completed(futures):
                completed.extend(future.result())
        except BaseException as error:
            failed(error)
            if not isinstance(error,Exception):raise
            raise failure[0]
    completed.sort(key=lambda row:row['trait_range'][0])
    if [row['trait_range'] for row in completed]!=[[offset,offset+width] for offset,width in tiles]:
        raise RuntimeError('Incomplete trait tile coverage')
    if before_publish is not None:
        before_publish()
    manifest=dict(format='torchgwas-binary-sumstats',version=2,layout='trait_tiles',
        shape=[int(n_variants),len(trait_names)],dtype='float32',byte_order='little',
        arrays=['beta','t_stat'] if store_beta else ['t_stat'],traits=list(trait_names),
        n_samples=int(n_samples),tiles=completed,
        **(dict(df=dict(layout='per_tile'),
                p_value='not stored; two-sided Student t using the matching tile df sidecar')
           if trait_df is None else
           dict(df=trait_df,p_value='not stored; two-sided Student t on t_stat with per-trait df')),
        excluded_convention='NaN marks an excluded variant in stored beta/t arrays',
        genotype_passes=len(tiles),reader_workers=reader_workers,
        scope='Each tile is a sequential variant-major store. No full-matrix assembly.',
        **(extra_manifest or {}))
    write_manifest(directory,manifest,fsync=fsync)
    return dict(directory=str(directory),cells=int(n_variants)*len(trait_names),
        payload_bytes=sum(row['write']['payload_bytes'] for row in completed),
        df_payload_bytes=sum(row['write'].get('df_payload_bytes',0) for row in completed),
        genotype_passes=len(tiles),trait_block=trait_block,devices=devices,
        reader_workers=reader_workers,tiles=completed,layout='trait_tiles')


class TiledSumstatsArray:
    """Lazy two-dimensional basic slicing over a completed trait-tiled store.

    Integer/slice indexing allocates only the requested result. np.asarray is
    an explicit request to materialize the complete array, as for a memmap.
    """
    def __init__(self,directory,manifest,field):
        self.directory=Path(directory)
        self.manifest=manifest
        self.field=field
        self.shape=tuple(manifest['shape'])
        self.dtype=np.dtype('<f4')
        self.ndim=2
        self.size=self.shape[0]*self.shape[1]

    def __len__(self):
        return self.shape[0]

    def __array__(self,dtype=None,copy=None):
        if copy is False:
            raise ValueError('A tiled store cannot expose a contiguous array without a copy')
        value=self[:,:]
        return value.astype(dtype,copy=False) if dtype is not None else value

    def __getitem__(self,key):
        if key is Ellipsis:
            key=(slice(None),slice(None))
        elif not isinstance(key,tuple):
            key=(key,slice(None))
        if len(key)!=2:
            raise IndexError('Tiled arrays require one or two basic indices')
        normalized=[];scalars=[]
        for value,size in zip(key,self.shape):
            if value is Ellipsis:
                value=slice(None)
            scalar=isinstance(value,(int,np.integer)) and not isinstance(value,(bool,np.bool_))
            if scalar:
                value=int(value)
                value=value+size if value<0 else value
                if not 0<=value<size:
                    raise IndexError('Tiled array index out of range')
                value=slice(value,value+1)
            if not isinstance(value,slice):
                raise IndexError('Tiled arrays support integer and slice indexing')
            normalized.append(value);scalars.append(scalar)
        rows,columns=normalized
        row_count=len(range(*rows.indices(self.shape[0])))
        columns=np.arange(*columns.indices(self.shape[1]),dtype=np.int64)
        result=np.empty((row_count,len(columns)),dtype=self.dtype)
        from .sumstats import open_binary_sumstats,open_binary_df
        for tile in self.manifest['tiles']:
            start,end=tile['trait_range']
            positions=np.flatnonzero((columns>=start)&(columns<end))
            if not len(positions) or not row_count:
                continue
            path=self.directory/tile['directory']
            if self.field=='df':
                values=open_binary_df(path)
                local=np.broadcast_to(values,(self.shape[0],end-start))[rows]
            else:
                beta,t_stat,_=open_binary_sumstats(path)
                local=(beta if self.field=='beta' else t_stat)[rows]
            result[:,positions]=local[:,columns[positions]-start]
        if scalars[0]:
            result=result[0]
        if scalars[1]:
            result=result[...,0]
        return result


def open_trait_tiled_sumstats(directory,manifest):
    """Validate tile coverage and child metadata before exposing lazy arrays."""
    directory=Path(directory)
    shape=manifest['shape']
    if manifest.get('version')!=2 or len(shape)!=2 or min(shape)<0 or len(manifest['traits'])!=shape[1]:
        raise ValueError('Invalid tiled sumstats manifest')
    offset=0
    for tile in manifest['tiles']:
        start,end=tile['trait_range']
        if start!=offset or not start<end<=shape[1]:
            raise ValueError('Trait tiles overlap or leave a gap')
        path=directory/tile['directory']
        if Path(tile['directory']).is_absolute() or not path.resolve().is_relative_to(directory.resolve()):
            raise ValueError('Tile path leaves store directory')
        child=read_manifest(path)
        if child.get('layout')=='trait_tiles' or child.get('format')!='torchgwas-binary-sumstats':
            raise ValueError('Invalid tile child format')
        if child['shape']!=[shape[0],end-start] or child['traits']!=manifest['traits'][start:end] or child['n_samples']!=manifest['n_samples']:
            raise ValueError('Tile metadata differs from root manifest')
        if set(child['arrays'])!=set(manifest['arrays']):
            raise ValueError('Tile stored fields differ')
        if isinstance(manifest.get('df'),list) and child['df']!=manifest['df'][start:end]:
            raise ValueError('Tile per-trait df differs from root manifest')
        for filename in child['arrays'].values():
            if (path/filename).stat().st_size!=shape[0]*(end-start)*4:
                raise ValueError('Tile payload length differs from manifest')
        from .sumstats import open_binary_df
        open_binary_df(path)
        offset=end
    if offset!=shape[1]:
        raise ValueError('Incomplete trait coverage')
    return (TiledSumstatsArray(directory,manifest,'beta') if 'beta' in manifest['arrays'] else None,
            TiledSumstatsArray(directory,manifest,'t_stat'),manifest)
