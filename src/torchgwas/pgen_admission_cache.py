"""Immutable, source-bound structural PGEN admission vectors.

Only exact locator-derived memory envelopes and non-LD base positions are
cached. This module never records timings, prices, or live capacity.
"""
import hashlib
import json
import os
from pathlib import Path
import struct
import tempfile
import time

import numpy as np

from .analytical_plan_cache import canonical, input_identity, input_is_stable
from .pgen_memory_layout import COMPACT_KIND

FORMAT = 'torchgwas-pgen-admission-v2'
MAX_ENTRY_BYTES = 256 << 20
MAX_ENTRIES = 8
MAX_MANIFEST_BYTES = 8192
_ARRAYS = (('fine_payload_bytes', '<u8'),
           ('fine_prefix_bytes', '<u8'),
           ('fine_extra_workspace_bytes', '<u8'),
           ('payload_cumulative_bytes', '<u8'),
           ('bases', '<u4'))


class PgenAdmissionCache:
    """Optional exact structural reuse; misses preserve ordinary admission."""

    def __init__(self, directory, identity, header, fine):
        self.directory = Path(directory).expanduser() / 'torchgwas-pgen-admission-v2'
        self.identity = identity
        self.header = header
        self.fine = fine
        source = Path(__file__).with_name('pgen_memory_layout.py')
        request = dict(format=FORMAT, input=identity, samples=header.sample_ct,
                       markers=header.variant_ct, fine=fine,
                       implementation_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
        self.key = hashlib.sha256(canonical(request)).hexdigest()
        self.path = self.directory / (self.key + '.bin')
        self.status = 'miss'
        self.write_status = 'not_attempted'
        self._pending = None
        self.retained_blob_bytes = 0
        self.publication_wall_seconds = 0.
        self.publication_cpu_seconds = 0.

    def _layout(self, arrays, payload):
        m = int(self.header.variant_ct)
        file_bytes = self.identity['bytes']
        return dict(kind=COMPACT_KIND, path=self.identity['path'],
                    samples=int(self.header.sample_ct), file_markers=m,
                    file_bytes=file_bytes, index_bytes=file_bytes-payload,
                    chunk_markers=self.fine, markers=m, variant_range=[0,m],
                    record_payload_bytes=payload,
                    fine_payload_bytes=arrays['fine_payload_bytes'],
                    fine_prefix_bytes=arrays['fine_prefix_bytes'],
                    fine_extra_workspace_bytes=arrays['fine_extra_workspace_bytes'],
                    payload_cumulative_bytes=arrays['payload_cumulative_bytes'],
                    logical_chunks=(m+self.fine-1)//self.fine,
                    scope='Full-file validated header vectors for exact fixed and shifted reader-memory maxima; no per-chunk Python records or genotype payload reads.')

    def _valid(self, arrays, payload, *, full=True, bases_required=True):
        count=(self.header.variant_ct+self.fine-1)//self.fine
        if (payload<1 or payload>self.identity['bytes'] or
            any(arrays[name].dtype != np.dtype(dtype) or arrays[name].ndim != 1 or
                len(arrays[name]) != (count+1 if name=='payload_cumulative_bytes' else count)
                for name,dtype in _ARRAYS[:-1])):
            return False
        if bases_required:
            bases=arrays['bases']
            if (bases.dtype != np.dtype('<u4') or bases.ndim != 1 or
                len(bases)<1 or len(bases)>self.header.variant_ct or
                int(bases[0])!=0 or int(bases[-1])>=self.header.variant_ct or
                (full and np.any(bases[1:]<=bases[:-1]))):
                return False
        cumulative=arrays['payload_cumulative_bytes']
        if (int(cumulative[0])!=0 or int(cumulative[-1])!=payload or
            (full and (np.any(cumulative[1:]-cumulative[:-1]!=arrays['fine_payload_bytes']) or
                       np.any(arrays['fine_payload_bytes']==0)))):
            return False
        return True

    def load(self):
        if not input_is_stable(self.identity):
            self.status='unstable_source'
            return None
        try:
            with self.path.open('rb') as stream:
                size=os.fstat(stream.fileno()).st_size
                if size>MAX_ENTRY_BYTES:
                    self.status='oversized';return None
                raw=stream.read(MAX_ENTRY_BYTES+1)
            if len(raw)>MAX_ENTRY_BYTES or len(raw)!=size or len(raw)<8:
                self.status='invalid';return None
            length=struct.unpack_from('<Q',raw)[0]
            if length>MAX_MANIFEST_BYTES or 8+length>len(raw):
                self.status='invalid';return None
            manifest=json.loads(raw[8:8+length])
            if (set(manifest)!={'format','key','payload_bytes','base_count','sha256'} or
                manifest['format']!=FORMAT or manifest['key']!=self.key or
                type(manifest['payload_bytes']) is not int or
                type(manifest['base_count']) is not int):
                self.status='invalid';return None
            count=(self.header.variant_ct+self.fine-1)//self.fine
            expected=(3*count+count+1)*8+manifest['base_count']*4
            body_start=8+length
            padding=(-body_start)%8
            if raw[body_start:body_start+padding]!=bytes(padding):
                self.status='invalid';return None
            body=memoryview(raw)[body_start+padding:]
            if (len(body)!=expected or
                hashlib.sha256(body).hexdigest()!=manifest['sha256']):
                self.status='invalid';return None
            arrays={};cursor=0
            for name,dtype in _ARRAYS:
                entries=manifest['base_count'] if name=='bases' else count+int(name=='payload_cumulative_bytes')
                amount=entries*np.dtype(dtype).itemsize
                arrays[name]=np.frombuffer(body,dtype=dtype,count=entries,offset=cursor)
                arrays[name].flags.writeable=False
                if not arrays[name].flags.aligned:
                    self.status='invalid';return None
                cursor+=amount
            if not self._valid(arrays,manifest['payload_bytes'],full=False) or input_identity(self.identity['path'])!=self.identity:
                self.status='invalid';return None
            try:os.utime(self.path,None)
            except OSError:pass
            self.status='hit'
            self.retained_blob_bytes=len(raw)
            return self._layout(arrays,manifest['payload_bytes']),arrays['bases']
        except FileNotFoundError:
            self.status='miss'
        except (OSError,ValueError,TypeError,KeyError,OverflowError,MemoryError):
            self.status='invalid'
        return None

    def stage(self, layout, bases):
        if self.status=='hit':
            return
        if (not input_is_stable(self.identity) or
            input_identity(self.identity['path'])!=self.identity or
            self.header.variant_ct>np.iinfo(np.uint32).max):
            self.write_status='unstable_source'
            return
        arrays={name:layout[name] for name,_ in _ARRAYS[:-1]}
        arrays['bases']=None if bases is None else np.asarray(bases,dtype='<u4')
        if not self._valid(arrays,layout['record_payload_bytes'],bases_required=bases is not None):
            self.write_status='invalid_source'
            return
        # A cold job keeps only fine-grid vectors until successful finish.
        # A cache hit already carries the base index for its first window.
        self._pending=(dict(arrays, bases=None),arrays['bases'],
                       int(layout['record_payload_bytes']))
        self.write_status='pending'

    def owned_extra_bytes(self, retained_bases_bytes):
        if self.status=='hit':return max(0,self.retained_blob_bytes-retained_bases_bytes)
        return self.pending_bytes

    @property
    def pending_bytes(self):
        if self._pending is None:return 0
        return sum(a.nbytes for a in self._pending[0].values() if a is not None)

    def publish(self, *, successful):
        if self._pending is None or not successful:
            if self._pending is not None:self.write_status='unsuccessful'
            self._pending=None
            return self.write_status
        temporary=None
        began=time.perf_counter();cpu=time.process_time()
        try:
            if input_identity(self.identity['path'])!=self.identity:
                self.write_status='source_changed';return self.write_status
            vectors,bases,payload=self._pending
            if bases is None:
                if self.header.variant_ct>np.iinfo(np.uint32).max:
                    self.write_status='oversized';return self.write_status
                bases=np.flatnonzero((self.header.vrtypes!=2)&(self.header.vrtypes!=3)).astype('<u4')
            arrays=dict(vectors,bases=bases)
            if not self._valid(arrays,payload):
                self.write_status='invalid_source';return self.write_status
            body=b''.join(np.asarray(arrays[name],dtype=dtype).tobytes(order='C')
                          for name,dtype in _ARRAYS)
            manifest=canonical(dict(format=FORMAT,key=self.key,payload_bytes=payload,
                                    base_count=len(bases),
                                    sha256=hashlib.sha256(body).hexdigest()))
            header=struct.pack('<Q',len(manifest))+manifest
            record=header+bytes((-len(header))%8)+body
            if len(manifest)>MAX_MANIFEST_BYTES or len(record)>MAX_ENTRY_BYTES:
                self.write_status='oversized';return self.write_status
            self.directory.mkdir(parents=True,exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=self.directory,prefix='.pending-',suffix='.bin',delete=False) as stream:
                temporary=Path(stream.name)
                stream.write(record);stream.flush();os.fsync(stream.fileno())
            os.replace(temporary,self.path);temporary=None
            entries=[]
            for path in self.directory.glob('*.bin'):
                if len(path.stem)==64:
                    try:entries.append((path.stat().st_mtime_ns,path))
                    except FileNotFoundError:pass
            for _,path in sorted(entries)[:max(0,len(entries)-MAX_ENTRIES)]:
                if path!=self.path:
                    try:path.unlink(missing_ok=True)
                    except OSError:pass
            self.write_status='stored'
        except (OSError,ValueError,TypeError,MemoryError):
            self.write_status='unavailable'
        finally:
            self.publication_wall_seconds=time.perf_counter()-began
            self.publication_cpu_seconds=time.process_time()-cpu
            self._pending=None
            if temporary is not None:
                try:temporary.unlink(missing_ok=True)
                except OSError:pass
        return self.write_status

    def audit(self):
        return dict(status=self.status,write_status=self.write_status,key=self.key,
                    pending_bytes=self.pending_bytes,retained_blob_bytes=self.retained_blob_bytes,
                    publication_wall_seconds=self.publication_wall_seconds,
                    publication_cpu_seconds=self.publication_cpu_seconds,
                    format=FORMAT,
                    scope='Source-bound exact structural admission only; no timing prices or live capacity.')
