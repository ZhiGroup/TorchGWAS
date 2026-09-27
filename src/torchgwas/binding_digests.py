"""Optional reuse of content digests under fresh filesystem identity checks.

This avoids rereading stable source/library bytes, not runtime settings. Like
input metadata binding, it trusts filesystem change metadata rather than
authenticating files against adversarial timestamp/inode manipulation.
"""
from copy import deepcopy
import hashlib
from pathlib import Path
import re
import time

from .analytical_plan_cache import input_identity
from .calibration_cache import CalibrationParameterCache,_digest


def _host_identity():
    # Keep device/inode identities local to one host and boot. A remounted or
    # rebooted environment must repopulate instead of trusting reused numbers.
    values=[Path(path).read_text().strip() for path in ('/etc/machine-id','/proc/sys/kernel/random/boot_id')]
    if not all(values):raise ValueError('Host and boot identity required for digest reuse')
    return hashlib.sha256('\n'.join(values).encode()).hexdigest()


def _stable(identity):
    return time.time_ns()-max(identity['mtime_ns'],identity['ctime_ns'])>=1_000_000_000


class BindingDigestCache:
    """Bounded digest groups, immutable on disk and published only on success.

    Every lookup checks resolved path, device/inode, size, mtime and ctime before
    and after use. Recent/future timestamps bypass reuse. Publication rechecks
    all touched files. Runtime hardware, settings and empirical ages never enter
    this cache; their current values still belong to the normal context reader.
    """
    def __init__(self,directory,*,max_groups=8,max_files=512):
        if any(type(v) is not int or v<1 for v in (max_groups,max_files)):
            raise ValueError('Positive binding digest limits required')
        self.cache=CalibrationParameterCache(Path(directory)/'binding-digests-v1')
        self.host_identity=_host_identity();self.max_groups=max_groups;self.max_files=max_files
        self._entries={};self._pending={};self._reference={};self._requested={};self._verification={}
        self._changed=False;self._closed=False
        self._disk_hits=self._memory_hits=self._misses=self._bypasses=self._hashed_files=self._hashed_bytes=0
        self._errors=[]

    def _error(self,error):
        if len(self._errors)<8:self._errors.append(type(error).__name__+': '+str(error))

    def digests(self,paths):
        from .detailed_calibration import sha256_file
        if self._closed:raise ValueError('Binding digest cache is closed')
        paths=[str(Path(path).absolute()) for path in paths]
        if not paths or len(paths)>self.max_files:raise ValueError('Bounded nonempty digest group required')
        identities=[input_identity(path) for path in paths]
        by_path={row['path']:row for row in identities}
        if len(by_path)!=len(paths):raise ValueError('Unique resolved digest paths required')
        if len(set(self._reference)|set(by_path))>self.max_files:
            raise ValueError('Binding digest file count exceeds budget')
        if len(set(self._requested)|set(paths))>self.max_files:
            raise ValueError('Binding digest requested path count exceeds budget')
        for path,row in zip(paths,identities):
            if self._requested.setdefault(path,row)!=row:self._changed=True
        for path,row in by_path.items():
            if self._reference.setdefault(path,row)!=row:self._changed=True
        dependencies=dict(protocol='torchgwas.binding_digests.v1',host_boot_sha256=self.host_identity,
            files=[by_path[path] for path in sorted(by_path)])
        key=_digest(dependencies);unstable=[row['path'] for row in identities if not _stable(row)]
        eligible=not unstable
        values=None
        if eligible and key in self._entries:
            self._memory_hits+=1;values=self._entries[key]
        elif eligible and len(self._entries)<self.max_groups:
            try:
                found=self.cache.lookup('source_work','binding_file_sha256.v1',dependencies=dependencies)
                if found['hit']:
                    candidate=found['record']['value']
                    if (not isinstance(candidate,dict) or set(candidate)!=set(by_path) or
                            any(not isinstance(v,str) or re.fullmatch('[0-9a-f]{64}',v) is None for v in candidate.values())):
                        raise ValueError('Invalid saved binding digests')
                    self._disk_hits+=1;values=candidate
            except (OSError,ValueError,TypeError,KeyError) as error:self._error(error)
        else:
            self._bypasses+=1;eligible=False
        if values is None:
            self._misses+=1;values={path:sha256_file(path) for path in by_path}
            self._hashed_files+=len(paths);self._hashed_bytes+=sum(row['bytes'] for row in identities)
            if eligible:self._pending[key]=(dependencies,deepcopy(values))
        # A path replacement or modification during read/lookup must not create
        # a binding assembled from two different source states.
        if any(input_identity(path)!=row for path,row in zip(paths,identities)):
            self._changed=True;self._pending.pop(key,None)
            raise ValueError('File identity changed during context digest binding')
        # Same-size writes within one filesystem clock tick can preserve stat
        # identity. Recheck bytes at completion for files initially too recent
        # to trust, keeping the first observed digest rather than renewing it.
        for path in unstable:
            if self._verification.setdefault(path,values[path])!=values[path]:self._changed=True
        if eligible:self._entries[key]=deepcopy(values)
        return [values[row['path']] for row in identities]

    def unchanged(self):
        from .detailed_calibration import sha256_file
        if self._closed or self._changed:return False
        try:
            if any(input_identity(path)!=row for path,row in self._requested.items()):return False
            for path,digest in self._verification.items():
                self._hashed_files+=1;self._hashed_bytes+=self._reference[path]['bytes']
                if sha256_file(path)!=digest:return False
            return all(input_identity(path)==row for path,row in self._requested.items())
        except OSError:return False

    def publish(self,*,successful):
        if type(successful) is not bool:raise ValueError('Boolean successful result required')
        if self._closed:raise ValueError('Binding digest cache is closed')
        if not successful or not self.unchanged():
            self._pending.clear()
            return dict(status='unsuccessful' if not successful else 'files_changed',stored=[])
        stored=[]
        for dependencies,values in self._pending.values():
            try:
                stored.append(self.cache.store('source_work','binding_file_sha256.v1',values,
                    dependencies=dependencies,provenance=dict(protocol='sha256-of-stable-file-bytes',
                        scope='Metadata-validated content digests only; no runtime settings or measured rates')))
            except (OSError,ValueError) as error:self._error(error)
        self._pending.clear()
        return dict(status='checked',stored=stored)

    def snapshot(self):
        return dict(disk_hits=self._disk_hits,memory_hits=self._memory_hits,misses=self._misses,
            bypasses=self._bypasses,hashed_files=self._hashed_files,hashed_bytes=self._hashed_bytes,
            groups=len(self._entries),pending=len(self._pending),tracked_files=len(self._reference),
            requested_paths=len(self._requested),content_rechecks=len(self._verification),
            max_groups=self.max_groups,max_files=self.max_files,changed=self._changed,closed=self._closed,
            errors=list(self._errors),scope='Optional stat-bound source/library digests, not cached runtime context. Fresh path/device/inode/size/mtime/ctime checks, one-second timestamp grace, host/boot binding and end-of-run revalidation. Default digest functions without this object still read all bytes.')

    def close(self):
        self._entries.clear();self._pending.clear();self._reference.clear()
        self._requested.clear();self._verification.clear();self._closed=True
