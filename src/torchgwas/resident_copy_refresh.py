"""Independent resident-copy service refreshed across useful output callbacks.

Only explicitly named CPU copy coefficients are replaced. No stage duration,
allocator rate, GPU/DRAM bandwidth or writer throughput is inferred here.
"""
from copy import deepcopy
from datetime import datetime,timezone
from pathlib import Path
import time

import numpy as np

from .calibration_cache import _age,_digest
from .cpu_service_refresh import CpuServiceRefresh
from .detailed_calibration import read_detailed_profile,sha256_file,write_detailed_profile
from .price_binding import _at,validate_price_bindings


NAME='resident_numpy_copy.v3'
ELEMENTS=8<<20
HOST_BYTES=32<<20  # two arrays, equality temporary, small probe bookkeeping
TARGETS=(('owned_result_copy_scenario','resident_cpu_seconds_per_byte'),
         ('process_units','numpy_copy_bytes'))


def validate_refresh_config(config):
    required={'cache_dir','profile_dir','binding_indexes','max_age_seconds',
              'expected_cpu_seconds','expected_wall_seconds'}
    if not isinstance(config,dict) or set(config)!=required:
        raise ValueError('Resident-copy refresh requires cache/profile directories, binding indexes, lifetime and cost forecasts')
    for field in ('cache_dir','profile_dir'):
        if not isinstance(config[field],str) or not config[field].strip():
            raise ValueError('Nonempty string '+field+' required')
    indexes=config['binding_indexes']
    if (not isinstance(indexes,list) or not 1<=len(indexes)<=128
            or any(type(i) is not int or i<0 for i in indexes) or len(set(indexes))!=len(indexes)):
        raise ValueError('Explicit unique resident-copy binding indexes required')
    for field in ('max_age_seconds','expected_cpu_seconds','expected_wall_seconds'):_age(config[field])


def _set(root,path,value):
    _at(root,path)  # Existing field only, never an implicit new price.
    _at(root,path[:-1])[path[-1]]=value


class ResidentCopyRefresh:
    """Lazy probe with exactly two check samples or seven refresh samples.

    Buffers and cache lookup begin inside the first charged productive step.
    Each advance takes at most two samples. Profile publication is immutable;
    matching cached evidence preserves original record and profile timestamps.
    """
    host_bytes=HOST_BYTES

    def __init__(self,profile,config,*,preserve_artifacts=()):
        validate_refresh_config(config)
        self.original=deepcopy(profile);self.config=deepcopy(config)
        self.indexes=tuple(config['binding_indexes']);self.profile=None
        self.preserve_artifacts=set(preserve_artifacts)
        bindings=profile.get('price_bindings',[])
        if not isinstance(bindings,list) or max(self.indexes)>=len(bindings):
            raise ValueError('Resident-copy refresh binding does not exist')
        self.targets=[]
        for index in self.indexes:
            binding=bindings[index]
            if binding.get('kind')!='cpu_capacity' or not isinstance(binding.get('targets'),list) or not binding['targets']:
                raise ValueError('Refresh requires a declared CPU coefficient')
            for target in binding['targets']:
                path=target.get('context_path')
                if (not isinstance(path,list) or len(path)!=5 or type(path[0]) is not int or path[1]!='profiles'
                        or tuple(path[3:]) not in TARGETS):
                    raise ValueError('Only resident NumPy copy coefficient targets may be refreshed')
                _at(profile['contexts'],path)
                self.targets.append(tuple(path))
        self.protocol=dict(operation='numpy.copyto',dtype='uint8',elements=ELEMENTS,repeats=7,loops=4,
            source='arange',destination='pre-touched',clock='thread_time',
            aggregation='median_repeat_cpu_seconds_per_byte',implementation_sha256=sha256_file(__file__))
        self.dependencies=dict(source_sha256=deepcopy(profile['source_sha256']),
            execution_context=deepcopy(profile['execution_context']),measurement_protocol=self.protocol)
        self.controller=CpuServiceRefresh(config['cache_dir'],NAME,dependencies=self.dependencies,
            work_units=ELEMENTS*4,max_age_seconds=config['max_age_seconds'])
        self.values=self.destination=None;self._closed=False;self._batches=[]
        self._published=None

    @property
    def pending(self):return self.profile is None

    @property
    def deferred_bindings(self):return self.indexes if self.pending else ()

    def sample(self):
        if self.values is None:
            self.values=np.arange(ELEMENTS,dtype=np.uint8)
            self.destination=np.empty_like(self.values)
            np.copyto(self.destination,self.values)
        observed=time.time();wall=time.perf_counter();cpu=time.thread_time()
        for _ in range(4):np.copyto(self.destination,self.values)
        cpu=time.thread_time()-cpu;wall=time.perf_counter()-wall
        if not np.array_equal(self.destination,self.values):raise ValueError('Resident copy probe verification failed')
        return dict(cpu_seconds=cpu,wall_seconds=wall,work_units=ELEMENTS*4,observed_unix_seconds=observed)

    def _bind(self,result):
        updated=deepcopy(self.original);artifact=str(Path(result['path']).resolve(strict=True))
        rate=result['record']['value']['cpu_seconds_per_unit']
        old_artifacts={updated['price_bindings'][index]['artifact'] for index in self.indexes}
        for index in self.indexes:
            targets=[]
            for target in updated['price_bindings'][index]['targets']:
                path=target['context_path'];_set(updated['contexts'],path,rate)
                targets.append(dict(context_path=list(path),value_path=['cpu_seconds_per_unit']))
            updated['price_bindings'][index]=dict(artifact=artifact,kind='cpu_capacity',name=NAME,
                dependencies=deepcopy(self.dependencies),targets=targets,max_age_seconds=self.config['max_age_seconds'])
        required={binding['artifact'] for binding in updated['price_bindings']}|self.preserve_artifacts
        for path in old_artifacts-required:updated['component_artifacts'].pop(path,None)
        updated['component_artifacts'][artifact]=sha256_file(artifact)
        if updated!=self.original:updated['bound_at_utc']=datetime.now(timezone.utc).isoformat()
        checked=validate_price_bindings(updated)
        path=Path(self.config['profile_dir']).expanduser().resolve() / (_digest(updated)+'.json')
        try:write_detailed_profile(updated,path)
        except FileExistsError:
            if read_detailed_profile(path)!=updated:raise ValueError('Immutable refreshed profile changed')
        # Publication itself must not renew or outlive any component evidence.
        validate_price_bindings(updated)
        self._published=dict(path=str(path),profile_sha256=_digest(updated),artifact_sha256=sha256_file(path),
            original_profile_sha256=_digest(self.original),price_evidence=checked)
        self.profile=updated

    def advance(self,*,written_events,issued_revision,validate=None):
        if self._closed:raise ValueError('Resident-copy refresh is closed')
        if not self.pending:return self.profile
        try:
            state=self.controller.advance(self.sample,validate=validate)
            self._batches.append(dict(written_events=written_events,issued_revision=issued_revision,
                state=state['state'],samples=len(state['samples'])))
            if state['state']=='ready':self._bind(state['result']);self.release_buffers()
            elif state['state'] not in ('checking','measuring'):
                raise ValueError('Resident-copy refresh did not produce a stable fresh window: '+state['state'])
            return self.profile
        except BaseException:
            self.release_buffers()
            raise

    def patch_tiles(self,tiles,context_name):
        if self.pending:raise ValueError('Pending measurements cannot price a forecast')
        context_index=next(i for i,c in enumerate(self.profile['contexts']) if c['name']==context_name)
        for path in self.targets:
            if path[0]!=context_index:continue
            for tile in tiles:
                if tile['device']==path[2]:_set(tile['profile'],path[3:],_at(self.profile['contexts'],path))

    def release_buffers(self):self.values=self.destination=None

    def close(self):
        self.release_buffers();self._closed=True

    def snapshot(self):
        return dict(controller=self.controller.snapshot(),batches=deepcopy(self._batches),
            published=deepcopy(self._published),pending=self.pending,closed=self._closed,
            buffers_released=self.values is None and self.destination is None,host_bytes=self.host_bytes,
            targets=[list(path) for path in self.targets],
            scope='Independent fixed-work resident copy CPU coefficient. No loaded-stage or resource-bandwidth estimate; unlisted prices remain unqualified.')
