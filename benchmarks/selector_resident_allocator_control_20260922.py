"""Independent resident primitives and held-out selectors with matched kernels.

Default, forwarding-handler and bounded resident-arena scopes are compared.
Resident observations are an experimental component boundary, not automatically
usable production prices. Input/output copies, pool preparation, restoration,
and returned-array destruction are outside the operation clock and recorded
separately where applicable. No GWAS inputs or timings are used.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import time
import numpy as np
import torch
import selector_allocator_control as meter
from direct_bounded_host_selection_prices_20260921 import measured_call
import direct_bounded_host_selection_prices_20260921 as meter_module
from torchgwas.detailed_calibration import source_identity,_numpy_core_context
from torchgwas.host_significance import ceil_float32,fill_predicate_mask,select_host_pairs
from torchgwas.native_host_predicate import context as native_context


def checked_snapshot():
    result=meter.snapshot()
    assert all(result[key]==0 for key in ['invalid_free','foreign_thread','failures','realloc_calls']),result
    return result


def call(function,mode,validate=None):
    meter.begin(mode)
    try:result,timing=measured_call(function)
    finally:meter.restore()
    during=checked_snapshot()
    if validate is not None:validate(result)
    cpu=time.thread_time();wall=time.perf_counter()
    del result
    release_cpu=time.thread_time()-cpu;release_wall=time.perf_counter()-wall
    after=checked_snapshot();assert after['live']==0,after
    timing.update(allocator=during,after_release=after,
        release_cpu_seconds=release_cpu,release_wall_seconds=release_wall,
        residual_cpu_seconds=timing['cpu_seconds']-(during['allocate_cpu_ns']+during['free_cpu_ns'])*1e-9)
    assert timing['residual_cpu_seconds']>=0
    return timing


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4)
    assert np.__version__=='2.2.6' and not np._core.multiarray._get_madvise_hugepage()
    assert np._core.multiarray.get_handler_name()=='default_allocator'
    source=source_identity();core=_numpy_core_context();native=native_context()
    files=[Path(__file__),Path(meter_module.__file__),Path(meter.__file__),Path(__file__).with_name('selector_allocator_control.c')]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    started=datetime.now(timezone.utc).isoformat();meter.configure(512<<20)
    report=dict(observation_started_at_utc=started,source_sha256=source,harness_sha256=hashes,
        context=dict(numpy_core=core,native=native,numpy=np.__version__,torch=torch.__version__,libc=list(platform.libc_ver()),
            affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),
            numpy_madvise_hugepage=False,page_bytes=os.sysconf('SC_PAGE_SIZE'),arena_bytes=512<<20,
            environment={key:os.getenv(key) for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS',
                'OMP_WAIT_POLICY','GOMP_SPINCOUNT','NUMPY_MADVISE_HUGEPAGE','MALLOC_MMAP_THRESHOLD_','MALLOC_TRIM_THRESHOLD_','GLIBC_TUNABLES']}),
        primitives=[],selectors=[],prices_published=0,prediction_complete=False,scope=__doc__)
    def save():(root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    rng=random.Random(9221401)
    b,k=256,4096
    values=np.linspace(-5.,5.,b*k,dtype=np.float32).reshape(b,k)
    limits=np.full((b,1),3.,np.float32);mask=np.empty((b,k),bool)
    masks={name:np.zeros((b,k),bool) for name in ['empty','sparse','dense']}
    masks['sparse'].ravel()[::64]=True;masks['dense'].fill(True)
    canonical=np.arange(b*k,dtype=np.int64);rows=canonical//k;flat_rows=canonical.copy()
    sparse=np.arange(0,b*k,64,dtype=np.int64);zero=canonical[:0]
    row_df=np.full((b,1),39971.,np.float32)
    critical=np.linspace(1.,9.,32768)[:,None]
    def predicate(backend,empty):
        os.environ['TORCHGWAS_HOST_PREDICATE']=backend
        a=values[:0] if empty else values;out=mask[:0] if empty else mask
        fill_predicate_mask(a,np.broadcast_to(limits[:0] if empty else limits,a.shape),out)
        return out
    def coordinates(a,width):
        out=np.empty_like(a);np.divmod(a,width,out=(a,out));return out
    def add(a):a+=8192;return a
    def item(name,pattern,units,function,reset=None):
        return dict(name=name,pattern=pattern,units=units,function=function,reset=reset)
    controls=[
        item('critical_round','empty',0,lambda:ceil_float32(critical[:0])),
        item('critical_round','rows',len(critical),lambda:ceil_float32(critical)),
        item('mask_allocate','fixed',0,lambda:np.empty((b,k),bool)),
        item('coordinate_divmod','empty',0,lambda:coordinates(zero,k)),
        item('coordinate_divmod','power_two',len(rows),lambda:coordinates(flat_rows,4096),lambda:np.copyto(flat_rows,canonical)),
        item('coordinate_divmod','non_power_two',len(rows),lambda:coordinates(flat_rows,4093),lambda:np.copyto(flat_rows,canonical)),
        item('df_gather_row','empty',0,lambda:row_df[:,0][zero]),
        item('df_gather_row','dense',len(rows),lambda:row_df[:,0][rows]),
        item('matrix_gather_flat','empty',0,lambda:values.reshape(-1)[zero]),
        item('matrix_gather_flat','dense',len(rows),lambda:values.reshape(-1)[canonical]),
        item('matrix_gather_flat','stride64',len(sparse),lambda:values.reshape(-1)[sparse]),
        item('inplace_index_add','empty',0,lambda:add(zero)),
        item('inplace_index_add','dense',len(rows),lambda:add(flat_rows),lambda:np.copyto(flat_rows,canonical))]
    for backend in ['numpy','native']:
        for empty in [True,False]:
            controls.append(item('predicate_'+backend,'empty' if empty else 'cells',0 if empty else values.size,
                lambda backend=backend,empty=empty:predicate(backend,empty)))
    for density in ['empty','sparse','dense']:
        for empty in [True,False]:
            controls.append(item('flatnonzero_'+density,'empty' if empty else 'cells',0 if empty else values.size,
                lambda density=density,empty=empty:np.flatnonzero(masks[density][:0] if empty else masks[density]).astype(np.int64,copy=False)))
    def fingerprint(result):
        if result is None:return None
        if isinstance(result,tuple):return [fingerprint(value) for value in result]
        return (result.shape,result.dtype.str,hashlib.sha256(memoryview(result)).hexdigest())
    for control in controls:
        if control['reset']:control['reset']()
        # Uninitialized mask allocation has no value equality contract.
        expected=None if control['name']=='mask_allocate' else fingerprint(control['function']())
        def validate(result):
            if control['name']=='mask_allocate':assert result.shape==(b,k) and result.dtype==np.bool_
            else:assert fingerprint(result)==expected
        for mode in [0,1,2,3]:
            if control['reset']:control['reset']()
            call(control['function'],mode,validate)
        for repeat in range(7):
            modes=[0,1,2,3];rng.shuffle(modes)
            for mode in modes:
                samples=[]
                for _ in range(4 if control['units'] else 8):
                    if control['reset']:control['reset']()
                    samples.append(call(control['function'],mode))
                report['primitives'].append(dict(primitive=control['name'],pattern=control['pattern'],units=control['units'],
                    repeat=repeat,mode=mode,samples=samples,
                    cpu_seconds=statistics.mean(s['cpu_seconds'] for s in samples),
                    residual_cpu_seconds=statistics.mean(s['residual_cpu_seconds'] for s in samples)))
        save();print(json.dumps(dict(primitive=control['name'],pattern=control['pattern'],complete=True)),flush=True)
    # Release reference arrays before the independent held-out panel.
    del controls,control,values,limits,mask,masks,canonical,rows,flat_rows,sparse,zero,row_df,critical
    for shape in [(256,4096),(1024,8193)]:
        tensors=[torch.empty(shape,dtype=torch.float32,pin_memory=True) for _ in range(4)]
        ring=[tensor.numpy() for tensor in tensors]
        df=np.full((shape[0],1),35000.,np.float32);cutoff=np.full((shape[0],1),2.,np.float64)
        for density in ['empty','sparse','dense']:
            for value in ring:
                value.fill(3. if density=='dense' else .25)
                if density=='sparse':value.ravel()[::31]=3.
            os.environ['TORCHGWAS_HOST_PREDICATE']='numpy'
            expected=fingerprint(select_host_pairs(None,ring[0],df,cutoff))
            count=expected[0][0][0]
            def validate_result(result):assert fingerprint(result)==expected
            for backend in ['numpy','native']:
                os.environ['TORCHGWAS_HOST_PREDICATE']=backend
                for mode in [0,1,2,3]:
                    call(lambda:select_host_pairs(None,ring[0],df,cutoff),mode,validate_result)
            for repeat in range(7):
                cases=[(backend,mode) for backend in ['numpy','native'] for mode in [0,1,2,3]];rng.shuffle(cases)
                for backend,mode in cases:
                    os.environ['TORCHGWAS_HOST_PREDICATE']=backend
                    observation=call(lambda:select_host_pairs(None,ring[repeat%4],df,cutoff),mode)
                    report['selectors'].append(dict(shape=shape,density=density,retained=count,backend=backend,
                        mode=mode,repeat=repeat,observation=observation))
            save();print(json.dumps(dict(shape=shape,density=density,complete=True)),flush=True)
        del tensors,ring,value,df,cutoff
    assert source==source_identity() and core==_numpy_core_context() and native==native_context()
    assert hashes=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    assert np._core.multiarray.get_handler_name()=='default_allocator' and meter.snapshot()['live']==0
    report['observation_finished_at_utc']=datetime.now(timezone.utc).isoformat();save()


if __name__=='__main__':main()
