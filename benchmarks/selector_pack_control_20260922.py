"""Matched prototype control: fuse selected gathers and coordinate conversion.

No production replacement or calculator prices are installed. Both paths use
the current production predicate and NumPy flatnonzero. Timings cover complete
selection including output allocation, with returned-array release separate.
"""
import argparse
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import time

from direct_bounded_host_selection_prices_20260921 import measured_call
import numpy as np
import torch
from torchgwas.detailed_calibration import source_identity, _numpy_core_context
from torchgwas.host_significance import ceil_float32, fill_predicate_mask, select_host_pairs
from torchgwas.native_host_predicate import context as predicate_context


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class PackControl:
    def __init__(self, library):
        self.library=ctypes.CDLL(str(library))
        self.call=self.library.selector_pack_control
        self.call.argtypes=[ctypes.c_void_p]*8+[ctypes.c_size_t,ctypes.c_int64,ctypes.c_int64,ctypes.c_ssize_t]
        self.call.restype=None

    def __call__(self,beta,values,row_df,critical,start=0):
        # Deliberately narrow experimental contract. Production fallback for
        # arbitrary strides/dtypes is outside this prototype.
        assert values.ndim==2 and values.dtype==np.float32 and values.flags.c_contiguous
        assert values.flags.aligned
        assert beta is None or (beta.shape==values.shape and beta.dtype==np.float32
            and beta.flags.c_contiguous and beta.flags.aligned)
        df=np.broadcast_to(row_df,values.shape)
        assert df.dtype==np.float32 and df.flags.aligned and df.strides[0]%4==0
        assert values.shape[1]<=1 or df.strides[1]==0
        assert 0<=start<=np.iinfo(np.int64).max-values.shape[0]
        keep=np.empty(values.shape,dtype=bool,order='C')
        limits=np.broadcast_to(ceil_float32(critical),values.shape)
        fill_predicate_mask(values,limits,keep)
        rows=np.flatnonzero(keep).astype(np.int64,copy=False)
        columns=np.empty_like(rows)
        selected_t=np.empty(rows.shape,np.float32)
        selected_df=np.empty(rows.shape,np.float32)
        selected_beta=None if beta is None else np.empty(rows.shape,np.float32)
        self.call(rows.ctypes.data,columns.ctypes.data,values.ctypes.data,
            None if beta is None else beta.ctypes.data,df.ctypes.data,
            selected_t.ctypes.data,None if selected_beta is None else selected_beta.ctypes.data,
            selected_df.ctypes.data,len(rows),values.shape[1],start,df.strides[0]//4)
        return rows,columns,selected_beta,selected_t,selected_df


def reference(beta,values,df,critical,start):
    keep=np.isfinite(values)&(np.abs(values).astype(np.float64)>=critical)
    rows,cols=np.nonzero(keep)
    return (rows.astype(np.int64)+start,cols.astype(np.int64),
        None if beta is None else beta[rows,cols],values[rows,cols],np.broadcast_to(df,values.shape)[rows,cols])


def validate(actual,expected,inputs):
    for got,want in zip(actual,expected):
        if want is None:assert got is None
        else:
            assert got.dtype==want.dtype and got.shape==want.shape
            assert np.array_equal(got,want,equal_nan=True)
            assert all(not np.shares_memory(got,value) for value in inputs)
            assert got.flags.writeable


def main(args):
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    cpp=Path(__file__).with_suffix('.cpp');library=root/'pack.so'
    flags=['-std=c++17','-O3','-fPIC','-shared','-fno-fast-math','-ffp-contract=off']
    command=[os.environ.get('CXX','c++'),*flags,str(cpp),'-o',str(library)]
    subprocess.run(command,check=True)
    compiled=dict(command=command,compiler=subprocess.check_output([command[0],'--version'],text=True).splitlines()[0],
        cpp_sha256=sha(cpp),binary_sha256=sha(library))
    fused=PackControl(library.resolve());os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4)
    assert os.environ['TORCHGWAS_HOST_PREDICATE']=='native'
    assert not np._core.multiarray._get_madvise_hugepage()
    source=source_identity();core=_numpy_core_context();native=predicate_context()
    began=datetime.now(timezone.utc).isoformat();records=[];summaries=[];checks=0
    randomizer=random.Random(9221631)
    # Boundary, empty-shape, scalar/reversed df and immutable-input controls.
    critical=np.array([0.,1e-50,1e-40,1.,np.nextafter(1.,2.),7.123456789,
        np.finfo(np.float32).max,1e40,np.inf,np.nan])[:,None]
    with np.errstate(over='ignore',invalid='ignore'):
        center=critical.astype(np.float32)
        boundary=np.concatenate([center,np.nextafter(center,np.float32(-np.inf)),
            np.nextafter(center,np.float32(np.inf)),-center],axis=1)
    for shape in [(0,7),(3,0),(1,1),(17,31),boundary.shape]:
        values=boundary.copy() if shape==boundary.shape else np.linspace(-7.,7.,np.prod(shape),dtype=np.float32).reshape(shape)
        cutoff=critical if shape==boundary.shape else np.asarray(1.23456789)
        for beta_present in [False,True]:
            beta=np.arange(values.size,dtype=np.float32).reshape(shape) if beta_present else None
            for row_df in [np.asarray(31.,np.float32),(np.arange(shape[0],dtype=np.float32)+20)[::-1,None]]:
                inputs=[values,row_df]+([] if beta is None else [beta])
                for value in inputs:value.flags.writeable=False
                expected=reference(beta,values,row_df,cutoff,19)
                validate(fused(beta,values,row_df,cutoff,19),expected,inputs);checks+=1
    for shape in [(256,4096),(1024,8193)]:
        for density in ['empty','sparse','dense']:
            values=np.full(shape,.5,np.float32)
            if density=='sparse':values.ravel()[::61]=3.;values.ravel()[17::61]=-3.
            elif density=='dense':values.fill(3.);values.ravel()[::2]=-3.
            values.ravel()[3::104729]=np.nan;values.ravel()[7::104729]=np.inf
            beta_values=np.linspace(-2.,2.,values.size,dtype=np.float32).reshape(shape)
            row_df=np.arange(shape[0],dtype=np.float32)[:,None]+500.
            cutoff=np.asarray(1.234567890123,np.float64)
            for beta_present in [False,True]:
                beta=beta_values if beta_present else None
                inputs=[values,row_df]+([] if beta is None else [beta])
                expected=reference(beta,values,row_df,cutoff,8192)
                functions={'current':lambda:select_host_pairs(beta,values,row_df,cutoff,8192),
                    'fused':lambda:fused(beta,values,row_df,cutoff,8192)}
                for function in functions.values():validate(function(),expected,inputs)
                for repeat in range(7):
                    order=list(functions);randomizer.shuffle(order)
                    for name in order:
                        result,observation=measured_call(functions[name])
                        validate(result,expected,inputs);checks+=1
                        wall=time.perf_counter();cpu=time.thread_time();del result
                        observation.update(release_cpu_seconds=time.thread_time()-cpu,
                            release_wall_seconds=time.perf_counter()-wall)
                        records.append(dict(shape=list(shape),density=density,return_beta=beta_present,
                            retained=len(expected[0]),repeat=repeat,implementation=name,observation=observation))
                subset=[r for r in records if r['shape']==list(shape) and r['density']==density and r['return_beta']==beta_present]
                summary=dict(shape=list(shape),density=density,return_beta=beta_present,retained=len(expected[0]),
                    medians={name:{key:statistics.median(r['observation'][key] for r in subset if r['implementation']==name)
                        for key in ['cpu_seconds','wall_seconds','release_cpu_seconds','minor_faults','major_faults']}
                        for name in functions})
                summary['paired_cpu_ratios']=[next(r['observation']['cpu_seconds'] for r in subset if r['repeat']==i and r['implementation']=='fused')/
                    next(r['observation']['cpu_seconds'] for r in subset if r['repeat']==i and r['implementation']=='current') for i in range(7)]
                summaries.append(summary);print(json.dumps(summary),flush=True)
                del expected
    assert source==source_identity() and core==_numpy_core_context() and native==predicate_context()
    report=dict(observation_started_at_utc=began,observation_finished_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=source,numpy_core=core,native_context=native,affinity=sorted(os.sched_getaffinity(0)),
        torch_threads=torch.get_num_threads(),numpy_madvise_hugepage=False,compiled=compiled,
        checks=checks,records=records,summaries=summaries,
        harness_sha256={str(p):sha(p) for p in [Path(__file__),Path(__file__).with_name('direct_bounded_host_selection_prices_20260921.py')]},
        scope=__doc__,prices_installed=False,production_changed=False)
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    main(parser.parse_args())
