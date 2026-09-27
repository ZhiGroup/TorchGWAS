"""Independent gather controls and held-out selector comparison, no GWAS fit."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import statistics
import time
import numpy as np
import torch
from torchgwas.host_significance import ceil_float32,fill_predicate_mask
from torchgwas.detailed_calibration import source_identity,_numpy_core_context


def legacy_select(beta,values,row_df,critical,start=0):
    # Exact production bounded-flat-v1 selection order, frozen for comparison.
    if values.dtype==np.float32:critical=ceil_float32(critical)
    keep=np.empty(values.shape,dtype=bool,order='C')
    limits=np.broadcast_to(critical,values.shape);fill_predicate_mask(values,limits,keep)
    rows=np.flatnonzero(keep).astype(np.int64,copy=False)
    columns=np.empty_like(rows)
    if values.shape[1]:np.divmod(rows,values.shape[1],out=(rows,columns))
    selected_df=np.broadcast_to(row_df,values.shape)[rows,columns]
    selected_beta=None if beta is None else beta[rows,columns]
    selected_t=values[rows,columns]
    rows+=start
    return rows,columns,selected_beta,selected_t,selected_df


def candidate_select(beta,values,row_df,critical,start=0):
    if values.dtype==np.float32:critical=ceil_float32(critical)
    keep=np.empty(values.shape,dtype=bool,order='C')
    limits=np.broadcast_to(critical,values.shape);fill_predicate_mask(values,limits,keep)
    rows=np.flatnonzero(keep).astype(np.int64,copy=False)
    flat_t=values.flags.c_contiguous
    flat_beta=beta is not None and beta.shape==values.shape and beta.flags.c_contiguous
    selected_beta=beta.reshape(-1)[rows] if flat_beta else None
    selected_t=values.reshape(-1)[rows] if flat_t else None
    columns=np.empty_like(rows)
    if values.shape[1]:np.divmod(rows,values.shape[1],out=(rows,columns))
    df=np.broadcast_to(row_df,values.shape)
    if not values.shape[1]:selected_df=np.empty(rows.shape,dtype=df.dtype)
    elif df.strides[1]==0:selected_df=df[:,0][rows]
    else:selected_df=df[rows,columns]
    if beta is not None and not flat_beta:selected_beta=beta[rows,columns]
    if not flat_t:selected_t=values[rows,columns]
    rows+=start
    return rows,columns,selected_beta,selected_t,selected_df


def exact(left,right):
    for a,b in zip(left,right):
        if a is None or b is None:assert a is b;continue
        assert a.shape==b.shape and a.dtype==b.dtype
        assert hashlib.sha256(memoryview(a).cast('B')).digest()==hashlib.sha256(memoryview(b).cast('B')).digest()


def measure(function):
    before=resource.getrusage(resource.RUSAGE_THREAD)
    start=time.perf_counter();cpu=time.thread_time()
    result=function()
    cpu=time.thread_time()-cpu;wall=time.perf_counter()-start
    after=resource.getrusage(resource.RUSAGE_THREAD)
    # Destruction stays outside the primitive/selector boundary in both arms.
    del result
    return dict(cpu_seconds=cpu,wall_seconds=wall,user_seconds=after.ru_utime-before.ru_utime,
        system_seconds=after.ru_stime-before.ru_stime,
        minor_faults=after.ru_minflt-before.ru_minflt,major_faults=after.ru_majflt-before.ru_majflt)


def main(args):
    directory=Path(args.out);directory.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(4);torch.set_num_interop_threads(1)
    assert not np._core.multiarray._get_madvise_hugepage()
    source=source_identity();started=datetime.now(timezone.utc).isoformat()
    rng=random.Random(922241);observations=[];selector=[]
    # Generic primitive controls, not a candidate layout grid. The 1M-cell
    # source and stride-64 survivor pattern match the independent collector.
    for shape in [(256,4096),(1024,8193)]:
        values=np.linspace(-3,3,np.prod(shape),dtype=np.float32).reshape(shape)
        df=np.arange(shape[0],dtype=np.float32)[:,None]+30000
        for pattern,flat in [('empty',np.empty(0,np.int64)),
                ('stride64',np.arange(0,values.size,64,dtype=np.int64)),
                ('dense',np.arange(values.size,dtype=np.int64))]:
            rows,columns=np.divmod(flat,shape[1])
            calls=dict(df_2d=lambda:np.broadcast_to(df,shape)[rows,columns],
                df_row=lambda:np.broadcast_to(df,shape)[:,0][rows],
                matrix_2d=lambda:values[rows,columns],matrix_flat=lambda:values.reshape(-1)[flat])
            exact([calls['df_2d'](),calls['matrix_2d']()],[calls['df_row'](),calls['matrix_flat']()])
            for function in calls.values():value=function();del value
            for repeat in range(7):
                order=list(calls);rng.shuffle(order)
                for name in order:
                    observations.append(dict(shape=list(shape),pattern=pattern,primitive=name,
                        retained=len(flat),repeat=repeat,**measure(calls[name])))
        del rows,columns,flat,calls,values,df
    (directory/'primitive_observations.json').write_text(json.dumps(observations,indent=2)+'\n')
    print('Independent primitives complete',flush=True)
    # Different survivor pattern in held-out whole-selector observations.
    for shape in [(256,4096),(1024,8193)]:
        tensors=[torch.empty(shape,dtype=torch.float32,pin_memory=True) for _ in range(4)]
        ring=[tensor.numpy() for tensor in tensors]
        for values in ring:
            values.fill(.25);values.ravel()[::31]=3.;values.ravel()[::71]=np.nan
            values.ravel()[1::127]=np.inf;values.ravel()[2::127]=-np.inf
        df=np.arange(shape[0],dtype=np.float32)[:,None]+30000
        for label,limit in [('sparse',2.),('dense',0.)]:
            critical=np.full((shape[0],1),limit,np.float64)
            for store_beta in [False,True]:
                beta=ring[1] if store_beta else None
                left=legacy_select(beta,ring[0],df,critical,17)
                right=candidate_select(beta,ring[0],df,critical,17)
                exact(left,right);retained=len(left[0]);del left,right
                for repeat in range(6):
                    order=['legacy','candidate'];rng.shuffle(order)
                    values=ring[repeat%4];beta=ring[(repeat+1)%4] if store_beta else None
                    for name in order:
                        function=legacy_select if name=='legacy' else candidate_select
                        selector.append(dict(shape=list(shape),occupancy=label,store_beta=store_beta,
                            retained=retained,implementation=name,repeat=repeat,
                            **measure(lambda:function(beta,values,df,critical,17))))
                print(json.dumps(dict(shape=shape,occupancy=label,store_beta=store_beta,
                    median_cpu_seconds={name:statistics.median(r['cpu_seconds'] for r in selector
                        if r['shape']==list(shape) and r['occupancy']==label and r['store_beta']==store_beta and r['implementation']==name)
                        for name in ['legacy','candidate']})),flush=True)
                (directory/'selector_observations.json').write_text(json.dumps(selector,indent=2)+'\n')
        del tensors,ring,values,beta,df,critical
    assert source_identity()==source
    report=dict(source_sha256=source,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        observation_started_at_utc=started,observation_finished_at_utc=datetime.now(timezone.utc).isoformat(),
        numpy_core=_numpy_core_context(),numpy_version=np.__version__,torch_version=torch.__version__,
        settings=dict(affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),
            torch_interop_threads=torch.get_num_interop_threads(),environment={k:os.getenv(k) for k in
            ['NUMPY_MADVISE_HUGEPAGE','OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS',
             'OMP_WAIT_POLICY','GOMP_SPINCOUNT','TORCHGWAS_HOST_PREDICATE']}),
        primitive_observations=observations,selector_observations=selector,outputs_exact=True,scope=__doc__)
    (directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    main(parser.parse_args())
