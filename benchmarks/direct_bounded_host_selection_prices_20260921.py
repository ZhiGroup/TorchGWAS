"""Fixed independent primitives for the bounded host selector, no GWAS input.

This records CPU observations and source logical traffic. It does not qualify
concurrent transfer, supply archive/storage prices, or enable automatic tuning.
"""
import argparse, ctypes, hashlib, json, os, platform, random, resource, statistics, time
from datetime import datetime, timezone
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
    NUMPY_MADVISE_HUGEPAGE='0', OMP_WAIT_POLICY='PASSIVE', GOMP_SPINCOUNT='0')
import numpy as np
import torch
from threadpoolctl import threadpool_info
from torchgwas.reduce import SignificantPairs
from torchgwas.host_significance import (ceil_float32, host_selector, NATIVE_HOST_SELECTOR,
    fill_predicate_mask, PREDICATE_MAX_CELLS)
from torchgwas.detailed_calibration import source_identity, _numpy_core_context
from torchgwas.numpy_nonzero_work import nonzero_protocol


class Mallinfo(ctypes.Structure):
    _fields_=[(name,ctypes.c_size_t) for name in
        ['arena','ordblks','smblks','hblks','hblkhd','usmblks','fsmblks','uordblks','fordblks','keepcost']]


_libc=ctypes.CDLL(None)
_libc.mallinfo2.restype=Mallinfo


def allocation_counters():
    value=_libc.mallinfo2()
    return {key:int(getattr(value,key)) for key in ['arena','hblks','hblkhd','uordblks']}


def measured_call(function):
    """Counters bracket the timer; no output destruction inside this boundary.

    mallinfo2 deltas are net live allocator state, not allocation event counts.
    Resource user/system counters are diagnostic and may have coarse resolution.
    """
    allocated_before=allocation_counters();before=resource.getrusage(resource.RUSAGE_THREAD)
    started=time.perf_counter();cpu=time.thread_time()
    result=function()
    cpu_seconds=time.thread_time()-cpu;wall_seconds=time.perf_counter()-started
    after=resource.getrusage(resource.RUSAGE_THREAD);allocated_after=allocation_counters()
    return result,dict(cpu_seconds=cpu_seconds,wall_seconds=wall_seconds,
        user_cpu_seconds=after.ru_utime-before.ru_utime,system_cpu_seconds=after.ru_stime-before.ru_stime,
        minor_faults=after.ru_minflt-before.ru_minflt,major_faults=after.ru_majflt-before.ru_majflt,
        allocator_before=allocated_before,allocator_after=allocated_after)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--out', required=True); args = parser.parse_args()
    root = Path(args.out); root.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, list(range(12,20))); torch.set_num_threads(4)
    assert not np._core.multiarray._get_madvise_hugepage()
    source = source_identity(); harness = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    numpy_core = _numpy_core_context()
    observation_started_at_utc = datetime.now(timezone.utc).isoformat()
    selector = host_selector()
    native = selector == NATIVE_HOST_SELECTOR
    predicate_context = None
    if native:
        from torchgwas.native_host_predicate import context
        predicate_context = context()
    b,k = 256,4096
    values = np.linspace(-5.,5.,b*k,dtype=np.float32).reshape(b,k)
    limits = np.full((b,1),3.,np.float32)
    mask = np.empty((b,k),bool)
    empty_mask = np.zeros((b,k),bool); dense_mask = np.ones_like(empty_mask)
    sparse_mask = empty_mask.copy(); sparse_mask.ravel()[::64] = True
    row_df = np.full((b,1),39971.,np.float32)
    df_bulk = np.full((32768,1),39971.,np.float32); df_empty = df_bulk[:0]
    critical = np.linspace(1.,9.,32768)[:,None]; critical_empty = critical[:0]
    sig = SignificantPairs(.05); sig.prepare_integer_df(40000,k); one = SignificantPairs(1.)
    # Production flatnonzero/divmod coordinates are separate contiguous arrays.
    # np.nonzero's coordinate views instead have interleaved strides and would
    # silently probe a different input layout for row gathers and index copies.
    canonical = np.arange(b*k,dtype=np.int64); rows = canonical//k
    zero = canonical[:0]; flat_rows = canonical.copy()
    sparse_flat = np.arange(0,b*k,64,dtype=np.int64)
    def predicate(a,lim,out):
        fill_predicate_mask(a, np.broadcast_to(lim,a.shape), out)
        return out
    def coordinates(a,width):
        result=np.empty_like(a); np.divmod(a,width,out=(a,result)); return result
    def add_inplace(a):
        a += 8192; return a
    def control(units,function,pattern,reset=None):
        return dict(units=units,function=function,pattern=pattern,reset=reset)
    controls={
        'critical_one': (29,[control(0,lambda:one.critical_abs_t(df_empty,k),'empty'),control(len(df_bulk),lambda:one.critical_abs_t(df_bulk,k),'integer')]),
        'critical_lookup': (138,[control(0,lambda:sig.critical_abs_t(df_empty,k),'empty'),control(len(df_bulk),lambda:sig.critical_abs_t(df_bulk,k),'integer')]),
        'critical_round': (50,[control(0,lambda:ceil_float32(critical_empty),'empty'),control(len(critical),lambda:ceil_float32(critical),'cutoff')]),
        'mask_allocate': (0,[control(0,lambda:np.empty((b,k),bool),'fixed_1M_uninitialized_allocation')]),
        ('predicate_native' if native else 'predicate_block'): (5 if native else 25,
            [control(0,lambda:predicate(values[:0],limits[:0],mask[:0]),'empty'),
             control(values.size,lambda:predicate(values,limits,mask),'fixed_1M_block')]),
        'flatnonzero_empty': (1,[control(0,lambda:np.flatnonzero(empty_mask[:0]).astype(np.int64,copy=False),'empty'),control(values.size,lambda:np.flatnonzero(empty_mask).astype(np.int64,copy=False),'all_false')]),
        'flatnonzero_sparse': (2,[control(0,lambda:np.flatnonzero(empty_mask[:0]).astype(np.int64,copy=False),'empty'),control(values.size,lambda:np.flatnonzero(sparse_mask).astype(np.int64,copy=False),'stride64')]),
        'flatnonzero_dense': (2,[control(0,lambda:np.flatnonzero(empty_mask[:0]).astype(np.int64,copy=False),'empty'),control(values.size,lambda:np.flatnonzero(dense_mask).astype(np.int64,copy=False),'all_true')]),
        'coordinate_divmod': (24,[control(0,lambda:coordinates(zero,k),'empty'),control(len(rows),lambda:coordinates(flat_rows,4096),'power_two',lambda:np.copyto(flat_rows,canonical)),control(len(rows),lambda:coordinates(flat_rows,4093),'non_power_two',lambda:np.copyto(flat_rows,canonical))]),
        'df_gather_row': (16,[control(0,lambda:np.broadcast_to(row_df,values.shape)[:,0][zero],'empty'),control(len(rows),lambda:np.broadcast_to(row_df,values.shape)[:,0][rows],'dense')]),
        'matrix_gather_flat': (16,[control(0,lambda:values.reshape(-1)[zero],'empty'),control(len(rows),lambda:values.reshape(-1)[canonical],'dense'),control(len(sparse_flat),lambda:values.reshape(-1)[sparse_flat],'stride64')]),
        'inplace_index_add': (16,[control(0,lambda:add_inplace(zero),'empty'),control(len(rows),lambda:add_inplace(flat_rows),'dense')]),
        'index_cast': (16,[control(0,lambda:zero.astype(np.int64),'empty'),control(len(rows),lambda:rows.astype(np.int64),'dense')]),
        'index_add': (16,[control(0,lambda:zero+8192,'empty'),control(len(rows),lambda:rows+8192,'dense')])}
    records=[]; rng=random.Random(9219179)
    for name,(traffic,variants) in controls.items():
        for item in variants:
            if item['reset']:item['reset']()
            result=item['function']();del result
        for repeat in range(9):
            order=list(range(len(variants)));rng.shuffle(order)
            for index in order:
                item=variants[index];samples=[]
                for call in range(4 if item['units'] else 32):
                    if item['reset']:item['reset']()
                    result,sample=measured_call(item['function'])
                    # Keep destruction separate from the historical price boundary.
                    started=time.perf_counter();cpu=time.thread_time()
                    del result
                    sample['release_cpu_seconds']=time.thread_time()-cpu
                    sample['release_wall_seconds']=time.perf_counter()-started
                    sample['allocator_released']=allocation_counters()
                    samples.append(sample)
                records.append(dict(primitive=name,control=index,repeat=repeat,units=item['units'],pattern=item['pattern'],
                    samples=samples,cpu_seconds=statistics.mean(r['cpu_seconds'] for r in samples),
                    wall_seconds=statistics.mean(r['wall_seconds'] for r in samples)))
        (root/'observations.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps(dict(primitive=name,complete=True)),flush=True)
    prices={};qualification=[]
    for name,(traffic,variants) in controls.items():
        medians=[statistics.median(r['cpu_seconds'] for r in records if r['primitive']==name and r['control']==index) for index in range(len(variants))]
        fixed=medians[0]
        rates=[(value-fixed)/item['units'] for value,item in zip(medians[1:],variants[1:])]
        if any(value<0 for value in rates):qualification.append(name+' bulk below fixed call: '+str(rates))
        prices[name]=dict(call_cpu_seconds=fixed,unit_cpu_seconds=max(rates) if rates else 0.,dram_bytes_per_unit=traffic)
    assert source == source_identity(), 'Source changed during probe'
    assert harness == hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    assert _numpy_core_context() == numpy_core
    if native: assert context() == predicate_context
    result=dict(host_selector=selector,predicate_max_cells=PREDICATE_MAX_CELLS,nonzero_protocol=nonzero_protocol(),prices=prices,
        observation_started_at_utc=observation_started_at_utc,
        observation_finished_at_utc=datetime.now(timezone.utc).isoformat(),
        predicate_context=predicate_context,
        observations=records,source_sha256=source,harness_sha256=harness,qualification_errors=qualification,
        input_layouts={name:dict(shape=list(value.shape),strides=list(value.strides),dtype=value.dtype.str)
            for name,value in [('values',values),('row_indices',rows),('flat_indices',canonical),('sparse_flat_indices',sparse_flat)]},
        primitive_differences_nonnegative=not qualification,concurrent_transfer_qualified=False,prediction_complete=False,
        context=dict(python=platform.python_version(),numpy=np.__version__,numpy_core=numpy_core,torch=torch.__version__,
            affinity=sorted(os.sched_getaffinity(0)),torch_threads=torch.get_num_threads(),
            torch_interop_threads=torch.get_num_interop_threads(),numpy_madvise_hugepage=bool(np._core.multiarray._get_madvise_hugepage()),
            environment={key:os.getenv(key) for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMPY_MADVISE_HUGEPAGE','OMP_WAIT_POLICY','GOMP_SPINCOUNT']},
            numerical_pools=threadpool_info()),
        allocation_observation_protocol='thread_faults_net_mallinfo2_separate_release_v1',
        scope='Fixed generic controls; maximum declared bulk-minus-empty CPU service per source unit, without interpolation, GWAS inputs or association timings. Timer overhead remains in fixed call costs. Output destruction and input restoration excluded from prices; destruction is observed separately. Mask allocation is a fixed 1M-cell reference. Logical traffic is source-derived. Per-call thread fault counters and net live allocator deltas expose reference memory state, but do not qualify transfer to different allocator history or concurrency. Archive/queue/storage prices must be supplied separately.')
    (root/'selection_prices.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if qualification:raise ValueError('; '.join(qualification))


if __name__ == '__main__':main()
