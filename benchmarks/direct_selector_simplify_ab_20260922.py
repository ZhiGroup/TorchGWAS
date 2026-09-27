"""Matched uninstrumented selector simplification; no custom CUDA kernels."""
import argparse,gc,hashlib,inspect,json,random,statistics,time
from pathlib import Path
import numpy as np
import torch
from torchgwas.reduce import device_significant_pairs
from torchgwas.geometry_collection import write_record
from torchgwas.detailed_calibration import sha256_file,source_identity
from direct_device_significance_primitive_measure_20260921 import telemetry
p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
root=Path(a.out);root.mkdir(parents=True,exist_ok=False);(root/'harness.py').write_bytes(Path(__file__).read_bytes())
torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.cuda.set_device(0)
prop=torch.cuda.get_device_properties(0);limit=512<<20
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,0)
if torch.cuda.mem_get_info(0)[0]<limit:raise ValueError('Insufficient free memory')
original=inspect.getsource(device_significant_pairs)
old_df='valid=(status[first:last]==0)&(df>0)&torch.isfinite(df)'
old_values='keep=torch.isfinite(values)&(values.abs()>=limits[:,None])&valid[:,None]'
assert original.count(old_df)==original.count(old_values)==1
proposal=original.replace(old_df,"valid=(status[first:last]==0)&(df>0)&(df<float('inf'))")
proposal=proposal.replace(old_values,"magnitude=values.abs()\n            keep=(magnitude<float('inf'))&(magnitude>=limits[:,None])&valid[:,None]\n            del magnitude")
ns=dict(torch=torch,np=np);exec(proposal,ns);alternative=ns['device_significant_pairs']
(root/'original.py').write_text(original);(root/'proposal.py').write_text(proposal)
source=source_identity();rng=random.Random(9220114);rows=[];states=[];checks=[]
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device_uuid=str(prop.uuid),
    compute_capability=[prop.major,prop.minor],sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
write_record(root/'protocol.json',dict(context=context,source_sha256=source,original_sha256=sha256_file(root/'original.py'),
    proposal_sha256=sha256_file(root/'proposal.py'),repetitions=9,scope='Same process, randomized paired plain selector API calls and final CUDA drain. No profiler is created or imported. Results are consumed and released each call. These are selector component timings, not full GWAS speedups.'))
for b,k,cells in [(13,7,17),(5,37,11),(257,4093,1<<20)]:
    for mode in ['empty','sparse','dense','invalid']:
        ids=np.arange(b*k,dtype=np.int64).reshape(b,k);beta=(ids%1009).astype(np.float32)
        values=np.zeros((b,k),np.float32)
        if mode=='sparse':values.flat[::17]=2.
        elif mode in ('dense','invalid'):values.fill(2.)
        status=np.zeros(b,np.uint8);df=np.full(b,38.,np.float32)
        if mode=='invalid':
            status[1::5]=1;status[2::5]=2;df[::7]=0.;df[3::17]=np.nan;df[4::17]=np.inf;df[5::17]=-np.inf
            values.flat[::13]=np.nan;values.flat[1::19]=np.inf;values.flat[2::19]=-np.inf
        table=np.ones(41,np.float32);table[0]=np.inf
        # Invalid df rows cannot index a CPU array; their effective limit is irrelevant.
        indices=np.zeros(b,np.int64);finite=np.isfinite(df);indices[finite]=np.clip(df[finite],0,40).astype(np.int64)
        ix=np.nonzero(np.isfinite(values)&(np.abs(values)>=table[indices,None])&(status[:,None]==0)&(df[:,None]>0)&finite[:,None])
        expected=[ix[0]+19,ix[1],beta[ix],values[ix],df[ix[0]]]
        tensors=[torch.as_tensor(v,device='cuda') for v in [beta,values,status,df,table]]
        for function in [device_significant_pairs,alternative]:
            output=list(function(*tensors,start=19,max_cells=cells));torch.cuda.synchronize()
            actual=[np.concatenate([part[i] for part in output]) for i in range(2,7)]
            order=np.lexsort((actual[1],actual[0]))
            for got,want in zip(actual,expected):np.testing.assert_array_equal(got[order],want)
            del output,actual,order
        loops=8 if b*k<1000 else 16
        arms=dict(original=device_significant_pairs,proposal=alternative)
        for function in arms.values():
            for _ in range(4):output=list(function(*tensors,start=19,max_cells=cells));del output
        torch.cuda.synchronize();gc.collect();states.append(dict(B=b,K=k,mode=mode,telemetry=telemetry()))
        for repeat in range(9):
            order=list(arms);rng.shuffle(order)
            for arm in order:
                torch.cuda.synchronize();t=time.perf_counter();cpu=time.thread_time()
                for _ in range(loops):output=list(arms[arm](*tensors,start=19,max_cells=cells));del output
                torch.cuda.synchronize()
                rows.append(dict(B=b,K=k,max_cells=cells,mode=mode,repeat=repeat,arm=arm,loops=loops,
                    seconds_per_call=(time.perf_counter()-t)/loops,cpu_seconds_per_call=(time.thread_time()-cpu)/loops))
        checks.append(dict(B=b,K=k,max_cells=cells,mode=mode,independent_arrays_equal=True))
        medians={arm:statistics.median(r['seconds_per_call'] for r in rows if (r['B'],r['K'],r['mode'],r['arm'])==(b,k,mode,arm)) for arm in arms}
        print(json.dumps(dict(B=b,K=k,mode=mode,median_seconds=medians,speedup=medians['original']/medians['proposal'])),flush=True)
        del tensors;gc.collect();torch.cuda.empty_cache()
assert source==source_identity()
write_record(root/'report.json',dict(protocol_sha256=sha256_file(root/'protocol.json'),rows=rows,checks=checks,telemetry=states))
