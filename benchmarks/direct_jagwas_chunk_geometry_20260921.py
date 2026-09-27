"""Untimed statistics/projection launches for one workload across chunk sizes."""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
import torch
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import _gpu_geometry,kernel_census,write_record
from torchgwas.jagwas_projection import JagwasReduction
from torchgwas.reduction_tensor_work import jagwas_tensor_work,jagwas_factor_memory_floor
from torchgwas.tensor_memory import eager_scan_memory

parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
parser.add_argument('--device',default='cuda:0');args=parser.parse_args()
root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
source=source_identity();device=args.device;n,k,c=2049,512,2;chunks=[1,128,256,512]
limit=2<<30
if torch.__version__.split('+')[0]!='2.5.1':raise ValueError('PyTorch 2.5.1 capture required')
torch.set_num_threads(2);torch.set_num_interop_threads(1)
torch.backends.cuda.matmul.allow_tf32=False;torch.cuda.set_device(device)
properties=torch.cuda.get_device_properties(device)
torch.cuda.set_per_process_memory_fraction(limit/properties.total_memory,device)
if torch.cuda.mem_get_info(device)[0]<limit:raise ValueError('Insufficient free memory for bounded census')
requests=[]
for b in chunks:
    ledger=eager_scan_memory(n,b,k,c,2,reduction='jagwas')
    bound=max(ledger['tensor_storage_budget'],jagwas_factor_memory_floor(n,k)['explicit_live_bytes'])
    if bound>limit//2:raise ValueError('Explicit probe request exceeds safety budget')
    requests.append(dict(device=device,shape=[n,b,k,c],device_budget_bytes=limit,source_memory_bound_bytes=bound))
statistics=_gpu_geometry(requests,root)['rows']
write_record(root/'statistics.json',dict(rows=statistics,durations_recorded=False,source_sha256=source))
rng=np.random.default_rng(9219156);y=rng.normal(size=(n,k)).astype(np.float32)
factor=JagwasReduction().prepare(y,device=device);projection=[]
for b in chunks:
    beta=torch.empty((b,k),dtype=torch.float32,device=device)
    stat=torch.randn((b,k),dtype=torch.float32,device=device)
    status=torch.zeros(b,dtype=torch.uint8,device=device)
    df=torch.full((b,),float(n-c-2),dtype=torch.float32,device=device)
    def run():return factor.reduce(beta,stat,status,df,1)
    warm=run();torch.cuda.synchronize(device);del warm
    with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,torch.profiler.ProfilerActivity.CUDA]) as profiler:
        value=run();torch.cuda.synchronize(device)
    del value
    with tempfile.TemporaryDirectory(dir=root,prefix='projection-trace-') as temporary:
        trace=Path(temporary)/'trace.json';profiler.export_chrome_trace(str(trace))
        kernels=kernel_census(json.loads(trace.read_text())['traceEvents'])
    row=dict(N=n,B=b,K=k,compute_dtype='float32',phase='reduce',kernels=kernels,
        tensor_work=jagwas_tensor_work(n,b,k,phase='reduce'),durations_recorded=False)
    projection.append(row);write_record(root/('projection_'+str(b)+'.json'),row)
    print(json.dumps(dict(B=b,statistics_kernels=len(statistics[chunks.index(b)]['kernels']),projection_kernels=len(kernels))),flush=True)
    del beta,stat,status,df,run,profiler;gc.collect();torch.cuda.empty_cache()
if source!=source_identity():raise ValueError('Source changed during capture')
write_record(root/'census.json',dict(source_sha256=source,benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    torch_version=torch.__version__,device=device,gpu=properties.name,compute_capability=[properties.major,properties.minor],
    allow_tf32=False,dimensions=dict(N=n,K=k,C=c),chunks=chunks,statistics=statistics,projection=projection,
    memory_probe_limit_bytes=limit,durations_recorded=False,
    scope='Exact launch geometry for one synthetic phenotype shape and several chunks, including a one-variant tail. No association timings, GPU capacity estimates, CPU prices or ranking validation.'))
