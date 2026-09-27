"""Fixed generic PyTorch integer operators; no custom CUDA code or GWAS input."""
import argparse
import gc
import json
import os
from pathlib import Path
import random
import statistics
import time
import torch
from torchgwas.detailed_calibration import sha256_file,source_identity
from torchgwas.geometry_collection import write_record
from direct_device_significance_primitive_measure_20260921 import telemetry

parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);parser.add_argument('--device',default='cuda:0')
args=parser.parse_args();root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
source=source_identity();(root/'harness.py').write_bytes(Path(__file__).read_bytes())
torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.cuda.set_device(args.device)
prop=torch.cuda.get_device_properties(args.device);limit=256<<20
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,args.device)
if torch.cuda.mem_get_info(args.device)[0]<limit:raise ValueError('Insufficient free memory')
graphs={};repeat_pairs=32
for size in [32,1<<20]:
    values=torch.arange(size,device=args.device,dtype=torch.int64)*17+2027
    quotient=torch.empty_like(values);remainder=torch.empty_like(values)
    def copy_pair():quotient.copy_(values);remainder.copy_(quotient)
    def add_pair():torch.add(values,17,out=quotient);torch.add(quotient,31,out=remainder)
    def divmod_pair():torch.div(values,17,rounding_mode='trunc',out=quotient);torch.fmod(quotient,31,out=remainder)
    for name,function in [('copy_pair',copy_pair),('add_pair',add_pair),('divmod_pair',divmod_pair)]:
        stream=torch.cuda.Stream();stream.wait_stream(torch.cuda.current_stream())
        with torch.cuda.stream(stream):
            for _ in range(4):function()
        torch.cuda.current_stream().wait_stream(stream);torch.cuda.synchronize()
        graph=torch.cuda.CUDAGraph()
        with torch.cuda.graph(graph):
            for _ in range(repeat_pairs):function()
        graph.replay();torch.cuda.synchronize()
        if name=='divmod_pair':
            expected=torch.fmod(torch.div(values,17,rounding_mode='trunc'),31)
            torch.testing.assert_close(remainder,expected,rtol=0,atol=0)
        graphs[(size,name)]=(graph,values,quotient,remainder)
plan=[];rng=random.Random(9211956)
for repeat in range(9):
    order=list(graphs);rng.shuffle(order);plan.append(dict(repeat=repeat,order=order))
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device=args.device,device_uuid=str(prop.uuid),
    name=prop.name,compute_capability=[prop.major,prop.minor],sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
write_record(root/'protocol.json',dict(context=context,source_sha256=source,harness_sha256=sha256_file(__file__),
    affinity=sorted(os.sched_getaffinity(0)),plan=plan,repeat_pairs=repeat_pairs,memory_limit_bytes=limit,
    scope='Fixed int64 quotient then remainder, add and copy pairs. CUDA graph uses existing PyTorch kernels and preallocated arrays. Events surround one replay of 32 pairs; kernel/device launch work remains included. No fitted scan time or subtraction of control durations.'))
observations=[];states=[];gc.collect()
for round_ in plan:
    states.append(dict(repeat=round_['repeat'],phase='before',telemetry=telemetry()))
    for size,name in round_['order']:
        graph=graphs[(size,name)][0]
        start=torch.cuda.Event(enable_timing=True);end=torch.cuda.Event(enable_timing=True)
        torch.cuda.synchronize();start.record();graph.replay();end.record();end.synchronize()
        seconds=start.elapsed_time(end)*1e-3
        observations.append(dict(repeat=round_['repeat'],cells=size,primitive=name,graph_seconds=seconds,
            seconds_per_pair=seconds/repeat_pairs,element_pairs_per_second=size*repeat_pairs/seconds))
    states.append(dict(repeat=round_['repeat'],phase='after',telemetry=telemetry()))
assert source==source_identity()
summaries=[]
for size,name in graphs:
    values=[r for r in observations if (r['cells'],r['primitive'])==(size,name)]
    summaries.append(dict(cells=size,primitive=name,
        seconds_per_pair=statistics.median(r['seconds_per_pair'] for r in values),
        element_pairs_per_second=statistics.median(r['element_pairs_per_second'] for r in values)))
write_record(root/'report.json',dict(protocol_sha256=sha256_file(root/'protocol.json'),observations=observations,summaries=summaries,
    telemetry=states,transfer_qualified=False,
    limitations=['Attained generic operator capacity includes memory and kernel-launch work',
        'PyTorch scalar-division kernels differ from the compiled coordinate-scatter implementation',
        'Cache residence, uniform operands and architecture can alter attainment'],
    scope='Independent resource observations. Transfer to coordinate scatter remains a held-out check, not a published device-selection speed claim.'))
for row in summaries:print(json.dumps(row),flush=True)
