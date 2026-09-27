"""Independent fixed generic CUDA resource measurements, without GWAS calls."""
import argparse,ctypes,gc,json,os,random,statistics,time
from pathlib import Path
import torch
from torchgwas.detailed_calibration import sha256_file
from torchgwas.geometry_collection import write_record
from direct_device_significance_primitive_measure_20260921 import telemetry

p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--device',default='cuda:0');a=p.parse_args()
root=Path(a.out);root.mkdir(parents=True,exist_ok=False);(root/'harness.py').write_bytes(Path(__file__).read_bytes())
torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.cuda.set_device(a.device)
torch.backends.cuda.matmul.allow_tf32=False
prop=torch.cuda.get_device_properties(a.device);limit=1<<30
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,a.device)
if torch.cuda.mem_get_info(a.device)[0]<limit:raise ValueError('Insufficient free memory')
driver=ctypes.CDLL('libcuda.so.1');assert driver.cuInit(0)==0
l2=ctypes.c_int();assert driver.cuDeviceGetAttribute(ctypes.byref(l2),38,torch.cuda.current_device())==0
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device=a.device,device_uuid=str(prop.uuid),
    name=prop.name,compute_capability=[prop.major,prop.minor],sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
plan=[dict(name='launch_fill',cells=1,nodes=256,work=256,unit='launches'),
    dict(name='l2_add',cells=l2.value//24,nodes=64,work=12*(l2.value//24)*64,unit='logical_bytes'),
    dict(name='hbm_add',cells=1<<24,nodes=16,work=12*(1<<24)*16,unit='logical_bytes'),
    dict(name='fp32_gemm',cells=4096,nodes=4,work=2*4096**3*4,unit='flops')]
write_record(root/'protocol.json',dict(context=context,plan=plan,l2_bytes=l2.value,memory_cap_bytes=limit,
    affinity=sorted(os.sched_getaffinity(0)),harness_sha256=sha256_file(__file__),
    scope='Fixed generic operations and existing PyTorch CUDA graphs. No GWAS shape, selector duration, timing-grid fit, or subtraction of controls.'))
rows=[];states=[]
for item in plan:
    n=item['cells'];name=item['name']
    if name=='fp32_gemm':
        x=torch.ones((n,n),device=a.device);y=torch.ones_like(x);z=torch.empty_like(x)
        def run():torch.mm(x,y,out=z)
    elif name.endswith('_add'):
        x=torch.ones(n,device=a.device);y=torch.ones_like(x);z=torch.empty_like(x)
        def run():torch.add(x,y,out=z)
    else:
        x=torch.zeros(1,device=a.device);y=z=x
        def run():x.fill_(1)
    warm=torch.cuda.Stream();warm.wait_stream(torch.cuda.current_stream())
    with torch.cuda.stream(warm):
        for _ in range(5):run()
    torch.cuda.current_stream().wait_stream(warm);torch.cuda.synchronize()
    graph=torch.cuda.CUDAGraph()
    with torch.cuda.graph(graph):
        for _ in range(item['nodes']):run()
    for _ in range(5):graph.replay()
    torch.cuda.synchronize();states.append(dict(phase=name+'_before',telemetry=telemetry()))
    for repeat in range(9):
        begin=torch.cuda.Event(enable_timing=True);end=torch.cuda.Event(enable_timing=True)
        begin.record();graph.replay();end.record();end.synchronize();seconds=begin.elapsed_time(end)*1e-3
        rows.append(dict(name=name,repeat=repeat,graph_seconds=seconds,work=item['work'],unit=item['unit'],
            rate_per_second=item['work']/seconds))
    states.append(dict(phase=name+'_after',telemetry=telemetry()))
    expected=4096 if name=='fp32_gemm' else 2 if name.endswith('_add') else 1
    assert bool(torch.all(z==expected).item())
    del graph,run,x,y,z,warm;gc.collect();torch.cuda.empty_cache()
def rate(name):return statistics.median(r['rate_per_second'] for r in rows if r['name']==name)
resources=dict(hbm_bytes_per_second=rate('hbm_add'),l2_bytes_per_second=rate('l2_add'),
    fp32_flops_per_second=rate('fp32_gemm'),kernel_launch_seconds=1/rate('launch_fill'),
    host_dispatch_cpu_seconds=0.,available_l2_bytes=l2.value,sm_count=prop.multi_processor_count)
write_record(root/'report.json',dict(context=context,protocol_sha256=sha256_file(root/'protocol.json'),rows=rows,
    resources=resources,telemetry=states,transfer_qualified=False,
    limitations=['Logical attained capacities under recorded load, not hardware counters',
        'FP32 GEMM issue capacity is an explicit proxy for scalar/integer operations',
        'Zero host service is a GPU-only scope and cannot price the full selector']))
print(json.dumps(dict(resources=resources,context=context)),flush=True)
