"""Paired profiler controls for generic eager and graph GPU execution."""
import argparse,gc,json,random,statistics,tempfile,time
from pathlib import Path
import torch
from torchgwas.detailed_calibration import sha256_file
from torchgwas.geometry_collection import write_record
from direct_device_significance_primitive_measure_20260921 import telemetry
p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
root=Path(a.out);root.mkdir(parents=True,exist_ok=False);(root/'harness.py').write_bytes(Path(__file__).read_bytes())
torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.cuda.set_device(0)
prop=torch.cuda.get_device_properties(0);limit=256<<20
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,0)
if torch.cuda.mem_get_info(0)[0]<limit:raise ValueError('Insufficient free memory')
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device_uuid=str(prop.uuid),
    compute_capability=[prop.major,prop.minor],sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
rows=[];states=[];rng=random.Random(9220110)
for name,cells,nodes in [('fill',1,256),('add',1<<20,32)]:
    x=torch.ones(cells,device='cuda');y=torch.ones_like(x);z=torch.empty_like(x)
    def op():
        if name=='fill':x.fill_(1)
        else:torch.add(x,y,out=z)
    stream=torch.cuda.Stream();stream.wait_stream(torch.cuda.current_stream())
    with torch.cuda.stream(stream):
        for _ in range(5):op()
    torch.cuda.current_stream().wait_stream(stream);torch.cuda.synchronize()
    graph=torch.cuda.CUDAGraph()
    with torch.cuda.graph(graph):
        for _ in range(nodes):op()
    for _ in range(5):graph.replay()
    torch.cuda.synchronize()
    for repeat in range(7):
        order=[(mode,measure) for mode in ['graph','eager'] for measure in ['plain','profile']];rng.shuffle(order)
        states.append(dict(name=name,repeat=repeat,telemetry=telemetry()))
        for mode,measure in order:
            gc.collect();torch.cuda.synchronize()
            begin=torch.cuda.Event(enable_timing=True);end=torch.cuda.Event(enable_timing=True)
            def run():
                begin.record()
                if mode=='graph':graph.replay()
                else:
                    for _ in range(nodes):op()
                end.record();end.synchronize()
            if measure=='profile':
                with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,torch.profiler.ProfilerActivity.CUDA]) as prof:run()
                with tempfile.TemporaryDirectory() as tmp:
                    path=Path(tmp)/'trace.json';prof.export_chrome_trace(str(path));events=json.loads(path.read_text())['traceEvents']
                kernels=[e for e in events if e.get('cat')=='kernel'];kernel_count=len(kernels)
                service=sum(e['dur']*1e-6 for e in kernels) if kernel_count==nodes else None
                del prof,events,kernels
            else:run();service=None;kernel_count=None
            rows.append(dict(name=name,cells=cells,nodes=nodes,repeat=repeat,mode=mode,instrument=measure,
                event_seconds=begin.elapsed_time(end)*1e-3,kernel_seconds=service,traced_kernel_count=kernel_count))
    del graph,op,x,y,z;gc.collect();torch.cuda.empty_cache()
write_record(root/'report.json',dict(context=context,observations=rows,telemetry=states,
    scope='Matched existing-PyTorch graph and eager generic operations. Event span and profiled kernel sum are separate. No correction factor is applied to any prediction.'))
for name in ['fill','add']:
    print(json.dumps(dict(name=name,summaries=[dict(mode=mode,instrument=instrument,
        event_seconds=statistics.median(r['event_seconds'] for r in rows if (r['name'],r['mode'],r['instrument'])==(name,mode,instrument)),
        kernel_seconds=statistics.median(r['kernel_seconds'] for r in rows if (r['name'],r['mode'],r['instrument'])==(name,mode,instrument)) if instrument=='profile' and all(r['kernel_seconds'] is not None for r in rows if (r['name'],r['mode'],r['instrument'])==(name,mode,instrument)) else None)
        for mode in ['graph','eager'] for instrument in ['plain','profile']])),flush=True)
