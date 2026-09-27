"""Freeze independent-resource predictions, then observe held-out selector kernels."""
import argparse,gc,json,random,statistics,tempfile,time
from pathlib import Path
import numpy as np
import torch
from torchgwas.detailed_calibration import sha256_file,source_identity
from torchgwas.device_significance_work import device_significant_tensor_work
from torchgwas.device_selection_gpu_work import selection_gpu_census,selection_gpu_work,selection_gpu_service
from torchgwas.tensor_service import DeviceService
from torchgwas.geometry_collection import kernel_census,write_record
from torchgwas.reduce import device_significant_pairs
from direct_device_significance_primitive_measure_20260921 import telemetry

p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--census',required=True)
p.add_argument('--resources',required=True);p.add_argument('--integer-capacity',required=True);a=p.parse_args()
root=Path(a.out);root.mkdir(parents=True,exist_ok=False);(root/'harness.py').write_bytes(Path(__file__).read_bytes())
census=json.loads(Path(a.census).read_text());capacities=json.loads(Path(a.resources).read_text())
integer=json.loads(Path(a.integer_capacity).read_text());context=capacities['context'];device=context['device']
integer_protocol=Path(a.integer_capacity).parent/'protocol.json'
assert integer['protocol_sha256']==sha256_file(integer_protocol)
for key in ['torch_version','cuda_runtime','device_uuid','compute_capability','sm_count','library_sha256']:
    assert context[key]==json.loads(integer_protocol.read_text())['context'][key]
torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.cuda.set_device(device)
prop=torch.cuda.get_device_properties(device);limit=512<<20
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,device)
if torch.cuda.mem_get_info(device)[0]<limit:raise ValueError('Insufficient free memory')
current=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device_uuid=str(prop.uuid),
    compute_capability=[prop.major,prop.minor],sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
for key,value in current.items():assert context[key]==value
source=source_identity();resources=DeviceService(**capacities['resources'])
integer_rate=next(r['element_pairs_per_second'] for r in integer['summaries'] if r['primitive']=='divmod_pair' and r['cells']==1<<20)
predictions=[]
for row in census['rows']:
    work=device_significant_tensor_work(row['N'],row['B'],row['K'],[b['retained'] for b in row['blocks']],max_cells=row['max_cells'])
    kernels=selection_gpu_census(work,census,current)
    ledger=selection_gpu_work(work,kernels,compute_capability=current['compute_capability'])
    predictions.append(dict(B=row['B'],K=row['K'],max_cells=row['max_cells'],mode=row['mode'],
        kernels=kernels,phases=[r['phase'] for r in ledger['kernels']],
        scenarios={mode:selection_gpu_service(ledger,resources,traffic_mode=mode,int64_divmod_per_second=integer_rate)
            for mode in ['logical_hbm','logical_l2']}))
# Immutable forecast is written before any observed selector runs in this process.
write_record(root/'predictions.json',dict(created=time.time(),context=current,source_sha256=source,
    census_sha256=sha256_file(a.census),resources_sha256=sha256_file(a.resources),integer_sha256=sha256_file(a.integer_capacity),
    harness_sha256=sha256_file(__file__),rows=predictions,
    scope='Held-out kernel-service transfer check only. No host barrier/dispatch/transfer wall-time prediction; cache scenarios are not bounds. No observations used to set prices.'))
frozen=sha256_file(root/'predictions.json');rng=random.Random(9220109);observations=[];states=[]
for index,row in enumerate(census['rows']):
    b,k,mode=row['B'],row['K'],row['mode'];ids=np.arange(b*k,dtype=np.int64).reshape(b,k)
    beta=(ids%1009).astype(np.float32);values=np.zeros((b,k),np.float32)
    if mode=='sparse':values.flat[::17]=2.
    elif mode in ('dense','invalid'):values.fill(2.)
    status=np.zeros(b,np.uint8);df=np.full(b,38.,np.float32)
    if mode=='invalid':
        status[1::5]=1;status[2::5]=2;df[::7]=0.;values.flat[::13]=np.nan
    table=np.ones(41,np.float32);table[0]=np.inf
    ix=np.nonzero(np.isfinite(values)&(np.abs(values)>=table[df.astype(np.int64),None])&(status[:,None]==0)&(df[:,None]>0))
    expected=[ix[0]+19,ix[1],beta[ix],values[ix],df[ix[0]]]
    tensors=[torch.as_tensor(v,device=device) for v in [beta,values,status,df,table]]
    def run():return list(device_significant_pairs(*tensors,start=19,max_cells=row['max_cells']))
    outputs=run();torch.cuda.synchronize(device)
    actual=[np.concatenate([part[i] for part in outputs]) for i in range(2,7)]
    order=np.lexsort((actual[1],actual[0]))
    for got,want in zip(actual,expected):np.testing.assert_array_equal(got[order],want)
    assert [len(part[2]) for part in outputs]==[r['retained'] for r in row['blocks']]
    del outputs,actual,expected,order
    for _ in range(3):outputs=run();torch.cuda.synchronize(device);del outputs
    states.append(dict(case=index,phase='before',telemetry=telemetry()))
    for repeat in range(5):
        modes=['plain','profile'];rng.shuffle(modes)
        for instrument in modes:
            gc.collect();torch.cuda.synchronize(device)
            def timed():
                t=time.perf_counter();cpu=time.thread_time();outputs=run()
                returned=time.perf_counter()-t;cpu_return=time.thread_time()-cpu
                torch.cuda.synchronize(device);complete=time.perf_counter()-t
                return outputs,returned,cpu_return,complete
            if instrument=='profile':
                with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU,torch.profiler.ProfilerActivity.CUDA]) as profiler:
                    outputs,returned,cpu_return,complete=timed()
                with tempfile.TemporaryDirectory() as tmp:
                    trace=Path(tmp)/'trace.json';profiler.export_chrome_trace(str(trace));events=json.loads(trace.read_text())['traceEvents']
                kernels=kernel_census(events)
                assert kernels==predictions[index]['kernels'], 'Changed compiled launch census'
                gpu=[e for e in events if e.get('cat')=='kernel']
                assert len(gpu)==len(kernels)
                times=[e['dur']*1e-6 for e in gpu]
                phase_times={phase:sum(t for t,p_ in zip(times,predictions[index]['phases']) if p_==phase)
                    for phase in set(predictions[index]['phases'])}
                del profiler,events,gpu
            else:
                outputs,returned,cpu_return,complete=timed()
                times=[];phase_times={}
            observations.append(dict(case=index,repeat=repeat,instrument=instrument,api_return_seconds=returned,
                api_return_cpu_seconds=cpu_return,complete_seconds=complete,kernel_seconds=times,phase_seconds=phase_times))
            del outputs
    states.append(dict(case=index,phase='after',telemetry=telemetry()))
    report=predictions[index]
    measured=statistics.median(sum(r['kernel_seconds']) for r in observations if r['case']==index and r['instrument']=='profile')
    print(json.dumps(dict(case=index,B=b,K=k,mode=mode,observed_gpu_seconds=measured,
        predicted_gpu_seconds={name:value['gpu_service_seconds'] for name,value in report['scenarios'].items()})),flush=True)
    write_record(root/f'observations_case_{index}.json',dict(predictions_sha256=frozen,observations=observations,telemetry=states,
        selector_arrays_equal=True,transfer_qualified=False,
        limitations=['Profiler kernel durations are not uninstrumented timing',
            'Paired plain/profile API intervals retain CPU and barriers; no subtraction or correction is fitted',
            'These synthetic selector cases do not establish GWAS candidate ranking']))
    del tensors,run;gc.collect();torch.cuda.empty_cache()
assert source==source_identity() and frozen==sha256_file(root/'predictions.json')

write_record(root/'observations.json',dict(predictions_sha256=frozen,observations=observations,telemetry=states,selector_arrays_equal=True,transfer_qualified=False))
