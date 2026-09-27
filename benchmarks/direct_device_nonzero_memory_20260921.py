"""Untimed installed-allocator requests for one bounded CUDA bool nonzero."""
import argparse
import gc
import json
from pathlib import Path
import torch
from torchgwas.detailed_calibration import sha256_file,source_identity
from torchgwas.geometry_collection import write_record
from torchgwas.selection_geometry import DEVICE_SELECTION_MAX_CELLS

parser=argparse.ArgumentParser()
parser.add_argument('--out',required=True)
parser.add_argument('--device',default='cuda:0')
parser.add_argument('--cells',nargs='+',type=int,default=[32,4093,1048576])
args=parser.parse_args()
if not args.cells or any(c<1 or c>DEVICE_SELECTION_MAX_CELLS for c in args.cells) or len(set(args.cells))!=len(args.cells):
    raise ValueError('Unique positive counts within the CUDA nonzero limit required')
root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
(root/'harness.py').write_bytes(Path(__file__).read_bytes())
torch.set_num_threads(2);torch.set_num_interop_threads(1)
torch.cuda.set_device(args.device)
prop=torch.cuda.get_device_properties(args.device)
# The dense control holds the mask and 16 bytes of coordinates per cell.
limit=max(256<<20,24*max(args.cells)+(64<<20))
torch.cuda.set_per_process_memory_fraction(limit/prop.total_memory,args.device)
if torch.cuda.mem_get_info(args.device)[0]<limit:raise ValueError('Insufficient free memory')
source=source_identity()
lib=Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device=args.device,
    device_uuid=str(prop.uuid),name=prop.name,compute_capability=[prop.major,prop.minor],
    sm_count=prop.multi_processor_count,library_sha256=sha256_file(lib),memory_limit_bytes=limit)
rows=[]
for cells in args.cells:
    for count in sorted({0,1,cells}):
        mask=torch.zeros((1,cells),device=args.device,dtype=torch.bool)
        if count:mask[0,:count]=True
        output=mask.nonzero(as_tuple=False)
        assert output.shape==(count,2)
        del output;torch.cuda.synchronize(args.device);gc.collect()
        baseline=torch.cuda.memory_allocated(args.device)
        torch.cuda.reset_peak_memory_stats(args.device)
        torch.cuda.memory._record_memory_history(enabled='all',context='all',stacks='all',max_entries=10000,device=args.device)
        output=mask.nonzero(as_tuple=False)
        torch.cuda.synchronize(args.device)
        peak=torch.cuda.max_memory_allocated(args.device)-baseline
        snapshot=torch.cuda.memory._snapshot()
        torch.cuda.memory._record_memory_history(enabled=None)
        trace=[{key:value for key,value in event.items() if key!='time_us'}
            for event in snapshot['device_traces'][torch.device(args.device).index]]
        allocations=[event for event in trace if event['action']=='alloc']
        frees=[event for event in trace if event['action']=='free_requested']
        output_ptr=output.data_ptr() if output.numel() else 0
        output_allocation_bytes=0
        for segment in snapshot['segments']:
            if segment['device']!=torch.device(args.device).index:continue
            cursor=segment['address']
            for block in segment['blocks']:
                address=block.get('address',cursor)
                if output_ptr and address==output_ptr:
                    assert block['state']=='active_allocated'
                    assert block['requested_size']==16*count
                    assert not output_allocation_bytes
                    output_allocation_bytes=block['size']
                cursor+=block['size']
        assert bool(output_allocation_bytes)==bool(count)
        row=dict(cells=cells,retained=count,shape=list(mask.shape),dtype=str(mask.dtype),
            output_ptr=output_ptr,output_bytes=16*count,output_allocation_bytes=output_allocation_bytes,extra_allocated_peak_bytes=peak,
            trace=trace,allocation_events=len(allocations),free_events=len(frees),
            duration_fields_recorded=False)
        rows.append(row);write_record(root/f'{cells}_{count}.json',row)
        print(json.dumps(dict(cells=cells,retained=count,peak_bytes=peak,output_allocation_bytes=output_allocation_bytes,
            allocation_sizes=[v['size'] for v in allocations],free_sizes=[v['size'] for v in frees])),flush=True)
        del output,mask,snapshot,trace,allocations,frees;gc.collect();torch.cuda.empty_cache()
assert source==source_identity()
write_record(root/'census.json',dict(context=context,source_sha256=source,rows=rows,
    duration_fields_recorded=False,
    scope='Exact installed allocator event requests and extra allocated peak for bounded contiguous 2D bool nonzero. Mask allocation excluded; coordinate output remains live. No timings, full scan, allocator-reservation claim or candidate-shape performance calibration.'))
