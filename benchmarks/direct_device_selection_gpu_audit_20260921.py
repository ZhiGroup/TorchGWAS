"""Audit every untimed selector launch against analytical CUB phase work."""
import argparse
import json
from pathlib import Path
import torch
from torchgwas.device_significance_work import device_significant_tensor_work
from torchgwas.device_selection_gpu_work import selection_gpu_census,selection_gpu_work
from torchgwas.detailed_calibration import sha256_file,source_identity
from torchgwas.geometry_collection import write_record

parser=argparse.ArgumentParser();parser.add_argument('--census',required=True);parser.add_argument('--out',required=True)
args=parser.parse_args();census=json.loads(Path(args.census).read_text())
device=census['context']['device'];prop=torch.cuda.get_device_properties(device)
context=dict(torch_version=torch.__version__,cuda_runtime=torch.version.cuda,device_uuid=str(prop.uuid),
    compute_capability=[prop.major,prop.minor],sm_count=prop.multi_processor_count,
    library_sha256=sha256_file(Path(torch.__file__).parent/'lib'/'libtorch_cuda.so'))
source=source_identity();reports=[]
for row in census['rows']:
    work=device_significant_tensor_work(row['N'],row['B'],row['K'],[b['retained'] for b in row['blocks']],max_cells=row['max_cells'])
    kernels=selection_gpu_census(work,census,context)
    ledger=selection_gpu_work(work,kernels,compute_capability=context['compute_capability'])
    assert ledger['kernel_count']==len(kernels)
    report=dict(N=row['N'],B=row['B'],K=row['K'],mode=row['mode'],kernel_count=ledger['kernel_count'],
        groups=[{phase:len(indices) for phase,indices in group.items()} for group in ledger['nonzero_groups']],
        logical_gpu_bytes=ledger['logical_bytes'],source_sha256=work['source_sha256'])
    reports.append(report);print(json.dumps(report),flush=True)
assert source==source_identity()
write_record(args.out,dict(census_sha256=sha256_file(args.census),source_sha256=source,
    context=context,rows=reports,durations_recorded=False,
    scope='Every compiled selector launch assigned once to its source phase; explicit source policy and context checks. No timing qualification.'))
