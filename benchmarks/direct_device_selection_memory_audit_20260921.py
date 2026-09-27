"""Reconcile exact CUB requests and a conservative selector tensor budget."""
import argparse
import json
from pathlib import Path
from torchgwas.nonzero_memory import device_nonzero_workspace,device_selection_memory
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record

parser=argparse.ArgumentParser()
parser.add_argument('--census',required=True)
parser.add_argument('--geometry',required=True)
parser.add_argument('--out',required=True)
args=parser.parse_args()
census=json.loads(Path(args.census).read_text())
geometry=json.loads(Path(args.geometry).read_text())
source=source_identity()
for name in ('reduce.py','selection_geometry.py'):
    assert geometry['source_sha256'][name]==source[name]
    assert census['source_sha256'][name]==source[name]
assert geometry['torch_version']==census['context']['torch_version']
assert geometry['device']==census['context']['name']
requests=[device_nonzero_workspace(census,cells,census['context']) for cells in sorted({r['cells'] for r in census['rows']})]
checks=[]
for row in geometry['rows']:
    if (row['B'],row['K'],row['max_cells'])!=(257,4093,1<<20):continue
    work=device_selection_memory(40,row['B'],row['K'],census=census,context=census['context'],max_cells=row['max_cells'])
    assert row['extra_allocated_bytes']<=work['selection_gpu_bytes']
    checks.append(dict(mode=row['mode'],observed_extra_bytes=row['extra_allocated_bytes'],selection_gpu_bytes=work['selection_gpu_bytes'],
        all_source_block_extents=[(b['rows'],b['traits'],b['cells']) for b in work['block_shapes']],
        single_block_tensor_bytes=work['single_block_tensor_bytes'],single_block_private_bytes=work['single_block_private_bytes']))
assert len(checks)==4
write_record(Path(args.out),dict(source_sha256=source,census_sha256=sha256_file(args.census),geometry_sha256=sha256_file(args.geometry),
    requests=requests,checks=checks,scope='Every captured workspace request reconciles across empty/single/dense occupancy; conservative block-local tensor budget exceeds the observed allocation peak in all four existing large selector cases. These observations do not establish an allocator-reserved whole-job bound.'))
print(json.dumps(dict(workspace_extents=len(requests),workspace_occupancy_cases=len(census['rows']),selector_peak_checks=checks)),flush=True)
