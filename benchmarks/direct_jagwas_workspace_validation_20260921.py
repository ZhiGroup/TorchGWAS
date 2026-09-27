"""Bind untimed vendor queries to the source memory composition."""
import argparse
import hashlib
import json
from pathlib import Path
from torchgwas.cusolver_memory import jagwas_factor_workspace
from torchgwas.detailed_calibration import source_identity
from torchgwas.tensor_memory import eager_memory_plan

parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
parser.add_argument('--census',required=True);args=parser.parse_args()
out=Path(args.out)
if out.exists():raise FileExistsError(out)
path=Path(args.census);census=json.loads(path.read_text())
assert census['benchmark_sha256']==hashlib.sha256(Path('benchmarks/direct_jagwas_workspace_20260921.py').read_bytes()).hexdigest()
source=source_identity()
for name in ['reduce.py','linear.py','native_scan.py','preprocess.py','api.py']:
    assert source[name]==census['source_sha256'][name]
profile=dict(torch_version=census['torch_version'],cuda_version=census['cuda_version'],
    preferred_linalg=census['preferred_linalg'],compute_capability=census['compute_capability'],
    gpu_name=census['gpu'],host=census['host'],cusolver_library=census['library'],
    # The captured A100 resource geometry. No throughput is used in this audit.
    sm_count=108,max_threads_per_sm=2048,cublas_workspace_config=None,cublas_handle_stream_pairs=1,
    jagwas_factor_workspace_census=census)
assert profile['compute_capability']==[8,0] and 'A100' in profile['gpu_name']
checks=[]
for row in census['rows']:
    k=row['traits'];query=jagwas_factor_workspace(census,k,profile)
    plan=eager_memory_plan(max(257,2*k+1),13,k,2,2,profile,reduction='jagwas')
    assert plan['factor_workspace']==query
    assert plan['factor_host_workspace_bytes']==row['host_workspace_bytes']
    assert plan['setup_bytes']>=4*max(257,2*k+1)*k+plan['factor_setup']['distinct_temporary_bytes']+query['device_rounded_bytes']+plan['cublas']['total_bytes']
    assert not plan['prediction_complete']
    checks.append(dict(traits=k,device_query_bytes=query['device_requested_bytes'],
        device_rounded_bytes=query['device_rounded_bytes'],host_query_bytes=query['host_requested_bytes'],
        model_setup_bytes=plan['setup_bytes'],model_device_bytes=plan['device_bytes']))
assert all(row['remaining_within_4096_bytes'] for row in census['observations'])
report=dict(checks=checks,source_sha256=source,census_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    audit_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),durations_recorded=False,prediction_complete=False,
    scope='Six exact vendor queries integrated as factor-preparation requests. No interpolation, association timings, full memory guarantee or throughput estimate.')
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(dict(queried_shapes=len(checks),all_requests_composed=True)),flush=True)
