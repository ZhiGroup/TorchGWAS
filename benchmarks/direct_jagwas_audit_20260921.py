"""Reconcile source-ledger changes with retained CUDA geometry and test artifacts."""
import ast
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record
from torchgwas.mechanistic_plan import _integer
from torchgwas.pinned_work import pinned_scan_work
from torchgwas.reduction_tensor_work import jagwas_tensor_work
from torchgwas.significant_host_work import indexed_part_work
from torchgwas.tensor_memory import eager_memory_plan

ROOT=Path('results/jagwas_source_audit_20260921');ROOT.mkdir(parents=True,exist_ok=True)
source=source_identity()
geometry=json.loads(Path('results/jagwas_geometry_20260921/census.json').read_text())
for name in ['reduce.py','linear.py','native_scan.py','preprocess.py','api.py','reduction_tensor_work.py']:
    if source[name]!=geometry['source_sha256'][name]:raise ValueError('Geometry execution source changed: '+name)


def normalize(steps):
    rows=[]
    for step in steps:
        if step['op'] in ['aten.detach.default','aten.lift_fresh.default']:continue
        if step['op']=='aten._to_copy.default' and step['inputs'][0].get('device_type')=='cpu':continue
        arrays=lambda key:[(v['shape'],v['dtype']) for v in step[key]]
        rows.append((step['op'],arrays('inputs'),arrays('outputs')))
    return rows


coverage=set();rows=[]
for row in geometry['rows']:
    key=(row['N'],row['B'],row['K'],row['compute_dtype'],row['phase'])
    assert key not in coverage;coverage.add(key)
    trace=jagwas_tensor_work(*key[:3],compute_dtype=key[3],phase=key[4])
    assert trace==row['tensor_work']
    assert normalize(trace['steps'])==normalize(row['observed_tensor_steps'])
    if row['phase']=='prepare':
        assert row['observed_extra_allocated_bytes']>=row['factor_capacity_floor']['explicit_live_bytes']
    rows.append(dict(shape=list(key[:3]),compute_dtype=key[3],phase=key[4],
        normalized_operations_match=True,kernel_count=len(row['kernels']),
        meta_dispatches=len(trace['steps']),observed_dispatches=len(row['observed_tensor_steps']),
        observed_extra_allocated_bytes=row['observed_extra_allocated_bytes']))
assert len(coverage)==12

# Moving NPZ accounting into the shared reduced-output module must preserve
# every old significant-part byte count. No association/runtime fit occurs.
prior=Path('/data484_4/zxie3/torchGWAS1.1/src/torchgwas/significant_host_work.py')
prior_bytes=prior.read_bytes();prior_sha=hashlib.sha256(prior_bytes).hexdigest()
assert prior_sha=='70d78a55ad51fa53c5720f5a5e155f67315701bba44da5605a6c2fc1f281aa2f'
node=next(n for n in ast.parse(prior_bytes).body if isinstance(n,ast.FunctionDef) and n.name=='indexed_part_work')
namespace=dict(io=io,np=np,_integer=_integer)
exec(compile(ast.Module(body=[node],type_ignores=[]),str(prior),'exec'),namespace)
compared=[]
for count in [0,1,37,4096,1<<20]:
    for store_beta in [False,True]:
        assert indexed_part_work(count,store_beta=store_beta)==namespace['indexed_part_work'](count,store_beta=store_beta)
        compared.append(dict(rows=count,store_beta=store_beta))
profile=dict(torch_version='2.5.1',sm_count=108,max_threads_per_sm=2048,
    compute_capability=[8,0],cublas_workspace_config=None,cublas_handle_stream_pairs=2)
memory=[]
for n,b,k,c in [(257,13,7,2),(2049,128,512,27),(4097,512,2048,27)]:
    plan=eager_memory_plan(n,b,k,c,3,profile,reduction='jagwas')
    pins=pinned_scan_work(n,b,k,3,reduction='jagwas')
    memory.append(dict(shape=[n,b,k,c],device=plan,pinned=pins))
xml=Path('results/jagwas_variant_checks_v8_20260921/tests.xml')
suite=ET.parse(xml).getroot().find('testsuite')
assert int(suite.attrib['failures'])==int(suite.attrib['errors'])==0
baseline=json.loads(Path('docs/jagwas_development_baseline_20260921.json').read_text())['source_sha256']
changed={name:dict(before=baseline.get(name),after=value) for name,value in source.items() if baseline.get(name)!=value}
assert source==source_identity()
write_record(ROOT/'audit.json',dict(source_sha256=source,changed_source=changed,
    geometry_source_sha256=geometry['source_sha256'],geometry_checks=rows,
    significant_npz_baseline_sha256=prior_sha,significant_npz_exact_replays=compared,
    memory_examples=memory,tests=dict(path=str(xml),sha256=hashlib.sha256(xml.read_bytes()).hexdigest(),
        **{key:suite.attrib[key] for key in ['tests','errors','failures','skipped']}),
    durations_recorded=False,prediction_complete=False,
    scope='Operation and storage audit only. Meta detach aliases and host-to-device upload are explicitly normalized against CUDA execution. No timing coefficients or readiness claim.'))
print(json.dumps(dict(geometry_cases=len(rows),significant_npz_exact_replays=len(compared),
    tests=suite.attrib['tests'],changed_source=list(changed))),flush=True)
