"""Untimed joint-ledger reconciliation and exact legacy default replays."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import types
import xml.etree.ElementTree as ET
from unittest.mock import patch
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record
from torchgwas.reduction_tensor_work import jagwas_tensor_work
from torchgwas.significant_schedule import significant_trait_schedule
from torchgwas.mechanistic_torch import torch_scan_work
from torchgwas.owned_result_work import owned_result_work
from torchgwas.native_control_work import native_control_work
from test_significant_schedule import tile,step
from test_mechanistic_shapes import fixture,component

parser=argparse.ArgumentParser()
parser.add_argument('--out',required=True)
parser.add_argument('--test-report',action='append',required=True)
parser.add_argument('--geometry',default='results/jagwas_geometry_20260921/census.json')
parser.add_argument('--chunk-geometry',default='tests/fixtures/jagwas_chunk_geometry.json')
args=parser.parse_args()
OUT=Path(args.out)
OUT.mkdir(parents=True,exist_ok=False)
source=source_identity()
baseline=json.loads(Path('docs/jagwas_development_baseline_20260921.json').read_text())['source_sha256']
geometry=json.loads(Path(args.geometry).read_text())
for name in ['reduce.py','linear.py','native_scan.py','preprocess.py','api.py']:
    assert source[name]==geometry['source_sha256'][name],name


def old_trace(value):
    # New fields describe host API grouping and tensor location. They must not
    # change the prior operation, storage identity, dtype or byte accounting.
    if isinstance(value,list):return [old_trace(v) for v in value]
    if isinstance(value,dict):
        drop={'host_calls','unpriced_host_calls','device','host_call_id','shape_only_inputs'}
        if 'op' in value:drop.add('phase')
        return {key:old_trace(v) for key,v in value.items() if key not in drop}
    return value


def operations(steps):
    result=[]
    for row in steps:
        if row['op'] in ['aten.detach.default','aten.lift_fresh.default']:continue
        if row['op']=='aten._to_copy.default' and row['inputs'][0].get('device_type')=='cpu':continue
        arrays=lambda key:[(v['shape'],v['dtype']) for v in row[key]]
        result.append((row['op'],arrays('inputs'),arrays('outputs')))
    return result


checks=[]
for row in geometry['rows']:
    trace=jagwas_tensor_work(row['N'],row['B'],row['K'],phase=row['phase'],compute_dtype=row['compute_dtype'])
    assert old_trace(trace)==old_trace(row['tensor_work'])
    assert operations(trace['steps'])==operations(row['observed_tensor_steps'])
    checks.append(dict(N=row['N'],B=row['B'],K=row['K'],phase=row['phase'],compute_dtype=row['compute_dtype'],
        operations_and_storage_unchanged=True,cuda_sequence_matches=True))
assert len(checks)==12

prior_root=Path('/data484_4/zxie3/torchGWAS1.1/src/torchgwas')
prior_hashes={}


def legacy(name):
    path=prior_root/(name+'.py');content=path.read_bytes()
    digest=hashlib.sha256(content).hexdigest()
    assert digest==baseline[name+'.py'],name+' frozen baseline changed'
    prior_hashes[name+'.py']=digest
    module=types.ModuleType('torchgwas._audit_baseline_'+name)
    module.__package__='torchgwas';module.__file__=str(path)
    sys.modules[module.__name__]=module
    exec(compile(content,str(path),'exec'),module.__dict__)
    return module


prior_schedule=legacy('significant_schedule').significant_trait_schedule
schedule_replays=[]
for backend in ['host','device']:
    for retained in [0,1]:
        for count in [1,2,4]:
            tiles=[tile(device='cuda:'+str(1+i%2),backend=backend,chunks=3,
                        kernel=1.+i,selection=2.,writer=5.+i,retained=retained) for i in range(count)]
            for depth in ([0] if count==1 else [1,3]):
                kwargs=dict(queue_depth=depth,shared_capacities={},
                    queue_service=dict(put=step(.2),get=step(.3)) if count>1 else None,finalize=step(.4))
                old=prior_schedule(copy.deepcopy(tiles),return_graph=True,**kwargs)
                new=significant_trait_schedule(copy.deepcopy(tiles),return_graph=True,**kwargs)
                assert old.__dict__==new.__dict__
                assert prior_schedule(copy.deepcopy(tiles),**kwargs)==significant_trait_schedule(copy.deepcopy(tiles),**kwargs)
                schedule_replays.append(dict(backend=backend,retained=retained,tiles=count,queue_depth=depth,
                    graph_demands_tokens_fifo_exact=True,solution_exact=True,nodes=len(new.nodes)))

prior_scan=legacy('mechanistic_torch');prior_scan.tensor_stage_service=component
scan_replays=[]
for k in [1,7,128]:
    data,profile=fixture(k=k,c=2)
    for copy_scenario in [False,True]:
        if copy_scenario:
            profile['owned_result_copy_scenario']=dict(resident_cpu_seconds_per_byte=1e-10,fresh_cpu_seconds_per_byte=8e-10,fresh_fraction=.5)
        with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
            assert torch_scan_work(data,profile)==prior_scan.torch_scan_work(data,profile)
        scan_replays.append(dict(K=k,copy_scenario=copy_scenario,all_fields_exact=True))
prior_owned=legacy('owned_result_work');prior_control=legacy('native_control_work')
for b in [1,32,1024]:
    for k in [1,7,2048]:assert owned_result_work(b,k)==prior_owned.owned_result_work(b,k)
for reused in [False,True]:assert native_control_work(reused)==prior_control.native_control_work(reused)

prior_tiling=legacy('trait_tiling_model')
from torchgwas.trait_tiling_model import trait_tiled_shape
from test_trait_tiling_model import candidate as trait_candidate
from test_pgen_native_reader import write_pgen
import tempfile
import numpy as np
trait_shape_replays=[]
with tempfile.TemporaryDirectory(prefix='joint_shape_audit_') as directory:
    path=Path(directory)/'source.pgen'
    write_pgen(path,np.zeros((10,32),np.uint8))
    for width in [1,2,5]:
        for count in [1,2]:
            if count>(5+width-1)//width:continue
            value=trait_candidate(path,width=width,count=count)
            assert trait_tiled_shape(value)==prior_tiling.trait_tiled_shape(value)
            trait_shape_replays.append(dict(traits=5,width=width,devices=count,all_fields_exact=True))
# The FP64 singleton/NSP extension must leave established FP32 services exact.
prior_tensor=legacy('tensor_service')
from torchgwas.tensor_service import tensor_stage_service,DeviceService,host_primitive_name
from torchgwas.tensor_work import eager_statistics_work
chunk_capture=json.loads(Path(args.chunk_geometry).read_text())
for name in ['reduce.py','linear.py','native_scan.py']:
    assert source[name]==chunk_capture['source_sha256'][name]
resources=dict(hbm_bytes_per_second=1e12,l2_bytes_per_second=2e12,fp32_flops_per_second=1e13,
    kernel_launch_seconds=1e-6,host_dispatch_cpu_seconds=1e-6,available_l2_bytes=40<<20,sm_count=108)
tensor_replays=[]
for row in chunk_capture['statistics']:
    if row['B']==1:continue  # Newly captured NSP path has no baseline support.
    work=eager_statistics_work(row['N'],row['B'],row['K'],row['C'],True)
    bank={host_primitive_name(call):1e-6 for call in work['host_calls']}
    old=prior_tensor.tensor_stage_service(work,prior_tensor.DeviceService(**resources),row['kernels'],host_primitives=bank)
    new=tensor_stage_service(work,DeviceService(**resources),row['kernels'],host_primitives=bank)
    assert old==new
    tensor_replays.append(dict(N=row['N'],B=row['B'],K=row['K'],kernels=len(row['kernels']),all_fields_exact=True))

tests=[]
for name in args.test_report:
    path=Path(name);suite=ET.parse(path).getroot().find('testsuite')
    assert int(suite.attrib['failures'])==int(suite.attrib['errors'])==0
    tests.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        **{key:suite.attrib[key] for key in ['tests','errors','failures','skipped']}))
assert source==source_identity()
write_record(OUT/'audit.json',dict(source_sha256=source,
    changed_source={name:dict(before=baseline.get(name),after=value) for name,value in source.items() if baseline.get(name)!=value},
    baseline_source_sha256=prior_hashes,geometry_checks=checks,
    geometry_artifacts={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in [args.geometry,args.chunk_geometry]},
    significant_schedule_exact_replays=schedule_replays,dense_scan_exact_replays=scan_replays,trait_shape_exact_replays=trait_shape_replays,
    tensor_service_exact_replays=tensor_replays,tests=tests,durations_recorded=False,prediction_complete=False,
    scope='Untimed operation/storage audit and exact default graph/work replays. Geometry execution code must match its capture; new host-call annotations cannot alter old tensor work. Synthetic services test accounting only, not runtime accuracy or autotune readiness.'))
print(json.dumps(dict(source_files=len(source),geometry_checks=len(checks),
    significant_graph_replays=len(schedule_replays),dense_scan_replays=len(scan_replays),trait_shape_replays=len(trait_shape_replays),tensor_service_replays=len(tensor_replays),tests=tests)),flush=True)
