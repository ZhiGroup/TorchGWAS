"""Validate independent joint primitives; never apply prices to a scan silently."""
import argparse
import hashlib
import json
from pathlib import Path
from torchgwas.detailed_calibration import source_identity
from torchgwas.gil_service import gil_service_prices
from torchgwas.result_service import result_finish_prices
from torchgwas.reduction_tensor_work import jagwas_tensor_work,jagwas_host_primitive_name

parser=argparse.ArgumentParser()
parser.add_argument('--dispatch',required=True);parser.add_argument('--finish',required=True)
parser.add_argument('--out',required=True)
args=parser.parse_args()
out=Path(args.out)
if out.exists():raise FileExistsError(out)
dispatch=json.loads(Path(args.dispatch).read_text());finish=json.loads(Path(args.finish).read_text())
errors=[];dispatch_prices={};finish_prices={}
required=set()
for dtype in ['float32','float64']:
    trace=jagwas_tensor_work(64,32,32,phase='reduce',compute_dtype=dtype)
    required.update(jagwas_host_primitive_name(call) for call in trace['host_calls'])
source=source_identity()
for title,probe in [('dispatch',dispatch),('finish',finish)]:
    for filename,digest in probe['source_sha256'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest()!=digest:
            raise ValueError(title+' primitive source changed: '+filename)
if dispatch.get('primitive_bank')!='jagwas' or dispatch.get('fixed_shape')!=[32,32]:
    raise ValueError('Wrong independent fixed joint bank')
for context in dispatch['results']:
    for worker in context['workers']:
        key=context['mode']+':'+str(worker['device'])
        try:
            value=gil_service_prices(dispatch,mode=context['mode'],device=worker['device'],
                torch_version=dispatch['torch_version'],cpu_affinity=dispatch['affinity'])
            if not required<=set(value['cpu_primitives']):raise ValueError('Incomplete joint API coverage')
            dispatch_prices[key]=value
        except ValueError as error:errors.append(dict(context=key,error=str(error)))
for context in finish['results']:
    count=context['worker_count']
    try:
        finish_prices[str(count)]=result_finish_prices(finish,workers=count,
            numpy_version=finish['numpy_version'],torch_version=finish['torch_version'],
            python_version=finish['python_version'],cpu_affinity=finish['affinity'],
            source_sha256=source['native_scan.py'],reduction='jagwas')
    except ValueError as error:errors.append(dict(context='finish:'+str(count),error=str(error)))
report=dict(source_sha256=source,dispatch_file_sha256=hashlib.sha256(Path(args.dispatch).read_bytes()).hexdigest(),
    finish_file_sha256=hashlib.sha256(Path(args.finish).read_bytes()).hexdigest(),
    host=dispatch['host'],affinity=dispatch['affinity'],devices=dispatch['devices'],
    dispatch_prices=dispatch_prices,finish_prices=finish_prices,errors=errors,
    measurement_controls_passed=not errors,prediction_complete=False,
    total_cpu_transfer_qualified=False,serial_transfer_qualified=False,
    qualification_required='Matched unhooked CPU and serial-partition controls; meter validity alone is insufficient',
    scope='Independent fixed32 API CPU and five-array finish observations in the captured host/runtime/load context. '
          'Prices are not transferred to another host or silently applied to scan predictions. GPU throughput, '
          'factorization, writer, status-size scaling and loaded-context transfer still require separate evidence.')
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(dict(measurement_controls_passed=not errors,dispatch_contexts=list(dispatch_prices),
    finish_contexts=list(finish_prices),errors=errors)),flush=True)
if errors:raise SystemExit(2)
