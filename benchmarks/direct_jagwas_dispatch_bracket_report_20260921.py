"""Validate raw dispatch controls and compare with an unhooked CPU bracket."""
import argparse
import hashlib
import json
from pathlib import Path
from torchgwas.gil_probe_qualification import compare_gil_probe_bracket

parser=argparse.ArgumentParser()
parser.add_argument('--root',required=True)
parser.add_argument('--maximum-relative-difference',required=True,type=float)
args=parser.parse_args();root=Path(args.root);out=root/'comparison.json'
if out.exists():raise FileExistsError(out)
paths=[root/(name+'.json') for name in ['before','audited','after']]
probes=[json.loads(path.read_text()) for path in paths]
for probe in probes:
    for filename,digest in probe['source_sha256'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest()!=digest:
            raise ValueError('Measurement source changed: '+filename)
try:
    report=compare_gil_probe_bracket(*probes,maximum_relative_difference=args.maximum_relative_difference)
    report['measurement_controls_passed']=True
except ValueError as error:
    report=dict(measurement_controls_passed=False,total_cpu_compatible=False,serial_transfer_qualified=False,
        prediction_complete=False,error=str(error),contexts={},
        scope='Rejected measurement evidence. No clipping, dropped repeats or transfer of its CPU/GIL prices.')
report['artifact_sha256']={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
report['qualification_source_sha256']={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [
    Path(__file__),Path('src/torchgwas/gil_probe_qualification.py'),Path('src/torchgwas/gil_service.py')]}
with out.open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(dict(measurement_controls_passed=report['measurement_controls_passed'],
    error=report.get('error'),total_cpu_compatible=report['total_cpu_compatible'],
    serial_transfer_qualified=False,contexts={name:dict(total_cpu_compatible=context['total_cpu_compatible'],
        incompatible_primitives=[key for key,value in context['primitives'].items() if not value['total_cpu_compatible']])
        for name,context in report['contexts'].items()})),flush=True)

if not report['measurement_controls_passed']:raise SystemExit(2)
