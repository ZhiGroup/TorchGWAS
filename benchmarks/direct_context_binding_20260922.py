"""Attribute context-binding startup in fresh processes; no GWAS capacity fit."""
import argparse
import cProfile
import json
import os
from pathlib import Path
import pstats
import subprocess
import sys
import time
from unittest.mock import patch

import torch
import torchgwas.detailed_calibration as binding
import torchgwas.gpu_identity as gpu_identity
# Match the public API's imports before the binding boundary in both arms.
from torchgwas.api import run_linear_gwas


def worker(args):
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    rows=[]
    original_hash=binding.sha256_file
    original_check=binding.subprocess.check_output
    for iteration in range(2):
        hashes=[];commands=[]
        def digest(path):
            started=time.perf_counter();value=original_hash(path)
            hashes.append(dict(path=str(path),bytes=Path(path).stat().st_size,
                seconds=time.perf_counter()-started))
            return value
        def check(command,*a,**kw):
            started=time.perf_counter();value=original_check(command,*a,**kw)
            commands.append(dict(command=command,seconds=time.perf_counter()-started))
            return value
        profiler=cProfile.Profile()
        native=gpu_identity._nvml_identity
        def query(uuids):
            if args.backend=='smi':raise OSError('benchmark subprocess arm')
            return native(uuids)
        with patch.object(binding,'sha256_file',digest),patch.object(binding.subprocess,'check_output',check),\
             patch.object(gpu_identity,'_nvml_identity',query):
            profiler.enable();started=time.perf_counter()
            sources=binding.source_identity();source_seconds=time.perf_counter()-started
            started=time.perf_counter()
            context=binding.execution_context(['cuda:0','cuda:1'],input_path=args.input,output_path=args.data)
            context_seconds=time.perf_counter()-started;profiler.disable()
        stats=pstats.Stats(profiler)
        functions=[dict(file=key[0],line=key[1],function=key[2],calls=value[1],
            self_seconds=value[2],cumulative_seconds=value[3]) for key,value in stats.stats.items()]
        rows.append(dict(iteration=iteration,source_seconds=source_seconds,context_seconds=context_seconds,
            hashes=hashes,commands=commands,functions=sorted(functions,key=lambda f:-f['cumulative_seconds'])[:45],
            source_sha256=sources,context=context))
    uuids=[g['uuid'] for g in context['devices'].values()]
    native=native(uuids);fallback=gpu_identity._smi_identity(uuids)
    assert native==fallback
    result=dict(rows=rows,backend=args.backend,nvml_equals_smi=True,source_equal=rows[0]['source_sha256']==rows[1]['source_sha256'],
        context_equal=rows[0]['context']==rows[1]['context'],affinity=sorted(os.sched_getaffinity(0)),
        script_sha256=original_hash(__file__),scope='Attributed binding wall time, imports excluded; profiler and hash instrumentation overhead included. Not a GWAS speedup measurement.')
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--data',required=True)
    p.add_argument('--out',required=True);p.add_argument('--worker',action='store_true')
    p.add_argument('--backend',choices=['smi','nvml'],default='nvml');args=p.parse_args()
    if args.worker:worker(args)
    else:
        directory=Path(args.out);directory.mkdir(parents=True,exist_ok=False)
        for index in range(3):
            for backend in (['smi','nvml'] if index%2==0 else ['nvml','smi']):
                subprocess.run([sys.executable,__file__,'--worker','--input',args.input,'--data',args.data,
                    '--backend',backend,'--out',str(directory/f'process{index}_{backend}.json')],check=True)
        reports=[json.loads(path.read_text()) for path in sorted(directory.glob('process*.json'))]
        assert all(row['context']==reports[0]['rows'][0]['context'] for report in reports for row in report['rows'])
        (directory/'report.json').write_text(json.dumps(reports,indent=2)+'\n')
        for i,report in enumerate(reports):
            for row in report['rows']:
                print(i,report['backend'],row['iteration'],'source',row['source_seconds'],'context',row['context_seconds'],
                    'hash',sum(h['seconds'] for h in row['hashes']),'commands',row['commands'],flush=True)
