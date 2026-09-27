"""Attribute context-binding cost before proposing any identity-cache change."""
import argparse
import cProfile
import json
from pathlib import Path
import pstats
import time
from unittest.mock import patch

import torch
from torchgwas.api import run_linear_gwas  # Match the public API's loaded libraries.
from torchgwas import detailed_calibration as calibration
from torchgwas import gpu_identity,pgen_native
import threadpoolctl


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    source=calibration.source_identity();records=[]
    for phase in ('first','repeat1','repeat2'):
        calls=[]
        def instrument(name,function):
            def wrapped(*a,**kw):
                wall=time.perf_counter();cpu=time.thread_time()
                try:return function(*a,**kw)
                finally:
                    row=dict(name=name,wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu)
                    if name=='hash':
                        path=Path(a[0]);row.update(path=str(path),bytes=path.stat().st_size)
                    if name=='subprocess':row['command']=a[0]
                    calls.append(row)
            return wrapped
        initialized=torch.cuda.is_initialized();wall=time.perf_counter();cpu=time.thread_time()
        with patch.object(calibration,'sha256_file',instrument('hash',calibration.sha256_file)),\
             patch.object(calibration.subprocess,'check_output',instrument('subprocess',calibration.subprocess.check_output)),\
             patch.object(gpu_identity,'physical_gpu_identity',instrument('physical_gpu_identity',gpu_identity.physical_gpu_identity)),\
             patch.object(torch.cuda,'get_device_properties',instrument('cuda_properties',torch.cuda.get_device_properties)),\
             patch.object(threadpoolctl,'threadpool_info',instrument('cpu_pool_discovery',threadpoolctl.threadpool_info)),\
             patch.object(pgen_native,'load_library',instrument('native_library_load',pgen_native.load_library)):
            current=calibration.source_identity()
            execution=calibration.execution_context(['cuda:0','cuda:1'],input_path=args.input,output_path=out)
        row=dict(phase=phase,cuda_initialized_before=initialized,wall_seconds=time.perf_counter()-wall,
            cpu_seconds=time.thread_time()-cpu,calls=calls,context=execution)
        assert current==source
        records.append(row)
        summary={name:dict(calls=sum(c['name']==name for c in calls),
            wall_seconds=sum(c['wall_seconds'] for c in calls if c['name']==name),
            cpu_seconds=sum(c['cpu_seconds'] for c in calls if c['name']==name)) for name in sorted({c['name'] for c in calls})}
        print(json.dumps(dict(phase=phase,wall_seconds=row['wall_seconds'],cpu_seconds=row['cpu_seconds'],components=summary)),flush=True)
    profile=cProfile.Profile();profile.enable()
    calibration.source_identity()
    calibration.execution_context(['cuda:0','cuda:1'],input_path=args.input,output_path=out)
    profile.disable();profile.dump_stats(out/'binding.prof')
    entries=[]
    for (filename,line,function),(primitive,total,self_time,cumulative,_) in pstats.Stats(profile).stats.items():
        entries.append(dict(file=filename,line=line,function=function,primitive_calls=primitive,calls=total,
            self_seconds=self_time,cumulative_seconds=cumulative))
    entries.sort(key=lambda row:row['cumulative_seconds'],reverse=True)
    assert calibration.source_identity()==source
    report=dict(source_sha256=source,script_sha256=calibration.sha256_file(__file__),records=records,
        profile_top=entries[:40],input_mount=calibration.storage_identity(args.input),
        package_mount=calibration.storage_identity(calibration.__file__),
        scope='Instrumented source/context binding only, with public API libraries imported. Includes first CUDA property initialization where indicated. Per-call timers and path stats add overhead; component totals may nest. cProfile is diagnostic, not an uninstrumented speed measurement. No GWAS or independent capacity measurement.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
