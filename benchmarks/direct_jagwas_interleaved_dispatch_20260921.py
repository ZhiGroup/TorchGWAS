"""Interleave independent fixed-API controls in two persistent Python processes.

Only one child measures at a time. Both complete an excluded full warm-up, so
interpreter/CUDA startup is outside every comparison. No association is run.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

THREAD_ENV=dict(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
    OMP_WAIT_POLICY='PASSIVE',GOMP_SPINCOUNT='0')
os.environ.update(THREAD_ENV)


def child(args):
    import threading
    from concurrent.futures import ThreadPoolExecutor
    import torch
    if args.worker=='audited':
        from direct_calculator_gil_audit import worker
        if os.environ.get('LD_AUDIT')!=args.library:raise ValueError('Exact audit library required')
    else:
        from direct_jagwas_host_cpu_control_20260921 import worker
        if os.environ.get('LD_AUDIT'):raise ValueError('Unhooked process required')
    os.sched_setaffinity(0,range(12,20));torch.set_num_threads(4)
    print(json.dumps(dict(ready=args.worker)),flush=True)
    for line in sys.stdin:
        request=json.loads(line);path=Path(request['out'])
        if path.exists():raise FileExistsError(path)
        path.parent.mkdir(parents=True,exist_ok=True);results=[]
        for mode,devices in [('single0',[0]),('threads',[0,2])]:
            barrier=threading.Barrier(len(devices))
            with ThreadPoolExecutor(max_workers=len(devices)) as pool:
                futures=[]
                for device in devices:
                    checkpoint=str(path.with_name(path.stem+'.'+mode+'.'+str(device)+'.json'))
                    if args.worker=='audited':
                        future=pool.submit(worker,device,barrier,args.library,'jagwas',checkpoint)
                    else:future=pool.submit(worker,device,barrier,args.library,checkpoint)
                    futures.append(future)
                results.append(dict(mode=mode,workers=[future.result() for future in futures]))
        paths=[Path(__file__),Path('benchmarks/direct_jagwas_host_primitives.py'),
            Path('benchmarks/direct_calculator_gil_audit.c'),Path(args.library),
            Path('benchmarks/direct_calculator_gil_audit.py' if args.worker=='audited'
                else 'benchmarks/direct_jagwas_host_cpu_control_20260921.py')]
        report=dict(results=results,host=os.uname().nodename,torch_version=torch.__version__,
            affinity=sorted(os.sched_getaffinity(0)),thread_environment=THREAD_ENV,
            primitive_bank='jagwas',fixed_shape=[32,32],
            devices={str(d):dict(name=torch.cuda.get_device_name(d),capability=list(torch.cuda.get_device_capability(d))) for d in [0,2]},
            source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            scope='Fixed32 APIs in a persistent interpreter; alternating separately hooked and unhooked processes. Existing raw observation and meter boundaries are unchanged. Idle companion CUDA contexts persist. Shared server and phase drift remain confounders.')
        if args.worker=='audited':
            report['context_verified']=all(w['context_verified'] for context in results for w in context['workers'])
        else:report['ld_audit']=None
        with path.open('x') as stream:json.dump(report,stream,indent=2)
        print(json.dumps(dict(complete=str(path))),flush=True)


def parent(args):
    root=Path(args.out)
    root.mkdir(parents=True,exist_ok=False)
    library=str(Path(args.library).resolve());children={};logs=[];order=[]
    try:
        for kind in ['unhooked','audited']:
            env=dict(os.environ);env.pop('LD_AUDIT',None)
            if kind=='audited':env['LD_AUDIT']=library
            log=(root/(kind+'.stderr.log')).open('w');logs.append(log)
            process=subprocess.Popen([sys.executable,__file__,'--worker',kind,'--library',library],
                stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log,text=True,env=env)
            children[kind]=process
            message=process.stdout.readline()
            if not message or json.loads(message).get('ready')!=kind:
                raise RuntimeError('Worker failed to start: '+kind)
        def run(kind,name):
            path=root/(name+'.json');process=children[kind]
            process.stdin.write(json.dumps(dict(out=str(path)))+'\n');process.stdin.flush()
            message=process.stdout.readline()
            if not message or json.loads(message).get('complete')!=str(path):
                raise RuntimeError('Worker failed: '+kind+' '+name+'; see retained stderr/checkpoints')
            order.append(dict(kind=kind,path=str(path)))
            print(json.dumps(dict(completed=name,kind=kind)),flush=True)
            return path
        run('unhooked','warm_unhooked');run('audited','warm_audited')
        control=run('unhooked','control_0')
        for index in range(args.rounds):
            audit=run('audited','audit_'+str(index))
            after=run('unhooked','control_'+str(index+1))
            triplet=root/('triplet_'+str(index));triplet.mkdir()
            for source,name in [(control,'before'),(audit,'audited'),(after,'after')]:
                shutil.copyfile(source,triplet/(name+'.json'))
            control=after
        (root/'order.json').write_text(json.dumps(order,indent=2))
    finally:
        for process in children.values():
            if process.stdin:process.stdin.close()
        for process in children.values():process.wait()
        for log in logs:log.close()
    # Comparison runs only after both measurement processes have exited.
    for index in range(args.rounds):
        subprocess.run([sys.executable,'benchmarks/direct_jagwas_dispatch_bracket_report_20260921.py',
            '--root',str(root/('triplet_'+str(index))),'--maximum-relative-difference','0.10'],check=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--library',required=True)
    parser.add_argument('--out');parser.add_argument('--rounds',type=int,default=3)
    parser.add_argument('--worker',choices=['audited','unhooked']);args=parser.parse_args()
    if args.worker:child(args)
    else:
        if not args.out or args.rounds<1:raise ValueError('Output root and positive rounds required')
        parent(args)


if __name__=='__main__':main()
