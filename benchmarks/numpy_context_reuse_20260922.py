"""Fresh-process runtime binding audit; no association or capacity measurement."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def write(path,value):
    with Path(path).open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def worker(args):
    import numpy as np
    import torch
    from torchgwas.api import run_linear_gwas  # Match public library imports.
    from torchgwas.binding_digests import BindingDigestCache
    from torchgwas import detailed_calibration as binding
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    cache=BindingDigestCache(args.cache)
    start=time.perf_counter();cpu=time.process_time()
    sources=binding.source_identity(digest_cache=cache)
    context=binding.execution_context(['cuda:1','cuda:2'],input_path=args.input,output_path=args.input,digest_cache=cache)
    seconds=time.perf_counter()-start;cpu_seconds=time.process_time()-cpu
    late=None
    if args.late_declaration:
        os.environ['NPY_DISABLE_CPU_FEATURES']='AVX2'
        late=binding.execution_context(['cuda:1','cuda:2'],input_path=args.input,output_path=args.input,digest_cache=cache)
        assert late['numpy_core']==context['numpy_core']
        assert late['environment']['NPY_DISABLE_CPU_FEATURES']=='AVX2'
    published=cache.publish(successful=True)
    assert published['status']=='checked'
    write(args.out,dict(context=context,late_context=late,source_sha256=sources,
        digest_state=cache.snapshot(),publication=published,binding_wall_seconds=seconds,
        binding_cpu_seconds=cpu_seconds,numpy_version=np.__version__,
        script_sha256=binding.sha256_file(__file__)))


def main(args):
    from torchgwas import detailed_calibration as binding
    directory=Path(args.out);directory.mkdir(parents=True,exist_ok=False)
    cache=directory/'cache'
    saved=None;reports=[]
    for mode in ['first','reuse','disabled']:
        environment=os.environ.copy()
        # Both feature controls are import-time inputs. Do not inherit a
        # parent's choice and mislabel the baseline.
        environment.pop('NPY_DISABLE_CPU_FEATURES',None)
        environment.pop('NPY_ENABLE_CPU_FEATURES',None)
        if mode=='disabled':environment['NPY_DISABLE_CPU_FEATURES']='AVX2'
        path=directory/(mode+'.json')
        command=[sys.executable,__file__,'--worker','--input',args.input,
                 '--cache',str(cache),'--out',str(path)]
        if mode=='reuse':command.append('--late-declaration')
        subprocess.run(command,env=environment,check=True)
        reports.append(json.loads(path.read_text()))
        artifacts={str(p):p.read_bytes() for p in cache.rglob('*.json')}
        if saved is None:saved=artifacts
        else:assert artifacts==saved,'Reusing digests changed an immutable record'
    first,reuse,disabled=reports
    assert first['context']==reuse['context']
    assert first['source_sha256']==reuse['source_sha256']==disabled['source_sha256']
    assert first['context']['numpy_core']['cpu_features']['AVX2'] is True
    assert disabled['context']['numpy_core']['cpu_features']['AVX2'] is False
    assert first['context']['numpy_core']['library_sha256']==disabled['context']['numpy_core']['library_sha256']
    assert all(row['digest_state']['hashed_files']==0 for row in reports[1:])
    assert all(not row['publication']['stored'] for row in reports[1:])
    profile=binding.bind_detailed_profile(
        [dict(name='identity-audit',devices=['cuda:1','cuda:2'],profiles={'cuda:1':{},'cuda:2':{}})],
        first['context'],sources=first['source_sha256'],
        component_artifacts={str(directory/'first.json'):binding.sha256_file(directory/'first.json')},
        limitations=['Identity audit only: no empirical component prices or runtime qualification.'])
    binding.validate_detailed_profile(profile,reuse['context'],sources=reuse['source_sha256'])
    rejections={}
    for label,current in [('disabled',disabled['context']),('late_declaration',reuse['late_context'])]:
        try:binding.validate_detailed_profile(profile,current,sources=first['source_sha256'])
        except ValueError as error:rejections[label]=str(error)
        else:raise AssertionError('Changed NumPy execution context accepted')
    assert 'numpy_core.cpu_features.AVX2' in rejections['disabled']
    assert 'environment.NPY_DISABLE_CPU_FEATURES' in rejections['late_declaration']
    report=dict(unchanged_context_reused=True,original_digest_records_unchanged=True,
        immutable_digest_records=len(saved),rejections=rejections,
        processes=[dict(mode=mode,binding_wall_seconds=row['binding_wall_seconds'],
            binding_cpu_seconds=row['binding_cpu_seconds'],digest_state=row['digest_state'])
            for mode,row in zip(['first','reuse','disabled'],reports)],
        source_sha256=first['source_sha256'],script_sha256=binding.sha256_file(__file__),
        scope='Fresh-process context and immutable-digest audit only. No capacity measurement or GWAS speedup claim.')
    write(directory/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256','processes']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True)
    parser.add_argument('--out',required=True);parser.add_argument('--cache')
    parser.add_argument('--worker',action='store_true');parser.add_argument('--late-declaration',action='store_true')
    args=parser.parse_args()
    if args.worker:worker(args)
    else:main(args)
