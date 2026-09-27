"""Compare exact context bindings with optional immutable file-digest reuse."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import torch
from torchgwas.api import run_linear_gwas  # Match the public API's loaded pools.
from torchgwas import detailed_calibration as calibration
from torchgwas.binding_digests import BindingDigestCache


def artifacts(directory):
    return {str(p):calibration.sha256_file(p) for p in directory.glob('**/*.json')}


def binding(args,*,cached):
    initialized=torch.cuda.is_initialized();wall=time.perf_counter();cpu=time.thread_time()
    cache=BindingDigestCache(args.cache) if cached else None
    source=calibration.source_identity(digest_cache=cache)
    context=calibration.execution_context(['cuda:0','cuda:1'],input_path=args.input,
        output_path=args.cache,digest_cache=cache)
    bound_wall=time.perf_counter()-wall;bound_cpu=time.thread_time()-cpu
    endwall=time.perf_counter();endcpu=time.thread_time()
    assert calibration.source_identity(digest_cache=cache)==source
    if cache is not None:assert cache.unchanged()
    recheck_wall=time.perf_counter()-endwall;recheck_cpu=time.thread_time()-endcpu
    endwall=time.perf_counter();endcpu=time.thread_time()
    publication=None;state=None
    if cache is not None:
        publication=cache.publish(successful=True);state=cache.snapshot();cache.close()
    publication_wall=time.perf_counter()-endwall;publication_cpu=time.thread_time()-endcpu
    return dict(cached=cached,cuda_initialized_before=initialized,
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu,
        binding_wall_seconds=bound_wall,binding_cpu_seconds=bound_cpu,
        recheck_wall_seconds=recheck_wall,recheck_cpu_seconds=recheck_cpu,
        publication_wall_seconds=publication_wall,publication_cpu_seconds=publication_cpu,
        state=state,publication=publication,source_sha256=source,context=context)


def worker(args):
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    if args.worker=='warm_pairs':
        warm=binding(args,cached=True);rows=[]
        for pair in range(6):
            for cached in ((False,True) if pair%2==0 else (True,False)):
                row=binding(args,cached=cached);row['pair']=pair
                assert row['source_sha256']==warm['source_sha256'] and row['context']==warm['context']
                if cached:
                    assert row['state']['disk_hits']==3 and row['state']['hashed_files']==0
                    assert not row['publication']['stored']
                rows.append(row)
        result=dict(warmup=warm,records=rows)
    else:result=binding(args,cached=args.worker!='control')
    (Path(args.out)/(args.worker+'.json')).write_text(json.dumps(result,indent=2)+'\n')


def main(args):
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    cache=Path(args.cache);cache.mkdir(parents=True,exist_ok=False)
    before=calibration.source_identity();reports={};saved=None
    for label in ('control','populate','reuse','warm_pairs'):
        subprocess.run([sys.executable,__file__,'--input',args.input,'--cache',args.cache,
            '--out',args.out,'--worker',label],check=True)
        row=json.loads((out/(label+'.json')).read_text());reports[label]=row
        if label=='warm_pairs':assert artifacts(cache)==saved;continue
        assert row['source_sha256']==before
        if label=='control':assert not artifacts(cache)
        else:
            assert row['context']==reports['control']['context']
            if label=='populate':
                assert row['state']['disk_hits']==0 and len(row['publication']['stored'])==3
                saved=artifacts(cache);assert len(saved)==3
            else:
                assert row['state']['disk_hits']==3 and row['state']['hashed_files']==0
                assert not row['publication']['stored'] and artifacts(cache)==saved
        print(json.dumps(dict(phase=label,wall_seconds=row['wall_seconds'],cpu_seconds=row['cpu_seconds'],
            state=row['state'])),flush=True)
    assert calibration.source_identity()==before
    report=dict(source_sha256=before,script_sha256=calibration.sha256_file(__file__),records=reports,
        artifacts=saved,mounts={name:calibration.storage_identity(path) for name,path in
            [('input',args.input),('cache',args.cache),('package',calibration.__file__)]},
        scope='Binding-only audit with public API imports. Separate-process control/population/reuse include first CUDA property initialization. Six alternating warmed pairs each include construction, source/context binding, final source/file checks and deferred publication. Exact hashes and contexts agree; saved records remain byte-identical. Timings exclude imports and scientific work; no GWAS throughput or model-accuracy claim.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--cache',required=True)
    p.add_argument('--out',required=True);p.add_argument('--worker');args=p.parse_args()
    worker(args) if args.worker else main(args)
