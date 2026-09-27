"""Fresh-process public-GWAS control of experimental CPU selected-result packing.

Fixed two-GPU layout, warmed server-local inputs, durable significant output.
No calculator fitting, price installation or production source replacement.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import sys
import time
from unittest.mock import patch

from selector_pack_control_20260922 import PackControl,sha
from native_host_predicate_gwas_20260922 import output
import numpy as np
import torch
from torchgwas import host_significance
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity,execution_context

ROOT=Path('results/selector_pack_gwas_control_20260922')
TARGET=Path('/data/zxie3/torchgwas_selector_pack_control_20260922')
SOURCE=Path('/data/zxie3/torchgwas_public_refresh_v3_20260922')
PHENO=Path('/data/zxie3/torchgwas_adaptive_candidate_fixture_v1_20260922')
PROTOTYPE=Path('results/selector_pack_control_20260922')
CASES={'empty_t':(1e-30,'t'),'sparse_beta_t':(.02,'beta+t'),
    'dense_t':(1.,'t'),'dense_beta_t':(1.,'beta+t')}


def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def read_control(path):
    began=time.perf_counter();count=0;digest=hashlib.sha256()
    with path.open('rb') as stream:
        while True:
            block=stream.read(1<<20)
            if not block:break
            count+=len(block);digest.update(block)
    return dict(bytes=count,wall_seconds=time.perf_counter()-began,sha256=digest.hexdigest(),
        scope='Sequential warm-cache read plus hashing, not a storage-bandwidth ceiling.')


def child(case,implementation,repeat):
    plan=json.loads((ROOT/'plan.json').read_text())
    assert source_identity()==plan['source_sha256']
    assert {p:sha(p) for p in plan['inputs']}==plan['inputs']
    assert {p:sha(p) for p in plan['harness_sha256']}==plan['harness_sha256']
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False;os.sched_setaffinity(0,list(range(12,20)))
    devices=['cuda:1','cuda:2']
    for device in devices:
        values=torch.ones((32,32),device=device);answer=values@values
        torch.cuda.synchronize(device);del values,answer
        torch.cuda.reset_peak_memory_stats(device)
    current=host_significance.select_host_pairs
    function=current if implementation=='current' else PackControl((PROTOTYPE/'pack.so').resolve())
    threshold,fields=CASES[case];label=f'{case}_{repeat}_{implementation}'
    y=np.load(PHENO/'phenotype.npy');cov=np.load(PHENO/'covariates.npy')
    before=read_control(SOURCE/'input.pgen')
    context=execution_context(devices,input_path=SOURCE/'input.pgen',output_path=TARGET)
    began=time.perf_counter();cpu=time.process_time()
    with patch.object(host_significance,'select_host_pairs',function):
        result=run_linear_gwas(str(SOURCE/'input.pgen'),y,cov,genotype_format='pgen',pgen_mode='hardcall',
            device='cuda:1',compute_dtype='float32',chunk_size=4096,reader_workers=4,prefetch_chunks=2,
            reduce='significant',significance_threshold=threshold,trait_block=256,trait_devices=devices,
            output_dir=TARGET/label,sumstats_format='binary',sumstats_fields=fields,
            sumstats_queue_depth=2,sumstats_fsync=True)
    cpu_seconds=time.process_time()-cpu;elapsed=time.perf_counter()-began
    after=read_control(SOURCE/'input.pgen')
    assert before['sha256']==after['sha256']==plan['inputs'][str(SOURCE/'input.pgen')]
    actual=output(TARGET/label/'sumstats')
    hashes={key:dict(shape=list(value.shape),dtype=value.dtype.str,
        sha256=hashlib.sha256(memoryview(value)).hexdigest()) for key,value in actual.items()}
    assert host_significance.select_host_pairs is current and source_identity()==plan['source_sha256']
    report=dict(case=case,implementation=implementation,repeat=repeat,pid=os.getpid(),
        observed_at_utc=datetime.now(timezone.utc).isoformat(),execution_context=context,
        api_seconds=elapsed,api_process_cpu_seconds=cpu_seconds,read_before=before,read_after=after,
        result_rows=result.run_metadata['n_result_rows'],writer=result.run_metadata['sumstats_write'],
        phase_breakdown=result.run_metadata.get('phase_breakdown'),output=str(TARGET/label),
        peak_allocated_bytes={d:torch.cuda.max_memory_allocated(d) for d in devices},output_hashes=hashes)
    save(ROOT/(label+'.json'),report)
    print(json.dumps({key:report[key] for key in ['case','implementation','repeat','api_seconds','result_rows']}),flush=True)


def main():
    ROOT.mkdir(parents=True,exist_ok=False);TARGET.mkdir(parents=True,exist_ok=False)
    control=json.loads((PROTOTYPE/'report.json').read_text())
    assert control['source_sha256']==source_identity()
    inputs=[SOURCE/name for name in ['input.pgen','input.pvar','input.psam']]+[PHENO/name for name in ['phenotype.npy','covariates.npy']]
    harness=[Path(__file__),Path(__file__).with_name('selector_pack_control_20260922.py'),
        Path(__file__).with_name('selector_pack_control_20260922.cpp'),
        Path(__file__).with_name('native_host_predicate_gwas_20260922.py'),
        Path(__file__).with_name('direct_bounded_host_selection_prices_20260921.py'),PROTOTYPE/'pack.so']
    assert sha(PROTOTYPE/'pack.so')==control['compiled']['binary_sha256']
    order=[];rng=random.Random(9221649)
    for repeat in range(3):
        cases=list(CASES);rng.shuffle(cases)
        for index,case in enumerate(cases):
            pair=['current','fused'] if (repeat+list(CASES).index(case))%2==0 else ['fused','current']
            order.extend(dict(case=case,implementation=name,repeat=repeat) for name in pair)
    plan=dict(source_sha256=source_identity(),inputs={str(p):sha(p) for p in inputs},
        harness_sha256={str(p):sha(p) for p in harness},order=order,
        shape=dict(samples=2049,markers=16385,traits=512,covariates=2),
        layout=dict(chunk_size=4096,trait_block=256,devices=['cuda:1','cuda:2'],reader_workers=4,prefetch_chunks=2,queue_depth=2),
        mount={str(p):subprocess.check_output(['findmnt','-T',str(p),'-n','-o','SOURCE,FSTYPE,TARGET'],text=True).strip() for p in inputs+[TARGET]},
        preregistration=dict(repeats=3,primary='Complete output-inclusive API time in fresh warmed CUDA processes.',
            secondary='Recorded scan-and-durable-write executor and API process CPU.',
            correctness='Every persisted selected field must hash identically across all six runs per case.',
            cache='Input hashes warm the same files; sequential read/API/read controls remain warm and are not an IO ceiling.',
            inference='Small fixed workload; retain every paired observation and do not fit calculator rates.'),scope=__doc__)
    assert all(v.split()==['/dev/md0','xfs','/data'] for v in plan['mount'].values()),plan['mount']
    save(ROOT/'plan.json',plan);runs=[];reference={}
    for entry in order:
        label=f"{entry['case']}_{entry['repeat']}_{entry['implementation']}"
        with (ROOT/(label+'.log')).open('w') as log:
            subprocess.run([sys.executable,__file__,'--child',entry['case'],entry['implementation'],str(entry['repeat'])],
                stdout=log,stderr=subprocess.STDOUT,check=True)
        row=json.loads((ROOT/(label+'.json')).read_text())
        key=entry['case'];identity=(row['result_rows'],row['output_hashes'])
        if key not in reference:reference[key]=identity
        assert reference[key]==identity,(key,entry)
        runs.append(row);print(json.dumps({k:row[k] for k in ['case','implementation','repeat','api_seconds','result_rows']}),flush=True)
        save(ROOT/'progress.json',dict(completed=len(runs),runs=runs))
    summaries=[]
    for case in CASES:
        rows=[r for r in runs if r['case']==case]
        ratios=[next(r['api_seconds'] for r in rows if r['repeat']==i and r['implementation']=='fused')/
            next(r['api_seconds'] for r in rows if r['repeat']==i and r['implementation']=='current') for i in range(3)]
        summaries.append(dict(case=case,paired_api_ratios=ratios,median_paired_api_ratio=statistics.median(ratios),
            median_api_seconds={name:statistics.median(r['api_seconds'] for r in rows if r['implementation']==name) for name in ['current','fused']}))
    save(ROOT/'report.json',dict(plan=plan,runs=runs,summaries=summaries,all_persisted_fields_exact=True,
        source_sha256=source_identity(),production_changed=False,prices_installed=False,scope=__doc__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--child',nargs=3)
    args=parser.parse_args()
    if args.child:child(args.child[0],args.child[1],int(args.child[2]))
    else:main()
