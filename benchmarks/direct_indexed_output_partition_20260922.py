"""Real reduced-output identity and cache audit over a nonzero source range."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import sys
import threading
import time
from unittest.mock import patch

import numpy as np
import torch
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.sumstats_indexed import open_indexed_sumstats
import torchgwas.sumstats_indexed as indexed


MODES=('significant_tiles','significant_empty','significant_full','jagwas')
FIRST,LAST,K=129,1154,512


def artifacts(path):return {str(p):sha256_file(p) for p in path.glob('**/*.json')}


def read_arrays(directory,mode):
    manifest,parts=open_indexed_sumstats(directory/'sumstats');parts=list(parts)
    fields=['variant_index','chi2'] if mode=='jagwas' else ['variant_index','trait_index','beta','t_stat','df']
    values={name:np.concatenate([part[name] for part in parts]) if parts else np.empty(0) for name in fields}
    vi=values['variant_index'];key=vi if mode=='jagwas' else vi*K+values['trait_index']
    assert len(np.unique(key))==len(key) and np.all((0<=vi)&(vi<LAST-FIRST))
    order=np.argsort(key)
    return {name:value[order] for name,value in values.items()}


def worker(args):
    mode,label=args.worker.split(':');fixture=Path(args.fixture);data=Path(args.output_data)
    directory=data/(mode+'_'+label);cache=data/'cache';events=[]
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    settings=dict(pgen_mode='hardcall',compute_dtype='float32',chunk_size=128,device='cuda:0',reader_workers=4,
        prefetch_chunks=2,sumstats_queue_depth=2,sumstats_fsync=True,output_dir=directory,variant_range=(FIRST,LAST))
    if mode=='jagwas':settings.update(reduce='jagwas',variant_devices=['cuda:0','cuda:1'])
    else:
        settings.update(reduce='significant',significance_threshold=1e-300 if mode=='significant_empty' else .01)
        if mode!='significant_full':settings.update(trait_block=128,trait_devices=['cuda:0','cuda:1'])
    if label!='control':
        settings['initial_calibration']=dict(cache_dir=str(cache),max_age_seconds=3600.,
            max_chunks_per_device=3,validation_chunks_per_device=2,warmup_chunks=1,stride=1,max_window_seconds=10.)
    before=artifacts(cache);original=indexed.write_indexed_sumstats
    def wrap(*a,**kw):
        callback=kw.get('on_chunk_written')
        def observe(event):
            events.append(event)
            if label=='failed' and len(events)==2:raise OSError('injected indexed callback failure')
            if callback is not None:callback(event)
        kw['on_chunk_written']=observe
        return original(*a,**kw)
    began=time.perf_counter()
    try:
        with patch.object(indexed,'write_indexed_sumstats',wrap):
            result=run_linear_gwas(fixture/'input.pgen',fixture/'phenotype.npy',fixture/'covariates.npy',**settings)
    except Exception as error:
        cause=error
        while cause is not None and not str(cause).startswith('injected '):cause=cause.__cause__
        if label!='failed' or cause is None:raise
        assert artifacts(cache)==before
        assert not [t.name for t in threading.enumerate() if t.name.startswith('torchgwas-')]
        report=dict(mode=mode,label=label,failure=str(error),source_sha256=source_identity(),
            records_unchanged=True,worker_threads_closed=True,events=len(events))
    else:
        assert label!='failed'
        elapsed=time.perf_counter()-began
        audit=result.run_metadata.get('initial_calibration');groups={}
        if audit:
            digests=audit['binding_digests']
            assert digests['publication']['status']=='checked' and not digests['state']['errors']
            if label=='second':
                assert digests['state']['disk_hits']==3 and digests['state']['hashed_files']==0
                assert not digests['publication']['stored']
            counts=audit['output_occupancy']
            assert counts['sample_status']=='sampled' and counts['unbound_events']==0
            assert len(counts['bins'])==(3 if mode=='significant_full' else 6)
            assert counts['status']==('published' if label=='first' else 'reused_without_renewal'),counts
            assert counts['cache_hit']==(label=='second')
            for event in events:
                part=event.partition;assert part is not None
                assert event.source_variant_range==(FIRST+event.start,FIRST+event.end)
                key=(part.device,part.variant_range,part.trait_range)
                groups.setdefault(key,[]).append(event.source_variant_range)
                if event.rows:
                    with np.load(directory/'sumstats'/event.part_file) as values:
                        assert len(values['variant_index'])==event.rows
                        assert np.all((event.start<=values['variant_index'])&(values['variant_index']<event.end))
                        if mode!='jagwas':assert np.all((part.trait_range[0]<=values['trait_index'])&(values['trait_index']<part.trait_range[1]))
                else:assert event.part_file is None and event.part_bytes==0
            for (device,span,traits),pieces in groups.items():
                pieces.sort();assert pieces[0][0]==span[0] and pieces[-1][1]==span[1]
                assert all(a[1]==b[0] for a,b in zip(pieces,pieces[1:]))
                if mode=='jagwas':assert traits==(0,K)
            expected_groups=2 if mode=='jagwas' else 1 if mode=='significant_full' else 4
            assert len(groups)==expected_groups
            assert sum((hi-lo)*(b-a) for _,(lo,hi),(a,b) in groups)==(LAST-FIRST)*K
            if label=='second':assert all(artifacts(cache).get(path)==digest for path,digest in before.items())
            assert json.loads((directory/'run.json').read_text())['initial_calibration']==audit
        if mode=='significant_empty':assert any(event.rows==0 for event in events)
        report=dict(mode=mode,label=label,source_sha256=source_identity(),api_wall_seconds=elapsed,
            first_output_seconds=min(event.completed for event in events)-began,
            rows=result.run_metadata['n_result_rows'],events=[asdict(event) for event in events],calibration=audit)
    (Path(args.out)/(mode+'_'+label+'.json')).write_text(json.dumps(report,indent=2)+'\n')


def main(args):
    fixture=Path(args.fixture);out=Path(args.out);data=Path(args.output_data)
    out.mkdir(parents=True,exist_ok=False);data.mkdir(parents=True,exist_ok=False)
    source=source_identity();inputs={name:sha256_file(fixture/name) for name in
        ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}
    mounts={name:storage_identity(path) for name,path in [('input',fixture),('output',data)]}
    assert all(row['fstype']=='xfs' and row['source']=='/dev/md0' for row in mounts.values())
    reports=[];checks={}
    for mode in MODES:
        for label in ('control','first','second'):
            subprocess.run([sys.executable,__file__,'--fixture',str(fixture),'--out',str(out),
                '--output-data',str(data),'--worker',mode+':'+label],check=True)
            row=json.loads((out/(mode+'_'+label+'.json')).read_text());assert row['source_sha256']==source
            reports.append(row)
            print(json.dumps(dict(mode=mode,label=label,seconds=row['api_wall_seconds'],rows=row['rows'])),flush=True)
        reference=read_arrays(data/(mode+'_control'),mode)
        for label in ('first','second'):
            actual=read_arrays(data/(mode+'_'+label),mode)
            for name,values in actual.items():
                if name in ('variant_index','trait_index','df'):np.testing.assert_array_equal(values,reference[name])
                else:np.testing.assert_allclose(values,reference[name],rtol=6e-5,atol=3e-4)
        checks[mode]=dict(output_equal=True,rows=len(reference['variant_index']))
    subprocess.run([sys.executable,__file__,'--fixture',str(fixture),'--out',str(out),
        '--output-data',str(data),'--worker','significant_tiles:failed'],check=True)
    failure=json.loads((out/'significant_tiles_failed.json').read_text())
    assert source_identity()==source and all(sha256_file(fixture/name)==digest for name,digest in inputs.items())
    (out/'report.json').write_text(json.dumps(dict(source_sha256=source,script_sha256=sha256_file(__file__),
        inputs=inputs,mounts=mounts,source_variant_range=[FIRST,LAST],records=reports,checks=checks,failure=failure,
        scope='Twelve fresh-process public API runs plus one failure case. Nonzero source range, empty and sparse significant output, multiple phenotype tiles and variant-sharded full-panel JAGWAS. Immutable output-count evidence, unchanged numerical results and cleanup; no model-accuracy or speedup claim. API times include binding/cache work and output, exclude imports. Control also records writer events but does not collect bound survivor evidence.'),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',required=True);p.add_argument('--out',required=True)
    p.add_argument('--output-data',required=True);p.add_argument('--worker');args=p.parse_args()
    worker(args) if args.worker else main(args)
