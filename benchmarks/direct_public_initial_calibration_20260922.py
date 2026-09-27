"""Fresh-process productive cache reuse and unchanged output on local storage.

This is a lifecycle/correctness audit, not a calibrated speedup experiment.
API wall time includes context binding; Python import time is excluded.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
from unittest.mock import patch
import numpy as np
import torch
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.geometry_collection import write_record
from torchgwas.sumstats import open_binary_sumstats,BinarySumstatsWriter,DenseWriteProgress
from direct_jagwas_bounded_execution_20260922 import read_result as read_jagwas
from direct_significant_bounded_execution_20260922 import read_result as read_significant


MODES=['full','full_tiles','full_variants','jagwas','significant_tiles']


def worker(args):
    fixture=Path(args.fixture);out=Path(args.out);data=Path(args.output_data)
    mode,label=args.worker.split(':');directory=data/(mode+'_'+label)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    kwargs=dict(pgen_mode='hardcall',chunk_size=128,device='cuda:0',reader_workers=4,
        prefetch_chunks=2,sumstats_queue_depth=2,sumstats_fsync=True,output_dir=directory)
    if mode.endswith('_tiles'):kwargs.update(trait_block=128,trait_devices=['cuda:0','cuda:1'])
    if mode=='jagwas':kwargs.update(reduce='jagwas',variant_devices=['cuda:0','cuda:1'])
    if mode=='full_variants':kwargs.update(variant_devices=['cuda:0','cuda:1'])
    if mode=='significant_tiles':kwargs.update(reduce='significant',significance_threshold=.01)
    if label!='control':
        kwargs['initial_calibration']=dict(cache_dir=str(data/'cache'),max_age_seconds=3600.,
            max_chunks_per_device=3,validation_chunks_per_device=2,warmup_chunks=1,stride=1,
            max_window_seconds=10.)
    events=[]
    import torchgwas.sumstats_indexed as indexed
    import torchgwas.api as api
    original=indexed.write_indexed_sumstats
    original_dense_init=BinarySumstatsWriter.__post_init__
    original_json=api.write_json
    fail=label in ('failed_writer','failed_metadata')
    before={str(p):sha256_file(p) for p in (data/'cache').glob('**/*.json')} if fail else None
    def written(event):
        events.append(event)
        if label=='failed_writer' and len(events)==(3 if mode.startswith('full') else 12):
            raise OSError('injected writer failure')
    def indexed_writer(*a,**kw):
        callback=kw.get('on_chunk_written')
        def deliver(event):
            written(event)
            if callback is not None:callback(event)
        kw['on_chunk_written']=deliver
        return original(*a,**kw)
    def dense_init(writer):
        callback=writer.on_write_progress
        def deliver(event):
            written(event)
            if callback is not None:callback(event)
        writer.on_write_progress=deliver
        return original_dense_init(writer)
    def json_writer(value,path):
        if label=='failed_metadata' and Path(path).name=='qc.json':raise OSError('injected metadata failure')
        return original_json(value,path)
    started=time.perf_counter()
    try:
        with patch.object(indexed,'write_indexed_sumstats',indexed_writer),patch.object(api,'write_json',json_writer),\
             patch.object(BinarySumstatsWriter,'__post_init__',dense_init):
            result=run_linear_gwas(fixture/'input.pgen',fixture/'phenotype.npy',fixture/'covariates.npy',**kwargs)
    except Exception as error:
        cause=error
        while cause is not None and not str(cause).startswith('injected '):cause=cause.__cause__
        if not fail or cause is None:raise
        import threading
        assert not [t.name for t in threading.enumerate() if t.name.startswith('torchgwas-')]
        assert before=={str(p):sha256_file(p) for p in (data/'cache').glob('**/*.json')}
        write_record(out/(mode+'_'+label+'.json'),dict(failure=str(error),
            original_records_unchanged=True,new_records_published=0,completed_parts=len(events),
            source_sha256=source_identity(),script_sha256=sha256_file(__file__)))
        return
    assert not fail,'Expected injected failure was not reached'
    elapsed=time.perf_counter()-started
    audit=result.run_metadata.get('initial_calibration')
    if audit:
        assert audit['successful'] and audit['publication_inputs_unchanged']
        assert not audit['unobserved_devices'] and not audit['cache_errors']
        assert len(audit['windows'])==(1 if mode=='full' else 2)
        for device,window in audit['windows'].items():
            assert not window['pending'] and 2<=len(window['observations'])<=3
            decision=window['refresh_decisions'][device]
            assert decision['cache_hit']==(label=='second'),decision
            assert not window['publication']['incomplete_devices']
        persisted=json.loads((directory/'run.json').read_text())
        assert persisted['initial_calibration']==audit
        progress=audit['output_progress']
        assert progress['events']==len(events) and progress['first_written_seconds'] is not None
        assert 0<=progress['first_written_seconds']<=progress['elapsed_to_report_seconds']<=elapsed
        if mode.startswith('full'):
            assert progress['dense_statistic_cells']==4097*512
            assert progress['dense_statistic_bytes']==4097*512*8
            assert progress['first_fsynced_part_seconds'] is None
        else:
            assert progress['first_fsynced_part_seconds'] is not None
    material=[e.completed for e in events if not isinstance(e,DenseWriteProgress) and e.part_file_fsynced]
    dense=[e.completed for e in events if isinstance(e,DenseWriteProgress)]
    row=dict(mode=mode,label=label,output=str(directory),api_wall_seconds=elapsed,
        first_fsynced_part_seconds=min(material)-started if material else None,
        first_dense_statistic_seconds=min(dense)-started if dense else None,
        first_output_progress_seconds=min(e.completed for e in events)-started if events else None,
        first_processed_chunk_seconds=min(e.completed for e in events if not isinstance(e,DenseWriteProgress))-started if material else None,
        indexed_events=len(events)-len(dense),dense_events=len(dense),rows=result.run_metadata['n_result_rows'],calibration=audit,
        source_sha256=source_identity())
    write_record(out/(mode+'_'+label+'.json'),row)


def arrays(directory,mode):
    if mode=='jagwas':return {'chi2':read_jagwas(directory,4097)[1]}
    if mode=='significant_tiles':return read_significant(directory,4097,512)[1]
    beta,tstat,manifest=open_binary_sumstats(directory/'sumstats')
    assert manifest['shape']==[4097,512]
    # Use bounded row slicing, including for lazy tiled readers.
    return {name:np.concatenate([np.asarray(value[i:i+128,:]) for i in range(0,4097,128)])
            for name,value in [('beta',beta),('t_stat',tstat)]}


def main(args):
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    data=Path(args.output_data);data.mkdir(parents=True,exist_ok=False)
    fixture=Path(args.fixture);sources=source_identity()
    inputs={name:sha256_file(fixture/name) for name in ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}
    mounts={key:storage_identity(path) for key,path in [('input',fixture),('output',data)]}
    assert all(row['fstype']=='xfs' and row['source']=='/dev/md0' for row in mounts.values()),mounts
    rows=[];checks={}
    for mode in MODES:
        records={}
        for label in ['control','first','second']:
            subprocess.run([sys.executable,__file__,'--fixture',str(fixture),'--out',str(out),
                '--output-data',str(data),'--worker',mode+':'+label],check=True)
            row=json.loads((out/(mode+'_'+label+'.json')).read_text())
            assert row['source_sha256']==sources
            rows.append(row)
            if label=='first':
                for device,window in row['calibration']['windows'].items():
                    saved=window['publication']['publications'][device]
                    records[device]=(saved,Path(saved['path']).read_bytes())
            if label=='second':
                for device,(saved,old) in records.items():
                    assert Path(saved['path']).read_bytes()==old
                    decision=row['calibration']['windows'][device]['refresh_decisions'][device]
                    assert decision['previous_record_sha256']==saved['record_sha256']
            print(json.dumps(dict(mode=mode,label=label,seconds=row['api_wall_seconds'],
                first_part_seconds=row['first_fsynced_part_seconds'])),flush=True)
        reference=arrays(data/(mode+'_control'),mode);deltas={}
        for label in ['first','second']:
            values=arrays(data/(mode+'_'+label),mode)
            assert values.keys()==reference.keys()
            for name,actual in values.items():
                if name in ('variant_index','trait_index','df'):np.testing.assert_array_equal(actual,reference[name])
                else:np.testing.assert_allclose(actual,reference[name],rtol=6e-5,atol=3e-4)
            deltas[label]={name:float(np.max(np.abs(actual-reference[name]))) for name,actual in values.items()}
        checks[mode]=dict(output_deltas=deltas,original_records_unchanged=True,
            baseline_observation_times={d:json.loads(old)['observed_unix_seconds'] for d,(_,old) in records.items()})
    failures={}
    for mode,label in [('jagwas','failed_writer'),('jagwas','failed_metadata'),('full','failed_writer')]:
        subprocess.run([sys.executable,__file__,'--fixture',str(fixture),'--out',str(out),
            '--output-data',str(data),'--worker',mode+':'+label],check=True)
        failures[mode+':'+label]=json.loads((out/(mode+'_'+label+'.json')).read_text())
    assert source_identity()==sources
    assert all(sha256_file(fixture/name)==value for name,value in inputs.items())
    write_record(out/'report.json',dict(source_sha256=sources,inputs=inputs,mounts=mounts,
        script_sha256=sha256_file(__file__),rows=rows,checks=checks,failures=failures,
        scope='Five output/partition modes, three fresh processes each. Both control and measured arms retain the same writer telemetry. Productive measurements only; cache validation may expand on drift. Timings include API setup and identity hashing but exclude imports. One run per condition cannot establish overhead or speedup. Dense first beta/t writes may precede df flushing and fsync; reduced output reports first part fsync. No automatic selection or hardware capacity inference.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--output-data',required=True);parser.add_argument('--worker')
    args=parser.parse_args()
    worker(args) if args.worker else main(args)
