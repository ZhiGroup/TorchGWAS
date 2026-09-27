"""Separate-process structural reuse audit; synthetic prices, no GWAS timing."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import time

from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.structural_tensor_cache import StructuralTensorWorkCache
from torchgwas.window_model import compare_prepared_windows
from test_jagwas_actual_candidate import actual_candidate,writer_prices
from test_significant_host_model import bank


MODES=('full','significant_empty','significant_sparse','significant_dense','jagwas')


def layout(path,header,b,mode):
    # Exact censuses in this helper are an offline fixture oracle. They are
    # outside the timer and never proposed as productive-job startup work.
    choice=actual_candidate(path,b,2);windows=[];counts=[]
    reduction=None if mode=='full' else 'jagwas' if mode=='jagwas' else 'significant'
    for i,tile in enumerate(choice['tiles']):
        w={key:deepcopy(tile[key]) for key in ('device','trait_range','data','profile')};w['issued_chunks']=4
        begin=2048*i if reduction=='jagwas' else 0;encoded=header.window(begin,begin+512,b)
        w['data'].update(markers=512,encoded=encoded);p=w['profile']
        p['decode_units']={key:1e-8 for row in encoded['chunks'] for key in row['source_units']}
        if reduction!='jagwas':
            w['trait_range']=[512*i,512*(i+1)];p.pop('reduction');p['result_ownership']='borrowed'
            p['result_finish_service']=dict(cpu_seconds=1e-6,serial_cpu_seconds=0.,baseline_copy_bytes=0,
                replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False)
            p['process_units']['bytearray_zero_bytes']=1e-10
            p['writeback_service'].update(submit_seconds=1e-6,wait_seconds=1e-6,fadvise_seconds=1e-6)
        counts.append([0 if mode=='significant_empty' else row['markers']//128 if mode=='significant_sparse'
            else row['markers']*512 if mode=='significant_dense' else row['markers'] for row in encoded['chunks']])
        windows.append(w)
    common=dict(total_traits=512 if reduction=='jagwas' else 1024,reduction=reduction,
        output=dict(block_bytes=None,queue_depth=2,store_beta=reduction!='jagwas',fsync=True),
        shared_capacities=choice['shared_capacities'],endpoint='upper',host_serial_fraction=.5,
        shared_storage_bytes_per_second=1e8,shared_links=[dict(devices=choice['devices'],h2d_bytes_per_second=1e8,d2h_bytes_per_second=1e8)])
    if reduction:common['prices']=writer_prices() if reduction=='jagwas' else bank()
    evidence=None
    if reduction:
        evidence=dict(input_identity=deepcopy(windows[0]['data']['encoded']['input_identity']),
            reduction=reduction,total_traits=common['total_traits'],significance_threshold=None,
            bins=[dict(variant_range=deepcopy(row['variant_range']),trait_range=deepcopy(w['trait_range']),retained=count)
                for w,values in zip(windows,counts) for row,count in zip(w['data']['encoded']['chunks'],values)])
    return dict(windows=windows,partition_axis='variant' if reduction=='jagwas' else 'trait'),common,evidence


def artifacts(directory):
    return {str(p.relative_to(directory)):sha256_file(p) for p in directory.glob('*/*.json')}


def worker(args):
    path=Path(args.input);header=PgenHeaderWork(path);prepared=[];sources=source_identity()
    for mode in MODES:
        before,common,evidence=layout(path,header,128,mode);after,_,_=layout(path,header,512,mode)
        prepared.append((mode,before,after,common,evidence))
    data=Path(args.cache_data);initial=artifacts(data)
    cpu=time.thread_time();wall=time.perf_counter()
    cache=None if args.worker=='control' else StructuralTensorWorkCache(data)
    initialize=dict(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu)
    rows=[]
    for mode,before,after,common,evidence in prepared:
        if cache is None:result=compare_prepared_windows(before,after,survivor_evidence=evidence,**common)
        else:
            with cache.activate():result=compare_prepared_windows(before,after,survivor_evidence=evidence,**common)
        rows.append(dict(mode=mode,comparison=result,cache=None if cache is None else cache.snapshot()))
        print(json.dumps(dict(phase=args.worker,mode=mode,cpu_seconds=result['calculation_cpu_seconds'],
            wall_seconds=result['calculation_wall_seconds'])),flush=True)
    cpu=time.thread_time();wall=time.perf_counter()
    publication=None if cache is None else cache.publish(successful=True)
    publication_cost=dict(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu)
    final=artifacts(data)
    if args.worker=='control':assert not final
    if args.worker=='populate':assert len(publication['stored'])==4 and len(final)==4
    if args.worker=='reuse':
        assert cache.snapshot()['disk_hits']==4 and cache.snapshot()['misses']==0
        assert publication['stored']==[] and final==initial
    assert source_identity()==sources
    result=dict(phase=args.worker,source_sha256=sources,script_sha256=sha256_file(__file__),
        initialization=initialize,publication_cost=publication_cost,publication=publication,
        artifacts_before=initial,artifacts_after=final,rows=rows,cache=None if cache is None else cache.snapshot())
    (Path(args.out)/(args.worker+'.json')).write_text(json.dumps(result,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--cache-data',required=True);parser.add_argument('--worker',choices=['control','populate','reuse'])
    args=parser.parse_args()
    if args.worker:return worker(args)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    data=Path(args.cache_data);data.mkdir(parents=True,exist_ok=False)
    input_hash=sha256_file(args.input);sources=source_identity()
    for phase in ('control','populate','reuse'):
        subprocess.run([sys.executable,__file__,'--input',args.input,'--out',args.out,'--cache-data',args.cache_data,'--worker',phase],check=True)
    reports={phase:json.loads((out/(phase+'.json')).read_text()) for phase in ('control','populate','reuse')}
    def prediction(row):
        value=deepcopy(row['comparison']);value.pop('calculation_wall_seconds');value.pop('calculation_cpu_seconds');return value
    for index in range(len(MODES)):
        expected=prediction(reports['control']['rows'][index])
        assert all(prediction(reports[phase]['rows'][index])==expected for phase in ('populate','reuse'))
    for path in data.glob('*/*.json'):
        record=json.loads(path.read_text());assert record['kind']=='source_work' and record['max_age_seconds'] is None
        assert 'observed_unix_seconds' not in record
        assert 'seconds' not in json.dumps(record['value'])
    assert sha256_file(args.input)==input_hash and source_identity()==sources
    summary=dict(source_sha256=sources,script_sha256=sha256_file(__file__),input_sha256=input_hash,
        helpers={name:sha256_file(name) for name in ['tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py',
            'tests/test_jagwas_scan_work.py','tests/test_significant_host_model.py','tests/fixtures/jagwas_chunk_geometry.json']},
        mounts=dict(input=storage_identity(args.input),cache=storage_identity(data)),reports=reports,
        scope='Three separate model processes, each comparing chunk sizes 128 and 512 over identical prepared windows in five output modes. Captured geometry and real header extents, synthetic independent prices. Exact model equality and immutable source-ledger reuse; no measured GWAS prediction accuracy or speedup. Timed comparison includes disk lookup and fresh pricing/solve, excludes fixture construction, header indexing and process imports. Initialization and publication are separately reported costs.')
    (out/'report.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
