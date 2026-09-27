"""Output/resource composition on real headers and captured GPU geometry.

All component prices are synthetic test controls, not a calibrated forecast.
No association run is timed. Real writer payload extents and model construction
cost are checked without claiming selection accuracy or end-to-end benefit.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import time
import numpy as np

from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.planning_session import PlanningWorkCache
from torchgwas.window_model import prepared_window_runtime,compare_prepared_windows
from test_jagwas_actual_candidate import actual_candidate,writer_prices
from test_significant_host_model import bank


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True);p.add_argument('--output-data',required=True)
    args=p.parse_args();path=Path(args.input);out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    data_out=Path(args.output_data);data_out.mkdir(parents=True,exist_ok=False)
    sources=source_identity();input_hash=sha256_file(path);header=PgenHeaderWork(path);rows=[];cache=PlanningWorkCache();layouts={};comparisons=[]
    with cache.activate():
        for b in (128,256,512):
            # This test helper's exact census is an offline fixture operation.
            # It is not a proposed production-startup action.
            choice=actual_candidate(path,b,2)
            for mode in ('full','significant_empty','significant_sparse','significant_dense','jagwas'):
                reduction=None if mode=='full' else 'jagwas' if mode=='jagwas' else 'significant'
                windows=[];retained=[];artifact_bytes=0
                output=dict(block_bytes=None,queue_depth=2,store_beta=mode!='jagwas',fsync=True)
                for i,tile in enumerate(choice['tiles']):
                    w={key:deepcopy(tile[key]) for key in ('device','trait_range','data','profile')};w['issued_chunks']=4
                    begin=2048*i if mode=='jagwas' else 0
                    window=header.window(begin,begin+512,b);w['data'].update(markers=512,encoded=window)
                    profile=w['profile'];profile['decode_units']={key:1e-8 for row in window['chunks'] for key in row['source_units']}
                    if reduction!='jagwas':
                        w['trait_range']=[i*512,(i+1)*512]
                        profile.pop('reduction');profile['result_ownership']='borrowed'
                        profile['result_finish_service']=dict(cpu_seconds=1e-6,serial_cpu_seconds=0.,baseline_copy_bytes=0,
                            replaces_fixed_finish_and_tensor_conversion=True,includes_ready_cuda_event=False)
                        profile['process_units']['bytearray_zero_bytes']=1e-10
                        profile['writeback_service'].update(submit_seconds=1e-6,wait_seconds=1e-6,fadvise_seconds=1e-6)
                    counts=[]
                    for j,row in enumerate(window['chunks']):
                        m=row['markers']
                        keep=0 if mode=='significant_empty' else m//128 if mode=='significant_sparse' else m*512 if mode=='significant_dense' else m
                        counts.append(keep)
                        if reduction is not None and keep:
                            if reduction=='jagwas':values=dict(variant_index=np.zeros(keep,dtype=np.int64),chi2=np.zeros(keep,dtype=np.float64))
                            else:values=dict(variant_index=np.zeros(keep,dtype=np.int64),trait_index=np.zeros(keep,dtype=np.int64),
                                beta=np.zeros(keep,dtype=np.float32),t_stat=np.zeros(keep,dtype=np.float32),df=np.zeros(keep,dtype=np.float32))
                            target=data_out/f'{mode}_{b}_{i}_{j}.npz'
                            np.savez(target,**values);artifact_bytes+=target.stat().st_size
                    windows.append(w);retained.append(counts)
                kwargs=dict(total_traits=512 if reduction=='jagwas' else 1024,reduction=reduction,output=output,
                    shared_capacities=choice['shared_capacities'],endpoint='upper',host_serial_fraction=.5,
                    shared_storage_bytes_per_second=1e8,shared_links=[dict(devices=choice['devices'],h2d_bytes_per_second=1e8,d2h_bytes_per_second=1e8)])
                if reduction:kwargs.update(prices=writer_prices() if reduction=='jagwas' else bank(),retained=retained)
                cpu=time.thread_time();wall=time.perf_counter()
                report=prepared_window_runtime(windows,**kwargs)
                wall=time.perf_counter()-wall;cpu=time.thread_time()-cpu
                if reduction:assert artifact_bytes==report['payload_bytes'],(mode,artifact_bytes,report['payload_bytes'])
                else:assert report['payload_bytes']==2*512*(8*512+4)
                rows.append(dict(chunk_size=b,mode=mode,planning_wall_seconds=wall,planning_cpu_seconds=cpu,
                    artifact_payload_bytes=artifact_bytes if reduction else None,model=report,cache=cache.snapshot()))
                layouts[b,mode]=(windows,kwargs,report)
                print(json.dumps(dict(chunk_size=b,mode=mode,wall_seconds=wall,cpu_seconds=cpu,nodes=report['graph_nodes'],payload=report['payload_bytes'])),flush=True)
    for mode in ('full','significant_empty','significant_sparse','significant_dense','jagwas'):
        before,kwargs,expected_before=layouts[128,mode];after,_,expected_after=layouts[512,mode]
        common=dict(kwargs);counts=common.pop('retained',None);reduction=common['reduction'];evidence=None
        if reduction:
            evidence=dict(input_identity=deepcopy(before[0]['data']['encoded']['input_identity']),reduction=reduction,
                total_traits=common['total_traits'],significance_threshold=common.get('significance_threshold'),
                bins=[dict(variant_range=deepcopy(chunk['variant_range']),trait_range=deepcopy(w['trait_range']),retained=count)
                    for w,values in zip(before,counts) for chunk,count in zip(w['data']['encoded']['chunks'],values)])
        axis='variant' if reduction=='jagwas' else 'trait';comparison_cache=PlanningWorkCache()
        with comparison_cache.activate():
            for repeat in range(2):
                result=compare_prepared_windows(dict(windows=before,partition_axis=axis),dict(windows=after,partition_axis=axis),
                    survivor_evidence=evidence,**common)
                assert result['baseline']['payload_bytes']==expected_before['payload_bytes']
                assert result['candidate']['payload_bytes']==expected_after['payload_bytes']
                assert result['baseline']['estimated_window_seconds']==expected_before['estimated_window_seconds']
                assert result['candidate']['estimated_window_seconds']==expected_after['estimated_window_seconds']
                comparisons.append(dict(mode=mode,repeat=repeat,model=result,cache=comparison_cache.snapshot()))
                print(json.dumps(dict(comparison=mode,repeat=repeat,wall_seconds=result['calculation_wall_seconds'],
                    cpu_seconds=result['calculation_cpu_seconds'])),flush=True)
        comparison_cache.close()
    assert source_identity()==sources and sha256_file(path)==input_hash
    report=dict(source_sha256=sources,script_sha256=sha256_file(__file__),input_sha256=input_hash,
        helpers={name:sha256_file(name) for name in ['tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py',
            'tests/test_jagwas_scan_work.py','tests/test_significant_host_model.py','tests/fixtures/jagwas_chunk_geometry.json']},
        mounts=dict(input=storage_identity(path),output=storage_identity(data_out)),rows=rows,comparisons=comparisons,
        scope='15 finite two-GPU prepared-window model scenarios and 10 equivalent-work comparisons on real PGEN header extents and captured shape geometry. Prices are synthetic. Part bytes compared with actual np.savez files; timing is model construction plus solution, not GWAS runtime. No measured prediction or autotuning qualification.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
