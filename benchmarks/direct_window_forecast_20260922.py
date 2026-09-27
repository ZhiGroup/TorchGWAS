"""Held-out analytical-window audit, not measured GWAS forecast qualification."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import time

from direct_structural_tensor_reuse_20260922 import layout,MODES
from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.planning_session import PlanningWorkCache
from torchgwas.window_model import compare_prepared_windows
from torchgwas.window_forecast import forecast_remaining_windows,remaining_forecast_payback


def extend(template,header,count,mode):
    result=deepcopy(template);bins=[]
    for window in result['windows']:
        start=window['data']['encoded']['variant_range'][0];b=window['profile']['chunk_markers']
        encoded=header.window(start,start+count,b,max_chunks=32);window['data'].update(markers=count,encoded=encoded)
        for row in encoded['chunks']:
            m=row['markers'];keep=0 if mode=='significant_empty' else m//128 if mode=='significant_sparse' else m*512 if mode=='significant_dense' else m
            bins.append(dict(variant_range=deepcopy(row['variant_range']),trait_range=deepcopy(window['trait_range']),retained=keep))
    return result,bins


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True)
    p.add_argument('--horizon-step',type=int,default=512)
    p.add_argument('--create-fixture',action='store_true')
    args=p.parse_args();path=Path(args.input);out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    if args.horizon_step not in (512,1024):raise ValueError('Bounded 512 or 1024 variant horizon step required')
    if args.create_fixture:
        import numpy as np
        from test_pgen_native_reader import write_pgen
        if path.exists():raise ValueError('Fixture creation requires a new path')
        path.parent.mkdir(parents=True,exist_ok=True)
        write_pgen(path,np.tile(np.arange(2049,dtype=np.uint32)%3,(8*args.horizon_step+1,1)).astype(np.uint8))
    source=source_identity();input_hash=sha256_file(path);header=PgenHeaderWork(path);records=[]
    for mode in MODES:
        baseline,common,evidence=layout(path,header,128,mode);alternative,_,_=layout(path,header,512,mode)
        if mode=='jagwas':
            for template in (baseline,alternative):
                for i,window in enumerate(template['windows']):
                    start=4*args.horizon_step*i
                    window['data']['encoded']=header.window(start,start+512,window['profile']['chunk_markers'])
        if mode=='full':common['output']['block_bytes']=1<<20
        common.update(max_source_chunks=64,model_identity=dict(source_sha256=source,protocol='synthetic-independent-prices-captured-geometry'))
        series=[];cache=PlanningWorkCache();started=time.perf_counter();cpu=time.thread_time()
        with cache.activate():
            for count in (args.horizon_step*i for i in (1,2,3,4)):
                before,bins=extend(baseline,header,count,mode);after,_=extend(alternative,header,count,mode)
                occupancy=None if evidence is None else dict(evidence,bins=bins)
                series.append(compare_prepared_windows(before,after,survivor_evidence=occupancy,**common))
        wall=time.perf_counter()-started;cpu=time.thread_time()-cpu
        units=4*args.horizon_step*1024  # same pair coverage for trait and JAGWAS variant layouts
        forecast=forecast_remaining_windows(series[:3],remaining_pairs=units,
            boundary_adjustments=dict(baseline=[0.,0.],candidate=[0.,0.]),relative_model_error=.05,
            max_slope_change=.1,max_extrapolation=2.,assumptions=dict(
                source_work='Declared hardcall fixture header work continues over this exact held-out extent.',
                output_occupancy='Same deliberately specified synthetic empty/sparse/dense survivor scenario.',
                resource_capacity='Same synthetic independent prices and declared capacities.',
                partition_balance='Both device windows contain the same marker count at all horizons.'))
        heldout={}
        for name in ('baseline','candidate'):
            expected=series[-1][name]['estimated_window_seconds'];f=forecast['forecasts'][name]
            heldout[name]=dict(modeled_seconds=expected,inside_scenario=f['lower_seconds']<=expected<=f['upper_seconds'],
                center_relative_error=((f['lower_seconds']+f['upper_seconds'])/2-expected)/expected)
        payback=remaining_forecast_payback(forecast,planning_seconds=sum(r['calculation_wall_seconds'] for r in series[:3]),
            switching_seconds=0.,publication_seconds=.06,reserve_seconds=.01)
        records.append(dict(mode=mode,windows=series,forecast=forecast,heldout=heldout,payback=payback,
            four_window_construction_cpu_seconds=cpu,four_window_construction_wall_seconds=wall,cache=cache.snapshot()))
        print(json.dumps(dict(mode=mode,status=forecast['status'],heldout=heldout,payback=payback)),flush=True)
        cache.close()
    assert source_identity()==source and sha256_file(path)==input_hash
    report=dict(source_sha256=source,script_sha256=sha256_file(__file__),input_sha256=input_hash,mount=storage_identity(path),
        horizon_step=args.horizon_step,
        helpers={name:sha256_file(name) for name in ['benchmarks/direct_structural_tensor_reuse_20260922.py',
            'tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py','tests/test_jagwas_scan_work.py',
            'tests/test_significant_host_model.py','tests/fixtures/jagwas_chunk_geometry.json']},records=records,
        scope='Three bounded analytical horizons predict a fourth held-out analytical horizon, using real header extents, captured geometry and synthetic prices. No measured GWAS prediction accuracy, future-source statistical guarantee or automatic tuning qualification. Publication/reserve costs in payback are explicit scenarios. Offline fixture profile construction includes exact census work outside timed calculations.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
