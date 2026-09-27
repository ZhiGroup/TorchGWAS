"""Matched calculator-only profiling after a saved structural reuse audit."""
import argparse
import cProfile
import io
import json
from pathlib import Path
import pstats
import time

from direct_structural_tensor_reuse_20260922 import layout
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.structural_tensor_cache import StructuralTensorWorkCache
from torchgwas.window_model import compare_prepared_windows


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--cache-data',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    path=Path(args.input);header=PgenHeaderWork(path);before,common,evidence=layout(path,header,128,'full');after,_,_=layout(path,header,512,'full')
    rows=[];expected=None
    for repeat in range(4):
        for reuse in ([False,True] if repeat%2==0 else [True,False]):
            wall=time.perf_counter();cpu=time.thread_time()
            cache=StructuralTensorWorkCache(args.cache_data) if reuse else None
            if cache is None:result=compare_prepared_windows(before,after,survivor_evidence=evidence,**common)
            else:
                with cache.activate():result=compare_prepared_windows(before,after,survivor_evidence=evidence,**common)
            row=dict(repeat=repeat,reuse=reuse,cpu_seconds=time.thread_time()-cpu,wall_seconds=time.perf_counter()-wall,
                cache=None if cache is None else cache.snapshot())
            predicted=(result['baseline']['estimated_window_seconds'],result['candidate']['estimated_window_seconds'])
            if expected is not None:assert predicted==expected
            expected=predicted
            if cache is not None:assert cache.snapshot()['disk_hits']==2 and cache.snapshot()['misses']==0
            rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='cache'}),flush=True)
    profiler=cProfile.Profile();profiler.enable()
    cache=StructuralTensorWorkCache(args.cache_data)
    with cache.activate():compare_prepared_windows(before,after,survivor_evidence=evidence,**common)
    profiler.disable();profiler.dump_stats(str(out/'reuse.prof'))
    text=io.StringIO();stats=pstats.Stats(profiler,stream=text).strip_dirs().sort_stats('cumulative')
    stats.print_stats(35);stats.print_stats('structural_tensor_cache');stats.print_stats('calibration_cache')
    (out/'profile.txt').write_text(text.getvalue());(out/'timings.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(text.getvalue(),flush=True)


if __name__=='__main__':main()
