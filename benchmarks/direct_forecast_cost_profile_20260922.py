"""Replay a verified productive frontier to attribute analytical planning cost."""
import argparse
import cProfile
from copy import deepcopy
import io
import json
from pathlib import Path
import pstats
import time

import torch
from test_jagwas_actual_candidate import writer_prices
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.planning_session import PlanningWorkCache
from torchgwas.productive_forecast import productive_window_proposal
from torchgwas.tensor_work import eager_statistics_work
from torchgwas.window_model import compare_prepared_windows,prepared_source_window


def replay(prior,cache,*,reuse=False):
    start=prior['admission'];state=prior['records'][1]['model_snapshot'];candidate=start['candidate']
    path=next(p for p in prior['inputs'] if p.endswith('/input.pgen'))
    stages=[]
    def timed(name,callback):
        wall=time.perf_counter();cpu=time.thread_time();value=callback()
        stages.append(dict(name=name,wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu))
        return value
    header=timed('header',lambda:PgenHeaderWork(path,max_cached_signatures=1024 if reuse else 0));comparisons=[]
    model=prior['records'][1]['comparisons'][0]['comparison_contract']['model_identity']
    common=dict(total_traits=512,reduction='jagwas',output=candidate['output'],shared_capacities=candidate['shared_capacities'],
        endpoint='upper',host_serial_fraction=.5,prices=writer_prices(),max_source_chunks=64,model_identity=model)
    with cache.activate():
        for count in (256,512,768):
            layouts=[];bins=[];wall=time.perf_counter();cpu=time.thread_time()
            for size in (state['current_chunk_size'],256):
                windows=[]
                for row,tile in zip(state['partitions'],candidate['tiles']):
                    lo=row['cursor']
                    if reuse:
                        w=prepared_source_window(tile,header,start=lo,stop=lo+count,chunk_markers=size,
                            issued_chunks=row['issued_chunks'],expected_input_identity=start['input_file_identity'],max_chunks=16)
                    else:
                        w={key:deepcopy(tile[key]) for key in ('device','trait_range','data','profile')}
                        w['issued_chunks']=row['issued_chunks'];w['profile']['chunk_markers']=size
                        w['data'].update(markers=count,encoded=header.window(lo,lo+count,size,max_chunks=16))
                    windows.append(w)
                    if len(layouts)==0:bins.append(dict(variant_range=[lo,lo+count],trait_range=[0,512],retained=count))
                layouts.append(dict(windows=windows,partition_axis='variant'))
            evidence=dict(input_identity=header.input_identity,reduction='jagwas',total_traits=512,significance_threshold=None,bins=bins)
            stages.append(dict(name='windows:'+str(count),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu))
            comparisons.append(timed('compare:'+str(count),lambda:compare_prepared_windows(*layouts,survivor_evidence=evidence,**common)))
    old=prior['records'][1]['steps'][0]['forecast_audit']
    proposal,audit=timed('forecast',lambda:productive_window_proposal(state,comparisons,chunk_sizes=[128,256],
        boundary_adjustments=dict(baseline=[0.,1.],candidate=[0.,1.]),relative_model_error=.05,
        max_slope_change=.1,max_extrapolation=8.,assumptions=old['assumptions']))
    assert audit==old and proposal==prior['records'][1]['steps'][0]['value']
    for current,expected in zip(comparisons,prior['records'][1]['comparisons']):
        for name in ('baseline','candidate'):
            assert current[name]==expected[name],name
    return dict(stages=stages,cache=cache.snapshot(),header_cache=header.cache_info(),proposal=proposal,audit=audit)


def main():
    p=argparse.ArgumentParser();p.add_argument('--execution',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    previous=Path(args.execution)/'report.json';prior=json.loads(previous.read_text());source=source_identity()
    assert all(sha256_file(path)==digest for path,digest in prior['inputs'].items())
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    # Match the prior audit's source-meta initialization before its live callback.
    eager_statistics_work(2049,128,512,2,True)
    rows=[]
    for phase in ('first','repeat1','repeat2'):
        cache=PlanningWorkCache();wall=time.perf_counter();cpu=time.thread_time()
        result=replay(prior,cache)
        result.update(phase=phase,wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu)
        cache.close();rows.append(result)
        print(json.dumps({key:result[key] for key in ('phase','wall_seconds','cpu_seconds','stages')}),flush=True)
    profiler=cProfile.Profile();cache=PlanningWorkCache();profiler.enable();replay(prior,cache);profiler.disable();cache.close()
    profiler.dump_stats(out/'planning.prof');stream=io.StringIO()
    stats=pstats.Stats(profiler,stream=stream).sort_stats('cumulative');stats.print_stats(55)
    (out/'profile.txt').write_text(stream.getvalue())
    top=[]
    for (filename,line,function),(primitive,total,self_time,cumulative,_) in stats.stats.items():
        top.append(dict(file=filename,line=line,function=function,calls=total,self_seconds=self_time,cumulative_seconds=cumulative))
    top.sort(key=lambda row:row['cumulative_seconds'],reverse=True)
    assert source_identity()==source and all(sha256_file(path)==digest for path,digest in prior['inputs'].items())
    report=dict(source_sha256=source,script_sha256=sha256_file(__file__),prior_report_sha256=sha256_file(previous),
        records=rows,profile_top=top[:80],
        scope='Calculator-only replay of an actual held productive source frontier. Exact model reports/proposal/audit equality to the prior run. Fresh planning cache per repetition, source-meta initialization before timings as in admission. Includes header/window construction and all three comparisons. cProfile attribution is instrumented and not a timing benchmark. Synthetic component prices; no association timing, prediction qualification or throughput claim.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
