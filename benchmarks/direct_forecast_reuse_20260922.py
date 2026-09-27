"""Matched calculator replay: bounded signature reuse and source-window copies."""
import argparse
import cProfile
import gc
import io
import json
from pathlib import Path
import pstats
import statistics
import time

import torch
from direct_forecast_cost_profile_20260922 import replay
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.planning_session import PlanningWorkCache
from torchgwas.tensor_work import eager_statistics_work


def main():
    p=argparse.ArgumentParser();p.add_argument('--execution',required=True);p.add_argument('--out',required=True)
    p.add_argument('--pairs',type=int,default=12)
    args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    previous=Path(args.execution)/'report.json';prior=json.loads(previous.read_text());source=source_identity()
    helpers={str(Path(__file__).with_name(name)):sha256_file(Path(__file__).with_name(name)) for name in
        ('direct_forecast_cost_profile_20260922.py',)}
    assert all(sha256_file(path)==digest for path,digest in prior['inputs'].items())
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    eager_statistics_work(2049,128,512,2,True)
    rows=[]
    if args.pairs<1:raise ValueError('Positive pair count required')
    for pair in range(args.pairs):
        for reuse in ((False,True) if pair%2==0 else (True,False)):
            collections=[];began_gc={}
            def observe_gc(phase,info):
                generation=info['generation']
                if phase=='start':began_gc[generation]=(time.perf_counter(),time.thread_time())
                else:
                    began,used=began_gc.pop(generation)
                    collections.append(dict(generation=generation,collected=info['collected'],
                        wall_seconds=time.perf_counter()-began,cpu_seconds=time.thread_time()-used))
            gc.callbacks.append(observe_gc)
            cache=PlanningWorkCache();wall=time.perf_counter();cpu=time.thread_time()
            try:result=replay(prior,cache,reuse=reuse)
            finally:gc.callbacks.remove(observe_gc)
            result.update(pair=pair,reuse=reuse,wall_seconds=time.perf_counter()-wall,cpu_seconds=time.thread_time()-cpu)
            result['gc_collections']=collections
            cache.close();rows.append(result)
            print(json.dumps({key:result[key] for key in ('pair','reuse','wall_seconds','cpu_seconds','header_cache','stages')}),flush=True)
    for reuse in (False,True):
        profiler=cProfile.Profile();cache=PlanningWorkCache();profiler.enable()
        replay(prior,cache,reuse=reuse);profiler.disable();cache.close()
        label='reuse' if reuse else 'control'
        profiler.dump_stats(out/(label+'.prof'));stream=io.StringIO()
        pstats.Stats(profiler,stream=stream).sort_stats('cumulative').print_stats(55)
        (out/(label+'.txt')).write_text(stream.getvalue())
    assert source_identity()==source and all(sha256_file(path)==digest for path,digest in prior['inputs'].items())
    assert all(sha256_file(path)==digest for path,digest in helpers.items())
    summary={str(reuse):dict(median_wall_seconds=statistics.median(r['wall_seconds'] for r in rows if r['reuse']==reuse),
        median_cpu_seconds=statistics.median(r['cpu_seconds'] for r in rows if r['reuse']==reuse)) for reuse in (False,True)}
    report=dict(source_sha256=source,script_sha256=sha256_file(__file__),helper_sha256=helpers,
        prior_report_sha256=sha256_file(previous),records=rows,summary=summary,
        scope='Alternating calculator-only replay of one actual held productive source frontier. Control disables signature reuse and copies source templates before replacing encoded data; reuse retains at most 1024 structural signatures and copies only retained template fields. New header and planning cache every replay; all constructor and three-horizon costs included. Garbage collection remains enabled with unchanged thresholds; callbacks report pauses without subtracting them from reported time. Exact model reports/proposal/audit equality to the prior run. Synthetic prices, no runtime qualification or GWAS throughput claim. Instrumented profiles are attribution only.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
