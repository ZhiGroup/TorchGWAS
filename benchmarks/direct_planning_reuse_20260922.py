"""Matched calculator reuse/cost gate audit; no production JIT claim."""
import argparse
from contextlib import nullcontext
import hashlib
import json
from pathlib import Path
import statistics
import time

from torchgwas.decoder_work import decoder_work
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.pgen_work_census import census
from torchgwas.planning_session import PlanningWorkCache,IncrementalPlanningBudget
from test_jagwas_actual_candidate import actual_candidate
from test_candidate_continuation import candidate_graph,future_choice


def digest(graph):
    return hashlib.sha256(json.dumps(graph.__dict__,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execution',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    prior_path=Path(args.execution)/'report.json';prior=json.loads(prior_path.read_text());identity=source_identity()
    for name,value in prior['inputs'].items():assert sha256_file(name)==value
    path=next(name for name in prior['inputs'] if name.endswith('.pgen'))
    source=census(path,128,include_chunks=True);choice=actual_candidate(path,512,2)
    for tile in choice['tiles']:
        tile['profile']['decode_units']={key:1e-9 for key in decoder_work(source,'torch_native_int8',restart_ld_bases=True)['source_units']}
    original=candidate_graph(choice);full=original.solve()
    at=min(value for name,value in full['end'].items() if name.endswith('jagwas:0:0:write:done'))
    checkpoint=original.checkpoint(at);remaining=original.resume(checkpoint)['remaining_seconds']
    cache=PlanningWorkCache();rows=[];expected={}
    # Alternate ordering so a single cold-first comparison cannot establish
    # the reported reuse effect. No whole-GWAS observations enter this audit.
    for repeat,order in enumerate(([128,256,512],[512,256,128],[256,128,512])):
        for size in order:
            for reuse in ([False,True] if repeat%2==0 else [True,False]):
                with cache.activate() if reuse else nullcontext():
                    began=time.perf_counter();cpu=time.process_time()
                    future,issued=future_choice(choice,source,checkpoint,size)
                    graph=candidate_graph(future)
                    build_cpu=time.process_time()-cpu;build_wall=time.perf_counter()-began
                    began=time.perf_counter();cpu=time.process_time()
                    solved=graph.resume(checkpoint)
                    resume_cpu=time.process_time()-cpu;resume_wall=time.perf_counter()-began
                checksum=digest(graph)
                if size not in expected:expected[size]=(checksum,solved['remaining_seconds'])
                assert (checksum,solved['remaining_seconds'])==expected[size]
                rows.append(dict(repeat=repeat,size=size,reused=reuse,issued_chunks=issued,
                    graph_sha256=checksum,remaining_seconds=solved['remaining_seconds'],
                    build_cpu_seconds=build_cpu,build_wall_seconds=build_wall,
                    resume_cpu_seconds=resume_cpu,resume_wall_seconds=resume_wall))
    costs={}
    for size in [128,256,512]:
        costs[str(size)]={}
        for reuse in [False,True]:
            selected=[row for row in rows if row['size']==size and row['reused']==reuse]
            costs[str(size)]['reused' if reuse else 'uncached']=dict(
                median_build_wall_seconds=statistics.median(row['build_wall_seconds'] for row in selected),
                median_total_wall_seconds=statistics.median(row['build_wall_seconds']+row['resume_wall_seconds'] for row in selected),
                median_cpu_seconds=statistics.median(row['build_cpu_seconds']+row['resume_cpu_seconds'] for row in selected))
    attempts=[]
    def forbidden():attempts.append(True);raise AssertionError('Cost gate invoked an unaffordable step')
    budget=IncrementalPlanningBudget(max_cpu_seconds=2.,max_window_seconds=10.)
    estimate=costs['128']['uncached']
    args=dict(remaining_seconds=remaining,expected_cpu_seconds=estimate['median_cpu_seconds'],
              expected_wall_seconds=estimate['median_total_wall_seconds'])
    early=budget.run_step(forbidden,**args)
    assert early['reason']=='no_useful_output_yet'
    budget.start_after_first_output()
    declined=budget.run_step(forbidden,**args)
    assert declined['reason']=='insufficient_remaining_horizon' and not attempts
    gate=budget.finish();reused=cache.snapshot();cache.close()
    assert source_identity()==identity
    write_record(out/'report.json',dict(source_sha256=identity,benchmark_sha256=sha256_file(__file__),
        helper_sha256={name:sha256_file(name) for name in ['tests/test_candidate_continuation.py',
            'tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py']},
        prior_execution_report_sha256=sha256_file(prior_path),inputs=prior['inputs'],dimensions=prior['dimensions'],
        baseline_remaining_seconds=remaining,rows=rows,summary=costs,cache=reused,
        startup_gate=early['reason'],short_horizon_gate=declined['reason'],budget=gate,
        scope='Alternating matched calculator evaluations on a retained real PGEN source with synthetic component prices. Graphs and continuation scores are identical with reuse. Cost-gate forecasts are audit controls; no live JIT execution, posterior calibration, end-to-end overhead or GWAS speedup is established.'))
    print(json.dumps(dict(summary=costs,cache=reused,short_horizon_gate=declined['reason'])),flush=True)


if __name__=='__main__':main()
