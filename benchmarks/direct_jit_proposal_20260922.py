"""Price one future-size alternative from a retained actual issue frontier.

Input and execution-prefix evidence are real. Component prices remain the
declared synthetic calculator controls. No production-benefit claim is made.
"""
import argparse
import json
from pathlib import Path
import time

from test_jagwas_actual_candidate import actual_candidate
from test_candidate_continuation import candidate_graph
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.jit_proposal import analytical_chunk_proposal
from torchgwas.pgen_work_census import census
from torchgwas.planning_session import PlanningWorkCache


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execution',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    evidence_path=Path(args.execution)/'report.json';evidence=json.loads(evidence_path.read_text())
    state=next(row for row in evidence['rows'] if row['label']=='scripted_future_switches')['first_snapshot']
    assert state['prefix_complete'] and state['finished'] is None and state['current_chunk_size']==128
    inputs=evidence['inputs']
    for name,digest in inputs.items():assert sha256_file(name)==digest
    path=next(Path(name) for name in inputs if name.endswith('/input.pgen'))
    candidate=actual_candidate(path,512,2);source=census(path,128,include_chunks=True)
    shifted_ld_controls={name:1e-9 for name in ['uleb1','uleb2','set_category',
        'difflist_group_absolute_id','difflist_category_extract','difflist_record_header']}
    for tile in candidate['tiles']:tile['profile']['decode_units'].update(shifted_ld_controls)
    sources=source_identity();cache=PlanningWorkCache();rows=[]
    for size in [512,256,128]:
        wall=time.perf_counter();cpu=time.thread_time()
        with cache.activate():
            proposal,audit=analytical_chunk_proposal(candidate,source_census=source,snapshot=state,
                chunk_sizes=[128,256,512],next_size=size,reduction='jagwas',graph_factory=candidate_graph)
        used=time.thread_time()-cpu;elapsed=time.perf_counter()-wall
        row=dict(proposal=proposal,audit=audit,planner_cpu_seconds=used,planner_wall_seconds=elapsed,
            improvement_exceeds_measured_planner_wall=audit['conditional_gain_floor_seconds']>elapsed)
        rows.append(row)
        print(json.dumps(dict(size=size,baseline_floor=proposal['baseline_seconds'],
            candidate_ceiling=proposal['candidate_seconds'],planning_wall=elapsed,
            conditional_gain_floor=audit['conditional_gain_floor_seconds'])),flush=True)
    for name,digest in inputs.items():assert sha256_file(name)==digest
    assert source_identity()==sources
    write_record(out/'report.json',dict(source_sha256=sources,benchmark_sha256=sha256_file(__file__),
        helper_sha256={name:sha256_file(name) for name in ['tests/test_candidate_continuation.py',
            'tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py']},
        execution_report_sha256=sha256_file(evidence_path),inputs=inputs,dimensions=evidence['dimensions'],
        source_mount=evidence['source_mount'],output_mount=evidence['output_mount'],
        actual_initial_snapshot=state,synthetic_shifted_ld_prices=shifted_ld_controls,rows=rows,cache=cache.snapshot(),
        scope='Offline source-calculator proposals from an actual held-prefix snapshot, using synthetic component prices. Exact reserved ranges and full phenotype JAGWAS are retained. Conditional model continuation bounds, not measured hardware intervals, Bayesian uncertainty, or live autotuning benefit.'))
    cache.close()


if __name__=='__main__':main()
