"""Counterfactual future chunks on a retained real source; synthetic services.

This measures calculator cost and checks model-state continuation. It does not
measure GWAS JIT overhead, infer hardware rates or execute an adaptive policy.
"""
import argparse
import copy
import json
from pathlib import Path
import time

from torchgwas.decoder_work import decoder_work
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.pgen_work_census import census
from test_jagwas_actual_candidate import actual_candidate
from test_candidate_continuation import candidate_graph,future_choice


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execution',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();previous=Path(args.execution)/'report.json';out=Path(args.out)
    out.mkdir(parents=True,exist_ok=False)
    executed=json.loads(previous.read_text());identity=source_identity()
    for name,digest in executed['inputs'].items():assert sha256_file(name)==digest
    path=next(name for name in executed['inputs'] if name.endswith('.pgen'))
    source=census(path,128,include_chunks=True);choice=actual_candidate(path,512,2)
    for tile in choice['tiles']:
        units=decoder_work(source,'torch_native_int8',restart_ld_bases=True)['source_units']
        tile['profile']['decode_units']={name:1e-9 for name in units}
    original=candidate_graph(choice);full=original.solve()
    at=min(value for name,value in full['end'].items() if name.endswith('jagwas:0:0:write:done'))
    state=original.checkpoint(at);saved=copy.deepcopy(state)
    baseline=original.resume(state);rows=[]
    for size in [128,256,512]:
        wall=time.perf_counter();cpu=time.process_time()
        future,issued=future_choice(choice,source,state,size)
        graph=candidate_graph(future)
        build_cpu=time.process_time()-cpu;build_wall=time.perf_counter()-wall
        wall=time.perf_counter();cpu=time.process_time()
        result=graph.resume(state)
        solve_cpu=time.process_time()-cpu;solve_wall=time.perf_counter()-wall
        assert all(result['start'][name]==value for name,value in state['start'].items())
        assert all(result['end'][name]==value for name,value in state['end'].items())
        if size==512:assert abs(result['seconds']-full['seconds'])<1e-10
        rows.append(dict(next_size=size,issued_chunks=issued,
            source_chunks=[len(tile['data']['encoded']['chunks']) for tile in future['tiles']],
            graph_nodes=len(graph.nodes),remaining_seconds=result['remaining_seconds'],
            forecast_gain_seconds=baseline['remaining_seconds']-result['remaining_seconds'],
            build_cpu_seconds=build_cpu,build_wall_seconds=build_wall,
            resume_cpu_seconds=solve_cpu,resume_wall_seconds=solve_wall))
    assert state==saved and source_identity()==identity
    write_record(out/'checkpoint.json',state)
    write_record(out/'report.json',dict(source_sha256=identity,benchmark_sha256=sha256_file(__file__),
        helper_sha256={name:sha256_file(name) for name in ['tests/test_candidate_continuation.py',
            'tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py']},
        prior_execution_report_sha256=sha256_file(previous),inputs=executed['inputs'],
        dimensions=executed['dimensions'],model_decision_seconds=at,
        completed_nodes=len(state['end']),running_nodes=len(state['remaining_service_seconds']),
        queued_consumers=sum(len(values) for values in state['fifo'].values()),
        baseline_remaining_seconds=baseline['remaining_seconds'],candidates=rows,
        prediction_validated=False,selection_validated=False,
        scope='Exact analytical scheduler continuation with fixed paid/issued prefix on the prior real native PGEN input. Synthetic component prices and modeled progress; not actual executor telemetry. Measured calculator CPU/wall costs are separate from forecast gains and are not end-to-end JIT overhead or speedup.'))
    print(json.dumps(dict(completed_nodes=len(state['end']),running_nodes=len(state['remaining_service_seconds']),
        candidates=rows)),flush=True)


if __name__=='__main__':main()
