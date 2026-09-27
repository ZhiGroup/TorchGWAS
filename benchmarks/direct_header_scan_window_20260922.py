"""Bounded real-header/captured-geometry arithmetic audit, synthetic prices.

No GWAS is timed or fitted. Exact payload censuses are offline comparison
oracles only; the timed window/scenario path reads index metadata alone.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import time

from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.mechanistic_torch import torch_scan_work,torch_scan_header_work
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.pgen_work_census import census
from torchgwas.planning_session import PlanningWorkCache
from test_jagwas_actual_candidate import actual_candidate


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();path=Path(args.input);out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    sources=source_identity();digest=sha256_file(path)
    began=time.perf_counter();header=PgenHeaderWork(path);index_seconds=time.perf_counter()-began
    records=[];cache=PlanningWorkCache()
    with cache.activate():
        for b in (128,256,512):
            # Full census in this test helper is deliberately outside measured
            # window work. It is NOT an admissible production startup path.
            tile=actual_candidate(path,b,1)['tiles'][0];profile=tile['profile']
            for start,stop in [(0,2*b),(2048,2048+2*b),(4096,4097)]:
                began=time.perf_counter();window=header.window(start,stop,b)
                header_seconds=time.perf_counter()-began
                profile['decode_units']={name:1e-8 for row in window['chunks'] for name in row['source_units']}
                data=deepcopy(tile['data']);data.update(markers=stop-start,encoded=window)
                scenarios=[];timings=[]
                for endpoint in ('lower','upper'):
                    wall=time.perf_counter();cpu=time.thread_time()
                    work=torch_scan_header_work(data,profile,endpoint=endpoint)
                    timings.append(dict(endpoint=endpoint,cpu_seconds=time.thread_time()-cpu,wall_seconds=time.perf_counter()-wall))
                    scenarios.append(work)
                exact_data=deepcopy(data);exact_data['encoded']=census(path,b,variant_range=(start,stop),include_chunks=True)
                exact=torch_scan_work(exact_data,profile)
                for low,real,high in zip(scenarios[0]['blocks'],exact['blocks'],scenarios[1]['blocks']):
                    for resource in ('cpu','dram'):
                        amounts=[x['decode_seconds']*x['decode_resources'][resource] for x in (low,real,high)]
                        assert amounts[0]-1e-8<=amounts[1]<=amounts[2]+1e-8
                    for key in set(real)-{'decode_seconds','decode_resources'}:assert low[key]==real[key]==high[key],key
                records.append(dict(chunk_size=b,range=[start,stop],header_seconds=header_seconds,
                    pricing=timings,source_work=window['structural_work'],
                    scan_cpu_work_interval=[w['cpu_work_seconds'] for w in scenarios],
                    exact_cpu_seconds=exact['cpu_work_seconds'],matches_exact_nondecode_work=True,
                    d2h_bytes=sum(x['d2h_bytes'] for x in exact['blocks']),cache=cache.snapshot()))
                print(json.dumps(records[-1]),flush=True)
    assert source_identity()==sources and sha256_file(path)==digest
    report=dict(source_sha256=sources,script_sha256=sha256_file(__file__),
        helper_sha256={p:sha256_file(p) for p in ['tests/test_jagwas_actual_candidate.py','tests/test_jagwas_candidate.py',
            'tests/test_jagwas_scan_work.py','tests/fixtures/jagwas_chunk_geometry.json']},
        input_sha256=digest,input_identity=header.input_identity,input_mount=storage_identity(path),
        index_construction_seconds=index_seconds,rows=records,
        scope='Real PGEN headers and captured tensor geometry, with synthetic independent prices for arithmetic validation. No pipeline runtime, measured capacity, output-inclusive forecast, future-source representativeness, or autotuning benefit is established. Index construction and exact-census oracle work are outside window pricing times.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
