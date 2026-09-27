"""Read-only full-header PGEN schedule costs on the large local-data source."""
import hashlib
import json
from pathlib import Path
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork,paired_schedule_source_difference


SOURCE=Path('/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen')
TARGET=Path('results/large_schedule_bounds_v3_20260922')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8<<20),b''):h.update(block)
    return h.hexdigest()


def main():
    TARGET.mkdir(parents=True,exist_ok=False)
    source_before=input_identity(SOURCE)
    began=time.perf_counter();cpu=time.process_time()
    header=PgenHeaderWork(SOURCE)
    parsed=dict(wall_seconds=time.perf_counter()-began,cpu_seconds=time.process_time()-cpu)
    assert header.input_identity==source_before
    rows=[];bounds=[]
    for size in (128,1024,4096):
        began=time.perf_counter();cpu=time.process_time()
        bound=header.schedule_bounds(0,int(header._header.variant_ct),size,
            max_records=int(header._header.variant_ct),max_signatures=65536,max_chunks=100000)
        bounds.append(bound)
        rows.append(dict(chunk_markers=size,chunk_count=bound['chunk_count'],
            ld_replay_count=bound['ld_replay_count'],read_bytes=bound['read_bytes'],
            decode_input_bytes=bound['decode_input_bytes'],
            primary_payload_bytes=bound['record_payload_bytes'],
            source_unit_intervals=bound['source_units'],
            wall_seconds=time.perf_counter()-began,cpu_seconds=time.process_time()-cpu,
            bounds_cache=header.bounds_cache_info()))
        print(json.dumps({key:rows[-1][key] for key in ('chunk_markers','chunk_count',
            'ld_replay_count','wall_seconds','cpu_seconds','read_bytes')}),flush=True)
    assert input_identity(SOURCE)==source_before
    source_file=Path(__file__).parents[1]/'src/torchgwas/pgen_work_bounds.py'
    paired=[]
    for before,after in zip(bounds,bounds[1:]):
        began=time.perf_counter();cpu=time.process_time()
        difference=paired_schedule_source_difference(before,after,{})
        paired.append(dict(baseline_chunk=before['chunk_markers'],candidate_chunk=after['chunk_markers'],
            wall_seconds=time.perf_counter()-began,cpu_seconds=time.process_time()-cpu,
            read_bytes_delta=difference['read_bytes_delta'],
            decode_input_bytes_delta=difference['decode_input_bytes_delta'],
            chunk_count_delta=difference['chunk_count_delta'],
            ld_replay_count_delta=difference['ld_replay_count_delta'],
            source_unit_delta_intervals=difference['source_unit_delta_intervals'],
            unpriced_source_units=difference['unpriced_source_units'],
            decoder_cpu_delta_upper_bound=difference['decoder_cpu_delta_upper_bound']))
    report=dict(source=source_before,header_parse=parsed,rows=rows,paired_source_differences=paired,
        implementation_sha256=sha(source_file),script_sha256=sha(__file__),
        scope=__doc__+' Conditional source work only; no payload read, component pricing, GPU/output model, or performance switch.')
    with (TARGET/'report.json').open('x') as out:json.dump(report,out,indent=2)


if __name__=='__main__':main()
