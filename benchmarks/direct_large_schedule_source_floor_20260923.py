"""Read-only whole-PGEN source-floor construction cost on the local data copy.

Synthetic primitive rates exercise conservation arithmetic only. Report no
GWAS runtime prediction or candidate ranking from these rates.
"""
import hashlib
import json
from pathlib import Path
import resource
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_work_bounds import PgenHeaderWork,native_schedule_source_floor


SOURCE=Path('/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen')
TARGET=Path('results/large_schedule_source_floor_v2_20260923')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    TARGET.mkdir(parents=True,exist_ok=False)
    before=input_identity(SOURCE)
    if before['bytes']!=20_838_552_600:
        raise ValueError('Large-source byte identity differs from the control')
    began,cpu=time.perf_counter(),time.process_time()
    header=PgenHeaderWork(SOURCE)
    opened=dict(wall_seconds=time.perf_counter()-began,
                cpu_seconds=time.process_time()-cpu)
    markers=int(header._header.variant_ct)
    if markers!=8_086_101:
        raise ValueError('Large-source variant count differs from the control')
    rows=[]
    for chunk,expected in ((1024,(7897,1903,20_826_445_834)),
                           (4096,(1975,423,20_820_180_640))):
        began,cpu=time.perf_counter(),time.process_time()
        schedule=header.schedule_bounds(0,markers,chunk,
            max_records=10_000_000,max_chunks=10_000,max_signatures=65_536)
        schedule_time=dict(wall_seconds=time.perf_counter()-began,
                           cpu_seconds=time.process_time()-cpu)
        if (schedule['chunk_count'],schedule['ld_replay_count'],schedule['read_bytes'])!=expected:
            raise ValueError('Whole-source schedule differs from the prior exact metadata control')
        prices={name:1e-8 for name,pair in schedule['source_units'].items() if pair[1]}
        profile=dict(decode_units=prices,cpu_fraction=.5,depth=4,decode_workers=4,
            cpu_available_cores=4.,shared_dram_bytes_per_second=1e11,
            read_bytes_per_second=1e10,
            input_read_cpu_prices=dict(cpu_seconds_per_byte=1e-10,
                                       cpu_seconds_per_call=1e-5))
        began,cpu=time.perf_counter(),time.process_time()
        floor=native_schedule_source_floor(schedule,profile,
            dict(cpu=4.,dram=1e11,input=1e10))
        priced_time=dict(wall_seconds=time.perf_counter()-began,
                         cpu_seconds=time.process_time()-cpu)
        if (floor['resource_work']['input_bytes']!=expected[2]
                or floor['chunk_count']!=expected[0]
                or floor['source_stage_floor_seconds'][0]>floor['source_stage_floor_seconds'][1]):
            raise ValueError('Aggregate source-floor arithmetic differs')
        rows.append(dict(chunk_markers=chunk,chunks=expected[0],ld_replays=expected[1],
                         indexed_read_bytes=expected[2],schedule=schedule_time,
                         source_floor=priced_time))
    if input_identity(SOURCE)!=before:
        raise ValueError('PGEN source changed during the diagnostic')
    report=dict(source=before,header_open=opened,rows=rows,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        source_model_sha256=digest(Path(__file__).parents[1]/'src/torchgwas/pgen_work_bounds.py'),
        script_sha256=digest(__file__),
        scope='Read-only local-data PGEN metadata, no genotype payload, GPU or GWAS. '
              'Primitive rates and capacities are synthetic arithmetic fixtures; '
              'only construction CPU/wall cost and conserved indexed work are measured. '
              'Source-stage floors are deliberately omitted from this report and '
              'cannot qualify a pipeline prediction or chunk switch.')
    with (TARGET/'report.json').open('x') as out:
        json.dump(report,out,indent=2)
    print(json.dumps(dict(rows=rows,header_open=opened,
                          peak_rss_kib=report['peak_rss_kib'])),flush=True)


if __name__=='__main__':
    main()
