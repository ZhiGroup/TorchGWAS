"""Preserve phase CPU/wall observations and paired instrumentation controls."""
import argparse
import json
from pathlib import Path
import statistics
from torchgwas.detailed_calibration import sha256_file
from torchgwas.geometry_collection import write_record

parser=argparse.ArgumentParser()
parser.add_argument('--root',required=True)
args=parser.parse_args()
root=Path(args.root)
report=json.loads((root/'report.json').read_text())
assert report['protocol_sha256']==sha256_file(root/'protocol.json')
phases=[];invocations=[]
for row in report['observations']:
    for trace in row['traces']:
        base={k:row[k] for k in ['repeat','scale','primitive','mode']}
        base['sample']=trace['sample']
        total={clock:(trace['end_'+clock]-trace['begin_'+clock])*1e-9 for clock in ['cpu','wall']}
        assert all(value>=0 for value in total.values())
        separated='return_cpu' in trace and 'return_wall' in trace
        api={clock:(trace.get('return_'+clock,trace['end_'+clock])-trace['begin_'+clock])*1e-9 for clock in ['cpu','wall']}
        release={clock:total[clock]-api[clock] for clock in ['cpu','wall']}
        assert all(api[clock]>=0 and release[clock]>=-1e-12 for clock in ['cpu','wall'])
        invocations.append(dict(base,whole=total,api=api,release=release,ownership_checkpoint=separated))
        if row['mode']!='record' or not row['primitive'].startswith(('copy_','nonzero_')):continue
        memcpy,wait=trace['events']
        assert [memcpy['kind'],wait['kind']]==[1,2]
        assert trace['count']==2 and trace['overflow']==0
        pieces={}
        for clock in ['cpu','wall']:
            stamps=[trace['begin_'+clock],memcpy['begin_'+clock],memcpy['end_'+clock],
                wait['begin_'+clock],wait['end_'+clock],trace.get('return_'+clock,trace['end_'+clock]),trace['end_'+clock]]
            assert all(a<=b for a,b in zip(stamps,stamps[1:]))
            pieces[clock]={name:(b-a)*1e-9 for name,a,b in zip(
                ['before_copy','copy_runtime','between_copy_and_wait','wait_runtime','after_wait','ownership_release'],stamps,stamps[1:])}
            assert abs(sum(pieces[clock].values())-total[clock])<1e-12
            pieces[clock]['nonzero_before_wait']=sum(pieces[clock][name] for name in ['before_copy','copy_runtime','between_copy_and_wait'])
            pieces[clock]['copy_and_wait']=sum(pieces[clock][name] for name in ['copy_runtime','between_copy_and_wait','wait_runtime'])
        phases.append(dict(base,whole=total,phase=pieces,bytes=memcpy['bytes'],
            copy_and_wait_cpu_fraction=pieces['cpu']['copy_and_wait']/pieces['wall']['copy_and_wait']))


def distribution(values):
    return dict(values=values,median=statistics.median(values),minimum=min(values),maximum=max(values))

summaries=[]
for scale,primitive in sorted({(r['scale'],r['primitive']) for r in invocations}):
    values=[r for r in invocations if (r['scale'],r['primitive'])==(scale,primitive)]
    repeats=sorted({r['repeat'] for r in values})
    medians={mode:{clock:[statistics.median(r['whole'][clock] for r in values if r['repeat']==repeat and r['mode']==mode)
        for repeat in repeats] for clock in ['cpu','wall']} for mode in ['plain','record']}
    summary=dict(scale=scale,primitive=primitive,repetitions=len(repeats),
        ownership_checkpoint=all(r['ownership_checkpoint'] for r in values),
        api_return_intervals={mode:{clock:distribution([statistics.median(r['api'][clock] for r in values if r['repeat']==repeat and r['mode']==mode)
            for repeat in repeats]) for clock in ['cpu','wall']} for mode in ['plain','record']},
        ownership_release_intervals={mode:{clock:distribution([statistics.median(r['release'][clock] for r in values if r['repeat']==repeat and r['mode']==mode)
            for repeat in repeats]) for clock in ['cpu','wall']} for mode in ['plain','record']},
        marker_intervals={mode:{clock:distribution(v) for clock,v in by_clock.items()} for mode,by_clock in medians.items()},
        paired_record_minus_plain={clock:distribution([a-b for a,b in zip(medians['record'][clock],medians['plain'][clock])]) for clock in ['cpu','wall']},
        paired_record_divided_by_plain={clock:distribution([a/b for a,b in zip(medians['record'][clock],medians['plain'][clock])]) for clock in ['cpu','wall']})
    records=[r for r in phases if (r['scale'],r['primitive'])==(scale,primitive)]
    if records:
        summary['phases']={clock:{name:distribution([statistics.median(r['phase'][clock][name] for r in records if r['repeat']==repeat)
            for repeat in repeats]) for name in records[0]['phase'][clock]} for clock in ['cpu','wall']}
        summary['copy_and_wait_cpu_fraction']=distribution([statistics.median(r['copy_and_wait_cpu_fraction'] for r in records if r['repeat']==repeat) for repeat in repeats])
    summaries.append(summary)
write_record(root/'phase_summary.json',dict(report_sha256=sha256_file(root/'report.json'),
    observations=phases,summaries=summaries,instrumentation_qualified=False,
    scope='Exact conservation of recorded marker and runtime intervals; medians within each eight-call repetition and across all nine repetitions. Negative paired differences are retained. No minimum clipping, regression, price publication, or GIL claim.',
    model_mapping={'nonzero':'before_wait CPU and after_wait CPU bracket the count barrier; CPU consumed inside cudaStreamSynchronize is waiting, not independent dispatch.',
        'copy':'before_copy and after_wait CPU bracket transport up to the API return checkpoint. cudaMemcpyAsync plus the gap and stream synchronization form the copy-and-wait interval; pageable runtime work may block.',
        'ownership':'When a return checkpoint exists, output destruction and return-marker exit/final-marker entry overhead remain in ownership_release. Legacy records do not separate these lifetimes.',
        'controls':'Both paths use identical boundary and return markers. Plain disables internal runtime clocks. Loop-control and ready-stream controls expose marker/driver overhead without subtracting it.'}))
for row in summaries:
    out={key:row[key] for key in ['scale','primitive']}
    out.update(plain_cpu_us=1e6*row['marker_intervals']['plain']['cpu']['median'],
        record_cpu_us=1e6*row['marker_intervals']['record']['cpu']['median'],
        paired_wall_ratio=row['paired_record_divided_by_plain']['wall']['median'])
    if 'phases' in row:
        out['phase_cpu_us']={key:1e6*value['median'] for key,value in row['phases']['cpu'].items()}
        out['copy_wait_cpu_fraction']=row['copy_and_wait_cpu_fraction']['median']
    print(json.dumps(out),flush=True)
