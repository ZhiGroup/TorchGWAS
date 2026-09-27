"""Read-only allocation evidence audit; never corrects or publishes prices."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics

from torchgwas.detailed_calibration import source_identity
from torchgwas.significant_host_work import host_significant_selection_work, host_selection_service


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    paths={backend:Path('results/selector_page_'+backend+'_prices_20260922/selection_prices.json')
        for backend in ['numpy','native']}
    paths['control']=Path('results/selector_first_touch_control_20260922/report.json')
    digests={key:hashlib.sha256(path.read_bytes()).hexdigest() for key,path in paths.items()}
    inputs={key:json.loads(path.read_text()) for key,path in paths.items()}
    control=inputs['control'];source=control['source_sha256']
    page=control['context']['page_bytes'];page_price=control['control_summary']['first_touch_cpu_seconds_per_page']
    assert control['observation_finished_at_utc'] and len(control['controls'])==7
    assert all(value['source_sha256']==source for value in inputs.values())
    rows=[];reference=[]
    for backend in ['numpy','native']:
        bank=inputs[backend]
        assert bank['context']['numpy_core']==control['context']['numpy_core']
        assert bank['allocation_observation_protocol']=='thread_faults_net_mallinfo2_separate_release_v1'
        for primitive,pattern in sorted({(r['primitive'],r['pattern']) for r in bank['observations']}):
            observations=[r for r in bank['observations'] if r['primitive']==primitive and r['pattern']==pattern]
            samples=[s for r in observations for s in r['samples']]
            assert all(s['major_faults']==0 for s in samples)
            reference.append(dict(backend=backend,primitive=primitive,pattern=pattern,
                median_cpu_seconds=statistics.median(r['cpu_seconds'] for r in observations),
                median_repeat_mean_minor_faults=statistics.median(statistics.mean(s['minor_faults'] for s in r['samples']) for r in observations),
                minor_faults=sum(s['minor_faults'] for s in samples),calls=len(samples),
                median_release_cpu_seconds=statistics.median(s['release_cpu_seconds'] for s in samples)))
        for shape in [(256,4096),(1024,8193)]:
            for density in ['sparse','dense']:
                selected=[r for r in control['selectors'] if r['backend']==backend and tuple(r['shape'])==shape and r['density']==density]
                assert len(selected)==7
                count=selected[0]['retained'];assert all(r['retained']==count for r in selected)
                os.environ['TORCHGWAS_HOST_PREDICATE']=backend
                work=host_significant_selection_work(*shape,count,return_beta=False)
                steps=host_selection_service(work,bank['prices'],cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)
                # select_host_pairs receives critical values, and its caller's
                # outer trait-coordinate rebasing is outside this control.
                predicted=sum(s['seconds']*s['resources']['cpu'] for s in steps[1:9])
                observed=statistics.median(r['observation']['cpu_seconds'] for r in selected)
                arrays=[r for r in work['allocations'] if r['lifetime']=='selected_result']
                assert sum(r['bytes'] for r in arrays)==selected[0]['selected_bytes']
                fresh_pages=sum(math.ceil(r['bytes']/page) for r in arrays)
                rows.append(dict(backend=backend,shape=shape,density=density,retained=count,
                    predicted_cpu_seconds=predicted,observed_cpu_seconds=observed,observed_over_predicted=observed/predicted,
                    median_system_cpu_seconds=statistics.median(r['observation']['system_cpu_seconds'] for r in selected),
                    median_minor_faults=statistics.median(r['observation']['minor_faults'] for r in selected),
                    median_release_cpu_seconds=statistics.median(r['release']['cpu_seconds'] for r in selected),
                    net_live_mapping_counts=[r['observation']['allocator_after']['hblks']-r['observation']['allocator_before']['hblks'] for r in selected],
                    net_live_mapping_bytes=[r['observation']['allocator_after']['hblkhd']-r['observation']['allocator_before']['hblkhd'] for r in selected],
                    selected_allocations=arrays,all_selected_pages_fresh_scenario=fresh_pages,
                    separate_all_fresh_cpu_scale_seconds=fresh_pages*page_price))
    assert digests=={key:hashlib.sha256(path.read_bytes()).hexdigest() for key,path in paths.items()}
    analysis_source=source_identity()
    report=dict(inputs={key:dict(path=str(path),sha256=digests[key],
        observation_started_at_utc=inputs[key]['observation_started_at_utc'],
        observation_finished_at_utc=inputs[key]['observation_finished_at_utc']) for key,path in paths.items()},
        observation_source_sha256=source,analysis_source_sha256=analysis_source,
        source_differences=[key for key in set(source)|set(analysis_source) if source.get(key)!=analysis_source.get(key)],
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        first_touch_control=control['control_summary'],reference=reference,cases=rows,
        artifacts_byte_unchanged=True,price_records_published=0,prediction_complete=False,selection_validated=False,
        scope='All-fresh page CPU scale is shown separately, never added to allocation-inclusive primitive prices. Net mallinfo2 deltas do not count transient allocations. Arena allocation does not guarantee resident pages. Source allocation extents do not certify future allocator routes, page state or release ownership.')
    (root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for row in rows:
        print(json.dumps({key:row[key] for key in ['backend','shape','density','predicted_cpu_seconds','observed_cpu_seconds','median_minor_faults','separate_all_fresh_cpu_scale_seconds']}),flush=True)


if __name__=='__main__':main()
