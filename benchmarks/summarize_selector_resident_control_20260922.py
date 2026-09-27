"""Read-only resident-price transfer and forwarding-handler control audit."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
from torchgwas.detailed_calibration import source_identity
from torchgwas.significant_host_work import host_significant_selection_work,host_selection_service
from torchgwas.cpu_service_refresh import CpuServiceRefresh


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--input',default='results/selector_resident_allocator_control_20260922/report.json')
    parser.add_argument('--reference-mode',type=int,choices=[2,3],default=2);args=parser.parse_args()
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    path=Path(args.input)
    digest=hashlib.sha256(path.read_bytes()).hexdigest();data=json.loads(path.read_text())
    assert data['observation_finished_at_utc'] and data['source_sha256']==source_identity()
    medians={};reference=[];modes=sorted({r['mode'] for r in data['primitives']})
    assert args.reference_mode in modes
    stability=CpuServiceRefresh(root/'not_published','resident_audit',dependencies=dict(measurement_protocol=dict(audit=True)),
        work_units=1,max_age_seconds=1)
    for primitive,pattern in sorted({(r['primitive'],r['pattern']) for r in data['primitives']}):
        for mode in modes:
            rows=[r for r in data['primitives'] if r['primitive']==primitive and r['pattern']==pattern and r['mode']==mode]
            assert len(rows)==7 and {r['repeat'] for r in rows}==set(range(7))
            medians[primitive,pattern,mode]=statistics.median(r['residual_cpu_seconds'] for r in rows)
            samples=[s for row in rows for s in row['samples']]
            reference.append(dict(primitive=primitive,pattern=pattern,mode=mode,
                median_repeat_mean_cpu_seconds=statistics.median(r['cpu_seconds'] for r in rows),
                median_repeat_mean_residual_cpu_seconds=medians[primitive,pattern,mode],
                max_minor_faults=max(s['minor_faults'] for s in samples),
                total_minor_faults=sum(s['minor_faults'] for s in samples),
                max_major_faults=max(s['major_faults'] for s in samples),calls=len(samples),
                max_faults_outside_prefault=max(s['minor_faults']-s['allocator'].get('prefault_minor_faults',0) for s in samples),
                stability=stability._window_stability([dict(cpu_seconds=r['residual_cpu_seconds']) for r in sorted(rows,key=lambda r:r['repeat'])])))
    rates={};fixed={}
    for primitive in sorted({r['primitive'] for r in data['primitives']}):
        patterns={r['pattern']:r['units'] for r in data['primitives'] if r['primitive']==primitive}
        fixed[primitive]=medians[primitive,'fixed' if primitive=='mask_allocate' else 'empty',args.reference_mode]
        rates[primitive]={pattern:(medians[primitive,pattern,args.reference_mode]-fixed[primitive])/units
            for pattern,units in patterns.items() if units}
        assert all(rate>=0 for rate in rates[primitive].values()),(primitive,rates[primitive])
    cases=[]
    for shape in [(256,4096),(1024,8193)]:
        for density in ['empty','sparse','dense']:
            for backend in ['numpy','native']:
                rows=[r for r in data['selectors'] if tuple(r['shape'])==shape and r['density']==density and r['backend']==backend]
                count=rows[0]['retained'];assert all(r['retained']==count for r in rows)
                observed={}
                for mode in modes:
                    values=[r['observation'] for r in rows if r['mode']==mode];assert len(values)==7
                    observed[mode]=dict(cpu_seconds=statistics.median(r['cpu_seconds'] for r in values),
                        residual_cpu_seconds=statistics.median(r['residual_cpu_seconds'] for r in values),
                        minor_faults=statistics.median(r['minor_faults'] for r in values),
                        max_minor_faults=max(r['minor_faults'] for r in values),
                        release_cpu_seconds=statistics.median(r['release_cpu_seconds'] for r in values),
                        internal_allocator_cpu_seconds=statistics.median((r['allocator']['allocate_cpu_ns']+r['allocator']['free_cpu_ns'])*1e-9 for r in values),
                        prefault_cpu_seconds=statistics.median(r['allocator'].get('prefault_cpu_ns',0)*1e-9 for r in values),
                        prefault_minor_faults=statistics.median(r['allocator'].get('prefault_minor_faults',0) for r in values),
                        max_faults_outside_prefault=max(r['minor_faults']-r['allocator'].get('prefault_minor_faults',0) for r in values),
                        max_arena_used_bytes=max(r['allocator']['arena_used_bytes'] for r in values))
                os.environ['TORCHGWAS_HOST_PREDICATE']=backend
                work=host_significant_selection_work(*shape,count,return_beta=False)
                predictions={};terms={}
                for policy in ['maximum_declared_pattern','source_known_dense']:
                    prices={}
                    for primitive,patterns in rates.items():
                        rate=max(patterns.values()) if patterns else 0.
                        if policy=='source_known_dense' and primitive=='matrix_gather_flat' and count==shape[0]*shape[1]:rate=patterns['dense']
                        # Do not interpolate timings across held-out sizes or fit
                        # rates to selectors. Both policies use fixed controls.
                        name='predicate_block' if primitive=='predicate_numpy' else primitive
                        prices[name]=dict(call_cpu_seconds=fixed[primitive],unit_cpu_seconds=rate,dram_bytes_per_unit=0.)
                    for name in ['critical_lookup','index_cast','index_add']:
                        prices[name]=dict(call_cpu_seconds=0.,unit_cpu_seconds=0.,dram_bytes_per_unit=0.)
                    steps=host_selection_service(work,prices,cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)
                    selected=steps[1:9]
                    predictions[policy]=sum(s['seconds']*s['resources']['cpu'] for s in selected)
                    terms[policy]=[s['seconds']*s['resources']['cpu'] for s in selected]
                cases.append(dict(shape=shape,density=density,backend=backend,retained=count,
                    observation=observed,resident_predictions=predictions,resident_prediction_terms=terms,
                    resident_observed_over_source_prediction=observed[args.reference_mode]['residual_cpu_seconds']/predictions['source_known_dense'],
                    passthrough_over_default=observed[1]['cpu_seconds']/observed[0]['cpu_seconds'],
                    source_policy='Dense payload-gather price only when every pair is retained; maximum declared pattern otherwise. All other primitives retain maximum declared pattern.'))
    assert digest==hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(input_path=str(path),input_sha256=digest,source_sha256=source_identity(),
        input_observation_started_at_utc=data['observation_started_at_utc'],
        input_observation_finished_at_utc=data['observation_finished_at_utc'],
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        reference_mode=args.reference_mode,resident_fixed_cpu_seconds=fixed,resident_unit_cpu_seconds_by_pattern=rates,reference=reference,cases=cases,
        reference_qualification_warnings=[dict(primitive=r['primitive'],pattern=r['pattern'],
            stable=r['stability']['stable'],max_faults_outside_prefault=r['max_faults_outside_prefault'])
            for r in reference if r['mode']==args.reference_mode and
            (not r['stability']['stable'] or r['max_faults_outside_prefault']>4 or r['max_major_faults'])],
        artifacts_byte_unchanged=True,price_records_published=0,prediction_complete=False,selection_validated=False,
        scope='Resident arena controls remove demand paging and default allocation/free service, while retaining NumPy numerical kernels. Callback subtraction excludes measured allocator bodies but not all wrapper bookkeeping. Pool placement/cache state differs from the default allocator. Reference-price transfer is evaluated against resident held-out selectors; default timing differences are not fitted as coefficients.')
    (root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for row in cases:
        print(json.dumps(dict(shape=row['shape'],density=row['density'],backend=row['backend'],
            predicted=row['resident_predictions']['source_known_dense'],resident_observed=row['observation'][args.reference_mode]['residual_cpu_seconds'],
            ratio=row['resident_observed_over_source_prediction'],passthrough_ratio=row['passthrough_over_default'])),flush=True)


if __name__=='__main__':main()
