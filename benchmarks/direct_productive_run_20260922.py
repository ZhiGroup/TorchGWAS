"""Actual two-GPU productive lifecycle, correctness and matched overhead audit.

The proposal values below are deliberately scripted controls, not calibrated
predictions or an autotuning policy. Fixed-size enabled/disabled pairs isolate
lifecycle overhead on a warm, small fixture; all early associations are kept.
"""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import time

import numpy as np
import torch

from test_jagwas_actual_candidate import actual_candidate
from direct_jagwas_bounded_execution_20260922 import read_result
from torchgwas.adaptive_candidate import adaptive_candidate_memory
from torchgwas.adaptive_chunks import AlignedChunkSizeControl
from torchgwas.api import load_genotype
from torchgwas.detailed_calibration import source_identity, sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.linear import linear_scan_multigpu
from torchgwas.pgen_work_census import census
from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.reduce import JagwasReduction
from torchgwas.sumstats_indexed import write_indexed_sumstats


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--execution',required=True)
    parser.add_argument('--out',required=True)
    parser.add_argument('--output-data',help='Separate remote local-storage directory for association files')
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--admit-start',action='store_true',help='Use the single-layout structural startup path')
    args=parser.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    if args.repeats<1:raise ValueError('At least one matched pair required')
    output_data=Path(args.output_data) if args.output_data else out
    if args.output_data:output_data.mkdir(parents=True,exist_ok=False)
    prior_path=Path(args.execution)/'report.json'
    prior=json.loads(prior_path.read_text())
    inputs=prior['inputs'];path=next(Path(p) for p in inputs if p.endswith('/input.pgen'))
    for name,digest in inputs.items():assert sha256_file(name)==digest
    fixture=path.parent
    y=np.load(fixture/'phenotype.npy',mmap_mode='r');cov=np.load(fixture/'covariates.npy')
    n,m,k,c=[prior['dimensions'][name] for name in ['N','M','K','C']]
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    source_hashes=source_identity()
    candidate=actual_candidate(path,512,2)
    profiles={}
    import os
    for device in candidate['devices']:
        p=torch.cuda.get_device_properties(device)
        profiles[device]=dict(torch_version=torch.__version__,sm_count=p.multi_processor_count,
            max_threads_per_sm=p.max_threads_per_multi_processor,compute_capability=[p.major,p.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
    startup=None
    from torchgwas.api import _available_host_bytes
    if args.admit_start:
        from torchgwas.adaptive_start import prepare_adaptive_start
        context=dict(name='two_gpu_fixture',devices=candidate['devices'],
            profiles={tile['device']:tile['profile'] for tile in candidate['tiles']},
            shared_capacities=candidate['shared_capacities'])
        workload=dict(genotype=str(path),samples=n,markers=m,traits=k,covariates=c,
            matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
        startup=prepare_adaptive_start(workload,context,chunk_sizes=[128,256,512],initial_size=128,
            partition_axis='variant',trait_block=None,reduction='jagwas',output=candidate['output'],
            cpu_workers=3,host_memory_bytes=_available_host_bytes(),
            device_memory_bytes={d:torch.cuda.mem_get_info(d)[0] for d in candidate['devices']},
            host_reserve_bytes=512<<20,device_reserve_bytes=256<<20,device_memory_profiles=profiles)
        candidate=startup['candidate'];fine=startup['source_layout'];admission=startup['memory']
        write_record(out/'startup.json',startup)
        assert startup['structural_work']['runtime_candidates_evaluated']==0
        assert startup['structural_work']['census_passes']==0
        assert startup['structural_work']['header_passes']==1
    else:
        fine=census(path,128,include_chunks=True)
        admission=adaptive_candidate_memory(candidate,chunk_sizes=[128,256,512],source_census=fine,
            reduction='jagwas',device_memory_profiles=profiles,host_reserve_bytes=512<<20,device_reserve_bytes=256<<20)
    assert not admission['missing_geometry']
    from torchgwas.api import _available_host_bytes
    assert admission['host_bytes']<_available_host_bytes()
    for device,value in admission['device_bytes'].items():assert value<torch.cuda.mem_get_info(device)[0]
    write_record(out/'admission.json',admission)
    partitions=[dict(id=str(i),device=tile['device'],variant_range=tile['variant_range'],trait_range=tile['trait_range'])
                for i,tile in enumerate(candidate['tiles'])]
    rows=[];reference=None

    def execute(label,enabled,*,switch=False):
        nonlocal reference
        directory=output_data/label
        for device in candidate['devices']:torch.cuda.synchronize(device)
        started=time.perf_counter()
        control=ProductiveTuningRun(partitions,chunk_sizes=[128,256,512],initial=128,
            max_issued_chunks=128,budget=IncrementalPlanningBudget(max_steps=2,
                max_cpu_seconds=.1,max_window_seconds=10.)) if enabled else AlignedChunkSizeControl([128,256,512],initial=128)
        if enabled:
            refused=control.planning_step(lambda snapshot: (_ for _ in ()).throw(AssertionError('startup planning')),
                remaining_seconds=10.,expected_cpu_seconds=.001,expected_wall_seconds=.002)
            assert not refused['evaluated'] and refused['reason']=='no_useful_output_yet'
        source=load_genotype(path,genotype_format='pgen',pgen_mode='hardcall',reader_workers=3)[0]
        chunks,basis=linear_scan_multigpu(source,y,cov,devices=candidate['devices'],chunk_size=512,
            reader_workers=3,prefetch_chunks=2,compute_dtype='float32',compute_p_values=False,
            ordered=False,shared_queue_depth=2,reduction_factory=JagwasReduction,
            _chunk_size_selector=control)
        delivered=[];events=[];steps=[];first_snapshot=None
        def checked():
            try:
                for item in chunks:
                    delivered.append([int(item[0]),int(item[1])]);yield item
            finally:chunks.close()
        def written(event):
            nonlocal first_snapshot
            events.append(event)
            if enabled:
                control.output_written(event)
                if first_snapshot is None:first_snapshot=control.snapshot()
                if switch and len(events) in (1,3):
                    size=256 if len(events)==1 else 512
                    # Lifecycle test only: estimates do not claim a speedup.
                    steps.append(control.planning_step(lambda snapshot:dict(chunk_size=size,
                        baseline_seconds=10.,candidate_seconds=5.),remaining_seconds=10.,
                        expected_cpu_seconds=.001,expected_wall_seconds=.002,expected_gain_seconds=5.))
        try:
            total,summary=write_indexed_sumstats(directory/'sumstats',[f'v{i}' for i in range(m)],
                [f't{i}' for i in range(k)],n,checked(),kind='jagwas',df=n-basis.shape[1]-2,
                chi2_df=k,store_beta=False,fsync=True,on_chunk_written=written)
        except BaseException:
            if enabled:control.finish(successful=False)
            raise
        snapshot=control.finish(successful=True) if enabled else None
        completed=time.perf_counter()
        _,values=read_result(directory,m)
        if reference is None:reference=values.copy()
        np.testing.assert_allclose(values,reference,rtol=6e-5,atol=3e-4)
        cursor=0
        for lo,hi in sorted(delivered):assert lo==cursor;cursor=hi
        assert cursor==m and total==m and len(set(map(tuple,delivered)))==len(delivered)
        assert events and all(event.part_file_fsynced for event in events)
        if enabled:
            ranges=sorted(span for partition in snapshot['partitions'] for span in partition['ranges'])
            assert ranges==sorted(delivered)
            assert sum(p['issued_chunks'] for p in snapshot['partitions'])==len(delivered)
            assert sum(p['issued_chunks'] for p in first_snapshot['partitions'])>1
            assert snapshot['part_bytes']==sum(event.part_bytes for event in events)
            assert snapshot['first_written']==events[0].completed
            assert snapshot['first_fsynced_part']==events[0].completed
            assert snapshot['finished']['successful'] and snapshot['cache']['closed']
        if switch:
            assert len(steps)==2 and all(step['applied'] for step in steps)
            assert {128,256,512}<=set(hi-lo for lo,hi in delivered)
        row=dict(label=label,enabled=enabled,scripted_switch=switch,
            time_to_first_written_seconds=events[0].completed-started,
            output_inclusive_seconds=completed-started,
            first_snapshot=first_snapshot,final_snapshot=snapshot,source_ranges=sorted(delivered),
            retained_variants=total,parts=len(events),part_bytes=sum(event.part_bytes for event in events),
            max_absolute_reference_difference=float(np.max(np.abs(values-reference))),steps=steps)
        rows.append(row)
        print(json.dumps({key:row[key] for key in ['label','time_to_first_written_seconds','output_inclusive_seconds','parts']}),flush=True)

    execute('warmup',False)
    for repeat in range(args.repeats):
        for enabled in ([False,True] if repeat%2==0 else [True,False]):
            execute(f'pair{repeat}_'+('tracked' if enabled else 'fixed'),enabled)
    execute('scripted_future_switches',True,switch=True)
    measured=[row for row in rows if row['label'].startswith('pair')]
    summary={label:{metric:statistics.median(row[metric] for row in measured if row['enabled']==enabled)
        for metric in ['time_to_first_written_seconds','output_inclusive_seconds']}
        for label,enabled in [('fixed',False),('tracked',True)]}
    paired=[]
    for repeat in range(args.repeats):
        pair={row['enabled']:row for row in measured if row['label'].startswith(f'pair{repeat}_')}
        paired.append(dict(repeat=repeat,**{metric:pair[True][metric]-pair[False][metric]
            for metric in ['time_to_first_written_seconds','output_inclusive_seconds']}))
    mounts={key:json.loads(subprocess.check_output(['findmnt','-J','-T',str(value),
        '-o','SOURCE,FSTYPE,TARGET'],text=True)) for key,value in [('input',path),('output',output_data)]}
    for name,digest in inputs.items():assert sha256_file(name)==digest
    assert source_identity()==source_hashes
    write_record(out/'report.json',dict(source_sha256=source_hashes,benchmark_sha256=sha256_file(__file__),
        helper_sha256={name:sha256_file(name) for name in ['tests/test_jagwas_actual_candidate.py',
            'tests/test_jagwas_candidate.py','benchmarks/direct_jagwas_bounded_execution_20260922.py']},
        prior_execution_report_sha256=sha256_file(prior_path),inputs=inputs,dimensions=prior['dimensions'],
        source_mount=mounts['input'],output_mount=mounts['output'],output_data_directory=str(output_data),
        rows=rows,matched_medians=summary,paired_added_seconds=paired,
        adaptive_start=None if startup is None else {key:startup[key] for key in
            ['admission_seconds','structural_work','api_kwargs','partitions','input_file_identity','capacity','initial_size']},
        scope='Warm small-fixture two-GPU JAGWAS lifecycle and exact coverage audit. First-output timestamps follow part-file fsync; final output includes manifest writing. Both matched arms retain the same minimal writer callback. Scripted proposal forecasts test plumbing only, not analytical prediction, Bayesian policy, or production autotuning benefit. Admission is outside matched scan timing; total public startup remains unqualified.'))
    print(json.dumps(dict(matched_medians=summary)),flush=True)


if __name__=='__main__':main()
