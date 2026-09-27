"""Real held-frontier JAGWAS model callback; synthetic prices, no speedup claim."""
import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import threading
import time

import numpy as np
import torch
from test_jagwas_scan_work import joint_fixture
from test_jagwas_actual_candidate import CAPTURE,writer_prices
from direct_jagwas_bounded_execution_20260922 import read_result
from torchgwas.adaptive_chunks import AlignedChunkSizeControl
from torchgwas.adaptive_start import prepare_adaptive_start
from torchgwas.api import load_genotype,_available_host_bytes
from torchgwas.detailed_calibration import source_identity,sha256_file,storage_identity
from torchgwas.linear import linear_scan_multigpu
from torchgwas.pgen_reader import read_header
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.planning_session import IncrementalPlanningBudget
from torchgwas.productive_run import ProductiveTuningRun
from torchgwas.reduce import JagwasReduction
from torchgwas.sumstats_indexed import write_indexed_sumstats,IndexedOutputPartition
from torchgwas.tensor_service import host_primitive_name
from torchgwas.tensor_work import eager_statistics_work
from torchgwas.window_model import compare_prepared_windows,prepared_source_window


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fixture',required=True)
    parser.add_argument('--out',required=True);parser.add_argument('--output-data',required=True)
    parser.add_argument('--copy-price-cache')
    parser.add_argument('--copy-price-max-age',type=float,default=300.)
    args=parser.parse_args();fixture=Path(args.fixture);out=Path(args.out);destination=Path(args.output_data)
    out.mkdir(parents=True,exist_ok=False);destination.mkdir(parents=True,exist_ok=False)
    path=fixture/'input.pgen';inputs={str(fixture/name):sha256_file(fixture/name) for name in
        ('input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy')}
    source=source_identity();torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    y=np.load(fixture/'phenotype.npy');cov=np.load(fixture/'covariates.npy');header=read_header(path)
    n,m,k=header.sample_ct,header.variant_ct,y.shape[1]
    assert (n,k,cov.shape[1])==(2049,512,2)
    devices=['cuda:0','cuda:1'];_,template=joint_fixture();caps=dict(cpu=3.,dram=1e9,input=1e8,output=1e8)
    profiles={};memory_profiles={};began=time.perf_counter()
    current_context=None
    if args.copy_price_cache:
        from torchgwas.detailed_calibration import execution_context
        current_context=execution_context(devices,input_path=path,output_path=destination)
    for index,device in enumerate(devices):
        profile=deepcopy(template)
        profile.update(result_ownership='owned',validate_range=True,event_wait_cpu_fraction=1.,
            decode_workers=2-index,cpu_available_cores=3.,chunk_markers=256,
            write_bytes_per_second=1e8,fsync_seconds=1e-6,
            writeback_service=dict(pagecache_seconds_per_byte=1e-9,storage_seconds_per_byte=1e-8),
            kernel_geometry=deepcopy(CAPTURE['statistics']),joint_kernel_geometry=deepcopy(CAPTURE['projection']),
            host_primitives={host_primitive_name(call):1e-6 for call in eager_statistics_work(n,128,k,2,True)['host_calls']})
        # Explicit synthetic decoder controls, not rates inferred from the job.
        profile['decode_units'].update({name:1e-8 for name in ['uleb1','uleb2','uleb3','uleb4','uleb5',
            'set_category','difflist_group_absolute_id','difflist_category_extract','difflist_record_header']})
        profiles[device]=profile;p=torch.cuda.get_device_properties(device)
        memory_profiles[device]=dict(torch_version=torch.__version__,sm_count=p.multi_processor_count,
            max_threads_per_sm=p.max_threads_per_multi_processor,compute_capability=[p.major,p.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
    start=prepare_adaptive_start(dict(genotype=str(path),samples=n,markers=m,traits=k,covariates=2,
        matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True),
        dict(name='synthetic-price-live-bridge',devices=devices,profiles=profiles,shared_capacities=caps),
        chunk_sizes=[128,256],initial_size=128,partition_axis='variant',trait_block=None,reduction='jagwas',
        output=dict(block_bytes=None,queue_depth=2,store_beta=False,fsync=True),cpu_workers=3,
        host_memory_bytes=_available_host_bytes(),device_memory_bytes={d:torch.cuda.mem_get_info(d)[0] for d in devices},
        host_reserve_bytes=512<<20,device_reserve_bytes=256<<20,device_memory_profiles=memory_profiles)
    preparation=time.perf_counter()-began
    assert start['structural_work']['census_passes']==0 and start['structural_work']['runtime_candidates_evaluated']==0
    assert not start['memory']['missing_geometry']
    options=dict(total_traits=k,reduction='jagwas',output=start['candidate']['output'],shared_capacities=caps,
        endpoint='upper',host_serial_fraction=.5,prices=writer_prices(),max_source_chunks=64,
        model_identity=dict(source_sha256=source,scope='Synthetic independent prices with captured untimed geometry.'))
    results=[];reference=None;component_evidence=None;bound_prices=None
    for label in ('control','forecast'):
        enabled=label=='forecast';began=time.perf_counter();events=[];delivered=[];steps=[];comparisons=[];snapshot_at_model=None
        probe=None;measurement_steps=[];model_attempted=False;measurement_finished=False;tuning_closed=False
        controller=ProductiveTuningRun(start['partitions'],chunk_sizes=[128,256],initial=128,
            budget=IncrementalPlanningBudget(max_steps=5 if args.copy_price_cache else 1,max_cpu_seconds=10.,max_window_seconds=30.)) if enabled else AlignedChunkSizeControl([128,256],initial=128)
        def measure(snapshot):
            nonlocal component_evidence,bound_prices,probe,measurement_finished
            from direct_initial_component_price_20260922 import CopyPriceProbe
            from torchgwas.calibration_cache import _digest
            contexts=[dict(name='productive-copy-price-audit',devices=devices,
                profiles={tile['device']:tile['profile'] for tile in start['candidate']['tiles']})]
            if probe is None:
                probe=CopyPriceProbe(args.copy_price_cache,contexts,source,current_context,
                    max_age_seconds=args.copy_price_max_age)
            bound_prices,component_evidence=probe.advance()
            measurement_steps.append(dict(written_parts=len(events),state=deepcopy(component_evidence)))
            measurement_finished=component_evidence['state'] not in ('checking','measuring')
            if bound_prices is not None:
                for tile in start['candidate']['tiles']:
                    tile['profile']=bound_prices['contexts'][0]['profiles'][tile['device']]
                options['model_identity']=dict(source_sha256=source,price_profile_sha256=_digest(bound_prices),
                    scope='One measured resident copy CPU coefficient; all other prices are synthetic controls.')
            return dict(chunk_size=snapshot['current_chunk_size'],baseline_seconds=0.,candidate_seconds=0.)
        def build(snapshot):
            nonlocal snapshot_at_model
            snapshot_at_model=snapshot;header_work=PgenHeaderWork(path)
            def priced_work(data,profile,**kwargs):
                from torchgwas.mechanistic_torch import torch_scan_header_work
                work=torch_scan_header_work(data,profile,**kwargs)
                coefficient=component_evidence['cpu_seconds_per_byte']
                assert profile['owned_result_copy_scenario']['fresh_fraction']==0.
                for shape,row in work['owned_result_work'].items():
                    service=row['service'];expected=service['additional_copy_bytes']*coefficient
                    assert expected>0 and service['additional_copy_cpu_seconds']==expected
                    component_evidence.setdefault('consumed_copy_prices',[]).append(dict(chunk_markers=shape,
                        additional_copy_bytes=service['additional_copy_bytes'],copy_cpu_seconds=expected,
                        measured_cpu_seconds_per_byte=coefficient))
                return work
            remaining=[row['variant_range'][1]-row['cursor'] for row in snapshot['partitions'] if row['cursor']<row['variant_range'][1]]
            if not remaining or min(remaining)<384:
                raise ValueError('Insufficient unissued source for three admitted audit horizons')
            unit=256 if min(remaining)>=768 else 128
            for count in (unit,2*unit,3*unit):
                layouts=[];bins=[]
                for size in (snapshot['current_chunk_size'],256):
                    windows=[]
                    for row,tile in zip(snapshot['partitions'],start['candidate']['tiles']):
                        if row['cursor']==row['variant_range'][1]:continue
                        lo=row['cursor'];assert lo+count<=row['variant_range'][1]
                        w=prepared_source_window(tile,header_work,start=lo,stop=lo+count,chunk_markers=size,
                            issued_chunks=row['issued_chunks'],expected_input_identity=start['input_file_identity'],max_chunks=16)
                        windows.append(w)
                        if len(layouts)==0:bins.append(dict(variant_range=[lo,lo+count],trait_range=[0,k],retained=count))
                    layouts.append(dict(windows=windows,partition_axis='variant'))
                evidence=dict(input_identity=header_work.input_identity,reduction='jagwas',total_traits=k,
                    significance_threshold=None,bins=bins)
                if bound_prices is None:
                    comparisons.append(compare_prepared_windows(*layouts,survivor_evidence=evidence,**options))
                else:
                    from unittest.mock import patch
                    with patch('torchgwas.window_model.torch_scan_header_work',side_effect=priced_work):
                        comparisons.append(compare_prepared_windows(*layouts,survivor_evidence=evidence,**options))
            return comparisons
        def written(event):
            nonlocal model_attempted,tuning_closed
            events.append(event)
            if enabled:
                controller.output_written(event)
                if tuning_closed:return
                if args.copy_price_cache and not measurement_finished:
                    result=controller.planning_step(measure,remaining_seconds=10.,expected_cpu_seconds=.01,
                        expected_wall_seconds=.02,switching_seconds=0.,publication_seconds=0.,reserve_seconds=.01)
                    steps.append(result)
                    tuning_closed=not result.get('usable_for_decision',False)
                    return
                if not model_attempted and (not args.copy_price_cache or bound_prices is not None):
                    model_attempted=True
                    steps.append(controller.forecast_step(build,price_profile=bound_prices,forecast_options=dict(
                        boundary_adjustments=dict(baseline=[0.,1.],candidate=[0.,1.]),
                        relative_model_error=.05,max_slope_change=.1,max_extrapolation=8.,
                        assumptions=dict(source_work='The bounded header-work scenario represents unissued source work.',
                            output_occupancy='Declared all-valid full-panel JAGWAS scenario.',
                            resource_capacity='Synthetic independent-price controls, not qualified runtime capacities.',
                            partition_balance='Exact remaining extents use the balanced core/envelope scenarios.')),
                        remaining_seconds=10.,expected_cpu_seconds=.01,expected_wall_seconds=.02,
                        switching_seconds=0.,publication_seconds=0.,reserve_seconds=.01))
        def partition_for_range(lo,hi):
            rows=[p for p in start['partitions'] if p['variant_range'][0]<=lo<hi<=p['variant_range'][1]]
            assert len(rows)==1;p=rows[0]
            return IndexedOutputPartition(p['device'],tuple(p['variant_range']),tuple(p['trait_range']))
        genotype=load_genotype(path,genotype_format='pgen',pgen_mode='hardcall',reader_workers=3)[0]
        chunks,basis=linear_scan_multigpu(genotype,y,cov,devices=devices,chunk_size=256,reader_workers=3,
            prefetch_chunks=2,compute_dtype='float32',compute_p_values=False,ordered=False,shared_queue_depth=2,
            reduction_factory=JagwasReduction,_chunk_size_selector=controller)
        def checked():
            try:
                for item in chunks:delivered.append([int(item[0]),int(item[1])]);yield item
            finally:chunks.close()
        try:
            total,_=write_indexed_sumstats(destination/label/'sumstats',[f'v{i}' for i in range(m)],
                [f't{i}' for i in range(k)],n,checked(),kind='jagwas',df=n-basis.shape[1]-2,chi2_df=k,
                store_beta=False,fsync=True,on_chunk_written=written,partition_for_range=partition_for_range)
        except BaseException:
            if enabled:controller.finish(successful=False)
            raise
        finally:
            if probe is not None:probe.close()
        final=controller.finish(successful=True) if enabled else None;elapsed=time.perf_counter()-began
        assert not [t.name for t in threading.enumerate() if t.name.startswith('torchgwas-')]
        _,values=read_result(destination/label,m)
        if reference is None:reference=values.copy()
        np.testing.assert_allclose(values,reference,rtol=6e-5,atol=3e-4)
        cursor=0
        for lo,hi in sorted(delivered):assert lo==cursor;cursor=hi
        assert cursor==m and total==m
        if enabled:
            assert len(steps)<=(5 if args.copy_price_cache else 1),steps
            if 'forecast_audit' in steps[-1]:
                assert len(comparisons)==3
                audit=steps[-1]['forecast_audit']
                assert audit['remaining_pairs']==sum((p['variant_range'][1]-p['cursor'])*k for p in snapshot_at_model['partitions'])
            assert all(not s['applied'] for s in steps if s.get('error'))
            assert sorted(span for p in final['partitions'] for span in p['ranges'])==sorted(delivered)
        row=dict(label=label,output_inclusive_seconds=elapsed,first_output_seconds=events[0].completed-began,
            rows=total,parts=len(events),steps=steps,comparisons=comparisons,model_snapshot=snapshot_at_model,
            final_snapshot=final,measurement_steps=measurement_steps,source_ranges=sorted(delivered),max_absolute_difference=float(np.max(np.abs(values-reference))))
        results.append(row);print(json.dumps({key:row[key] for key in ['label','output_inclusive_seconds','first_output_seconds','rows']}),flush=True)
    assert source_identity()==source and all(sha256_file(p)==digest for p,digest in inputs.items())
    report=dict(source_sha256=source,script_sha256=sha256_file(__file__),inputs=inputs,records=results,
        helper_sha256={p:sha256_file(p) for p in ['tests/test_jagwas_scan_work.py','tests/test_jagwas_actual_candidate.py',
            'tests/test_jagwas_candidate.py','tests/test_jagwas_tensor_service.py','tests/test_mechanistic_shapes.py',
            'tests/test_jagwas_writer_service.py','tests/fixtures/jagwas_chunk_geometry.json',
            'benchmarks/direct_jagwas_bounded_execution_20260922.py']},
        admission_seconds=preparation,admission=start,component_evidence=component_evidence,
        price_profile=bound_prices,mounts={name:storage_identity(value) for name,value in [('input',path),('output',destination)]},
        scope='Actual two-GPU native PGEN/JAGWAS source frontier, prepared-window calculation and controller integration. Bounded header construction and all three comparisons run after the first fsynced part, with the issue frontier held and cost charged. Structural admission precedes scan without payload census or runtime ranking. Synthetic prices, explicit 0..1 second in-flight boundary scenario, and generous audit-only budget; no measured prediction accuracy, production default budget compliance or automatic-tuning speedup claim. Per-execution times include output, exclude shared preparation/admission and imports.')
    if args.copy_price_cache:
        report['scope']+=' Optional independent-copy measurement/cache checks run at most two samples per charged productive step. Refresh uses seven samples across four useful-output callbacks; consistent cache reuse uses two fresh samples. Price-bound forecasting follows a later written part. All costs enter cumulative tuning accounting. The measured coefficient does not qualify other synthetic prices.'
        report['helper_sha256']['benchmarks/direct_initial_component_price_20260922.py']=sha256_file('benchmarks/direct_initial_component_price_20260922.py')
        assert execution_context(devices,input_path=path,output_path=destination)==current_context
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
