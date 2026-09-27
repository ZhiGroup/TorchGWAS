"""Public deferred-calculator audit with real output, not a speedup benchmark.

Only resident-copy CPU service is measured. Other component prices remain
explicit synthetic controls, so this does not qualify production decisions.
"""
import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import time

import numpy as np
import torch

from direct_initial_component_price_20260922 import CopyPriceProbe
from direct_jagwas_bounded_execution_20260922 import read_result
from test_jagwas_actual_candidate import CAPTURE,writer_prices
from test_jagwas_preparation import profile as preparation_profile
from test_jagwas_scan_work import joint_fixture
from torchgwas.api import run_linear_gwas
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import execution_context,source_identity,sha256_file,write_detailed_profile
from torchgwas.tensor_service import host_primitive_name
from torchgwas.tensor_work import eager_statistics_work


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fixture',required=True)
    parser.add_argument('--out',required=True);parser.add_argument('--output-data',required=True)
    parser.add_argument('--structural-cache-dir')
    parser.add_argument('--stage-observations',action='store_true')
    parser.add_argument('--planning-cost-history',action='store_true')
    parser.add_argument('--background-planning',action='store_true')
    args=parser.parse_args();fixture=Path(args.fixture);out=Path(args.out);destination=Path(args.output_data)
    out.mkdir(parents=True,exist_ok=False);destination.mkdir(parents=True,exist_ok=False)
    preparation_started=time.perf_counter()
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    n,m,k=2049,4097,512;devices=['cuda:0','cuda:1'];source=source_identity()
    y=np.load(fixture/'phenotype.npy');cov=np.load(fixture/'covariates.npy')
    _,template=joint_fixture();profiles={};memory={};caps=dict(cpu=3.,dram=1e9,input=1e8,output=1e8)
    for index,device in enumerate(devices):
        entry=preparation_profile(template);props=torch.cuda.get_device_properties(device)
        entry.update(return_beta=True,result_ownership='owned',validate_range=True,event_wait_cpu_fraction=1.,
            decode_workers=2-index,cpu_available_cores=3.,chunk_markers=256,
            write_bytes_per_second=1e8,fsync_seconds=1e-6,
            writeback_service=dict(pagecache_seconds_per_byte=1e-9,storage_seconds_per_byte=1e-8),
            kernel_geometry=deepcopy(CAPTURE['statistics']),joint_kernel_geometry=deepcopy(CAPTURE['projection']),
            host_primitives={host_primitive_name(call):1e-6 for call in eager_statistics_work(n,128,k,2,True)['host_calls']})
        entry['gpu_resources']['sm_count']=props.multi_processor_count
        entry['decode_units'].update({name:1e-8 for name in ['uleb1','uleb2','uleb3','uleb4','uleb5',
            'set_category','difflist_group_absolute_id','difflist_category_extract','difflist_record_header']})
        profiles[device]=entry
        memory[device]=dict(torch_version=torch.__version__,sm_count=props.multi_processor_count,
            max_threads_per_sm=props.max_threads_per_multi_processor,compute_capability=[props.major,props.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
    context=dict(name='public-jit-audit',devices=devices,profiles=profiles,shared_capacities=caps)
    runtime=execution_context(devices,input_path=fixture/'input.pgen',output_path=destination)
    probe=CopyPriceProbe(out/'measurements',[context],source,runtime,max_age_seconds=3600.)
    profile=None
    try:
        for _ in range(4):
            profile,evidence=probe.advance()
            if profile is not None:break
    finally:probe.close()
    assert profile is not None,evidence
    record=CalibrationParameterCache(out/'measurements').store('cpu_capacity','jagwas_components',
        dict(writer_prices=writer_prices(),preparation_services={'unused_during_prepared_windows':{}}),
        dependencies=dict(source_sha256=source,execution_context=runtime),
        provenance={'scope':'Synthetic indexed-writer controls. No observed GWAS timing.'},
        observed_unix_seconds=time.time(),max_age_seconds=3600.)
    profile['component_artifacts'][str(Path(record['path']).resolve())]=sha256_file(record['path'])
    write_detailed_profile(profile,out/'profile.json')
    config=dict(bounds=dict(chunks=[128,256]),qc_trait_block=128,jagwas_services=record['path'],
        joint=dict(cpu_workers=3,host_memory_bytes=4<<30,device_memory_bytes={d:2<<30 for d in devices},
            host_reserve_bytes=512<<20,device_reserve_bytes=256<<20,device_memory_profiles=memory,
            host_scenarios={'specified':dict(host_serial_fraction=.5,host_serial_policy='fluid')},
            occupancy_scenarios={'all':'dense'}),
        initial_chunks=dict(context=context['name'],chunk_size=128,partition_axis='variant',trait_block=None,
            window_markers=[256,512,768],budget=dict(max_steps=1,max_cpu_seconds=10.,max_window_seconds=30.),
            cost_forecasts=dict(remaining_seconds=10.,expected_cpu_seconds=.01,expected_wall_seconds=.02,
                switching_seconds=0.,publication_seconds=0.,reserve_seconds=.01),
            forecast_options=dict(boundary_adjustments=dict(baseline=[0.,1.],candidate=[0.,1.]),
                relative_model_error=.05,max_slope_change=.1,max_extrapolation=8.,
                assumptions=dict(source_work='Fixed synthetic fixture; header-bound source work scenario.',
                    output_occupancy='All valid JAGWAS rows retained.',resource_capacity='One measured copy coefficient; other rates synthetic.',
                    partition_balance='Balanced core/envelope scenarios over exact unissued extents.'))))
    if args.structural_cache_dir:
        config['initial_chunks']['structural_cache_dir']=args.structural_cache_dir
    if args.background_planning:
        config['initial_chunks']['background_planning']=True
    if args.stage_observations:
        config['initial_chunks']['stage_observations']=dict(max_chunks_per_device=2,
            warmup_chunks=1,stride=1,max_window_seconds=20.,cuda_events=True,
            measurement_reserve_seconds=.01)
    if args.planning_cost_history:
        config['initial_chunks']['planning_cost_history']=dict(
            cache_dir=str(out/'planning_cost_history'),max_age_seconds=3600.,
            publication_seconds=.01)
    (out/'config.json').write_text(json.dumps(config,indent=2))
    preparation_seconds=time.perf_counter()-preparation_started
    original_profile=(out/'profile.json').read_bytes()
    original_measurement=Path(evidence['result']['path']).read_bytes()
    prior_binding=None
    results=[];truth=None
    for label in ('control','deferred','reuse'):
        options=(dict(chunk_size=128,device='cuda:0',variant_devices=devices,reader_workers=3,prefetch_chunks=2)
            if label=='control' else dict(autotune_profile=out/'profile.json',autotune_config=config))
        began=time.perf_counter()
        result=run_linear_gwas(fixture/'input.pgen',y,cov,pgen_mode='hardcall',compute_dtype='float32',
            reduce='jagwas',output_dir=destination/label,sumstats_fields='t',sumstats_queue_depth=2,**options)
        elapsed=time.perf_counter()-began
        _,values=read_result(destination/label,m)
        if truth is None:truth=values.copy()
        np.testing.assert_allclose(values,truth,rtol=6e-5,atol=3e-4)
        audit=result.run_metadata.get('autotune')
        if audit is not None:
            if args.structural_cache_dir:
                if label=='deferred':
                    assert audit['selected']['memory']['retained_index_bases_bytes']==0
                else:
                    assert audit['selected']['memory']['retained_index_bases_bytes']>0
            else:
                assert audit['selected']['memory']['retained_index_bases_bytes']==0
            assert audit['candidates_evaluated']==0 and audit['structural_work']['census_passes']==0
            assert audit['productive']['finished']['successful'] is True
            if args.planning_cost_history:
                history=audit['productive']['planning_cost_history']
                assert history['lookup']['hit']==(label=='reuse')
                assert history['publication']['status']=='stored'
                assert len(history['observations'])==1
                if label=='reuse':
                    assert history['prior']['observed_unix_seconds']==prior_cost_observed
                    assert audit['productive']['planning']['steps'][0]['expected_wall_seconds']>=prior_cost_wall
                else:
                    prior_cost_observed=history['observations'][0]['observed_unix_seconds']
                    prior_cost_wall=history['observations'][0]['wall_seconds']
            if args.stage_observations:
                stage=audit['productive']['stage_sample']
                assert not stage['pending'] and len(stage['observations'])==4
                assert {row['device'] for row in stage['observations']}==set(devices)
            if args.structural_cache_dir:
                expected='miss' if label=='deferred' else 'hit'
                assert audit['structural_work']['admission_cache']['status']==expected
                assert audit['productive']['admission_cache']['write_status']==('stored' if label=='deferred' else 'not_attempted')
            assert audit['productive']['planning']['steps']
            assert audit['productive']['decisions'][0]['error'] is None
            assert len(audit['productive']['forecast_attempts'])>=1
            assert all(p['cursor']==p['variant_range'][1] for p in audit['productive']['partitions'])
            bound=audit['productive']['forecast_attempts'][0]['scenarios'][0]['price_evidence']['bindings'][0]
            assert bound['record_sha256']==evidence['result']['record_sha256']
            assert bound['observed_unix_seconds']==evidence['result']['record']['observed_unix_seconds']
            if prior_binding is not None:
                assert bound['age_seconds']>prior_binding['age_seconds']
            prior_binding=bound
        assert (out/'profile.json').read_bytes()==original_profile
        assert Path(evidence['result']['path']).read_bytes()==original_measurement
        first=None if audit is None else audit['productive']['first_written']-began
        results.append(dict(label=label,api_seconds=elapsed,first_written_after_api_seconds=first,
            autotune=audit,max_absolute_difference=float(np.max(np.abs(values-truth)))))
        print('COMPLETE',label,elapsed,flush=True)
    assert source_identity()==source
    (out/'report.json').write_text(json.dumps(dict(results=results,source_sha256=source,
        measured_copy_evidence=evidence,profile=profile,config=config,
        preparation_seconds=preparation_seconds,benchmark_sha256=sha256_file(__file__),
        scope=__doc__,inputs={str(fixture/name):sha256_file(fixture/name) for name in
            ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}),indent=2))


if __name__=='__main__':main()
