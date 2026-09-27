"""Real public JIT copy refresh/reuse/expiry with unchanged JAGWAS output.

Only resident-copy CPU service is measured. Other prices are explicit synthetic
controls. This audit verifies lifecycle and costs, not tuning accuracy/speedup.
"""
import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import time

import numpy as np
import torch

from direct_jagwas_bounded_execution_20260922 import read_result
from test_jagwas_actual_candidate import CAPTURE,writer_prices
from test_jagwas_preparation import profile as preparation_profile
from test_jagwas_scan_work import joint_fixture
from test_pgen_native_reader import write_pgen
from torchgwas.api import run_linear_gwas
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import (bind_detailed_profile,execution_context,
    source_identity,sha256_file,write_detailed_profile,read_detailed_profile)
from torchgwas.tensor_service import host_primitive_name
from torchgwas.tensor_work import eager_statistics_work


def main():
    parser=argparse.ArgumentParser()
    for name in ('fixture','out','data'):parser.add_argument('--'+name,required=True)
    parser.add_argument('--digest-comparison',action='store_true',
        help='Run three counterbalanced full-byte/cached-digest pairs after a refresh.')
    args=parser.parse_args();fixture=Path(args.fixture);out=Path(args.out);data=Path(args.data)
    out.mkdir(parents=True,exist_ok=False);data.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    n,m,k=2049,16385,512;devices=['cuda:0','cuda:1'];source=source_identity()
    y=np.load(fixture/'phenotype.npy');c=np.load(fixture/'covariates.npy')
    calls=np.tile((np.arange(n)%3).astype(np.uint8),(m,1))
    for row in range(m):calls[row,row%n]=(calls[row,row%n]+1)%3
    path=data/'input.pgen';write_pgen(path,calls);del calls
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
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
        entry['owned_result_copy_scenario']['resident_cpu_seconds_per_byte']=1e-10
        profiles[device]=entry
        memory[device]=dict(torch_version=torch.__version__,sm_count=props.multi_processor_count,
            max_threads_per_sm=props.max_threads_per_multi_processor,compute_capability=[props.major,props.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
    context=dict(name='public-refresh-audit',devices=devices,profiles=profiles,shared_capacities=caps)
    execution=execution_context(devices,input_path=path,output_path=data)
    cache=CalibrationParameterCache(out/'seed')
    deps=dict(source_sha256=source,execution_context=execution,measurement_protocol={'operation':'synthetic seed control'})
    seed=cache.store('cpu_capacity','seed-copy',dict(rate=1e-10),dependencies=deps,
        provenance={'scope':'Synthetic coefficient, intentionally expired by a short caller lifetime'},
        observed_unix_seconds=time.time(),max_age_seconds=3600.)
    writer=cache.store('cpu_capacity','jagwas_components',
        dict(writer_prices=writer_prices(),preparation_services={'unused_during_prepared_windows':{}}),
        dependencies=dict(source_sha256=source,execution_context=execution),
        provenance={'scope':'Synthetic indexed-writer controls'},observed_unix_seconds=time.time(),max_age_seconds=3600.)
    binding=dict(artifact=seed['path'],kind='cpu_capacity',name='seed-copy',dependencies=deps,max_age_seconds=None,
        targets=[dict(context_path=[0,'profiles',d,'owned_result_copy_scenario','resident_cpu_seconds_per_byte'],
                      value_path=['rate']) for d in devices])
    profile=bind_detailed_profile([context],execution,sources=source,limitations=[__doc__],
        component_artifacts={r['path']:sha256_file(r['path']) for r in (seed,writer)},price_bindings=[binding])
    profile['price_bindings'][0]['max_age_seconds']=1e-6
    initial=out/'initial_profile.json';write_detailed_profile(profile,initial)
    config=dict(bounds=dict(chunks=[128,256]),qc_trait_block=128,jagwas_services=writer['path'],
        joint=dict(cpu_workers=3,host_memory_bytes=4<<30,device_memory_bytes={d:2<<30 for d in devices},
            host_reserve_bytes=512<<20,device_reserve_bytes=256<<20,device_memory_profiles=memory,
            host_scenarios={'specified':dict(host_serial_fraction=.5,host_serial_policy='fluid')},
            occupancy_scenarios={'all':'dense'}),
        initial_chunks=dict(context=context['name'],chunk_size=128,partition_axis='variant',trait_block=None,
            window_markers=[256,512,768],budget=dict(max_steps=5,max_cpu_seconds=10.,max_window_seconds=60.),
            cost_forecasts=dict(remaining_seconds=100.,expected_cpu_seconds=.01,expected_wall_seconds=.02,
                switching_seconds=0.,publication_seconds=0.,reserve_seconds=.01),
            resident_copy_refresh=dict(cache_dir=str(out/'parameters'),profile_dir=str(out/'profiles'),
                binding_indexes=[0],max_age_seconds=30.,expected_cpu_seconds=.04,expected_wall_seconds=.1),
            forecast_options=dict(boundary_adjustments=dict(baseline=[0.,1.],candidate=[0.,1.]),
                relative_model_error=.05,max_slope_change=.1,max_extrapolation=32.,
                assumptions=dict(source_work='Header-bound native PGEN synthetic input.',output_occupancy='All valid joint rows.',
                    resource_capacity='Resident copy measured during useful chunks; other rates synthetic.',
                    partition_balance='Exact unissued partition extents with core/envelope scenarios.'))))
    if args.digest_comparison:
        config['initial_chunks']['resident_copy_refresh']['max_age_seconds']=3600.
    (out/'config.json').write_text(json.dumps(config,indent=2))
    original_files={str(initial):sha256_file(initial),seed['path']:sha256_file(seed['path'])}
    report=dict(scope=__doc__,source_sha256=source,benchmark_sha256=sha256_file(__file__),config=config,
        input_sha256={str(p):sha256_file(p) for p in [path,path.with_suffix('.pvar'),path.with_suffix('.psam'),fixture/'phenotype.npy',fixture/'covariates.npy']},runs=[],
        coverage=dict(measured=False,reused=False,drift=False,expiry=False,rejected_unstable=False))
    active=initial;truth=None;latest_record=None
    try:
        labels=(('control','refresh','check_strict_1','check_cached_1','check_cached_2','check_strict_2',
            'check_strict_3','check_cached_3') if args.digest_comparison else
            ('control','refresh','check','check_again','expired'))
        for label in labels:
            if label=='expired' and latest_record is not None:
                wait=max(0.,latest_record['record']['observed_unix_seconds']+latest_record['record']['max_age_seconds']-time.time()+.05)
                print('WAIT original measurement expiry',wait,flush=True)
                time.sleep(wait)
            run_config=deepcopy(config)
            run_config['initial_chunks']['reuse_binding_digests']='strict' not in label
            options=(dict(chunk_size=128,device=devices[0],variant_devices=devices,reader_workers=3,prefetch_chunks=2)
                if label=='control' else dict(autotune_profile=active,autotune_config=run_config))
            wall=time.time();began=time.perf_counter()
            result=run_linear_gwas(path,y,c,pgen_mode='hardcall',compute_dtype='float32',reduce='jagwas',
                output_dir=data/label,sumstats_fields='t',sumstats_queue_depth=2,**options)
            elapsed=time.perf_counter()-began;_,values=read_result(data/label,m)
            if truth is None:truth=values.copy()
            np.testing.assert_allclose(values,truth,rtol=6e-5,atol=3e-4)
            audit=result.run_metadata.get('autotune')
            row=dict(label=label,api_seconds=elapsed,max_absolute_difference=float(np.max(np.abs(values-truth))),autotune=audit)
            report['runs'].append(row)
            if audit is not None:
                state=audit['productive'];refresh=state['resident_copy_refresh'];observed=refresh['controller']['result']
                assert state['finished']['successful'] and refresh['buffers_released']
                assert all(p['cursor']==p['variant_range'][1] for p in state['partitions'])
                row['first_written_after_api_seconds']=state['first_written']-began
                if refresh['pending']:
                    # A noisy window is a valid fail-closed outcome. Keep this
                    # evidence instead of widening thresholds or hiding a run.
                    assert refresh['published'] is None and state['current_chunk_size']==128
                    assert all(not decision.get('applied',False) for decision in state['decisions'])
                    row['outcome']=refresh['controller']['state']
                    report['coverage']['rejected_unstable']|=row['outcome']=='unstable_measurement'
                    assert original_files=={p:sha256_file(p) for p in original_files}
                    print('COMPLETE',label,elapsed,row['outcome'],flush=True)
                    continue
                assert observed is not None
                assert state['forecast_attempts'] and not state['decisions'][-1].get('error')
                assert audit['selected']['memory']['resident_copy_refresh_bytes']==32<<20
                row['outcome']=observed['status']
                if label.startswith('check') and observed['status']=='reused_original':
                    report['coverage']['reused']=True
                    assert observed['record_sha256']==latest_record['record_sha256']
                    assert observed['record']['observed_unix_seconds']==latest_record['record']['observed_unix_seconds']
                    assert len(refresh['controller']['samples'])==2
                else:
                    assert [b['samples'] for b in refresh['batches']]==[2,4,6,7]
                    assert observed['status']=='measured_and_published'
                    report['coverage']['measured']=True
                    assert observed['record']['observed_unix_seconds']>=wall+row['first_written_after_api_seconds']-.01
                    if latest_record is not None:assert observed['record_sha256']!=latest_record['record_sha256']
                    if label.startswith('check') and refresh['controller']['check'] is not None:
                        assert refresh['controller']['check']['status'] in ('drift','unstable_check','expired_or_invalid_during_check')
                        report['coverage']['drift']|=refresh['controller']['check']['status']=='drift'
                    if label=='expired' and latest_record is not None:report['coverage']['expiry']=True
                latest_record=deepcopy(observed)
                active=Path(refresh['published']['path'])
                assert read_detailed_profile(active)['source_sha256']==source
                for artifact in (str(active),observed['path']):original_files.setdefault(artifact,sha256_file(artifact))
                assert original_files=={p:sha256_file(p) for p in original_files}
            print('COMPLETE',label,elapsed,flush=True)
        assert source_identity()==source
        report['complete']=True
    finally:
        report['preserved_artifacts']=original_files
        (out/'report.json').write_text(json.dumps(report,indent=2))


if __name__=='__main__':main()
