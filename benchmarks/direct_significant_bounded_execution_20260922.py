"""Numerical audit of bounded significant-pairs proposals, including both tails.

All service prices below are synthetic controls, not measured capacities.
Fresh CUDA captures retain launch geometry only. This verifies execution
wiring and statistical invariance, not runtime ranking or throughput.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np
from scipy.special import stdtrit
import torch
from test_pgen_native_reader import write_pgen
from test_mechanistic_shapes import fixture as scan_fixture
from test_significant_host_model import bank as selector_fixture
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import _gpu_geometry, write_record
from torchgwas.native_control_work import native_control_work
from torchgwas.significant_candidate_space import prepare_significant_host_candidates, bounded_significant_host_plan
from torchgwas.significant_host_work import significant_host_memory
from torchgwas.sumstats_indexed import open_indexed_sumstats
from torchgwas.tensor_work import eager_statistics_work
from torchgwas.tensor_service import host_primitive_name


def synthetic_profile(n, c, workers, device):
    _, profile = scan_fixture(c=c)
    profile.update(depth=2, decode_workers=workers, validate_range=True,
        event_wait_cpu_fraction=1., result_ownership='borrowed',
        write_bytes_per_second=1e8, fsync_seconds=1e-6,
        pin_cpu_seconds_per_page=1e-6, pin_driver_seconds_per_page=1e-6,
        pin_cached_cpu_seconds_per_call=1e-7,
        writeback_service=dict(pagecache_seconds_per_byte=1e-9, storage_seconds_per_byte=1e-8,
            submit_seconds=1e-6, wait_seconds=1e-6, fadvise_seconds=1e-6),
        control_primitives={name:1e-6 for counts in native_control_work().values() for name in counts},
        result_finish_service=dict(cpu_seconds=1e-6, serial_cpu_seconds=0., baseline_copy_bytes=0,
            replaces_fixed_finish_and_tensor_conversion=True, includes_ready_cuda_event=False),
        setup_primitives={name:dict(reference_shape=[max(32,c+3),1,c],cpu_seconds=1e-6,non_cpu_seconds=1e-6)
            for name in ['residual_common','residual_block','design_common','design_block']},
        kernel_geometry=[])
    profile['gpu_resources'].update(gpu_fraction=1.,sm_count=torch.cuda.get_device_properties(device).multi_processor_count)
    profile['process_units'].update(bytearray_zero_bytes=1e-10,covariate_basis_work=1e-10)
    profile['decode_units']['expand_int8_tail_sample']=1e-9
    return profile


def reference(calls, y, covariates, alpha):
    # Independent FP64 Frisch-Waugh-Lovell OLS, with no torchGWAS statistics.
    x=np.column_stack([np.ones(len(y)),covariates.astype(np.float64)])
    yr=y.astype(np.float64);yr-=x@np.linalg.lstsq(x,yr,rcond=None)[0]
    g=calls.T.astype(np.float64);g-=x@np.linalg.lstsq(x,g,rcond=None)[0]
    ss=np.sum(g*g,axis=0)[:,None]
    score=g.T@yr
    beta=score/ss
    df=len(y)-x.shape[1]-1
    t=score/np.sqrt(ss*(np.sum(yr*yr,axis=0)[None,:]-score*beta)/df)
    critical=float(stdtrit(df,1-alpha/2))
    vi,ti=np.nonzero(np.abs(t)>=critical)
    # API beta is per residual phenotype population standard deviation. The
    # independent raw-phenotype OLS t statistic and pair identities are invariant.
    return dict(variant_index=vi,trait_index=ti,beta=beta[vi,ti]/np.std(yr,axis=0)[ti],t_stat=t[vi,ti],
        df=np.full(len(vi),df,dtype=np.int32)),float(np.min(np.abs(np.abs(t)-critical)))


def read_result(path, markers, traits):
    manifest,parts=open_indexed_sumstats(path/'sumstats')
    fields=['variant_index','trait_index','beta','t_stat','df']
    collected={field:[] for field in fields}
    for part in parts:
        assert set(part)==set(fields)
        for field in fields:collected[field].append(np.asarray(part[field]).copy())
    result={field:np.concatenate(rows) for field,rows in collected.items()}
    vi,ti=result['variant_index'],result['trait_index']
    assert np.all((vi>=0)&(vi<markers)) and np.all((ti>=0)&(ti<traits))
    keys=vi*traits+ti
    assert len(np.unique(keys))==len(keys)
    order=np.argsort(keys)
    return manifest,{field:values[order] for field,values in result.items()}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();fixture=Path(args.fixture);out=Path(args.out)
    fixture.mkdir(parents=True,exist_ok=False);out.mkdir(parents=True,exist_ok=False)
    n,m,k,c=2049,1025,4097,2
    seed=9222207;alpha=.001;rng=np.random.default_rng(seed)
    calls=rng.integers(0,3,size=(m,n),dtype=np.uint8)
    y=rng.normal(size=(n,k)).astype(np.float32)
    covariates=rng.normal(size=(n,c)).astype(np.float32)
    y[:,:4]+=.5*calls[:4].T;y[:,-1]+=.5*calls[-1]
    path=fixture/'input.pgen';write_pgen(path,calls)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(
        f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    np.save(fixture/'phenotype.npy',y);np.save(fixture/'covariates.npy',covariates)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    sources=source_identity();contexts=[]
    for count in (1,2):
        devices=[f'cuda:{i}' for i in range(count)]
        contexts.append(dict(name=f'gpus{count}',devices=devices,
            profiles={d:synthetic_profile(n,c,2//count,d) for d in devices},
            shared_capacities=dict(cpu=2.,dram=1e9,input=1e8,output=1e8)))
    workload=dict(genotype=str(path),samples=n,markers=m,traits=k,covariates=c,
        matching_sample_order=True,complete_phenotypes=True,phenotype_c_contiguous=True)
    bounds=dict(chunks=[128,256],trait_blocks=[512,1024,k],max_candidates=12)
    output=dict(block_bytes=None,queue_depth=1,store_beta=True,fsync=True)
    space=prepare_significant_host_candidates(workload,contexts,output=output,**bounds)
    reserve=64<<20
    memories=[significant_host_memory(candidate,device_reserve_bytes=reserve) for candidate in space['candidates']]
    small=max(max(memory['device_bytes'].values()) for candidate,memory in zip(space['candidates'],memories) if candidate['trait_block']!=k)
    large=min(max(memory['device_bytes'].values()) for candidate,memory in zip(space['candidates'],memories) if candidate['trait_block']==k)
    assert small<large
    limit=(small+large)//2
    # A deliberately constrained budget forces actual phenotype partitioning.
    # An equal physical allocator cap is applied during capture and execution.
    for d in contexts[-1]['devices']:
        assert torch.cuda.mem_get_info(d)[0]>limit
        torch.cuda.set_per_process_memory_fraction(limit/torch.cuda.get_device_properties(d).total_memory,d)
    requests={}
    for candidate,memory in zip(space['candidates'],memories):
        if any(value>limit for value in memory['device_bytes'].values()):continue
        for tile in candidate['tiles']:
            width=tile['trait_range'][1]-tile['trait_range'][0]
            chunk=tile['profile']['chunk_markers'];device=tile['device']
            for b in sorted({chunk,m%chunk}):
                key=(device,n,b,width,c)
                requests[key]=dict(device=device,shape=[n,b,width,c],
                    source_memory_bound_bytes=memory['device_bytes'][device],device_budget_bytes=limit)
    print(json.dumps(dict(stage='geometry',unique_shapes=len(requests),device_budget_bytes=limit,
        largest_tiled_bound_bytes=small,smallest_full_panel_bound_bytes=large)),flush=True)
    captured=_gpu_geometry(list(requests.values()),out);write_record(out/'geometry.json',captured)
    for context in contexts:
        for device,profile in context['profiles'].items():
            rows=[{name:value for name,value in row.items() if name!='device'} for row in captured['rows'] if row['device']==device]
            profile['kernel_geometry']=copy.deepcopy(rows)
            for row in rows:
                work=eager_statistics_work(n,row['B'],row['K'],c,True)
                profile['host_primitives'].update({host_primitive_name(call):1e-6 for call in work['host_calls']})
    joint=dict(occupancy_scenarios={'none':'empty','all':'dense'},
        host_scenarios={'fluid':dict(host_serial_fraction=0.)},cpu_workers=2,
        host_memory_bytes=4<<30,host_reserve_bytes=512<<20,
        device_memory_bytes={d:limit for d in contexts[-1]['devices']},device_reserve_bytes=reserve)
    spec=dict(workload=workload,contexts=contexts,bounds=bounds,output=output,joint=joint,
        prices=selector_fixture(),significance_threshold=alpha)
    write_record(out/'request.json',dict(model='detailed_significant_host_space',**spec))
    print(json.dumps(dict(stage='planning')),flush=True)
    plan=bounded_significant_host_plan(**spec);write_record(out/'plan.json',plan)
    assert len(plan['candidates'])==8 and len(plan['rejected'])==2
    assert all(row['reason']=='device_memory' for row in plan['rejected'])
    assert all(row['trait_block']<k for row in plan['candidates'])
    truth,margin=reference(calls,y,covariates,alpha);np.savez(out/'reference.npz',**truth)
    assert len(truth['variant_index'])>100
    assert ((truth['variant_index']==m-1)&(truth['trait_index']==k-1)).any()
    reports=[];baseline=None
    for row in plan['candidates']:
        index=row['candidate_index'];directory=out/('candidate_'+str(index))
        options=dict(row['api_kwargs']);options.pop('genotype')
        for name,value in row['required_environment'].items():
            if os.getenv(name)!=value:raise ValueError('Executor environment differs: '+name)
        print(json.dumps(dict(stage='execute',candidate_index=index,chunk_size=row['chunk_size'],
            trait_block=row['trait_block'],devices=row['devices'])),flush=True)
        result=run_linear_gwas(path,fixture/'phenotype.npy',fixture/'covariates.npy',output_dir=directory,**options)
        manifest,values=read_result(directory,m,k)
        assert manifest['shape']==[m,k]
        for field in ['variant_index','trait_index','df']:
            np.testing.assert_array_equal(values[field],truth[field])
        for field in ['beta','t_stat']:
            np.testing.assert_allclose(values[field],truth[field],rtol=3e-5,atol=3e-5)
        if baseline is None:baseline=values
        for field in ['beta','t_stat']:
            np.testing.assert_allclose(values[field],baseline[field],rtol=3e-5,atol=3e-5)
        assert result.run_metadata['trait_devices']==row['devices']
        check=dict(candidate_index=index,chunk_size=row['chunk_size'],trait_block=row['trait_block'],
            devices=row['devices'],retained_pairs=len(values['variant_index']),
            exact_pair_set_matches=True,both_axis_tail_pair_retained=True,
            max_absolute_t_error=float(np.max(np.abs(values['t_stat']-truth['t_stat']))),
            max_absolute_beta_error=float(np.max(np.abs(values['beta']-truth['beta']))),
            max_cross_candidate_t_error=float(np.max(np.abs(values['t_stat']-baseline['t_stat']))),
            pair_sha256=hashlib.sha256(np.column_stack([values['variant_index'],values['trait_index']]).tobytes()).hexdigest())
        reports.append(check);write_record(out/('check_'+str(index)+'.json'),check)
        print(json.dumps(check),flush=True)
    assert source_identity()==sources
    write_record(out/'report.json',dict(source_sha256=sources,
        benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dimensions=dict(N=n,M=m,K=k,C=c),seed=seed,significance_threshold=alpha,
        reference_threshold_margin=margin,
        source_mount=json.loads(subprocess.check_output(['findmnt','-J','-T',str(path),'-o','SOURCE,FSTYPE,TARGET'],text=True)),
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),selected_candidate=plan['selected']['candidate_index'],
        cuda_allocator_cap_per_device=limit,largest_tiled_bound_bytes=small,
        smallest_full_panel_bound_bytes=large,rejected_full_panel_candidates=len(plan['rejected']),candidates=reports,
        scope='Numerical execution of every feasible bounded significant-host configuration with fresh duration-free geometry, nonempty output and explicit memory pressure. Synthetic prices test planning wiring only; no throughput, ranking accuracy or capacity qualification.'))


if __name__=='__main__':main()
