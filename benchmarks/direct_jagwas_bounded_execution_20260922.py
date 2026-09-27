"""Execute every bounded JAGWAS proposal and independently check its output.

Synthetic prices exercise planning/execution wiring only. This is not a timing
ranking or production capacity experiment. Inputs and metadata live on /data.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import torch
from test_jagwas_candidate_space import specification
from test_pgen_native_reader import write_pgen
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity
from torchgwas.geometry_collection import write_record
from torchgwas.jagwas_candidate_space import bounded_jagwas_plan
from torchgwas.sumstats_indexed import open_indexed_sumstats


def reference(calls,y,c,indices):
    x=np.column_stack([np.ones(len(y)),c.astype(np.float64)])
    values=y.astype(np.float64)
    residual=values-x@np.linalg.lstsq(x,values,rcond=None)[0]
    correlation=np.corrcoef(residual.T)
    statistics=[]
    for index in indices:
        design=np.column_stack([x,calls[index].astype(np.float64)])
        coefficients=np.linalg.lstsq(design,values,rcond=None)[0]
        errors=values-design@coefficients
        variance=np.sum(errors*errors,axis=0)/(len(values)-np.linalg.matrix_rank(design))
        standard_error=np.sqrt(variance*np.linalg.inv(design.T@design)[-1,-1])
        statistics.append(coefficients[-1]/standard_error)
    t=np.asarray(statistics)
    return np.sum(t*np.linalg.solve(correlation,t.T).T,axis=1)


def read_result(path,markers):
    manifest,parts=open_indexed_sumstats(path/'sumstats')
    result=np.full(markers,np.nan,dtype=np.float64)
    seen=np.zeros(markers,dtype=bool)
    for part in parts:
        assert set(part)=={'variant_index','chi2'}
        indices=part['variant_index']
        assert len(np.unique(indices))==len(indices) and not seen[indices].any()
        result[indices]=part['chi2'];seen[indices]=True
    assert seen.all() and np.isfinite(result).all()
    return manifest,result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args()
    fixture=Path(args.fixture);out=Path(args.out)
    fixture.mkdir(parents=True,exist_ok=False);out.mkdir(parents=True,exist_ok=False)
    n,m,k=2049,1025,512
    rng=np.random.default_rng(9221703)
    calls=rng.integers(0,3,size=(m,n),dtype=np.uint8)
    y=rng.normal(size=(n,k)).astype(np.float32)
    covariates=rng.normal(size=(n,2)).astype(np.float32)
    y[:,:4]+=.25*calls[:4].T
    path=fixture/'input.pgen';write_pgen(path,calls)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(
        f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    np.save(fixture/'phenotype.npy',y);np.save(fixture/'covariates.npy',covariates)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    devices=['cuda:0','cuda:1'];limit=2<<30
    for device in devices:
        props=torch.cuda.get_device_properties(device)
        torch.cuda.set_per_process_memory_fraction(limit/props.total_memory,device)
        if torch.cuda.mem_get_info(device)[0]<limit:raise ValueError('Insufficient live memory for audit')
    sources=source_identity();spec=specification(path)
    # Reserves are explicit and the physical allocator also caps each GPU.
    spec['joint'].update(host_memory_bytes=4<<30,host_reserve_bytes=512<<20,
        device_memory_bytes={d:limit for d in devices},device_reserve_bytes=256<<20)
    write_record(out/'request.json',dict(model='detailed_jagwas_space',**spec))
    plan=bounded_jagwas_plan(**spec);write_record(out/'plan.json',plan)
    indices=np.unique(np.r_[np.arange(4),rng.choice(m,29,replace=False),m-1])
    truth=reference(calls,y,covariates,indices)
    np.savez(out/'independent_reference.npz',indices=indices,chi2=truth)
    reports=[];baseline=None
    for candidate in plan['candidates']:
        index=candidate['candidate_index'];directory=out/('candidate_'+str(index))
        options=dict(candidate['api_kwargs']);options.pop('genotype')
        for key,value in candidate['required_environment'].items():
            if os.getenv(key)!=value:raise ValueError('Executor environment differs: '+key)
        # Count scenarios in planning do not guess occupancy from alpha. This
        # execution explicitly retains every valid statistic for verification.
        result=run_linear_gwas(path,fixture/'phenotype.npy',fixture/'covariates.npy',
            output_dir=directory,significance_threshold=1.,**options)
        manifest,values=read_result(directory,m)
        assert manifest['shape']==[m,k] and manifest['df']==k
        np.testing.assert_allclose(values[indices],truth,rtol=3e-4,atol=3e-4)
        if baseline is None:baseline=values.copy()
        np.testing.assert_allclose(values,baseline,rtol=6e-5,atol=3e-4)
        assert result.run_metadata['variant_devices']==candidate['devices']
        rows=dict(candidate_index=index,chunk_size=candidate['chunk_size'],devices=candidate['devices'],
            all_variants_retained=True,independent_cells=len(indices),
            max_absolute_reference_error=float(np.max(np.abs(values[indices]-truth))),
            max_relative_reference_error=float(np.max(np.abs(values[indices]-truth)/np.abs(truth))),
            max_absolute_cross_candidate_error=float(np.max(np.abs(values-baseline))),
            output_sha256=hashlib.sha256(values.tobytes()).hexdigest())
        reports.append(rows);write_record(out/('check_'+str(index)+'.json'),rows)
        print(json.dumps(rows),flush=True)
    assert source_identity()==sources
    write_record(out/'report.json',dict(source_sha256=sources,
        benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dimensions=dict(N=n,M=m,K=k,C=2),seed=9221703,
        source_mount=json.loads(subprocess.check_output(['findmnt','-J','-T',str(path),'-o','SOURCE,FSTYPE,TARGET'],text=True)),
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),selected_candidate=plan['selected']['candidate_index'],
        candidates=reports,cuda_allocator_cap_per_device=limit,
        scope='Numerical execution of all six emitted JAGWAS configurations on synthetic native PGEN input. Synthetic component prices and retained captured geometry test model wiring only; no throughput, ranking accuracy or current-source calibration qualification.'))


if __name__=='__main__':main()
