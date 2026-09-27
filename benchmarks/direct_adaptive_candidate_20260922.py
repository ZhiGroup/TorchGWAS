"""Execute an admitted shared-ring candidate through bounded early switches.

The switching sequence is a correctness control, not a fitted tuning policy.
Source prices/compiled geometry are synthetic or retained test fixtures.
"""
import argparse
import hashlib
import json
from pathlib import Path
import os
import subprocess
import threading
import numpy as np
import torch
from test_pgen_native_reader import write_pgen_mixed
from test_jagwas_actual_candidate import actual_candidate
from direct_jagwas_bounded_execution_20260922 import reference,read_result
from torchgwas.adaptive_candidate import adaptive_candidate_memory
from torchgwas.adaptive_chunks import AlignedChunkSizeControl,InitialChunkMeasurements
from torchgwas.api import load_genotype,run_linear_gwas
from torchgwas.detailed_calibration import source_identity,sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.linear import linear_scan_multigpu
from torchgwas.pgen_work_census import census
from torchgwas.reduce import JagwasReduction
from torchgwas.sumstats_indexed import write_indexed_sumstats


class Window(InitialChunkMeasurements):
    """Only the first six completed sampled chunks choose sizes in this audit."""
    def __init__(self,control):
        super().__init__(['cuda:0','cuda:1'],max_chunks_per_device=8,warmup_chunks=0,
                         stride=1,max_window_seconds=60.)
        self.control=control;self.count=0;self.transitions=[];self.decision_lock=threading.Lock()

    def __call__(self,row):
        super().__call__(row)
        with self.decision_lock:
            self.count+=1
            if self.count in (2,6):
                size=256 if self.count==2 else 512
                self.control.set_size(size)
                self.transitions.append(dict(completed_samples=self.count,size=size,
                    after_device=row.device,after_range=[row.start,row.end]))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();fixture=Path(args.fixture);out=Path(args.out)
    fixture.mkdir(parents=True,exist_ok=False);out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    n,m,k,c=2049,4097,512,2
    rng=np.random.default_rng(9223311)
    calls=rng.integers(0,3,size=(m,n),dtype=np.uint8)
    forms=[0 if i%512==0 else 2 for i in range(m)]
    path=fixture/'input.pgen';write_pgen_mixed(path,calls,forms)
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    y=rng.normal(size=(n,k)).astype(np.float32);cov=rng.normal(size=(n,c)).astype(np.float32)
    y[:,:4]+=.25*calls[:4].T
    np.save(fixture/'phenotype.npy',y);np.save(fixture/'covariates.npy',cov)
    inputs={str(fixture/name):sha256_file(fixture/name) for name in ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}
    sources=source_identity();candidate=actual_candidate(path,512,2);profiles={};cap=2<<30
    for device in candidate['devices']:
        p=torch.cuda.get_device_properties(device)
        profiles[device]=dict(torch_version=torch.__version__,sm_count=p.multi_processor_count,
            max_threads_per_sm=p.max_threads_per_multi_processor,compute_capability=[p.major,p.minor],
            cublas_workspace_config=os.getenv('CUBLAS_WORKSPACE_CONFIG'),cublas_handle_stream_pairs=2)
        assert torch.cuda.mem_get_info(device)[0]>cap
        torch.cuda.set_per_process_memory_fraction(cap/p.total_memory,device)
    fine=census(path,128,include_chunks=True)
    memory=adaptive_candidate_memory(candidate,chunk_sizes=[128,256,512],source_census=fine,
        reduction='jagwas',device_memory_profiles=profiles,host_reserve_bytes=512<<20,device_reserve_bytes=256<<20)
    assert not memory['missing_geometry'] and all(value<cap for value in memory['device_bytes'].values())
    assert memory['host_bytes']<4<<30
    write_record(out/'admission.json',memory)
    print(json.dumps(dict(phase='admitted',device_bytes=memory['device_bytes'],
                         decoder_extra=memory['decoder_extra_bytes_by_device'])),flush=True)
    run_linear_gwas(path,y,cov,pgen_mode='hardcall',reduce='jagwas',compute_dtype='float32',
        chunk_size=512,variant_devices=candidate['devices'],reader_workers=3,prefetch_chunks=2,
        sumstats_queue_depth=2,output_dir=out/'fixed')
    _,fixed=read_result(out/'fixed',m)
    control=AlignedChunkSizeControl([128,256,512],initial=128);window=Window(control)
    source=load_genotype(path,genotype_format='pgen',pgen_mode='hardcall',reader_workers=3)[0]
    for device in candidate['devices']:torch.cuda.reset_peak_memory_stats(device)
    chunks,basis=linear_scan_multigpu(source,y,cov,devices=candidate['devices'],chunk_size=512,
        reader_workers=3,prefetch_chunks=2,compute_dtype='float32',compute_p_values=False,
        ordered=False,shared_queue_depth=2,reduction_factory=JagwasReduction,
        _chunk_size_selector=control,_chunk_observer=window)
    delivered=[]
    def checked():
        for row in chunks:
            delivered.append((int(row[0]),int(row[1])))
            yield row
    write_indexed_sumstats(out/'adaptive'/'sumstats',[f'v{i}' for i in range(m)],[f't{i}' for i in range(k)],
        n,checked(),kind='jagwas',df=n-basis.shape[1]-2,chi2_df=k,store_beta=False,fsync=True)
    _,values=read_result(out/'adaptive',m);snapshot=window.snapshot()
    cursor=0
    for lo,hi in sorted(delivered):
        assert lo==cursor;cursor=hi
    assert cursor==m and len(set(delivered))==len(delivered)
    assert set(hi-lo for lo,hi in delivered)>={128,256,512,1}
    assert [row['size'] for row in window.transitions]==[256,512]
    assert not snapshot['pending'] and len(snapshot['observations'])<=16
    np.testing.assert_allclose(values,fixed,rtol=6e-5,atol=3e-4)
    indices=np.unique(np.r_[np.arange(4),rng.choice(m,29,replace=False),m-1])
    truth=reference(calls,y,cov,indices)
    np.testing.assert_allclose(values[indices],truth,rtol=3e-4,atol=3e-4)
    for path,digest in inputs.items():assert sha256_file(path)==digest
    assert source_identity()==sources
    peaks={d:torch.cuda.max_memory_allocated(d) for d in candidate['devices']}
    assert all(peaks[d]<=memory['device_bytes'][d] for d in peaks)
    write_record(out/'report.json',dict(source_sha256=sources,benchmark_sha256=sha256_file(__file__),
        dimensions=dict(N=n,M=m,K=k,C=c),inputs=inputs,
        source_mount=json.loads(subprocess.check_output(['findmnt','-J','-T',str(fixture/'input.pgen'),'-o','SOURCE,FSTYPE,TARGET'],text=True)),
        transitions=window.transitions,observations=snapshot,source_ranges=sorted(delivered),
        actual_shapes=sorted({hi-lo for lo,hi in delivered}),retained_variants=m,independent_variants=len(indices),
        max_absolute_fixed_difference=float(np.max(np.abs(values-fixed))),
        max_absolute_reference_error=float(np.max(np.abs(values[indices]-truth))),
        device_peak_allocated_bytes=peaks,modeled_device_bytes=memory['device_bytes'],
        output_sha256=hashlib.sha256(values.tobytes()).hexdigest(),
        scope='Fixed-ring memory and numerical audit of bounded early chunk transitions, full-panel JAGWAS and native LD replay. The transition sequence is a test control, not autotuning/ranking qualification; source accounting plus explicit reserves is not a general allocator guarantee.'))
    print(json.dumps(dict(phase='complete',retained_variants=m,shapes=sorted({hi-lo for lo,hi in delivered}),
                         transitions=window.transitions)),flush=True)


if __name__=='__main__':main()
