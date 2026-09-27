"""Two-job audit of bounded real-chunk measurements and immutable cache reuse.

Synthetic numerical fixture in native PGEN format. This is a correctness and
cache-lifecycle audit, not a workload-ranking or throughput benchmark.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import torch

from torchgwas.adaptive_chunks import InitialChunkMeasurements
from torchgwas.api import load_genotype
from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import execution_context, source_identity, sha256_file
from torchgwas.geometry_collection import write_record
from torchgwas.linear import linear_scan_multigpu, linear_scan
from torchgwas.reduce import JagwasReduction
from torchgwas.reduction_tensor_work import jagwas_tensor_work
from test_pgen_native_reader import write_pgen


def prepare(root):
    root.mkdir(parents=True, exist_ok=True)
    if (root/'manifest.json').exists():
        raise ValueError('Use an existing fixture without --prepare')
    n, m, k = 2049, 32768, 129
    rng = np.random.default_rng(92311)
    categories = rng.integers(0, 3, (m, n), dtype=np.uint8)
    y = rng.normal(size=(n, k)).astype(np.float32)
    cov = rng.normal(size=(n, 3)).astype(np.float32)
    path = root/'input.pgen'
    write_pgen(path, categories)
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n'+''.join(f'1\t{i+1}\tv{i}\tA\tC\n' for i in range(m)))
    path.with_suffix('.psam').write_text('#IID\n'+''.join(f's{i}\n' for i in range(n)))
    np.save(root/'phenotype.npy', y); np.save(root/'covariates.npy', cov)
    indices = np.linspace(0, m-1, 41, dtype=np.int64)
    np.savez(root/'reference_cells.npz', indices=indices, calls=categories[indices].T)
    write_record(root/'manifest.json', dict(samples=n, markers=m, traits=k, covariates=3,
        input_sha256=sha256_file(path), seed=92311,
        scope='Synthetic native-PGEN correctness fixture, not production-scale timing evidence'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--cache', required=True)
    parser.add_argument('--prepare', action='store_true')
    args=parser.parse_args()
    root=Path(args.fixture); out=Path(args.out); out.mkdir(parents=True, exist_ok=True)
    if args.prepare: prepare(root)
    fixture=json.loads((root/'manifest.json').read_text())
    path=root/'input.pgen'
    if sha256_file(path)!=fixture['input_sha256']: raise ValueError('Fixture changed')
    devices=['cuda:0','cuda:1']; chunk=128
    context=execution_context(devices,input_path=path,output_path=out)
    cache=CalibrationParameterCache(args.cache)
    # A duration-free source ledger is reusable until its source/shape changes.
    source=source_identity()
    structural_dependencies=dict(source_sha256={'reduce.py':source['reduce.py'],
        'reduction_tensor_work.py':source['reduction_tensor_work.py']},
        shape=[fixture['samples'],chunk,fixture['traits']],compute_dtype='float32',
        torch_version=context['torch_version'],protocol='jagwas_meta_reduce_v1')
    structural=cache.lookup('source_work','jagwas_reduce',dependencies=structural_dependencies)
    structural_publication=None
    if not structural['hit']:
        work=jagwas_tensor_work(fixture['samples'],chunk,fixture['traits'],phase='reduce')
        structural_publication=cache.store('source_work','jagwas_reduce',work,
            dependencies=structural_dependencies,provenance=dict(job=str(out),source='duration-free meta execution'))
    dependencies=dict(source_sha256=source,execution_context=context,
        workload={k:fixture[k] for k in ('samples','markers','traits','covariates')},
        input_sha256=fixture['input_sha256'],devices=devices,chunk_size=chunk,
        reader_workers=4,depth=2,reduction='jagwas',protocol='early_chunk_cuda_events_v1')
    previous=cache.lookup('stage_observations','initial_chunks',dependencies=dependencies,max_age_seconds=300.)
    # Even a fresh empirical hit gets early-job observations. These are stored
    # separately from the original fuller window rather than replacing it.
    window=InitialChunkMeasurements(devices,max_chunks_per_device=2 if previous['hit'] else 8,
        warmup_chunks=2,stride=4,max_window_seconds=60.)
    y=np.load(root/'phenotype.npy',mmap_mode='r'); cov=np.load(root/'covariates.npy')
    source_geno=load_genotype(path,genotype_format='pgen',pgen_mode='hardcall',reader_workers=4)[0]
    values=np.empty(fixture['markers'],np.float64); seen=np.zeros(len(values),bool)
    chunks,_=linear_scan_multigpu(source_geno,y,cov,devices=devices,chunk_size=chunk,
        reader_workers=4,prefetch_chunks=2,compute_p_values=False,ordered=False,
        shared_queue_depth=2,reduction_factory=JagwasReduction,_chunk_observer=window)
    count=0
    for start,end,_beta,chi2,_p,_index in chunks:
        if seen[start:end].any(): raise AssertionError('Duplicate variant result')
        values[start:end]=chi2[:,0]; seen[start:end]=True; count+=1
    if not seen.all(): raise AssertionError('Missing variants')
    snapshot=window.snapshot()
    expected=2 if previous['hit'] else 8
    if snapshot['pending'] or len(snapshot['observations'])!=expected*2:
        raise AssertionError('Incomplete expected early-window measurements')
    if any(row['cuda'] is None for row in snapshot['observations']):
        raise AssertionError('Missing device measurements')
    # Validate after the scan, outside any reported component spans.
    with np.load(root/'reference_cells.npz') as cells:
        reference=linear_scan(cells['calls'].astype(np.float64),y,cov,device='cpu',compute_dtype='float64')
        x=np.column_stack([np.ones(len(y)),cov.astype(np.float64)])
        yr=y-x@np.linalg.lstsq(x,y,rcond=None)[0]
        correlation=np.corrcoef(yr.T)
        truth=np.einsum('ij,ij->i',reference[1],np.linalg.solve(correlation,reference[1].T).T)
        error=float(np.max(np.abs(values[cells['indices']]-truth)))
        np.testing.assert_allclose(values[cells['indices']],truth,rtol=3e-4,atol=3e-5)
    np.save(out/'all_joint_statistics.npy',values)
    if execution_context(devices,input_path=path,output_path=out)!=context:
        raise ValueError('Execution context changed during scan')
    if previous['hit']:
        publication=cache.store('stage_observations','initial_chunks_validation',snapshot,
            dependencies=dependencies,provenance=dict(job=str(out),previous_record_sha256=previous['record_sha256']),
            max_age_seconds=300.)
    else:
        publication=window.publish(cache,dependencies=dependencies,provenance=dict(job=str(out)),max_age_seconds=300.)
    write_record(out/'report.json',dict(schema='torchgwas.initial_chunk_audit.v1',fixture=fixture,
        context=context,source_sha256=source,source_cache_hit=structural['hit'],
        source_publication=structural_publication,previous_measurements_hit=previous['hit'],
        previous_measurement_record=previous.get('record_sha256'),publication=publication,
        measurements=snapshot,source_chunks=count,all_variants_retained=True,
        independent_reference_cells=41,max_absolute_chi2_error=error,
        output_sha256=sha256_file(out/'all_joint_statistics.npy'),
        scope='Real CUDA/PGEN execution of a synthetic numerical fixture. First job collects 8 chunks per GPU; cache-hit job retains structural work and takes 2 new observations per GPU. No timing-rank, tuning-readiness, hardware-rate or durable-write-speed claim.'))
    print(json.dumps(dict(source_cache_hit=structural['hit'],previous_measurements_hit=previous['hit'],
        measurements=len(snapshot['observations']),source_chunks=count,variants=len(values),
        max_absolute_chi2_error=error,output_sha256=sha256_file(out/'all_joint_statistics.npy'))))


if __name__=='__main__': main()
