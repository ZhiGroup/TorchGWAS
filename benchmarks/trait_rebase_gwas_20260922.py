"""Persisted public-GWAS parity for old versus single-copy trait rebasing."""
import hashlib
import json
import os
from pathlib import Path
import time
from unittest.mock import patch
import numpy as np
import torch
from torchgwas import api
from torchgwas.detailed_calibration import source_identity,execution_context
from native_host_predicate_gwas_20260922 import output,sha
from legacy_trait_rebase_20260922 import legacy_trait_blocked_significant_chunks,ORIGINAL_FUNCTION_SHA256


def main():
    fixture=Path('/data/zxie3/torchgwas_adaptive_candidate_fixture_v1_20260922')
    root=Path('results/trait_rebase_gwas_20260922');root.mkdir(parents=True,exist_ok=False)
    target=Path('/data/zxie3/torchgwas_trait_rebase_20260922');target.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    source=source_identity();inputs={name:sha(fixture/name) for name in
        ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}
    y=np.load(fixture/'phenotype.npy');covariates=np.load(fixture/'covariates.npy')
    current=api._trait_blocked_significant_chunks;runs=[]
    cases=[('one_gpu_tiled_sparse',.02,'beta+t',['cuda:1']),
           ('two_gpu_tiled_sparse',.02,'t',['cuda:1','cuda:2']),
           ('two_gpu_tiled_dense',1.,'beta+t',['cuda:1','cuda:2']),
           ('two_gpu_tiled_empty',1e-30,'t',['cuda:1','cuda:2'])]
    for i,(case,threshold,fields,devices) in enumerate(cases):
        for backend,selection,predicate in [('host_numpy','host','numpy'),('host_native','host','native'),('device','device','numpy')]:
            os.environ['TORCHGWAS_SIGNIFICANCE_BACKEND']=selection;os.environ['TORCHGWAS_HOST_PREDICATE']=predicate
            context=execution_context(devices,input_path=fixture/'input.pgen',output_path=target)
            reference=None
            for name,function in ([('legacy',legacy_trait_blocked_significant_chunks),('current',current)] if i%2==0
                                  else [('current',current),('legacy',legacy_trait_blocked_significant_chunks)]):
                directory=target/(case+'_'+backend+'_'+name);started=time.perf_counter()
                with patch.object(api,'_trait_blocked_significant_chunks',function):
                    result=api.run_linear_gwas(str(fixture/'input.pgen'),y,covariates,
                        genotype_format='pgen',device='cuda:1',chunk_size=128,reader_workers=3,prefetch_chunks=3,
                        reduce='significant',significance_threshold=threshold,output_dir=directory,
                        sumstats_format='binary',sumstats_fields=fields,sumstats_queue_depth=2,sumstats_fsync=True,
                        trait_block=193,trait_devices=devices)
                elapsed=time.perf_counter()-started;data=output(directory/'sumstats')
                if reference is None:reference=data
                else:
                    assert set(data)==set(reference)
                    for field in data:
                        assert data[field].shape==reference[field].shape and data[field].dtype==reference[field].dtype
                        assert data[field].tobytes()==reference[field].tobytes(),field
                row=dict(case=case,backend=backend,implementation=name,devices=devices,threshold=threshold,
                    fields=fields,output_dir=str(directory),execution_context=context,api_seconds=elapsed,
                    result_rows=result.run_metadata['n_result_rows'],
                    output_hashes={field:hashlib.sha256(memoryview(value)).hexdigest() for field,value in data.items()})
                runs.append(row)
                print(json.dumps({key:row[key] for key in ['case','backend','implementation','result_rows','api_seconds']}),flush=True)
    assert source==source_identity() and inputs=={name:sha(fixture/name) for name in inputs}
    report=dict(runs=runs,all_persisted_fields_exact=True,source_sha256=source,inputs=inputs,
        legacy_function_sha256=ORIGINAL_FUNCTION_SHA256,
        harness_sha256={name:sha(Path(__file__).parent/name) for name in
            ['trait_rebase_gwas_20260922.py','legacy_trait_rebase_20260922.py','native_host_predicate_gwas_20260922.py']},
        scope='Twenty-four public scans establish before/after persisted output parity across host and device selection, tiled one/two-GPU execution, empty/sparse/dense output and t/beta+t. One observation per condition is not a speedup estimate.')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
