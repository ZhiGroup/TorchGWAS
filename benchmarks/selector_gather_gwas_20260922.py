"""Persisted output parity for gather v1/v2 and both host predicates."""
import hashlib
import json
import os
from pathlib import Path
import time
from unittest.mock import patch
import numpy as np
import torch
from torchgwas.api import run_linear_gwas
from torchgwas.detailed_calibration import source_identity,execution_context
from torchgwas import host_significance
from native_host_predicate_gwas_20260922 import output,sha
from selector_gather_control_20260922 import legacy_select,exact


def main():
    fixture=Path('/data/zxie3/torchgwas_adaptive_candidate_fixture_v1_20260922')
    root=Path('results/selector_gather_gwas_20260922')
    target=Path('/data/zxie3/torchgwas_selector_gather_20260922')
    root.mkdir(parents=True,exist_ok=False);target.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    source=source_identity();inputs={n:sha(fixture/n) for n in
        ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']}
    y=np.load(fixture/'phenotype.npy');covariates=np.load(fixture/'covariates.npy')
    cases=[('single_sparse',.02,'beta+t',{}),
        ('tiled_sparse',.02,'t',dict(trait_block=193,trait_devices=['cuda:1','cuda:2'])),
        ('tiled_dense',1.,'beta+t',dict(trait_block=193,trait_devices=['cuda:1','cuda:2'])),
        ('tiled_empty',1e-30,'t',dict(trait_block=193,trait_devices=['cuda:1','cuda:2']))]
    current=host_significance.select_host_pairs;reports=[]
    for i,(label,threshold,fields,options) in enumerate(cases):
        reference=None
        for backend in ['numpy','native']:
            os.environ['TORCHGWAS_HOST_PREDICATE']=backend
            context=execution_context(['cuda:1','cuda:2'],input_path=fixture/'input.pgen',output_path=target)
            for implementation in (['legacy','current'] if i%2==0 else ['current','legacy']):
                directory=target/(label+'_'+backend+'_'+implementation)
                function=current if implementation=='current' else legacy_select
                started=time.perf_counter()
                with patch.object(host_significance,'select_host_pairs',function):
                    result=run_linear_gwas(str(fixture/'input.pgen'),y,covariates,
                        genotype_format='pgen',device='cuda:1',chunk_size=128,
                        reader_workers=3,prefetch_chunks=3,reduce='significant',
                        significance_threshold=threshold,output_dir=directory,
                        sumstats_format='binary',sumstats_fields=fields,
                        sumstats_queue_depth=2,sumstats_fsync=True,**options)
                elapsed=time.perf_counter()-started;data=output(directory/'sumstats')
                if reference is None:reference=data
                else:
                    assert set(data)==set(reference)
                    exact([data[key] for key in sorted(data)],[reference[key] for key in sorted(reference)])
                row=dict(case=label,backend=backend,implementation=implementation,
                    threshold=threshold,fields=fields,output_dir=str(directory),
                    context=context,api_seconds=elapsed,result_rows=result.run_metadata['n_result_rows'],
                    output_hashes={k:hashlib.sha256(memoryview(v).cast('B')).hexdigest() for k,v in data.items()})
                reports.append(row)
                print(json.dumps({k:row[k] for k in ['case','backend','implementation','result_rows','api_seconds']}),flush=True)
    assert source_identity()==source
    assert inputs=={n:sha(fixture/n) for n in inputs}
    report=dict(runs=reports,all_persisted_fields_exact=True,source_sha256=source,inputs=inputs,
        harness_sha256={name:sha(Path(__file__).parent/name) for name in
            ['selector_gather_gwas_20260922.py','selector_gather_control_20260922.py','native_host_predicate_gwas_20260922.py']},
        scope='Sixteen public scans compare identical statistics and durable output fields. Single observations per configuration do not establish GWAS speedup.')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
