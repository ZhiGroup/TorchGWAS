"""Public output equivalence and measured payload audit; no speedup claim."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np
import torch

from torchgwas.api import run_linear_gwas
import torchgwas.native_scan as native
from torchgwas.pinned_work import pinned_scan_work
from torchgwas.result_service import result_finish_prices
from torchgwas.sumstats import open_binary_sumstats,open_binary_df
from torchgwas.sumstats_indexed import open_indexed_sumstats


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_output(path,significant):
    if significant:
        manifest,parts=open_indexed_sumstats(path)
        parts=list(parts)
        assert parts and manifest['rows']>0
        assert all('beta' not in part for part in parts)
        data={name:np.concatenate([p[name] for p in parts]) for name in ('variant_index','trait_index','t_stat','df')}
        order=np.lexsort((data['trait_index'],data['variant_index']))
        assert len(np.unique(np.column_stack((data['variant_index'],data['trait_index'])),axis=0))==len(order)
        return {name:array[order] for name,array in data.items()}
    beta,t,_=open_binary_sumstats(path)
    assert beta is None and not list(path.rglob('beta.f32'))
    return dict(t_stat=np.asarray(t),df=np.asarray(open_binary_df(path)))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fixture',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--output-data',required=True);parser.add_argument('--finish-probe',required=True)
    args=parser.parse_args();fixture=Path(args.fixture);out=Path(args.out);target=Path(args.output_data)
    out.mkdir(parents=True,exist_ok=False);target.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    y=np.load(fixture/'phenotype.npy');c=np.load(fixture/'covariates.npy')
    original=native.dosage_cuda_iterator
    rows=[];equivalence=[]
    scenarios=[('dense_single',{}),('dense_variants',dict(variant_devices=['cuda:1','cuda:2'])),
        ('dense_traits',dict(trait_block=193,trait_devices=['cuda:1','cuda:2'])),
        ('significant_traits',dict(trait_block=193,trait_devices=['cuda:1','cuda:2'],
                                  reduce='significant',significance_threshold=.02))]
    for index,(name,options) in enumerate(scenarios):
        reference=None
        for omit in ([False,True] if index%2==0 else [True,False]):
            traces=[];first=[];started=time.perf_counter()
            def traced(source,phenotype,basis,chunk_size,device,reader_workers=None,
                       prefetch_chunks=None,compute_p_values=True,**kw):
                assert kw.get('return_beta') is False,'Public t-only dispatch did not omit beta'
                kw['return_beta']=not omit
                chunks=original(source,phenotype,basis,chunk_size,device,reader_workers,
                                prefetch_chunks,compute_p_values,**kw)
                try:
                    for item in chunks:
                        if not first:first.append(time.perf_counter()-started)
                        assert (item[2] is None)==omit
                        yield item
                finally:
                    chunks.close()
                span=kw.get('variant_range') or (0,source.shape[1]);m=span[1]-span[0];k=phenotype.shape[1]
                profile=dict(source._last_scan_profile)
                expected=m*(4*k*(1+int(not omit))+5)
                assert profile['result_payload_bytes']==expected
                traces.append(dict(device=str(device),variants=m,traits=k,variant_range=list(span),
                    payload_bytes=expected,profile=profile,
                    pinned_work=pinned_scan_work(source.shape[0],chunk_size,k,prefetch_chunks,return_beta=not omit)))
            directory=target/(name+('_omit' if omit else '_baseline'))
            with patch.object(native,'dosage_cuda_iterator',traced):
                result=run_linear_gwas(str(fixture/'input.pgen'),y,c,genotype_format='pgen',
                    device='cuda:1',chunk_size=128,reader_workers=3,prefetch_chunks=3,
                    output_dir=directory,sumstats_format='binary',sumstats_fields='t',
                    sumstats_block_bytes=1<<20,sumstats_queue_depth=2,sumstats_fsync=True,**options)
            elapsed=time.perf_counter()-started
            assert traces
            values=read_output(directory/'sumstats',name.startswith('significant'))
            if reference is None:reference={field:array.copy() for field,array in values.items()}
            else:
                assert set(reference)==set(values)
                for field in values:np.testing.assert_array_equal(values[field],reference[field])
                equivalence.append(dict(scenario=name,fields=list(values),identical=True))
            rows.append(dict(scenario=name,omit_beta=omit,api_seconds=elapsed,
                first_scientific_chunk_seconds=first[0],traces=traces,
                payload_bytes=sum(t['payload_bytes'] for t in traces),
                pinned_requested_bytes=sum(t['pinned_work']['requested_bytes'] for t in traces),
                output_dir=str(directory),sumstats_summary=result.run_metadata['sumstats_write']))
            print('COMPLETE',name,omit,elapsed,flush=True)
    for name,_ in scenarios:
        old,new=[next(row for row in rows if row['scenario']==name and row['omit_beta']==omit) for omit in [False,True]]
        expected=4*sum(t['variants']*t['traits'] for t in old['traces'])
        assert old['payload_bytes']-new['payload_bytes']==expected
    probe=json.loads(Path(args.finish_probe).read_text());finish={}
    for workers in (1,4):
        try:
            finish[str(workers)]=dict(status='accepted',service=result_finish_prices(probe,workers=workers,
                numpy_version=np.__version__,torch_version=torch.__version__,python_version=__import__('sys').version,
                cpu_affinity=list(range(12,20)),source_sha256=digest('src/torchgwas/native_scan.py'),return_beta=False))
        except ValueError as error:finish[str(workers)]=dict(status='rejected',reason=str(error))
    report=dict(runs=rows,equivalence=equivalence,finish_validation=finish,
        source_sha256={str(p):digest(p) for p in sorted(Path('src/torchgwas').glob('*.py'))},
        inputs={str(fixture/n):digest(fixture/n) for n in ['input.pgen','input.pvar','input.psam','phenotype.npy','covariates.npy']},
        scope='Synthetic PGEN fixture, N=2049 M=4097 K=512; real public CUDA execution and durable output on local /data. Baseline forces beta transfer while keeping t-only writer. One pair per mode with alternating order; no controlled performance inference. First chunk means computed output available, not durable file publication. Statistics kernels are unchanged.')
    (out/'report.json').write_text(json.dumps(report,indent=2))


if __name__=='__main__':main()
