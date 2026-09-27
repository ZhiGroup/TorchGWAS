"""Compare independently measured v2 primitives with held-out selectors."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
from torchgwas.detailed_calibration import source_identity
from torchgwas.numpy_nonzero_work import validate_host_price_protocol
from torchgwas.significant_host_work import host_significant_selection_work,host_selection_service


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(args):
    root=Path(args.out);root.mkdir(parents=True,exist_ok=False)
    held_path=Path(args.held_out)
    held=json.loads(held_path.read_text());originals={str(held_path):sha(held_path)}
    source=source_identity();assert held['source_sha256']==source
    results=[];observed_dates={}
    for backend in ['numpy','native']:
        path=Path(getattr(args,backend+'_prices'))
        bank=json.loads(path.read_text());originals[str(path)]=sha(path)
        assert bank['source_sha256']==source
        assert bank['primitive_differences_nonnegative']
        for key in ['row_indices','flat_indices','sparse_flat_indices']:
            assert bank['input_layouts'][key]['strides']==[8]
            assert bank['input_layouts'][key]['dtype']=='<i8'
        assert bank['context']['affinity']==held['affinity']
        assert bank['context']['numpy']==held['numpy']
        for key in ['environment','torch_threads','numpy_madvise_hugepage','numpy_core']:
            assert bank['context'][key]==held['settings'][key],key
        if backend=='native':assert bank['predicate_context']==held['native_context']
        os.environ['TORCHGWAS_HOST_PREDICATE']=backend
        validate_host_price_protocol(bank)
        observed_dates[backend]={key:bank[key] for key in
            ['observation_started_at_utc','observation_finished_at_utc']}
        grouped={}
        for row in held['records']:
            if row['backend']!=backend:continue
            key=(tuple(row['shape']),row['density'],row['retained'],row['operation'])
            grouped.setdefault(key,[]).append(row)
        for (shape,density,retained,operation),observations in grouped.items():
            assert len(observations)==8
            work=host_significant_selection_work(*shape,retained,return_beta=False)
            service=host_selection_service(work,bank['prices'],cpu_fraction=1.,dram_bytes_per_second=1e30,host_serial_fraction=0.)
            selected=service[3:4] if operation=='predicate' else service[1:9]
            prediction=sum(s['seconds']*s['resources']['cpu'] for s in selected)
            median=statistics.median(r['cpu_seconds'] for r in observations)
            row=dict(backend=backend,shape=list(shape),density=density,retained=retained,
                operation=operation,predicted_cpu_seconds=prediction,observed_median_cpu_seconds=median,
                observed_min_cpu_seconds=min(r['cpu_seconds'] for r in observations),
                observed_max_cpu_seconds=max(r['cpu_seconds'] for r in observations),
                predicted_over_observed=prediction/median)
            results.append(row)
            if operation=='selector':print(json.dumps(row),flush=True)
    assert all(sha(path)==digest for path,digest in originals.items())
    (root/'report.json').write_text(json.dumps(dict(results=results,source_sha256=source,
        original_artifacts=originals,original_artifacts_unchanged=True,observation_dates=observed_dates,
        fitted_to_held_out=False,runtime_prediction_validated=False,script_sha256=sha(__file__),
        scope='Independent fixed controls transferred to a different survivor pattern and larger buffers. CPU service only; held-out values do not alter prices or original dates. No whole-GWAS fit or automatic profile installation.'),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='results/selector_gather_price_transfer_v2_20260922')
    parser.add_argument('--held-out',default='results/selector_gather_heldout_v2_20260922/report.json')
    parser.add_argument('--numpy-prices',default='results/selector_gather_numpy_prices_v2_20260922/selection_prices.json')
    parser.add_argument('--native-prices',default='results/selector_gather_native_prices_v2_20260922/selection_prices.json')
    main(parser.parse_args())
