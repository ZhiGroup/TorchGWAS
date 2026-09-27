"""Derive branch-specific prices from original independent observations.

No new timing, observation-date renewal, whole-GWAS fit or profile installation.
The saved measurements remain byte-identical. The new artifact distinguishes
the original measurement source from the current analytical implementation.
"""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import statistics
from torchgwas.detailed_calibration import source_identity
from torchgwas.numpy_nonzero_work import nonzero_protocol, validate_host_price_protocol
from torchgwas.significant_host_work import host_significant_selection_work, host_selection_service


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root=Path('results/nonzero_branch_repricing_20260922');root.mkdir(parents=True,exist_ok=False)
    original_path=Path('results/native_host_predicate_prices_20260922/selection_prices.json')
    held_path=Path('results/native_host_predicate_matched_v2_20260922/report.json')
    original_hash=sha(original_path);held_hash=sha(held_path)
    original=json.loads(original_path.read_text());held=json.loads(held_path.read_text())
    assert original['source_sha256']==held['source_sha256']
    current=source_identity()
    for name in ['host_significance.py','native_host_predicate.py','_host_predicate.cpp']:
        assert current[name]==original['source_sha256'][name], 'Selector implementation changed'
    assert original['context']['numpy']=='2.2.6'
    assert original['predicate_context']==held['native_context']
    assert original['context']['affinity']==held['affinity']
    for key in ['environment','torch_threads','numpy_madvise_hugepage']:
        assert original['context'][key]==held['settings'][key]
    observations=[r for r in original['observations'] if r['primitive']=='flatnonzero_nonempty']
    fixed_rows=[r for r in observations if r['pattern']=='empty']
    assert len(fixed_rows)==9 and all(r['units']==0 for r in fixed_rows)
    fixed=statistics.median(r['cpu_seconds'] for r in fixed_rows)
    derived=deepcopy(original);derived['prices'].pop('flatnonzero_nonempty')
    branches=[]
    for regime,pattern in [('sparse','stride64'),('dense','all_true')]:
        rows=[r for r in observations if r['pattern']==pattern]
        assert len(rows)==9 and len({r['units'] for r in rows})==1
        units=rows[0]['units'];median=statistics.median(r['cpu_seconds'] for r in rows)
        rate=(median-fixed)/units
        assert rate>=0
        derived['prices']['flatnonzero_'+regime]=dict(call_cpu_seconds=fixed,
            unit_cpu_seconds=rate,dram_bytes_per_unit=2)
        branches.append(dict(regime=regime,pattern=pattern,reference_cells=units,
            median_cpu_seconds=median,call_cpu_seconds=fixed,unit_cpu_seconds=rate,
            original_observation_indices=[i for i,r in enumerate(original['observations'])
                if r['primitive']=='flatnonzero_nonempty' and r['pattern'] in ('empty',pattern)]))
    derived.update(nonzero_protocol=nonzero_protocol(),calculator_source_sha256=current,
        derived_from=dict(path=str(original_path),sha256=original_hash),
        derivation='Split existing independent empty/all-true/stride64 observations by the NumPy source branch; original observation dates preserved.',
        prediction_complete=False,concurrent_transfer_qualified=False)
    assert derived['observation_started_at_utc']==original['observation_started_at_utc']
    assert derived['observation_finished_at_utc']==original['observation_finished_at_utc']
    (root/'derived_prices.json').write_text(json.dumps(derived,indent=2)+'\n')
    os.environ['TORCHGWAS_HOST_PREDICATE']='native';validate_host_price_protocol(derived)
    old_report=json.loads(Path('results/native_selector_price_transfer_20260922.json').read_text())
    results=[]
    for old in old_report['results']:
        if old['operation']!='selector':continue
        b,k=old['shape'];retained=old['retained']
        work=host_significant_selection_work(b,k,retained,return_beta=False)
        steps=host_selection_service(work,derived['prices'],cpu_fraction=1.,
            dram_bytes_per_second=1e30,host_serial_fraction=0.)
        # Supplied critical cutoffs and no outer phenotype rebase in this probe.
        predicted=sum(s['seconds']*s['resources']['cpu'] for s in steps[1:9])
        observed=old['observed_median_cpu_seconds']
        row=dict(shape=old['shape'],density=old['density'],retained=retained,
            primitive=work['nonzero_primitive'],original_prediction=old['predicted_cpu_seconds'],
            new_prediction=predicted,observed_median=observed,predicted_over_observed=predicted/observed)
        results.append(row);print(json.dumps(row),flush=True)
    assert sha(original_path)==original_hash and sha(held_path)==held_hash
    (root/'report.json').write_text(json.dumps(dict(branches=branches,results=results,
        original_price_sha256=original_hash,held_out_sha256=held_hash,
        original_observation_started_at_utc=original['observation_started_at_utc'],
        original_observation_finished_at_utc=original['observation_finished_at_utc'],
        original_artifacts_unchanged=True,fitted_to_held_out=False,runtime_prediction_validated=False,
        harness_sha256=sha(__file__),scope=__doc__),indent=2)+'\n')


if __name__=='__main__':main()
