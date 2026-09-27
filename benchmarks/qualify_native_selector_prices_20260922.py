"""Held-out component-price transfer audit; never fit or republish a price."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import statistics


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    price_path=Path('results/native_host_predicate_prices_20260922/selection_prices.json')
    observation_path=Path('results/native_host_predicate_matched_v2_20260922/report.json')
    bank=json.loads(price_path.read_text()); held=json.loads(observation_path.read_text())
    assert bank['host_selector']=='native_row_flat_v1'
    assert bank['source_sha256']==held['source_sha256']
    assert bank['predicate_context']==held['native_context']
    assert bank['context']['affinity']==held['affinity']
    assert bank['context']['environment']==held['settings']['environment']
    assert bank['context']['torch_threads']==held['settings']['torch_threads']
    assert bank['context']['numpy_madvise_hugepage']==held['settings']['numpy_madvise_hugepage']
    groups=defaultdict(list)
    for row in held['records']:
        if row['backend']=='native': groups[tuple(row['shape']),row['density'],row['operation']].append(row)
    results=[]
    for (shape,density,operation),rows in sorted(groups.items()):
        b,k=shape; cells=b*k; retained=rows[0]['retained']
        assert all(row['retained']==retained for row in rows)
        terms=[('predicate_native',1,cells)]
        if operation=='selector':
            # Exactly select_host_pairs: cutoffs are supplied, beta is omitted,
            # and the outer phenotype-tile coordinate adjustment is not timed.
            terms=[('critical_round',1,b),('mask_allocate',1,0),*terms,
                ('flatnonzero_nonempty' if retained else 'flatnonzero_empty',1,cells),
                ('coordinate_divmod',1,retained),('df_gather',1,retained),
                ('matrix_gather',1,retained),('inplace_index_add',1,retained)]
        priced=[dict(primitive=name,calls=calls,units=units,
            cpu_seconds=calls*bank['prices'][name]['call_cpu_seconds']+units*bank['prices'][name]['unit_cpu_seconds'])
            for name,calls,units in terms]
        predicted=sum(row['cpu_seconds'] for row in priced)
        observed=statistics.median(row['cpu_seconds'] for row in rows)
        results.append(dict(shape=shape,density=density,operation=operation,retained=retained,
            predicted_cpu_seconds=predicted,observed_median_cpu_seconds=observed,
            predicted_over_observed=predicted/observed,
            observed_min_cpu_seconds=min(row['cpu_seconds'] for row in rows),
            observed_max_cpu_seconds=max(row['cpu_seconds'] for row in rows),
            min_minor_faults=min(row['minor_faults'] for row in rows),
            max_minor_faults=max(row['minor_faults'] for row in rows),terms=priced))
    report=dict(results=results,price_artifact=dict(path=str(price_path),sha256=sha(price_path),
        original_observation_started_at_utc=bank['observation_started_at_utc'],
        original_observation_finished_at_utc=bank['observation_finished_at_utc']),
        held_out_artifact=dict(path=str(observation_path),sha256=sha(observation_path)),
        harness_sha256=sha(__file__),source_sha256=bank['source_sha256'],
        fitted_parameters=False,price_records_modified=False,runtime_prediction_validated=False,
        scope=__doc__)
    with Path('results/native_selector_price_transfer_20260922.json').open('x') as stream:
        json.dump(report,stream,indent=2);stream.write('\n')
    for row in results:
        print(json.dumps({k:v for k,v in row.items() if k!='terms'}),flush=True)


if __name__=='__main__': main()
