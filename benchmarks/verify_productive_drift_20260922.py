"""Verify immutable records, source consumption and completed real audit jobs."""
import json
from pathlib import Path

from torchgwas.calibration_cache import _digest
from torchgwas.detailed_calibration import sha256_file, source_identity


def main():
    sources = source_identity(); reports = []; rows = []
    for name in ('first','reuse','expiry','followup'):
        path = Path(f'results/productive_drift_{name}_v1_20260922/report.json')
        report = json.loads(path.read_text()); reports.append(report)
        assert report['source_sha256']==sources
        assert report['script_sha256']==sha256_file('benchmarks/direct_productive_forecast_20260922.py')
        assert all(sha256_file(path)==digest for path,digest in report['helper_sha256'].items())
        assert all(sha256_file(path)==digest for path,digest in report['inputs'].items())
        evidence = report['component_evidence']; result = evidence['result']
        assert evidence['state']=='ready'
        record = result['record']; artifact = Path(result['path'])
        assert json.loads(artifact.read_text())==record
        assert _digest(record)==artifact.stem==result['record_sha256']
        binding = evidence['evidence']['bindings'][0]
        assert sha256_file(artifact)==binding['artifact_sha256']
        assert record['observed_unix_seconds']==min(sample['observed_unix_seconds'] for sample in record['value']['samples'])
        assert len(record['value']['samples'])==7
        for execution in report['records']:
            assert execution['rows']==4097 and execution['parts']==33
            assert execution['max_absolute_difference']==0.
        execution = report['records'][1]; steps = execution['steps']; measurements = execution['measurement_steps']
        counts = [len(step['state']['samples']) for step in measurements]
        assert counts==([2] if result['status']=='reused_original' else [2,4,6,7])
        assert [step['written_parts'] for step in measurements]==list(range(1,len(measurements)+1))
        assert len(steps)==len(measurements)+1
        assert all(step['evaluated'] and not step['applied'] and not step.get('error') for step in steps)
        for previous,current in zip(measurements,measurements[1:]):
            before=previous['state']['samples'];assert current['state']['samples'][:len(before)]==before
        audit = steps[-1]['forecast_audit']; price = audit['price_evidence']
        assert price['verified_targets']==2 and audit['forecast_status']=='unstable_marginal_cost'
        assert price['bindings'][0]['record_sha256']==result['record_sha256']
        assert price['bindings'][0]['age_seconds']>=binding['age_seconds']
        used = evidence['consumed_copy_prices']; assert len(used)==12
        for service in used:
            assert service['copy_cpu_seconds']>0
            assert service['measured_cpu_seconds_per_byte']==record['value']['cpu_seconds_per_unit']
            assert service['copy_cpu_seconds']==service['additional_copy_bytes']*record['value']['cpu_seconds_per_unit']
        assert abs(steps[-1]['total_tuning_cost_seconds']-sum(step['wall_seconds'] for step in steps)-.01)<1e-10
        rows.append(dict(job=name,record_sha256=artifact.stem,artifact_sha256=sha256_file(artifact),
            status=result['status'],lookup_reason=evidence['lookup_reason'],check=evidence['check'],
            samples_per_callback=[counts[0]]+[b-a for a,b in zip(counts,counts[1:])],
            measurement_wall_seconds=sum(step['wall_seconds'] for step in steps[:-1]),
            forecast_wall_seconds=steps[-1]['wall_seconds'],total_tuning_cost_seconds=steps[-1]['total_tuning_cost_seconds'],
            binding_age_seconds=binding['age_seconds'],decision_age_seconds=price['bindings'][0]['age_seconds'],
            cpu_seconds_per_byte=evidence['cpu_seconds_per_byte'],source_consumption_checks=len(used),
            preparation_seconds=report['admission_seconds'],
            output_inclusive_seconds=[row['output_inclusive_seconds'] for row in report['records']],
            first_written_seconds=[row['first_output_seconds'] for row in report['records']]))
    assert rows[0]['status']=='measured_and_published'
    assert rows[1]['status']=='reused_original' and rows[1]['check']['status']=='consistent'
    assert reports[0]['component_evidence']['result']['record']==reports[1]['component_evidence']['result']['record']
    assert rows[0]['record_sha256']==rows[1]['record_sha256']
    assert rows[0]['artifact_sha256']==rows[1]['artifact_sha256']
    assert rows[2]['status']=='measured_and_published' and rows[2]['lookup_reason']=='expired'
    assert rows[2]['record_sha256']!=rows[0]['record_sha256']
    result=dict(source_files_verified=len(sources),all_package_script_helper_input_hashes_verified=True,
        original_artifact_unchanged=True,real_executions_verified=8,rows=rows,
        scope='Evidence lifecycle and model wiring, not prediction accuracy or throughput benefit. Controls run first; warm states differ.')
    Path('results/productive_drift_verification_v1_20260922.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
