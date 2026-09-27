"""Audit saved independent primitive windows without renewing measurements."""
from collections import Counter
import hashlib
import json
from pathlib import Path

from torchgwas.cpu_service_refresh import CpuServiceRefresh
from torchgwas.detailed_calibration import source_identity


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root=Path('results/cpu_window_stability_audit_20260922')
    root.mkdir(parents=True,exist_ok=False)
    paths=[Path('results')/name/'selection_prices.json' for name in [
        'selector_gather_numpy_prices_20260922','selector_gather_native_prices_20260922',
        'selector_gather_numpy_prices_v2_20260922','selector_gather_native_prices_v2_20260922']]
    before={str(path):sha(path) for path in paths}
    inspector=CpuServiceRefresh(root/'unused-cache','historical-window-audit',
        dependencies=dict(measurement_protocol={'operation':'read-only historical inspection'}),
        work_units=1,max_age_seconds=3600.)
    windows=[];original_dates={}
    for path in paths:
        artifact=json.loads(path.read_text())
        original_dates[str(path)]={key:artifact[key] for key in
            ['observation_started_at_utc','observation_finished_at_utc']}
        groups={}
        for row in artifact['observations']:
            groups.setdefault((row['primitive'],row['control']),[]).append(row)
        for (primitive,control),rows in groups.items():
            rows=sorted(rows,key=lambda r:r['repeat'])
            assert [r['repeat'] for r in rows]==list(range(9))
            # Fixed selection rule, declared before looking at outcomes: the
            # first seven repeat means, matching the controller window length.
            selected=rows[:7]
            assert len({r['units'] for r in selected})==1
            spread=inspector._spread(selected)
            stability=inspector._window_stability(selected)
            windows.append(dict(artifact=str(path),primitive=primitive,control=control,
                units=selected[0]['units'],pattern=selected[0]['pattern'],repeats=list(range(7)),
                cpu_seconds=[r['cpu_seconds'] for r in selected],
                previous_spread_check_passed=spread['stable'],window_stability=stability,
                newly_rejected=spread['stable'] and not stability['stable']))
    assert before=={str(path):sha(path) for path in paths}
    assert not (root/'unused-cache').exists()
    summary=dict(windows=len(windows),previously_accepted=sum(w['previous_spread_check_passed'] for w in windows),
        accepted_after_temporal_check=sum(w['window_stability']['stable'] for w in windows),
        newly_rejected=sum(w['newly_rejected'] for w in windows),
        newly_rejected_primitives=dict(Counter(w['primitive'] for w in windows if w['newly_rejected'])))
    report=dict(summary=summary,windows=windows,original_artifacts=before,
        original_dates=original_dates,original_artifacts_unchanged=True,measurements_published=0,
        source_sha256=source_identity(),script_sha256=sha(__file__),policy=inspector.snapshot()['policy'],
        scope='Historical audit of the first seven per-call repeat means from independent primitive controls. No new timings, no fitted coefficients, no measurement-age renewal and no statistical confidence guarantee. Earlier non-v2 gather controls retain their known index-layout limitation.')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    for window in windows:
        if window['newly_rejected']:
            print(json.dumps({key:window[key] for key in ['artifact','primitive','control','cpu_seconds','window_stability']}),flush=True)


if __name__=='__main__':main()
