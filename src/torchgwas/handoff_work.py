"""Conditional elapsed handoff service from independent Python wait probes."""
import math
import statistics


def handoff_wakeup_prices(probe, *, pairs, python_version, cpu_affinity):
    if probe['python_version']!=python_version or probe['affinity']!=cpu_affinity:
        raise ValueError('Handoff runtime or CPU context differs from the probe')
    if isinstance(pairs,bool) or pairs not in (1,2):raise ValueError('Unsupported measured handoff context')
    result={}
    for kind in ('queue','future'):
        rows=[row for row in probe['rows'] if row['kind']==kind and row['pairs']==pairs]
        if len(rows)<3:raise ValueError('Insufficient handoff repeats')
        values=[]
        for row in rows:
            observations=row['observations']
            if not observations or any(item.get('blocked_verified') is not True for item in observations):
                raise ValueError('Blocked handoff context was not verified')
            latency=[item['publication_to_return_seconds'] for item in observations]
            ready=row['ready']
            numbers=latency+[ready['producer_cpu_seconds'],ready['receiver_cpu_seconds']]
            if any(isinstance(value,bool) or not math.isfinite(value) or value<0 for value in numbers):
                raise ValueError('Invalid independent handoff service')
            # Queue put/get CPU are already in the source control ledger.
            # Future publication CPU is not, so only remove ready-result CPU.
            baseline=ready['receiver_cpu_seconds']+(ready['producer_cpu_seconds'] if kind=='queue' else 0.)
            values.append(max(0.,statistics.fmean(latency)-baseline))
        result[kind]=statistics.median(values)
    return dict(wakeup_seconds=result,pairs=pairs,
                scope='Median of per-repeat mean blocked handoffs less independently measured ready work. Applied only when source dependency timestamps imply blocking. Extra wake/park CPU is not separately represented; live-context changes remain unpriced. No scan observations used.')
