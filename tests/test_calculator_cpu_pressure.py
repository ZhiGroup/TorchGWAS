"""CPU telemetry invariants, independent of CUDA and system load."""
import importlib.util
from pathlib import Path
import pytest
spec=importlib.util.spec_from_file_location('cpu_pressure',Path(__file__).resolve().parents[1]/'benchmarks/direct_calculator_cpu_pressure.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def test_guest_ticks_are_not_double_counted():
    row=module.parse_affinity_cpu_stat('cpu1 100 20 30 40 5 2 3 4 70 8\n',[1])
    assert sum(row['cpu1'].values())==204
    assert 'guest' not in row['cpu1']

@pytest.mark.parametrize('text',['cpu1 1 2 3','cpu2 1 2 3 4 5 6 7 8','cpu1 -1 2 3 4 5 6 7 8'])
def test_incomplete_or_invalid_observations_rejected(text):
    with pytest.raises(ValueError):module.parse_affinity_cpu_stat(text,[1])

def test_interval_separates_idle_iowait_and_busy():
    ticks=module.parse_affinity_cpu_stat('cpu1 0 0 0 0 0 0 0 0',[1])
    before=dict(affinity=[1],ticks_per_second=100,monotonic_seconds=1,cpu_ticks=ticks)
    after=dict(before,monotonic_seconds=2,cpu_ticks=module.parse_affinity_cpu_stat('cpu1 20 0 10 60 5 1 2 2',[1]))
    result=module.affinity_cpu_delta(before,after)
    assert result['total_cpu_seconds']==pytest.approx(1)
    assert result['busy_cpu_seconds']==pytest.approx(.35)
    assert module.affinity_cpu_delta(before,dict(after,affinity=[2]))['status']=='incomparable'
    assert module.affinity_cpu_delta(after,dict(before,monotonic_seconds=3))['status']=='incomparable'
