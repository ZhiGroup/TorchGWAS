"""Transfer prices qualify only on held-out agreement and bind exactly."""
from copy import deepcopy

import pytest

from torchgwas.calibration_cache import CalibrationParameterCache
from torchgwas.detailed_calibration import bind_detailed_profile, sha256_file
from torchgwas.layout_transfer_links import transfer_link_loads
from torchgwas.price_binding import validate_price_bindings
from torchgwas.transfer_calibration import (KIND, NAME, apply_transfer_prices, declared_groups,
    summarize_transfer_groups, transfer_price_targets)


RATES = {'cuda:0': 10e9, 'cuda:1': 12e9, 'cuda:0,cuda:1': 15e9}
PAIR = [['cuda:0'], ['cuda:1'], ['cuda:0', 'cuda:1']]
TRIO_RATES = {'cuda:0': 10e9, 'cuda:1': 12e9, 'cuda:2': 11e9, 'cuda:0,cuda:1': 15e9,
              'cuda:0,cuda:1,cuda:2': 20e9}
TRIO = [['cuda:0'], ['cuda:1'], ['cuda:2'], ['cuda:0', 'cuda:1'], ['cuda:0', 'cuda:1', 'cuda:2']]


def raw(rates=RATES, rounds=6, jitter=0.02, busy_round=None, nodes=None, groups=PAIR):
    devices = [g[0] for g in groups if len(g) == 1]
    records = []
    for index in range(rounds):
        scale = 1 + jitter*((index % 3)-1)
        observations = [dict(direction=d, devices=g,
                             bytes_per_second=rates[','.join(g)]*scale*(0.9 if d == 'd2h' else 1))
                        for d in ('h2d', 'd2h') for g in groups]
        records.append(dict(round=index, half='calibration' if index % 2 == 0 else 'holdout',
                            started_unix_seconds=1000.+index, observations=observations,
                            busy_selected_gpus=['uuid0'] if index == busy_round else []))
    placement = {d: dict(sampled_pages=8, total_pages=8, status_counts=nodes or {'0': 8})
                 for d in devices}
    return dict(devices=devices, groups=deepcopy(groups), uuids={d: 'u'+d[5:] for d in devices},
                placement_before=placement, placement_after=deepcopy(placement), records=records)


def test_groups_are_singles_links_then_whole_set():
    assert declared_groups(['cuda:0', 'cuda:1', 'cuda:2'], [('cuda:0', 'cuda:1')]) == [
        ('cuda:0',), ('cuda:1',), ('cuda:2',), ('cuda:0', 'cuda:1'), ('cuda:0', 'cuda:1', 'cuda:2')]
    assert declared_groups(['cuda:3']) == [('cuda:3',)]
    for links in ([('cuda:0',)], [('cuda:0', 'cuda:9')], [('cuda:0', 'cuda:1'), ('cuda:1', 'cuda:0')]):
        with pytest.raises(ValueError):
            declared_groups(['cuda:0', 'cuda:1'], links)


def test_stable_rounds_qualify_with_ordered_scenarios():
    summary = summarize_transfer_groups(raw())
    assert summary['qualified'] and not summary['failures']
    scenarios = summary['value']['scenarios']
    assert (scenarios['low']['shared_transfer_capacities']['h2d'] <
            scenarios['median']['shared_transfer_capacities']['h2d'] <
            scenarios['high']['shared_transfer_capacities']['h2d'])
    assert scenarios['median']['per_device']['cuda:1']['h2d'] == pytest.approx(12e9)
    assert scenarios['median']['shared_links'] == []


@pytest.mark.parametrize('change,check', [
    (dict(busy_round=3), 'selected_gpu_idle'),
    (dict(nodes={'0': 4, '1': 4}), 'pinned_single_numa_node'),
    (dict(nodes={'-2': 8}), 'pinned_single_numa_node'),
    (dict(rates={'cuda:0': 10e9, 'cuda:1': 12e9, 'cuda:0,cuda:1': 30e9}), 'group_not_superadditive'),
])
def test_violations_block_publication(change, check):
    summary = summarize_transfer_groups(raw(**change))
    assert not summary['qualified'] and 'value' not in summary
    assert check in {row['check'] for row in summary['failures']}


def test_holdout_drift_blocks_publication():
    data = raw()
    for record in data['records']:
        if record['half'] == 'holdout':
            for row in record['observations']:
                row['bytes_per_second'] *= 0.7
    summary = summarize_transfer_groups(data)
    assert not summary['qualified']
    assert {row['check'] for row in summary['failures']} == {'holdout'}


def test_tight_span_tolerates_small_holdout_offset_but_not_busy_start():
    data = raw(jitter=0.001)
    for record in data['records']:
        if record['half'] == 'holdout':
            for row in record['observations']:
                row['bytes_per_second'] *= 0.995
    assert summarize_transfer_groups(data)['qualified']
    data['initial_busy_selected_gpus'] = ['u0']
    summary = summarize_transfer_groups(data)
    assert [row['check'] for row in summary['failures']] == ['selected_gpu_idle_at_start']


def test_expected_node_is_enforced():
    summary = summarize_transfer_groups(raw(), expected_nodes={'cuda:0': 0, 'cuda:1': 1})
    assert [row['device'] for row in summary['failures']] == ['cuda:1', 'cuda:1']


def test_applied_scenario_binds_to_immutable_record(tmp_path, monkeypatch):
    now = [2000.]
    monkeypatch.setattr('time.time', lambda: now[0])
    value = summarize_transfer_groups(raw(TRIO_RATES, groups=TRIO))['value']
    assert [link['devices'] for link in value['scenarios']['high']['shared_links']] == [['cuda:0', 'cuda:1']]
    source = {'m.py': 'x'}
    execution = dict(devices={d: dict(uuid='u'+d[5:]) for d in value['devices']}, affinity=[0])
    deps = dict(source_sha256=source, execution_context=execution,
                measurement_protocol=dict(operation='pinned'))
    record = CalibrationParameterCache(tmp_path).store(KIND, NAME, value, dependencies=deps,
        provenance=dict(job='t'), max_age_seconds=60., observed_unix_seconds=1990.)
    context = dict(name='trio', devices=['cuda:0', 'cuda:1', 'cuda:2'],
                   profiles={'cuda:0': dict(h2d_bytes_per_second=1.), 'cuda:1': dict(), 'cuda:2': dict()})
    with pytest.raises(ValueError, match='exactly'):
        apply_transfer_prices(dict(context, devices=['cuda:0']), value, 'high')
    applied = apply_transfer_prices(context, value, 'high')
    assert context['profiles']['cuda:0']['h2d_bytes_per_second'] == 1.
    assert applied['profiles']['cuda:0']['h2d_bytes_per_second'] == value['scenarios']['high']['per_device']['cuda:0']['h2d']
    binding = dict(artifact=record['path'], kind=KIND, name=NAME, dependencies=deps,
                   max_age_seconds=None, targets=transfer_price_targets(0, value, 'high'))
    profile = bind_detailed_profile([applied], execution, sources=source, limitations=['test'],
        component_artifacts={record['path']: sha256_file(record['path'])}, price_bindings=[binding])
    evidence = validate_price_bindings(profile)
    assert evidence['status'] == 'declared_targets_verified' and evidence['verified_targets'] == 10
    loads = transfer_link_loads({'cuda:0': 10, 'cuda:1': 10, 'cuda:2': 10}, applied['shared_links'],
                                direction='h2d')
    assert loads['links'][0]['bytes'] == 20
    tampered = deepcopy(profile)
    tampered['contexts'][0]['shared_transfer_capacities']['h2d'] *= 2
    with pytest.raises(ValueError):
        validate_price_bindings(tampered)
    now[0] = 2051.
    with pytest.raises(ValueError, match='[Ee]xpired'):
        validate_price_bindings(profile)


def test_attach_keeps_existing_bindings_and_refuses_rebinding(tmp_path, monkeypatch):
    from torchgwas.transfer_calibration import attach_transfer_prices
    monkeypatch.setattr('time.time', lambda: 2000.)
    value = summarize_transfer_groups(raw(TRIO_RATES, groups=TRIO))['value']
    source = {'m.py': 'x'}
    execution = dict(devices={d: dict(uuid='u'+d[5:]) for d in value['devices']}, affinity=[0])
    deps = dict(source_sha256=source, execution_context=execution,
                measurement_protocol=dict(operation='independent'))
    cache = CalibrationParameterCache(tmp_path)
    cpu = cache.store('cpu_capacity', 'copy', dict(rate=3.), dependencies=deps,
                      provenance=dict(job='cpu'), max_age_seconds=60., observed_unix_seconds=1990.)
    transfer = cache.store(KIND, NAME, value, dependencies=deps, provenance=dict(job='t'),
                           max_age_seconds=60., observed_unix_seconds=1995.)
    contexts = [dict(name='trio', devices=value['devices'],
                     profiles={d: dict(rate=3.) for d in value['devices']})]
    cpu_binding = dict(artifact=cpu['path'], kind='cpu_capacity', name='copy', dependencies=deps,
        max_age_seconds=None, targets=[dict(context_path=[0, 'profiles', 'cuda:0', 'rate'],
                                            value_path=['rate'])])
    profile = bind_detailed_profile(contexts, execution, sources=source, limitations=['test'],
        component_artifacts={cpu['path']: sha256_file(cpu['path'])}, price_bindings=[cpu_binding])
    attached = attach_transfer_prices(profile, transfer['path'], context_name='trio', scenario='low')
    assert len(attached['price_bindings']) == 2
    assert validate_price_bindings(attached)['verified_targets'] == 11
    assert attached['contexts'][0]['shared_transfer_capacities'] == value['scenarios']['low']['shared_transfer_capacities']
    with pytest.raises(ValueError):
        attach_transfer_prices(attached, transfer['path'], context_name='trio', scenario='low')
    with pytest.raises(ValueError, match='context'):
        attach_transfer_prices(dict(profile, execution_context=dict(execution, affinity=[1])),
                               transfer['path'], context_name='trio', scenario='low')
    with pytest.raises(ValueError, match='named'):
        attach_transfer_prices(profile, transfer['path'], context_name='other', scenario='low')


def test_link_scenario_targets_cover_every_link():
    value = dict(devices=['cuda:0'], scenarios={s: dict(per_device={'cuda:0': dict(h2d=1., d2h=1.)},
        shared_transfer_capacities=dict(h2d=1., d2h=1.),
        shared_links=[dict(devices=['cuda:0'], h2d_bytes_per_second=1., d2h_bytes_per_second=1.)])
        for s in ('low', 'median', 'high')})
    paths = [t['context_path'] for t in transfer_price_targets(3, value, 'low')]
    assert [3, 'shared_links', 0, 'd2h_bytes_per_second'] in paths and len(paths) == 6
    with pytest.raises(ValueError):
        transfer_price_targets(0, value, 'best')


def test_real_gpu_measurement_records_placement_and_rounds():
    torch = pytest.importorskip('torch')
    if not torch.cuda.is_available():
        pytest.skip('CUDA required')
    from torchgwas.transfer_calibration import measure_transfer_groups, pinned_page_nodes
    host = torch.empty(1 << 20, dtype=torch.uint8, pin_memory=True)
    host.fill_(1)
    nodes = pinned_page_nodes(host.data_ptr(), host.numel())
    assert nodes['sampled_pages'] > 0 and all(int(k) >= 0 for k in nodes['status_counts'])
    data = measure_transfer_groups(['cuda:0'], size_bytes=1 << 20, copies=1, samples=1,
                                   warmups=0, rounds=4, idle_utilization_limit=100)
    assert [r['half'] for r in data['records']] == ['calibration', 'holdout']*2
    assert all(o['bytes_per_second'] > 0 for r in data['records'] for o in r['observations'])
