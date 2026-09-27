"""No staged work price is qualified by a matching profile name alone."""
from copy import deepcopy

import pytest

from torchgwas.productive_staged_work_binding import audit_staged_work_price_binding


def case():
    writer=dict(cpu_fraction=.5,executor_cpu_seconds=.01,
        process_units=dict(bytearray_zero_bytes=.001,numpy_copy_bytes=.002),
        writeback_service=dict(pagecache_seconds_per_byte=.001,
            storage_seconds_per_byte=.002,submit_seconds=.003,
            wait_seconds=.004,fadvise_seconds=.005),fsync_seconds=.006,
        writer_copy_service=dict(cpu_seconds_per_byte=.007,
                                 cpu_seconds_per_call=.008))
    active=dict(writer,h2d_bytes_per_second=800.,d2h_bytes_per_second=400.,
        shared_dram_bytes_per_second=1e9,
        gpu_resources=dict(fp32_flops_per_second=100000.,
                           fp64_flops_per_second=10000.,
                           kernel_launch_seconds=.0001,gpu_fraction=1.))
    profile=dict(contexts=[dict(name='current',profiles={'cuda:0':active},
        shared_capacities=dict(output=200.),
        shared_transfer_capacities=dict(h2d=1000.,d2h=500.))],
        price_bindings=[])
    candidate=dict(id='baseline',partitions=[dict(id='all',device='cuda:0')],
        compute_options=dict(shared_h2d_bytes_per_second=1000.,
            per_device_h2d_bytes_per_second={'cuda:0':800.},
            peak_fp32_flops_per_second={'cuda:0':100000.}),
        output_options=dict(shared_d2h_bytes_per_second=500.,
            per_device_d2h_bytes_per_second={'cuda:0':400.},
            output_bytes_per_second=200.),
        mode_service_options=dict(writer_profiles={'cuda:0':deepcopy(writer)}))
    return profile,candidate


def test_matching_values_still_need_declarations_and_gpu_shape_service():
    profile,candidate=case()
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert report['status']=='work_prices_match_but_unbound'
    assert report['unbound_price_leaves']==report['matched_price_leaves']
    assert report['missing_active_prices']==['gpu_shape_service']
    paths=[[0,'profiles','cuda:0'],[0,'shared_capacities'],
           [0,'shared_transfer_capacities']]
    profile['price_bindings']=[dict(targets=[dict(context_path=path)
        for path in paths])]
    declared=audit_staged_work_price_binding(profile,'current',[candidate])
    assert declared['declared_price_leaves']==declared['matched_price_leaves']
    assert declared['status']=='work_prices_match_but_unbound'
    assert declared['missing_active_prices']==['gpu_shape_service']


@pytest.mark.parametrize('path',[
    ('compute_options','per_device_h2d_bytes_per_second','cuda:0'),
    ('compute_options','peak_fp32_flops_per_second','cuda:0'),
    ('output_options','per_device_d2h_bytes_per_second','cuda:0'),
    ('mode_service_options','writer_profiles','cuda:0','writer_copy_service',
     'cpu_seconds_per_call')])
def test_changed_active_work_rate_is_rejected(path):
    profile,candidate=case()
    value=candidate
    for key in path[:-1]:value=value[key]
    value[path[-1]]*=2
    with pytest.raises(ValueError,match='differs from active profile'):
        audit_staged_work_price_binding(profile,'current',[candidate])


def test_shared_scenario_cannot_exceed_context_and_missing_bus_is_visible():
    profile,candidate=case()
    candidate['output_options']['output_bytes_per_second']=100.
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert report['matched_price_leaves']>0
    candidate['output_options']['output_bytes_per_second']=300.
    with pytest.raises(ValueError,match='exceeds active context'):
        audit_staged_work_price_binding(profile,'current',[candidate])
    profile,candidate=case()
    del profile['contexts'][0]['shared_transfer_capacities']['h2d']
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert 'shared_h2d_bytes_per_second' in report['missing_active_prices']
    profile,candidate=case()
    context=profile['contexts'][0]
    context['shared_capacities'].update(context.pop('shared_transfer_capacities'))
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert {'shared_h2d_bytes_per_second','shared_d2h_bytes_per_second'} <= set(
        report['missing_active_prices'])
    profile,candidate=case()
    profile['contexts'][0]['shared_links']=[dict(devices=['cuda:0'],
        h2d_bytes_per_second=1000.,d2h_bytes_per_second=500.)]
    with pytest.raises(ValueError,match='omits active shared transfer links'):
        audit_staged_work_price_binding(profile,'current',[candidate])


def test_multigpu_without_declared_link_topology_is_unbound():
    profile,candidate=case()
    context=profile['contexts'][0]
    context['profiles']['cuda:1']=deepcopy(context['profiles']['cuda:0'])
    candidate['partitions'].append(dict(id='other',device='cuda:1'))
    candidate['compute_options']['per_device_h2d_bytes_per_second']['cuda:1']=800.
    candidate['compute_options']['peak_fp32_flops_per_second']['cuda:1']=100000.
    candidate['output_options']['per_device_d2h_bytes_per_second']['cuda:1']=400.
    candidate['mode_service_options']['writer_profiles']['cuda:1']=deepcopy(
        candidate['mode_service_options']['writer_profiles']['cuda:0'])
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert 'shared_link_topology' in report['missing_active_prices']


@pytest.mark.parametrize('mode', ['significant','jagwas','device_significant'])
def test_reduced_primitive_banks_remain_explicitly_unbound(mode):
    profile,candidate=case()
    active=profile['contexts'][0]['profiles']['cuda:0']
    if mode=='significant':
        candidate['output_options']['significant_backend']='host'
        candidate['mode_service_options']=dict(
            selection_prices={'placeholder':1.},
            selection_profiles={'cuda:0':deepcopy(active)},
            archive_prices={'placeholder':2.},
            archive_profiles={'cuda:0':deepcopy(active)})
    elif mode=='jagwas':
        candidate['output_options']['jagwas_writer_fsync']=True
        candidate['compute_options']['peak_fp64_flops_per_second']={'cuda:0':10000.}
        candidate['mode_service_options']=dict(
            selection_prices={'placeholder':1.},
            cpu_fraction_by_device={'cuda:0':.5},
            archive_price={'placeholder':2.},
            archive_profiles={'cuda:0':deepcopy(active)})
    else:
        candidate['output_options']['significant_backend']='device'
        candidate['mode_service_options']=dict(
            launch_profiles={'cuda:0':dict(kernel_launch_seconds=.0001,
                gpu_fraction=1.)},count_transfer_prices={'placeholder':1.},
            archive_prices={'placeholder':2.},
            archive_profiles={'cuda:0':deepcopy(active)})
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert report['status']=='work_prices_match_but_unbound'
    assert any(item.startswith('reduced_service.')
               for item in report['missing_active_prices'])


def test_changed_reduced_archive_profile_is_rejected():
    profile,candidate=case()
    active=deepcopy(profile['contexts'][0]['profiles']['cuda:0'])
    active['fsync_seconds']*=2
    candidate['output_options']['significant_backend']='host'
    candidate['mode_service_options']=dict(
        selection_prices={},selection_profiles={'cuda:0':deepcopy(active)},
        archive_prices={},archive_profiles={'cuda:0':active})
    with pytest.raises(ValueError,match='differs from active profile'):
        audit_staged_work_price_binding(profile,'current',[candidate])


def test_shape_profile_must_match_active_fixed_prices_and_geometry():
    profile,candidate=case()
    active=profile['contexts'][0]['profiles']['cuda:0']
    active['kernel_geometry']=[dict(N=129,B=4,K=5,C=3,validate_range=True)]
    candidate['shape_profiles']={'cuda:0':deepcopy(active)}
    report=audit_staged_work_price_binding(profile,'current',[candidate])
    assert report['missing_active_prices']==[
        'gpu_shape_service.cuda:0.host_primitives']
    assert [0,'profiles','cuda:0','gpu_resources','fp32_flops_per_second'] in (
        report['unbound_paths'])
    candidate['shape_profiles']['cuda:0']['gpu_resources'][
        'fp32_flops_per_second']*=2
    with pytest.raises(ValueError,match='GPU shape prices differ'):
        audit_staged_work_price_binding(profile,'current',[candidate])
    candidate['shape_profiles']={'cuda:0':deepcopy(active)}
    candidate['shape_profiles']['cuda:0']['kernel_geometry'][0]['B']=8
    with pytest.raises(ValueError,match='GPU shape geometry differs'):
        audit_staged_work_price_binding(profile,'current',[candidate])


def test_shape_prices_can_be_bound_to_immutable_targets():
    profile,candidate=case()
    active=profile['contexts'][0]['profiles']['cuda:0']
    active['host_primitives']=dict(synthetic_dispatch_seconds=.001)
    active['kernel_geometry']=[dict(N=129,B=4,K=5,C=3,
                                   validate_range=True,kernels=[])]
    candidate['shape_profiles']={'cuda:0':deepcopy(active)}
    unbound=audit_staged_work_price_binding(profile,'current',[candidate])
    assert unbound['missing_active_prices']==[]
    assert [0,'profiles','cuda:0','host_primitives',
            'synthetic_dispatch_seconds'] in unbound['unbound_paths']
    profile['price_bindings']=[dict(targets=[dict(context_path=path)
        for path in ([0,'profiles','cuda:0'],[0,'shared_capacities'],
                     [0,'shared_transfer_capacities'])])]
    bound=audit_staged_work_price_binding(profile,'current',[candidate])
    assert bound['status']=='declared_work_prices_verified'
    assert bound['declared_price_leaves']==bound['matched_price_leaves']


@pytest.mark.parametrize('mode', ['jagwas','significant'])
def test_reduced_banks_match_the_original_record_exactly(mode):
    profile,candidate=case()
    active=deepcopy(profile['contexts'][0]['profiles']['cuda:0'])
    if mode=='jagwas':
        candidate['output_options']['jagwas_writer_fsync']=True
        candidate['compute_options']['peak_fp64_flops_per_second']={'cuda:0':10000.}
        candidate['mode_service_options']=dict(selection_prices={'primitive':.01},
            archive_price={'call':.02},cpu_fraction_by_device={'cuda:0':.5},
            archive_profiles={'cuda:0':active})
        value=dict(writer_prices=dict(prices={'primitive':.01},
                                      archive={'call':.02}))
    else:
        candidate['output_options']['significant_backend']='host'
        value=dict(prices={'primitive':.01},archive={'call':.02},
                   protocol='synthetic')
        candidate['mode_service_options']=dict(selection_prices=deepcopy(value),
            archive_prices={'call':.02},selection_profiles={'cuda:0':active},
            archive_profiles={'cuda:0':active})
    evidence=dict(record=dict(value=value),record_sha256='a'*64,
                  artifact_sha256='b'*64)
    report=audit_staged_work_price_binding(profile,'current',[candidate],
        reduction_price_evidence=evidence)
    assert len(report['external_price_fields_bound'])==2
    assert report['external_record_sha256']=='a'*64
    assert report['external_artifact_sha256']=='b'*64
    incomplete=dict(evidence,artifact_sha256=None)
    unbound=audit_staged_work_price_binding(profile,'current',[candidate],
        reduction_price_evidence=incomplete)
    assert not unbound['external_price_fields_bound']
    changed=deepcopy(candidate)
    changed['mode_service_options']['selection_prices']={'other':1.}
    with pytest.raises(ValueError,match='differs from bound reduction record'):
        audit_staged_work_price_binding(profile,'current',[changed],
            reduction_price_evidence=evidence)


def test_device_significant_archive_binds_but_count_latency_remains_unbound():
    profile,candidate=case()
    active=deepcopy(profile['contexts'][0]['profiles']['cuda:0'])
    candidate['output_options']['significant_backend']='device'
    candidate['mode_service_options']=dict(
        launch_profiles={'cuda:0':dict(kernel_launch_seconds=.0001,
            gpu_fraction=1.)},count_transfer_prices={'cuda:0':dict(latency_seconds=.01)},
        archive_prices={'True':dict(call_cpu_seconds=.02)},
        archive_profiles={'cuda:0':active})
    evidence=dict(record=dict(value=dict(
        archive={'True':dict(call_cpu_seconds=.02)})),
        record_sha256='a'*64,artifact_sha256='b'*64)
    report=audit_staged_work_price_binding(profile,'current',[candidate],
        reduction_price_evidence=evidence)
    assert report['external_price_fields_bound']==['significant.archive_prices']
    assert 'reduced_service.count_transfer_prices' in report['missing_active_prices']
