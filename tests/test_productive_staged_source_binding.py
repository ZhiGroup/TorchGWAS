"""A live staged screen must not silently borrow arbitrary source rates."""
from copy import deepcopy

import pytest

from torchgwas.productive_staged_source_binding import audit_staged_source_price_binding


def case():
    source=dict(decode_units=dict(copy_packed_byte=.01,uleb1=.02),
        cpu_fraction=.5,depth=2,decode_workers=2,cpu_available_cores=4.,
        shared_dram_bytes_per_second=1e9,read_bytes_per_second=1e8)
    context=dict(name='current',shared_capacities=dict(cpu=4.,dram=1e9,
        input=1e8,output=2e8),profiles={'cuda:0':source})
    profile=dict(contexts=[context],price_bindings=[])
    candidate=dict(id='baseline',partitions=[dict(id='only',device='cuda:0')],
        source_profiles={'only':deepcopy(source)})
    capacities=dict(cpu=2.,dram=5e8,input=5e7)
    return profile,candidate,capacities


def test_matching_source_values_are_not_automatically_measurement_bound():
    profile,candidate,capacities=case()
    report=audit_staged_source_price_binding(profile,'current',[candidate],capacities)
    assert report['status']=='source_prices_match_but_unbound'
    assert report['matched_price_leaves']==report['unbound_price_leaves']
    assert report['declared_price_leaves']==0
    assert any(path[-1]=='copy_packed_byte' for path in report['unbound_paths'])


def test_declared_parent_targets_cover_all_source_leaves():
    profile,candidate,capacities=case()
    paths=[[0,'profiles','cuda:0',field] for field in candidate[
        'source_profiles']['only']]
    paths += [[0,'shared_capacities',field] for field in ('cpu','dram','input')]
    profile['price_bindings']=[dict(targets=[dict(context_path=path)
        for path in paths])]
    report=audit_staged_source_price_binding(profile,'current',[candidate],capacities)
    assert report['status']=='declared_source_prices_verified'
    assert report['unbound_price_leaves']==0
    assert report['declared_price_leaves']==report['matched_price_leaves']


def test_changed_or_unrecorded_source_price_cannot_bind():
    profile,candidate,capacities=case()
    candidate['source_profiles']['only']['decode_units']['uleb1']*=2
    with pytest.raises(ValueError,match='differs from active profile'):
        audit_staged_source_price_binding(profile,'current',[candidate],capacities)
    profile,candidate,capacities=case()
    candidate['source_profiles']['only']['input_read_cpu_prices']=dict(
        cpu_seconds_per_byte=.01,cpu_seconds_per_call=.02)
    with pytest.raises(ValueError,match='input_read_cpu_prices'):
        audit_staged_source_price_binding(profile,'current',[candidate],capacities)
    profile,candidate,capacities=case()
    profile['contexts'][0]['profiles']['cuda:0']['input_read_cpu_prices']=dict(
        cpu_seconds_per_byte=.01,cpu_seconds_per_call=.02)
    with pytest.raises(ValueError,match='omits active buffered-read'):
        audit_staged_source_price_binding(profile,'current',[candidate],capacities)
    profile,candidate,capacities=case()
    capacities['cpu']=5.
    with pytest.raises(ValueError,match='exceeds active context'):
        audit_staged_source_price_binding(profile,'current',[candidate],capacities)
