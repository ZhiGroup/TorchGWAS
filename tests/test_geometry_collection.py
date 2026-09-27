"""Coverage, probe bounds, provenance and publication contracts."""
import copy
import json
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from test_trait_candidate_space import spec
from test_trait_tiling_model import component
from test_pageable_host_service import inputs
from torchgwas.setup_work import setup_work
from torchgwas import geometry_collection as g


def inventory(value,**limits):
    return g.geometry_requirements(**value,**limits)


def empty_kernels(value):
    for context in value['contexts']:
        for profile in context['profiles'].values():profile['kernel_geometry']=[]


def add_host(value):
    for context in value['contexts']:
        for profile in context['profiles'].values():
            data=inputs(setup_work(32,5))
            data['geometry']['arrays']={'numpy':{},'torch':{}}
            data['geometry']['source_sha256']={}
            profile['pageable_host_service']=data


def test_collect_only_missing_exact_shapes_including_variant_tails(tmp_path):
    value=spec(tmp_path);value['bounds']['partition_axes']=['trait','variant']
    empty_kernels(value);before=copy.deepcopy(value)
    result=inventory(value)
    identities=[(r['device'],*r['shape']) for r in result['kernel_requests']]
    assert len(identities)==len(set(identities))
    assert ('cuda:0',32,4,5,8) in identities
    assert ('cuda:1',32,2,5,8) in identities
    assert all(r['contexts'] and r['candidate_indices'] for r in result['kernel_requests'])
    for request in result['kernel_requests']:
        n,b,k,c=request['shape']
        for context in value['contexts']:
            if context['name'] in request['contexts']:
                context['profiles'][request['device']]['kernel_geometry'].append(
                    dict(N=n,B=b,K=k,C=c,validate_range=True,kernels=[]))
    assert not inventory(value)['kernel_requests']
    assert before['contexts'][0]['profiles']['cuda:0']['kernel_geometry']==[]


def test_memory_rejection_precedes_probe_and_host_extent_dedup(tmp_path):
    value=spec(tmp_path);empty_kernels(value);add_host(value)
    value['joint']['device_memory_bytes']['cuda:1']=1
    result=inventory(value)
    assert result['memory_rejected']
    assert {r['device'] for r in result['kernel_requests']}=={'cuda:0'}
    ids=[(r['kind'],r['bytes']) for r in result['host_requests']]
    assert len(ids)==len(set(ids)) and ids
    assert result['host_touch_bytes']==4*sum(r['bytes'] for r in result['host_requests'])
    assert all(t[1]=='cuda:0' for r in result['host_requests'] for t in r['targets'])
    value['joint']['device_memory_bytes']['cuda:0']=1
    with pytest.raises(ValueError,match='No memory-feasible'):inventory(value)


@pytest.mark.parametrize('limits,match',[
    ({'max_kernel_shapes':1},'max_kernel_shapes'),({'max_device_probe_bytes':1},'max_device_probe_bytes'),
    ({'max_host_extents':1},'max_host_extents'),({'max_host_probe_bytes':1},'max_host_probe_bytes'),
    ({'max_host_touch_bytes':1},'max_host_touch_bytes'),({'max_kernel_shapes':True},'max_kernel_shapes')])
def test_all_probe_limits_checked_before_any_collection(tmp_path,limits,match):
    value=spec(tmp_path);empty_kernels(value);add_host(value)
    with patch.object(g,'_run_worker') as worker:
        with pytest.raises(ValueError,match=match):inventory(value,**limits)
        worker.assert_not_called()


def test_ambiguous_geometry_is_not_silently_replaced(tmp_path):
    value=spec(tmp_path)
    bank=value['contexts'][0]['profiles']['cuda:0']['kernel_geometry'];bank.append(copy.deepcopy(bank[0]))
    with pytest.raises(ValueError,match='Ambiguous'):inventory(value)


def test_kernel_census_keeps_work_but_strips_every_timing_field():
    events=[dict(cat='cpu',name='discard',dur=10),dict(cat='kernel',name='gemm',ts=15,dur=3,
        args={'grid':[2,1,1],'block':[128,1,1],'registers per thread':32,'shared memory':0,'duration':5,'stream':7})]
    rows=g.kernel_census(events)
    assert rows==[dict(name='gemm',geometry={'grid':[2,1,1],'block':[128,1,1],'registers per thread':32,'shared memory':0})]
    for bad in [[],[dict(cat='kernel',name='x',args={'grid':[0,1,1],'block':[1,1,1]})]]:
        with pytest.raises(ValueError):g.kernel_census(bad)


def test_atomic_records_never_replace_existing_or_publish_nan(tmp_path):
    path=tmp_path/'record.json';g.write_record(path,{'old':1})
    with pytest.raises(FileExistsError):g.write_record(path,{'new':2})
    assert json.loads(path.read_text())=={'old':1}
    with pytest.raises(ValueError):g.write_record(tmp_path/'bad.json',{'x':float('nan')})
    assert sorted(p.name for p in tmp_path.iterdir())==['record.json']


def controller(value):
    return SimpleNamespace(profile=dict(contexts=copy.deepcopy(value['contexts']),component_artifacts={},limitations=[]),
        config={key:value[key] for key in ['bounds','joint']},input_path=value['workload']['genotype'],
        output_path=str(value['workload']['genotype'])+'.out',context={'devices':{}},devices=['cuda:0','cuda:1'])


def test_inspection_creates_nothing_and_partial_failure_publishes_no_profile(tmp_path):
    value=spec(tmp_path);empty_kernels(value);control=controller(value);root=tmp_path/'capture'
    with patch('torchgwas.detailed_autotune.DetailedAutotune',return_value=control),patch.object(g,'_run_worker',side_effect=RuntimeError('worker failed')) as worker:
        result=g.complete_profile_geometry({},value['workload'],{},output=value['output'],output_path=control.output_path,collection_dir=root,collect=False)
        assert result['kernel_requests'] and not root.exists();worker.assert_not_called()
        with pytest.raises(RuntimeError,match='worker failed'):
            g.complete_profile_geometry({},value['workload'],{},output=value['output'],output_path=control.output_path,collection_dir=root)
        assert (root/'coverage.json').is_file() and not (root/'profile.json').exists()


@pytest.mark.parametrize('bound_prices',[False,True])
def test_completed_profile_preserves_prices_and_reuses_exact_planner(tmp_path,bound_prices):
    value=spec(tmp_path);empty_kernels(value);control=controller(value)
    if bound_prices:control.profile['price_bindings']=[dict(record='original measurement declaration')]
    before=copy.deepcopy(control.profile);root=tmp_path/'capture'
    def worker(kind,requests,directory,*args):
        assert kind=='gpu'
        rows=[dict(device=r['device'],**dict(zip(['N','B','K','C'],r['shape'])),validate_range=True,
                   kernels=[dict(name='synthetic',geometry={'grid':[1,1,1],'block':[1,1,1]})]) for r in requests]
        path=directory/'gpu_geometry.json';g.write_record(path,dict(rows=rows));return dict(rows=rows),path
    def bind(contexts,current,**kwargs):
        assert g._price_contexts(contexts)==g._price_contexts(before['contexts'])
        assert kwargs['price_bindings']==before.get('price_bindings')
        return dict(contexts=contexts,execution_context=current,**kwargs)
    with patch('torchgwas.detailed_autotune.DetailedAutotune',return_value=control),patch.object(g,'_run_worker',side_effect=worker),patch.object(g,'execution_context',return_value=control.context),patch.object(g,'validate_detailed_profile'),patch.object(g,'bind_detailed_profile',side_effect=bind),patch.object(g,'write_detailed_profile',side_effect=lambda p,path:g.write_record(path,p)),patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        result=g.complete_profile_geometry({},value['workload'],{},output=value['output'],output_path=control.output_path,collection_dir=root)
    assert result['component_prices_unchanged'] and not result['durations_retained']
    assert result['kernel_shapes_collected']>0 and (root/'manifest.json').is_file()
    assert control.profile==before
    completed=json.loads((root/'profile.json').read_text())['contexts']
    assert g._price_contexts(completed)==g._price_contexts(before['contexts'])
    assert not inventory(dict(value,contexts=completed))['kernel_requests']


def test_host_artifact_rejects_hidden_duration_fields():
    record=dict(inputs(setup_work(32,5))['geometry'],source_sha256={},scope='test')
    g.validate_host_capture(record)
    for level in ['root','row','observation']:
        bad=copy.deepcopy(record)
        target=bad if level=='root' else next(iter(bad['arrays']['torch'].values()))
        if level=='observation':target=target['observations'][0]
        target['seconds']=1
        with pytest.raises(ValueError,match='Unexpected allocator'):g.validate_host_capture(bad)

@pytest.mark.parametrize('bad',['duration','duplicate','missing','unknown_path'])
def test_bad_or_unsupported_capture_never_publishes_profile(tmp_path,bad):
    value=spec(tmp_path);empty_kernels(value);control=controller(value);root=tmp_path/'capture'
    def worker(kind,requests,directory,*args):
        rows=[dict(device=r['device'],**dict(zip(['N','B','K','C'],r['shape'])),validate_range=True,
                   kernels=[dict(name='unknown',geometry={'grid':[1,1,1],'block':[1,1,1]})]) for r in requests]
        if bad=='duration':rows[0]['kernels'][0]['geometry']['duration']=.01
        if bad=='duplicate':rows.append(copy.deepcopy(rows[0]))
        if bad=='missing':rows.pop()
        path=directory/'gpu_geometry.json';g.write_record(path,dict(rows=rows));return dict(rows=rows),path
    with patch('torchgwas.detailed_autotune.DetailedAutotune',return_value=control),patch.object(g,'_run_worker',side_effect=worker):
        with pytest.raises(ValueError):
            g.complete_profile_geometry({},value['workload'],{},output=value['output'],output_path=control.output_path,collection_dir=root)
    assert not (root/'profile.json').exists()
    assert not (root/'manifest.json').exists()
