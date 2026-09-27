"""Prepared-window composition keeps real writer scopes and shared resources."""
from copy import deepcopy
from unittest.mock import patch
import pytest

from test_trait_tiling_model import candidate,input_path
from test_mechanistic_shapes import component
from test_significant_host_model import bank
from torchgwas.pgen_work_bounds import PgenHeaderWork
from torchgwas.window_model import prepared_window_runtime,compare_prepared_windows
from torchgwas.significant_host_work import indexed_part_work


def configuration(path,*,reduction=None,count=2,traits=5,width=2):
    choice=candidate(path,width=width,count=count,traits=traits,block_bytes=None)
    header=PgenHeaderWork(path);windows=[]
    for tile in choice['tiles']:
        window=deepcopy(tile);window['issued_chunks']=4
        window['data']['encoded']=header.window(0,10,4)
        windows.append(window)
    options=dict(total_traits=traits,reduction=reduction,output=choice['output'],
        shared_capacities=choice['shared_capacities'],endpoint='upper',host_serial_fraction=.5)
    if reduction:
        options.update(prices=bank(),retained=[[r['markers']*w['data']['traits_analyzed'] for r in w['data']['encoded']['chunks']] for w in windows])
    return windows,options


def evaluate(windows,options,**kw):
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return prepared_window_runtime(windows,**dict(options,**kw))


@pytest.mark.parametrize('mode',[None,'significant'])
def test_vectorized_source_window_preserves_output_graph(input_path,mode):
    windows,options=configuration(input_path,reduction=mode)
    baseline=evaluate(windows,options,return_graph=True)
    header=PgenHeaderWork(input_path)
    altered=deepcopy(windows)
    for row in altered:
        row['data']['encoded']=header.vectorized_window(0,10,4)
    vectorized=evaluate(altered,options,return_graph=True)
    assert vectorized.nodes==baseline.nodes
    assert vectorized.demands==baseline.demands
    assert vectorized.token_capacities==baseline.token_capacities
    assert vectorized.solve()['seconds']==baseline.solve()['seconds']


def test_window_selector_prices_follow_explicit_runtime_choice(input_path, monkeypatch):
    windows, options = configuration(input_path, reduction='significant')
    monkeypatch.setenv('TORCHGWAS_HOST_PREDICATE', 'native')
    with pytest.raises(ValueError, match='host selector'):
        evaluate(windows, options)
    options['prices']['host_selector'] = 'native_row_flat_v2'
    options['prices']['prices']['predicate_native'] = options['prices']['prices'].pop('predicate_block')
    report = evaluate(windows, options)
    assert report['estimated_window_seconds'] > 0


@pytest.mark.parametrize('count',[1,2])
@pytest.mark.parametrize('block_bytes',[None,48])
def test_dense_payload_df_and_sequential_device_windows(input_path,count,block_bytes):
    windows,options=configuration(input_path,count=count);options['output']['block_bytes']=block_bytes
    original=deepcopy(windows);graph=evaluate(windows,options,return_graph=True);solved=graph.solve()
    report=evaluate(windows,options)
    assert windows==original
    assert report['payload_bytes']==4*10*(3*5+3)
    assert sum(r['fsync_calls'] for r in report['windows'])==12  # beta, t, -log10 P, df per window
    assert report['source_chunks']==9 and report['estimated_window_seconds']==solved['seconds']
    assert solved['start'][f'window:{count}:submit_decode:0']>=solved['end']['window:0:complete']
    for i in range(3):
        assert solved['end'][f'window:{i}:writer:df:fsync']<=solved['seconds']
    assert not any(':prepare:' in name for name in graph.nodes)


@pytest.mark.parametrize('count',[1,2])
@pytest.mark.parametrize('empty',[False,True])
def test_significant_parts_share_writer_and_empty_keeps_dense_transfers(input_path,count,empty):
    windows,options=configuration(input_path,reduction='significant',count=count)
    if empty:options['retained']=[[0]*3 for _ in windows]
    report=evaluate(windows,options);graph=evaluate(windows,options,return_graph=True);solved=graph.solve()
    expected=sum(indexed_part_work(v)['file_bytes'] for counts in options['retained'] for v in counts)
    assert report['payload_bytes']==expected
    assert sum(r['parts'] for r in report['windows'])==(0 if empty else 9)
    assert sum(r['d2h_bytes'] for r in report['windows'])==550
    assert ('significant:writer' in graph.token_capacities)==(count>1)
    assert solved['seconds']==solved['end']['significant:finalize:done']
    assert all(r['fsync_calls']==r['parts'] for r in report['windows'])


def test_explicit_occupancy_not_inferred_from_significance_threshold(input_path):
    windows,options=configuration(input_path,reduction='significant')
    a=evaluate(windows,options,significance_threshold=1e-8)
    b=evaluate(windows,options,significance_threshold=1.)
    assert a['payload_bytes']==b['payload_bytes']
    with pytest.raises(ValueError,match='survivors'):evaluate(windows,options,retained=None)


@pytest.mark.parametrize('mode',[None,'significant'])
def test_writer_and_selection_services_match_existing_exact_calculator(input_path,mode):
    from torchgwas.trait_tiling_model import torch_trait_tiled_runtime
    from torchgwas.significant_host_model import significant_host_runtime
    windows,options=configuration(input_path,reduction=mode)
    for w in windows:w['issued_chunks']=0
    choice=candidate(input_path,block_bytes=None)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        old=(torch_trait_tiled_runtime(choice,host_serial_fraction=.5,return_graph=True) if mode is None else
            significant_host_runtime(choice,bank(),occupancy=options['retained'],host_serial_fraction=.5,return_graph=True))
    new=evaluate(windows,options,return_graph=True);compared=0
    for name,(seconds,_) in old.nodes.items():
        if mode is None:
            if ':writer:' not in name:continue
            first,rest=name.split(':',1);mapped='window:'+first.removeprefix('tile')+':'+rest
        else:
            if ':significant:' not in name or not any(s in name for s in (':select:',':write:')):continue
            mapped=name
        assert new.nodes[mapped][0]==seconds
        assert new.demands.get(mapped,{})==old.demands.get(name,{})
        compared+=1
    assert compared>30


def test_shared_storage_link_work_and_capacity_pressure(input_path):
    windows,options=configuration(input_path)
    links=[dict(devices=['cuda:0','cuda:1'],h2d_bytes_per_second=1e8,d2h_bytes_per_second=1e8)]
    a=evaluate(windows,options,shared_links=links,shared_storage_bytes_per_second=1e8)
    resources=a['resource_balance']['resources']
    assert resources['storage']['work']==pytest.approx(resources['input']['work']+resources['output']['work'])
    for direction in ('h2d','d2h'):
        assert resources['link:0:'+direction]['work']==pytest.approx(sum(resources[d+':'+direction]['work'] for d in ['cuda:0','cuda:1']))
    caps=deepcopy(options['shared_capacities']);caps['output']=1.
    original=deepcopy(windows)
    b=evaluate(windows,dict(options,shared_capacities=caps))
    assert windows==original
    assert b['estimated_window_seconds']>=b['payload_bytes']
    assert b['estimated_window_seconds']>a['estimated_window_seconds']
    assert b['resource_balance']['largest_resource_bounds']==['output']


@pytest.mark.parametrize('resource',['cpu','dram','input','output'])
def test_conditional_shared_capacity_preserves_immutable_window_prices(input_path,resource):
    windows,options=configuration(input_path,count=1)
    original=deepcopy(windows)
    caps=dict(options['shared_capacities']);caps[resource]*=.5
    graph=evaluate(windows,dict(options,shared_capacities=caps),return_graph=True)
    assert graph.capacities[resource]==caps[resource]
    assert windows==original
    compared=compare(windows,deepcopy(windows),dict(options,shared_capacities=caps))
    contract=compared['comparison_contract']
    assert contract['shared_capacities']==caps
    assert contract['baseline']['windows'][0]['profile_sha256']==contract['candidate']['windows'][0]['profile_sha256']
    assert compared['estimated_window_gain_seconds']==0
    caps[resource]=options['shared_capacities'][resource]*1.01
    with patch('torchgwas.window_model.torch_scan_header_work',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError,match='exceeds window profile capacity'):
            prepared_window_runtime(windows,**dict(options,shared_capacities=caps))


def test_full_candidate_accepts_only_lower_conditional_capacity(input_path):
    from torchgwas.trait_tiling_model import trait_tiled_shape
    choice=candidate(input_path)
    original=deepcopy(choice)
    assert trait_tiled_shape(choice)==trait_tiled_shape(original)
    choice['shared_capacities']['cpu']*=.5
    assert trait_tiled_shape(choice)==trait_tiled_shape(original)
    assert [tile['profile'] for tile in choice['tiles']]==[tile['profile'] for tile in original['tiles']]
    choice['shared_capacities']['cpu']=original['shared_capacities']['cpu']*1.01
    with pytest.raises(ValueError,match='exceeds tile profile'):
        trait_tiled_shape(choice)


@pytest.mark.parametrize('change',[
    lambda w,o:o.update(max_source_chunks=2),lambda w,o:o.update(max_records=10),
    lambda w,o:o.update(max_windows=1),lambda w,o:o['output'].update(block_bytes=1),
    lambda w,o:w[1]['data']['encoded'].update(input_identity={}),
    lambda w,o:w[0]['profile'].update(write_bytes_per_second=123.),
    lambda w,o:w[1].update(trait_range=[0,2]),
    lambda w,o:o.update(host_serial_fraction=float('nan'))])
def test_preflight_rejects_budget_binding_and_overlap_before_pricing(input_path,change):
    windows,options=configuration(input_path);change(windows,options)
    with patch('torchgwas.window_model.torch_scan_header_work',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError):prepared_window_runtime(windows,**options)


def test_huge_coalesced_writeback_is_rejected_before_graph_expansion(input_path):
    windows,options=configuration(input_path,count=1,traits=1,width=1)
    windows[0]['trait_range']=[0,10**9];windows[0]['data']['traits_analyzed']=10**9
    options.update(total_traits=10**9);options['output']['block_bytes']=10**12
    with patch('torchgwas.window_model.torch_scan_header_work',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError,match='writeback'):prepared_window_runtime(windows,**options)


def test_real_jagwas_projection_narrow_output_and_no_trait_partition(tmp_path):
    import numpy as np
    from test_pgen_native_reader import write_pgen
    from test_jagwas_actual_candidate import actual_candidate,writer_prices
    path=tmp_path/'joint.pgen';write_pgen(path,(np.arange(1025*2049,dtype=np.uint32).reshape(1025,2049)%3).astype(np.uint8))
    choice=actual_candidate(path,128,2);header=PgenHeaderWork(path);windows=[];counts=[]
    for tile in choice['tiles']:
        window={key:deepcopy(tile[key]) for key in ('device','trait_range','data','profile')};window['issued_chunks']=4
        start=tile['variant_range'][0];stop=start+128
        window['data'].update(markers=128,encoded=header.window(start,stop,128));windows.append(window);counts.append([128])
    options=dict(total_traits=512,reduction='jagwas',output=choice['output'],shared_capacities=choice['shared_capacities'],
        endpoint='upper',host_serial_fraction=.5,prices=writer_prices(),retained=counts)
    result=prepared_window_runtime(windows,**options)
    assert sum(r['d2h_bytes'] for r in result['windows'])==17*256
    assert sum(r['parts'] for r in result['windows'])==2
    assert result['payload_bytes']>256*12
    lower_caps=dict(options['shared_capacities']);lower_caps['cpu']*=.5
    conditional=prepared_window_runtime(windows,**dict(options,shared_capacities=lower_caps))
    assert conditional['payload_bytes']==result['payload_bytes']
    assert [w['trait_range'] for w in windows]==[[0,512],[0,512]]
    evidence=survivors(windows,'jagwas',512,lambda w,r:r['markers'])
    common=dict(options);common.pop('retained')
    comparison=compare_prepared_windows(dict(windows=windows,partition_axis='variant'),
        dict(windows=deepcopy(windows),partition_axis='variant'),survivor_evidence=evidence,**common)
    assert comparison['estimated_window_gain_seconds']==0
    assert sum(r['retained'] for r in comparison['candidate']['windows'])==256
    windows[0]['trait_range']=[0,256];windows[0]['data']['traits_analyzed']=256
    with pytest.raises(ValueError,match='complete phenotype'):prepared_window_runtime(windows,**options)


def compare(before,after,options,**kw):
    common=dict(options);common.pop('retained',None)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        return compare_prepared_windows(dict(windows=before,partition_axis='trait'),
            dict(windows=after,partition_axis='trait'),**dict(common,**kw))


def rechunk(path,windows,size):
    result=deepcopy(windows);header=PgenHeaderWork(path)
    for w in result:
        lo,hi=w['data']['encoded']['variant_range']
        w['data']['encoded']=header.window(lo,hi,size);w['profile']['chunk_markers']=size
        # Durations are patched by these composition tests. Real geometry is
        # exercised separately by the JAGWAS case and remote arithmetic audit.
        w['profile']['kernel_geometry']=[dict(N=32,B=r['markers'],K=w['data']['traits_analyzed'],C=8,
            validate_range=True,kernels=[]) for r in w['data']['encoded']['chunks']]
    return result


def survivors(windows,reduction,traits,count,threshold=None):
    return dict(input_identity=deepcopy(windows[0]['data']['encoded']['input_identity']),reduction=reduction,
        total_traits=traits,significance_threshold=threshold,bins=[dict(variant_range=deepcopy(r['variant_range']),
            trait_range=deepcopy(w['trait_range']),retained=count(w,r)) for w in windows for r in w['data']['encoded']['chunks']])


def test_equal_work_compares_different_tiles_devices_and_chunkening(input_path):
    before,options=configuration(input_path,count=1,traits=4,width=2)
    after,_=configuration(input_path,count=2,traits=4,width=1);after=rechunk(input_path,after,8)
    original=deepcopy((before,after,options));result=compare(before,after,options)
    assert (before,after,options)==original
    assert result['coverage']==[dict(variant_range=[0,10],trait_ranges=[[0,4]])]
    assert result['baseline']['source_chunks']==6 and result['candidate']['source_chunks']==8
    # Same associations, with one df stream per phenotype tile.
    assert result['baseline']['payload_bytes']==4*10*(3*4+2)
    assert result['candidate']['payload_bytes']==4*10*(3*4+4)
    assert result['calculation_wall_seconds']>0 and result['calculation_cpu_seconds']>0
    assert not result['selection_validated']


def test_equal_dense_work_can_compare_trait_and_variant_partitions(input_path):
    before,options=configuration(input_path)
    full,_=configuration(input_path,traits=5,width=5);header=PgenHeaderWork(input_path);after=[]
    for i,(lo,hi) in enumerate(((0,4),(4,10))):
        w=deepcopy(full[0]);w['device']=f'cuda:{i}'
        w['data'].update(markers=hi-lo,encoded=header.window(lo,hi,4));after.append(w)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        report=compare_prepared_windows(dict(windows=before,partition_axis='trait'),
            dict(windows=after,partition_axis='variant'),**options)
    assert report['candidate']['payload_bytes']==4*10*(3*5+1)
    assert report['coverage']==[dict(variant_range=[0,10],trait_ranges=[[0,5]])]
    after[1]['device']='cuda:0'
    with pytest.raises(ValueError,match='one window'):evaluate(after,options,partition_axis='variant')


@pytest.mark.parametrize('axis',['',False,'variant','unknown'])
def test_significant_rejects_unsupported_partition_axis(input_path,axis):
    windows,options=configuration(input_path,reduction='significant')
    with pytest.raises(ValueError,match='Partition axis'):evaluate(windows,options,partition_axis=axis)


def test_trait_axis_cannot_reuse_one_trait_range_across_variant_windows(input_path):
    windows,options=configuration(input_path,traits=4)
    header=PgenHeaderWork(input_path)
    for i,w in enumerate(windows):
        w['trait_range']=[0,2];w['data'].update(markers=4,encoded=header.window(4*i,4*i+4,4))
    with pytest.raises(ValueError,match='Trait partitions overlap'):evaluate(windows,options)


@pytest.mark.parametrize('mutation',[
    lambda w:w.pop(),lambda w:w.append(deepcopy(w[0])),
    lambda w:w[0].update(trait_range=[2,0]),lambda w:w[0].update(trait_range=[False,2]),
    lambda w:w[0]['data'].update(samples=33),
    lambda w:w[0]['data']['encoded'].update(input_identity={})])
def test_comparison_refuses_coverage_and_binding_changes_before_pricing(input_path,mutation):
    before,options=configuration(input_path);after=deepcopy(before);mutation(after)
    with patch('torchgwas.window_model.prepared_window_runtime',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError):compare(before,after,options)


def test_comparison_coalesces_exact_survivors_and_prices_part_headers(input_path):
    before,options=configuration(input_path,reduction='significant',traits=4)
    after=rechunk(input_path,before,8)
    evidence=survivors(before,'significant',4,lambda w,r:1,threshold=.05)
    result=compare(before,after,options,survivor_evidence=evidence,significance_threshold=.05)
    assert sum(r['retained'] for r in result['baseline']['windows'])==6
    assert sum(r['retained'] for r in result['candidate']['windows'])==6
    assert result['baseline']['payload_bytes']==6*indexed_part_work(1)['file_bytes']
    assert result['candidate']['payload_bytes']==2*(indexed_part_work(2)['file_bytes']+indexed_part_work(1)['file_bytes'])
    assert result['candidate']['payload_bytes']<result['baseline']['payload_bytes']


@pytest.mark.parametrize('keep',[0,40])
def test_empty_and_full_evidence_can_split_across_new_chunks_and_tiles(input_path,keep):
    before,options=configuration(input_path,reduction='significant',traits=4,width=4)
    after,_=configuration(input_path,reduction='significant',traits=4,width=1)
    evidence=survivors(before,'significant',4,lambda w,r:0)
    evidence['bins']=[dict(variant_range=[0,10],trait_range=(0,4),retained=keep)]
    result=compare(before,after,options,survivor_evidence=evidence)
    assert sum(r['retained'] for r in result['candidate']['windows'])==keep
    assert sum(r['retained'] for r in result['baseline']['windows'])==keep


def test_partial_survivor_bin_cannot_be_guessed_when_chunk_splits_it(input_path):
    before,options=configuration(input_path,reduction='significant',traits=4)
    after=rechunk(input_path,before,2);evidence=survivors(before,'significant',4,lambda w,r:1)
    with patch('torchgwas.window_model.prepared_window_runtime',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError,match='finer evidence'):compare(before,after,options,survivor_evidence=evidence)


@pytest.mark.parametrize('mutation',[
    lambda e:e.update(input_identity={}),lambda e:e.update(reduction='jagwas'),
    lambda e:e.update(total_traits=5),lambda e:e.update(significance_threshold=.05),
    lambda e:e['bins'].pop(),lambda e:e['bins'].append(deepcopy(e['bins'][0]))])
def test_changed_survivor_evidence_is_not_reused(input_path,mutation):
    before,options=configuration(input_path,reduction='significant',traits=4)
    evidence=survivors(before,'significant',4,lambda w,r:1);mutation(evidence)
    with patch('torchgwas.window_model.prepared_window_runtime',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError):compare(before,deepcopy(before),options,survivor_evidence=evidence)


@pytest.mark.parametrize('budget',[dict(max_windows=1),dict(max_source_chunks=2),dict(max_survivor_bins=1)])
def test_comparison_budget_rejected_before_pricing(input_path,budget):
    before,options=configuration(input_path,reduction='significant',traits=4)
    evidence=survivors(before,'significant',4,lambda w,r:1)
    with patch('torchgwas.window_model.prepared_window_runtime',side_effect=AssertionError('preflight first')):
        with pytest.raises(ValueError):compare(before,deepcopy(before),options,survivor_evidence=evidence,**budget)
