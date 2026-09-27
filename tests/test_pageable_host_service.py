import copy
from unittest.mock import patch

import pytest

from torchgwas.pageable_host_service import pageable_host_prices,setup_host_allocations,setup_host_service
from torchgwas.setup_work import setup_work,setup_service
from torchgwas.trait_tiling_model import _prepare_graph


def primitive():
    size=64<<20;pages=size//4096
    rows=[]
    for operation in ['numpy_copy','pageable_d2h']:
        for repeat in range(7):
            rows.append(dict(operation=operation,repeat=repeat,bytes=size,
                mapping=dict(heap=False,mapping_bytes=size+4096),
                allocation=dict(seconds=.001,cpu_seconds=.001),
                first=dict(seconds=.09,cpu_seconds=.08,minor_faults=pages,major_faults=0),
                warm=dict(seconds=.03,cpu_seconds=.02,minor_faults=0,major_faults=0),
                release=dict(seconds=.01,cpu_seconds=.01),signed_first_less_warm_cpu=.06,signed_first_less_warm_faults=pages))
    return dict(bytes=size,page_bytes=4096,numpy_version='n',torch_version='t',python_version='p',libc=['glibc','g'],
        affinity=[1,2],torch_threads=4,numpy_madvise_hugepage=False,allocator_environment={},source_sha256='s',
        contexts=[dict(devices=['cuda:0'],workers=[dict(device='cuda:0',rows=rows)])])


def prices(probe=None):
    return pageable_host_prices(probe or primitive(),devices=['cuda:0'],device='cuda:0',numpy_version='n',torch_version='t',python_version='p',libc=['glibc','g'],cpu_affinity=[1,2],source_sha256='s')


def inputs(work,route='mmap'):
    p=prices();arrays={'numpy':{},'torch':{}}
    for a in setup_host_allocations(work):
        size=a['bytes'];mapped=(size//4096+1)*4096
        obs=[dict(repeat=i,route=route,allocate_count=1 if route=='mmap' else 0,
            release_count=-1 if route=='mmap' else 0,allocate_bytes=mapped if route=='mmap' else 0,
            release_bytes=-mapped if route=='mmap' else 0) for i in range(4)]
        arrays[a['kind']][str(size)]=dict(route=route,mapped_bytes=mapped if route=='mmap' else None,observations=obs)
    return dict(prices=p,geometry=dict(p['context'],arrays=arrays,durations_recorded=False,page_bytes=4096),arena_fresh_fraction=0.)


def profile(work):
    return dict(setup_primitives={p:dict(reference_shape=[32,1,8],cpu_seconds=.001,non_cpu_seconds=.002)
        for p in ['residual_common','residual_block','design_common','design_block']},cpu_fraction=.5,
        gpu_resources=dict(gpu_fraction=1.,hbm_bytes_per_second=1e10,fp32_flops_per_second=1e12),
        h2d_bytes_per_second=1e10,d2h_bytes_per_second=1e10,process_units=dict(numpy_copy_bytes=99.,covariate_basis_work=0.),
        pin_cpu_seconds_per_page=0.,pin_driver_seconds_per_page=0.,pin_cached_cpu_seconds_per_call=0.,
        pageable_host_service=inputs(work))


def test_price_context_complete_repeats_and_disturbed_controls_retained():
    p=primitive();p['contexts'][0]['workers'][0]['rows'][0]['first']['minor_faults']+=99
    p['contexts'][0]['workers'][0]['rows'][0]['signed_first_less_warm_faults']+=99
    result=prices(p)
    assert len(result['control_warnings'])==1
    assert len(result['prices']['numpy']['repeat_prices'])==7
    assert result['prices']['numpy']['first_touch_cpu_seconds_per_page']==pytest.approx(.06/16384)
    for mutate in [lambda x:x.update(source_sha256='wrong'),lambda x:x.update(numpy_madvise_hugepage=True),
                   lambda x:x['contexts'][0]['workers'][0]['rows'].pop(),
                   lambda x:x['contexts'][0]['workers'][0]['rows'][0]['first'].update(major_faults=1),
                   lambda x:x['contexts'][0]['workers'][0]['rows'][0].update(signed_first_less_warm_cpu=.01)]:
        bad=copy.deepcopy(p);mutate(bad)
        with pytest.raises(ValueError):prices(bad)


def test_source_allocation_lifetimes_and_multiblock_conservation():
    single=setup_host_allocations(setup_work(4096,4096,input_contiguous=False))
    assert [(a['kind'],a['release_phase']) for a in single]==[('numpy',1),('torch','after_scan')]
    full=setup_host_allocations(setup_work(4096,4096,input_contiguous=True))
    assert len(full)==1 and full[0]['kind']=='torch'
    work=setup_work(32768,40001,input_contiguous=True)
    arrays=setup_host_allocations(work)
    assembled=next(a for a in arrays if a['name']=='assembled_result')
    assert assembled['release_phase']=='after_scan'
    assert sum(t['bytes'] for t in assembled['touches'])==4*32768*40001
    assert all(a['release_phase']!='after_scan' for a in arrays if a['kind']=='torch')
    service=setup_host_service(work,inputs(work))
    assert sum(r['touch_cpu_seconds'] for r in service['phases'])==pytest.approx(sum(a['touch_pages']*service_inputs_price(a['kind']) for a in service['allocations']))
    assert service['cleanup']['serial_cpu_seconds']==service['cleanup']['cpu_seconds']>0


def service_inputs_price(kind):return prices()['prices'][kind]['first_touch_cpu_seconds_per_page']


def test_allocator_context_geometry_and_explicit_arena_uncertainty():
    work=setup_work(4096,2048);x=inputs(work,'arena')
    old=copy.deepcopy(x);zero=setup_host_service(work,x)
    assert sum(r['touch_cpu_seconds'] for r in zero['phases'])==0
    x['arena_fresh_fraction']=1.;fresh=setup_host_service(work,x)
    assert sum(r['touch_cpu_seconds'] for r in fresh['phases'])>0
    assert fresh['cleanup']==dict(cpu_seconds=0.,serial_cpu_seconds=0.)
    for value in [None,True,float('nan'),-1,2]:
        bad=copy.deepcopy(x);bad['arena_fresh_fraction']=value
        with pytest.raises(ValueError):setup_host_service(work,bad)
    for mutate in [lambda g:g.update(durations_recorded=True),lambda g:g.update(affinity=[0]),
        lambda g:g['arrays']['torch'][str(4*4096*2048)]['observations'].pop(),
        lambda g:g['arrays']['torch'][str(4*4096*2048)].update(route='mmap')]:
        bad=copy.deepcopy(old);mutate(bad['geometry'])
        with pytest.raises(ValueError):setup_host_service(work,bad)
    bad=inputs(work);bad['geometry']['arrays']['torch']={}
    with pytest.raises(ValueError,match='Missing pageable geometry'):setup_host_service(work,bad)
    bad=inputs(work);row=bad['geometry']['arrays']['torch'][str(4*4096*2048)]
    row['observations'][0]['allocate_bytes']+=4096;row['observations'][0]['release_bytes']-=4096;row['mapped_bytes']=None
    with pytest.raises(ValueError,match='Unstable pageable mapped extent'):setup_host_service(work,bad)


def test_transfer_not_counted_twice_and_graph_conserves_cpu_gil_dram():
    work=setup_work(4096,4096,input_contiguous=False);p=profile(work);result=setup_service(work,p)
    row=result['phases'][1];q=p['cpu_fraction'];extra=4*4096*4096-128
    warm=extra*p['pageable_host_service']['prices']['prices']['torch']['warm_cpu_seconds_per_byte']/q
    assert row['transfer_seconds']==pytest.approx(extra/p['h2d_bytes_per_second']+max(extra/p['d2h_bytes_per_second'],warm))
    assert row['host_copy_seconds']==pytest.approx(.02/q)
    graph,report=_prepare_graph(work,p,1.,0,0)
    for i,cost in enumerate(report['phases']):
        seconds=cost['seconds'];r=graph.demands['phase:'+str(i)]
        assert r['cpu']*seconds==pytest.approx(q*(cost['fixed_cpu_seconds']+cost['host_copy_seconds']+cost['host_page_seconds']+cost['pageable_d2h_cpu_seconds']))
        assert r['host_serial']*seconds==pytest.approx(q*(cost['fixed_cpu_seconds']+cost['host_page_serial_seconds']))
        assert r['dram']*seconds==pytest.approx(cost['host_dram_bytes'])
    assert report['cleanup']['cpu_seconds']>0
    assert report['cleanup']['serial_cpu_seconds']==0  # Tensor storage destruction releases GIL.
    assert result['seconds']==sum(r['seconds'] for r in result['phases'])  # Cleanup belongs after scan.


def test_small_numpy_first_touch_is_held_and_torch_is_detached():
    work=setup_work(32,1,input_contiguous=False);result=setup_host_service(work,inputs(work))
    row=result['phases'][1]
    assert row['serial_cpu_seconds']==pytest.approx(row['touch_cpu_seconds']+row['allocation_cpu_seconds']+row['release_cpu_seconds'])
    assert result['cleanup']['serial_cpu_seconds']==0


def test_cleanup_after_generator_drain_before_writer_close_and_next_tile(tmp_path):
    import numpy as np
    from test_trait_tiling_model import candidate,component
    from test_pgen_native_reader import write_pgen
    from torchgwas.trait_tiling_model import torch_trait_tiled_runtime
    path=tmp_path/'input.pgen';write_pgen(path,np.arange(320,dtype=np.uint8).reshape(10,32)%4)
    c=candidate(path,count=1)
    for tile in c['tiles']:
        work=setup_work(32,tile['data']['traits_analyzed'])
        tile['profile']['pageable_host_service']=inputs(work)
    before=copy.deepcopy(c)
    with patch('torchgwas.mechanistic_torch.tensor_stage_service',side_effect=component):
        graph=torch_trait_tiled_runtime(c,host_serial_fraction=1.,pageable_arena_fresh_fraction=1.,return_graph=True)
    assert c==before
    result=graph.solve()
    for i in range(3):
        prefix=f'tile{i}:';start=result['start'][prefix+'prepare:host_cleanup'];end=result['end'][prefix+'prepare:host_cleanup']
        assert start>=result['end'][prefix+'consume:2']
        assert start>=result['end'][prefix+'release:2']
        assert result['start'][prefix+'writer:close:start']>=end
        if i<2:assert result['start'][f'tile{i+1}:prepare:covariate_basis']>=end
