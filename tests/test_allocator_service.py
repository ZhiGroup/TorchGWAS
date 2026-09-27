import copy
import pytest
from torchgwas.allocator_service import allocator_service_prices,owned_allocator_service,validate_allocator_geometry
from torchgwas.owned_result_work import owned_result_work


def probe():
    extent=1<<26;page=4096
    def sample(cpu,detached=0.,intervals=0):
        return dict(cpu_seconds=cpu,detached_cpu_seconds=detached,detached_intervals=intervals,balance_errors=0)
    worker=dict(worker=0,context_verified=True,controls={'native':sample(10.,9.5,1),'python':sample(1.)},meter_rows=[],rows=[])
    for repeat in range(3):
        for recording in [False,True]:
            worker['meter_rows'].append(dict(sample(20. if recording else 10.,10. if recording else 0.,20000 if recording else 0),repeat=repeat,recording=recording,intervals=20000))
            for iteration in range(2):
                for size in [0,extent]:
                    worker['rows'].append(dict(worker=0,repeat=repeat,iteration=iteration,recording=recording,bytes=size,empty_cpu_seconds=.1,
                        allocate=sample((3. if size else 1.)+.1),release=sample((7. if size else 1.)+.1),
                        copy=sample(10.,9.5 if recording and size else 0.,1 if recording and size else 0),
                        mapping=dict(allocate_count=int(bool(size)),allocate_bytes=extent+page if size else 0,
                            release_count=-int(bool(size)),release_bytes=-extent-page if size else 0)))
    return dict(context_verified=True,numpy_version='2.2.6',libc=['glibc','2.35'],affinity=[12,13],numpy_madvise_hugepage=False,
        bytes=extent,page_bytes=page,results=[dict(worker_count=1,workers=[worker])])


def prices(value):
    return allocator_service_prices(value,workers=1,numpy_version='2.2.6',libc=['glibc','2.35'],cpu_affinity=[12,13])


def test_independent_prices_remove_tiny_baseline_and_keep_expensive_samples():
    value=probe();result=prices(value)
    assert result['allocate_cpu_seconds_per_mmap']==pytest.approx(2.)
    assert result['release_cpu_seconds_per_mapped_page']==pytest.approx(6./16385)
    # One costly call per repeat must contribute to additive service.
    for row in value['results'][0]['workers'][0]['rows']:
        if row['bytes'] and row['recording'] and row['iteration']==1:row['release']['cpu_seconds']+=8.
    assert prices(value)['release_cpu_seconds_per_mapped_page']==pytest.approx(10./16385)


@pytest.mark.parametrize('fault',['context','missing','duplicate','mapping','held','copy'])
def test_bad_allocator_evidence_rejected(fault):
    value=probe();worker=value['results'][0]['workers'][0]
    if fault=='context':value['numpy_madvise_hugepage']=True
    if fault=='missing':worker['rows'].pop()
    if fault=='duplicate':worker['rows'].append(copy.deepcopy(worker['rows'][0]))
    if fault=='mapping':worker['rows'][1]['mapping']['release_bytes']=0
    if fault=='held':worker['rows'][5]['release'].update(detached_intervals=1,detached_cpu_seconds=.1)
    if fault=='copy':worker['rows'][5]['copy']['detached_intervals']=0
    with pytest.raises(ValueError):prices(value)


def test_array_roles_conserve_mapped_service_without_copy_double_counting():
    w=owned_result_work(1024,8192)
    geometry={'array_bytes':{str(size):dict(route='mmap' if size>=1<<25 else 'arena',mapped_bytes=size+4096 if size>=1<<25 else None) for size in w['array_bytes'].values()}}
    r=owned_allocator_service(w,prices(probe()),geometry)
    assert r['allocate_cpu_seconds']==pytest.approx(4.)
    assert r['worker_release_cpu_seconds']==0.
    assert r['consumer_release_cpu_seconds']==pytest.approx(2*8193*6/16385)
    assert r['arrays']['beta']['release_owner']=='consumer'
    assert r['arrays']['status']['release_owner']=='finish_worker'
    assert r['unpriced_terms']
    del geometry['array_bytes']['1024']
    with pytest.raises(ValueError,match='Missing'):owned_allocator_service(w,prices(probe()),geometry)


def test_noisy_control_difference_is_retained_but_negative_service_is_rejected():
    value=probe()
    for row in value['results'][0]['workers'][0]['rows']:
        if row['repeat']==0 and not row['bytes']:
            row['allocate']['cpu_seconds']+=5.
    result=prices(value)['observations']['allocate']
    assert result['repeat_mean_seconds']==pytest.approx([-3.,2.,2.])
    assert result['paired_repeat_mean_seconds'][0]['recording_disabled']==pytest.approx(-3.)
    assert result['paired_relative_deltas'][0] is None
    assert result['additional_cpu_seconds']==pytest.approx(2.)
    for row in value['results'][0]['workers'][0]['rows']:
        if row['repeat']==1 and not row['bytes']:
            row['allocate']['cpu_seconds']+=5.
    with pytest.raises(ValueError,match='Negative aggregate'):prices(value)


def test_allocator_geometry_is_a_verified_route_not_a_threshold_guess():
    observations=[dict(repeat=i,route='mmap',allocate_count=1,release_count=-1,allocate_bytes=8192,release_bytes=-8192) for i in range(4)]
    geometry=dict(numpy_version='2.2.6',libc=['glibc','2.35'],affinity=[12,13],numpy_madvise_hugepage=False,
        allocator_handler='default_allocator',allocator_environment={},array_bytes={'4096':dict(route='mmap',mapped_bytes=8192,observations=observations)})
    kwargs=dict(numpy_version='2.2.6',libc=['glibc','2.35'],cpu_affinity=[12,13],allocator_environment={})
    assert validate_allocator_geometry(geometry,**kwargs) is geometry
    geometry['array_bytes']['4096']['route']='arena'
    with pytest.raises(ValueError,match='summary'):validate_allocator_geometry(geometry,**kwargs)
