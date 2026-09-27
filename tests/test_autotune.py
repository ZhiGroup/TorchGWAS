"""Joint search invariants, shared bottlenecks, tails and bounded delivery."""
from dataclasses import replace
import threading
import time
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pytest

from torchgwas.autotune import joint_plan
from torchgwas.pipeline_model import Workload, InputProfile, Hardware
from torchgwas.linear import linear_scan_multigpu, multigpu_variant_ranges


def profile(**changes):
    w = Workload(10000, 1000, 128, 2)
    p = InputProfile('test', 2500000, 250, 4, cpu_decode_core_seconds_per_variant=1e-7,
                     host_staging_copies=0, direct_native_fill=True)
    h = Hardware(1e9, 1e10, 1e11, 1e10, 1e11, 2**32, 2**30, 8,
                 d2h_bytes_per_second=1e10)
    options = dict(devices=('cuda:0', 'cuda:1'), chunks=(1000,), workers=(2,), depths=(2,),
                   tie_fraction=0, **changes)
    return w, p, h, options


def test_shared_storage_does_not_double_with_gpus():
    w, p, h, options = profile()
    h = replace(h, disk_bytes_per_second=1e4, gpu_flops_per_second=1e15)
    one = joint_plan(w, p, h, **options, device_sets=[['cuda:0']])['selected']
    two = joint_plan(w, p, h, **options, device_sets=[['cuda:0','cuda:1']])['selected']
    assert one['resource_seconds']['storage'] == two['resource_seconds']['storage']
    assert two['resource_floor_seconds'] == one['resource_floor_seconds']


def test_compute_bound_selects_two_and_produces_executable_worker_budget():
    w, p, h, options = profile()
    result = joint_plan(w, p, h, **options)
    assert len(result['selected']['devices']) == 2
    assert result['selected']['scan_kwargs']['reader_workers'] == 4
    assert result['selected']['scan_kwargs']['ordered'] is False
    assert result['predicted_runtime_seconds'] is None


def test_trait_tail_and_repeated_reads():
    w, p, h, options = profile(reduction='significant', trait_blocks=(50,))
    result = joint_plan(w, p, h, **options, device_sets=[['cuda:0','cuda:1']])['selected']
    assert [job['traits'] for job in result['assignments']] == [50,50,28]
    assert result['genotype_passes'] == 3
    assert result['resource_seconds']['storage'] == pytest.approx(3*p.stored_bytes/h.disk_bytes_per_second)
    assert result['scan_kwargs']['trait_block'] == 50


def test_cpu_and_host_budgets_are_shared():
    w, p, h, options = profile()
    result = joint_plan(w,p,replace(h,cpu_workers=3),**options)
    assert result['rejected']['shared_cpu_budget'] > 0
    single = joint_plan(w,p,h,**options,device_sets=[['cuda:0']])['selected']['host_bytes']
    result = joint_plan(w,p,replace(h,host_memory_bytes=single+1),**options)
    assert result['rejected']['host_memory'] > 0
    assert len(result['selected']['devices']) == 1


def test_common_pcie_uplink_caps_aggregate():
    w,p,h,options = profile()
    result = joint_plan(w,p,h,**options,device_sets=[['cuda:0','cuda:1']],
        shared_links=[dict(devices=['cuda:0','cuda:1'],h2d_bytes_per_second=1000,d2h_bytes_per_second=1000)])['selected']
    assert result['resource_seconds']['link0:h2d'] == (w.variants*p.transfer_bytes_per_variant + 2*4*w.samples*(w.traits+w.covariates+1))/1000


def test_bounds_and_unsupported_modes():
    w,p,h,options = profile()
    with pytest.raises(ValueError,match='max_candidates'):
        joint_plan(w,p,h,**options,max_candidates=1)
    with pytest.raises(ValueError,match='JAGWAS'):
        joint_plan(w,p,h,**options,reduction='jagwas')
    with pytest.raises(ValueError,match='no feasible'):
        joint_plan(w,p,h,**options,device_reserve_bytes=h.device_memory_bytes)


def test_ranges_cover_tail_without_empty_shards():
    assert multigpu_variant_ranges(11,3,3) == [(0,6),(6,9),(9,11)]
    assert multigpu_variant_ranges(3,8,4) == [(0,3)]


def fake_preprocess(y,c,**kwargs):
    return y,c,np.full(y.shape[1],y.shape[0])


def fixture_scan(fake, ordered=False, **kwargs):
    source=SimpleNamespace(shape=(10,20))
    with patch('torchgwas.linear.residualize_and_standardize',side_effect=fake_preprocess), \
         patch('torchgwas.linear.linear_scan_streaming_chunks',side_effect=fake):
        iterator,q=linear_scan_multigpu(source,np.ones((10,2)),np.ones((10,1)),
            devices=['cuda:0','cuda:1'],chunk_size=1,ordered=ordered,result_queue_depth=1,**kwargs)
        yield iterator,q


def test_completion_order_releases_later_shard_while_first_waits():
    later_finished=threading.Event()
    def fake(source,y,c,variant_range,**kwargs):
        def generate():
            start,end=variant_range
            if start==0:
                assert later_finished.wait(2), 'later shard blocked behind first shard'
            for i in range(start,end):
                yield (i,i+1,np.array([i]),None,None)
            if start:
                later_finished.set()
        return generate(),c
    for iterator,q in fixture_scan(fake):
        rows=list(iterator)
    assert sorted(row[0] for row in rows)==list(range(20))
    assert rows[0][0]==10
    assert q is not None


def test_ordered_contract_and_close_cancel_workers():
    closed=[]
    def fake(source,y,c,variant_range,**kwargs):
        def generate():
            try:
                for i in range(*variant_range):
                    yield (i,i+1,None,None,None)
            finally:
                closed.append(variant_range)
        return generate(),c
    for iterator,q in fixture_scan(fake,ordered=True):
        assert [row[0] for row in iterator]==list(range(20))
    closed.clear()
    for iterator,q in fixture_scan(fake):
        next(iterator)
        iterator.close()
    assert len(closed)==2
    assert not any(t.name.startswith('torchgwas-shard-') for t in threading.enumerate())


def test_worker_failure_cancels_full_queue_peers():
    def fake(source,y,c,variant_range,**kwargs):
        def generate():
            if variant_range[0]:
                raise RuntimeError('decode failed')
            for i in range(*variant_range):
                yield (i,i+1,None,None,None)
        return generate(),c
    for iterator,q in fixture_scan(fake):
        with pytest.raises(RuntimeError,match='decode failed'):
            list(iterator)
    assert not any(t.name.startswith('torchgwas-shard-') for t in threading.enumerate())

@pytest.mark.parametrize('missing', [False, True])
def test_real_two_gpu_bed_matches_single_with_partial_chunks(tmp_path, missing):
    import torch
    from test_statistics import _write_bed
    from torchgwas.bed import PlinkBedGenotype
    from torchgwas.linear import linear_scan_streaming_chunks
    if not torch.cuda.is_available() or torch.cuda.device_count() < 2:
        pytest.skip('two CUDA devices required')
    rng=np.random.default_rng(739)
    n,m,k=257,73,5
    genotype=rng.integers(0,3,size=(n,m)).astype(float)
    phenotype=rng.normal(size=(n,k))
    covariates=rng.normal(size=(n,2))
    if missing:
        phenotype[:11,0]=np.nan
        phenotype[::9,2]=np.nan
    bed=_write_bed(tmp_path/'joint',genotype)
    source=PlinkBedGenotype(bed,reader_workers=2,prefetch_chunks=2)
    def collect(iterator):
        b=np.full((m,k),np.nan);t=b.copy();p=b.copy();coverage=np.zeros(m,int)
        for start,end,beta,stat,prob in iterator:
            b[start:end]=beta;t[start:end]=stat;p[start:end]=prob;coverage[start:end]+=1
        assert (coverage==1).all()
        return b,t,p
    one,_=linear_scan_streaming_chunks(source,phenotype,covariates,chunk_size=8,device='cuda:0')
    expected=collect(one)
    for ordered in (False,True):
        many,_=linear_scan_multigpu(source,phenotype,covariates,chunk_size=8,
            devices=['cuda:0','cuda:1'],reader_workers=4,prefetch_chunks=2,ordered=ordered)
        for a,b in zip(collect(many),expected):
            np.testing.assert_allclose(a,b,rtol=5e-5,atol=3e-6,equal_nan=True)


def test_significance_workers_cancel_when_consumer_closes():
    from torchgwas.api import _trait_blocked_significant_chunks
    closed=[]
    def scan_once(offset,width,device):
        try:
            for i in range(100):
                yield (i,i+1,np.array([i]),np.array([0]),np.array([1.]),np.array([2.]),np.array([8.]))
        finally:
            closed.append(device)
    with patch('torchgwas.linear._significant_pairs_iterator', side_effect=lambda chunks,*args: chunks):
        iterator=_trait_blocked_significant_chunks(scan_once,None,4,2,8,['cuda:0','cuda:1'])
        next(iterator)
        iterator.close()
    assert closed
    assert not any(t.name.startswith('torchgwas-sigshard-') for t in threading.enumerate())


def test_direct_native_rounded_allocations_can_reject_host_budget():
    w,p,h,options=profile()
    options['device_sets']=[['cuda:0']]
    selected=joint_plan(w,p,h,**options)['selected']
    from torchgwas.pinned_work import pinned_scan_work
    pins=pinned_scan_work(w.samples,1000,w.traits,2,transfer_bytes_per_variant=250)
    extra=pins['allocator_bytes']-pins['requested_bytes']
    assert extra>0
    # A budget above the old logical accounting but below the actual rounded
    # request accounting must reject this sole candidate.
    budget=int(selected['host_bytes']-extra/2)
    with pytest.raises(ValueError,match='no feasible'):
        joint_plan(w,p,replace(h,host_memory_bytes=budget),**options)


def test_preprocessing_peak_can_reject_an_otherwise_small_scan_ring():
    w,p,h,options=profile()
    w=replace(w,traits=4096,covariates=8)
    options.update(chunks=(8,),device_sets=[['cuda:0']],modes=('variant',))
    from torchgwas.setup_work import setup_memory
    bound=setup_memory(w.samples,w.traits,w.covariates)['device_live_bytes_upper']
    with pytest.raises(ValueError,match='no feasible'):
        joint_plan(w,p,replace(h,device_memory_bytes=bound-1),**options)
    chosen=joint_plan(w,p,replace(h,device_memory_bytes=bound+1),**options)['selected']
    assert chosen['device_bytes']['cuda:0']>=bound


def test_multigpu_propagates_reader_budget_to_native_limit():
    budgets=[]
    def fake(source,y,c,variant_range,reader_workers,_reader_worker_limit,**kwargs):
        budgets.append((reader_workers,_reader_worker_limit))
        return iter(()),None
    for iterator,q in fixture_scan(fake,reader_workers=4):list(iterator)
    assert budgets==[(2,2),(2,2)]


def test_native_worker_limit_preserves_single_scan_preference():
    from torchgwas.native_scan import resolve_reader_workers
    source=SimpleNamespace(decode_workers=8)
    assert resolve_reader_workers(source,4)==8
    assert resolve_reader_workers(source,2,2)==2
    assert source.decode_workers==8
    assert resolve_reader_workers(SimpleNamespace(),None,3)==3
    assert resolve_reader_workers(SimpleNamespace(decode_workers=1),2,2)==1
    for limit in [0,-1,True,1.5]:
        with pytest.raises(ValueError,match='positive integer'):
            resolve_reader_workers(source,2,limit)


def test_multigpu_preprocess_uses_first_selected_device_once():
    import torch
    source=SimpleNamespace(shape=(10,20))
    with patch('torchgwas.linear.choose_device',return_value=torch.device('cuda:2')), \
         patch('torchgwas.linear.residualize_and_standardize',side_effect=fake_preprocess) as prep, \
         patch('torchgwas.linear.linear_scan_streaming_chunks',return_value=(iter(()),None)):
        iterator,_=linear_scan_multigpu(source,np.ones((10,2)),np.ones((10,1)),
            devices=['cuda:2','cuda:3'],chunk_size=1,reader_workers=4,ordered=False)
        list(iterator)
        prep.assert_called_once()
        assert prep.call_args.kwargs['device']==torch.device('cuda:2')


def test_multigpu_reader_remainder_is_not_lost():
    budgets=[]
    def fake(source,y,c,variant_range,reader_workers,_reader_worker_limit,**kwargs):
        budgets.append((variant_range[0],reader_workers,_reader_worker_limit))
        return iter(()),None
    for iterator,q in fixture_scan(fake,reader_workers=5):list(iterator)
    assert sorted(budgets)==[(0,3,3),(10,2,2)]


def test_multigpu_collapsed_shards_use_actual_device_count_and_cap():
    source=SimpleNamespace(shape=(10,20))
    with patch('torchgwas.linear.linear_scan_streaming_chunks',return_value=(iter(()),None)) as scan:
        iterator,_=linear_scan_multigpu(source,np.ones((10,2)),None,
            devices=['cuda:0','cuda:1','cuda:2'],chunk_size=64,reader_workers=1)
        list(iterator)
        scan.assert_called_once()
        assert scan.call_args.kwargs['device']=='cuda:0'
        assert scan.call_args.kwargs['reader_workers']==1
        assert scan.call_args.kwargs['_reader_worker_limit']==1


@pytest.mark.parametrize('workers',[0,-1,True,2.5])
def test_multigpu_rejects_noninteger_reader_budgets(workers):
    with pytest.raises(ValueError,match='positive integer'):
        linear_scan_multigpu(SimpleNamespace(shape=(10,20)),np.ones((10,2)),None,
            devices=['cuda:0','cuda:1'],reader_workers=workers)
