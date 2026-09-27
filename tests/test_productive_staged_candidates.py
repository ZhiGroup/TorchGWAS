"""First-chunk candidates use active prices and the exact unissued suffix."""
from copy import deepcopy

import pytest

from test_productive_source_floor import setup
from test_productive_staged_work_binding import case
from torchgwas.productive_staged_candidates import fixed_partition_chunk_candidates
from torchgwas.productive_staged_source_binding import audit_staged_source_price_binding
from torchgwas.productive_staged_work_binding import audit_staged_work_price_binding


def inputs(tmp_path,mode):
    _,source,frontier=setup(tmp_path,reduction=mode)
    profile,_=case()
    context=profile['contexts'][0]
    active=context['profiles']['cuda:0']
    active.update(source)
    active['host_primitives']=dict(synthetic_dispatch_seconds=.001)
    active['kernel_geometry']=[dict(N=129,B=b,K=5,C=3,kernels=[])
                               for b in (2,4)]
    context['shared_capacities'].update(cpu=1.,dram=5e7,input=5e6)
    output=dict(block_bytes=32,queue_depth=2,store_beta=True,fsync=True)
    writer=dict(block_bytes=32,queue_depth=2,borrow_chunks=True,fsync=True,
                writeback_bytes=0,sync_file_range=False,store_variant_df=False)
    return profile,context,frontier,output,writer


@pytest.mark.parametrize('mode', [None,'significant','jagwas'])
def test_chunk_candidates_preserve_frontier_and_bind_active_prices(tmp_path,mode):
    profile,context,frontier,output,writer=inputs(tmp_path,mode)
    if mode=='jagwas':
        context['profiles']['cuda:0']['gpu_resources']['fp64_flops_per_second']=10000.
        prices=dict(writer_prices=dict(prices={'primitive':.01},
                                       archive={'call':.02}))
    elif mode=='significant':
        prices=dict(prices={'primitive':.01},archive={'call':.02})
    else:
        prices=None
    rows=fixed_partition_chunk_candidates(frontier,context,
        chunk_sizes=[4,2],partition_axis=('variant' if mode=='jagwas' else 'trait'),
        covariate_rank=3,output=output,
        writer_options=(writer if mode is None else None),
        reduction_prices=prices,
        significance_threshold=(.05 if mode=='significant' else None))
    assert [row['chunk_markers'] for row in rows]==[4,2]
    assert all(row['partitions']==frontier['rectangles'] for row in rows)
    assert rows[0]['source_profiles']['original']['decode_units']==(
        context['profiles']['cuda:0']['decode_units'])
    assert rows[0]['compute_options']['shared_h2d_bytes_per_second']==1000.
    assert rows[0]['output_options']['shared_d2h_bytes_per_second']==500.
    assert rows[0]['shape_profiles']['cuda:0']['chunk_markers']==4
    assert rows[1]['shape_profiles']['cuda:0']['chunk_markers']==2
    assert audit_staged_source_price_binding(profile,'current',rows,
        dict(cpu=1.,dram=5e7,input=5e6))['candidate_count']==2
    evidence=(None if prices is None else dict(record={'value':prices},
        record_sha256='a'*64,artifact_sha256='b'*64))
    assert audit_staged_work_price_binding(profile,'current',rows,
        reduction_price_evidence=evidence)['candidate_count']==2
    if mode=='jagwas':
        assert all(row['partitions'][0]['trait_range']==[0,5] for row in rows)
        assert all(row['partition_axis']=='variant' for row in rows)
    if mode=='significant':
        assert rows[0]['output_options']['significant_backend']=='host'


def test_candidate_factory_requires_shared_transfer_measurement(tmp_path):
    _,context,frontier,output,writer=inputs(tmp_path,None)
    del context['shared_transfer_capacities']
    with pytest.raises(ValueError,match='transfer capacities'):
        fixed_partition_chunk_candidates(frontier,context,chunk_sizes=[4],
            partition_axis='trait',covariate_rank=3,output=output,
            writer_options=writer)


def test_candidate_factory_rejects_a_different_writer_and_duplicate_chunk(tmp_path):
    _,context,frontier,output,writer=inputs(tmp_path,None)
    mismatch=deepcopy(writer);mismatch['block_bytes']=64
    with pytest.raises(ValueError,match='Actual dense writer'):
        fixed_partition_chunk_candidates(frontier,context,chunk_sizes=[4],
            partition_axis='trait',covariate_rank=3,output=output,
            writer_options=mismatch)
    with pytest.raises(ValueError,match='Bounded active'):
        fixed_partition_chunk_candidates(frontier,context,chunk_sizes=[4,4],
            partition_axis='trait',covariate_rank=3,output=output,
            writer_options=writer)


def test_constructed_dense_candidates_run_the_real_post_output_screen(
        tmp_path,monkeypatch):
    from test_layout_gpu_shape_service import _component
    from test_productive_staged_screen import staged
    from torchgwas.productive_staged_screen import productive_staged_partial_screen
    header,_,frontier=setup(tmp_path)
    second=tmp_path/'second';second.mkdir()
    _,context,_,output,writer=inputs(second,None)
    active=context['profiles']['cuda:0']
    active['reduction']=None
    active['gpu_resources'].update(hbm_bytes_per_second=1e12,
        l2_bytes_per_second=2e12,host_dispatch_cpu_seconds=2e-6,
        available_l2_bytes=1<<20,sm_count=100)
    active['kernel_geometry']=[dict(N=129,B=b,K=5,C=3,kernels=[])
                               for b in (2,4)]
    rows=fixed_partition_chunk_candidates(frontier,context,
        chunk_sizes=[2,4],partition_axis='trait',covariate_rank=3,
        output=output,writer_options=writer)
    with monkeypatch.context() as patch:
        patch.setattr('torchgwas.mechanistic_torch._shape_component',_component)
        report=productive_staged_partial_screen(frontier,staged(header),rows,
            shared_source_capacities=dict(cpu=1.,dram=5e7,input=5e6))
    assert report['evaluated_candidates']==2
    assert all(row['partial']['envelope']['coverage']['exact_coverage']
               for row in report['candidates'])
    assert all(row['partial']['gpu_shape_service']['distinct_shapes']>=1
               for row in report['candidates'])
