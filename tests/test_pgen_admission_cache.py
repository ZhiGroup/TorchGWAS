"""Source-bound structural admission reuse without renewing performance prices."""
import os
from unittest.mock import patch

import pytest

from torchgwas.adaptive_start import prepare_adaptive_start
from test_adaptive_start import options
from test_trait_candidate_space import spec as dense_spec


@pytest.mark.parametrize('reduction',[None,'significant'])
def test_cached_compact_admission_matches_uncached_and_invalidates(tmp_path,reduction):
    spec=dense_spec(tmp_path);context=spec['contexts'][1]
    output=dict(spec['output'],block_bytes=None) if reduction else spec['output']
    kw=options(context,reduction=reduction)
    if reduction:kw['significance_threshold']=.01
    baseline=prepare_adaptive_start(spec['workload'],context,output=output,**kw,
                                    _compact_memory=True)
    caches=[];indexes=[]
    with patch('torchgwas.pgen_admission_cache.input_is_stable',return_value=True):
        first=prepare_adaptive_start(spec['workload'],context,output=output,**kw,
            _compact_memory=True,_compact_cache_dir=tmp_path,
            _compact_cache_receiver=caches.append,_index_receiver=indexes.append)
        assert len(caches)==1 and caches[0].status=='miss'
        assert caches[0].pending_bytes>0
        assert not caches[0].path.exists()
        assert caches[0].publish(successful=True)=='stored'
        assert caches[0].path.exists()
        with patch('torchgwas.pgen_memory_layout.memory_layout',side_effect=AssertionError('structural vectors rebuilt')):
            second=prepare_adaptive_start(spec['workload'],context,output=output,**kw,
                _compact_memory=True,_compact_cache_dir=tmp_path,
                _compact_cache_receiver=caches.append,_index_receiver=indexes.append)
        assert caches[-1].status=='hit' and caches[-1].pending_bytes==0
        assert caches[-1].retained_blob_bytes==caches[-1].path.stat().st_size
        assert caches[-1].owned_extra_bytes(indexes[-1][2].nbytes)>0
        assert second['memory']==baseline['memory']==first['memory']
        assert second['partitions']==baseline['partitions']
        assert second['api_kwargs']==baseline['api_kwargs']
        assert indexes[0][2].tolist()==indexes[1][2].tolist()
        assert not indexes[1][2].flags.writeable
        assert indexes[1][2].flags.aligned
        caches[-1].path.write_bytes(b'corrupt')
        again=prepare_adaptive_start(spec['workload'],context,output=output,**kw,
            _compact_memory=True,_compact_cache_dir=tmp_path,
            _compact_cache_receiver=caches.append,_index_receiver=indexes.append)
        assert caches[-1].status=='invalid'
        assert again['memory']==baseline['memory']
        assert caches[-1].publish(successful=False)=='unsuccessful'
        stat=os.stat(spec['workload']['genotype'])
        os.utime(spec['workload']['genotype'],ns=(stat.st_atime_ns,stat.st_mtime_ns+2_000_000_000))
        changed=prepare_adaptive_start(spec['workload'],context,output=output,**kw,
            _compact_memory=True,_compact_cache_dir=tmp_path,
            _compact_cache_receiver=caches.append,_index_receiver=indexes.append)
        assert caches[-1].status=='miss'
        assert changed['memory']==baseline['memory']


def test_unstable_source_never_loads_or_publishes(tmp_path):
    spec=dense_spec(tmp_path);context=spec['contexts'][1];caches=[]
    with patch('torchgwas.pgen_admission_cache.input_is_stable',return_value=False):
        prepare_adaptive_start(spec['workload'],context,output=spec['output'],
            **options(context),_compact_memory=True,_compact_cache_dir=tmp_path,
            _compact_cache_receiver=caches.append,_index_receiver=lambda _:None)
    assert len(caches)==1 and caches[0].status=='unstable_source'
    assert caches[0].pending_bytes==0
    assert caches[0].publish(successful=True)=='unstable_source'
    assert not caches[0].path.exists()


def test_cold_deferred_admission_publishes_bases_after_success(tmp_path):
    spec=dense_spec(tmp_path);context=spec['contexts'][1]
    kw=options(context)
    baseline=prepare_adaptive_start(spec['workload'],context,output=spec['output'],
                                    **kw,_compact_memory=True)
    caches=[];indexes=[]
    with patch('torchgwas.pgen_admission_cache.input_is_stable',return_value=True):
        cold=prepare_adaptive_start(spec['workload'],context,output=spec['output'],
            **kw,_compact_memory=True,_defer_compact_bases=True,
            _compact_cache_dir=tmp_path,_compact_cache_receiver=caches.append,
            _index_receiver=indexes.append)
        assert not indexes
        assert caches[-1].status=='miss' and caches[-1].pending_bytes>0
        assert cold['memory']==baseline['memory']
        assert caches[-1].publish(successful=True)=='stored'
        assert caches[-1].publication_cpu_seconds>=0
        with patch('torchgwas.pgen_memory_layout.memory_layout',
                   side_effect=AssertionError('cold memory layout rebuilt')):
            hot=prepare_adaptive_start(spec['workload'],context,output=spec['output'],
                **kw,_compact_memory=True,_defer_compact_bases=True,
                _compact_cache_dir=tmp_path,_compact_cache_receiver=caches.append,
                _index_receiver=indexes.append)
        assert caches[-1].status=='hit' and len(indexes)==1
        assert hot['memory']==baseline['memory']
        assert not indexes[0][2].flags.writeable
        assert indexes[0][2].flags.aligned
