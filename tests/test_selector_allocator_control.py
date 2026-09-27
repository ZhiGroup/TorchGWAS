"""Safety/ownership checks for the separately built experiment-only allocator."""
import gc
import numpy as np
import pytest
meter=pytest.importorskip('selector_allocator_control',reason='Experiment-only allocator must be built explicitly')


@pytest.fixture(scope='module',autouse=True)
def arena():
    meter.configure(64<<20)


@pytest.mark.parametrize('mode',[0,1,2,3])
@pytest.mark.parametrize('size',[0,4097,1<<20])
def test_values_and_handler_lifetime(mode,size):
    expected=np.arange(size,dtype=np.float32)
    meter.begin(mode)
    try:
        result=np.arange(size,dtype=np.float32)
        view=result[:]
    finally:meter.restore()
    assert np._core.multiarray.get_handler_name()=='default_allocator'
    np.testing.assert_array_equal(result,expected)
    if mode:
        assert meter.snapshot()['live']==1
        with pytest.raises(ValueError,match='live arrays'):meter.begin(mode)
    del result
    if mode:assert meter.snapshot()['live']==1
    del view
    gc.collect()
    snapshot=meter.snapshot()
    assert snapshot['live']==snapshot['invalid_free']==snapshot['foreign_thread']==0
    if mode:assert snapshot['malloc_calls']+snapshot['calloc_calls']==snapshot['free_calls']==1


def test_calloc_and_independent_allocations():
    meter.begin(2)
    try:
        first=np.zeros(1309,dtype=np.float64)
        second=np.empty(1309,dtype=np.float64);second.fill(91.)
    finally:meter.restore()
    assert not np.shares_memory(first,second)
    assert np.count_nonzero(first)==0 and np.all(second==91.)
    del first,second
    assert meter.snapshot()['live']==0


def test_arena_exhaustion_restores_without_reusing_live_storage():
    meter.begin(2)
    try:
        first=np.arange(13,dtype=np.int64)
        with pytest.raises(MemoryError):np.empty(128<<20,dtype=np.uint8)
    finally:meter.restore()
    np.testing.assert_array_equal(first,np.arange(13,dtype=np.int64))
    assert meter.snapshot()['failures']==1
    del first
    assert meter.snapshot()['live']==0


def test_realloc_is_rejected_and_old_array_survives():
    meter.begin(2)
    try:
        value=np.arange(13,dtype=np.int64)
        with pytest.raises(MemoryError):value.resize(100000,refcheck=False)
    finally:meter.restore()
    np.testing.assert_array_equal(value,np.arange(13,dtype=np.int64))
    assert meter.snapshot()['realloc_calls']==1
    del value
    assert meter.snapshot()['live']==0


def test_nested_scope_is_rejected():
    meter.begin(2)
    try:
        with pytest.raises(ValueError):meter.begin(2)
    finally:meter.restore()
    with pytest.raises(ValueError):meter.restore()


def test_default_prefault_has_independent_cpu_and_fault_accounting():
    meter.begin(3)
    try:
        value=np.empty(64<<20,dtype=np.uint8)
        stats=meter.snapshot()
    finally:meter.restore()
    assert stats['arena_used_bytes']==0
    assert stats['prefault_cpu_ns']>0 and stats['allocate_cpu_ns']>=stats['prefault_cpu_ns']
    assert stats['prefault_major_faults']==0
    # Only page starts were initialized; the rest remains np.empty storage.
    assert value[0]==0
    del value
    assert meter.snapshot()['live']==0
