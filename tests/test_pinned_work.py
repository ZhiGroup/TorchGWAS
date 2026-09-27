import pytest
from torchgwas.pinned_work import pinned_scan_work

def test_each_distinct_allocation_is_rounded():
    r=pinned_scan_work(10000,1000,129,3)
    assert r['requested_bytes']==3*(10000000+2*516000+1000+4000)
    assert r['allocator_bytes']==3*((1<<24)+2*(1<<19)+1024+4096)
    assert r['allocation_count']==15
    assert r['allocation_pages']==3*(4096+128+128+1+1)

def test_small_allocations_do_not_share_a_page_in_service_accounting():
    r=pinned_scan_work(32,8,3,2)
    assert r['allocation_pages']==10
    assert r['allocator_bytes']>=r['requested_bytes']

@pytest.mark.parametrize('bad',[0,-1,True,1.5])
def test_invalid_dimensions(bad):
    with pytest.raises(ValueError):pinned_scan_work(bad,8,3,2)