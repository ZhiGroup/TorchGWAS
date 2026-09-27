"""Consumer-visible ownership, independent of CUDA or result-array size."""
import threading
import weakref
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pytest
from torchgwas.linear import linear_scan_multigpu


@pytest.mark.parametrize('retain_first',[False,True])
def test_discarded_multigpu_results_are_released_before_waiting_for_more(retain_first):
    release_next=threading.Event();freed={i:threading.Event() for i in (0,2)}
    paused={i:threading.Event() for i in (0,2)};references=[]
    def item(index):
        array=np.full(4,index,dtype=np.float32)
        if index in freed:references.append(weakref.ref(array,lambda ref,i=index:freed[i].set()))
        return index,index+1,array,None,None
    def fake_scan(source,y,c,variant_range,**kwargs):
        lo,hi=variant_range
        def generate():
            yield item(lo)
            paused[lo].set()
            if not release_next.wait(5):raise RuntimeError('test producer not released')
            yield item(lo+1)
        return generate(),c
    def preprocess(y,c,**kwargs):return y,c,np.full(y.shape[1],y.shape[0])
    with patch('torchgwas.linear.choose_device',return_value=SimpleNamespace(type='cpu')), \
         patch('torchgwas.linear.residualize_and_standardize',side_effect=preprocess), \
         patch('torchgwas.linear.linear_scan_streaming_chunks',side_effect=fake_scan):
        iterator,_=linear_scan_multigpu(SimpleNamespace(shape=(10,4)),np.ones((10,1)),None,
            devices=['cuda:0','cuda:1'],chunk_size=1,ordered=False,result_queue_depth=1)
        errors=[];thread=None
        try:
            first=next(iterator);first_index=first[0]
            retained=first if retain_first else None
            del first
            second=next(iterator);del second
            assert all(event.wait(5) for event in paused.values())
            def advance():
                try:next(iterator)
                except BaseException as error:errors.append(error)
            thread=threading.Thread(target=advance);thread.start()
            assert all(event.wait(2) for index,event in freed.items() if not retain_first or index!=first_index), 'discarded arrays retained while waiting for next chunk'
            if retain_first:
                assert not freed[first_index].is_set()
                np.testing.assert_array_equal(retained[2],np.full(4,first_index,dtype=np.float32))
                retained=None
                assert freed[first_index].wait(2)
        finally:
            release_next.set()
            if thread is not None:thread.join(5);assert not thread.is_alive()
            iterator.close()
        assert not errors
