"""A memory-policy change must not reuse an old independent price context."""
from copy import deepcopy
import ctypes
import errno
import json
import subprocess
import sys

import pytest

from torchgwas import numa_context as numa
from torchgwas import detailed_calibration as calibration
from test_detailed_calibration import fixture as profile_fixture


class Policy:
    def __init__(self):self.mode=0;self.nodes=[];self.allowed=[0,1,130];self.failure=0;self.calls=0
    def numa_num_possible_nodes(self):return 1024
    def get_mempolicy(self,mode,mask,possible,address,flags):
        self.calls+=1
        if self.failure:ctypes.set_errno(self.failure);return -1
        if mode is not None:ctypes.cast(mode,ctypes.POINTER(ctypes.c_int))[0]=self.mode
        bits=8*ctypes.sizeof(ctypes.c_ulong)
        for node in self.allowed if flags==4 else self.nodes:mask[node//bits]|=1<<(node%bits)
        return 0


@pytest.fixture
def policy(monkeypatch):
    value=Policy();monkeypatch.setattr(numa,'_library',lambda:value)
    monkeypatch.setattr(numa,'_automatic_balancing',lambda:1)
    return value


def test_policy_reads_are_fresh_and_preserve_high_mask_bits(policy):
    before=numa.memory_policy_context()
    assert before['thread_policy_mode']==0 and before['thread_policy_nodes']==[]
    assert before['allowed_nodes']==[0,1,130]
    policy.mode=2|32768;policy.nodes=[130];policy.allowed=[1,130]
    after=numa.memory_policy_context()
    assert after['thread_policy_mode']==32770 and after['thread_policy_nodes']==[130]
    assert after['allowed_nodes']==[1,130] and policy.calls==6
    assert before['allowed_nodes']==[0,1,130]


@pytest.mark.parametrize('field,value',[
    ('thread_policy_mode',2),('thread_policy_mode',2|8192),('thread_policy_nodes',[1]),
    ('allowed_nodes',[1]),('automatic_balancing',0),('automatic_balancing',None)])
def test_changed_policy_invalidates_saved_profile_without_modifying_it(policy,tmp_path,field,value):
    profile,execution,sources,_=profile_fixture(tmp_path)
    execution['numa_policy']=numa.memory_policy_context();profile['execution_context']=deepcopy(execution)
    saved=json.dumps(profile,sort_keys=True)
    execution['numa_policy'][field]=value
    with pytest.raises(ValueError,match='numa_policy.'+field):
        calibration.validate_detailed_profile(profile,execution,sources=sources)
    assert json.dumps(profile,sort_keys=True)==saved


def test_legacy_profile_is_not_silently_given_a_current_memory_policy(policy,tmp_path):
    profile,execution,sources,_=profile_fixture(tmp_path)
    execution['numa_policy']=numa.memory_policy_context()
    with pytest.raises(ValueError,match='context.numa_policy'):
        calibration.validate_detailed_profile(profile,execution,sources=sources)
    assert 'numa_policy' not in profile['execution_context']


@pytest.mark.parametrize('failure',[errno.EPERM,errno.ENOSYS,errno.EINVAL])
def test_unreadable_policy_fails_closed(policy,failure):
    policy.failure=failure
    with pytest.raises(ValueError,match='get_mempolicy failed'):numa.memory_policy_context()


def test_empty_allowed_mask_fails_closed(policy):
    policy.allowed=[]
    with pytest.raises(ValueError,match='allowed-node mask is empty'):numa.memory_policy_context()


def test_policy_changes_during_capture_are_rejected(policy,monkeypatch):
    original=policy.get_mempolicy
    def changed(*args):
        result=original(*args)
        if policy.calls==2:policy.mode=2;policy.nodes=[1]
        return result
    monkeypatch.setattr(policy,'get_mempolicy',changed)
    with pytest.raises(ValueError,match='changed during'):numa.memory_policy_context()


def test_balancing_changes_during_capture_are_rejected(policy,monkeypatch):
    values=iter([1,0]);monkeypatch.setattr(numa,'_automatic_balancing',lambda:next(values))
    with pytest.raises(ValueError,match='changed during'):numa.memory_policy_context()


def test_real_thread_policy_change_is_private_and_inherited_by_new_worker():
    original=numa.memory_policy_context()
    # Set policy only in this fresh child. No system/cpuset settings are changed.
    command=r'''
import ctypes,json,threading
from torchgwas.numa_context import memory_policy_context,_library
before=memory_policy_context();node=before['allowed_nodes'][0];lib=_library()
possible=lib.numa_num_possible_nodes();bits=8*ctypes.sizeof(ctypes.c_ulong)
mask=(ctypes.c_ulong*((possible+bits-1)//bits))();mask[node//bits]|=1<<(node%bits)
setter=lib.set_mempolicy
setter.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_ulong),ctypes.c_ulong]
setter.restype=ctypes.c_long
assert setter(2,mask,possible)==0,ctypes.get_errno()
after=memory_policy_context();worker=[]
thread=threading.Thread(target=lambda:worker.append(memory_policy_context()))
thread.start();thread.join()
assert after['thread_policy_mode']==2 and after['thread_policy_nodes']==[node]
assert worker==[after]
assert setter(2|8192,mask,possible)==0,ctypes.get_errno()
balanced=memory_policy_context()
assert balanced['thread_policy_mode']==8194 and balanced['thread_policy_nodes']==[node]
print(json.dumps(dict(before=before,after=after,worker=worker,balanced=balanced)))
'''
    result=json.loads(subprocess.check_output([sys.executable,'-c',command],text=True))
    assert result['after']['thread_policy_mode']==2
    assert numa.memory_policy_context()==original
