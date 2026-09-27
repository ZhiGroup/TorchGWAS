"""Read the calling thread's NUMA allocation policy without changing it.

Thread policy is inherited by newly created workers. It does not describe
pre-existing threads, per-VMA mbind overrides, residency or live capacity.
"""
import ctypes
from functools import lru_cache
import os
from pathlib import Path
import sys


@lru_cache(maxsize=1)
def _library():
    try:
        library=ctypes.CDLL('libnuma.so.1',use_errno=True)
        library.numa_num_possible_nodes.argtypes=[]
        library.numa_num_possible_nodes.restype=ctypes.c_int
        library.get_mempolicy.argtypes=[ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_ulong),
            ctypes.c_ulong,ctypes.c_void_p,ctypes.c_ulong]
        library.get_mempolicy.restype=ctypes.c_long
        return library
    except (OSError,AttributeError) as error:
        raise ValueError('libnuma get_mempolicy is required to bind NUMA context') from error


def _automatic_balancing():
    try:text=Path('/proc/sys/kernel/numa_balancing').read_text().strip()
    except FileNotFoundError:return None
    except OSError as error:raise ValueError('Automatic NUMA-balancing setting unavailable') from error
    try:value=int(text)
    except ValueError as error:raise ValueError('Invalid automatic NUMA-balancing setting') from error
    if value<0:raise ValueError('Invalid automatic NUMA-balancing setting')
    return value


def _policy(library,possible,flags):
    bits=8*ctypes.sizeof(ctypes.c_ulong)
    mask=(ctypes.c_ulong*((possible+bits-1)//bits))()
    mode=ctypes.c_int()
    ctypes.set_errno(0)
    result=library.get_mempolicy(ctypes.byref(mode) if not flags else None,mask,possible,None,flags)
    if result!=0:
        error=ctypes.get_errno()
        raise ValueError('NUMA get_mempolicy failed: '+os.strerror(error))
    nodes=[node for node in range(possible) if mask[node//bits] & (1<<(node%bits))]
    if not flags and mode.value<0:raise ValueError('Invalid NUMA policy mode')
    return dict(mode=mode.value,nodes=nodes) if not flags else nodes


def memory_policy_context():
    """Fresh policy and allowed-node reads; only the library handle is cached."""
    if sys.platform!='linux':raise ValueError('NUMA context currently requires Linux')
    library=_library();possible=library.numa_num_possible_nodes()
    if type(possible) is not int or not 1<=possible<=8*os.sysconf('SC_PAGE_SIZE'):
        raise ValueError('Invalid NUMA possible-node mask size')
    balancing=_automatic_balancing()
    policy=_policy(library,possible,0)
    allowed=_policy(library,possible,4)  # MPOL_F_MEMS_ALLOWED, not a policy change.
    if not allowed:raise ValueError('NUMA allowed-node mask is empty')
    if policy!=_policy(library,possible,0) or balancing!=_automatic_balancing():
        raise ValueError('NUMA policy changed during context capture')
    return dict(thread_policy_mode=policy['mode'],thread_policy_nodes=policy['nodes'],
        allowed_nodes=allowed,automatic_balancing=balancing,
        scope='Calling thread default policy and allowed nodes; excludes other existing threads, VMA policies, residency and available capacity.')
