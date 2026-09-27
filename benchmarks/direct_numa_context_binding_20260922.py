"""Actual context-binding checks under private thread memory-policy changes.

The rate is a synthetic test control. This collects no performance coefficient,
does not alter kernel settings, and restores the calling thread's original policy.
"""
import ctypes
import json
from pathlib import Path
import time

from torchgwas.detailed_calibration import (execution_context,source_identity,sha256_file,
    bind_detailed_profile,write_detailed_profile,validate_detailed_profile)
from torchgwas.numa_context import _library,memory_policy_context


def main():
    root=Path('results/numa_context_binding_20260922/audit');root.mkdir(parents=True,exist_ok=False)
    sources=source_identity();harness=sha256_file(__file__)
    library=_library();possible=library.numa_num_possible_nodes();bits=8*ctypes.sizeof(ctypes.c_ulong)
    setter=library.set_mempolicy
    setter.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_ulong),ctypes.c_ulong]
    setter.restype=ctypes.c_long
    def set_policy(mode,nodes):
        mask=(ctypes.c_ulong*((possible+bits-1)//bits))()
        for node in nodes:mask[node//bits]|=1<<(node%bits)
        if setter(mode,mask if nodes else None,possible if nodes else 0)!=0:
            raise OSError(ctypes.get_errno(),'Private thread memory-policy change failed')
    def capture():
        return execution_context(['cuda:1'],input_path='/data/zxie3/torchgwas_public_refresh_v3_20260922/input.pgen',
            output_path=root)
    original=memory_policy_context();before=capture();node=original['allowed_nodes'][0]
    artifact=root/'synthetic.json';artifact.write_text(json.dumps(dict(value=1.,scope='Synthetic control; never a hardware rate.')))
    contexts=[dict(name='one',devices=['cuda:1'],profiles={'cuda:1':dict(primitive=1.)},shared_capacities=dict(cpu=1.))]
    profile=bind_detailed_profile(contexts,before,component_artifacts={str(artifact.resolve()):sha256_file(artifact)},
        limitations=['Context validation only; all numeric prices are synthetic controls.'],sources=sources)
    path=root/'original_profile.json';write_detailed_profile(profile,path);saved=path.read_bytes()
    validate_detailed_profile(profile,before,sources=sources)
    rows=[]
    try:
        for mode in [2,8194]:
            set_policy(mode,[node]);current=capture()
            assert current['numa_policy']['thread_policy_mode']==mode
            assert current['numa_policy']['thread_policy_nodes']==[node]
            try:validate_detailed_profile(profile,current,sources=sources)
            except ValueError as error:
                assert 'numa_policy' in str(error);rejection=str(error)
            else:raise AssertionError('A changed NUMA policy reused the old price context')
            assert path.read_bytes()==saved
            rows.append(dict(requested_mode=mode,execution_context=current,rejection=rejection))
    finally:set_policy(original['thread_policy_mode'],original['thread_policy_nodes'])
    after=capture();assert after==before
    validate_detailed_profile(profile,after,sources=sources)
    assert source_identity()==sources and sha256_file(__file__)==harness and path.read_bytes()==saved
    report=dict(source_sha256=sources,harness_sha256=harness,original_context=before,rows=rows,
        original_profile_sha256=sha256_file(path),original_profile_bytes_unchanged=True,
        original_policy_restored=True,old_profile_accepted_after_restore=True,
        observed_unix_seconds=time.time(),prices_measured=False,scope=__doc__)
    with (root/'report.json').open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps(dict(rejected_modes=[r['requested_mode'] for r in rows],original_policy_restored=True,
        original_profile_bytes_unchanged=True,source_files=len(sources))),flush=True)


if __name__=='__main__':main()
