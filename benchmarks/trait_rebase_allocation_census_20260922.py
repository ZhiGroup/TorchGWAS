"""Allocation counts at NumPy's temporary-elision boundary; no timing output."""
import hashlib
import json
from pathlib import Path
import numpy as np
import selector_allocator_control as meter
from torchgwas.api import _trait_blocked_significant_chunks
from torchgwas.detailed_calibration import source_identity,_numpy_core_context
from legacy_trait_rebase_20260922 import legacy_trait_blocked_significant_chunks,ORIGINAL_FUNCTION_SHA256


def main():
    root=Path('results/trait_rebase_allocation_census_20260922');root.mkdir(parents=True,exist_ok=False)
    source=source_identity();core=_numpy_core_context();records=[]
    for count in [0,7,1024,32767,32768,32769,1<<20,1024*8193]:
        incoming=np.zeros(count,np.int64);incoming.flags.writeable=False
        values=np.zeros(count,np.float32);empty=np.empty(0,np.int64);empty_values=np.empty(0,np.float32)
        def scan(offset,width,device):
            yield (0,1,empty,empty,None,empty_values,empty_values) if offset==0 else (0,1,incoming,incoming,None,values,values)
        functions={'legacy':legacy_trait_blocked_significant_chunks,'current':_trait_blocked_significant_chunks}
        for repeat in range(4):
            for name in (list(functions) if repeat%2==0 else list(reversed(functions))):
                meter.begin(1)
                try:result=list(functions[name](scan,None,8194,4097,17))
                finally:meter.restore()
                before=meter.snapshot()
                assert len(result)==2 and result[1][3].shape==incoming.shape
                assert not np.shares_memory(result[1][3],incoming)
                assert np.all(result[1][3]==4097) and np.all(incoming==0)
                assert before['live']==2
                del result
                after=meter.snapshot();assert after['live']==0
                for stats in [before,after]:
                    assert all(stats[key]==0 for key in ['failures','invalid_free','foreign_thread','realloc_calls'])
                fields=['malloc_calls','calloc_calls','free_calls','allocate_bytes','free_bytes','live']
                records.append(dict(count=count,repeat=repeat,implementation=name,
                    before_release={key:before[key] for key in fields},after_release={key:after[key] for key in fields}))
    summary=[]
    for count in sorted({r['count'] for r in records}):
        allocated={}
        for name in ['legacy','current']:
            values={r['before_release']['allocate_bytes'] for r in records if r['count']==count and r['implementation']==name}
            assert len(values)==1,(count,name,values);allocated[name]=values.pop()
        summary.append(dict(count=count,allocated_bytes=allocated,saved_allocation_bytes=allocated['legacy']-allocated['current']))
    assert source==source_identity() and core==_numpy_core_context()
    report=dict(records=records,summary=summary,duration_fields_retained=False,source_sha256=source,numpy_core=core,
        legacy_function_sha256=ORIGINAL_FUNCTION_SHA256,
        harness_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),Path(meter.__file__),Path(__file__).with_name('legacy_trait_rebase_20260922.py')]},
        scope='Four counterbalanced forwarding-allocator observations per count. Callback timings are not retained or used. Source boundary is NPY_MIN_ELIDE_BYTES=256KiB in NumPy 2.2.6 non-debug builds; observed allocator evidence is runtime-specific. Explicit new rebase does not depend on elision.')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    for row in summary:print(json.dumps(row),flush=True)


if __name__=='__main__':main()
