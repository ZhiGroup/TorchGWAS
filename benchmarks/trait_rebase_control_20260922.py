"""Paired source rebase CPU/allocation controls, independent of GWAS."""
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import time
import numpy as np
import torch
import selector_allocator_control as meter
from torchgwas.api import _trait_blocked_significant_chunks
from torchgwas.detailed_calibration import source_identity,_numpy_core_context
from legacy_trait_rebase_20260922 import legacy_trait_blocked_significant_chunks,ORIGINAL_FUNCTION_SHA256
from selector_resident_allocator_control_20260922 import call


def main():
    root=Path('results/trait_rebase_control_20260922');root.mkdir(parents=True,exist_ok=False)
    os.sched_setaffinity(0,list(range(12,20)));torch.set_num_threads(4)
    assert not np._core.multiarray._get_madvise_hugepage()
    source=source_identity();core=_numpy_core_context();records=[];rng=random.Random(9221521)
    started=datetime.now(timezone.utc).isoformat()
    for rows in [0,1<<20,1024*8193]:
        incoming=np.arange(rows,dtype=np.int64);incoming%=4097;incoming.flags.writeable=False
        variants=np.zeros(rows,np.int64);values=np.ones(rows,np.float32)
        empty=np.empty(0,np.int64);empty_values=np.empty(0,np.float32)
        def scan(offset,width,device):
            if offset==0:yield (0,1,empty,empty,None,empty_values,empty_values)
            else:yield (0,1,variants,incoming,None,values,values)
        expected=incoming+4097
        def validate(result):
            assert len(result)==2
            out=result[1][3];assert out.dtype==np.int64 and out.shape==incoming.shape
            assert not np.shares_memory(out,incoming) and out.flags.writeable
            assert out.tobytes()==expected.tobytes()
        functions={'legacy':legacy_trait_blocked_significant_chunks,'current':_trait_blocked_significant_chunks}
        for function in functions.values():
            for mode in [0,1]:call(lambda:list(function(scan,None,8194,4097,17)),mode,validate)
        for repeat in range(7):
            order=[(name,mode) for name in functions for mode in [0,1]];rng.shuffle(order)
            for name,mode in order:
                observation=call(lambda:list(functions[name](scan,None,8194,4097,17)),mode)
                records.append(dict(rows=rows,repeat=repeat,implementation=name,mode=mode,observation=observation))
        print(json.dumps(dict(rows=rows,medians={name:statistics.median(r['observation']['cpu_seconds'] for r in records
            if r['rows']==rows and r['mode']==0 and r['implementation']==name) for name in functions})),flush=True)
    assert source==source_identity() and core==_numpy_core_context() and meter.snapshot()['live']==0
    report=dict(observation_started_at_utc=started,observation_finished_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=source,numpy_core=core,affinity=sorted(os.sched_getaffinity(0)),
        torch_threads=torch.get_num_threads(),numpy_madvise_hugepage=False,legacy_function_sha256=ORIGINAL_FUNCTION_SHA256,
        records=records,harness_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),Path(meter.__file__),Path(__file__).with_name('legacy_trait_rebase_20260922.py'),
             Path(__file__).with_name('selector_resident_allocator_control_20260922.py'),
             Path(__file__).with_name('direct_bounded_host_selection_prices_20260921.py')]},
        scope='One empty first tile and one borrowed/read-only selected-index tile. Timed source emit/list completion excludes validation, input creation and returned-array destruction. Default timing and separately forwarding-wrapped allocation evidence, seven counterbalanced repeats. Not a whole-GWAS speedup or a reusable price bank.')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
