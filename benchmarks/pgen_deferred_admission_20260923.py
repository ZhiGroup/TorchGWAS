"""Ordered large-PGEN exact admission with deferred full base indexing."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_admission_cache import PgenAdmissionCache
from torchgwas.pgen_memory_layout import (memory_layout,compact_rechunk_memory_layout,
    compact_shifted_reader_envelope)
from torchgwas.pgen_reader import read_header,ld_safe_start
from torchgwas.pgen_work_bounds import PgenHeaderWork


def timed(call):
    wall=time.perf_counter();cpu=time.process_time()
    value=call()
    return value,dict(wall_seconds=time.perf_counter()-wall,
                      cpu_seconds=time.process_time()-cpu)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--path',required=True);p.add_argument('--report',required=True)
    p.add_argument('--fine',type=int,default=128)
    p.add_argument('--cache-dir')
    args=p.parse_args();path=Path(args.path).resolve(strict=True)
    header,header_time=timed(lambda:read_header(path))
    observations=[];results={}
    for name,defer in [('full',False),('deferred',True),('deferred',True),('full',False)]:
        layout,seconds=timed(lambda:memory_layout(path,args.fine,header=header,
                                                   compact=True,_defer_bases=defer))
        observations.append(dict(name=name,**seconds))
        results.setdefault(name,layout)
    full=results['full'];deferred=results['deferred']
    for field in ('fine_payload_bytes','fine_prefix_bytes','fine_extra_workspace_bytes',
                  'payload_cumulative_bytes'):
        np.testing.assert_array_equal(full[field],deferred[field])
    m=header.variant_ct
    envelopes=[]
    for capacity in (args.fine,8*args.fine):
        fixed=compact_rechunk_memory_layout(full,capacity)
        shifted=compact_shifted_reader_envelope(full,(0,m),capacity)
        assert fixed==compact_rechunk_memory_layout(deferred,capacity)
        assert shifted==compact_shifted_reader_envelope(deferred,(0,m),capacity)
        envelopes.append(dict(capacity=capacity,fixed_read_bytes=
            fixed['chunks'][0]['record_payload_bytes'],shifted_read_bytes=shifted[0],
            shifted_scratch_bytes=shifted[1]))
    identity=input_identity(path)
    worker,worker_time=timed(lambda:PgenHeaderWork(path,_prepared_header=(identity,header)))
    assert worker._bases is None
    probe=worker.window(0,min(m,1024),128,max_records=2048,max_chunks=8)
    cache_evidence=None
    if args.cache_dir:
        cache=PgenAdmissionCache(args.cache_dir,identity,header,args.fine)
        cache.stage(deferred,None)
        pending=cache.pending_bytes
        status,publish_time=timed(lambda:cache.publish(successful=True))
        assert status=='stored'
        reuse=PgenAdmissionCache(args.cache_dir,identity,header,args.fine)
        hit,load_time=timed(reuse.load)
        assert hit is not None and reuse.status=='hit'
        cached,bases=hit
        for key in ('fine_payload_bytes','fine_prefix_bytes','fine_extra_workspace_bytes',
                    'payload_cumulative_bytes'):
            np.testing.assert_array_equal(cached[key],deferred[key])
        for start in np.linspace(0,m-1,100,dtype=np.int64):
            selected=int(bases[np.searchsorted(bases,bases.dtype.type(start),side='right')-1])
            assert selected==ld_safe_start(header.vrtypes,int(start))
        cache_evidence=dict(pending_vector_bytes=pending,
            publication=publish_time,publication_audit=cache.audit(),
            load=load_time,load_audit=reuse.audit(),artifact_bytes=reuse.path.stat().st_size,
            sampled_predecessors_equal=True)
    report=dict(input=identity,samples=header.sample_ct,markers=m,
        initial_header=header_time,observations=observations,
        productive_header_constructor=worker_time,
        productive_probe_chunks=len(probe['chunks']),
        fine_vector_bytes=sum(full[key].nbytes for key in ('fine_payload_bytes',
            'fine_prefix_bytes','fine_extra_workspace_bytes','payload_cumulative_bytes')),
        vectors_equal=True,envelopes=envelopes,cache=cache_evidence,
        scope='Ordered structural component diagnostic with the same prepared header; no end-to-end speedup or decoder price claim.')
    dest=Path(args.report);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':main()
