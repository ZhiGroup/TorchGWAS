"""Measure a real source-bound structural admission miss and reuse."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_admission_cache import PgenAdmissionCache
from torchgwas.pgen_memory_layout import (
    memory_layout, compact_rechunk_memory_layout, compact_shifted_reader_envelope)
from torchgwas.pgen_reader import read_header


def measured(call):
    w=time.perf_counter();c=time.process_time()
    value=call()
    return value,dict(wall_seconds=time.perf_counter()-w,cpu_seconds=time.process_time()-c)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--path',required=True)
    parser.add_argument('--cache-dir',required=True)
    parser.add_argument('--report',required=True)
    parser.add_argument('--fine',type=int,default=128)
    args=parser.parse_args()
    path=Path(args.path).resolve(strict=True)
    identity=input_identity(path)
    header,header_time=measured(lambda:read_header(path))
    bases=[]
    layout,build_time=measured(lambda:memory_layout(path,args.fine,header=header,
        compact=True,_validated_bases_receiver=bases.append))
    base=bases[0].astype('uint32')
    base.flags.writeable=False
    cache=PgenAdmissionCache(args.cache_dir,identity,header,args.fine)
    before=cache.audit()
    _,stage_time=measured(lambda:cache.stage(layout,base))
    pending=cache.pending_bytes
    publication,publish_time=measured(lambda:cache.publish(successful=True))
    reuse=PgenAdmissionCache(args.cache_dir,identity,header,args.fine)
    hit,load_time=measured(reuse.load)
    if hit is None:raise RuntimeError('Structural cache miss after publication: '+repr(reuse.audit()))
    cached,cached_bases=hit
    for key in ('fine_payload_bytes','fine_prefix_bytes','fine_extra_workspace_bytes',
                'payload_cumulative_bytes'):
        if not np.array_equal(layout[key],cached[key]):
            raise RuntimeError('Cached vector differs: '+key)
    if not np.array_equal(base,cached_bases):
        raise RuntimeError('Cached base index differs')
    m=header.variant_ct
    capacities=[args.fine,8*args.fine]
    envelopes=[]
    for capacity in capacities:
        regular=compact_rechunk_memory_layout(layout,capacity)
        replay=compact_rechunk_memory_layout(cached,capacity)
        shifted=compact_shifted_reader_envelope(layout,(0,m),capacity)
        shifted_replay=compact_shifted_reader_envelope(cached,(0,m),capacity)
        if regular!=replay or shifted!=shifted_replay:
            raise RuntimeError('Cached admission envelope differs')
        envelopes.append(dict(capacity=capacity,fixed_reader_bytes=
            regular['chunks'][0]['record_payload_bytes'],
            shifted_reader_bytes=shifted[0],shifted_scratch_bytes=shifted[1]))
    report=dict(input=identity,samples=header.sample_ct,markers=m,fine=args.fine,
        logical_chunks=layout['logical_chunks'],retained_base_bytes=base.nbytes,
        pending_vector_bytes=pending,cache_file_bytes=reuse.path.stat().st_size,
        header=header_time,uncached_compact_build=build_time,cache_stage=stage_time,
        cache_publish=publish_time,cache_load=load_time,
        before=before,after=reuse.audit(),publication=publication,
        vectors_equal=True,bases_equal=True,envelopes=envelopes,
        scope='Ordered single-process source-index, structural-build, publication and hit diagnostic; no GWAS end-to-end speedup claim.')
    report_path=Path(args.report);report_path.parent.mkdir(parents=True,exist_ok=True)
    report_path.write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':main()
