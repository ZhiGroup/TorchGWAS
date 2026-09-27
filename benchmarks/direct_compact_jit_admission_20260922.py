"""Compare exact full and compact PGEN memory envelopes, without payload reads."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from torchgwas.adaptive_candidate import _reader_envelope
from torchgwas.analytical_plan_cache import input_identity
from torchgwas.decoder_work import native_reader_workspace
from torchgwas.linear import multigpu_variant_ranges
from torchgwas.pgen_memory_layout import (memory_layout,rechunk_memory_layout,
    compact_rechunk_memory_layout,compact_shifted_reader_envelope)
from torchgwas.pgen_reader import read_header


def digest(path):
    with open(path,'rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pgen',required=True)
    parser.add_argument('--fine',required=True,type=int)
    parser.add_argument('--capacity',required=True,type=int)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    if args.fine<1 or args.capacity<args.fine or args.capacity%args.fine:
        raise ValueError('Capacity must be a multiple of the positive fine grid')
    path=Path(args.pgen).resolve(strict=True);identity=input_identity(path)
    out=Path(args.out)
    if out.exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    def timed(name,fn):
        wall=time.perf_counter();cpu=time.process_time();value=fn()
        rows.append(dict(name=name,wall_seconds=time.perf_counter()-wall,
            cpu_seconds=time.process_time()-cpu))
        return value
    header=timed('read_header',lambda:read_header(path))
    compact=timed('compact_fine',lambda:memory_layout(path,args.fine,header=header,compact=True))
    full=timed('full_fine',lambda:memory_layout(path,args.fine,header=header))
    spans=[(0,header.variant_ct)]
    shards=multigpu_variant_ranges(header.variant_ct,args.capacity,2)
    if len(shards)==2:spans.extend(shards)
    comparisons=[]
    for lo,hi in spans:
        def full_bounds():
            fixed=rechunk_memory_layout(full,args.capacity,(lo,hi))
            shifted=_reader_envelope(full['chunks'][lo//args.fine:(hi+args.fine-1)//args.fine],
                capacity=args.capacity,fine=args.fine)
            return dict(fixed=native_reader_workspace(fixed['chunks']),shifted=shifted,
                payload=fixed['record_payload_bytes'])
        def compact_bounds():
            fixed=compact_rechunk_memory_layout(compact,args.capacity,(lo,hi))
            shifted=compact_shifted_reader_envelope(compact,(lo,hi),args.capacity)
            return dict(fixed=native_reader_workspace(fixed['chunks']),shifted=shifted,
                payload=fixed['record_payload_bytes'])
        old=timed(f'full_envelope_{lo}_{hi}',full_bounds)
        new=timed(f'compact_envelope_{lo}_{hi}',compact_bounds)
        if old!=new:raise AssertionError(f'Memory envelope differs for {lo}:{hi}: {old} != {new}')
        comparisons.append(dict(variant_range=[lo,hi],fixed_read_bytes=old['fixed'][0],
            fixed_scratch_bytes=old['fixed'][1],shifted_read_bytes=old['shifted'][0],
            shifted_scratch_bytes=old['shifted'][1],payload_bytes=old['payload']))
    if input_identity(path)!=identity:raise ValueError('PGEN changed during comparison')
    source={name:digest(Path(__file__).resolve().parents[1]/'src'/'torchgwas'/name)
        for name in ('pgen_memory_layout.py','adaptive_candidate.py','trait_candidate_space.py')}
    report=dict(scope=__doc__,input_identity=identity,samples=header.sample_ct,
        variants=header.variant_ct,fine=args.fine,capacity=args.capacity,
        logical_fine_chunks=len(full['chunks']),compact_vector_bytes=sum(
            compact[key].nbytes for key in ('fine_payload_bytes','fine_prefix_bytes',
                'fine_extra_workspace_bytes','payload_cumulative_bytes')),
        rows=rows,comparisons=comparisons,script_sha256=digest(__file__),source_sha256=source,
        caveat='One ordered process. Shared load, cache and allocator state are uncontrolled; component cost is not a full-job speedup.')
    with out.open('x') as handle:json.dump(report,handle,indent=2)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
