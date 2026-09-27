"""Ordered PGEN metadata components for JIT admission, not a GWAS speedup test."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_memory_layout import memory_layout
from torchgwas.pgen_reader import read_header
from torchgwas.pgen_work_bounds import PgenHeaderWork


def digest(path):
    with open(path,'rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pgen',required=True)
    parser.add_argument('--chunk',required=True,type=int)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    path=Path(args.pgen).resolve(strict=True)
    out=Path(args.out)
    if out.exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    def timed(name,fn):
        wall=time.perf_counter();cpu=time.process_time();value=fn()
        rows.append(dict(name=name,wall_seconds=time.perf_counter()-wall,
            cpu_seconds=time.process_time()-cpu))
        return value
    identity=input_identity(path)
    header=timed('read_header',lambda:read_header(path))
    bases=[]
    layout=timed('memory_layout',lambda:memory_layout(path,args.chunk,header=header,
        _validated_bases_receiver=bases.append))
    if len(bases)!=1:raise ValueError('Expected one certified LD-base array')
    compact_bases=bases[0].astype('uint32')
    compact_bases.flags.writeable=False
    prepared=(identity,header);certificate=(identity,header,compact_bases)
    certified=timed('certified_work_constructor',lambda:PgenHeaderWork(path,
        _prepared_header=prepared,_prepared_index=certificate))
    fallback=timed('revalidating_work_constructor',lambda:PgenHeaderWork(path,
        _prepared_header=prepared))
    lo=min(args.chunk,header.variant_ct-1)
    hi=min(lo+2*args.chunk,header.variant_ct)
    if lo>=hi:lo=0;hi=header.variant_ct
    if hi-lo>65536:hi=lo+65536
    if certified.window(lo,hi,args.chunk)!=fallback.window(lo,hi,args.chunk):
        raise AssertionError('Certified and revalidated source windows differ')
    if input_identity(path)!=identity:raise ValueError('Input changed during benchmark')
    source={name:digest(Path(__file__).resolve().parents[1]/'src'/'torchgwas'/name)
        for name in ('pgen_reader.py','pgen_memory_layout.py','pgen_work_bounds.py')}
    report=dict(scope=__doc__,input_identity=identity,samples=header.sample_ct,
        variants=header.variant_ct,chunk_markers=args.chunk,fine_chunks=len(layout['chunks']),
        retained_base_bytes=int(compact_bases.nbytes),rows=rows,
        validated_window=[lo,hi],script_sha256=digest(__file__),source_sha256=source,
        caveat='One ordered process on server-local input; page cache and shared load are uncontrolled. Not a matched cold test or whole-GWAS speedup.')
    with out.open('x') as handle:json.dump(report,handle,indent=2)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
