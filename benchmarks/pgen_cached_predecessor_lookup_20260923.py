"""Measure dtype-matched cached PGEN predecessor lookup on a real index."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.pgen_admission_cache import PgenAdmissionCache
from torchgwas.pgen_reader import read_header


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--path',required=True);p.add_argument('--cache-dir',required=True)
    p.add_argument('--report',required=True);p.add_argument('--fine',type=int,default=128)
    args=p.parse_args();path=Path(args.path).resolve(strict=True)
    header=read_header(path)
    cache=PgenAdmissionCache(args.cache_dir,input_identity(path),header,args.fine)
    hit=cache.load()
    if hit is None:raise RuntimeError('Expected aligned structural cache hit')
    bases=hit[1];assert bases.flags.aligned and not bases.flags.writeable
    probe=header.variant_ct//2;typed=bases.dtype.type(probe)
    start=time.process_time()
    for _ in range(1000):
        typed_index=np.searchsorted(bases,typed,side='right')
    typed_cpu=time.process_time()-start
    start=time.process_time()
    plain_index=np.searchsorted(bases,probe,side='right')
    plain_cpu=time.process_time()-start
    assert typed_index==plain_index
    report=dict(input=input_identity(path),base_count=len(bases),base_dtype=str(bases.dtype),
        aligned=bases.flags.aligned,typed_1000_cpu_seconds=typed_cpu,
        python_int_1_cpu_seconds=plain_cpu,same_result=True,cache=cache.audit(),
        scope='One ordered diagnostic under shared load; dtype and alignment behavior, not GWAS end-to-end timing.')
    dest=Path(args.report);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':main()
