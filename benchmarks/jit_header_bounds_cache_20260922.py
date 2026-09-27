"""Matched header-window replay on the real public two-shard JAGWAS fixture.

This isolates bounded source metadata construction. It does not time a GWAS,
GPU work, output delivery, or the analytical schedule comparison.
"""
import argparse
import json
from pathlib import Path
import statistics
import time

from torchgwas.analytical_plan_cache import input_identity
from torchgwas.detailed_calibration import sha256_file,source_identity
from torchgwas.pgen_reader import read_header
from torchgwas.pgen_work_bounds import PgenHeaderWork


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--trials',type=int,default=8)
    args=parser.parse_args()
    if not 2<=args.trials<=32:raise ValueError('Two to 32 paired trials required')
    path=Path(args.input);header=read_header(path);identity=input_identity(path)
    starts=(512,2688);horizons=(256,512,768);sizes=(128,256)
    if header.variant_ct!=4097 or header.sample_ct!=2049:
        raise ValueError('This audit requires the public two-shard fixture')
    results=[];reference=None
    for trial in range(args.trials):
        for limit in ((0,64) if trial%2==0 else (64,0)):
            index=PgenHeaderWork(path,_prepared_header=(identity,header),max_cached_bounds=limit)
            began=time.perf_counter()
            windows=[index.window(start,start+horizon,size,max_chunks=32,max_records=65536)
                for horizon in horizons for size in sizes for start in starts]
            elapsed=time.perf_counter()-began
            if reference is None:reference=windows
            elif windows!=reference:raise AssertionError('Cached header windows changed source work')
            results.append(dict(trial=trial,cache_limit=limit,window_seconds=elapsed,
                bounds_cache=index.bounds_cache_info()))
    if input_identity(path)!=identity:raise ValueError('Input changed during matched replay')
    medians={str(limit):statistics.median(row['window_seconds'] for row in results
        if row['cache_limit']==limit) for limit in (0,64)}
    report=dict(input=str(path),input_identity=identity,input_sha256=sha256_file(path),
        source_sha256=source_identity(),benchmark_sha256=sha256_file(__file__),
        sequence=dict(starts=starts,horizons=horizons,chunk_sizes=sizes),
        results=results,median_seconds=medians,scope=__doc__)
    target=Path(args.out);target.parent.mkdir(parents=True,exist_ok=False)
    target.write_text(json.dumps(report,indent=2))
    print(json.dumps(dict(median_seconds=medians,cache_info=results[-1]['bounds_cache'])))


if __name__=='__main__':main()
