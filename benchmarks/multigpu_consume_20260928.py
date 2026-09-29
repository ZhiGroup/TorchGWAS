"""Does the writer or the scan limit K = 512 on several GPUs? Consume the shards' chunks, write nothing.

The same variant-sharded min-p scan the API runs (linear_scan_multigpu,
chunk 4096, 4 decoders per GPU, depth 4, the default statistics), on 1, 2
and 4 of `--devices`, its chunks consumed and dropped on the main thread.
Against the API's executor seconds for the same layout, the gap is what
writing costs; if one GPU and four take the same time here too, the limit
is upstream of the writer.

    python benchmarks/multigpu_consume_20260928.py --data /data/zxie3/torchgwas_bench/full_scale_k512_20260926 \\
        --devices cuda:4 cuda:5 cuda:6 cuda:7
"""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--variants', type=int, default=None)
    args = parser.parse_args()
    os.environ.setdefault('TORCHGWAS_PGEN_BACKEND', 'native')
    # Each shard's own scan profile (fetch, GPU compute, result waits).
    os.environ['TORCHGWAS_SCAN_PROFILE'] = '1'
    from torchgwas import sumstats_tiled
    views = []
    original_init = sumstats_tiled.ScanSourceView.__init__

    def recording_init(self, source):
        original_init(self, source)
        views.append(self)
    sumstats_tiled.ScanSourceView.__init__ = recording_init
    from torchgwas.linear import linear_scan_multigpu
    from torchgwas.min_p import MinPReduction
    from torchgwas.pgen import PgenGenotype
    # The API's source: 16 decoders (load_genotype's reader_workers).
    source = PgenGenotype(args.data/'input.pgen', mode='hardcall', metadata_cache_dir=args.data/'metadata_cache',
                          reader_workers=16)
    phenotype = np.load(args.data/'phenotype.npy').astype(np.float32)
    covariates = np.load(args.data/'covariates.npy').astype(np.float32)
    span = None if args.variants is None else (0, args.variants)
    for count in (1, 2, 4):
        devices = args.devices[:count]
        started = time.perf_counter()
        iterator, _ = linear_scan_multigpu(
            source, phenotype, covariates, devices=devices, chunk_size=4096, reader_workers=4 * count,
            prefetch_chunks=4, compute_p_values=False, ordered=False, reduction_factory=MinPReduction,
            compute_log10_p=True, log10_p_dtype='float32', variant_range=span, borrow_results=True)
        chunks = 0
        waited = 0.0
        mark = time.perf_counter()
        for _ in iterator:
            now = time.perf_counter()
            waited += now - mark
            chunks += 1
            mark = time.perf_counter()
        seconds = time.perf_counter() - started
        for view in views:
            profile = view.__dict__.get('_last_scan_profile', {})
            print(json.dumps({key: round(float(profile.get(key, 0)), 3) for key in (
                'fetch_seconds', 'copy_wait_seconds', 'result_wait_seconds', 'gpu_compute_milliseconds',
                'gpu_result_milliseconds', 'chunks')}), flush=True)
        views.clear()
        print(json.dumps(dict(devices=count, chunks=chunks, seconds=round(seconds, 2),
                              consumer_wait_seconds=round(waited, 2),
                              backend=getattr(source, '_last_scan_profile', {}).get('statistics_backend'),
                              encoding=source.native_encoding)), flush=True)


if __name__ == '__main__':
    main()
