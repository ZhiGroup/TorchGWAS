"""Host decode cost of a PGEN sample selection: reader subset against whole rows.

A dropped subject (missing_phenotype='drop_subject') is a sample selection.
The native reader decodes a selection by expanding every sample, gathering
the kept columns and mapping them through the call table (three passes, two
full-width temporaries); the whole cohort is one pass straight into the
output. Whole rows plus a device-side missing mask keep the one-pass decode
and pay |dropped| / n more transfer and GPU rows instead.

Per configuration, one thread fills `--chunks` chunks of `--chunk` variants
(int8 hard calls, the default transport) and reports seconds per chunk,
after one warm pass over the same range so the page cache is equal.

    python benchmarks/pgen_subset_decode_20260928.py --pgen /data/zxie3/torchgwas_bench/full_scale_k512_20260926/input.pgen
"""
import argparse
import json
import os
import time

import numpy as np


def fill_seconds(source, starts, chunk):
    with source.native_reader_session() as read_into:
        out = np.empty((chunk, source.native_row_width), dtype=source.native_transfer_dtype)
        for start in starts:  # warm the page cache and the thread's reader
            read_into(start, start + chunk, out)
        began = time.perf_counter()
        for start in starts:
            read_into(start, start + chunk, out)
        return (time.perf_counter() - began) / len(starts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pgen', required=True)
    parser.add_argument('--chunk', type=int, default=4096)
    parser.add_argument('--chunks', type=int, default=24)
    parser.add_argument('--drop', type=float, nargs='+', default=[0.01, 0.05, 0.5])
    args = parser.parse_args()
    os.environ.setdefault('TORCHGWAS_PGEN_BACKEND', 'native')
    os.environ['TORCHGWAS_PGEN_PACKED'] = '0'
    from torchgwas.pgen import PgenGenotype
    source = PgenGenotype(args.pgen, mode='hardcall', reader_workers=1)
    n = source.shape[0]
    starts = [1_000_000 + i * args.chunk for i in range(args.chunks)]
    ids = np.asarray(source.sample_ids, dtype=object)
    rng = np.random.default_rng(20260928)
    rows = [dict(config='all_samples', kept=n, seconds_per_chunk=round(fill_seconds(source, starts, args.chunk), 5))]
    for fraction in args.drop:
        kept = np.sort(rng.choice(n, size=int(round(n * (1 - fraction))), replace=False))
        source.select_samples(ids[kept])
        rows.append(dict(config=f'subset_drop_{fraction}', kept=int(kept.size),
                         seconds_per_chunk=round(fill_seconds(source, starts, args.chunk), 5)))
    source.select_samples(None)
    for row in rows:
        row['ns_per_call'] = round(1e9 * row['seconds_per_chunk'] / (args.chunk * n), 3)
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
