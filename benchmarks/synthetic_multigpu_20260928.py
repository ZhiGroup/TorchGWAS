"""Does the scan pipeline itself serialize across GPUs? A source with no decode and no disk.

A real PGEN of `--samples` samples and a handful of variants, presented as
`--variants` variants whose packed rows are memcpy'd from a random pool
(PgenGenotype.native_reader_session overridden): no page cache, no decoder,
no writer. linear_scan_multigpu (the API's variant shards: min-p, chunk
4096, 4 fill workers per GPU, depth 4, the default statistics) on 1, 2, ...
of `--devices`, its chunks consumed and dropped. If more GPUs do not finish
sooner, the serial stage is the scan's own per-chunk host work.

    python benchmarks/synthetic_multigpu_20260928.py --devices cuda:1 cuda:2 cuda:7
"""
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tests'))


def synthetic_source(directory, samples, variants, traits_seed=20260928):
    from test_pgen_native_reader import write_pgen
    from torchgwas.pgen import PgenGenotype
    rng = np.random.default_rng(traits_seed)
    path = Path(directory) / 'tiny.pgen'
    write_pgen(path, rng.integers(0, 3, size=(8, samples)).astype(np.uint8))
    path.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
        f'1\t{i + 1}\tv{i}\tA\tC\n' for i in range(8)))
    path.with_suffix('.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(samples)))

    class Synthetic(PgenGenotype):
        @contextmanager
        def native_reader_session(self):
            pool = self._pool

            def read_into(start, end, out):
                rows = end - start
                first = start % (pool.shape[0] - rows)
                np.copyto(out, pool[first:first + rows])
                return out
            yield read_into

    source = Synthetic(path, mode='hardcall', reader_workers=16)
    source._raw_n_variants = variants
    pool_rows = 16_384
    if source.native_encoding == 'pgen_2bit':
        source._pool = rng.integers(0, 256, size=(pool_rows, source.native_row_width), dtype=np.uint8)
        # Code 3 is missing; clearing a bit per pair keeps most calls observed.
        source._pool &= np.uint8(0b10111011)
    else:  # int8 transport (TORCHGWAS_PGEN_PACKED=0): 0/1/2 calls, -9 missing
        pool = rng.integers(0, 3, size=(pool_rows, source.native_row_width), dtype=np.int8)
        pool[rng.random(pool.shape) < 0.01] = -9
        source._pool = pool
    return source


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--devices', nargs='+', required=True)
    parser.add_argument('--samples', type=int, default=22_250)
    parser.add_argument('--variants', type=int, default=400_000)
    parser.add_argument('--traits', type=int, default=512)
    parser.add_argument('--mode', default='min-p', choices=['min-p', 'min-p-no-tail', 'full'])
    parser.add_argument('--counts', type=int, nargs='+', default=None)
    parser.add_argument('--profile', action='store_true', help='print every shard profile (TORCHGWAS_SCAN_PROFILE)')
    args = parser.parse_args()
    os.environ.setdefault('TORCHGWAS_PGEN_BACKEND', 'native')
    if args.profile:
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
    rng = np.random.default_rng(1)
    phenotype = rng.normal(size=(args.samples, args.traits)).astype(np.float32)
    covariates = rng.normal(size=(args.samples, 10)).astype(np.float32)
    with tempfile.TemporaryDirectory() as directory:
        source = synthetic_source(directory, args.samples, args.variants)
        from torchgwas.reduce import VariantReduction
        extra = {'min-p': dict(reduction_factory=MinPReduction, compute_log10_p=True, log10_p_dtype='float32'),
                 # The same |t| ranking with no -log10 P tail at all.
                 'min-p-no-tail': dict(reduction_factory=lambda: VariantReduction('min-p')),
                 'full': dict(compute_log10_p=True, log10_p_dtype='float32')}[args.mode]
        for count in (args.counts or range(1, len(args.devices) + 1)):
            for attempt in ('warm', 'timed'):
                span = (0, 40_960 if attempt == 'warm' else args.variants)
                started = time.perf_counter()
                iterator, _ = linear_scan_multigpu(
                    source, phenotype, covariates, devices=args.devices[:count], chunk_size=4096,
                    reader_workers=4 * count, prefetch_chunks=4, compute_p_values=False, ordered=False,
                    variant_range=span, borrow_results=True, **extra)
                chunks = sum(1 for _ in iterator)
                seconds = time.perf_counter() - started
            if args.profile:
                keys = ('fetch_seconds', 'copy_wait_seconds', 'result_wait_seconds', 'gpu_compute_milliseconds',
                        'gpu_result_milliseconds', 'chunks')
                for view in views[-count:]:
                    profile = view.__dict__.get('_last_scan_profile', {})
                    print(json.dumps({key: round(float(profile.get(key, 0)), 3) for key in keys}), flush=True)
            views.clear()
            if args.profile:
                import torch
                for device in args.devices[:count]:
                    stats = torch.cuda.memory_stats(device)
                    print(json.dumps(dict(device=device, **{key: stats.get(key, 0) for key in (
                        'num_alloc_retries', 'num_device_alloc', 'num_device_free', 'num_sync_all_streams',
                        'num_ooms')}, reserved_gb=round(stats.get('reserved_bytes.all.peak', 0) / 2**30, 2))),
                          flush=True)
            print(json.dumps(dict(devices=count, mode=args.mode, chunks=chunks, seconds=round(seconds, 3),
                                  us_per_variant=round(1e6 * seconds / args.variants, 3),
                                  backend=getattr(source, '_last_scan_profile', {}).get('statistics_backend'))),
                  flush=True)


if __name__ == '__main__':
    main()
