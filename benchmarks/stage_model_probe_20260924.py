"""Per-chunk stage timings across chunk sizes, for the online tuning model.

Usage: python stage_model_probe_20260924.py DATA OUT.json [--traits 8192]
       [--sizes 256 512 1024 2048 4096] [--block 6] [--readers 4] [--depth 4]

One single-GPU significant-pair scan over DATA (from
empirical_layout_bench_20260923.py make). The chunk size changes every
--block chunks to a size drawn at random (seeded), so load bursts are not
aligned with sizes. Every chunk after two warmup chunks is measured by
InitialChunkMeasurements: reader wall, reader CPU and scheduler wait (from
schedstat), and CUDA spans for H2D, conversion and statistics. All output is
real work. The JSON keeps every observation; the printout checks the model

    stage_time(c) = a + b*c

per stage: the fit across all sizes, and a prediction of each size's chunk
period from a fit on two sizes only.
"""
import argparse
import json
import os
import random
import socket
import time
from pathlib import Path

import numpy as np

ENV = dict(TORCHGWAS_PGEN_BACKEND='native', TORCHGWAS_PGEN_PACKED='0', TORCHGWAS_NATIVE_STATS='0',
           TORCHGWAS_SCAN_PROFILE='0', NUMPY_MADVISE_HUGEPAGE='0')


def fit(sizes, values):
    """Least squares values = a + b*sizes; returns a, b, R^2."""
    x = np.asarray(sizes, dtype=float)
    y = np.asarray(values, dtype=float)
    design = np.column_stack([np.ones_like(x), x])
    (a, b), *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y-(a+b*x)
    total = ((y-y.mean())**2).sum()
    return float(a), float(b), float(1-(residual**2).sum()/total) if total > 0 else float('nan')


def analyse(rows, readers):
    by_size = {}
    for row in rows:
        by_size.setdefault(row['size'], []).append(row)
    sizes = sorted(by_size)
    stages = ('read_wall', 'read_cpu', 'read_wait', 'h2d', 'statistics', 'consumer')
    print('size  chunks  period_ms  rows_per_s  ' + '  '.join(f'{s}_ms' for s in stages))
    observed = {}
    for size in sizes:
        group = by_size[size]
        periods = [r['period'] for r in group if r['period'] is not None]
        observed[size] = float(np.median(periods)) if periods else None
        medians = [np.nanmedian([np.nan if r[s] is None else r[s] for r in group])*1e3 for s in stages]
        print(f'{size:5d} {len(group):6d}  {observed[size]*1e3 if observed[size] else float("nan"):9.2f}  '
              f'{size/observed[size] if observed[size] else float("nan"):10.0f}  '
              + '  '.join(f'{m:9.3f}' for m in medians))
    print('\nfit over all sizes: stage = a + b*c   (a in ms, b in us/variant)')
    models = {}
    for stage in stages:
        pairs = [(r['size'], r[stage]) for r in rows if r[stage] is not None]
        if len(pairs) < 3:
            continue
        a, b, r2 = fit(*zip(*pairs))
        models[stage] = (a, b)
        print(f'  {stage:10s} a={a*1e3:8.3f}  b={b*1e6:8.3f}  R2={r2:6.3f}')
    # Predict each size's period from two sizes: bottleneck of the stages.
    train = sizes[1:3] if len(sizes) >= 3 else sizes[:2]
    print(f'\nperiod predicted from sizes {train} only (bottleneck of read/{readers} readers, h2d, statistics, consumer):')
    two = {}
    for stage in ('read_wall', 'h2d', 'statistics', 'consumer'):
        pairs = [(r['size'], r[stage]) for r in rows if r[stage] is not None and r['size'] in train]
        if len(pairs) >= 2:
            two[stage] = fit(*zip(*pairs))[:2]
    for size in sizes:
        terms = {stage: (a+b*size)/(readers if stage == 'read_wall' else 1) for stage, (a, b) in two.items()}
        bottleneck = max(terms, key=terms.get)
        predicted = terms[bottleneck]
        seen = observed[size]
        error = (predicted/seen-1)*100 if seen else float('nan')
        print(f'  {size:5d}: predicted {predicted*1e3:8.2f} ms ({bottleneck}), observed {seen*1e3 if seen else float("nan"):8.2f} ms, '
              f'error {error:+6.1f}%{"  (trained)" if size in train else ""}')
    return models


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--traits', type=int, default=8192)
    parser.add_argument('--sizes', type=int, nargs='+', default=[256, 512, 1024, 2048, 4096])
    parser.add_argument('--block', type=int, default=6)
    parser.add_argument('--readers', type=int, default=4)
    parser.add_argument('--depth', type=int, default=4)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--threshold', type=float, default=1e-5)
    parser.add_argument('--seed', type=int, default=20260924)
    parser.add_argument('--analyse-only', action='store_true')
    args = parser.parse_args()
    if args.analyse_only:
        record = json.loads(args.out.read_text())
        analyse(record['chunks'], record['readers'])
        return
    os.environ.update(ENV)
    from torchgwas.adaptive_chunks import AlignedChunkSizeControl, InitialChunkMeasurements
    from torchgwas.io import load_genotype
    from torchgwas.linear import linear_scan_streaming_chunks
    from torchgwas.reduce import SignificantPairs
    genotype = load_genotype(str(args.data/'input.pgen'), genotype_format='pgen', pgen_mode='hardcall',
                             reader_workers=args.readers)[0]
    y = np.ascontiguousarray(np.load(args.data/'phenotype.npy', mmap_mode='r')[:, :args.traits], dtype=np.float32)
    covariates = np.load(args.data/'covariates.npy').astype(np.float32)
    sizes = sorted(args.sizes)
    control = AlignedChunkSizeControl(sizes, initial=sizes[0])
    rng = random.Random(args.seed)
    calls = [0]

    def selector(start, stop, capacity):
        if calls[0] % args.block == 0:
            control.set_size(rng.choice(sizes))
        calls[0] += 1
        return control(start, stop, capacity)

    window = InitialChunkMeasurements([args.device], max_chunks_per_device=10**6, warmup_chunks=2, stride=1,
                                      max_window_seconds=1e6, cuda_events=True)
    started = time.perf_counter()
    chunks, _ = linear_scan_streaming_chunks(
        genotype, y, covariates, chunk_size=max(sizes), device=args.device, compute_dtype='float32',
        reader_workers=args.readers, prefetch_chunks=args.depth, compute_p_values=False,
        significance=SignificantPairs(args.threshold), significance_n_traits=y.shape[1], return_df=True,
        _chunk_size_selector=selector, _chunk_observer=window)
    pairs = sum(1 for _ in chunks)
    wall = time.perf_counter()-started
    observations = sorted(window.snapshot()['observations'], key=lambda row: row['completed'])
    rows, previous = [], None
    for row in observations:
        size = row['end']-row['start']
        cuda = row.get('cuda') or {}
        period = (row['completed']-previous['completed']
                  if previous is not None and previous['end']-previous['start'] == size else None)
        rows.append(dict(start=row['start'], size=size, completed=row['completed'], period=period,
                         read_wall=row['read_finished']-row['read_started'], read_cpu=row['read_cpu_seconds'],
                         read_wait=row['read_runnable_wait_seconds'], h2d=cuda.get('h2d'),
                         conversion=cuda.get('conversion'), statistics=cuda.get('statistics_and_reduction'),
                         consumer=row['consumer_seconds']))
        previous = row
    record = dict(host=socket.gethostname(), device=args.device, samples=int(y.shape[0]), traits=int(y.shape[1]),
                  variants=int(genotype.shape[1]), sizes=sizes, block=args.block, readers=args.readers,
                  depth=args.depth, wall_seconds=wall, selected_pair_chunks=pairs,
                  load_average=os.getloadavg(), chunks=rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record))
    print(f"{record['host']} N={record['samples']} K={record['traits']} M={record['variants']} "
          f"wall {wall:.1f} s, load {record['load_average'][0]:.0f}, {len(rows)} measured chunks")
    analyse(rows, args.readers)


if __name__ == '__main__':
    main()
