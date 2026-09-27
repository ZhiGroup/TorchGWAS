"""Does a JAGWAS join stage cost throughput?  Experiment only, not a product path.

JAGWAS in torchGWAS runs only when the full panel and its factor fit every GPU
(variant shards). This asks the what-if: split the panel over G GPUs and join
each chunk's t as its own pipeline stage. Does the extra stage hide behind the
scan, or does it cost throughput?

Synthetic and decode-free, so only the GPU side is compared. The genotype is a
pinned uint8 chunk copied to each GPU that needs it, and the factor is a random
lower triangle (timing does not depend on its values).

  vshards        each GPU takes every G-th chunk with the full panel: t = X Y,
                 then the block-triangular FP64 projection, inline (as the product does)
  vshards_stream same, with the projection on a second stream
  split_inline   GPU b holds traits [s_b, e_b). Every GPU scores every chunk (t_b = X Y_b).
                 GPU b pulls t_g for g < b and computes sum_g t_g L^-1[b, g]' and ||.||^2,
                 on the scoring stream
  split_stream   same, with the join on a second stream, issued one chunk behind
  *_noproj       scoring only (plus the split's peer pulls), which is the stage-1 floor

Tile bounds: 'equal' K_b = K/G; 'triangle', an equal area of the L^-1 triangle per GPU (FP64
projection balance); 'balanced', equal a*K_b + p*area_b, with the scoring cost per column (a) and the
projection cost per unit area (p) measured on the first GPU.

    $PY benchmarks/jagwas_join_stage_bench_20260925.py --devices cuda:0,cuda:1 --traits 8192 --out results/jagwas_join_20260925/HOST.jsonl
"""
import argparse
import json
import math
import threading
import time
from pathlib import Path

import torch

from torchgwas.jagwas_projection import triangular_blocks

MODES = ('vshards', 'vshards_stream', 'vshards_noproj', 'split_inline', 'split_stream', 'split_noproj')


class Board:
    """Cross-thread handoff of CUDA events keyed by chunk."""

    def __init__(self, timeout=600.0):
        self.cond = threading.Condition()
        self.items = {}
        self.failed = None
        self.timeout = timeout

    def put(self, key, value):
        with self.cond:
            self.items[key] = value
            self.cond.notify_all()

    def get(self, key):
        deadline = time.monotonic() + self.timeout
        with self.cond:
            while key not in self.items:
                if self.failed is not None:
                    raise RuntimeError('peer failed') from self.failed
                left = deadline - time.monotonic()
                if left <= 0:
                    raise TimeoutError(key)
                self.cond.wait(left)
            return self.items[key]

    def fail(self, error):
        with self.cond:
            self.failed = error
            self.cond.notify_all()


def measured_costs(device, *, samples, traits, chunk, repeats=5):
    """Seconds per trait column (scoring) and per unit of L^-1 triangle area (projection)."""
    x = torch.randint(0, 3, (chunk, samples), dtype=torch.uint8, device=device)
    panel = torch.randn(samples, traits, device=device)
    factor = torch.randn(traits, traits, dtype=torch.float64, device=device).tril_()
    blocks = [(e, factor[s:e, :e]) for s, e in triangular_blocks(traits)]

    def score():
        return x.float() @ panel

    def project(t):
        t = t.double()
        return sum(((t[:, :e] @ rows.T) ** 2).sum(dim=1) for e, rows in blocks)

    def median(fn):
        fn()
        torch.cuda.synchronize(device)
        values = []
        for _ in range(repeats):
            start = time.perf_counter()
            fn()
            torch.cuda.synchronize(device)
            values.append(time.perf_counter() - start)
        return sorted(values)[len(values) // 2]

    scores = score()
    per_column = median(score) / traits
    per_area = median(lambda: project(scores)) / (traits * traits / 2)
    return per_column, per_area


def tile_bounds(traits, count, kind, costs=None):
    if kind == 'balanced':
        # GPU b costs a*K_b + p*(e_b^2 - s_b^2)/2; bisect the common cost.
        a, p = costs

        def edges(target):
            out, start = [0], 0
            for _ in range(count):
                lo, hi = start, traits
                while lo < hi:
                    mid = (lo + hi + 1) // 2
                    if a * (mid - start) + p * (mid * mid - start * start) / 2 <= target:
                        lo = mid
                    else:
                        hi = mid - 1
                start = lo
                out.append(start)
            return out

        low, high = 0.0, a * traits + p * traits * traits / 2
        for _ in range(60):
            middle = (low + high) / 2
            low, high = (middle, high) if edges(middle)[-1] < traits else (low, middle)
        found = edges(high)
        found[-1] = traits
        return list(zip(found[:-1], found[1:]))
    if kind == 'equal':
        width = math.ceil(traits / count)
        edges = [min(traits, b * width) for b in range(count + 1)]
    else:
        # Equal triangle area per GPU (the diagonal block is sub-blocked too):
        # cumulative work to e is e^2 / 2, so edges at K sqrt(b/G).
        edges = [round(traits * math.sqrt(b / count)) for b in range(count + 1)]
    return list(zip(edges[:-1], edges[1:]))


def build(mode, devices, *, samples, traits, chunk, tiles, depth, costs=None):
    split = mode.startswith('split')
    bounds = tile_bounds(traits, len(devices), tiles, costs) if split else [(0, traits)] * len(devices)
    state = []
    for b, device in enumerate(devices):
        device = torch.device(device)
        start, end = bounds[b]
        width = end - start
        gen = torch.Generator(device=device).manual_seed(100 + b)
        panel = torch.randn(samples, width, device=device, generator=gen) / math.sqrt(samples)
        # Rows [start, end) of a random lower-triangular L^-1, columns [0, end).
        factor = torch.randn(width, end, dtype=torch.float64, device=device, generator=gen)
        factor[:, start:] = factor[:, start:].tril()
        factor[:, start:].diagonal().abs_().add_(1.0)
        state.append(dict(
            device=device, bounds=(start, end), panel=panel, factor=factor,
            blocks=[(e, factor[s:e, :e]) for s, e in triangular_blocks(traits)] if not split else None,
            genotype=[torch.empty((chunk, samples), dtype=torch.uint8, device=device) for _ in range(depth)],
            scores=[torch.empty((chunk, width), device=device) for _ in range(depth)],
            received={g: torch.empty((chunk, bounds[g][1] - bounds[g][0]), device=device) for g in range(b)}
            if split else {},
            partial=[torch.empty((chunk, 1), dtype=torch.float64, pin_memory=True) for _ in range(depth)],
            compute=torch.cuda.Stream(device), join=torch.cuda.Stream(device)))
    return state, bounds


def pipeline(mode, state, host, chunks, depth):
    count = len(state)
    split = mode.startswith('split')
    project = not mode.endswith('noproj')
    second = mode in ('split_stream', 'vshards_stream')
    lag = 1 if mode == 'split_stream' else 0
    board = Board()
    errors = []
    barrier = threading.Barrier(count + 1)

    def worker(b):
        try:
            me = state[b]
            torch.cuda.set_device(me['device'])
            compute = me['compute']
            join = me['join'] if second else compute
            mine = [i for i in range(chunks) if split or i % count == b]
            barrier.wait()

            done = []

            def score(position, i):
                # Split GPUs take every chunk, so position == i and peers agree on slots.
                slot = position % depth
                if second and position >= depth:
                    # The own join must be done reading this slot.
                    compute.wait_event(done[position - depth])
                if split and i >= depth:
                    # Peers must have pulled chunk i - depth before its slot is reused.
                    for peer in range(b + 1, count):
                        compute.wait_event(board.get(('used', i - depth, b, peer)))
                with torch.cuda.stream(compute):
                    me['genotype'][slot].copy_(host[i % len(host)], non_blocking=True)
                    torch.matmul(me['genotype'][slot].float(), me['panel'], out=me['scores'][slot])
                    event = torch.cuda.Event()
                    event.record(compute)
                if split:
                    board.put(('t', i, b), event)
                return event

            def reduce(position, i, event):
                slot = position % depth
                with torch.cuda.stream(join):
                    if split:
                        start, width = me['bounds'][0], me['bounds'][1] - me['bounds'][0]
                        if project:
                            acc = torch.zeros((me['scores'][slot].shape[0], width), dtype=torch.float64,
                                              device=me['device'])
                        for g in range(b + 1):
                            if g == b:
                                join.wait_event(event)
                                source = me['scores'][slot]
                            else:
                                join.wait_event(board.get(('t', i, g)))
                                source = me['received'][g]
                                source.copy_(state[g]['scores'][slot], non_blocking=True)
                                used = torch.cuda.Event()
                                used.record(join)
                                board.put(('used', i, g, b), used)
                            if not project:
                                continue
                            wide = source.double()
                            if g == b:
                                # The diagonal block is itself lower triangular: row sub-blocks.
                                for r0, r1 in triangular_blocks(width):
                                    acc[:, r0:r1].addmm_(wide[:, :r1], me['factor'][r0:r1, start:start + r1].T)
                            else:
                                s, e = state[g]['bounds']
                                acc.addmm_(wide, me['factor'][:, s:e].T)
                        if project:
                            me['partial'][slot].copy_((acc * acc).sum(dim=1, keepdim=True), non_blocking=True)
                    elif project:
                        join.wait_event(event)
                        t = me['scores'][slot].double()
                        statistic = None
                        for e, rows in me['blocks']:
                            projected = t[:, :e] @ rows.T
                            piece = (projected * projected).sum(dim=1, keepdim=True)
                            statistic = piece if statistic is None else statistic + piece
                        me['partial'][slot].copy_(statistic, non_blocking=True)
                    finished = torch.cuda.Event()
                    finished.record(join)
                    done.append(finished)

            pending = []
            for position, i in enumerate(mine):
                pending.append((position, i, score(position, i)))
                if len(pending) > lag:
                    reduce(*pending.pop(0))
            for item in pending:
                reduce(*item)
            torch.cuda.synchronize(me['device'])
        except BaseException as error:
            errors.append(error)
            board.fail(error)
            barrier.abort()

    threads = [threading.Thread(target=worker, args=(b,), daemon=True) for b in range(count)]
    for thread in threads:
        thread.start()
    barrier.wait()
    start = time.perf_counter()
    for thread in threads:
        thread.join()
    seconds = time.perf_counter() - start
    if errors:
        raise errors[0]
    return seconds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--devices', required=True)
    parser.add_argument('--traits', default='8192')
    parser.add_argument('--samples', type=int, default=16384)
    parser.add_argument('--chunk', type=int, default=2048)
    parser.add_argument('--chunks', type=int, default=48)
    parser.add_argument('--modes', default=','.join(MODES))
    parser.add_argument('--tiles', default='equal,triangle,balanced')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--depth', type=int, default=3)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    devices = args.devices.split(',')
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generator = torch.Generator().manual_seed(1)
    host = [torch.randint(0, 3, (args.chunk, args.samples), dtype=torch.uint8, generator=generator).pin_memory()
            for _ in range(2)]
    with out.open('a') as handle:
        for traits in map(int, args.traits.split(',')):
            costs = measured_costs(torch.device(devices[0]), samples=args.samples, traits=traits, chunk=args.chunk)
            for mode in args.modes.split(','):
                for tiles in (args.tiles.split(',') if mode.startswith('split') else ['-']):
                    state, bounds = build(mode, devices, samples=args.samples, traits=traits, chunk=args.chunk,
                                          tiles=tiles, depth=args.depth, costs=costs)
                    pipeline(mode, state, host, 2 * len(devices) + args.depth, args.depth)
                    times = sorted(pipeline(mode, state, host, args.chunks, args.depth) for _ in range(args.repeats))
                    seconds = times[len(times) // 2]
                    row = dict(device=torch.cuda.get_device_name(devices[0]), devices=len(devices), mode=mode,
                               tiles=tiles, traits=traits, samples=args.samples, chunk=args.chunk,
                               chunks=args.chunks, seconds=seconds, times=times,
                               variants_per_second=args.chunks * args.chunk / seconds, bounds=bounds,
                               seconds_per_column=costs[0], seconds_per_area=costs[1])
                    print(json.dumps(row), flush=True)
                    handle.write(json.dumps(row) + '\n')
                    del state
                    torch.cuda.empty_cache()


if __name__ == '__main__':
    main()
