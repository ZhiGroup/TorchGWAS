"""Cost of CompleteCasePlan.correct per chunk, and of its gather, on one GPU.

The full-scale missing panel's shape: 22,250 samples, K = 512, 64 traits
each missing ~445 samples (their own). Times correct() on a synthetic
(4096, n) chunk with CUDA events, and the gather alone in both
orientations: columns of the (chunk, n) calls, or rows of their transpose.

    python benchmarks/complete_case_gather_probe_20260927.py cuda:4
"""
import json
import sys

import numpy as np
import torch

from torchgwas.complete_case import CompleteCasePlan


def timed(fn, repeats=5):
    fn()
    torch.cuda.synchronize()
    start, stop = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
    start.record()
    for _ in range(repeats):
        fn()
    stop.record()
    stop.synchronize()
    return start.elapsed_time(stop) / repeats


def main():
    device = torch.device(sys.argv[1])
    torch.cuda.set_device(device)
    rng = np.random.default_rng(0)
    n, k, chunk, rank = 22_250, 512, 4096, 10
    observed = np.ones((n, k), dtype=bool)
    for trait in range(64):
        observed[rng.choice(n, 445, replace=False), trait] = False
    q, _ = np.linalg.qr(rng.normal(size=(n, rank)) - 0.0)
    q -= q.mean(0)
    q, _ = np.linalg.qr(q)
    plan = CompleteCasePlan(observed, q)
    centered = torch.randn((chunk, n), device=device)
    rows = torch.as_tensor(np.concatenate(plan.rows), device=device)
    zg = torch.randn((chunk, rank + 1), device=device)
    products = torch.randn((chunk, k), device=device)
    beta, t = torch.randn((chunk, k), device=device), torch.randn((chunk, k), device=device)
    ss = torch.full((k,), float(n), device=device)
    variant_df = torch.full((chunk,), float(n - rank - 2), device=device)
    sums = (centered * centered).sum(1)
    result = dict(
        missing_cells=int(rows.numel()),
        correct_ms=timed(lambda: plan.correct(centered, sums, zg, products, ss, variant_df, beta, t)),
        gather_columns_ms=timed(lambda: centered.index_select(1, rows)),
        transpose_ms=timed(lambda: centered.t().contiguous()),
    )
    transposed = centered.t().contiguous()
    result['gather_rows_of_transpose_ms'] = timed(lambda: transposed.index_select(0, rows))
    print(json.dumps({key: round(value, 3) if isinstance(value, float) else value for key, value in result.items()}))


if __name__ == '__main__' and len(sys.argv) == 2:
    main()


def profile(device_name):
    """Top CUDA kernels of one correct() call (torch.profiler)."""
    import torch.profiler as tp
    device = torch.device(device_name)
    torch.cuda.set_device(device)
    rng = np.random.default_rng(0)
    n, k, chunk, rank = 22_250, 512, 4096, 10
    observed = np.ones((n, k), dtype=bool)
    for trait in range(64):
        observed[rng.choice(n, 445, replace=False), trait] = False
    q, _ = np.linalg.qr(rng.normal(size=(n, rank)))
    q -= q.mean(0)
    q, _ = np.linalg.qr(q)
    plan = CompleteCasePlan(observed, q)
    centered = torch.randn((chunk, n), device=device)
    args = (centered, (centered * centered).sum(1), torch.randn((chunk, rank + 1), device=device),
            torch.randn((chunk, k), device=device), torch.full((k,), float(n), device=device),
            torch.full((chunk,), float(n - rank - 2), device=device),
            torch.randn((chunk, k), device=device), torch.randn((chunk, k), device=device))
    plan.correct(*args)
    torch.cuda.synchronize()
    with tp.profile(activities=[tp.ProfilerActivity.CUDA]) as prof:
        plan.correct(*args)
        torch.cuda.synchronize()
    print(prof.key_averages().table(sort_by='cuda_time_total', row_limit=14))


if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'profile':
    profile(sys.argv[1])
