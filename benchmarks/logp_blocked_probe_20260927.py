"""Probe: one continued fraction per cell, compiled in blocks of iterations.

Each cell needs either the direct fraction K(a, b, x) or the reflected
K(b, a, 1 - x), never both, so the parameters are swapped per cell and one
fraction is evaluated. The 40 iterations are compiled as a block of BLOCK
iterations reused 40 / BLOCK times, which keeps the traced graph (and so the
per-process compile) small while each call is still one fused kernel.

    python benchmarks/logp_blocked_probe_20260927.py --device cuda:0 --block 8
"""
import argparse
import json
import math
import time

import numpy as np
import torch
from scipy import special

from torchgwas.tails import upper_tail_log10_from_t_torch

TINY = 1e-300
ITERATIONS = 40
LN10 = math.log(10.0)


def _guard(v):
    return torch.where(v.abs() < TINY, TINY, v)


def prologue(t, df):
    t = t.double().abs()
    df = df.double()
    a = df * 0.5
    t2 = t * t
    s = df + t2
    x = df / s
    y = t2 / s
    reflect = x >= (a + 1.0) / (a + 2.5)
    # Per cell: the direct fraction K(a, 1/2, x) or the reflected K(1/2, a, y).
    A = torch.where(reflect, 0.5, a)
    B = torch.where(reflect, a, 0.5)
    X = torch.where(reflect, y, x)
    d = 1.0 / _guard(1.0 - (A + B) * X / (A + 1.0))
    return A, B, X, torch.ones_like(X), d, d.clone(), x, y, reflect


def block(A, B, X, c, d, h, first):
    qab, qap, qam = A + B, A + 1.0, A - 1.0
    for offset in range(BLOCK):
        m = first + offset
        m2 = 2.0 * m
        step = m * (B - m) * X / ((qam + m2) * (A + m2))
        d = 1.0 / _guard(1.0 + step * d)
        c = _guard(1.0 + step / c)
        h = h * d * c
        step = -(A + m) * (qab + m) * X / ((A + m2) * (qap + m2))
        d = 1.0 / _guard(1.0 + step * d)
        c = _guard(1.0 + step / c)
        h = h * d * c
    return c, d, h


def epilogue(df, x, y, reflect, h):
    a = df.double() * 0.5
    log_beta = torch.lgamma(a) + math.lgamma(0.5) - torch.lgamma(a + 0.5)
    log_direct = (a * torch.log(x.clamp_min(TINY)) + 0.5 * torch.log(y.clamp_min(TINY))
                  - torch.log(a) - log_beta + torch.log(h))
    upper = torch.exp(-log_beta + 0.5 * torch.log(y.clamp_min(TINY)) + a * torch.log1p(-y)) * h / 0.5
    log_reflect = torch.log((1.0 - upper).clamp(TINY, 1.0))
    return -torch.where(reflect, log_reflect, log_direct) / LN10


def neg_log10_p(t, df, stages=None):
    pro, blk, epi = stages or (prologue, block, epilogue)
    A, B, X, c, d, h, x, y, reflect = pro(t, df)
    for first in range(1, ITERATIONS + 1, BLOCK):
        c, d, h = blk(A, B, X, c, d, h, torch.tensor(float(first), dtype=torch.float64, device=X.device))
    return epi(df, x, y, reflect, h)


def timed(function, *args, repeats=5):
    function(*args)
    torch.cuda.synchronize()
    times = []
    for _ in range(repeats):
        started = time.perf_counter()
        function(*args)
        torch.cuda.synchronize()
        times.append(time.perf_counter() - started)
    return 1e3 * min(times)


def main():
    global BLOCK
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--block', type=int, default=8)
    args = parser.parse_args()
    BLOCK = args.block
    device = torch.device(args.device)
    torch.cuda.set_device(device)
    for df_value in (22238.0, 500.0, 30.0, 3.0):
        grid = np.concatenate([np.linspace(0, 8, 4001), np.linspace(8, 37, 2000)])
        reference = -np.log10(2 * special.stdtr(df_value, -grid))
        finite = np.isfinite(reference)
        value = neg_log10_p(torch.as_tensor(grid, device=device)[:, None],
                            torch.full((grid.size, 1), df_value, device=device)).cpu().numpy()[:, 0]
        print(json.dumps(dict(df=df_value, max_abs_err=float(np.max(np.abs(value - reference)[finite])))), flush=True)
    started = time.perf_counter()
    stages = tuple(torch.compile(f, dynamic=True) for f in (prologue, block, epilogue))
    rng = np.random.default_rng(20260927)
    for rows, traits in ((1024, 512), (1024, 8192), (4096, 2048)):
        t = torch.as_tensor(rng.standard_t(22238.0, size=(rows, traits)).astype(np.float32), device=device)
        t.view(-1)[::997] *= 12.0
        df = torch.full((rows, 1), 22238.0, device=device)
        df[::13] -= 3.0
        first = time.perf_counter()
        value = neg_log10_p(t, df, stages)
        torch.cuda.synchronize()
        first_call = time.perf_counter() - first
        torch.cuda.reset_peak_memory_stats(device)
        base = torch.cuda.memory_allocated(device)
        ms = timed(neg_log10_p, t, df, stages)
        peak = torch.cuda.max_memory_allocated(device) - base
        release = upper_tail_log10_from_t_torch(t, df)
        print(json.dumps(dict(block=BLOCK, shape=[rows, traits], compiled_ms=ms, peak_bytes=int(peak),
                              first_call_seconds=first_call, since_compile_start=time.perf_counter() - started,
                              max_abs_vs_release=float((value - release).abs().max()))), flush=True)


if __name__ == '__main__':
    main()
