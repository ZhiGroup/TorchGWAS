"""Probe: a fusion-friendly exact -log10 P against tails.upper_tail_log10_from_t_torch.

Same identity and continued fraction as the release (tails.py), rearranged so
the compiler can fuse it: no host synchronisation, per-variant terms computed
once per row, 1 - x formed as t^2 / (df + t^2), both branches evaluated and
selected elementwise. Times eager and torch.compile forms on one scan chunk and
checks them against scipy's stdtr where it is representable.

    python benchmarks/logp_fused_probe_20260927.py --device cuda:0
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


def _cf(a, b, x):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = torch.ones_like(x)
    d = 1.0 - qab * x / qap
    d = 1.0 / torch.where(d.abs() < TINY, TINY, d)
    h = d
    for m in range(1, ITERATIONS + 1):
        m2 = 2 * m
        step = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + step * d
        d = 1.0 / torch.where(d.abs() < TINY, TINY, d)
        c = 1.0 + step / c
        c = torch.where(c.abs() < TINY, TINY, c)
        h = h * d * c
        step = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + step * d
        d = 1.0 / torch.where(d.abs() < TINY, TINY, d)
        c = 1.0 + step / c
        c = torch.where(c.abs() < TINY, TINY, c)
        h = h * d * c
    return h


def neg_log10_p(t, df):
    """t (rows, traits) any float dtype; df (rows, 1) or broadcastable."""
    t = t.double().abs()
    df = df.double()
    a = df * 0.5
    b = 0.5
    t2 = t * t
    s = df + t2
    x = df / s
    y = t2 / s
    log_beta = torch.lgamma(a) + math.lgamma(b) - torch.lgamma(a + b)
    reflect = x >= (a + 1.0) / (a + b + 2.0)
    # Direct: log I_x(a, b); the fraction is fed a harmless 0.5 where reflected.
    xd = torch.where(reflect, 0.5, x)
    log_direct = (a * torch.log(x.clamp_min(TINY)) + b * torch.log(y.clamp_min(TINY))
                  - torch.log(a) - log_beta + torch.log(_cf(a, b, xd)))
    # Reflected: I_x(a, b) = 1 - I_y(b, a), O(1) here.
    yr = torch.where(reflect, y, 0.5)
    upper = torch.exp(-log_beta + b * torch.log(yr.clamp_min(TINY)) + a * torch.log1p(-yr)) * _cf(b, a, yr) / b
    log_reflect = torch.log((1.0 - upper).clamp(TINY, 1.0))
    return -torch.where(reflect, log_reflect, log_direct) / LN10


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
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', default='cuda:0')
    args = parser.parse_args()
    device = torch.device(args.device)
    torch.cuda.set_device(device)
    rng = np.random.default_rng(20260927)
    for df_value in (22238.0, 500.0, 30.0):
        grid = np.concatenate([np.linspace(0, 8, 4001), np.linspace(8, 37, 2000)])
        reference = -np.log10(2 * special.stdtr(df_value, -grid))
        t = torch.as_tensor(grid, device=device)[:, None]
        df = torch.full((grid.size, 1), df_value, device=device, dtype=torch.float64)
        value = neg_log10_p(t, df).cpu().numpy()[:, 0]
        print(json.dumps(dict(df=df_value, max_abs_err=float(np.max(np.abs(value - reference))))), flush=True)
    compiled = torch.compile(neg_log10_p, dynamic=True)
    for rows, traits in ((1024, 8192), (4096, 2048), (1024, 512)):
        t = torch.as_tensor(rng.standard_t(22238.0, size=(rows, traits)).astype(np.float32), device=device)
        t.view(-1)[::997] *= 12.0
        df = torch.full((rows, 1), 22238.0, device=device)
        df[::13] -= 3.0
        release = upper_tail_log10_from_t_torch(t, df)
        started = time.perf_counter()
        fused = compiled(t, df)
        torch.cuda.synchronize()
        compile_seconds = time.perf_counter() - started
        torch.cuda.reset_peak_memory_stats(device)
        base = torch.cuda.memory_allocated(device)
        compiled_ms = timed(compiled, t, df)
        peak = torch.cuda.max_memory_allocated(device) - base
        print(json.dumps(dict(shape=[rows, traits], release_ms=timed(upper_tail_log10_from_t_torch, t, df),
                              eager_rearranged_ms=timed(neg_log10_p, t, df), compiled_ms=compiled_ms,
                              compiled_peak_bytes=int(peak), first_call_seconds=compile_seconds,
                              max_abs_vs_release=float((fused - release).abs().max()))), flush=True)


if __name__ == '__main__':
    main()
