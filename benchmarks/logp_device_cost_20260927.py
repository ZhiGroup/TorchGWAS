"""What exact -log10 P costs on the device, per scan chunk.

The public release computes it with tails.upper_tail_log10_from_t_torch: 40
continued-fraction iterations of eager FP64 elementwise operations, each
materialising a chunk-sized temporary. This times that function, and the same
function under torch.compile (one fused kernel), on scan-shaped inputs: FP32 t
and one df per variant, as the native scan produces them.

    python benchmarks/logp_device_cost_20260927.py --device cuda:0 --out results/logp_device_cost_20260927.json
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch

from torchgwas.tails import upper_tail_log10_from_t, upper_tail_log10_from_t_torch


def timed(function, *args, repeats=5):
    function(*args)
    torch.cuda.synchronize()
    times = []
    for _ in range(repeats):
        started = time.perf_counter()
        function(*args)
        torch.cuda.synchronize()
        times.append(time.perf_counter() - started)
    return min(times), sorted(times)[len(times) // 2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--shapes', nargs='+', default=['1024x8192', '4096x2048', '1024x512'])
    parser.add_argument('--df', type=float, default=22238.0)
    args = parser.parse_args()
    device = torch.device(args.device)
    torch.cuda.set_device(device)
    rng = np.random.default_rng(20260927)
    compiled = torch.compile(upper_tail_log10_from_t_torch, dynamic=False)
    rows = []
    for shape in args.shapes:
        rows_, traits = (int(v) for v in shape.split('x'))
        t = rng.standard_t(args.df, size=(rows_, traits)).astype(np.float32)
        t.flat[::997] *= 12.0  # a sprinkle of strong associations
        df = np.full((rows_, 1), args.df, dtype=np.float32)
        df[::13] -= 3.0  # per-variant df, as missing calls produce
        t_d = torch.as_tensor(t, device=device)
        df_d = torch.as_tensor(df, device=device)
        torch.cuda.reset_peak_memory_stats(device)
        base = torch.cuda.memory_allocated(device)
        eager_min, eager_median = timed(upper_tail_log10_from_t_torch, t_d, df_d)
        eager_peak = torch.cuda.max_memory_allocated(device) - base
        started = time.perf_counter()
        compiled(t_d, df_d)
        torch.cuda.synchronize()
        compile_seconds = time.perf_counter() - started
        torch.cuda.reset_peak_memory_stats(device)
        fused_min, fused_median = timed(compiled, t_d, df_d)
        fused_peak = torch.cuda.max_memory_allocated(device) - base
        eager = upper_tail_log10_from_t_torch(t_d, df_d).cpu().numpy()
        fused = compiled(t_d, df_d).cpu().numpy()
        sample = rng.choice(t.size, size=min(t.size, 200_000), replace=False)
        reference = upper_tail_log10_from_t(t.reshape(-1)[sample], np.broadcast_to(df, t.shape).reshape(-1)[sample])
        row = dict(shape=[rows_, traits], cells=rows_ * traits,
                   eager_ms=1e3 * eager_min, eager_median_ms=1e3 * eager_median, eager_peak_bytes=int(eager_peak),
                   compiled_ms=1e3 * fused_min, compiled_median_ms=1e3 * fused_median,
                   compiled_peak_bytes=int(fused_peak), compile_seconds=compile_seconds,
                   eager_vs_compiled_max_abs=float(np.nanmax(np.abs(eager - fused))),
                   eager_vs_scipy_max_rel=float(np.nanmax(np.abs(eager.reshape(-1)[sample] - reference)
                                                          / np.maximum(np.abs(reference), 1e-12))),
                   stored_float32_max_rel=float(np.nanmax(np.abs(eager.reshape(-1)[sample].astype(np.float32)
                                                                 - reference) / np.maximum(np.abs(reference), 1e-12))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    record = dict(device=torch.cuda.get_device_name(device), torch=torch.__version__, df=args.df, rows=rows,
                  scope='Device time of exact -log10 P for one scan chunk (FP32 t, per-variant df), eager and '
                        'torch.compile, with accuracy against the float64 numpy reference on a sample.')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')


if __name__ == '__main__':
    main()
