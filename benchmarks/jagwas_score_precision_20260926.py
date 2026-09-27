"""Precision of the scan's score z (FP32 GEMM, t formula, score transform) on a GPU.

The JAGWAS cutoff keeps the longest pivot-ordered trait set whose null
rounding error in T, 2 eps_z sqrt(tr R_S^-1), stays within a target; eps_z is
the rms rounding of each z (unit null variance). This measures it: the scan's
FP32 path (linear._linear_stats on centred dosages and the processed panel,
then z = t / sqrt(1 + t^2 / df)) against the same arithmetic in FP64 on the
same FP32 inputs. Null variants give eps_z; injected strong effects give the
relative rounding at large |z|. sqrt(df) r has error ~ u * O(1) independent
of N (|g||y| ~ N), which the two sample counts check.
    jagwas_score_precision_20260926.py [--device cuda:0]
"""
import argparse
import json

import numpy as np
import torch


def statistics(g, design, k, df, dtype):
    g = g.to(dtype); design = design.to(dtype)
    products = g @ design
    gy, gc = products[:, :k], products[:, k:]
    ss = torch.clamp((g * g).sum(1) - (gc * gc).sum(1), min=1e-12)
    beta = gy / ss[:, None]
    explained = gy * gy / ss[:, None]
    phenotype_ss = (design[:, :k] * design[:, :k]).sum(0)
    residual = torch.clamp(phenotype_ss[None, :] - explained, min=1e-12)
    t = beta / torch.sqrt(residual / df / ss[:, None])
    return t * torch.rsqrt(t.square() / df + 1.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', default='cuda:0')
    args = parser.parse_args()
    torch.backends.cuda.matmul.allow_tf32 = False
    device = torch.device(args.device)
    rng = np.random.default_rng(20260926)
    for n in (8000, 35365):
        k, c, chunk = 512, 27, 4096
        covariates = np.column_stack([np.ones(n), rng.standard_normal((n, c - 1))])
        q, _ = np.linalg.qr(covariates)
        panel = rng.standard_normal((n, k)) @ rng.standard_normal((k, k)) * 0.1 + rng.standard_normal((n, k))
        panel -= q @ (q.T @ panel)
        panel /= panel.std(0)
        design = torch.as_tensor(np.column_stack([panel, q]), dtype=torch.float32, device=device)
        df = float(n - c - 1)
        rows = {}
        for effect in (0.0, 0.05, 0.15):
            maf = rng.uniform(0.05, 0.5, chunk)
            dosage = rng.binomial(2, maf, (n, chunk)).astype(np.float64)
            dosage += effect * panel[:, :1] / 0.5  # a real effect on trait 0 when effect > 0
            dosage -= dosage.mean(0)
            g = torch.as_tensor(dosage.T, dtype=torch.float32, device=device)
            z32 = statistics(g, design, k, df, torch.float32).double()
            z64 = statistics(g, design, k, df, torch.float64)
            error = (z32 - z64).cpu().numpy()
            scale = np.maximum(np.abs(z64.cpu().numpy()), 1.0)
            rows[effect] = dict(max_abs_z=float(np.abs(z64.cpu().numpy()).max()),
                                rms_abs=float(np.sqrt(np.mean(error ** 2))),
                                p999_abs=float(np.quantile(np.abs(error), 0.999)),
                                rms_relative=float(np.sqrt(np.mean((error / scale) ** 2))),
                                max_relative=float(np.abs(error / scale).max()))
        print(json.dumps(dict(gpu=torch.cuda.get_device_name(device), samples=n, traits=k,
                              unit_roundoff=2.0 ** -24, rows=rows)), flush=True)


if __name__ == '__main__':
    main()
