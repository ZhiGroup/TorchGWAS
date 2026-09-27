"""Per-chunk JAGWAS stage times: is computing z from r (skipping t) worth an API change?

CUDA-event times per stage for one chunk on one GPU, full-scale sample count:
  gemm        genotype x [phenotype | covariates] products (common to both paths)
  t_tail      linear._linear_stats after the products: beta, se, t, validity
  score_on_t  reduction prologue now: guards, z = t / sqrt(1 + t^2 / df), FP64 cast
  z_from_r    the alternative: z = sqrt(df) gy / sqrt(ss_g ss_y), guards, FP64 cast
  projection  block-triangular L^-1 z, square, column sum (FP64)
    jagwas_score_path_bench_20260926.py [--device cuda:0]
"""
import argparse
import json

import torch

from torchgwas.jagwas_blocks import triangular_blocks


def timed(fn, repeats=20):
    for _ in range(3):
        fn()
    torch.cuda.synchronize()
    start, end = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
    start.record()
    for _ in range(repeats):
        fn()
    end.record(); torch.cuda.synchronize()
    return start.elapsed_time(end) / repeats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--samples', type=int, default=35365)
    args = parser.parse_args()
    torch.backends.cuda.matmul.allow_tf32 = False
    device = torch.device(args.device)
    n, c = args.samples, 28
    df = float(n - c - 1)
    for chunk in (2048, 4096):
        for k in (512, 2048, 8192):
            g = torch.randn(chunk, n, device=device)
            design = torch.randn(n, k + c, device=device)
            ss_y = (design[:, :k] ** 2).sum(0)
            factor = torch.linalg.cholesky(torch.eye(k, dtype=torch.float64, device=device) * 2
                                           + torch.full((k, k), 0.1, dtype=torch.float64, device=device))
            inverse = torch.linalg.inv(factor)
            blocks = [(s, e, inverse[s:e, :e]) for s, e in triangular_blocks(k)]
            products = g @ design
            gy, gc = products[:, :k], products[:, k:]
            residual_ss = (g * g).sum(1) - (gc * gc).sum(1)
            valid = residual_ss > 1e-12
            variant_df = torch.full_like(residual_ss, df)

            def gemm():
                return g @ design

            def t_tail():
                safe_df = torch.clamp(variant_df, min=1.0); safe_ss = torch.clamp(residual_ss, min=1e-12)
                beta = gy / safe_ss[:, None]
                explained = gy * gy / safe_ss[:, None]
                residual = torch.clamp(ss_y[None, :] - explained, min=1e-12)
                se = torch.sqrt(residual / safe_df[:, None] / safe_ss[:, None])
                t = beta / se
                return torch.where(valid[:, None], beta, torch.zeros_like(beta)), torch.where(valid[:, None], t, torch.zeros_like(t))

            _, t = t_tail()

            def score_on_t():
                scores = torch.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0)
                bad = ~torch.isfinite(t).all(dim=1)
                scores.mul_(scores.square().div_(variant_df.unsqueeze(1)).add_(1.0).rsqrt_())
                return scores.double(), bad

            def z_from_r():
                scale = torch.rsqrt(torch.clamp(residual_ss, min=1e-12)) * df ** 0.5
                z = gy * scale[:, None] * torch.rsqrt(ss_y)[None, :]
                bad = ~torch.isfinite(z).all(dim=1)
                return torch.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0).double(), bad

            z64, _ = score_on_t()

            def projection():
                transposed = z64.T
                out = torch.empty((k, chunk), dtype=torch.float64, device=device)
                for s, e, rows in blocks:
                    torch.mm(rows, transposed[:e], out=out[s:e])
                return out.square_().sum(dim=0)

            times = {name: timed(fn) for name, fn in (('gemm', gemm), ('t_tail', t_tail), ('score_on_t', score_on_t),
                                                      ('z_from_r', z_from_r), ('projection', projection))}
            now = times['gemm'] + times['t_tail'] + times['score_on_t'] + times['projection']
            alternative = times['gemm'] + times['z_from_r'] + times['projection']
            print(json.dumps(dict(gpu=torch.cuda.get_device_name(device), chunk=chunk, traits=k,
                                  ms={key: round(value, 3) for key, value in times.items()},
                                  current_ms=round(now, 3), z_from_r_ms=round(alternative, 3),
                                  saving_fraction=round(1 - alternative / now, 4))), flush=True)
            del g, design, products, gy, gc, t, z64, inverse, factor, blocks
            torch.cuda.empty_cache()


if __name__ == '__main__':
    main()
