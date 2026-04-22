from __future__ import annotations

import argparse
import json
import math
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import stats

from torchgwas.kernels import linear_chunk_kernel
from torchgwas.preprocess import residualize_and_standardize


def _rss_mb() -> float:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return float(value) / 1024.0


def _run_worker(args: argparse.Namespace) -> int:
    rng = np.random.default_rng(args.seed)
    phenotype = rng.normal(size=(args.n_samples, args.n_traits)).astype(np.float64)
    covariates = rng.normal(size=(args.n_samples, args.n_covariates)).astype(np.float64)
    pheno_proc, _ = residualize_and_standardize(phenotype, covariates)

    if args.device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA worker requested but torch.cuda.is_available() is False")
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()

    torch_device = torch.device(args.device)
    pheno_t = torch.as_tensor(pheno_proc, dtype=torch.float64, device=torch_device)
    total_chunks = math.ceil(args.n_markers_total / args.chunk_size)
    timed_chunks = 0
    timed_seconds = 0.0
    max_cpu_rss = _rss_mb()

    for chunk_index in range(args.warmup_chunks + args.profile_chunks):
        current_chunk = min(args.chunk_size, args.n_markers_total)
        genotype_chunk = rng.normal(size=(args.n_samples, current_chunk)).astype(np.float64)
        geno_t = torch.as_tensor(genotype_chunk, dtype=torch.float64, device=torch_device)

        start = time.perf_counter()
        corr_t, t_stat_t = linear_chunk_kernel(geno_t, pheno_t)
        if args.device == "cuda":
            torch.cuda.synchronize()
        t_stat = t_stat_t.detach().cpu().numpy()
        _ = 2.0 * stats.t.sf(np.abs(t_stat), df=args.n_samples - 2)
        elapsed = time.perf_counter() - start

        max_cpu_rss = max(max_cpu_rss, _rss_mb())
        if chunk_index >= args.warmup_chunks:
            timed_chunks += 1
            timed_seconds += elapsed

    if args.device == "cuda":
        gpu_peak_allocated_mb = torch.cuda.max_memory_allocated() / (1024.0**2)
        gpu_peak_reserved_mb = torch.cuda.max_memory_reserved() / (1024.0**2)
    else:
        gpu_peak_allocated_mb = 0.0
        gpu_peak_reserved_mb = 0.0

    avg_chunk_seconds = timed_seconds / max(timed_chunks, 1)
    payload = {
        "analysis": "linear_streaming_profile",
        "device": args.device,
        "n_samples": args.n_samples,
        "n_markers_total": args.n_markers_total,
        "n_traits": args.n_traits,
        "n_covariates": args.n_covariates,
        "chunk_size": args.chunk_size,
        "warmup_chunks": args.warmup_chunks,
        "profile_chunks": args.profile_chunks,
        "cpu_peak_rss_mb": round(max_cpu_rss, 3),
        "gpu_peak_allocated_mb": round(gpu_peak_allocated_mb, 3),
        "gpu_peak_reserved_mb": round(gpu_peak_reserved_mb, 3),
        "avg_chunk_seconds": round(avg_chunk_seconds, 6),
        "estimated_full_runtime_seconds": round(avg_chunk_seconds * total_chunks, 6),
        "provenance": "simulated_streaming_large_scale_profile",
    }
    print(json.dumps(payload))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=None)
    parser.add_argument("--devices", nargs="+", default=["cpu", "cuda"])
    parser.add_argument("--n-samples", type=int, default=20000)
    parser.add_argument("--n-markers-total", type=int, default=1000000)
    parser.add_argument("--n-traits", type=int, default=2000)
    parser.add_argument("--n-covariates", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--warmup-chunks", type=int, default=1)
    parser.add_argument("--profile-chunks", type=int, default=2)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()

    if args.worker:
        return _run_worker(args)

    rows: list[dict] = []
    script = Path(__file__).resolve()
    for device in args.devices:
        if device == "cuda" and not torch.cuda.is_available():
            continue
        cmd = [
            sys.executable,
            str(script),
            "--worker",
            "--device",
            device,
            "--n-samples",
            str(args.n_samples),
            "--n-markers-total",
            str(args.n_markers_total),
            "--n-traits",
            str(args.n_traits),
            "--n-covariates",
            str(args.n_covariates),
            "--chunk-size",
            str(args.chunk_size),
            "--warmup-chunks",
            str(args.warmup_chunks),
            "--profile-chunks",
            str(args.profile_chunks),
            "--seed",
            str(args.seed),
        ]
        env = dict(os.environ)
        root = script.resolve().parents[1]
        env["PYTHONPATH"] = str(root / "src")
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        rows.append(json.loads(completed.stdout.strip().splitlines()[-1]))

    if args.output is not None:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    else:
        print(pd.DataFrame(rows).to_csv(sep="\t", index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
