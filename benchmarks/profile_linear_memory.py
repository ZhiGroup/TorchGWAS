from __future__ import annotations

import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from torchgwas.api import run_linear_gwas


def _parse_shape(token: str) -> tuple[int, int, int]:
    parts = token.lower().split("x")
    if len(parts) != 3:
        raise ValueError(f"invalid shape token {token}; expected n_samplesxn_markersxn_traits")
    return int(parts[0]), int(parts[1]), int(parts[2])


def _rss_mb() -> float:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return float(value) / 1024.0


def _run_worker(args: argparse.Namespace) -> int:
    rng = np.random.default_rng(args.seed)
    genotype = rng.normal(size=(args.n_samples, args.n_markers)).astype(np.float64)
    phenotype = rng.normal(size=(args.n_samples, args.n_traits)).astype(np.float64)
    covariates = rng.normal(size=(args.n_samples, args.n_covariates)).astype(np.float64)

    if args.device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA worker requested but torch.cuda.is_available() is False")
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()

    started = time.perf_counter()
    result = run_linear_gwas(
        genotype=genotype,
        phenotype=phenotype,
        covariates=covariates,
        chunk_size=args.chunk_size,
        device=args.device,
    )
    if args.device == "cuda":
        torch.cuda.synchronize()
        gpu_peak_allocated_mb = torch.cuda.max_memory_allocated() / (1024.0**2)
        gpu_peak_reserved_mb = torch.cuda.max_memory_reserved() / (1024.0**2)
    else:
        gpu_peak_allocated_mb = 0.0
        gpu_peak_reserved_mb = 0.0
    runtime_seconds = round(time.perf_counter() - started, 6)

    payload = {
        "analysis": "linear",
        "device": args.device,
        "n_samples": args.n_samples,
        "n_markers": args.n_markers,
        "n_traits": args.n_traits,
        "n_covariates": args.n_covariates,
        "chunk_size": args.chunk_size,
        "cpu_peak_rss_mb": round(_rss_mb(), 3),
        "gpu_peak_allocated_mb": round(gpu_peak_allocated_mb, 3),
        "gpu_peak_reserved_mb": round(gpu_peak_reserved_mb, 3),
        "runtime_seconds": runtime_seconds,
        "n_result_rows": len(result.table),
        "provenance": "simulated_profile_linear_memory",
    }
    print(json.dumps(payload))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=None)
    parser.add_argument("--shapes", nargs="+", default=["1024x4096x8", "2048x8192x8"])
    parser.add_argument("--devices", nargs="+", default=["cpu", "cuda"])
    parser.add_argument("--n-covariates", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--n-samples", type=int, default=None)
    parser.add_argument("--n-markers", type=int, default=None)
    parser.add_argument("--n-traits", type=int, default=None)
    args = parser.parse_args()

    if args.worker:
        return _run_worker(args)

    rows: list[dict] = []
    script = Path(__file__).resolve()
    for shape in args.shapes:
        n_samples, n_markers, n_traits = _parse_shape(shape)
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
                str(n_samples),
                "--n-markers",
                str(n_markers),
                "--n-traits",
                str(n_traits),
                "--n-covariates",
                str(args.n_covariates),
                "--chunk-size",
                str(args.chunk_size),
                "--seed",
                str(args.seed),
            ]
            env = dict(__import__("os").environ)
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
