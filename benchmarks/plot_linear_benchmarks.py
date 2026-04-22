from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _device_color(device: str) -> str:
    return {"cpu": "#1f4e79", "cuda": "#c45a1a"}.get(device, "#666666")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results/benchmarks")
    parser.add_argument("--output-dir", default="results/figures")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    small = pd.read_table(results_dir / "linear_memory_profile.tsv")
    large = pd.read_table(results_dir / "linear_streaming_profile_large.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    for device, subdf in small.groupby("device"):
        x = subdf["n_samples"].astype(str) + "x" + subdf["n_markers"].astype(str) + "x" + subdf["n_traits"].astype(str)
        ax.plot(x, subdf["cpu_peak_rss_mb"], marker="o", linewidth=2, label=device.upper(), color=_device_color(device))
    ax.set_title("Small/Medium Linear GWAS CPU Peak RSS")
    ax.set_ylabel("Peak RSS (MB)")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    for device, subdf in small.groupby("device"):
        x = subdf["n_samples"].astype(str) + "x" + subdf["n_markers"].astype(str) + "x" + subdf["n_traits"].astype(str)
        ax.plot(x, subdf["runtime_seconds"], marker="o", linewidth=2, label=device.upper(), color=_device_color(device))
    ax.set_title("Small/Medium Linear GWAS Runtime")
    ax.set_ylabel("Runtime (s)")
    ax.tick_params(axis="x", rotation=20)

    ax = axes[1, 0]
    memory_cols = ["cpu_peak_rss_mb", "gpu_peak_allocated_mb", "gpu_peak_reserved_mb"]
    labels = ["CPU RSS", "GPU alloc", "GPU reserved"]
    width = 0.22
    xpos = [0, 1]
    for offset, (col, label) in enumerate(zip(memory_cols, labels)):
        vals = [large.loc[large["device"] == "cpu", col].iloc[0], large.loc[large["device"] == "cuda", col].iloc[0]]
        ax.bar([x + (offset - 1) * width for x in xpos], vals, width=width, label=label)
    ax.set_xticks(xpos)
    ax.set_xticklabels(["CPU run", "CUDA run"])
    ax.set_title("Large Streaming Profile Peak Memory")
    ax.set_ylabel("Memory (MB)")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    runtime_cols = ["avg_chunk_seconds", "estimated_full_runtime_seconds"]
    labels = ["Avg chunk time", "Estimated full runtime"]
    width = 0.28
    for offset, (col, label) in enumerate(zip(runtime_cols, labels)):
        vals = [large.loc[large["device"] == "cpu", col].iloc[0], large.loc[large["device"] == "cuda", col].iloc[0]]
        ax.bar([x + (offset - 0.5) * width for x in xpos], vals, width=width, label=label)
    ax.set_xticks(xpos)
    ax.set_xticklabels(["CPU run", "CUDA run"])
    ax.set_title("Large Streaming Profile Runtime")
    ax.set_ylabel("Seconds")
    ax.legend(frameon=False)

    fig.suptitle("TorchGWAS Linear Benchmark Summary", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_dir / "linear_benchmark_summary.png", dpi=220)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
