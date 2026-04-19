from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results/benchmarks")
    parser.add_argument("--output-dir", default="results/figures")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    accuracy = pd.read_table(results_dir / "accuracy_plink.tsv")
    runtime = pd.read_table(results_dir / "runtime_scaling.tsv")
    multi = pd.read_table(results_dir / "multivariate_case.tsv")

    plt.figure(figsize=(6, 4))
    plt.scatter(accuracy["trait_index"], accuracy["pearson_r"], s=16)
    plt.xlabel("Trait index")
    plt.ylabel("Pearson r vs PLINK2")
    plt.ylim(0.0, 1.01)
    plt.tight_layout()
    plt.savefig(output_dir / "figure2_accuracy.png", dpi=200)
    plt.close()

    plt.figure(figsize=(6, 4))
    for tool, subdf in runtime.groupby("tool"):
        plt.plot(subdf["n_traits"], subdf["runtime_seconds"], marker="o", label=tool)
    plt.xlabel("Number of traits")
    plt.ylabel("Runtime (seconds)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "figure3_runtime.png", dpi=200)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.bar(multi["marker_id"], multi["-log10_p"])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("-log10 p")
    plt.tight_layout()
    plt.savefig(output_dir / "figure3_multivariate_case.png", dpi=200)
    plt.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

