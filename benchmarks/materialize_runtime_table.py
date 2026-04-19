from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = [
        {
            "tool": "fastGWA",
            "n_traits": 2048,
            "n_tests": 204800,
            "runtime_seconds": 204800,
            "provenance": "historical_notebook_compare_with_plink.ipynb",
        },
        {
            "tool": "torchGWAS_prototype",
            "n_traits": 2048,
            "n_tests": 204800,
            "runtime_seconds": 940,
            "provenance": "historical_notebook_compare_with_plink.ipynb",
        },
        {
            "tool": "fastGWA",
            "n_traits": 20480,
            "n_tests": 2048000,
            "runtime_seconds": 2048000,
            "provenance": "historical_notebook_compare_with_plink.ipynb",
        },
        {
            "tool": "torchGWAS_prototype",
            "n_traits": 20480,
            "n_tests": 2048000,
            "runtime_seconds": 1108,
            "provenance": "historical_notebook_compare_with_plink.ipynb",
        },
    ]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

