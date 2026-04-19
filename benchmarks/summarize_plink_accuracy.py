from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_plink_additive(path: Path) -> np.ndarray:
    values = []
    with path.open() as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            if fields[6] == "ADD":
                try:
                    p = float(fields[11])
                    value = -np.log10(max(p, np.finfo(float).tiny))
                except ValueError:
                    value = np.nan
                values.append(value)
    return np.asarray(values, dtype=np.float64)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast-gwas-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--trait-indices", default="0,17,63,127")
    args = parser.parse_args()

    root = Path(args.fast_gwas_root)
    logp = np.load(root / "logp.npy", mmap_mode="r")
    rows = []
    trait_indices = [int(token) for token in args.trait_indices.split(",") if token.strip()]
    for trait_idx in trait_indices:
        if trait_idx >= logp.shape[1]:
            continue
        plink_path = root / f"plink2.QT{trait_idx}.glm.linear"
        if not plink_path.exists():
            continue
        plink = load_plink_additive(plink_path)
        torch_vals = np.asarray(logp[: len(plink), trait_idx], dtype=np.float64)
        valid = np.isfinite(torch_vals) & np.isfinite(plink)
        if valid.sum() == 0:
            continue
        torch_vals = torch_vals[valid]
        plink = plink[valid]
        corr = float(np.corrcoef(torch_vals, plink)[0, 1])
        mae = float(np.mean(np.abs(torch_vals - plink)))
        max_abs = float(np.max(np.abs(torch_vals - plink)))
        rows.append(
            {
                "trait_index": trait_idx,
                "n_markers": int(valid.sum()),
                "pearson_r": corr,
                "mae_log10p": mae,
                "max_abs_diff_log10p": max_abs,
                "source_file": str(plink_path.name),
            }
        )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
