from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    root = Path(__file__).resolve().parents[1] / "examples" / "toy"
    root.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(42)
    n_samples = 64
    n_markers = 12
    n_covars = 4
    n_traits = 3

    covar = rng.normal(size=(n_samples, n_covars))
    genotype = rng.integers(0, 3, size=(n_samples, n_markers)).astype(np.float64)
    signal = genotype[:, [2, 7, 10]].astype(np.float64)
    noise = rng.normal(scale=0.5, size=(n_samples, n_traits))
    pheno = np.column_stack(
        [
            1.2 * signal[:, 0] - 0.5 * covar[:, 0] + noise[:, 0],
            -0.8 * signal[:, 1] + 0.7 * covar[:, 1] + noise[:, 1],
            0.9 * signal[:, 2] + 0.4 * signal[:, 1] - 0.4 * covar[:, 2] + noise[:, 2],
        ]
    )

    np.save(root / "genotype.npy", genotype)
    np.save(root / "covar.npy", covar)
    np.save(root / "pheno.npy", pheno)
    marker_ids = [f"rsToy{i:03d}" for i in range(n_markers)]
    sample_ids = [f"sample_{i:03d}" for i in range(n_samples)]
    (root / "markers.tsv").write_text("\n".join(marker_ids) + "\n")
    (root / "samples.tsv").write_text("\n".join(sample_ids) + "\n")
    pheno_table = pd.DataFrame(pheno, columns=["trait_0", "trait_1", "trait_2"])
    pheno_table.insert(0, "IID", sample_ids)
    pheno_table.insert(0, "FID", sample_ids)
    pheno_table.to_csv(root / "pheno.tsv", sep="\t", index=False)
    covar_table = pd.DataFrame(covar, columns=["covar_0", "covar_1", "covar_2", "covar_3"])
    covar_table.insert(0, "IID", sample_ids)
    covar_table.insert(0, "FID", sample_ids)
    covar_table.to_csv(root / "covar.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
