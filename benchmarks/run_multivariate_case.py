from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from torchgwas.api import run_multivariate_gwas
from torchgwas.datasets import get_toy_dataset_paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    toy = get_toy_dataset_paths()
    result = run_multivariate_gwas(
        genotype=toy["genotype"],
        phenotype=toy["phenotype"],
        covariates=toy["covariates"],
        marker_ids=toy["marker_ids"],
        sample_ids=toy["sample_ids"],
        chunk_size=4,
    )
    df = pd.DataFrame(result.table).sort_values("-log10_p", ascending=False).head(10)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

