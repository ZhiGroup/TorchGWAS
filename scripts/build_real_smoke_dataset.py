from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genotype-npy", default="/data4012/zxie3/fast_GWAS/genotype.npy")
    parser.add_argument("--marker-npy", default="/data4012/zxie3/fast_GWAS/markers.npy")
    parser.add_argument("--sample-list", default="/data4012/zxie3/fast_GWAS/subset")
    parser.add_argument("--pheno-table", default="real_case/t1_128/t1_pheno.tsv")
    parser.add_argument("--covar-table", default="real_case/t1_128/t1_covar.tsv")
    parser.add_argument("--output-dir", default="real_case/t1_128_smoke")
    parser.add_argument("--n-samples", type=int, default=512)
    parser.add_argument("--n-markers", type=int, default=1024)
    parser.add_argument("--n-traits", type=int, default=8)
    parser.add_argument("--use-source-markers", action="store_true")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    sample_ids = [line.strip() for line in Path(args.sample_list).read_text().splitlines() if line.strip()]
    pheno = pd.read_table(args.pheno_table)
    covar = pd.read_table(args.covar_table)
    pheno["IID"] = pheno["IID"].astype(str)
    covar["IID"] = covar["IID"].astype(str)
    pheno_ids = set(pheno["IID"])
    covar_ids = set(covar["IID"])
    common = [iid for iid in sample_ids if iid in pheno_ids and iid in covar_ids]
    selected = common[: args.n_samples]
    if len(selected) < args.n_samples:
        raise ValueError(f"only found {len(selected)} overlapping samples, fewer than requested {args.n_samples}")

    genotype_mmap = np.load(args.genotype_npy, mmap_mode="r")
    genotype = np.asarray(genotype_mmap[: args.n_samples, : args.n_markers], dtype=np.float64)

    if args.use_source_markers:
        markers = np.load(args.marker_npy, allow_pickle=True)[: args.n_markers]
        marker_names = [str(v) for v in markers.tolist()]
    else:
        marker_names = [f"real_marker_{i}" for i in range(args.n_markers)]
    pheno = pheno.drop_duplicates(subset=["IID"]).set_index("IID")
    covar = covar.drop_duplicates(subset=["IID"]).set_index("IID")
    trait_columns = [col for col in pheno.columns if col.startswith("QT")][: args.n_traits]
    covariate_columns = [col for col in covar.columns if col not in {"FID", "IID"}]

    pheno_out = pheno.loc[selected, ["FID", *trait_columns]].copy()
    pheno_out.insert(0, "IID", pheno_out.index.astype(str))
    pheno_out = pheno_out[["FID", "IID", *trait_columns]].reset_index(drop=True)

    covar_out = covar.loc[selected, ["FID", *covariate_columns]].copy()
    covar_out.insert(0, "IID", covar_out.index.astype(str))
    covar_out = covar_out[["FID", "IID", *covariate_columns]].reset_index(drop=True)

    np.save(out / "genotype.npy", genotype)
    pheno_out.to_csv(out / "pheno.tsv", sep="\t", index=False)
    covar_out.to_csv(out / "covar.tsv", sep="\t", index=False)
    (out / "markers.tsv").write_text("\n".join(marker_names) + "\n")
    (out / "samples.tsv").write_text("\n".join(selected) + "\n")
    (out / "trait_columns.txt").write_text("\n".join(trait_columns) + "\n")
    (out / "covariate_columns.txt").write_text("\n".join(covariate_columns) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "n_samples": len(selected),
                "n_markers": int(args.n_markers),
                "n_traits": len(trait_columns),
                "n_covariates": len(covariate_columns),
                "source_genotype": args.genotype_npy,
                "source_pheno": args.pheno_table,
                "source_covar": args.covar_table,
            },
            indent=2,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
