from __future__ import annotations

import argparse
import json
from pathlib import Path

from .api import run_linear_gwas, run_multivariate_gwas
from .datasets import get_toy_dataset_paths
from .io import DiskBackedGenotype, align_table_to_samples, load_array, load_genotype
from .preprocess import prepare_inputs_for_prep, residualize_and_standardize
from .utils import mkdir, write_json


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="torchgwas")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prep", help="Validate and preprocess phenotype/covariates")
    prep.add_argument("--genotype", required=True)
    prep.add_argument("--genotype-format", default="auto", choices=["auto", "npy", "plink", "bgen"])
    prep.add_argument("--phenotype", default=None)
    prep.add_argument("--phenotype-table", default=None)
    prep.add_argument("--covariates", default=None)
    prep.add_argument("--covariates-table", default=None)
    prep.add_argument("--trait-columns", default=None)
    prep.add_argument("--covariate-columns", default=None)
    prep.add_argument("--sample-id-column", default="IID")
    prep.add_argument("--sample-ids", default=None)
    prep.add_argument("--sample-file", default=None)
    prep.add_argument("--genotype-cache-dir", default=None)
    prep.add_argument("--bim", default=None)
    prep.add_argument("--fam", default=None)
    prep.add_argument("--plink2-binary", default=None)
    prep.add_argument("--output-dir", required=True)

    linear = subparsers.add_parser("linear", help="Run linear GWAS")
    linear.add_argument("--genotype", required=True)
    linear.add_argument("--genotype-format", default="auto", choices=["auto", "npy", "plink", "bgen"])
    linear.add_argument("--phenotype", default=None)
    linear.add_argument("--phenotype-table", default=None)
    linear.add_argument("--covariates", default=None)
    linear.add_argument("--covariates-table", default=None)
    linear.add_argument("--trait-columns", default=None)
    linear.add_argument("--covariate-columns", default=None)
    linear.add_argument("--sample-id-column", default="IID")
    linear.add_argument("--marker-ids", default=None)
    linear.add_argument("--sample-ids", default=None)
    linear.add_argument("--sample-file", default=None)
    linear.add_argument("--genotype-cache-dir", default=None)
    linear.add_argument("--bim", default=None)
    linear.add_argument("--fam", default=None)
    linear.add_argument("--plink2-binary", default=None)
    linear.add_argument("--chunk-size", type=int, default=None)
    linear.add_argument("--output-dir", required=True)

    multi = subparsers.add_parser("multi", help="Run multivariate GWAS")
    multi.add_argument("--genotype", required=True)
    multi.add_argument("--genotype-format", default="auto", choices=["auto", "npy", "plink", "bgen"])
    multi.add_argument("--phenotype", default=None)
    multi.add_argument("--phenotype-table", default=None)
    multi.add_argument("--covariates", default=None)
    multi.add_argument("--covariates-table", default=None)
    multi.add_argument("--trait-columns", default=None)
    multi.add_argument("--covariate-columns", default=None)
    multi.add_argument("--sample-id-column", default="IID")
    multi.add_argument("--marker-ids", default=None)
    multi.add_argument("--sample-ids", default=None)
    multi.add_argument("--sample-file", default=None)
    multi.add_argument("--genotype-cache-dir", default=None)
    multi.add_argument("--bim", default=None)
    multi.add_argument("--fam", default=None)
    multi.add_argument("--plink2-binary", default=None)
    multi.add_argument("--chunk-size", type=int, default=None)
    multi.add_argument("--ridge", type=float, default=1e-6)
    multi.add_argument("--output-dir", required=True)

    demo = subparsers.add_parser("demo", help="Run the bundled toy example")
    demo.add_argument("--output-dir", required=True)
    return parser


def _run_prep(args) -> int:
    genotype, sample_ids, _, genotype_meta = load_genotype(
        args.genotype,
        genotype_format=args.genotype_format,
        bim=args.bim,
        fam=args.fam,
        sample_file=args.sample_file,
        genotype_cache_dir=args.genotype_cache_dir,
        plink2_binary=args.plink2_binary,
    )
    if sample_ids is None and args.sample_ids is not None:
        sample_ids = load_array(args.sample_ids) if str(args.sample_ids).endswith(".npy") else None
        if sample_ids is None:
            from .io import load_vector

            sample_ids = load_vector(args.sample_ids)
    trait_columns = None if args.trait_columns is None else [token.strip() for token in args.trait_columns.split(",") if token.strip()]
    covariate_columns = None if args.covariate_columns is None else [token.strip() for token in args.covariate_columns.split(",") if token.strip()]
    if args.phenotype_table is not None:
        phenotype, trait_columns = align_table_to_samples(
            args.phenotype_table, sample_ids=sample_ids, value_columns=trait_columns, sample_id_column=args.sample_id_column
        )
    else:
        phenotype = load_array(args.phenotype)
    if args.covariates_table is not None:
        covariates, covariate_columns = align_table_to_samples(
            args.covariates_table, sample_ids=sample_ids, value_columns=covariate_columns, sample_id_column=args.sample_id_column
        )
    else:
        covariates = None if args.covariates is None else load_array(args.covariates)
    genotype_array = genotype.genotype if isinstance(genotype, DiskBackedGenotype) else genotype
    phenotype, covariates, qc = prepare_inputs_for_prep(genotype_array, phenotype, covariates)
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    out = mkdir(args.output_dir)
    import numpy as np

    np.save(out / "phenotype_processed.npy", pheno_proc)
    if q_matrix is not None:
        np.save(out / "covariate_q.npy", q_matrix)
    if covariates is not None:
        np.save(out / "covariates_aligned.npy", covariates)
    if sample_ids is not None:
        (out / "samples.tsv").write_text("\n".join([str(v) for v in sample_ids]) + "\n")
    write_json(qc, out / "qc.json")
    write_json(
        {
            **genotype_meta,
            "sample_id_column": args.sample_id_column,
            "trait_columns": trait_columns,
            "covariate_columns": covariate_columns,
        },
        out / "prep.json",
    )
    return 0


def _run_linear(args) -> int:
    trait_columns = None if args.trait_columns is None else [token.strip() for token in args.trait_columns.split(",") if token.strip()]
    covariate_columns = None if args.covariate_columns is None else [token.strip() for token in args.covariate_columns.split(",") if token.strip()]
    run_linear_gwas(
        genotype=args.genotype,
        phenotype=args.phenotype,
        covariates=args.covariates,
        phenotype_table=args.phenotype_table,
        covariates_table=args.covariates_table,
        trait_columns=trait_columns,
        covariate_columns=covariate_columns,
        genotype_format=args.genotype_format,
        sample_file=args.sample_file,
        genotype_cache_dir=args.genotype_cache_dir,
        bim=args.bim,
        fam=args.fam,
        plink2_binary=args.plink2_binary,
        sample_id_column=args.sample_id_column,
        marker_ids=args.marker_ids,
        sample_ids=args.sample_ids,
        chunk_size=args.chunk_size,
        output_dir=args.output_dir,
    )
    return 0


def _run_multi(args) -> int:
    trait_columns = None if args.trait_columns is None else [token.strip() for token in args.trait_columns.split(",") if token.strip()]
    covariate_columns = None if args.covariate_columns is None else [token.strip() for token in args.covariate_columns.split(",") if token.strip()]
    run_multivariate_gwas(
        genotype=args.genotype,
        phenotype=args.phenotype,
        covariates=args.covariates,
        phenotype_table=args.phenotype_table,
        covariates_table=args.covariates_table,
        trait_columns=trait_columns,
        covariate_columns=covariate_columns,
        genotype_format=args.genotype_format,
        sample_file=args.sample_file,
        genotype_cache_dir=args.genotype_cache_dir,
        bim=args.bim,
        fam=args.fam,
        plink2_binary=args.plink2_binary,
        sample_id_column=args.sample_id_column,
        marker_ids=args.marker_ids,
        sample_ids=args.sample_ids,
        chunk_size=args.chunk_size,
        output_dir=args.output_dir,
        ridge=args.ridge,
    )
    return 0


def _run_demo(args) -> int:
    toy = get_toy_dataset_paths()
    out = mkdir(args.output_dir)
    linear_out = out / "linear"
    multi_out = out / "multi"
    linear_result = run_linear_gwas(
        genotype=toy["genotype"],
        phenotype=toy["phenotype"],
        covariates=toy["covariates"],
        marker_ids=toy["marker_ids"],
        sample_ids=toy["sample_ids"],
        chunk_size=4,
        output_dir=linear_out,
    )
    multi_result = run_multivariate_gwas(
        genotype=toy["genotype"],
        phenotype=toy["phenotype"],
        covariates=toy["covariates"],
        marker_ids=toy["marker_ids"],
        sample_ids=toy["sample_ids"],
        chunk_size=4,
        output_dir=multi_out,
    )
    summary = {
        "linear_rows": len(linear_result.table),
        "multi_rows": len(multi_result.table),
        "linear_top_hit": max(linear_result.table, key=lambda row: row["-log10_p"]),
        "multi_top_hit": max(multi_result.table, key=lambda row: row["-log10_p"]),
    }
    write_json(summary, out / "run_summary.json")
    return 0


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()
    if args.command == "prep":
        return _run_prep(args)
    if args.command == "linear":
        return _run_linear(args)
    if args.command == "multi":
        return _run_multi(args)
    if args.command == "demo":
        return _run_demo(args)
    raise ValueError(f"unknown command: {args.command}")
