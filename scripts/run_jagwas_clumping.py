#!/usr/bin/env python3
"""Run current TorchGWAS JAGWAS reduction and FUMA-style locus clumping.

The scan writes one chi-square statistic per variant with K degrees of freedom.
This wrapper converts that indexed output to the harmonized table expected by
the lab's validated local-clumping implementation, then runs locus clumping.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import special

from torchgwas.api import run_linear_gwas
from torchgwas.sumstats_indexed import open_indexed_sumstats


MHC_START = 29_614_758
MHC_END = 33_170_276
HARMONIZED_FORMAT = "fuma_clump_harmonized_v1"


def parse_variant_range(value: str | None) -> tuple[int, int] | None:
    if value is None:
        return None
    try:
        start_text, end_text = value.split(":", 1)
        start, end = int(start_text), int(end_text)
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError("variant range must be START:END") from error
    if start < 0 or end <= start:
        raise argparse.ArgumentTypeError("variant range must satisfy 0 <= START < END")
    return start, end


def load_maf(path: Path, variant_range: tuple[int, int] | None) -> np.ndarray:
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            if "maf" not in archive.files:
                raise ValueError(f"{path} has no 'maf' array")
            maf = np.asarray(archive["maf"], dtype=np.float32)
    else:
        maf = np.asarray(np.load(path, mmap_mode="r", allow_pickle=False), dtype=np.float32)
    if variant_range is not None:
        maf = maf[slice(*variant_range)]
    return maf


def collect_jagwas(sumstats_dir: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    manifest, parts = open_indexed_sumstats(sumstats_dir)
    if manifest.get("kind") != "jagwas":
        raise ValueError(f"expected jagwas output, got {manifest.get('kind')!r}")
    indices = []
    statistics = []
    for part in parts:
        indices.append(np.asarray(part["variant_index"], dtype=np.int64))
        statistics.append(np.asarray(part["chi2"], dtype=np.float64))
    if not indices:
        return manifest, np.empty(0, np.int64), np.empty(0, np.float64)
    return manifest, np.concatenate(indices), np.concatenate(statistics)


def harmonized_table(
    sumstats_dir: Path,
    maf_path: Path,
    *,
    gwas_p: float,
    maf_min: float,
    exclude_mhc: bool,
    variant_range: tuple[int, int] | None,
) -> tuple[pd.DataFrame, dict]:
    manifest, variant_index, chi2 = collect_jagwas(sumstats_dir)
    degrees_of_freedom = int(manifest["df"])
    if degrees_of_freedom < 1:
        raise ValueError("JAGWAS degrees of freedom must be positive")
    metadata_path = sumstats_dir / "variant_metadata.npz"
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)
    with np.load(metadata_path, allow_pickle=False) as metadata:
        required = {"chromosome", "position", "effect_allele", "other_allele"}
        missing = required.difference(metadata.files)
        if missing:
            raise ValueError(f"variant metadata is missing {sorted(missing)}")
        chromosome = np.asarray(metadata["chromosome"], dtype=str)
        position = np.asarray(metadata["position"], dtype=np.int64)
        effect = np.asarray(metadata["effect_allele"], dtype=str)
        other = np.asarray(metadata["other_allele"], dtype=str)
    marker_ids = np.load(sumstats_dir / "variant_ids.npy", allow_pickle=False).astype(str)
    maf = load_maf(maf_path, variant_range)
    n_variants = len(marker_ids)
    metadata_arrays = (chromosome, position, effect, other)
    if any(len(array) != n_variants for array in metadata_arrays):
        if variant_range is None:
            raise ValueError("marker and variant-metadata arrays are not aligned")
        start, end = variant_range
        if not all(len(array) >= end for array in metadata_arrays):
            raise ValueError("full variant metadata does not cover the requested range")
        chromosome = chromosome[start:end]
        position = position[start:end]
        effect = effect[start:end]
        other = other[start:end]
    if not all(len(array) == n_variants for array in (chromosome, position, effect, other, maf)):
        raise ValueError("marker, metadata, and MAF arrays are not aligned")
    if len(variant_index) and (variant_index.min() < 0 or variant_index.max() >= n_variants):
        raise ValueError("indexed JAGWAS output refers outside the marker table")

    # P(chi2_K >= T) = Q(K/2, T/2).  gammaincc is the regularized upper
    # incomplete gamma.  Clamp only exact floating-point underflow so the
    # strongest valid statistics are retained by the clumper's P > 0 filter.
    p_value = special.gammaincc(degrees_of_freedom / 2.0, chi2 / 2.0)
    p_value = np.maximum(p_value, np.nextafter(np.float64(0), np.float64(1)))

    chrom_selected = chromosome[variant_index]
    position_selected = position[variant_index]
    effect_selected = np.char.upper(effect[variant_index])
    other_selected = np.char.upper(other[variant_index])
    marker_selected = marker_ids[variant_index]
    maf_selected = maf[variant_index]
    chrom_numeric = pd.to_numeric(pd.Series(chrom_selected), errors="coerce").to_numpy()
    acgt = np.array(["A", "C", "G", "T"])
    keep = np.isfinite(p_value) & (p_value > 0.0) & (p_value <= gwas_p)
    keep &= np.isfinite(maf_selected) & (maf_selected >= maf_min)
    keep &= np.isfinite(chrom_numeric) & (chrom_numeric >= 1) & (chrom_numeric <= 22)
    keep &= np.isin(effect_selected, acgt) & np.isin(other_selected, acgt)
    keep &= np.char.startswith(marker_selected, "rs")
    if exclude_mhc:
        keep &= ~(
            (chrom_numeric == 6)
            & (position_selected >= MHC_START)
            & (position_selected <= MHC_END)
        )

    selected = np.flatnonzero(keep)
    a1 = effect_selected[selected]
    a2 = other_selected[selected]
    # NumPy has no minimum/maximum loop for Unicode arrays.  The comparison
    # itself is vectorized, so select the alphabetically ordered pair directly.
    a1_first = a1 <= a2
    allele_lo = np.where(a1_first, a1, a2)
    allele_hi = np.where(a1_first, a2, a1)
    chromosome_text = chrom_numeric[selected].astype(np.int16).astype(str)
    table = pd.DataFrame(
        {
            "CHR": chromosome_text,
            "SNP": marker_selected[selected],
            "POS": position_selected[selected],
            "A1": a1,
            "A2": a2,
            "N": np.int32(manifest["n_samples"]),
            "AF1": maf_selected[selected].astype(np.float32),
            "P": p_value[selected],
            "uniqID": np.char.add(
                np.char.add(
                    np.char.add(np.char.add(chromosome_text, ":"), position_selected[selected].astype(str)),
                    np.char.add(":", allele_lo),
                ),
                np.char.add(":", allele_hi),
            ),
        }
    )
    before_dedup = len(table)
    table = table.sort_values("P").drop_duplicates("uniqID")
    table = table.sort_values(
        ["CHR", "POS"],
        key=lambda column: column.astype(int) if column.name == "CHR" else column,
    ).reset_index(drop=True)
    report = {
        "jagwas_df": degrees_of_freedom,
        "jagwas_rows": int(len(chi2)),
        "rows_before_dedup": int(before_dedup),
        "harmonized_rows": int(len(table)),
        "gwas_p": float(gwas_p),
        "maf_min": float(maf_min),
        "exclude_mhc": bool(exclude_mhc),
    }
    return table, report


def write_harmonized_npz(table: pd.DataFrame, path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: dict[str, np.ndarray] = {}
    string_columns = []
    for name in table.columns:
        values = table[name].to_numpy()
        if values.dtype.kind in {"O", "U"}:
            values = values.astype("S")
            string_columns.append(name)
        columns[name] = values
    np.savez(
        path,
        _format=np.array(HARMONIZED_FORMAT),
        _gwasP=np.array(report["gwas_p"], dtype=np.float64),
        _n_input=np.array(report["jagwas_rows"], dtype=np.int64),
        _str_cols=np.asarray(string_columns, dtype="S"),
        **columns,
    )


def run_clumping(
    harmonized_path: Path,
    output_dir: Path,
    *,
    clumping_dir: Path,
    ld_dir: Path,
    lead_p: float,
    gwas_p: float,
) -> None:
    sys.path.insert(0, str(clumping_dir))
    import fuma_clump

    output_dir.mkdir(parents=True, exist_ok=True)
    fuma_clump.run(
        str(harmonized_path),
        str(output_dir),
        str(ld_dir),
        lead_p=lead_p,
        gwas_p=gwas_p,
        no_mhc=False,  # MHC filtering was already applied during harmonization.
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genotype", type=Path, required=True)
    parser.add_argument(
        "--genotype-format",
        choices=("auto", "npy", "plink", "bgen", "pgen", "zstd"),
        default="auto",
    )
    parser.add_argument("--sample-file", type=Path)
    parser.add_argument(
        "--sample-ids", type=Path,
        help="requested sample IDs in phenotype row order; used to subset BGEN",
    )
    parser.add_argument(
        "--bgen-decode-backend", choices=("auto", "cpu", "gpu"), default="auto",
    )
    parser.add_argument("--phenotype", type=Path, required=True)
    parser.add_argument("--covariates", type=Path)
    parser.add_argument("--maf", type=Path, help="aligned .npy, or .npz containing 'maf'")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--compute-dtype", choices=("float32", "float64"), default="float32")
    parser.add_argument("--chunk-size", type=int)
    parser.add_argument("--reader-workers", type=int, default=24)
    parser.add_argument("--prefetch-chunks", type=int, default=4)
    parser.add_argument("--variant-range", type=parse_variant_range)
    parser.add_argument("--gwas-p", type=float, default=0.05)
    parser.add_argument("--lead-p", type=float, default=5e-8)
    parser.add_argument("--maf-min", type=float, default=0.01)
    parser.add_argument("--include-mhc", action="store_true")
    parser.add_argument("--reuse-scan", action="store_true")
    parser.add_argument(
        "--clumping-dir", type=Path,
        default=Path("/data484_4/zxie3/local_clumping"),
    )
    parser.add_argument(
        "--ld-dir", type=Path,
        default=Path("/data/zxie3/ld_tables"),
    )
    args = parser.parse_args()
    if not (0.0 < args.gwas_p <= 1.0 and 0.0 < args.lead_p <= args.gwas_p):
        parser.error("require 0 < lead-p <= gwas-p <= 1")
    if not (0.0 <= args.maf_min <= 0.5):
        parser.error("maf-min must be in [0, 0.5]")
    default_maf = Path(f"{args.genotype}.maf.npy")
    maf_path = args.maf or default_maf
    for path in (
        args.genotype, args.phenotype, maf_path,
        args.clumping_dir / "fuma_clump.py", args.ld_dir,
        *(() if args.covariates is None else (args.covariates,)),
        *(() if args.sample_file is None else (args.sample_file,)),
        *(() if args.sample_ids is None else (args.sample_ids,)),
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    started = time.perf_counter()
    scan_dir = args.output_dir / "torchgwas"
    sumstats_dir = scan_dir / "sumstats"
    if args.reuse_scan:
        if not (sumstats_dir / "manifest.json").is_file():
            raise FileNotFoundError(sumstats_dir / "manifest.json")
    else:
        run_linear_gwas(
            genotype=args.genotype,
            genotype_format=args.genotype_format,
            phenotype=args.phenotype,
            covariates=args.covariates,
            sample_file=args.sample_file,
            sample_ids=args.sample_ids,
            bgen_decode_backend=args.bgen_decode_backend,
            device=args.device,
            compute_dtype=args.compute_dtype,
            chunk_size=args.chunk_size,
            reader_workers=args.reader_workers,
            prefetch_chunks=args.prefetch_chunks,
            variant_range=args.variant_range,
            reduce="jagwas",
            sumstats_fields="t",
            output_dir=scan_dir,
        )
    scan_done = time.perf_counter()

    table, report = harmonized_table(
        sumstats_dir,
        maf_path,
        gwas_p=args.gwas_p,
        maf_min=args.maf_min,
        exclude_mhc=not args.include_mhc,
        variant_range=args.variant_range,
    )
    harmonized_path = args.output_dir / "jagwas_clump_input.npz"
    write_harmonized_npz(table, harmonized_path, report)
    harmonized_done = time.perf_counter()
    run_clumping(
        harmonized_path,
        args.output_dir / "loci",
        clumping_dir=args.clumping_dir,
        ld_dir=args.ld_dir,
        lead_p=args.lead_p,
        gwas_p=args.gwas_p,
    )
    clumping_done = time.perf_counter()
    report.update(
        {
            "genotype": str(args.genotype),
            "genotype_format": args.genotype_format,
            "phenotype": str(args.phenotype),
            "covariates": None if args.covariates is None else str(args.covariates),
            "sample_file": None if args.sample_file is None else str(args.sample_file),
            "sample_ids": None if args.sample_ids is None else str(args.sample_ids),
            "maf": str(maf_path),
            "variant_range": args.variant_range,
            "torchgwas_seconds": scan_done - started,
            "harmonization_seconds": harmonized_done - scan_done,
            "clumping_seconds": clumping_done - harmonized_done,
            "total_seconds": clumping_done - started,
        }
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pipeline.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
