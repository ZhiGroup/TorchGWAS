#!/usr/bin/env python3
"""Run current TorchGWAS JAGWAS reduction and FUMA-style locus clumping.

The scan writes one chi-square statistic per variant with r degrees of freedom,
r the number of traits the JAGWAS rank cutoff keeps (K for a well-conditioned
panel; see torchgwas.jagwas_projection).
This wrapper converts that indexed output to the harmonized table expected by
the lab's validated local-clumping implementation, then runs locus clumping.

--phenotype-group NAME=PATH (repeatable, instead of --phenotype) runs one
joint test per group from a single genotype pass over the groups' concatenated
phenotypes; each group's loci and pipeline.json go to OUTPUT_DIR/NAME.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import json
import multiprocessing as mp
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy import special

from torchgwas.api import run_linear_gwas
from torchgwas.io import load_genotype
from torchgwas.preprocess import residualize_and_standardize
from torchgwas.sumstats_indexed import open_indexed_sumstats

from jagwas_streaming import (
    GroupedClumpers,
    PipelinedClumper,
    build_clumping_cache,
    validate_clumping_cache,
)


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
    excluded: int = 0,
) -> tuple[pd.DataFrame, dict]:
    manifest, variant_index, chi2 = collect_jagwas(sumstats_dir)
    if chi2.ndim != 1:
        raise ValueError("this scan has JAGWAS groups; harmonize each group's column")
    variants = variant_table(sumstats_dir, maf_path, variant_range)
    return harmonize(
        variants, variant_index, chi2,
        degrees_of_freedom=manifest["df"], n_samples=manifest["n_samples"] - excluded,
        gwas_p=gwas_p, maf_min=maf_min, exclude_mhc=exclude_mhc,
    )


def variant_table(
    sumstats_dir: Path, maf_path: Path, variant_range: tuple[int, int] | None,
) -> tuple[np.ndarray, ...]:
    """(marker IDs, chromosome, position, effect, other, MAF), aligned to the scan."""
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
    return marker_ids, chromosome, position, effect, other, maf


def harmonize(
    variants: tuple[np.ndarray, ...],
    variant_index: np.ndarray,
    chi2: np.ndarray,
    *,
    degrees_of_freedom: int,
    n_samples: int,
    gwas_p: float,
    maf_min: float,
    exclude_mhc: bool,
) -> tuple[pd.DataFrame, dict]:
    marker_ids, chromosome, position, effect, other, maf = variants
    n_variants = len(marker_ids)
    degrees_of_freedom = int(degrees_of_freedom)
    if degrees_of_freedom < 1:
        raise ValueError("JAGWAS degrees of freedom must be positive")
    if len(variant_index) and (variant_index.min() < 0 or variant_index.max() >= n_variants):
        raise ValueError("indexed JAGWAS output refers outside the marker table")

    # P(chi2_K >= T) = Q(K/2, T/2).  gammaincc is the regularized upper
    # incomplete gamma.  Clamp only exact floating-point underflow so the
    # strongest valid statistics are retained by the clumper's P > 0 filter.
    p_value = special.gammaincc(degrees_of_freedom / 2.0, chi2 / 2.0)
    p_value = np.maximum(p_value, np.nextafter(np.float64(0), np.float64(1)))
    # The P filter first (the filters are a conjunction, and order is kept),
    # so the string work below touches ~gwas_p of the rows rather than all.
    candidate = np.isfinite(p_value) & (p_value > 0.0) & (p_value <= gwas_p)
    variant_index, p_value = variant_index[candidate], p_value[candidate]

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
            "N": np.int32(n_samples),
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


def ensure_clumping_cache(
    args, maf_path: Path, cache_dir: Path, n_samples: int,
) -> tuple[dict, float]:
    """Validate an existing cache or build it once from the genotype metadata."""
    started = time.perf_counter()
    if (cache_dir / "manifest.json").is_file():
        # Opening the genotype merely to learn M would defeat the point of a
        # reusable cache. Its manifest records the aligned full variant count;
        # run_linear_gwas will independently reject an incompatible input.
        manifest = json.loads((cache_dir / "manifest.json").read_text())
        manifest = validate_clumping_cache(
            cache_dir,
            n_variants=int(manifest["n_variants"]),
            maf_min=args.maf_min,
            exclude_mhc=not args.include_mhc,
        )
        return manifest, time.perf_counter() - started

    requested_sample_ids = (
        None
        if args.sample_ids is None
        else np.load(args.sample_ids, allow_pickle=False)
    )
    genotype, _sample_ids, marker_ids, _metadata = load_genotype(
        args.genotype,
        genotype_format=args.genotype_format,
        sample_file=args.sample_file,
        selected_sample_ids=requested_sample_ids,
        genotype_cache_dir=args.genotype_cache_dir,
        reader_workers=args.reader_workers,
        zstd_read_workers=args.zstd_read_workers,
        prefetch_chunks=args.prefetch_chunks,
        bgen_decode_backend=args.bgen_decode_backend,
    )
    try:
        if int(genotype.shape[0]) != n_samples:
            raise ValueError(
                f"genotype has {genotype.shape[0]} selected samples but phenotype "
                f"has {n_samples} rows"
            )
        variant_metadata = getattr(genotype, "variant_metadata", None)
        if variant_metadata is None or marker_ids is None:
            raise ValueError("genotype input does not expose variant metadata")
        manifest = build_clumping_cache(
            cache_dir,
            marker_ids=marker_ids,
            variant_metadata=variant_metadata,
            maf=load_maf(maf_path, None),
            maf_min=args.maf_min,
            exclude_mhc=not args.include_mhc,
            mhc_start=MHC_START,
            mhc_end=MHC_END,
        )
    finally:
        close = getattr(genotype, "close", None)
        if close is not None:
            close()
    return manifest, time.perf_counter() - started


def common_report(args, maf_path: Path, phenotype: Path | None = None, excluded: int | None = None) -> dict:
    phenotype = args.phenotype if phenotype is None else phenotype
    return {
        "genotype": str(args.genotype),
        "genotype_format": args.genotype_format,
        "phenotype": None if phenotype is None else str(phenotype),
        "covariates": None if args.covariates is None else str(args.covariates),
        "sample_file": None if args.sample_file is None else str(args.sample_file),
        "sample_ids": None if args.sample_ids is None else str(args.sample_ids),
        "maf": str(maf_path),
        "variant_range": args.variant_range,
        "jagwas_rcond": args.jagwas_rcond,
        "jagwas_min_residual": args.jagwas_min_residual,
        "phenotype_outlier_sd": args.phenotype_outlier_sd,
        "excluded_samples": excluded,
        "phenotype_screen_seconds": getattr(args, "screen_seconds", None),
    }


def scanned_span(variant_range, full_variant_count: int) -> tuple[int, int]:
    """(offset, count) of the scanned variants within the clumping cache's axis."""
    if variant_range is None:
        return 0, full_variant_count
    start, end = variant_range
    if end > full_variant_count:
        raise ValueError("variant range exceeds the clumping cache")
    return start, end - start


def write_report(directory: Path, report: dict) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "pipeline.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )


GROUP_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
RESERVED_GROUP_NAMES = {"torchgwas", "clump-cache", "loci"}


def parse_phenotype_group(value: str) -> tuple[str, Path]:
    name, separator, path = value.partition("=")
    if not separator or not path or not GROUP_NAME.fullmatch(name) or name in RESERVED_GROUP_NAMES:
        raise argparse.ArgumentTypeError(
            "phenotype group must be NAME=PATH, NAME of letters, digits, '.', '_' "
            f"or '-' and not one of {sorted(RESERVED_GROUP_NAMES)}"
        )
    return name, Path(path)


def parse_group_rcond(value: str) -> tuple[str, float]:
    name, separator, text = value.partition("=")
    try:
        number = float(text)
    except ValueError:
        number = None
    if not separator or number is None or not 0.0 < number < 1.0:
        raise argparse.ArgumentTypeError("a group cutoff must be NAME=VALUE with 0 < VALUE < 1")
    return name, number


def load_panels(paths) -> dict:
    """Each phenotype file read once, several at a time: they sit on a network drive.

    On the shared server one 18 MB panel took 6 s to read while a scan was
    loading the file server, and each was read twice (screen, then scan).
    """
    unique = list(dict.fromkeys(paths))
    with ThreadPoolExecutor(max_workers=min(8, len(unique))) as pool:
        return dict(zip(unique, pool.map(lambda path: np.load(path, allow_pickle=False), unique)))


def group_panel(groups, panels=None) -> tuple[np.ndarray, list[str], list[tuple[str, list[int]]]]:
    """The groups' phenotypes side by side, their trait names, and each group's columns.

    A trait is named from a NAME.traits.txt sidecar (one name per line) when
    it lists every column, else trait_<i>, prefixed with its group.
    """
    arrays, names, columns = [], [], []
    for name, path in groups:
        values = np.load(path, allow_pickle=False) if panels is None else panels[path]
        if values.ndim != 2:
            raise ValueError(f"phenotype group {name} is not a two-dimensional array")
        if arrays and values.shape[0] != arrays[0].shape[0]:
            raise ValueError(
                f"phenotype group {name} has {values.shape[0]} rows but the first "
                f"group has {arrays[0].shape[0]}"
            )
        sidecar = path.with_suffix(".traits.txt")
        labels = (
            [line.strip() for line in sidecar.read_text().splitlines() if line.strip()]
            if sidecar.is_file()
            else []
        )
        if len(labels) != values.shape[1]:
            labels = [f"trait_{index}" for index in range(values.shape[1])]
        offset = sum(array.shape[1] for array in arrays)
        columns.append((name, list(range(offset, offset + values.shape[1]))))
        names.extend(f"{name}:{label}" for label in labels)
        arrays.append(values)
    return np.concatenate(arrays, axis=1), names, columns


DEFAULT_OUTLIER_SD = 5.0
DEFAULT_RCOND = 1e-3


def phenotype_outliers(values, covariates, threshold: float) -> np.ndarray:
    """Rows of one panel beyond `threshold` SD in any covariate-residualised, standardised trait.

    The pipeline sets those rows of that panel alone to missing, and TorchGWAS's
    phenotype-missingness path mean-imputes them: the genotype, and every other
    group, keep the sample. A group's whole row goes, not the one value. On
    near-collinear panels, masking or clipping just the extreme values broke
    the traits' linear relations for those samples and made the low-variance
    directions heavier-tailed (median kurtosis 105 -> 726 at 1e-4..1e-3 of the
    largest eigenvalue on fourier_PE_L4_xyz); dropping the row made them
    Gaussian (-> 9), and with 1e-3..1e-2 from 46 to 0.6.
    """
    standardized, _ = residualize_and_standardize(np.asarray(values, dtype=np.float64), covariates)
    return (np.abs(standardized) > threshold).any(axis=1)


def resolve_defaults(args) -> None:
    """Fill the pipeline's standard QC: 5-SD phenotype rows dropped, eigen truncation at 1e-3.

    0 switches either off. --reuse-scan scans nothing, so neither applies to it.
    --jagwas-min-residual drops traits instead of eigen-directions.

    Why eigen truncation rather than dropping traits: on the collinear imaging
    panels, a one-ulp FP32 perturbation of the phenotype (what another device
    or kernel rounds differently) changed which traits greedy pivoting kept,
    df by 1-2, and T by up to 7% (rounding cutoff) or 12% (VIF <= 100): those
    panels hold near-exact ties. The eigen-directions above 1e-3 of the
    largest kept the same df, and T within 1e-7.
    """
    if args.phenotype_outlier_sd is None:
        args.phenotype_outlier_sd = None if args.reuse_scan else DEFAULT_OUTLIER_SD
    elif args.phenotype_outlier_sd == 0:
        args.phenotype_outlier_sd = None
    elif args.phenotype_outlier_sd < 0:
        raise ValueError("phenotype-outlier-sd must be positive, or 0 to keep every row")
    elif args.reuse_scan:
        raise ValueError("phenotype-outlier-sd changes the scan; it cannot apply to --reuse-scan")
    if args.jagwas_min_residual == 0:
        args.jagwas_min_residual = None
    elif args.jagwas_min_residual is not None and not 0.0 < args.jagwas_min_residual < 1.0:
        raise ValueError("jagwas-min-residual must be in (0, 1)")
    if args.jagwas_rcond is None:
        args.jagwas_rcond = (None if args.reuse_scan or args.jagwas_min_residual is not None
                             else DEFAULT_RCOND)
    elif args.jagwas_rcond == 0:
        args.jagwas_rcond = 0.0  # TorchGWAS's rounding cutoff over traits
    elif not 0.0 < args.jagwas_rcond < 1.0:
        raise ValueError("jagwas-rcond must be in (0, 1), or 0 for the rounding cutoff alone")
    elif args.jagwas_min_residual is not None:
        raise ValueError("pass --jagwas-rcond or --jagwas-min-residual, not both")
    if args.jagwas_rcond == 0 and args.jagwas_min_residual is not None:
        args.jagwas_rcond = None


def write_excluded(ids, directory: Path, flagged: np.ndarray) -> None:
    """The IDs (row numbers without --sample-ids) whose phenotype row was dropped."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "excluded_samples.txt").write_text("".join(f"{sample}\n" for sample in ids[flagged]))


def excluded_count(args, name) -> int | None:
    flagged = args.outlier_rows.get(name)
    return None if flagged is None else int(flagged.sum())


def recorded_excluded(args, name, scan_dir: Path) -> int:
    """This run's dropped rows, or with --reuse-scan those the scanning run recorded beside its scan."""
    count = excluded_count(args, name)
    if count is not None:
        return count
    record = scan_dir.resolve().parent / ("" if name is None else name) / "excluded_samples.txt"
    return len(record.read_text().splitlines()) if record.is_file() else 0


def masked_panel(args, groups):
    """group_panel with each group's own outlier rows set to missing in its columns only."""
    panel, trait_names, columns = group_panel(groups, args.panels)
    for name, group_columns in columns:
        flagged = args.outlier_rows.get(name)
        if flagged is not None and flagged.any():
            if not np.issubdtype(panel.dtype, np.floating):
                panel = panel.astype(np.float64)
            panel[np.ix_(flagged, group_columns)] = np.nan
    return panel, trait_names, columns


def single_phenotype(args):
    """--phenotype as a path, or with its outlier rows set to missing."""
    flagged = args.outlier_rows.get(None)
    if flagged is None or not flagged.any():
        return args.phenotype
    values = np.array(args.panels[args.phenotype], dtype=np.float64)
    values[flagged] = np.nan
    return values


def group_cutoffs(args, columns):
    """jagwas_groups entries, with a group's --group-rcond or --group-min-residual where given."""
    entries = []
    for name, group_columns in columns:
        cutoff = {}
        if name in args.group_rcond:
            cutoff["rcond"] = args.group_rcond[name]
        if name in args.group_min_residual:
            cutoff["min_residual"] = args.group_min_residual[name]
        entries.append((name, group_columns, cutoff or None))
    return entries


def scan_arguments(args) -> dict:
    """run_linear_gwas arguments shared by every mode."""
    return dict(
        genotype=args.genotype,
        genotype_format=args.genotype_format,
        covariates=args.covariates,
        sample_file=args.sample_file,
        sample_ids=args.sample_ids,
        jagwas_rcond=args.jagwas_rcond,
        jagwas_min_residual=args.jagwas_min_residual,
        bgen_decode_backend=args.bgen_decode_backend,
        genotype_cache_dir=args.genotype_cache_dir,
        device=args.device,
        compute_dtype=args.compute_dtype,
        chunk_size=args.chunk_size,
        reader_workers=args.reader_workers,
        zstd_read_workers=args.zstd_read_workers,
        prefetch_chunks=args.prefetch_chunks,
        variant_range=args.variant_range,
        reduce="jagwas",
    )


def run_groups_overlapped(args, groups, maf_path: Path, started: float) -> dict:
    """One genotype pass for every group, with one overlapped clumping worker per group."""
    all_samples = int(args.panels[groups[0][1]].shape[0])
    cache_dir = args.clump_cache or (args.output_dir / "clump-cache")
    cache_manifest, cache_seconds = ensure_clumping_cache(args, maf_path, cache_dir, all_samples)
    variant_offset, scan_variant_count = scanned_span(
        args.variant_range, int(cache_manifest["n_variants"])
    )
    # The workers fork here, before this process holds the panel or CUDA.
    workers = []
    try:
        for name, _path in groups:
            workers.append(
                PipelinedClumper(
                    cache_dir=cache_dir,
                    output_dir=args.output_dir / name / "loci",
                    clumping_dir=args.clumping_dir,
                    ld_dir=args.ld_dir,
                    n_samples=all_samples - (excluded_count(args, name) or 0),
                    n_variants=scan_variant_count,
                    gwas_p=args.gwas_p,
                    lead_p=args.lead_p,
                    variant_offset=variant_offset,
                )
            )
    except BaseException:
        for worker in workers:
            worker.abort()
        raise
    clumpers = GroupedClumpers(workers)
    scan_started = time.perf_counter()
    try:
        panel, trait_names, columns = masked_panel(args, groups)
        gwas_result = run_linear_gwas(
            **scan_arguments(args),
            phenotype=panel,
            trait_columns=trait_names,
            jagwas_groups=group_cutoffs(args, columns),
            sumstats_format="none",
            result_chunk_callback=clumpers.consume,
            jagwas_rank_callback=clumpers.set_rank,
            output_dir=args.output_dir / "torchgwas",
        )
        del panel
        observed_variant_count = int(gwas_result.run_metadata["genotype_shape"][1])
        if observed_variant_count != scan_variant_count:
            raise ValueError(
                f"clumping cache expects {scan_variant_count} scanned "
                f"variants but TorchGWAS reported {observed_variant_count}"
            )
        scan_done = time.perf_counter()
        collected = clumpers.finish()
    except BaseException:
        clumpers.abort()
        raise
    shared = {
        "pipeline_mode": "overlapped_groups",
        "groups_in_scan": len(groups),
        "clump_cache": str(cache_dir),
        "clump_cache_seconds": cache_seconds,
        "clump_cache_eligible_variants": int(cache_manifest["eligible_variants"]),
        "scanned_variants": scan_variant_count,
        "torchgwas_seconds": scan_done - scan_started,
    }
    summaries = []
    for (name, path), (_name, group_columns), worker, (worker_report, done) in zip(
        groups, columns, workers, collected
    ):
        report = {
            **common_report(args, maf_path, path, excluded_count(args, name)),
            **worker_report,
            **shared,
            "group": name,
            "traits": len(group_columns),
            "jagwas_df": worker.degrees_of_freedom,
            "jagwas_rank": worker.rank_report,
            "jagwas_rows": int(worker_report["finite_jagwas_rows"]),
            "post_scan_wait_seconds": done - scan_done,
            "total_seconds_excluding_cache": done - scan_started,
            "total_seconds": done - started,
        }
        write_report(args.output_dir / name, report)
        summaries.append({key: report[key] for key in (
            "group", "traits", "jagwas_df", "excluded_samples", "harmonized_rows", "ok")})
    finished = max(done for _report, done in collected)
    return {
        **common_report(args, maf_path),
        **shared,
        "phenotype_groups": {name: str(path) for name, path in groups},
        "groups": summaries,
        "post_scan_wait_seconds": finished - scan_done,
        "total_seconds_excluding_cache": finished - scan_started,
        "total_seconds": finished - started,
    }


def run_groups_sequential(args, groups, maf_path: Path, started: float) -> dict:
    """One indexed grouped scan (or --reuse-scan), then each group's harmonization and clumping.

    The variant table is loaded once. Clumping runs in spawned processes,
    overlapping this process's harmonization of the next group.
    """
    sumstats_dir = args.output_dir / "torchgwas" / "sumstats"
    if args.reuse_scan:
        if not (sumstats_dir / "manifest.json").is_file():
            raise FileNotFoundError(sumstats_dir / "manifest.json")
    else:
        panel, trait_names, columns = masked_panel(args, groups)
        run_linear_gwas(
            **scan_arguments(args),
            phenotype=panel,
            trait_columns=trait_names,
            jagwas_groups=group_cutoffs(args, columns),
            sumstats_fields="t",
            output_dir=args.output_dir / "torchgwas",
        )
        del panel
    scan_done = time.perf_counter()
    manifest, variant_index, chi2 = collect_jagwas(sumstats_dir)
    scanned = manifest.get("groups")
    if scanned is None:
        raise ValueError(f"{sumstats_dir} is not a grouped JAGWAS scan")
    absent = [name for name, _path in groups if name not in scanned]
    if absent:
        raise ValueError(f"{sumstats_dir} has no JAGWAS group(s) {absent}")
    chi2 = chi2.reshape(len(chi2), len(scanned))
    variants = variant_table(sumstats_dir, maf_path, args.variant_range)
    workers = args.clump_workers or min(8, len(groups))
    pending = []
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
        for name, path in groups:
            column = scanned.index(name)
            group_started = time.perf_counter()
            finite = np.isfinite(chi2[:, column])
            excluded = recorded_excluded(args, name, sumstats_dir.parent)
            table, report = harmonize(
                variants, variant_index[finite], chi2[finite, column],
                degrees_of_freedom=manifest["df"][column],
                n_samples=manifest["n_samples"] - excluded,
                gwas_p=args.gwas_p, maf_min=args.maf_min, exclude_mhc=not args.include_mhc,
            )
            harmonized_path = args.output_dir / name / "jagwas_clump_input.npz"
            write_harmonized_npz(table, harmonized_path, report)
            report.update(
                {
                    **common_report(args, maf_path, path, excluded),
                    "group": name,
                    "jagwas_rank": manifest["jagwas_rank"][column],
                    "harmonization_seconds": time.perf_counter() - group_started,
                }
            )
            future = pool.submit(
                run_clumping,
                harmonized_path,
                args.output_dir / name / "loci",
                clumping_dir=args.clumping_dir,
                ld_dir=args.ld_dir,
                lead_p=args.lead_p,
                gwas_p=args.gwas_p,
            )
            pending.append((name, future, report))
        summaries = []
        for name, future, report in pending:
            future.result()
            done = time.perf_counter()
            report.update(
                {
                    "ok": True,
                    "pipeline_mode": "sequential_groups",
                    "groups_in_scan": len(scanned),
                    "torchgwas_seconds": scan_done - started,
                    "total_seconds": done - started,
                }
            )
            write_report(args.output_dir / name, report)
            summaries.append({key: report[key] for key in (
                "group", "jagwas_df", "excluded_samples", "harmonized_rows", "ok")})
    finished = time.perf_counter()
    return {
        **common_report(args, maf_path),
        "pipeline_mode": "sequential_groups",
        "phenotype_groups": {name: str(path) for name, path in groups},
        "groups": summaries,
        "torchgwas_seconds": scan_done - started,
        "total_seconds": finished - started,
    }


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
    parser.add_argument("--phenotype", type=Path)
    parser.add_argument(
        "--phenotype-group",
        type=parse_phenotype_group,
        action="append",
        metavar="NAME=PATH",
        help=(
            "repeatable, instead of --phenotype: one joint test per group from a "
            "single genotype pass; each group's loci and pipeline.json go to "
            "OUTPUT_DIR/NAME"
        ),
    )
    parser.add_argument(
        "--clump-workers",
        type=int,
        help="concurrent group clumping processes with --sequential/--reuse-scan (default min(8, groups))",
    )
    parser.add_argument(
        "--jagwas-rcond",
        type=float,
        help=(
            f"eigen truncation: keep the trait correlation's eigen-directions above "
            f"RCOND x the largest eigenvalue (default {DEFAULT_RCOND:g} for a new scan; "
            f"0 keeps only the rounding cutoff over traits)"
        ),
    )
    parser.add_argument(
        "--group-rcond",
        type=parse_group_rcond,
        action="append",
        default=[],
        metavar="NAME=RCOND",
        help="repeatable: one group's eigen truncation, overriding the global cutoff",
    )
    parser.add_argument(
        "--jagwas-min-residual",
        type=float,
        help=(
            "drop traits instead of eigen-directions: while less than this fraction "
            "of a trait's variance is its own given the traits kept before it, i.e. "
            "VIF above 1/VALUE (replaces the default eigen truncation)"
        ),
    )
    parser.add_argument(
        "--group-min-residual",
        type=parse_group_rcond,
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="repeatable: one group's trait-dropping threshold, overriding the global cutoff",
    )
    parser.add_argument(
        "--phenotype-outlier-sd",
        type=float,
        help=(
            f"set a panel's rows to missing where any of its covariate-residualised traits "
            f"is beyond this many SD; per group, genotypes untouched (default "
            f"{DEFAULT_OUTLIER_SD:g} for a new scan; 0 keeps every row)"
        ),
    )
    parser.add_argument("--covariates", type=Path)
    parser.add_argument("--maf", type=Path, help="aligned .npy, or .npz containing 'maf'")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--compute-dtype", choices=("float32", "float64"), default="float32")
    parser.add_argument("--chunk-size", type=int)
    parser.add_argument("--reader-workers", type=int, default=24)
    parser.add_argument("--zstd-read-workers", type=int)
    parser.add_argument("--prefetch-chunks", type=int, default=4)
    parser.add_argument("--genotype-cache-dir", type=Path)
    parser.add_argument("--variant-range", type=parse_variant_range)
    parser.add_argument("--gwas-p", type=float, default=0.05)
    parser.add_argument("--lead-p", type=float, default=5e-8)
    parser.add_argument("--maf-min", type=float, default=0.01)
    parser.add_argument("--include-mhc", action="store_true")
    parser.add_argument("--reuse-scan", action="store_true")
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="write indexed JAGWAS output, then harmonize and clump sequentially",
    )
    parser.add_argument(
        "--clump-cache",
        type=Path,
        help="reusable compact variant cache; defaults inside output-dir",
    )
    parser.add_argument(
        "--clumping-dir", type=Path,
        default=Path("/path/to/local_clumping"),
    )
    parser.add_argument(
        "--ld-dir", type=Path,
        default=Path("/path/to/ld_tables"),
    )
    args = parser.parse_args()
    if not (0.0 < args.gwas_p <= 1.0 and 0.0 < args.lead_p <= args.gwas_p):
        parser.error("require 0 < lead-p <= gwas-p <= 1")
    if not (0.0 <= args.maf_min <= 0.5):
        parser.error("maf-min must be in [0, 0.5]")
    if (args.phenotype is None) == (args.phenotype_group is None):
        parser.error("pass either --phenotype or --phenotype-group")
    if args.phenotype_group is not None:
        group_names = [name for name, _path in args.phenotype_group]
        if len(set(group_names)) != len(group_names):
            parser.error("phenotype group names must be unique")
    if args.clump_workers is not None and args.clump_workers < 1:
        parser.error("clump-workers must be positive")
    try:
        resolve_defaults(args)
    except ValueError as error:
        parser.error(str(error))
    groups = {name for name, _path in args.phenotype_group or ()}
    for option in ("group_rcond", "group_min_residual"):
        entries = dict(getattr(args, option))
        if len(entries) != len(getattr(args, option)):
            parser.error(f"each group takes one --{option.replace('_', '-')}")
        if set(entries) - groups:
            parser.error(f"--{option.replace('_', '-')} names no --phenotype-group: {sorted(set(entries) - groups)}")
        setattr(args, option, entries)
    if set(args.group_rcond) & set(args.group_min_residual):
        parser.error("a group takes --group-rcond or --group-min-residual, not both")
    default_maf = Path(f"{args.genotype}.maf.npy")
    maf_path = args.maf or default_maf
    if args.genotype_format in {"auto", "zstd"} and not args.genotype.exists():
        zstd_paths = tuple(
            Path(f"{args.genotype}{suffix}")
            for suffix in (
                ".zst", ".idx.npz", ".samples.tsv", ".variants.tsv",
                ".complete.json",
            )
        )
        missing_zstd_paths = [path for path in zstd_paths if not path.exists()]
        if missing_zstd_paths:
            raise FileNotFoundError(missing_zstd_paths[0])
    elif not args.genotype.exists():
        raise FileNotFoundError(args.genotype)

    phenotype_paths = (
        (args.phenotype,)
        if args.phenotype_group is None
        else tuple(path for _name, path in args.phenotype_group)
    )
    for path in (
        *phenotype_paths, maf_path,
        args.clumping_dir / "fuma_clump.py", args.ld_dir,
        *(() if args.covariates is None else (args.covariates,)),
        *(() if args.sample_file is None else (args.sample_file,)),
        *(() if args.sample_ids is None else (args.sample_ids,)),
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    started = time.perf_counter()
    args.panels = {} if args.reuse_scan else load_panels(phenotype_paths)
    # Each panel's own outlier rows, set to missing in that panel only.
    args.outlier_rows = {}
    if args.phenotype_outlier_sd is not None:
        screen_started = time.perf_counter()
        covariates = None if args.covariates is None else np.load(args.covariates, allow_pickle=False)
        panels = [(None, args.phenotype)] if args.phenotype_group is None else list(args.phenotype_group)
        with ThreadPoolExecutor(max_workers=min(8, len(args.panels))) as pool:
            flags_by_path = dict(zip(args.panels, pool.map(
                lambda values: phenotype_outliers(values, covariates, args.phenotype_outlier_sd),
                args.panels.values())))
        ids = (np.arange(len(next(iter(flags_by_path.values())))) if args.sample_ids is None
               else np.load(args.sample_ids, allow_pickle=False))
        for name, path in panels:
            flagged = args.outlier_rows[name] = flags_by_path[path]
            write_excluded(ids, args.output_dir if name is None else args.output_dir / name, flagged)
            print(f"{name or 'phenotype'}: {int(flagged.sum())} of {len(flagged)} phenotype rows beyond "
                  f"{args.phenotype_outlier_sd:g} SD set to missing", flush=True)
        args.screen_seconds = time.perf_counter() - screen_started
    if args.phenotype_group is not None:
        run = (
            run_groups_sequential
            if args.reuse_scan or args.sequential
            else run_groups_overlapped
        )
        report = run(args, args.phenotype_group, maf_path, started)
        write_report(args.output_dir, report)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    scan_dir = args.output_dir / "torchgwas"
    sumstats_dir = scan_dir / "sumstats"
    if args.reuse_scan or args.sequential:
        if args.reuse_scan:
            if not (sumstats_dir / "manifest.json").is_file():
                raise FileNotFoundError(sumstats_dir / "manifest.json")
        else:
            run_linear_gwas(
                **scan_arguments(args),
                phenotype=single_phenotype(args),
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
            excluded=recorded_excluded(args, None, scan_dir),
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
                **common_report(args, maf_path, None, excluded_count(args, None)),
                "pipeline_mode": "sequential",
                "torchgwas_seconds": scan_done - started,
                "harmonization_seconds": harmonized_done - scan_done,
                "clumping_seconds": clumping_done - harmonized_done,
                "total_seconds": clumping_done - started,
            }
        )
    else:
        phenotype_shape = np.load(
            args.phenotype, mmap_mode="r", allow_pickle=False
        ).shape
        if len(phenotype_shape) != 2:
            raise ValueError("phenotype must be a two-dimensional NumPy array")
        all_samples = int(phenotype_shape[0])
        n_samples = all_samples - (excluded_count(args, None) or 0)
        cache_dir = args.clump_cache or (args.output_dir / "clump-cache")
        cache_manifest, cache_seconds = ensure_clumping_cache(
            args, maf_path, cache_dir, all_samples
        )
        variant_offset, scan_variant_count = scanned_span(
            args.variant_range, int(cache_manifest["n_variants"])
        )
        clumper = PipelinedClumper(
            cache_dir=cache_dir,
            output_dir=args.output_dir / "loci",
            clumping_dir=args.clumping_dir,
            ld_dir=args.ld_dir,
            n_samples=n_samples,
            n_variants=scan_variant_count,
            gwas_p=args.gwas_p,
            lead_p=args.lead_p,
            variant_offset=variant_offset,
        )
        scan_started = time.perf_counter()
        try:
            gwas_result = run_linear_gwas(
                **scan_arguments(args),
                phenotype=single_phenotype(args),
                sumstats_format="none",
                result_chunk_callback=clumper.consume,
                jagwas_rank_callback=clumper.set_rank,
                output_dir=scan_dir,
            )
            observed_variant_count = int(
                gwas_result.run_metadata["genotype_shape"][1]
            )
            if observed_variant_count != scan_variant_count:
                raise ValueError(
                    f"clumping cache expects {scan_variant_count} scanned "
                    f"variants but TorchGWAS reported {observed_variant_count}"
                )
            scan_done = time.perf_counter()
            pipeline_report = clumper.finish()
        except BaseException:
            clumper.abort()
            raise
        finished = time.perf_counter()
        report = {
            **common_report(args, maf_path, None, excluded_count(args, None)),
            **pipeline_report,
            "pipeline_mode": "overlapped",
            "clump_cache": str(cache_dir),
            "clump_cache_seconds": cache_seconds,
            "clump_cache_eligible_variants": int(
                cache_manifest["eligible_variants"]
            ),
            "jagwas_df": clumper.degrees_of_freedom,
            "jagwas_rank": clumper.rank_report,
            "jagwas_rows": int(pipeline_report["finite_jagwas_rows"]),
            "scanned_variants": scan_variant_count,
            "torchgwas_seconds": scan_done - scan_started,
            "post_scan_wait_seconds": finished - scan_done,
            "total_seconds_excluding_cache": finished - scan_started,
            "total_seconds": finished - started,
        }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pipeline.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
