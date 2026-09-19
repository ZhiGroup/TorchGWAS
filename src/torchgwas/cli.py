from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .api import run_linear_gwas
from .datasets import get_toy_dataset_paths
from .io import (
    align_table_to_samples,
    load_array,
    load_bgen_genotype,
    load_genotype,
    load_pgen_genotype,
    load_vector,
)
from .preprocess import prepare_inputs_for_prep, residualize_and_standardize
from .pgen import (
    DEFAULT_PGEN_COMPRESSION_WORKERS,
    DEFAULT_PGEN_DECODE_BATCH_SIZE,
    DEFAULT_PGEN_DECODE_WORKERS,
)
from .streaming import ChunkedGenotype
from .utils import mkdir, write_json


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="torchgwas")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prep", help="Validate and preprocess phenotype/covariates")
    prep.add_argument("--genotype", required=True)
    prep.add_argument("--genotype-format", default="auto", choices=["auto", "plink", "bgen", "pgen", "zstd"])
    prep.add_argument("--phenotype", default=None)
    prep.add_argument("--phenotype-table", default=None)
    prep.add_argument("--covariates", default=None)
    prep.add_argument("--covariates-table", default=None)
    prep.add_argument("--trait-columns", default=None)
    prep.add_argument("--covariate-columns", default=None)
    prep.add_argument("--sample-id-column", default="IID")
    prep.add_argument("--sample-ids", default=None)
    prep.add_argument("--sample-file", default=None)
    prep.add_argument("--bgen-decode-backend", choices=["auto", "cpu", "gpu"], default="auto")
    prep.add_argument("--pvar", default=None)
    prep.add_argument("--psam", default=None)
    prep.add_argument("--pgen-mode", default="auto", choices=["auto", "hardcall", "dosage"])
    prep.add_argument("--pgen-decode-workers", type=int, default=None,
                        help="PGEN decode readers; defaults to --reader-workers")
    prep.add_argument("--pgen-decode-batch-size", type=int, default=DEFAULT_PGEN_DECODE_BATCH_SIZE)
    prep.add_argument("--pgen-compression-workers", type=int, default=DEFAULT_PGEN_COMPRESSION_WORKERS)
    prep.add_argument("--genotype-cache-dir", default=None)
    prep.add_argument("--bim", default=None)
    prep.add_argument("--fam", default=None)
    prep.add_argument("--reader-workers", type=int, default=24)
    prep.add_argument("--prefetch-chunks", type=int, default=4)
    prep.add_argument("--output-dir", required=True)

    linear = subparsers.add_parser("linear", help="Run linear GWAS")
    linear.add_argument("--genotype", required=True)
    linear.add_argument(
        "--hardcall-store", default=None,
        help="read genotype bytes from a zstd hard-call store instead of the "
             ".bed. The .bim/.fam are still used, so this is a transport swap "
             "only. Measured 2.15x faster at K=1 and no faster at K=512, "
             "where the write and the GEMM set the floor.")
    linear.add_argument("--genotype-format", default="auto", choices=["auto", "plink", "bgen", "pgen", "zstd"])
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
    linear.add_argument("--bgen-decode-backend", choices=["auto", "cpu", "gpu"], default="auto")
    linear.add_argument("--pvar", default=None)
    linear.add_argument("--psam", default=None)
    linear.add_argument("--pgen-mode", default="auto", choices=["auto", "hardcall", "dosage"])
    linear.add_argument("--pgen-decode-workers", type=int, default=None,
                        help="PGEN decode readers; defaults to --reader-workers")
    linear.add_argument("--pgen-decode-batch-size", type=int, default=DEFAULT_PGEN_DECODE_BATCH_SIZE)
    linear.add_argument("--pgen-compression-workers", type=int, default=DEFAULT_PGEN_COMPRESSION_WORKERS)
    linear.add_argument("--genotype-cache-dir", default=None)
    linear.add_argument("--bim", default=None)
    linear.add_argument("--fam", default=None)
    linear.add_argument("--reader-workers", type=int, default=None)
    linear.add_argument("--prefetch-chunks", type=int, default=None)
    linear.add_argument("--compute-dtype", default="auto", choices=["auto", "float32", "float64"])
    linear.add_argument("--device", default="auto")
    linear.add_argument("--chunk-size", type=int, default=None)
    linear.add_argument("--pipeline-profile", default=None,
                        help="JSON with explicit hardware and input costs for the shared calculator")
    linear.add_argument("--p-value-threshold", type=float, default=None)
    # TWO MODES, and the CLI previously offered NEITHER of them: it exposed
    # only the internal top-k spellings, so the reduction a user actually wants
    # was unreachable from the command line while four pieces of machinery
    # were. `significant` is the default answer; `jagwas` needs the whole
    # phenotype resident and refuses when it will not fit.
    linear.add_argument(
        "--reduce",
        choices=("significant", "jagwas"),
        default=None,
        help="reduce across traits on the device instead of writing the full "
             "variant x trait table: 'significant' keeps the pairs clearing "
             "the threshold (default 5e-8/K), 'jagwas' computes the quadratic "
             "form over the whole trait set",
    )
    linear.add_argument(
        "--significance-threshold", type=float, default=None,
        help="p-value threshold for --reduce significant; default 5e-8/K, "
             "which already carries the correction for the effective number "
             "of independent common variants",
    )
    linear.add_argument(
        "--trait-block", type=int, default=None,
        help="process the traits in blocks of this many columns, so the device "
             "never holds the full (samples x traits) design. Requires "
             "--reduce. Costs one extra pass over the genotypes per block",
    )
    linear.add_argument(
        "--sumstats-format",
        default="binary",
        choices=["binary", "none"],
        help=(
            "binary writes tiled float32 beta/t_stat/-log10(P) "
            "arrays under sumstats/ (12 bytes per marker-trait cell); "
            "none runs the scan and discards results, isolating scan cost"
        ),
    )
    linear.add_argument(
        "--sumstats-fields",
        default="beta+t",
        choices=["beta+t", "t"],
        help=(
            "binary store contents: 'beta+t' (default, 12 bytes per cell) keeps "
            "effect size, statistic and -log10(P); 't' writes t and -log10(P) at "
            "8 bytes per cell but is screening-only, with no effect size or standard error, so "
            "it cannot be meta-analysed"
        ),
    )
    linear.add_argument(
        "--sumstats-block-bytes",
        type=int,
        default=None,
        help="binary write block size in bytes (default 16 MiB)",
    )
    linear.add_argument(
        "--sumstats-queue-depth",
        type=int,
        default=None,
        help="binary write blocks in flight per array (default 3)",
    )
    linear.add_argument(
        "--no-sumstats-fsync",
        dest="sumstats_fsync",
        action="store_false",
        help="skip the final fsync; reported write time then excludes durability",
    )
    linear.add_argument(
        "--no-sumstats-variant-ids",
        dest="sumstats_variant_ids",
        action="store_false",
        help="skip writing variant_ids.txt next to the binary store",
    )
    linear.add_argument("--output-dir", required=True)

    demo = subparsers.add_parser("demo", help="Run the bundled toy example")
    demo.add_argument("--output-dir", required=True)

    convert = subparsers.add_parser(
        "convert-pgen", help="Convert PLINK 2 PGEN to the native zstd scan store"
    )
    convert.add_argument("--genotype", required=True)
    convert.add_argument("--pvar", default=None)
    convert.add_argument("--psam", default=None)
    convert.add_argument("--sample-ids", default=None)
    convert.add_argument("--cache-dir", required=True)
    convert.add_argument("--pgen-mode", default="auto", choices=["auto", "hardcall", "dosage"])
    convert.add_argument("--decode-workers", type=int, default=DEFAULT_PGEN_DECODE_WORKERS)
    convert.add_argument("--decode-batch-size", type=int, default=DEFAULT_PGEN_DECODE_BATCH_SIZE)
    convert.add_argument("--compression-workers", type=int, default=DEFAULT_PGEN_COMPRESSION_WORKERS)
    convert.add_argument("--zstd-chunk-size", type=int, default=2500)
    convert.add_argument("--zstd-level", type=int, default=15)
    convert.add_argument("--reader-workers", type=int, default=4)
    convert.add_argument("--prefetch-chunks", type=int, default=4)
    convert.add_argument("--output-json", default=None)
    store = subparsers.add_parser(
        "build-hardcall-store",
        help="Compress a PLINK .bed into a zstd hard-call store")
    store.add_argument("--genotype", required=True,
                       help="a .bed, or the triplet prefix")
    store.add_argument("--out", required=True, help="store prefix")
    store.add_argument("--frame-variants", type=int, default=2048,
                       help="variants per independently decodable frame; a "
                            "chunk that straddles frames pays for the whole "
                            "covering frames, so prefer a divisor of the "
                            "scan chunk (default 2048)")
    store.add_argument("--level", type=int, default=3,
                       help="zstd level. 3 measured 6.4x; level 10 reaches "
                            "6.7x for 5.4x the compression cost (default 3)")
    store.add_argument("--workers", type=int, default=8)
    store.add_argument("--verify", type=int, default=20,
                       help="random ranges to read back and compare against "
                            "the .bed, 0 to skip (default 20)")
    convert_bgen = subparsers.add_parser("convert-bgen", help="Optionally cache BGEN as a reusable zstd store")
    convert_bgen.add_argument("--genotype", required=True)
    convert_bgen.add_argument("--sample-file", default=None)
    convert_bgen.add_argument("--cache-dir", required=True)
    convert_bgen.add_argument("--reader-workers", type=int, default=4)
    return parser


def _run_build_hardcall_store(args) -> int:
    """Compress a PLINK `.bed` into a zstd hard-call store, and verify it.

    Verification is ON by default and is not ceremony: the store substitutes
    for the genotype bytes of every scan that uses it, so a silent corruption
    would produce plausible association results for the wrong markers. Reading
    random ranges back and comparing against the `.bed` costs seconds against
    an encode measured at 28 s for a full genome.
    """
    import os
    import time

    import numpy as np

    from .bed import resolve_plink_triplet
    from .hardcall_store import HardcallStore, encode_hardcall_store

    bed_path, _bim_path, fam_path = resolve_plink_triplet(args.genotype)
    for path in (bed_path, fam_path):
        if not path.is_file():
            raise SystemExit(f"missing {path}")
    with fam_path.open("rb") as fam_handle:
        samples = sum(1 for _ in fam_handle)
    stride = (samples + 3) // 4
    variants = (bed_path.stat().st_size - 3) // stride

    handle = os.open(str(bed_path), os.O_RDONLY)
    try:
        header = os.pread(handle, 3, 0)
        if header[:2] != b"\x6c\x1b":
            raise SystemExit(f"not a PLINK 1 bed file: {bed_path}")
        if header[2] != 1:
            raise SystemExit("sample-major bed is not supported")

        def read_packed(start: int, end: int) -> np.ndarray:
            want = (end - start) * stride
            out = bytearray(want)
            view = memoryview(out)
            got = 0
            while got < want:
                amount = os.preadv(handle, [view[got:]],
                                   3 + start * stride + got)
                if amount == 0:
                    raise OSError(f"unexpected EOF in {bed_path}")
                got += amount
            return np.frombuffer(bytes(out), dtype=np.uint8).reshape(
                end - start, stride)

        print(f"{bed_path}: {variants:,} variants x {samples:,} samples")
        started = time.perf_counter()
        manifest = encode_hardcall_store(
            read_packed, variants, samples, args.out,
            frame_variants=args.frame_variants, level=args.level,
            workers=args.workers)
        elapsed = time.perf_counter() - started
        print(f"{manifest['ratio']:.3f}x  "
              f"{manifest['raw_bytes'] / 1e9:.2f} -> "
              f"{manifest['compressed_bytes'] / 1e9:.2f} GB in {elapsed:,.1f} s")

        if args.verify > 0:
            rng = np.random.default_rng(0)
            checked = 0
            with HardcallStore(args.out) as store:
                for _ in range(args.verify):
                    start = int(rng.integers(0, variants))
                    end = int(min(variants, start + rng.integers(1, 20_000)))
                    if not np.array_equal(store.read_packed(start, end),
                                          read_packed(start, end)):
                        raise SystemExit(
                            f"MISMATCH in [{start}, {end}) -- store not written")
                    checked += end - start
                for start, end in ((0, 1), (variants - 1, variants)):
                    if not np.array_equal(store.read_packed(start, end),
                                          read_packed(start, end)):
                        raise SystemExit(f"MISMATCH at edge [{start}, {end})")
            print(f"verified {checked:,} variants byte-identical, plus edges")
    finally:
        os.close(handle)
    return 0


def _run_prep(args) -> int:
    requested_sample_ids = None
    if args.sample_ids is not None:
        requested_sample_ids = (
            load_array(args.sample_ids)
            if str(args.sample_ids).lower().endswith(".npy")
            else load_vector(args.sample_ids)
        )
    genotype, sample_ids, _, genotype_meta = load_genotype(
        args.genotype,
        genotype_format=args.genotype_format,
        bim=args.bim,
        fam=args.fam,
        sample_file=args.sample_file,
        bgen_decode_backend=args.bgen_decode_backend,
        pvar=args.pvar,
        psam=args.psam,
        selected_sample_ids=requested_sample_ids,
        genotype_cache_dir=args.genotype_cache_dir,
        reader_workers=args.reader_workers,
        prefetch_chunks=args.prefetch_chunks,
        pgen_mode=args.pgen_mode,
        pgen_decode_workers=args.pgen_decode_workers,
        pgen_decode_batch_size=args.pgen_decode_batch_size,
        pgen_compression_workers=args.pgen_compression_workers,
    )
    if args.sample_ids is not None:
        if isinstance(genotype, ChunkedGenotype) and hasattr(genotype, "select_samples"):
            genotype.select_samples(requested_sample_ids)
            sample_ids = genotype.sample_ids
        elif sample_ids is None:
            sample_ids = requested_sample_ids
        elif not np.array_equal(
            np.asarray(sample_ids, dtype=str),
            np.asarray(requested_sample_ids, dtype=str),
        ):
            raise ValueError(
                f"{type(genotype).__name__} does not support sample subsetting; "
                "pre-align or convert the genotype input first"
            )
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
    genotype_array = genotype.genotype if isinstance(genotype, ChunkedGenotype) else genotype
    phenotype, covariates, qc = prepare_inputs_for_prep(
        genotype_array,
        phenotype,
        covariates,
        genotype_chunk_size=getattr(genotype_array, "preferred_chunk_size", None),
    )
    pheno_proc, q_matrix = residualize_and_standardize(phenotype, covariates)
    out = mkdir(args.output_dir)
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
        bgen_decode_backend=args.bgen_decode_backend,
        pvar=args.pvar,
        psam=args.psam,
        pgen_mode=args.pgen_mode,
        pgen_decode_workers=args.pgen_decode_workers,
        pgen_decode_batch_size=args.pgen_decode_batch_size,
        pgen_compression_workers=args.pgen_compression_workers,
        genotype_cache_dir=args.genotype_cache_dir,
        hardcall_store=args.hardcall_store,
        bim=args.bim,
        fam=args.fam,
        reader_workers=args.reader_workers,
        prefetch_chunks=args.prefetch_chunks,
        sample_id_column=args.sample_id_column,
        marker_ids=args.marker_ids,
        sample_ids=args.sample_ids,
        compute_dtype=args.compute_dtype,
        device=args.device,
        chunk_size=args.chunk_size,
        pipeline_profile=args.pipeline_profile,
        p_value_threshold=args.p_value_threshold,
        reduce=args.reduce,
        significance_threshold=args.significance_threshold,
        trait_block=args.trait_block,
        sumstats_format=args.sumstats_format,
        sumstats_block_bytes=args.sumstats_block_bytes,
        sumstats_queue_depth=args.sumstats_queue_depth,
        sumstats_fsync=args.sumstats_fsync,
        sumstats_variant_ids=args.sumstats_variant_ids,
        sumstats_fields=args.sumstats_fields,
        output_dir=args.output_dir,
    )
    return 0


def _run_convert_pgen(args) -> int:
    selected_sample_ids = None
    if args.sample_ids is not None:
        selected_sample_ids = (
            load_array(args.sample_ids)
            if str(args.sample_ids).lower().endswith(".npy")
            else load_vector(args.sample_ids)
        )
    genotype, _, _ = load_pgen_genotype(
        args.genotype,
        pvar=args.pvar,
        psam=args.psam,
        selected_sample_ids=selected_sample_ids,
        cache_dir=args.cache_dir,
        reader_workers=args.reader_workers,
        prefetch_chunks=args.prefetch_chunks,
        pgen_mode=args.pgen_mode,
        decode_workers=args.decode_workers,
        decode_batch_size=args.decode_batch_size,
        zstd_chunk_size=args.zstd_chunk_size,
        zstd_level=args.zstd_level,
        compression_workers=args.compression_workers,
    )
    manifest = json.loads(genotype.manifest_path.read_text())
    payload = {"cache_prefix": str(genotype.prefix), **manifest}
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_json is None:
        print(rendered, end="")
    else:
        Path(args.output_json).write_text(rendered)
    return 0


def _run_demo(args) -> int:
    toy = get_toy_dataset_paths()
    out = mkdir(args.output_dir)
    linear_out = out / "linear"
    linear_result = run_linear_gwas(
        genotype=np.load(toy["genotype"], allow_pickle=False),
        phenotype=np.load(toy["phenotype"], allow_pickle=False),
        covariates=np.load(toy["covariates"], allow_pickle=False),
        marker_ids=toy["marker_ids"],
        sample_ids=toy["sample_ids"],
        chunk_size=4,
        output_dir=linear_out,
    )
    summary = {
        "linear_rows": len(linear_result.table),
        "linear_top_hit": max(linear_result.table, key=lambda row: row["-log10_p"]),
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
    if args.command == "demo":
        return _run_demo(args)
    if args.command == "convert-pgen":
        return _run_convert_pgen(args)
    if args.command == "build-hardcall-store":
        return _run_build_hardcall_store(args)
    if args.command == "convert-bgen":
        genotype, _, _ = load_bgen_genotype(args.genotype, sample_file=args.sample_file,
            cache_dir=args.cache_dir, reader_workers=args.reader_workers)
        print(genotype.manifest_path.read_text())
        genotype.close()
        return 0
    raise ValueError(f"unknown command: {args.command}")
