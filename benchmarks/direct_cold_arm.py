"""One cold full-file run of one input format. One arm, one process.

A process per arm, not a loop inside one process, for two reasons that have
both bitten this project: the decoder library is memoized on first use, so two
builds of it cannot be compared inside one interpreter, and repeated scans in
one process used to degrade by 17.9x until the pinned-buffer lifetime was
fixed. A fresh process removes both questions from the measurement.

Prints one line of JSON so a driver can collect it without parsing prose.

**An output directory is required even when nothing is written.** Without one
`run_linear_gwas` returns the whole result instead of streaming it, which for
8.93M variants and 128 traits is not a benchmark of anything -- the first
attempt reached 146 GB of resident memory and was still going seven minutes
after the last byte of the input had been read.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

POSIX_FADV_DONTNEED = 4


def drop_cache(path: str) -> int:
    """Evict a file and its sidecars from the page cache.

    Returns the number of files evicted, so a caller can tell "cold" from
    "the sidecar was not where I thought it was".
    """
    libc = ctypes.CDLL("libc.so.6", use_errno=True)
    libc.posix_fadvise.argtypes = [ctypes.c_int, ctypes.c_int64,
                                   ctypes.c_int64, ctypes.c_int]
    stem, _ = os.path.splitext(path)
    candidates = [path, stem + ".pvar", stem + ".psam", stem + ".bim",
                  stem + ".fam", path + ".bgi", stem + ".zst", stem + ".json"]
    evicted = 0
    for candidate in candidates:
        if not os.path.exists(candidate):
            continue
        fd = os.open(candidate, os.O_RDONLY)
        try:
            os.fsync(fd)
            if libc.posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED) == 0:
                evicted += 1
        finally:
            os.close(fd)
    return evicted


def free_gpu() -> str:
    used = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used",
         "--format=csv,noheader,nounits"],
        capture_output=True, text=True).stdout.strip()
    busy = {int(line.split(",")[0]) for line in used.split("\n")
            if int(line.split(",")[1]) > 100}
    for index in range(8):
        if index not in busy:
            return f"cuda:{index}"
    raise SystemExit("every GPU on this host is already in use")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sumstats-format", default="none",
                        choices=("none", "binary", "tsv"))
    parser.add_argument("--pgen-mode", default=None,
                        choices=(None, "auto", "hardcall", "dosage"))
    parser.add_argument("--genotype-cache-dir", default=None)
    parser.add_argument("--library", default=None,
                        help="TORCHGWAS_PGEN_LIBRARY for this arm only")
    parser.add_argument("--zstd-read-workers", type=int, default=None)
    parser.add_argument("--hardcall-store", default=None,
                        help="prefix of a zstd hard-call store to read "
                             "instead of the .bed. Metadata still comes "
                             "from the .bim/.fam, so this is a transport "
                             "swap and nothing else.")
    parser.add_argument("--traits", type=int, default=128)
    parser.add_argument("--covariates", type=int, default=27)
    # No default: an explicit chunk_size OVERRIDES _default_chunk_variants,
    # which derives the chunk from free device memory. Forcing 20,000 here --
    # a value tuned at 22,250 samples -- is what OOM'd an 80 GB card at
    # 35,365 samples on the float32 dosage transport, 16x the bytes per sample
    # of the packed path. Leave it unset unless a run is deliberately pinning
    # the geometry.
    parser.add_argument("--chunk-size", type=int, default=None)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--prefetch", type=int, default=32)
    parser.add_argument("--variants", type=int, default=None,
                        help="stop after this many variants; whole file if unset")
    parser.add_argument("--samples", type=int, default=None,
                        help="subset to this many samples, for an N sweep; "
                             "the id lookup happens before the timer starts")
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--autotune", action="store_true",
                        help="autotune=True: the planner and tuner choose GPUs, layout, chunk and readers "
                             "(device, --chunk-size, --workers and --prefetch are not passed to the scan)")
    parser.add_argument("--keep-cache", action="store_true",
                        help="do not evict the input from the page cache first "
                             "(a warm arm, to separate storage from decode)")
    args = parser.parse_args()

    if args.library:
        os.environ["TORCHGWAS_PGEN_LIBRARY"] = args.library

    # A zstd store is named by a prefix, not by a file: `load_genotype` strips
    # a trailing `.zst` and finds the rest, but `getsize` and the cache
    # eviction below need a real path. Resolve it here so a caller can pass
    # either spelling.
    if not os.path.exists(args.path) and os.path.exists(args.path + ".zst"):
        args.path = args.path + ".zst"

    from torchgwas.api import run_linear_gwas
    from torchgwas.io import load_genotype

    evicted = 0 if args.keep_cache else drop_cache(args.path)
    # A store holds the bytes actually read, so evicting only the `.bed`
    # would leave the store warm and hand it a free win over the format
    # it is being compared against. This is the whole measurement.
    if args.hardcall_store and not args.keep_cache:
        evicted += drop_cache(args.hardcall_store + '.zhc')
    shutil.rmtree(args.output_dir, ignore_errors=True)
    load_before = float(open("/proc/loadavg").read().split()[0])

    loader_kwargs = dict(reader_workers=args.workers,
                         prefetch_chunks=args.prefetch)
    if args.pgen_mode:
        loader_kwargs["pgen_mode"] = args.pgen_mode
    if args.genotype_cache_dir:
        loader_kwargs["genotype_cache_dir"] = args.genotype_cache_dir
    if args.zstd_read_workers is not None:
        loader_kwargs["zstd_read_workers"] = args.zstd_read_workers
    if args.hardcall_store:
        loader_kwargs["hardcall_store"] = args.hardcall_store

    # Peak memory is a reported quantity, not a footnote: a tool that is fast
    # only because it took 90 GB of an 80 GB card is not fast, and the
    # calculator predicts `device_buffer_bytes` and `host_buffer_bytes` that
    # nothing has ever been checked against.
    import resource

    import torch

    # An N sweep needs a sample subset, and picking one needs the ids -- so
    # open once to read them, untimed, then time the real open. Charging the
    # lookup to the measurement would make every subset arm look slower than
    # the whole-cohort arm for a reason that has nothing to do with N.
    if args.samples is not None:
        probe, probe_ids, _pm, _pmeta = load_genotype(args.path, **loader_kwargs)
        chosen = [str(name) for name in np.asarray(probe_ids, dtype=str)[:args.samples]]
        del probe
        loader_kwargs["selected_sample_ids"] = chosen
        drop_cache(args.path)

    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()

    started = time.perf_counter()
    source, _, _, _ = load_genotype(args.path, **loader_kwargs)
    opened = time.perf_counter() - started

    rng = np.random.default_rng(args.seed)
    phenotypes = rng.normal(size=(source.shape[0], args.traits)).astype(np.float32)
    covariates = rng.normal(size=(source.shape[0], args.covariates))
    scan_kwargs = {}
    if args.variants is not None:
        scan_kwargs["variant_range"] = (0, args.variants)
    if args.autotune:
        scan_kwargs["autotune"] = True
    else:
        scan_kwargs.update(device=free_gpu(), chunk_size=args.chunk_size,
                           reader_workers=args.workers, prefetch_chunks=args.prefetch)
    result = run_linear_gwas(source, phenotypes, covariates=covariates,
                             compute_dtype="float32",
                             output_dir=args.output_dir,
                             sumstats_format=args.sumstats_format, **scan_kwargs)
    wall = time.perf_counter() - started
    # Did the variant range actually apply? An arm that ignores it scans the
    # whole file and returns a timing 45x too large for the M it claims -- and
    # nothing else in this record would show it, because `variants` below is
    # the file's variant count, not the scanned one. The zstd loader's range
    # path in particular used to raise TypeError and is newly fixed.
    rows_written = None
    metadata = getattr(result, "run_metadata", None)
    if isinstance(metadata, dict):
        rows_written = metadata.get("n_result_rows")

    output_bytes = 0
    for root, _dirs, files in os.walk(args.output_dir):
        for name in files:
            output_bytes += os.path.getsize(os.path.join(root, name))

    # Report the path actually taken, not the one the flags asked for. Two
    # nested opt-ins decide it (TORCHGWAS_NATIVE_STATS and
    # TORCHGWAS_PGEN_PACKED), a source silently keeps the int8 transport when
    # it is not eligible, and a whole afternoon of PGEN measurements turned out
    # to have been taken on a transport nobody intended.
    print(json.dumps({
        "label": args.label,
        "autotune": bool(args.autotune),
        "autotune_layout": ((metadata or {}).get("autotune") or {}).get("layout") if args.autotune else None,
        "autotune_chunk": {k: v for k, v in (((metadata or {}).get("autotune") or {}).get("chunk") or {}).items()
                           if k in ("choice", "state", "reason", "reprobes", "decided_fraction")} if args.autotune else None,
        "devices_used": (metadata or {}).get("variant_devices") or (metadata or {}).get("trait_devices")
                        or (metadata or {}).get("device_used"),
        "sumstats_format": args.sumstats_format,
        "encoding": getattr(source, "native_encoding", None),
        "row_bytes": (getattr(source, "native_row_width", 0)
                      * getattr(getattr(source, "native_transfer_dtype", None),
                                "itemsize", 1)) or None,
        "statistics_backend": getattr(source, "_last_scan_profile", {}).get(
            "statistics_backend"),
        # The decoder actually used. A GPU library that fails to load (ABI
        # mismatch, unbuilt) falls back to the CPU decoder under the default
        # 'auto' -- 30x slower and indistinguishable from a slow GPU arm
        # without this field.
        "decode_backend": getattr(source, "backend_used", None),
        "wall": round(wall, 3),
        "scan": round(wall - opened, 3),
        "open": round(opened, 3),
        "traits": int(args.traits),
        # The sample and variant counts the run ACTUALLY used, not the ones
        # requested. bed and zstd used to accept a `--samples` subset and scan
        # the whole cohort anyway, so an N sweep could report four identical
        # runs as four different values of N; a sweep that does not print what
        # it got cannot detect that happening again.
        "samples": int(source.shape[0]),
        "samples_requested": args.samples,
        "variants_requested": args.variants,
        "variants_available": int(source.shape[1]),
        "variants_scanned_expected": int(args.variants or source.shape[1]),
        "n_result_rows": rows_written,
        # rows / K should equal the M asked for, give or take the handful of
        # monomorphic or all-missing variants QC drops. Far above it means the
        # range was ignored.
        "rows_per_trait": (None if not rows_written
                           else rows_written / max(1, int(args.traits))),
        "chunk_size": args.chunk_size,
        # The geometry the scan RESOLVED, which is what a prediction must be
        # checked against. `chunk_size` above is only what was requested, and
        # it is normally None because the planner chooses.
        "resolved_chunk": getattr(source, "_last_scan_profile", {}).get(
            "chunk_variants"),
        "resolved_depth": getattr(source, "_last_scan_profile", {}).get("depth"),
        "transfer_bytes_per_variant": getattr(
            source, "_last_scan_profile", {}).get("transfer_bytes_per_variant"),
        "decode_on_gpu": getattr(source, "_last_scan_profile", {}).get(
            "decode_on_gpu"),
        # `max_memory_allocated` is what the scan asked for; `reserved` is what
        # the caching allocator actually held, and that is the number deciding
        # whether anything else fits on the card. ru_maxrss is host peak, KB on
        # Linux. The calculator predicts `device_buffer_bytes` and
        # `host_buffer_bytes`; nothing had ever been checked against them.
        "peak_gpu_allocated_bytes": (int(torch.cuda.max_memory_allocated())
                                     if torch.cuda.is_available() else None),
        "peak_gpu_reserved_bytes": (int(torch.cuda.max_memory_reserved())
                                    if torch.cuda.is_available() else None),
        "peak_host_rss_bytes": int(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024,
        "samples": int(source.shape[0]),
        "variants": int(source.shape[1]),
        "input_bytes": os.path.getsize(args.path),
        "output_bytes": output_bytes,
        "evicted_files": evicted,
        "load_before": load_before,
        "load_after": float(open("/proc/loadavg").read().split()[0]),
        # Only when TORCHGWAS_BGEN_PROFILE=1 asked for it: the decoder's
        # per-stage CUDA-event and host-stage accounting for this scan.
        **({"decode_profile": getattr(source, "last_decode_profile", None)}
           if os.environ.get("TORCHGWAS_BGEN_PROFILE", "0") != "0" else {}),
        # Likewise the consumer's accounting (fetch wait, GPU compute ms)
        # under TORCHGWAS_SCAN_PROFILE=1.
        **({"scan_profile": getattr(source, "_last_scan_profile", None)}
           if os.environ.get("TORCHGWAS_SCAN_PROFILE", "0") != "0" else {}),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
