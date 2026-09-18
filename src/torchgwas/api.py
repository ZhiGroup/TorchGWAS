from __future__ import annotations

import json
import time
import warnings
from pathlib import Path


import numpy as np
import torch

from .io import align_table_to_samples, load_array, load_genotype, load_vector
from .linear import linear_scan, linear_scan_streaming, linear_scan_streaming_chunks
from .pgen import (
    DEFAULT_PGEN_COMPRESSION_WORKERS,
    DEFAULT_PGEN_DECODE_BATCH_SIZE,
    DEFAULT_PGEN_DECODE_WORKERS,
)
from .preprocess import prepare_inputs, prepare_inputs_for_prep
from .streaming import ChunkedGenotype, _resolve_variant_range
from .types import GWASResult
from .utils import choose_device, elapsed, mkdir, timestamp, write_json


DEFAULT_SUMSTATS_BLOCK_BYTES = 16 << 20
DEFAULT_SUMSTATS_QUEUE_DEPTH = 3


def _phase_breakdown(entered, prep_started, prep_done, write_started, finished):
    """Attribute the run's wall clock to named phases.

    `runtime_seconds` and `sumstats_write.scan_and_write_seconds` were the only
    two clocks a run reported, so everything between them was invisible. On the
    600,000-trait stress run that gap was 1,804 of 4,183 seconds -- 43% of the
    wall, spent on residualising an 85 GB phenotype before a single variant was
    read, and nothing in the output said so.

    Phases are reported only when their boundaries were actually reached; a run
    that returns an in-memory table never sets `write_started`, and reporting
    that as zero seconds would be a lie rather than a gap.
    """
    phases: dict[str, float] = {}
    if prep_started is not None:
        phases["open_and_resolve"] = prep_started - entered
    if prep_started is not None and prep_done is not None:
        phases["input_qc"] = prep_done - prep_started
    if prep_done is not None and write_started is not None:
        # Residualisation, trait blocking and scan construction. Everything
        # between the inputs being clean and the first byte being read.
        phases["scan_setup"] = write_started - prep_done
    if write_started is not None:
        phases["scan_and_write"] = finished - write_started
    phases["total"] = finished - entered
    return phases or None



def _available_host_bytes() -> int:
    """Host memory a pinned ring may claim, or 0 when it cannot be determined.

    MemAvailable, not MemTotal: pinned pages cannot be swapped or reclaimed,
    so the budget is what is free right now on a box that is shared with other
    jobs, not what the machine was built with.
    """
    try:
        with open("/proc/meminfo", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return 0


def _coerce_array_or_path(value):
    if value is None:
        return None
    if isinstance(value, (str, Path)):
        return load_array(value)
    return np.asarray(value)


def _coerce_vector_or_path(value):
    if value is None:
        return None
    if isinstance(value, (str, Path)):
        if Path(value).suffix.lower() == ".npy":
            return np.asarray(np.load(value, allow_pickle=False))
        return load_vector(value)
    return np.asarray(value)


def _prefer_vector(primary, fallback):
    return primary if primary is not None else fallback


def _select_chunked_samples(genotype, genotype_sample_ids, requested_sample_ids):
    if requested_sample_ids is None or genotype_sample_ids is None:
        return genotype, genotype_sample_ids, requested_sample_ids
    requested = np.asarray([str(value) for value in requested_sample_ids], dtype=object)
    stored = np.asarray([str(value) for value in genotype_sample_ids], dtype=object)
    if hasattr(genotype, "select_samples"):
        genotype.select_samples(requested)
        selected = np.asarray(genotype.sample_ids, dtype=object)
        return genotype, selected, selected
    if not np.array_equal(requested, stored):
        raise ValueError(
            f"{type(genotype).__name__} does not support sample subsetting; "
            "pre-align or convert the genotype input first"
        )
    return genotype, genotype_sample_ids, stored


def _resolve_linear_compute_dtype(genotype, compute_dtype: str) -> str:
    if compute_dtype != "auto":
        return compute_dtype
    return "float32" if isinstance(genotype, ChunkedGenotype) else "float64"






_QUOTE_TRIGGERS = ("\t", "\n", "\r", '"')










def _accumulate(accumulated, order, reduction, chunk, offset, guard=None):
    """Merge one reduced chunk into the running per-variant accumulator."""
    start, end, beta, t, p, index = chunk

    def as_tensor(array):
        return (None if array is None
                else torch.from_numpy(np.ascontiguousarray(array)))

    incoming = (as_tensor(beta), as_tensor(t), as_tensor(p), as_tensor(index),
                # The status word is not part of the chunk contract and does not
                # need to be: an invalid variant arrives as NaN in every column,
                # and NaN is scrubbed to -inf before the merge picks, so it can
                # never win a row.
                torch.zeros(end - start, dtype=torch.uint8),
                torch.zeros(end - start))
    key = (start, end)
    if guard is None:
        if key not in accumulated:
            order.append(key)
        accumulated[key] = reduction.merge(accumulated.get(key), incoming, offset)
        return
    with guard:
        if key not in accumulated:
            order.append(key)
        accumulated[key] = reduction.merge(accumulated.get(key), incoming, offset)


def _trait_blocked_significant_chunks(scan_once, significance, n_traits,
                                      trait_block, df, devices=None):
    """Significant (variant, trait) pairs for all K traits, block by block.

    `reduce='significant'` cannot use `_trait_blocked_reduced_chunks`: that
    one merges a `VariantReduction`'s one-row-per-variant result and expects
    a sixth tuple element carrying the winning trait index. Significance has
    no winner and no merge -- it emits every pair over a threshold, so the
    number of rows per chunk varies and blocks simply concatenate. Routing
    it through the reduced driver raised `not enough values to unpack
    (expected 6, got 5)` the moment blocking first engaged.

    **The threshold is Bonferroni over the WHOLE trait axis, never the
    block.** A block is an implementation detail chosen from the available
    memory, and if it changed which pairs were significant then the results
    would depend on the size of the card they ran on. `n_traits` here is
    always the full K.

    Trait indices come back numbered within their block and are rebased onto
    the full axis before they are yielded, which is the only bookkeeping the
    concatenation needs.
    """
    from .linear import _significant_pairs_iterator

    if trait_block < 1:
        raise ValueError("trait_block must be positive")
    blocks = [(offset, min(trait_block, n_traits - offset))
              for offset in range(0, n_traits, trait_block)]

    def emit(offset, width, device):
        pairs = _significant_pairs_iterator(
            scan_once(offset, width, device), significance, n_traits, df)
        for start, end, variant_index, trait_index, beta, t_stat, row_df in pairs:
            yield (start, end, variant_index,
                   trait_index.astype(np.int64) + offset,
                   beta, t_stat, row_df)

    if not devices or len(devices) <= 1:
        for offset, width in blocks:
            yield from emit(offset, width, None)
        return

    # One thread per device, each taking its own blocks. Rows are handed
    # back through a bounded queue rather than accumulated: the whole point
    # of this mode is that the output is small but unbounded in principle,
    # so it streams to the writer instead of being held.
    import queue
    import threading

    results: queue.Queue = queue.Queue(maxsize=4 * len(devices))
    failures: list[BaseException] = []
    done = object()

    def run(device, assigned):
        try:
            for offset, width in assigned:
                for item in emit(offset, width, device):
                    results.put(item)
        except BaseException as exc:  # noqa: BLE001 - re-raised below
            failures.append(exc)
        finally:
            results.put(done)

    workers = []
    for index, device in enumerate(devices):
        assigned = blocks[index::len(devices)]
        if not assigned:
            continue
        thread = threading.Thread(target=run, args=(device, assigned),
                                  name=f"torchgwas-sigshard-{index}",
                                  daemon=True)
        thread.start()
        workers.append(thread)
    remaining = len(workers)
    while remaining:
        item = results.get()
        if item is done:
            remaining -= 1
            continue
        yield item
    for thread in workers:
        thread.join()
    if failures:
        raise failures[0]

def _trait_blocked_reduced_chunks(scan_once, reduction, n_traits, trait_block,
                                  devices=None):
    """Reduced chunks for all K traits, processing the traits in blocks.

    **Why this exists.** The reduction fixed the output side of a high-trait
    scan; this fixes the input side. The residualised phenotype is resident on
    the device for the whole scan, and at 22,250 samples by 2.085M traits that
    is **185 GB** against an 80 GB card. Blocking the traits is the only way
    such a scan runs at all.

    **It is exact, not an approximation.** Residualisation is column-separable
    (the projection is `basis @ (basis.T @ values)` and the standardisation is
    per column), and the residual degrees of freedom are per *variant*, so a
    trait block computes precisely what the full matrix would. Merging keeps
    the true top k because the top k of a union is contained in the union of
    the two top-k sets.

    **Loop order: blocks outside chunks is the right way round, and the
    intuition that says otherwise is wrong.** This reads the genotypes once per
    block, which looks wasteful -- 47 passes for the 2M-trait voxel workload --
    so the obvious improvement is to put the block loop *inside* the chunk loop
    and read the genotypes once in total. That is worse, by a lot, because it
    re-uploads the **design** once per chunk instead, and at high K the design is
    the large object: 22,250 x 2.085M float32 is 185 GB, and 1,091 chunks of it
    is **202 TB against 2.5 TB** for this order. At K = 32,768 it is still 8.1x
    worse. Writing out the traffic:

        outer (this):  B * M * samples * geno_bytes  +  samples * K * 4
        inner:         M * samples * geno_bytes      +  C * samples * K * 4

    with B blocks and C chunks. For two-bit genotypes the inner form only wins
    when `block < chunk / 16` -- under 512 traits at a chunk of 8,192 -- and
    blocks that small give GEMMs too narrow to be worth having.

    **Sharding the blocks across GPUs.** Pass `devices` and each device takes
    its own trait blocks, concurrently, merging into the same accumulator. This
    is the axis `linear_scan_multigpu` cannot reach: that shards by *variant
    range* and hands every shard the same phenotype, so each device still needs
    the full K-wide design -- at 185 GB, variant sharding does not help the
    memory wall at all. Trait sharding divides both the memory and the block
    loop. The merge is safe to run from several threads because taking the top
    k of a union does not depend on the order the union was assembled in; a lock
    serialises only the accumulator update itself.

    `scan_once(offset, width, device)` runs one reduced scan over
    `phenotype[:, offset:offset + width]` and yields its
    `(start, end, beta, t, p, index)` chunks.
    """
    if trait_block < 1:
        raise ValueError("trait_block must be positive")
    accumulated: dict[tuple[int, int], tuple] = {}
    order: list[tuple[int, int]] = []
    blocks = [(offset, min(trait_block, n_traits - offset))
              for offset in range(0, n_traits, trait_block)]

    if devices and len(devices) > 1:
        import threading

        guard = threading.Lock()
        failures: list[BaseException] = []

        def run(device, assigned):
            try:
                for offset, width in assigned:
                    for chunk in scan_once(offset, width, device):
                        _accumulate(accumulated, order, reduction, chunk,
                                    offset, guard)
            except BaseException as exc:  # noqa: BLE001 - re-raised below
                with guard:
                    failures.append(exc)

        workers = []
        for index, device in enumerate(devices):
            assigned = blocks[index::len(devices)]
            if not assigned:
                continue
            thread = threading.Thread(target=run, args=(device, assigned),
                                      name=f"torchgwas-traitshard-{index}",
                                      daemon=True)
            thread.start()
            workers.append(thread)
        for thread in workers:
            thread.join()
        if failures:
            raise failures[0]
        for key in sorted(order):
            start, end = key
            beta, t, p, index, _status, _df = accumulated[key]
            yield (start, end, beta.numpy(), t.numpy(),
                   None if p is None else p.numpy(), index.numpy())
        return

    for offset, width in blocks:
        for chunk in scan_once(offset, width, None):
            _accumulate(accumulated, order, reduction, chunk, offset)
    for key in order:
        start, end = key
        beta, t, p, index, _status, _df = accumulated[key]
        yield (start, end, beta.numpy(), t.numpy(),
               None if p is None else p.numpy(), index.numpy())













def _write_linear_binary_streaming(
    directory: str | Path,
    marker_names: list[str],
    trait_names: list[str],
    n_samples: int,
    chunk_iterator,
    n_variants: int,
    df: int | list[int],
    block_bytes: int | None,
    queue_depth: int,
    fsync: bool,
    write_variant_ids: bool,
    store_beta: bool = True,
    extra_manifest: dict | None = None,
    borrow_results: bool = False,
) -> tuple[int, dict]:
    """Stream beta/t_stat/-log10(P) chunks into a binary directory.

    Returns the cell count and a timing summary. The iterator is consumed in
    variant order; storage back-pressure reaches the scan through the writer's
    bounded block ring.
    """
    from .sumstats import BinarySumstatsWriter

    directory = Path(directory)
    writer = BinarySumstatsWriter(
        directory=directory,
        n_variants=n_variants,
        trait_names=trait_names,
        n_samples=n_samples,
        df=df,
        block_bytes=block_bytes,
        queue_depth=queue_depth,
        fsync=fsync,
        store_beta=store_beta,
        extra_manifest=extra_manifest or {},
        # False when the scan is handing out views into its result ring: the
        # writer must then copy into staging rather than queueing the caller's
        # array to its writer thread, because that array is refilled as soon as
        # the next chunk is asked for. The two settings do the same number of
        # memcpys of the result stream -- the scan's copy simply moves into the
        # writer's staging buffer, which is reused rather than freshly
        # allocated per chunk.
        borrow_chunks=not borrow_results,
    )
    try:
        for start, end, beta_chunk, t_chunk, _p_chunk, logp_chunk in chunk_iterator:
            writer.write_chunk(start, end, beta_chunk, t_chunk, logp_chunk)
    except BaseException:
        writer.abort()
        raise
    summary = writer.close()
    if write_variant_ids:
        started = time.perf_counter()
        written = _write_variant_ids(directory / "variant_ids.txt",
                                     marker_names[:n_variants])
        summary["variant_id_seconds"] = time.perf_counter() - started
        summary["variant_id_bytes"] = written
    return summary["cells"], summary


def _write_variant_ids(path, marker_names) -> int:
    """Write one marker id per line, and do it without rebuilding every string.

    **This sidecar is the fixed cost of writing sumstats at all.** Measured on
    the full genome, enabling binary sumstats at K=1 costs **+3.72 s** for a
    0.17 GB payload that should take 0.035 s at the measured write rate. The
    difference is here: 8,931,083 ids, flat in K and proportional to M.

    Timed on the real marker list (8,931,083 rows, 96.4 MB out):

        "\\n".join(map(str, names))        3.96 s   <- what this was
        "\\n".join(names)                  3.46 s
        chunked writes, no big string     3.14 s
        bytes join + binary write         2.24 s   <- what this is
        np.savetxt                        7.02 s

    `marker_ids` is already a numpy `<U16` array, so `map(str, ...)` built 8.9
    million Python strings that already existed in another form. Joining bytes
    and writing binary skips both that and the encode of a 96 MB `str`.

    Falls back to the text path for anything the byte path cannot represent --
    a non-ASCII id is unusual but not impossible, and silently mangling marker
    names to save a second would be a very bad trade.
    """
    import numpy as np

    try:
        blob = b"\n".join(
            np.asarray(marker_names).astype("S").tolist()) + b"\n"
    except (UnicodeEncodeError, ValueError, TypeError):
        # Non-ASCII, ragged, or not array-like: correctness wins.
        #
        # **UTF-8 explicitly.** The original line was `write_text(payload)`,
        # which encodes with the locale default -- ASCII on the lab host. So a
        # single non-ASCII marker id would raise `UnicodeEncodeError` there and
        # kill the run AFTER the whole scan had finished, at the last step, on
        # some machines and not others. That bug predates this fast path; the
        # test for the fallback is what exposed it.
        payload = "\n".join(map(str, marker_names)) + "\n"
        encoded = payload.encode("utf-8")
        Path(path).write_bytes(encoded)
        return len(encoded)
    Path(path).write_bytes(blob)
    return len(blob)


def _drain_linear_chunks(chunk_iterator) -> tuple[int, dict]:
    """Consume the scan without serializing results.

    This isolates scan cost from output cost so both can be reported; it is a
    measurement mode, not an analysis mode, and produces no result table.
    """
    cells = 0
    started = time.perf_counter()
    # Tuple-length agnostic: a reduced scan yields a sixth element, and this
    # mode has to work there too -- measuring scan cost apart from output cost
    # is exactly how the reduction's benefit is separated from the writer's.
    for chunk in chunk_iterator:
        if len(chunk) == 7:  # selected variant/trait pairs
            cells += len(chunk[2])
        else:
            cells += int(np.asarray(chunk[3]).size)
    return cells, {"discarded": True, "scan_seconds": time.perf_counter() - started}


def run_linear_gwas(
    genotype,
    phenotype,
    covariates=None,
    phenotype_table=None,
    covariates_table=None,
    trait_columns: list[str] | None = None,
    covariate_columns: list[str] | None = None,
    genotype_format: str = "auto",
    sample_file: str | Path | None = None,
    bgen_decode_backend: str = "auto",
    # Read genotype bytes from a zstd hard-call store instead of the `.bed`.
    # The `.bim`/`.fam` are still used, so this is a transport swap and cannot
    # change which variants or samples are analysed. `load_genotype` REFUSES
    # it for any non-PLINK format rather than ignoring it.
    hardcall_store: str | None = None,
    pvar: str | Path | None = None,
    psam: str | Path | None = None,
    pgen_mode: str = "auto",
    pgen_decode_workers: int | None = None,
    pgen_decode_batch_size: int = DEFAULT_PGEN_DECODE_BATCH_SIZE,
    pgen_compression_workers: int = DEFAULT_PGEN_COMPRESSION_WORKERS,
    genotype_cache_dir: str | Path | None = None,
    bim: str | Path | None = None,
    fam: str | Path | None = None,
    plink2_binary: str | Path | None = None,
    reader_workers: int | None = None,
    zstd_read_workers: int | None = None,
    prefetch_chunks: int | None = None,
    sample_id_column: str = "IID",
    sample_ids=None,
    marker_ids=None,
    device: str = "auto",
    compute_dtype: str = "auto",
    chunk_size: int | None = None,
    topk_per_trait: int | None = None,
    p_value_threshold: float | None = None,
    reduce: str | None = None,
    significance_threshold: float | None = None,
    reduce_top_k: int | None = None,
    # Internal, for the machinery's own correctness tests. See the guard in the
    # body: `reduce=` accepts only 'significant' and 'jagwas', and the
    # per-variant top-k reductions they are built on are reached through here
    # so their tests can drive the real pipeline without reopening the public
    # surface.
    _internal_reduction=None,
    trait_block: int | None = None,
    trait_devices=None,
    return_beta: bool = True,
    return_se: bool = True,
    return_t: bool = True,
    output_dir: str | Path | None = None,
    variant_range: tuple[int, int] | None = None,
    pipeline_profile: dict | str | Path | None = None,
    sumstats_format: str = "binary",
    sumstats_block_bytes: int | None = None,
    sumstats_queue_depth: int | None = None,
    sumstats_fsync: bool = True,
    sumstats_variant_ids: bool = True,
    sumstats_fields: str = "beta+t",
) -> GWASResult:
    if sumstats_format not in {"binary", "none"}:
        raise ValueError("sumstats_format must be 'binary' or 'none'; TSV output has been removed")
    sumstats_summary: dict = {}
    requested_reader_workers = reader_workers
    requested_prefetch_chunks = prefetch_chunks
    reader_workers = 24 if reader_workers is None else reader_workers
    prefetch_chunks = 4 if prefetch_chunks is None else prefetch_chunks
    start = timestamp()
    # Phase boundaries, left as None on the paths that never reach them so a
    # missing phase reads as absent rather than as zero seconds.
    _phase_entered = time.perf_counter()
    _phase_prep_started = _phase_prep_done = write_started = None
    resolved_device = choose_device(device)
    genotype_meta = {}
    requested_sample_ids = _coerce_vector_or_path(sample_ids)
    if isinstance(genotype, (str, Path)):
        genotype, geno_sample_ids, geno_marker_ids, genotype_meta = load_genotype(
            genotype,
            genotype_format=genotype_format,
            bim=bim,
            fam=fam,
            sample_file=sample_file,
            bgen_decode_backend=bgen_decode_backend,
            hardcall_store=hardcall_store,
            pvar=pvar,
            psam=psam,
            selected_sample_ids=requested_sample_ids,
            genotype_cache_dir=genotype_cache_dir,
            plink2_binary=plink2_binary,
            reader_workers=reader_workers,
            zstd_read_workers=zstd_read_workers,
            prefetch_chunks=prefetch_chunks,
            pgen_mode=pgen_mode,
            pgen_decode_workers=pgen_decode_workers,
            pgen_decode_batch_size=pgen_decode_batch_size,
            pgen_compression_workers=pgen_compression_workers,
        )
    elif isinstance(genotype, ChunkedGenotype):
        geno_sample_ids, geno_marker_ids = genotype.sample_ids, genotype.marker_ids
    else:
        if pipeline_profile is not None:
            raise ValueError("pipeline profiles require a streaming genotype source")
        genotype = np.asarray(genotype)
        geno_sample_ids, geno_marker_ids = None, None
    if isinstance(genotype, ChunkedGenotype):
        genotype, geno_sample_ids, requested_sample_ids = _select_chunked_samples(
            genotype,
            geno_sample_ids,
            requested_sample_ids,
        )
    marker_ids = _prefer_vector(_coerce_vector_or_path(marker_ids), geno_marker_ids)
    sample_ids = _prefer_vector(requested_sample_ids, geno_sample_ids)
    variant_metadata = getattr(genotype, "variant_metadata", None)
    if phenotype_table is not None:
        if sample_ids is None:
            raise ValueError("tabular phenotype input requires genotype sample IDs")
        phenotype, trait_columns = align_table_to_samples(
            phenotype_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=trait_columns, sample_id_column=sample_id_column
        )
    else:
        phenotype = _coerce_array_or_path(phenotype)
    if covariates_table is not None:
        if sample_ids is None:
            raise ValueError("tabular covariate input requires genotype sample IDs")
        covariates, covariate_columns = align_table_to_samples(
            covariates_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=covariate_columns, sample_id_column=sample_id_column
        )
    else:
        covariates = _coerce_array_or_path(covariates)
    resolved_compute_dtype = _resolve_linear_compute_dtype(genotype, compute_dtype)
    if sumstats_fields not in {"beta+t", "t"}:
        raise ValueError(
            f"sumstats_fields must be 'beta+t' or 't', got {sumstats_fields!r}"
        )
    if sumstats_fields == "t" and sumstats_format != "binary":
        raise ValueError(
            "sumstats_fields='t' applies only to sumstats_format='binary'"
        )
    if sumstats_format == "none" and (topk_per_trait is not None or p_value_threshold is not None):
        raise ValueError("row selection requires binary output")
    if topk_per_trait is not None and topk_per_trait <= 0:
        raise ValueError("topk_per_trait must be positive")
    if p_value_threshold is not None and not (0.0 < p_value_threshold <= 1.0):
        raise ValueError("p_value_threshold must be in (0, 1]")
    # A reduction changes what the scan produces, not merely which rows are
    # kept, so every combination that assumes a full (variant x trait) matrix is
    # refused here rather than silently producing a narrower one.
    reduction = None
    significance = None
    jagwas = None
    if reduce == "jagwas":
        from .reduce import JagwasReduction

        if output_dir is None:
            raise ValueError(
                "reduce='jagwas' requires output_dir: it is a streaming "
                "reduction and the in-memory path returns the full matrix")
        if reduce_top_k is not None:
            raise ValueError("reduce_top_k does not apply to 'jagwas'")
        if trait_block is not None:
            raise ValueError(
                "jagwas cannot be trait-blocked: the statistic is a quadratic "
                "form over the whole trait correlation, so a block of traits "
                "does not carry enough information to be merged")
        # Which makes jagwas feasible only while the WHOLE phenotype fits --
        # it needs the full K-wide residualised matrix and a K x K correlation
        # at once, and neither can be tiled. Say so here, with the arithmetic,
        # rather than letting a voxel-scale K fail somewhere inside the scan:
        # at K = 2,085,000 the correlation alone is 17 TB.
        traits = int(np.asarray(phenotype).shape[1])
        samples = int(np.asarray(phenotype).shape[0])
        correlation_bytes = traits * traits * 4
        residual_bytes = samples * traits * 4
        try:
            import torch as _torch

            budget = (int(_torch.cuda.mem_get_info(_torch.device(device))[0])
                      if str(device).startswith("cuda") and _torch.cuda.is_available()
                      else None)
        except Exception:  # noqa: BLE001 - a probe failure must not decide this
            budget = None
        if budget is not None and correlation_bytes + residual_bytes > budget:
            raise ValueError(
                f"reduce='jagwas' needs the whole phenotype at once: a "
                f"{traits} x {traits} trait correlation is "
                f"{correlation_bytes / 1e9:.1f} GB and the residualised "
                f"phenotype is {residual_bytes / 1e9:.1f} GB, against "
                f"{budget / 1e9:.1f} GB free on {device}. The statistic is a "
                f"quadratic form over all traits, so it cannot be tiled -- use "
                f"reduce='significant', which streams and is what scales to "
                f"this many traits")
        jagwas = JagwasReduction()
        reduction = jagwas
        reduce = None
    if reduce == "significant":
        from .reduce import SignificantPairs

        if output_dir is None:
            raise ValueError(
                "reduce='significant' requires output_dir: it is a streaming "
                "selection, and the in-memory path returns the full matrix it "
                "exists to avoid")
        if reduce_top_k is not None:
            raise ValueError("reduce_top_k does not apply to 'significant'")
        significance = SignificantPairs(significance_threshold)
        # Run the proven fast path and select **after** it, not inside it.
        # Two earlier designs are recorded here because each was wrong in its
        # own way. Routing significance to the generic device path measured
        # **307x slower** (1362.75 s against 4.44 s on the same 200,000-variant
        # scan). Keeping the top `k` traits per variant and thresholding
        # afterwards was fast but **not correct**: exact only while fewer than
        # `k` traits at a variant clear the threshold, and a real association
        # in a correlated trait set can clear hundreds at once, so the file
        # quietly lost rows with a warning standing in for a guarantee.
        # There is now no `k` anywhere in this mode.
        reduce = None
    if _internal_reduction is not None:
        # INTERNAL ONLY, and deliberately awkward to reach.
        #
        # The per-variant top-k machinery still has to be correct: significant
        # pairs is built on it, and the trait-blocked paths merge through it.
        # Its correctness tests are worth keeping and are only meaningful when
        # driven through the REAL pipeline -- "blocked matches unblocked" says
        # nothing if it bypasses the scan. So the tests hand a pre-built
        # `VariantReduction` in here rather than going through `reduce=`, which
        # stays closed to callers.
        if reduce is not None:
            raise ValueError("pass reduce= or _internal_reduction=, not both")
        if output_dir is None:
            raise ValueError("_internal_reduction requires output_dir")
        reduction = _internal_reduction
    elif reduce is not None:
        # THE REDUCTION IS TWO MODES: 'significant' and 'jagwas'. Both are
        # handled above and set `reduce` back to None, so anything still set
        # here is one of the per-variant top-k spellings ('max-abs-t',
        # 'max-t2', 'min-p', 'top-k'). Those remain as MACHINERY --
        # SignificantPairs is built on top-k, and `VariantReduction` is still
        # used internally by the trait-blocked paths -- but they are not a
        # user-facing answer and are no longer accepted here.
        #
        # This is a standing decision, and every stress run today broke it by
        # passing reduce='max-abs-t'. Refusing it at the entry point is what
        # stops that recurring: a decision that lives only in a document gets
        # violated by whoever did not reread the document, which this time
        # was me.
        raise ValueError(
            f"reduce={reduce!r} is not a user-facing mode. The reduction has "
            "two: reduce='significant' (significant pairs; the threshold "
            "defaults to 5e-8/K) and reduce='jagwas'. The per-variant top-k "
            "spellings are internal machinery that significant pairs is built "
            "on -- if you want one trait per variant, that is a significance "
            "threshold, not a separate mode.")
    if reduce_top_k is not None:
        raise ValueError(
            "reduce_top_k does not apply: the user-facing reductions are "
            "'significant' and 'jagwas', and neither takes a k")
    if trait_block is not None:
        # Trait blocking only makes sense with a reduction: without one the
        # scan emits every (variant, trait) cell anyway, so nothing is saved
        # and the blocks would just be re-read genotypes.
        if reduction is None and significance is None:
            raise ValueError(
                "trait_block requires reduce: blocking the traits only helps "
                "when the result is reduced across them, otherwise every cell "
                "is emitted regardless and the blocks cost extra passes")
        if int(trait_block) < 1:
            raise ValueError("trait_block must be positive")
    elif ((reduction is not None or significance is not None)
          and jagwas is None and str(device).startswith("cuda")):
        # `significance is not None` matters as much as `reduction`, and
        # leaving it out meant the ONE mode the voxel stress test uses was
        # the one mode that never blocked. `reduce='significant'` builds a
        # `SignificantPairs` in its own variable and sets `reduction` back to
        # None, so this branch simply did not fire: a 600,000-voxel run put a
        # 79.4 GB design on GPU 0 and left the other seven cards idle, with
        # no error, because the design happened to fit in 85 GB. It would
        # not have fit at 2,085,000 voxels.
        # `jagwas is None` matters: jagwas is a quadratic form over the whole
        # trait correlation and cannot be blocked at all, so auto-blocking it
        # would silently produce a statistic that cannot be merged. Its own
        # feasibility check above refuses it when the phenotype does not fit.
        #
        # Choose the trait geometry rather than making the caller do it. At
        # voxel scale K is the wall, not the variant count: the design matrix
        # alone is `n_samples * K * 4`, which at 33,417 subjects and 2,085,000
        # voxels is 279 GB and fits on no device at any chunk size. The block
        # width comes from the same `device_ring_bytes` model that sizes the
        # chunk, so one formula governs both axes -- and when the whole matrix
        # already fits, the block is K and nothing changes.
        try:
            import torch as _torch

            from .pipeline_model import auto_trait_block

            if _torch.cuda.is_available():
                free_bytes, _total = _torch.cuda.mem_get_info(_torch.device(device))
                width = getattr(genotype, "native_row_width", 0)
                if (getattr(genotype, "native_encoding", None)
                        in {"pgen_2bit", "plink_2bit"} and width):
                    per_variant = float(width)
                elif hasattr(genotype, "iter_packed_chunks"):
                    per_variant = float((int(genotype.shape[0]) + 3) // 4)
                else:
                    per_variant = float(genotype.shape[0]) * 4.0
                derived = auto_trait_block(
                    n_samples=int(genotype.shape[0]),
                    n_traits=int(phenotype.shape[1]),
                    covariate_rank=int(0 if covariates is None
                                       else np.asarray(covariates).shape[1]),
                    chunk_variants=int(chunk_size or 4096),
                    depth=max(1, int(prefetch_chunks or 4)),
                    transfer_bytes_per_variant=per_variant,
                    device_memory_bytes=int(free_bytes),
                    decode_on_gpu=hasattr(genotype, "iter_device_chunks"),
                    # The host is the ceiling that actually binds at voxel K.
                    # Every ring slot is pinned on both sides, and the result
                    # ring is chunk-by-block rather than samples-by-block, so
                    # the device-only bisection chose blocks needing twice the
                    # host's memory in pinned -- therefore unswappable -- pages.
                    host_memory_bytes=_available_host_bytes(),
                    reduction_width=(None if reduction is None else
                                     reduction.resolved_width(
                                         int(phenotype.shape[1]))),
                    trait_devices=(len(trait_devices) if trait_devices
                                   else max(_torch.cuda.device_count(), 1)))
                # Largest-that-fits minimises the number of blocks, which is
                # the wrong objective once the blocks run concurrently: at
                # 2,085,000 voxels the widest fitting block is 518,821, which
                # is 5 blocks on an 8-card host -- three cards idle and a
                # 1.6x longer wall clock than spreading the same work eight
                # ways. Narrow the block until there is at least one per card.
                cards = (len(trait_devices) if trait_devices
                         else max(_torch.cuda.device_count(), 1))
                # ONLY when blocking is already required. Narrowing a block
                # that already fits turns a one-pass scan into one pass per
                # card, and each pass re-reads the whole genotype: measured
                # at 150,000 voxels, which fit a single card, spreading over
                # eight made the scan 129.1 s -> ~254 s. Sharding pays only
                # for work that had to be split anyway.
                if cards > 1 and int(derived) < int(phenotype.shape[1]):
                    spread = -(-int(phenotype.shape[1]) // cards)
                    derived = min(int(derived), spread)
                if derived < int(phenotype.shape[1]):
                    trait_block = int(derived)
                    # Every visible GPU unless the caller named some. Trait
                    # blocks are independent and each carries its own slice of
                    # the design, so this is the axis that shards without
                    # duplicating the K-wide matrix onto every device.
                    if not trait_devices and _torch.cuda.device_count() > 1:
                        trait_devices = [f"cuda:{index}" for index
                                         in range(_torch.cuda.device_count())]
        except Exception as error:  # noqa: BLE001 - never fail a scan over a size choice
            # Falling back to 'no blocking' is a safe default only when the
            # whole trait matrix genuinely fits. When it does not, swallowing
            # this silently puts an 80 GB design on one card and leaves the
            # other seven idle -- which is exactly what it did, invisibly,
            # until a run was watched with nvidia-smi. Say something.
            warnings.warn(
                f'automatic trait blocking failed ({type(error).__name__}: {error}); the scan will use the whole trait axis on one device',
                RuntimeWarning, stacklevel=2)
    if trait_devices and trait_block is None:
        raise ValueError(
            "trait_devices requires trait_block: the devices are given trait "
            "blocks to work on, so there must be blocks to give")

    # PREFLIGHT. Refuse an impossible plan here, before a byte is read.
    #
    # Supplementary Methods S3.1 lists as a limitation that "users must
    # currently split phenotype groups manually if the processed matrix does
    # not fit". Until now the unreduced path had no feasibility check at all:
    # an over-large trait panel ran until CUDA raised out-of-memory, after the
    # genotype read had begun, naming neither the cause nor a setting that
    # works. The calculator can answer this in milliseconds, so it should.
    #
    # Only when nothing has already rescued the plan -- a caller-supplied or
    # auto-derived `trait_block` means the traits are being split, which is
    # the fix this would otherwise recommend.
    if (trait_block is None and str(device).startswith("cuda")
            and isinstance(genotype, ChunkedGenotype)):
        try:
            import torch as _torch

            from .preflight import require_fit

            if _torch.cuda.is_available():
                free_bytes, _total = _torch.cuda.mem_get_info(
                    _torch.device(device))
                width = getattr(genotype, "native_row_width", 0)
                if (getattr(genotype, "native_encoding", None)
                        in {"pgen_2bit", "plink_2bit"} and width):
                    per_variant = float(width)
                elif hasattr(genotype, "iter_packed_chunks"):
                    per_variant = float((int(genotype.shape[0]) + 3) // 4)
                else:
                    per_variant = float(genotype.shape[0]) * 4.0
                require_fit(
                    chunk_variants=int(chunk_size or 4096),
                    depth=int(prefetch_chunks or 4),
                    n_samples=int(genotype.shape[0]),
                    n_traits=int(phenotype.shape[1]),
                    covariate_rank=int(0 if covariates is None
                                       else np.asarray(covariates).shape[1]),
                    transfer_bytes_per_variant=per_variant,
                    device_memory_bytes=float(free_bytes),
                    reduced=(reduction is not None or significance is not None))
        except ImportError:
            pass
        # PlanTooLarge is NOT caught: refusing early with a workable setting is
        # the entire point, and swallowing it would restore the opaque OOM this
        # replaces. Any OTHER failure in the size estimate must not sink a scan
        # that might well have run, so it warns and continues.
        except Exception as error:  # noqa: BLE001
            from .preflight import PlanTooLarge
            if isinstance(error, PlanTooLarge):
                raise
            warnings.warn(
                f"preflight size check failed ({type(error).__name__}: "
                f"{error}); running without it",
                RuntimeWarning, stacklevel=2)
    if isinstance(genotype, ChunkedGenotype):
        fused_bed_qc = (
            resolved_device.type == "cuda"
            and resolved_compute_dtype == "float32"
            and (hasattr(genotype, "iter_packed_chunks") or getattr(genotype, "supports_fused_qc", False))
        )
        effective_chunk_size = (
            chunk_size
            or (
                getattr(genotype, "preferred_gpu_chunk_size", None)
                if fused_bed_qc
                else None
            )
            or getattr(genotype, "preferred_chunk_size", None)
            or min(genotype.shape[1], 4096)
            or 1
        )
        # Work in the precision the scan actually computes in. The float32
        # CUDA path ends at `torch.as_tensor(pheno_proc, torch.float32)`, so
        # promoting to float64 here builds a full-size copy that is downcast
        # again and discarded: 160 GB at 600,000 voxels, 557 GB at 2,085,000.
        # Phase clocks. The 600,000-trait stress run reported 4,183 s total
        # against 2,379 s of `scan_and_write_seconds`, leaving 43% of the wall
        # attributed to nothing at all -- and at 2,085,000 voxels that unnamed
        # 30 minutes becomes an unnamed hour and three quarters. One
        # perf_counter per boundary costs nothing and turns "somewhere in
        # setup" into a number that can be argued with.
        _phase_prep_started = time.perf_counter()
        phenotype, covariates, qc = prepare_inputs_for_prep(
            genotype.genotype,
            phenotype,
            covariates,
            genotype_chunk_size=effective_chunk_size,
            validate_genotype=not fused_bed_qc,
            dtype=(np.float32 if resolved_compute_dtype == 'float32'
                   else np.float64),
        )
        _phase_prep_done = time.perf_counter()
        if pipeline_profile is not None:
            if resolved_device.type != "cuda" or resolved_compute_dtype != "float32":
                raise ValueError("pipeline profiles currently describe the float32 CUDA scan")
            from .pipeline_model import apply_pipeline_profile
            profile = (json.loads(Path(pipeline_profile).read_text())
                       if isinstance(pipeline_profile, (str, Path)) else dict(pipeline_profile))
            profile["candidates"] = dict(profile.get("candidates", {}))
            if chunk_size is not None:
                profile["candidates"]["chunk_variants"] = [chunk_size]
            if requested_reader_workers is not None:
                profile["candidates"]["workers"] = [requested_reader_workers]
            if requested_prefetch_chunks is not None:
                profile["candidates"]["depths"] = [requested_prefetch_chunks]
            if hasattr(genotype, "resolve_decode_backend"):
                genotype.resolve_decode_backend(resolved_device)
            planning = apply_pipeline_profile(genotype, phenotype, covariates, profile)
            chunk_size = planning["scan_kwargs"]["chunk_size"]
            reader_workers = planning["scan_kwargs"]["reader_workers"]
            prefetch_chunks = planning["scan_kwargs"]["prefetch_chunks"]
            genotype_meta["pipeline_plan"] = planning
            qc["genotype_qc_chunk_size"] = chunk_size
        marker_names = (
            [f"marker_{i}" for i in range(genotype.shape[1])]
            if marker_ids is None
            else marker_ids[: genotype.shape[1]]
        )
        trait_names = trait_columns or [f"trait_{i}" for i in range(phenotype.shape[1])]
        genotype_shape = list(genotype.shape)
        # A ranged scan reports on its range, not on the file. The variant count
        # and the marker names both have to be narrowed here, or the sumstats
        # writer is told to expect every variant and the names it writes belong
        # to the wrong rows -- which would be wrong output rather than an error.
        if variant_range is not None:
            first, last = _resolve_variant_range(variant_range, genotype_shape[1])
            genotype_shape[1] = last - first
            if marker_names is not None:
                marker_names = marker_names[first:last]
        if output_dir is not None:
            # Who may borrow the result ring instead of being handed copies.
            # Audited individually, because getting this wrong corrupts output
            # silently rather than raising:
            #   safe   `_drain_linear_chunks`        reads shapes, retains nothing
            # Discard-only scans may borrow result buffers. Binary writers
            # retain owned arrays until their background writes complete.
            borrow_results = (sumstats_format == "none"
                              and trait_block is None)

            dense_binary = (
                sumstats_format == "binary"
                and significance is None
                and jagwas is None
                and reduction is None
                and topk_per_trait is None
                and p_value_threshold is None
            )

            def _scan(trait_slice=None, device=None):
                return linear_scan_streaming_chunks(
                    genotype,
                    phenotype if trait_slice is None else phenotype[:, trait_slice],
                    covariates,
                    chunk_size=chunk_size,
                    device=str(device or resolved_device),
                    compute_dtype=resolved_compute_dtype,
                    reader_workers=reader_workers,
                    prefetch_chunks=prefetch_chunks,
                    # `significance is None` is not a detail: it was the
                    # single largest cost in a voxel-scale run. This flag
                    # makes the scan run `scipy.special.stdtr` over the FULL
                    # chunk x K t-matrix, on one host core, for every chunk --
                    # 153.6 million p-values per chunk at K = 150,000 -- and
                    # the significance path then throws every one of them
                    # away: `_significant_pairs_iterator` takes that chunk as
                    # `_p_chunk` and ignores it, and the writer computes
                    # -log10(P) only for pairs that clear the threshold. It
                    # presented as a multi-GPU problem (one core pinned, all
                    # cards idle, disk idle, worse as K grows, unaffected by
                    # running one process per card) and it was not one.
                    compute_p_values=False,
                    compute_log10_p=dense_binary,
                    variant_range=variant_range,
                    reduction=reduction,
                    borrow_results=borrow_results,
                )

            blocked_significance = False
            if trait_block is None:
                chunk_iterator, q_matrix = _scan()
            else:
                # One block's design is resident at a time, so the device never
                # holds the full (samples x K) matrix.
                #
                # `q_matrix` comes from the *covariates*, not the phenotype, so
                # it is identical for every block and is derived directly rather
                # than by starting a scan and throwing it away -- a discarded
                # iterator would leave its reader threads and pinned ring
                # allocated for the life of the run.
                from .preprocess import _covariate_basis

                q_matrix = (None if covariates is None or covariates.shape[1] == 0
                            else _covariate_basis(covariates))

                def _blocked(offset, width, device=None):
                    iterator, _ = _scan(slice(offset, offset + width), device)
                    return iterator

                blocked_significance = significance is not None
                if blocked_significance:
                    chunk_iterator = _trait_blocked_significant_chunks(
                        _blocked, significance, int(phenotype.shape[1]),
                        int(trait_block),
                        genotype_shape[0] - (0 if q_matrix is None
                                             else q_matrix.shape[1]) - 2,
                        devices=trait_devices)
                else:
                    chunk_iterator = _trait_blocked_reduced_chunks(
                        _blocked, reduction, int(phenotype.shape[1]),
                        int(trait_block), devices=trait_devices)
            if variant_range is not None:
                # The scan yields absolute variant indices; `marker_names` and
                # `n_variants` above were narrowed to the range. Rebase the
                # chunk bounds so every consumer indexes the same way -- without
                # this the writers index sliced names with absolute positions,
                # which the test caught as an IndexError at the range's end and
                # which would otherwise have written the wrong marker per row.
                offset = _resolve_variant_range(variant_range, genotype.shape[1])[0]

                def _rebased(source):
                    # Tuple-length agnostic: a reduced scan carries a sixth
                    # element (the trait index) that must survive rebasing.
                    # A trait-blocked SIGNIFICANCE chunk is different again:
                    # it is seven long and its third element is an array of
                    # absolute variant indices, which needs the same shift as
                    # the bounds. Shifting only the bounds there would name
                    # the wrong marker on every row.
                    for chunk in source:
                        if len(chunk) == 7:
                            yield (chunk[0] - offset, chunk[1] - offset,
                                   chunk[2] - offset, *chunk[3:])
                        else:
                            yield (chunk[0] - offset, chunk[1] - offset,
                                   *chunk[2:])

                chunk_iterator = _rebased(chunk_iterator)
            out = mkdir(output_dir)
            residual_df = (
                genotype_shape[0] - (0 if q_matrix is None else q_matrix.shape[1]) - 2
            )
            phenotype_observed_counts = np.asarray(
                qc.get("phenotype_observed_counts", [genotype_shape[0]] * len(trait_names)),
                dtype=np.int64,
            )
            trait_df = phenotype_observed_counts - (0 if q_matrix is None else q_matrix.shape[1]) - 2
            manifest_df = (
                int(residual_df)
                if np.all(trait_df == residual_df)
                else trait_df.astype(int).tolist()
            )
            if significance is not None and not blocked_significance:
                from .linear import _significant_pairs_iterator

                # Placed here, after the residual df exists and after any
                # variant-range rebasing, so the chunk bounds this filter
                # reports are the file's own numbering.
                chunk_iterator = _significant_pairs_iterator(
                    chunk_iterator, significance, len(trait_names),
                    residual_df)
            write_started = time.perf_counter()
            if sumstats_format == "none":
                n_rows, sumstats_summary = _drain_linear_chunks(chunk_iterator)
            elif (significance is not None or jagwas is not None or reduction is not None
                  or topk_per_trait is not None or p_value_threshold is not None):
                from .sumstats_indexed import write_indexed_sumstats
                kind = ("significant" if significance is not None else "jagwas" if jagwas is not None
                        else "reduced" if reduction is not None else "filtered")
                n_rows, sumstats_summary = write_indexed_sumstats(
                    out / "sumstats", marker_names, trait_names, genotype_shape[0], chunk_iterator,
                    kind=kind, df=residual_df,
                    chi2_df=jagwas.degrees_of_freedom if jagwas is not None else None,
                    topk_per_trait=topk_per_trait,p_value_threshold=p_value_threshold,
                    variant_metadata=variant_metadata,fsync=sumstats_fsync,
                    store_beta=sumstats_fields != "t")
            else:
                n_rows, sumstats_summary = _write_linear_binary_streaming(
                    out / "sumstats",
                    marker_names=marker_names,
                    trait_names=trait_names,
                    n_samples=genotype_shape[0],
                    chunk_iterator=chunk_iterator,
                    n_variants=genotype_shape[1],
                    df=manifest_df,
                    block_bytes=sumstats_block_bytes,
                    queue_depth=sumstats_queue_depth or DEFAULT_SUMSTATS_QUEUE_DEPTH,
                    fsync=sumstats_fsync,
                    write_variant_ids=sumstats_variant_ids,
                    store_beta=sumstats_fields != "t",
                    extra_manifest={"genotype_format": genotype_meta.get("format")},
                    borrow_results=borrow_results,
                )
            sumstats_summary["scan_and_write_seconds"] = time.perf_counter() - write_started
            table: list[dict] = []
            p_value = None
            beta = None
            t_stat = None
        else:
            # Refused rather than ignored. `linear_scan_streaming` has no
            # variant range, and silently scanning the whole file when a subset
            # was asked for is precisely the failure that cost this project a
            # day: it produces a plausible number for the wrong extent and
            # raises nothing.
            if variant_range is not None:
                raise ValueError(
                    "variant_range requires output_dir: the in-memory result "
                    "path scans the whole file and would silently ignore it")
            beta, t_stat, logp, q_matrix = linear_scan_streaming(
                genotype,
                phenotype,
                covariates,
                chunk_size=chunk_size,
                device=str(resolved_device),
                compute_dtype=resolved_compute_dtype,
                reader_workers=reader_workers,
                prefetch_chunks=prefetch_chunks,
                return_log10_p=True,
            )
            table = []
            for marker_index, marker_name in enumerate(marker_names):
                for trait_index, trait_name in enumerate(trait_names):
                    row = {
                        "marker_id": marker_name,
                        "trait": trait_name,
                        "n": int(genotype_shape[0]),
                        "-log10_p": float(logp[marker_index, trait_index]),
                    }
                    if return_beta:
                        row["beta"] = float(beta[marker_index, trait_index])
                    if return_t:
                        row["t_stat"] = float(t_stat[marker_index, trait_index])
                    if return_se:
                        denom = abs(t_stat[marker_index, trait_index])
                        row["se"] = float(abs(beta[marker_index, trait_index]) / denom) if denom > 0 else float("nan")
                    if variant_metadata is not None:
                        for field in ("chromosome", "position", "effect_allele", "other_allele"):
                            value = variant_metadata[field][marker_index]
                            row[field] = int(value) if field == "position" else str(value)
                    table.append(row)
            n_rows = int(beta.shape[0] * beta.shape[1])
        scan_exclusions = getattr(genotype, "_last_scan_exclusion_counts", None)
        if scan_exclusions is not None:
            qc["genotype_exclusion_counts"] = {
                key: int(value) for key, value in scan_exclusions.items()
            }
            qc["n_variants_excluded"] = int(sum(scan_exclusions.values()))
    else:
        genotype, phenotype, covariates, qc = prepare_inputs(genotype, phenotype, covariates)
        beta, t_stat, logp, q_matrix = linear_scan(
            genotype,
            phenotype,
            covariates,
            chunk_size=chunk_size,
            device=str(resolved_device),
            compute_dtype=resolved_compute_dtype,
            return_log10_p=True,
        )
        genotype_shape = list(genotype.shape)
        # Zero-variance columns are dropped by prepare_inputs, so labels follow
        # the retained indices rather than a prefix of the original order.
        kept_markers = qc.get("genotype_kept_column_indices")
        kept_traits = qc.get("phenotype_kept_column_indices")
        if marker_ids is None:
            marker_names = [f"marker_{i}" for i in (kept_markers or range(beta.shape[0]))]
        else:
            selected = kept_markers if kept_markers is not None else range(beta.shape[0])
            marker_names = [str(marker_ids[i]) for i in selected]
        if trait_columns is None:
            trait_names = [f"trait_{i}" for i in (kept_traits or range(beta.shape[1]))]
        elif kept_traits is not None:
            trait_names = [trait_columns[i] for i in kept_traits]
        else:
            trait_names = list(trait_columns)
        table = []
        for marker_index, marker_name in enumerate(marker_names):
            for trait_index, trait_name in enumerate(trait_names):
                row = {
                    "marker_id": marker_name,
                    "trait": trait_name,
                    "n": int(genotype_shape[0]),
                    "-log10_p": float(logp[marker_index, trait_index]),
                }
                if return_beta:
                    row["beta"] = float(beta[marker_index, trait_index])
                if return_t:
                    row["t_stat"] = float(t_stat[marker_index, trait_index])
                if return_se:
                    denom = abs(t_stat[marker_index, trait_index])
                    row["se"] = float(abs(beta[marker_index, trait_index]) / denom) if denom > 0 else float("nan")
                table.append(row)
        n_rows = len(table)

    analysis_chunk_size = int(
        chunk_size
        or (
            getattr(genotype, "preferred_gpu_chunk_size", None)
            if resolved_device.type == "cuda"
            else None
        )
        or getattr(genotype, "preferred_chunk_size", None)
        or min(genotype_shape[1], 4096)
    )
    # The covariate SVD truncates rank-deficient input, and until now it did so
    # silently: the only visible trace was a `df` that disagreed with the number
    # of covariate columns supplied, which nothing in the output explained.
    # Collinear covariates are common enough (a redundant principal component, a
    # dummy set including its own reference level) that a run must say when it
    # happened. plink2 reports this; so do we now.
    covariate_rank_used = 0 if q_matrix is None else int(q_matrix.shape[1])
    covariate_components_dropped = (
        None if covariates is None
        else int(covariates.shape[1]) - covariate_rank_used)
    if covariate_components_dropped:
        warnings.warn(
            f"{covariate_components_dropped} of {covariates.shape[1]} covariate "
            f"columns are linearly dependent and were dropped; the scan used "
            f"rank {covariate_rank_used}, so the residual degrees of freedom are "
            f"{covariate_components_dropped} higher than the column count implies",
            RuntimeWarning,
            stacklevel=2,
        )
    run_metadata = {
        "analysis": "linear",
        "chunk_size": analysis_chunk_size,
        "device_requested": device,
        "device_used": str(resolved_device),
        "gpu_name": (
            torch.cuda.get_device_name(resolved_device)
            if resolved_device.type == "cuda"
            else None
        ),
        "compute_dtype_requested": compute_dtype,
        "compute_dtype_used": resolved_compute_dtype,
        "topk_per_trait": topk_per_trait,
        "p_value_threshold": p_value_threshold,
        "reduce": ("jagwas" if jagwas is not None
                   else "significant" if significance is not None
                   else (None if reduction is None else reduction.mode)),
        "significance_threshold": (
            None if significance is None
            else significance.resolved_threshold(len(trait_names))),
        "reduce_top_k": None if reduction is None else reduction.width,
        # Recorded because it changes how many passes the run made over the
        # genotypes, which is the first thing to check against a wall time.
        "trait_block": None if trait_block is None else int(trait_block),
        # WHICH devices the blocks ran on, not just how wide the blocks were.
        # The 600,000-trait stress run recorded `trait_block: 75000` and
        # nothing else, so "did it use eight cards or one?" could only be
        # argued from the fact that 75,000 is 600,000/8 -- an inference, when
        # the run already knew the answer and simply never wrote it down.
        "trait_devices": (None if not trait_devices
                          else [str(d) for d in trait_devices]),
        "genotype_shape": genotype_shape,
        "phenotype_shape": list(phenotype.shape),
        "covariate_shape": None if covariates is None else list(covariates.shape),
        # The rank actually used, not the number of columns supplied. The SVD
        # truncates rank-deficient covariates, and the retained rank is what
        # sets `df` and so every p-value -- so a run whose covariates were
        # collinear must say so rather than leaving the reader to infer it from
        # a df that disagrees with the column count.
        "covariate_rank": covariate_rank_used,
        "covariate_components_dropped": covariate_components_dropped,
        "has_sample_ids": bool(sample_ids is not None),
        "has_marker_ids": bool(marker_ids is not None),
        "sample_id_column": sample_id_column,
        "trait_columns": trait_names,
        "covariate_columns": covariate_columns,
        "q_matrix_shape": None if q_matrix is None else list(q_matrix.shape),
        "results_streamed": bool(isinstance(genotype, ChunkedGenotype) and output_dir is not None),
        "reader_workers": int(reader_workers) if isinstance(genotype, ChunkedGenotype) else None,
        "zstd_read_workers": (
            int(zstd_read_workers)
            if isinstance(genotype, ChunkedGenotype) and zstd_read_workers is not None
            else None
        ),
        "prefetch_chunks": int(prefetch_chunks) if isinstance(genotype, ChunkedGenotype) else None,
        "n_result_rows": int(n_rows),
        "sumstats_format": sumstats_format,
        "sumstats_fields": sumstats_fields,
        "sumstats_write": sumstats_summary or None,
        "phase_seconds": _phase_breakdown(
            _phase_entered, _phase_prep_started, _phase_prep_done,
            write_started, time.perf_counter()),
        "runtime_seconds": elapsed(start),
        "version": "0.1.0",
    }
    if hasattr(genotype, 'backend_used'):
        genotype_meta['decode_backend_used'] = genotype.backend_used
        if hasattr(genotype, 'backend_reason'):
            genotype_meta['decode_backend_reason'] = genotype.backend_reason
    run_metadata.update(genotype_meta)
    result = GWASResult(table=table, run_metadata=run_metadata, qc_summary=qc)
    if output_dir is not None:
        out = mkdir(output_dir)
        streamed = isinstance(genotype, ChunkedGenotype) and run_metadata["results_streamed"]
        if not streamed:
            if sumstats_format == "binary" and (topk_per_trait is not None or p_value_threshold is not None):
                from .sumstats_indexed import write_indexed_sumstats
                written, sumstats_summary = write_indexed_sumstats(
                    out / "sumstats", marker_names, trait_names, genotype_shape[0],
                    iter([(0, len(marker_names), beta, t_stat, None)]), kind="filtered",
                    df=genotype_shape[0] - covariate_rank_used - 2,
                    topk_per_trait=topk_per_trait,p_value_threshold=p_value_threshold,
                    variant_metadata=variant_metadata,fsync=sumstats_fsync,
                    store_beta=sumstats_fields != "t")
                run_metadata["n_result_rows"] = written
                run_metadata["sumstats_write"] = sumstats_summary
            elif sumstats_format == "binary":
                n_variants = int(np.shape(beta)[0])
                _, sumstats_summary = _write_linear_binary_streaming(
                    out / "sumstats",
                    marker_names=marker_names,
                    trait_names=trait_names,
                    n_samples=genotype_shape[0],
                    chunk_iterator=iter([(0, n_variants, beta, t_stat, None, logp)]),
                    n_variants=n_variants,
                    df=genotype_shape[0]
                    - (0 if q_matrix is None else q_matrix.shape[1])
                    - 2,
                    block_bytes=sumstats_block_bytes,
                    queue_depth=sumstats_queue_depth or DEFAULT_SUMSTATS_QUEUE_DEPTH,
                    fsync=sumstats_fsync,
                    write_variant_ids=sumstats_variant_ids,
                    store_beta=sumstats_fields != "t",
                )
                run_metadata["sumstats_write"] = sumstats_summary
        write_json(run_metadata, out / "run.json")
        write_json(qc, out / "qc.json")
    return result
