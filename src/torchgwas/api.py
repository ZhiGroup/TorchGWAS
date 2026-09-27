from __future__ import annotations

import json
import time
import warnings
from pathlib import Path


import numpy as np
import torch

from .io import align_table_to_samples, load_array, load_genotype, load_vector
from .executor_timing import streaming_timing
from .linear import linear_scan, linear_scan_streaming, linear_scan_streaming_chunks
from .pgen import (
    DEFAULT_PGEN_COMPRESSION_WORKERS,
    DEFAULT_PGEN_DECODE_BATCH_SIZE,
    DEFAULT_PGEN_DECODE_WORKERS,
)
from .preprocess import prepare_inputs, prepare_inputs_for_prep
from .streaming import ChunkedGenotype, _resolve_variant_range
from .types import GWASResult
from .utils import (choose_device, elapsed, mkdir,
                    timestamp, upper_tail_log10, write_json)


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


def _coerce_array_or_path(value, *, mmap=False):
    if value is None:
        return None
    if isinstance(value, (str, Path)):
        if mmap and Path(value).suffix.lower()=='.npy':
            return np.load(value,mmap_mode='r',allow_pickle=False)
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
                                      trait_block, df, devices=None, queue_depth=None,
                                      partition_context=None):
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
    if partition_context is not None and not callable(partition_context):
        raise ValueError('Callable indexed producer binding required')
    blocks = [(offset, min(trait_block, n_traits - offset))
              for offset in range(0, n_traits, trait_block)]

    def emit(offset, width, device):
        partition=None if partition_context is None else partition_context(offset,width,device)
        source = scan_once(offset, width, device)
        pairs = _significant_pairs_iterator(source, significance, n_traits, df)
        try:
            for start, end, variant_index, trait_index, beta, t_stat, row_df in pairs:
                # Keep an owned copy: selected inputs may be views or reused by
                # their producer. Rebase it in place to avoid a second array.
                rebased_trait_index=trait_index.astype(np.int64)
                rebased_trait_index+=offset
                item=(start,end,variant_index,rebased_trait_index,beta,t_stat,row_df)
                if partition is not None:
                    from .sumstats_indexed import PartitionedIndexedChunk
                    item=PartitionedIndexedChunk(item,partition)
                yield item
        finally:
            try:
                if hasattr(pairs, "close"):
                    pairs.close()
            finally:
                if hasattr(source, "close"):
                    source.close()
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

    if queue_depth is None:
        queue_depth = 4 * len(devices)
    if isinstance(queue_depth, bool) or not isinstance(queue_depth, int) or queue_depth < 1:
        raise ValueError('queue_depth must be a positive integer')
    results: queue.Queue = queue.Queue(maxsize=queue_depth)
    failures: list[BaseException] = []
    done = object()
    stop = threading.Event()

    def publish(item):
        while not stop.is_set():
            try:
                results.put(item, timeout=0.05)
                return True
            except queue.Full:
                pass
        return False

    def run(device, assigned):
        iterator = None
        try:
            for offset, width in assigned:
                if stop.is_set():
                    break
                iterator = emit(offset, width, device)
                for item in iterator:
                    if not publish(item):
                        return
                iterator.close()
                iterator = None
        except BaseException as exc:
            failures.append(exc)
            stop.set()
        finally:
            try:
                if iterator is not None:
                    iterator.close()
            except BaseException as exc:
                failures.append(exc)
                stop.set()
            publish(done)

    workers = []
    for index, device in enumerate(devices):
        assigned = blocks[index::len(devices)]
        if not assigned:
            continue
        thread = threading.Thread(target=run, args=(device, assigned),
                                  name=f"torchgwas-sigshard-{index}", daemon=True)
        thread.start()
        workers.append(thread)
    remaining = len(workers)
    try:
        while remaining:
            if failures:
                raise failures[0]
            try:
                item = results.get(timeout=0.05)
            except queue.Empty:
                continue
            if item is done:
                remaining -= 1
            else:
                yield item
    finally:
        stop.set()
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
    store_variant_df: bool = False,
    on_write_progress=None,
) -> tuple[int, dict]:
    """Stream beta/t_stat chunks into a binary sumstats directory.

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
        store_variant_df=store_variant_df,
        on_write_progress=on_write_progress,
    )
    try:
        for chunk in chunk_iterator:
            # (start, end, beta, t, p, -log10 P[, df]): the scan computed
            # -log10 P on the device (compute_log10_p).
            start, end, beta_chunk, t_chunk, _p_chunk, logp_chunk = chunk[:6]
            if store_variant_df and len(chunk)!=7:
                raise ValueError('Per-variant df output requires scan df metadata')
            writer.write_chunk(start, end, beta_chunk, t_chunk, logp_chunk,
                               variant_df=chunk[6] if store_variant_df else None)
    except BaseException:
        writer.abort()
        raise
    finally:
        if hasattr(chunk_iterator,'close'):
            chunk_iterator.close()
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


from .initial_chunk_autotune import productive_api_lifecycle
from .host_pages import numpy_hugepage_advice, without_numpy_hugepages


@without_numpy_hugepages
@productive_api_lifecycle
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
    sumstats_variant_ids: bool = False,
    sumstats_fields: str = "beta+t",
    autotune_profile: dict | str | Path | None = None,
    autotune_config: dict | str | Path | None = None,
    variant_devices=None,
    initial_calibration: dict | None = None,
    autotune: bool | str | None = None,
    autotune_options: dict | None = None,
    # reduce='jagwas' only, one or neither. rcond (default 1e-3,
    # TORCHGWAS_JAGWAS_RCOND): eigen truncation, keeping R's eigen-directions
    # above rcond x the largest eigenvalue; 0 selects the rounding cutoff over
    # traits. min_residual: drop traits while less than that fraction of a
    # trait's variance is its own given the traits kept before it (VIF >
    # 1 / min_residual). The planner prices TORCHGWAS_JAGWAS_RCOND's method.
    jagwas_rcond: float | None = None,
    jagwas_min_residual: float | None = None,
    # reduce='jagwas' only: [(name, phenotype column indices[, cutoff]), ...] or
    # a dict, for one joint test per group from a single genotype pass (one chi2
    # column and one df per group; jagwas_projection.JagwasGroups). A cutoff is
    # a dict with rcond or min_residual, or a bare rcond; jagwas_rcond and
    # jagwas_min_residual apply to groups without their own.
    jagwas_groups=None,
    # Opt-in: mask phenotype values beyond this many SD of their
    # covariate-residualised trait before the scan (preprocess.
    # mask_phenotype_outliers). reduce='jagwas' drops the sample's whole
    # panel row; per-trait scans drop the value.
    phenotype_outlier_sd: float | None = None,
) -> GWASResult:
    """Run associations, optionally saving bounded early productive measurements.

    initial_calibration={'cache_dir': ...} enables dependency-bound observation
    reuse for explicit native hardcall PGEN/CUDA FP32 configurations. Supply
    chunk_size; this option neither benchmarks nor searches before useful work,
    and does not automatically change the configuration. See the calibration
    guide for window budgets and the distinction between spans and capacities.

    autotune=True (or 'empirical') picks GPUs, phenotype tiles or variant
    shards and reader workers at startup, then tunes the chunk size from the
    job's own first chunks (see empirical_autotune). Explicit trait_block,
    trait_devices, variant_devices, reader_workers or a cuda:N device are kept.
    autotune_options may set chunk_sizes, devices, min_tile_traits,
    max_tile_traits, cpus_per_device (fixed CPU cap instead of the measured
    decode demand), shard_setup_seconds (instead of the measured per-GPU
    setup; 0 ignores setup), min_job_seconds, shared_decode, gpu_fanout,
    warmup_fraction, trial_fraction, repeats, min_gain, max_utilization,
    min_free_bytes, tuner, probe_chunks and split.
    """
    _api_entered=time.perf_counter()
    if jagwas_groups is not None and reduce != "jagwas":
        raise TypeError("jagwas_groups needs reduce='jagwas'")
    # JAGWAS has one joint statistic over the complete retained phenotype
    # panel. Only variants may be partitioned, including across GPUs.
    if reduce == "jagwas" and (trait_block is not None or trait_devices is not None):
        raise ValueError(
            "jagwas cannot be trait-blocked: trait_block and trait_devices "
            "are unsupported. The full phenotype panel and joint-test state "
            "must fit on every active device. Use variant_devices to shard "
            "variants; reduce='significant' supports phenotype partitioning.")
    if variant_devices is not None:
        if (not isinstance(variant_devices,(list,tuple)) or not variant_devices
            or len(set(map(str,variant_devices)))!=len(variant_devices)):
            raise ValueError('variant_devices must be a nonempty unique device list')
        variant_devices=list(map(str,variant_devices))
        if (trait_block is not None or trait_devices is not None or reduce not in (None, "jagwas", "significant")
            or _internal_reduction is not None or pipeline_profile is not None
            or output_dir is None or sumstats_format!='binary'
            or p_value_threshold is not None):
            raise ValueError('variant_devices requires full, significant or jagwas binary output without trait tiling, row filters or coarse profiles')
    empirical = None
    empirical_tuner = None
    if autotune not in (None, False):
        if autotune not in (True, 'empirical'):
            raise ValueError("autotune must be True or 'empirical'")
        if (autotune_profile is not None or autotune_config is not None
                or initial_calibration is not None or pipeline_profile is not None):
            raise ValueError("autotune='empirical' cannot be combined with autotune_profile/"
                             "autotune_config, initial_calibration or pipeline_profile")
        if output_dir is None:
            raise ValueError("autotune='empirical' requires output_dir: it tunes the streaming writer path")
        options = dict(autotune_options or {})
        unknown = set(options)-{'chunk_sizes', 'devices', 'min_tile_traits', 'max_tile_traits', 'cpus_per_device',
                                'min_job_seconds', 'shared_decode', 'gpu_fanout',
                                'warmup_fraction', 'trial_fraction', 'repeats', 'min_gain',
                                'max_utilization', 'min_free_bytes', 'tuner', 'probe_chunks', 'split',
                                'shard_setup_seconds'}
        if unknown:
            raise ValueError(f'Unknown autotune_options: {sorted(unknown)}')
        empirical = dict(options=options)
    elif autotune_options is not None:
        raise ValueError('autotune_options requires autotune=True')
    autotuner = None
    productive = None
    if (autotune_profile is None) != (autotune_config is None):
        raise ValueError('Supply autotune_profile and autotune_config together')
    if autotune_profile is not None:
        from .detailed_autotune import DetailedAutotune, validate_autotune_request
        validate_autotune_request(genotype=genotype, phenotype=phenotype,
                                  output_dir=output_dir, options=locals())
    if sumstats_format not in {"binary", "none"}:
        raise ValueError("sumstats_format must be 'binary' or 'none'; TSV output has been removed")
    calibrator = None
    if initial_calibration is not None:
        from .run_calibration import RunCalibration,validate_initial_calibration
        from .pgen import resolve_pgen_triplet
        validate_initial_calibration(initial_calibration,genotype=genotype,options=locals())
        input_pgen,input_pvar,input_psam=resolve_pgen_triplet(genotype,pvar=pvar,
            psam=psam if psam is not None else sample_file)
        calibrator=RunCalibration(initial_calibration,output_path=output_dir,
            started_perf_counter=_api_entered,
            inputs=dict(genotype=input_pgen,pvar=input_pvar,psam=input_psam,
                phenotype=phenotype,covariates=covariates,sample_ids=sample_ids,marker_ids=marker_ids))
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
    autotune_basis = None
    if autotune_profile is not None:
        autotuner = DetailedAutotune(autotune_profile, autotune_config,
                                    input_path=genotype, output_path=output_dir,
                                    reduction=reduce, significance_threshold=significance_threshold)
        # A bounded QC width is needed before the execution tile is selected.
        if reduce != 'jagwas':
            trait_block = autotuner.qc_trait_block
    resolved_device = choose_device(autotuner.devices[0] if autotuner else variant_devices[0] if variant_devices else device)
    # Dense binary output stores -log10 P computed on each scan device
    # (tails.neg_log10_p_device); compiling it for a GPU takes seconds, so it
    # starts here and overlaps input loading and preparation.
    if (output_dir is not None and sumstats_format == "binary" and reduce is None
            and significance_threshold is None and p_value_threshold is None
            and _internal_reduction is None):
        from .tails import prepare_device_tail_async
        prepare_device_tail_async(autotuner.devices if autotuner else
                                  variant_devices or trait_devices or [resolved_device])
    genotype_meta = {}
    requested_sample_ids = _coerce_vector_or_path(sample_ids)
    # Streaming stores record this input instead of rewriting its variant IDs
    # (variant_source.py); these are the options needed to reopen it.
    genotype_input = genotype if isinstance(genotype, (str, Path)) else None
    genotype_load = dict(pvar=pvar, psam=psam, bim=bim, fam=fam, sample_file=sample_file)
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
        phenotype = _coerce_array_or_path(phenotype, mmap=reduce in ("significant", "jagwas") or ((trait_block is not None or variant_devices is not None) and reduce is None))
    if covariates_table is not None:
        if sample_ids is None:
            raise ValueError("tabular covariate input requires genotype sample IDs")
        covariates, covariate_columns = align_table_to_samples(
            covariates_table, sample_ids=np.asarray(sample_ids, dtype=object), value_columns=covariate_columns, sample_id_column=sample_id_column
        )
    else:
        covariates = _coerce_array_or_path(covariates)
    outlier_rows = None
    if phenotype_outlier_sd is not None:
        from .preprocess import mask_phenotype_outliers
        phenotype, outlier_rows = mask_phenotype_outliers(
            phenotype, None if covariates is None else np.asarray(covariates), float(phenotype_outlier_sd),
            whole_rows=reduce == "jagwas")
    if empirical is not None and not isinstance(genotype, ChunkedGenotype):
        # An in-memory genotype array runs one in-memory scan: nothing to tune.
        empirical.update(sizes=None, layout=dict(why=['in-memory genotype array: no streaming scan to tune']),
                         devices_seen=None, tuner_disabled_reason='in-memory genotype array')
        empirical['in_memory'] = True
    if empirical is not None and not empirical.get('in_memory'):
        from .empirical_autotune import (DEFAULT_CHUNK_SIZES, eligible_gpus, gpu_free_bytes,
                                         plan_layout, usable_cpus)
        import threading as _threading
        import torch as _torch
        # The idle-CPU sample (0.25 s) overlaps the nvidia-smi device query.
        _cpu_sample = {}
        _cpu_thread = _threading.Thread(target=lambda: _cpu_sample.setdefault('value', usable_cpus()), daemon=True)
        _cpu_thread.start()
        options = empirical['options']
        first, last = _resolve_variant_range(variant_range, int(genotype.shape[1]))
        sizes = sorted(set(int(s) for s in (options.get('chunk_sizes')
                       or ([chunk_size] if chunk_size is not None else DEFAULT_CHUNK_SIZES))))
        # A framed source decodes whole frames: default candidates become
        # multiples of its frame when chunks start on a frame boundary.
        frame = getattr(genotype, 'chunk_alignment_variants', None)
        if frame and not options.get('chunk_sizes') and chunk_size is None and first % int(frame) == 0:
            from .empirical_autotune import frame_aligned_sizes
            sizes = frame_aligned_sizes(sizes, frame)
        # Keep sizes a job can use; the smallest stays so alignment holds.
        sizes = [s for s in sizes if s <= max(last-first, 1)] or sizes[:1]
        explicit_layout = trait_block is not None or trait_devices is not None or variant_devices is not None
        seen = None
        if options.get('devices') is not None:
            devices = [str(d) for d in options['devices']]
        elif explicit_layout:
            devices = [str(d) for d in (variant_devices or trait_devices or [choose_device(device)])]
        elif str(device).startswith('cuda:'):
            devices = [str(device)]
        else:
            devices, seen = eligible_gpus(max_utilization=options.get('max_utilization', 20),
                                          min_free_bytes=options.get('min_free_bytes', 4 << 30))
            if not devices:
                raise ValueError('autotune found no idle GPU (no foreign process, low utilization, free '
                                 'memory); pass autotune_options={"devices": [...]} to choose explicitly')
        layout = dict(explicit=True, devices=devices)
        if not explicit_layout:
            width = getattr(genotype, "native_row_width", 0)
            if getattr(genotype, "native_encoding", None) in {"pgen_2bit", "plink_2bit"} and width:
                per_variant = float(width)
            elif hasattr(genotype, "iter_packed_chunks"):
                per_variant = float((int(genotype.shape[0]) + 3) // 4)
            else:
                per_variant = float(genotype.shape[0]) * 4.0
            mode = 'jagwas' if reduce == 'jagwas' else 'significant' if reduce == 'significant' else 'full'
            _cpu_thread.join()
            cpus, cpu_detail = _cpu_sample['value']
            unfiltered = (sumstats_format == 'binary'
                          and p_value_threshold is None and _internal_reduction is None)
            # Measured per-variant GPU time and per-shard setup size the
            # shard count; with decode CPU per variant they also cap GPUs on
            # a busy host (cpus_per_device overrides that cap).
            # The same decode/GPU ratio sizes the readers per GPU, so it is
            # measured on one GPU too (~0.1 s); the setup probe needs a second.
            from .empirical_autotune import (decode_cpu_seconds_per_variant, gpu_seconds_per_variant,
                                             shard_setup_seconds)
            decode_cpu = setup = None
            # A grouped JAGWAS job projects and factors each group alone.
            group_sizes = None
            if mode == 'jagwas' and jagwas_groups is not None:
                from .jagwas_projection import JagwasGroups
                group_sizes = [len(columns) for columns in JagwasGroups(jagwas_groups).columns]
            gpu_per_variant = gpu_seconds_per_variant(devices[0], mode=mode, n_samples=int(genotype.shape[0]),
                                                      n_traits=int(np.shape(phenotype)[1]),
                                                      group_sizes=group_sizes)
            if len(devices) > 1:
                setup = (float(options['shard_setup_seconds']) if options.get('shard_setup_seconds') is not None
                         else shard_setup_seconds(devices[1]))  # a GPU the job uses (at least two are kept)
            if options.get('cpus_per_device') is None:
                decode_cpu = decode_cpu_seconds_per_variant(genotype, first, last)
            output_rates = None
            if mode == 'full' and unfiltered and output_dir is not None and len(devices) > 1:
                # Dense output: each variant shard brings its own writer (output_write_rates).
                if options.get('output_rates') is not None:
                    output_rates = dict(options['output_rates'])
                else:
                    # The probe writes at most 5% of the job's own output; below
                    # 16 MB per writer it cannot resolve a rate and is skipped.
                    width = int(np.shape(phenotype)[1])
                    # beta (optional), t and -log10 P per cell, one df per variant.
                    job_bytes = (last - first) * (width * (12 if sumstats_fields != 't' else 8) + 4)
                    probe = min(256 << 20, int(0.05 * job_bytes / (1 + len(devices))))
                    if probe >= 16 << 20:
                        from .empirical_autotune import output_write_rates
                        output_rates = output_write_rates(mkdir(output_dir), n_traits=width,
                                                          store_beta=sumstats_fields != 't', writers=len(devices),
                                                          probe_bytes=probe)
            layout = plan_layout(
                mode=mode, n_samples=int(genotype.shape[0]), n_traits=int(np.shape(phenotype)[1]),
                covariate_rank=int(0 if covariates is None else np.shape(covariates)[1]),
                n_variants=last-first, devices=devices, cpus=cpus, capacity=max(sizes), chunk_sizes=sizes,
                depth=int(prefetch_chunks), transfer_bytes_per_variant=per_variant,
                device_free_bytes=min(gpu_free_bytes(devices, seen).values()),
                host_free_bytes=_available_host_bytes(),
                min_tile_traits=int(options.get('min_tile_traits', 2048)),
                max_tile_traits=int(options.get('max_tile_traits', 32768)),
                cpus_per_device=(int(options['cpus_per_device']) if options.get('cpus_per_device') is not None
                                 else None),
                decode_cpu_per_variant=decode_cpu, gpu_seconds_per_variant=gpu_per_variant,
                shard_setup_seconds=setup, cpu_cores=cpu_detail.get('affinity'),
                cpu_load=cpu_detail.get('load_1min'), allow_partitions=unfiltered, output_rates=output_rates,
                group_sizes=group_sizes)
            layout['cpus'] = cpu_detail
            sizes = layout.get('chunk_sizes') or sizes  # sizes whose rings do not fit are dropped
            layout['per_variant_transfer_bytes'] = per_variant
            split = options.get('split', 'auto')
            if split not in ('auto', 'variants', 'traits'):
                raise ValueError("autotune_options['split'] must be 'auto', 'variants' or 'traits'")
            tiled = layout.get('trait_devices') or []
            if (mode == 'significant' and len(tiled) > 1 and split != 'traits'
                    and layout.get('fit_traits', 0) >= int(np.shape(phenotype)[1])):
                # The panel fits every GPU, so either split works; price data
                # movement with measured bandwidth (layout_pricing).
                from .empirical_autotune import best_shard_count
                from .layout_pricing import choose_split, measure_transfer, peak_fp32_flops, split_costs
                # Tiles keep the tile rule's GPU count; shards take the count
                # that balances per-GPU setup against the work, over every GPU
                # the CPU supply allows (full-scale H100: 2 shards 347 s,
                # 4 shards 106 s). Each layout is priced at its own count.
                considered = list(layout.get('devices_considered') or tiled)
                shard_count = len(tiled)
                if gpu_per_variant and setup:
                    shard_count = max(2, best_shard_count((last - first) * gpu_per_variant, setup, len(considered)))
                shard_devices = considered[:shard_count]
                prices = measure_transfer(list(dict.fromkeys(shard_devices + list(tiled))))
                common = dict(n_samples=int(genotype.shape[0]), n_traits=int(np.shape(phenotype)[1]),
                              n_variants=last-first, genotype_bytes_per_variant=per_variant, **prices)
                costs = dict(
                    variant_shards=split_costs(devices=shard_devices, flops={d: peak_fp32_flops(d) for d in shard_devices},
                                               **common)['variant_shards'],
                    phenotype_tiles=split_costs(devices=list(tiled), flops={d: peak_fp32_flops(d) for d in tiled},
                                                **common)['phenotype_tiles'])
                chosen = choose_split(costs) if split == 'auto' else split
                layout['split_pricing'] = dict(prices=prices, costs=costs, choice=chosen, requested=split,
                                               shard_devices=shard_devices, tile_devices=list(tiled))
                if chosen == 'variants':
                    # reader_workers is a total that shards divide; keep each GPU's readers.
                    per_device = max(1, int(layout['reader_workers']) // max(1, len(tiled)))
                    layout.update(variant_devices=shard_devices, trait_devices=None, trait_block=None,
                                  reader_workers=per_device * len(shard_devices))
                    layout['why'].append(
                        f"variant shards over {len(shard_devices)} GPUs: priced {costs['variant_shards']['total']:.2f} s "
                        f"vs {costs['phenotype_tiles']['total']:.2f} s for {len(tiled)} phenotype tiles")
                else:
                    layout['why'].append(
                        f"phenotype tiles kept: priced {costs['phenotype_tiles']['total']:.2f} s "
                        f"vs {costs['variant_shards']['total']:.2f} s for variant shards")
            trait_block, trait_devices, variant_devices = (layout['trait_block'], layout['trait_devices'],
                                                           layout['variant_devices'])
            device = layout['device']
            if requested_reader_workers is None:
                reader_workers = int(layout['reader_workers'])
                if hasattr(genotype, 'decode_workers'):
                    genotype.decode_workers = reader_workers
            if requested_prefetch_chunks is None:
                prefetch_chunks = int(layout['prefetch_chunks'])
        resolved_device = choose_device(variant_devices[0] if variant_devices else
                                        trait_devices[0] if trait_devices else device)
        chunk_size = max(sizes)
        empirical.update(sizes=sizes, layout=layout, devices_seen=seen)
    resolved_compute_dtype = _resolve_linear_compute_dtype(genotype, compute_dtype)
    if calibrator is not None and resolved_compute_dtype!='float32':
        raise ValueError('initial_calibration requires resolved float32 computation')
    if autotuner is not None:
        autotuner.validate_inputs(genotype, phenotype, covariates)
    if sumstats_fields not in {"beta+t", "t"}:
        raise ValueError(
            f"sumstats_fields must be 'beta+t' or 't', got {sumstats_fields!r}"
        )
    if sumstats_fields == "t" and sumstats_format != "binary":
        raise ValueError(
            "sumstats_fields='t' applies only to sumstats_format='binary'"
        )
    if sumstats_format == "none" and (p_value_threshold is not None):
        raise ValueError("row selection requires binary output")
    if p_value_threshold is not None and not (0.0 < p_value_threshold <= 1.0):
        raise ValueError("p_value_threshold must be in (0, 1]")
    # A reduction changes what the scan produces, not merely which rows are
    # kept, so every combination that assumes a full (variant x trait) matrix is
    # refused here rather than silently producing a narrower one.
    if (jagwas_rcond is not None or jagwas_min_residual is not None) and reduce != "jagwas":
        raise TypeError("jagwas_rcond and jagwas_min_residual need reduce='jagwas'")
    reduction = None
    significance = None
    jagwas = None
    if reduce == "jagwas":
        if output_dir is None:
            raise ValueError(
                "reduce='jagwas' requires output_dir: it is a streaming "
                "reduction and the in-memory path returns the full matrix")
        if reduce_top_k is not None:
            raise ValueError("reduce_top_k does not apply to 'jagwas'")
        from .jagwas_projection import JagwasGroups, JagwasReduction
        jagwas = (JagwasReduction(rcond=jagwas_rcond, min_residual=jagwas_min_residual)
                  if jagwas_groups is None else
                  JagwasGroups(jagwas_groups, rcond=jagwas_rcond, min_residual=jagwas_min_residual))
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
        # Use the native pipeline. The optional device selector transfers only
        # passing pairs; the host fallback filters complete result chunks.
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
    full_tiled_output = trait_block is not None and reduction is None and significance is None
    full_variant_output = variant_devices is not None and jagwas is None and significance is None
    jagwas_variant_output = variant_devices is not None and jagwas is not None
    # Significant pairs are selected per (variant, phenotype) cell: variant
    # shards need no cross-shard state, only the shared indexed writer.
    significant_variant_output = variant_devices is not None and significance is not None
    if variant_devices is not None and not isinstance(genotype,ChunkedGenotype):
        raise ValueError('variant_devices requires a streaming genotype source')
    if full_tiled_output and pipeline_profile is not None:
        raise ValueError('The coarse pipeline profile does not model full-output trait tiling')
    if trait_block is not None:
        if full_tiled_output and (output_dir is None or sumstats_format!='binary'
                                  or p_value_threshold is not None
                                  or not isinstance(genotype,ChunkedGenotype)):
            raise ValueError(
                'full-output trait_block requires streaming genotype, output_dir and unfiltered binary output')
        if isinstance(trait_block,bool) or not isinstance(trait_block,(int,np.integer)):
            raise ValueError('trait_block must be a positive integer')
        if int(trait_block) < 1:
            raise ValueError("trait_block must be positive")
    elif ((reduction is not None or significance is not None)
          and jagwas is None and str(device).startswith("cuda") and empirical is None):
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
                    reduced=(reduction is not None or significance is not None),
                    # Dense binary output stages -log10 P as well.
                    compute_log10_p=(output_dir is not None and sumstats_format == "binary"
                                     and reduction is None and significance is None
                                     and p_value_threshold is None))
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
        import os as _os_env
        phenotype, covariates, qc = prepare_inputs_for_prep(
            genotype.genotype,
            phenotype,
            covariates,
            genotype_chunk_size=effective_chunk_size,
            validate_genotype=not fused_bed_qc,
            dtype=(np.float32 if resolved_compute_dtype == 'float32'
                   else np.float64),
            phenotype_block_size=(autotuner.qc_trait_block if autotuner is not None else int(trait_block) if full_tiled_output or (significance is not None and trait_block is not None) else min(4096,phenotype.shape[1]) if full_variant_output or jagwas is not None or significance is not None else None),
            qc_device=(str(resolved_device) if getattr(resolved_device, 'type', None) == 'cuda'
                       and _os_env.environ.get('TORCHGWAS_PHENOTYPE_QC', 'device') == 'device' else None),
        )
        _phase_prep_done = time.perf_counter()
        if jagwas is not None:
            # Missing values are mean-imputed like full output's; JAGWAS then
            # takes the imputed panel's t with the common df (linear.py).
            if phenotype.shape[1] > phenotype.shape[0]-1:
                raise ValueError('jagwas trait count exceeds the residual phenotype rank')
            # Blocked QC preserves source precision and lazy column subsets.
            # The joint preprocessing and factor must use the scan precision.
            phenotype = np.asarray(phenotype, dtype=np.float32 if resolved_compute_dtype == 'float32' else np.float64)
        if significance is not None and not qc['phenotype_missing_cells']:
            significance.prepare_integer_df(int(phenotype.shape[0]), int(phenotype.shape[1]))
        if autotuner is not None:
            settings, audit, autotune_basis = autotuner.select(genotype, phenotype, covariates, qc,
                output=dict(block_bytes=sumstats_block_bytes,
                    queue_depth=DEFAULT_SUMSTATS_QUEUE_DEPTH if sumstats_queue_depth is None else sumstats_queue_depth,
                    store_beta=jagwas is None and sumstats_fields!='t', fsync=sumstats_fsync))
            chunk_size = settings['chunk_size']
            trait_block = settings.get('trait_block')
            trait_devices = settings.get('trait_devices')
            variant_devices = settings.get('variant_devices')
            jagwas_variant_output = variant_devices is not None and jagwas is not None
            full_variant_output = variant_devices is not None and jagwas is None and significance is None
            significant_variant_output = variant_devices is not None and significance is not None
            full_tiled_output = trait_block is not None and reduction is None and significance is None
            effective_chunk_size = chunk_size
            reader_workers = settings['reader_workers']
            # The path-backed source was opened before the planner chose its
            # reader budget. Apply that choice to its native decode preference
            # as well, including choices above the loader's initial default.
            if hasattr(genotype, 'decode_workers'):
                genotype.decode_workers = reader_workers
            prefetch_chunks = settings['prefetch_chunks']
            resolved_device = choose_device(settings['device'])
            qc['genotype_qc_chunk_size'] = chunk_size
            qc['phenotype_qc_trait_block'] = autotuner.qc_trait_block
            genotype_meta['autotune'] = audit
            productive = getattr(autotuner,'productive',None)
        if jagwas is not None:
            # Check the actual selected devices, after live planner admission.
            # Every active worker retains the complete phenotype factor.
            from .reduction_tensor_work import require_jagwas_factor_capacity
            factor_devices = [str(resolved_device)]
            if variant_devices is not None:
                from .linear import multigpu_variant_ranges
                first, last = _resolve_variant_range(variant_range, int(genotype.shape[1]))
                active = len(multigpu_variant_ranges(last-first, int(effective_chunk_size), len(variant_devices))) if last > first else 1
                factor_devices = variant_devices[:max(1,active)]
            if jagwas_groups is None:
                require_jagwas_factor_capacity(int(phenotype.shape[0]), int(phenotype.shape[1]),
                    factor_devices, compute_dtype=resolved_compute_dtype,
                    method='eigen' if jagwas.rcond is not None else 'rounding')
            else:
                sizes = [len(columns) for columns in jagwas.columns]
                largest = jagwas.reductions[int(np.argmax(sizes))]
                require_jagwas_factor_capacity(int(phenotype.shape[0]), int(phenotype.shape[1]),
                    factor_devices, compute_dtype=resolved_compute_dtype,
                    method='eigen' if largest.rcond is not None else 'rounding', group_sizes=sizes)
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
        original_trait_names = trait_columns or [f"trait_{i}" for i in range(qc['phenotype_columns_input'])]
        trait_names = ([original_trait_names[i] for i in qc['phenotype_kept_column_indices']]
                       if 'phenotype_kept_column_indices' in qc else original_trait_names)
        if jagwas_groups is not None and 'phenotype_kept_column_indices' in qc:
            # Group columns index the input panel. QC can drop traits (constant
            # or empty columns), so map each group onto the scanned panel;
            # unmapped indices would silently join the wrong traits.
            position = {int(original): i for i, original in enumerate(qc['phenotype_kept_column_indices'])}
            remapped = []
            for name, columns, cutoff in zip(jagwas.names, jagwas.columns, jagwas.cutoffs):
                kept = [position[int(column)] for column in columns if int(column) in position]
                if len(kept) < len(columns):
                    warnings.warn(f"jagwas group {name}: {len(columns) - len(kept)} of {len(columns)} "
                                  f"traits removed by phenotype QC", UserWarning, stacklevel=2)
                if not kept:
                    raise ValueError(f"jagwas group {name} has no trait left after phenotype QC")
                remapped.append((name, kept, cutoff))
            jagwas = reduction = JagwasGroups(remapped)
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
            if variant_metadata is not None:
                variant_metadata = {key: np.asarray(value)[first:last] for key, value in variant_metadata.items()}
        from .variant_source import variant_source_record
        variant_source = variant_source_record(
            genotype=genotype_input,
            genotype_format=genotype_meta.get('genotype_format') or genotype_format,
            marker_ids=None if marker_ids is None else marker_names,
            variant_offset=0 if variant_range is None else _resolve_variant_range(variant_range, int(genotype.shape[1]))[0],
            n_variants=genotype_shape[1], load=genotype_load)
        # A genotype passed as an object has no path to record: keep its IDs.
        embed_variant_ids = bool(sumstats_variant_ids) or variant_source is None or variant_source['genotype'] is None
        if calibrator is not None:
            calibration_devices=variant_devices or trait_devices or [str(resolved_device)]
            calibrator.prepare(genotype,devices=calibration_devices,request=dict(
                capacity=int(chunk_size),prefetch_chunks=int(prefetch_chunks),reader_workers=int(reader_workers),
                genotype_shape=list(genotype.shape),phenotype_shape=list(phenotype.shape),
                covariate_shape=None if covariates is None else list(covariates.shape),
                variant_range=list(_resolve_variant_range(variant_range,int(genotype.shape[1]))),
                trait_block=trait_block,trait_devices=trait_devices,variant_devices=variant_devices,
                reduction='jagwas' if jagwas is not None else 'significant' if significance is not None else 'full',
                significance_threshold=None if significance is None else significance.resolved_threshold(len(trait_names)),
                p_value_threshold=p_value_threshold,
                sumstats_fields=sumstats_fields,sumstats_block_bytes=sumstats_block_bytes,
                sumstats_queue_depth=sumstats_queue_depth,sumstats_fsync=sumstats_fsync,
                sumstats_variant_ids=sumstats_variant_ids,pgen_decode_workers=pgen_decode_workers,
                pgen_decode_batch_size=pgen_decode_batch_size,trait_names=list(trait_names),
                covariate_columns=covariate_columns))
        if empirical is not None:
            from .empirical_autotune import EmpiricalChunkTuner
            from .adaptive_chunks import validate_chunk_control
            passes = 1 if trait_block is None else -(-int(phenotype.shape[1])//int(trait_block))
            options, sizes = empirical['options'], empirical['sizes']
            concurrent = len(variant_devices or trait_devices or [resolved_device])
            empirical['tuner_disabled_reason'] = None
            from .adaptive_chunks import chunk_control_path
            empirical['control_path'] = chunk_control_path(
                genotype, choose_device((variant_devices or trait_devices or [str(resolved_device)])[0]))
            if len(sizes) > 1:
                try:
                    tuner_kind = options.get('tuner', 'model')
                    if tuner_kind == 'model':
                        # Per-chunk samples, load-adjusted, re-planned on drift.
                        from .model_autotune import ModelChunkTuner
                        candidate = ModelChunkTuner(
                            sizes, total_rows=int(genotype_shape[1])*passes,
                            warmup_fraction=float(options.get('warmup_fraction', 0.02)),
                            max_probe_share=float(options.get('trial_fraction', 0.3)),
                            probe_chunks=int(options.get('probe_chunks', 4)),
                            margin=float(options.get('min_gain', 0.03)),
                            depth=int(prefetch_chunks), concurrent=concurrent,
                            min_job_seconds=float(options.get('min_job_seconds', 20.0)))
                    elif tuner_kind == 'segments':
                        candidate = EmpiricalChunkTuner(
                            sizes, total_rows=int(genotype_shape[1])*passes,
                            warmup_fraction=float(options.get('warmup_fraction', 0.03)),
                            trial_fraction=float(options.get('trial_fraction', 0.25)),
                            repeats=int(options.get('repeats', 2)), min_gain=float(options.get('min_gain', 0.02)),
                            depth=int(prefetch_chunks), concurrent=concurrent,
                            min_job_seconds=float(options.get('min_job_seconds', 20.0)))
                    else:
                        raise ValueError("autotune_options['tuner'] must be 'model' or 'segments'")
                    for tuned_device in (variant_devices or trait_devices or [str(resolved_device)]):
                        validate_chunk_control(genotype, choose_device(tuned_device), resolved_compute_dtype,
                                               int(chunk_size), candidate.control, candidate)
                    empirical_tuner = candidate
                except ValueError as error:
                    empirical['tuner_disabled_reason'] = str(error)
            else:
                empirical['tuner_disabled_reason'] = 'single chunk size'
            if empirical_tuner is None:
                # Without switchable chunks, use one mid-range size rather
                # than the largest candidate the ring was sized for.
                chunk_size = effective_chunk_size = sizes[(len(sizes)-1)//2]
                qc['genotype_qc_chunk_size'] = chunk_size
        output_progress = (productive.output_written if productive is not None else
                           None if calibrator is None else calibrator.output_written)
        # One decode pass shared by concurrently scanning phenotype tiles
        # (shared_decode.py): on for autotuned runs, or TORCHGWAS_SHARED_DECODE=1.
        import os as _os
        import threading
        # Significant pairs are selected on the GPU unless the threshold passes a
        # large share of cells (reduce.default_significance_backend; layout_profile
        # 2026-09-24: never slower, 10x faster on a CPU-loaded host). An explicit
        # TORCHGWAS_SIGNIFICANCE_BACKEND still wins.
        significance_backend = None
        if significance is not None:
            from .significance_backend import default_significance_backend
            significance_backend = _os.environ.get('TORCHGWAS_SIGNIFICANCE_BACKEND') or default_significance_backend(
                significance, int(np.shape(phenotype)[1]))
            genotype_meta['significance_backend'] = significance_backend
            if 'TORCHGWAS_SIGNIFICANCE_BACKEND' in _os.environ:
                significance_backend = None  # linear reads the variable itself
        shared_state = dict(hub=None, enabled=False)
        # More phenotype tiles than GPUs means rounds, and every round decodes
        # the genotype again. TORCHGWAS_GENOTYPE_CACHE=1 decodes once into
        # host memory instead. Opt-in: with hardcall PGEN (cheap decode,
        # GPU-bound rounds) it measured no gain and cost the cache fills
        # (docs/autotune_design_20260924.md); it is for expensive decoders.
        tile_source = genotype
        if (trait_block is not None and productive is None and calibrator is None
                and _os.environ.get('TORCHGWAS_GENOTYPE_CACHE', '0') == '1'):
            tile_count = -(-int(phenotype.shape[1])//int(trait_block))
            if tile_count > len(trait_devices or [resolved_device]):
                from .genotype_cache import CachedFillSource, fill_cache_for
                fill_cache = fill_cache_for(genotype, _resolve_variant_range(variant_range, int(genotype.shape[1])))
                if fill_cache is not None:
                    tile_source = CachedFillSource(genotype, fill_cache)
                    shared_state['genotype_cache'] = fill_cache
        if trait_block is not None and productive is None and calibrator is None and (
                (empirical is not None and empirical['options'].get('shared_decode', True))
                or _os.environ.get('TORCHGWAS_SHARED_DECODE', '0') == '1'):
            from .shared_decode import shared_decode_eligible
            shared_devices = trait_devices or [str(resolved_device)]
            shared_tiles = -(-int(phenotype.shape[1])//int(trait_block))
            shared_state['enabled'] = shared_decode_eligible(
                genotype, devices=shared_devices, tiles=shared_tiles, compute_dtype=resolved_compute_dtype)
            shared_state.update(tiles=shared_tiles, lock=threading.Lock(), fanout=None)
            if shared_state['enabled']:
                # GPU fan-out: one PCIe copy to the first tile GPU, then peer
                # copies. 'auto' and 'pcie' copy from host to every GPU. Fan-out
                # doubles raw delivery on lab-a100 pairs that share a PCIe
                # uplink (probe: 6.7 versus 3.3 GB/s each), yet end to end it
                # tied at 2 tiles and lost at 4 there and on lab-h100: host
                # transfer was not the bottleneck. 'uplink' fans out only over
                # NVLink between GPUs sharing an uplink; 'nvlink' whenever
                # NVLink connects them; 'peer' always.
                from .shared_decode import fanout_root, nvlink_root, root_slots
                fanout_mode = ((empirical['options'].get('gpu_fanout') if empirical is not None else None)
                               or _os.environ.get('TORCHGWAS_GPU_FANOUT', 'auto'))
                if fanout_mode not in ('auto', 'pcie', 'uplink', 'nvlink', 'peer'):
                    raise ValueError("gpu_fanout must be 'auto', 'pcie', 'uplink', 'nvlink' or 'peer'")
                tile_devices = shared_devices[:shared_tiles]
                root = (str(tile_devices[0]) if fanout_mode == 'peer' else
                        nvlink_root(tile_devices) if fanout_mode == 'nvlink' else
                        fanout_root(tile_devices) if fanout_mode == 'uplink' else None)
                shared_state.update(fanout_mode=fanout_mode, fanout_reason=None if root else 'not selected')
                if root is not None and empirical is not None:
                    # The root keeps an extra ring of chunks beside its own tile's rings.
                    from .empirical_autotune import gpu_free_bytes
                    from .pipeline_model import device_ring_bytes
                    per_variant = (empirical.get('layout') or {}).get('per_variant_transfer_bytes',
                                                                      float(genotype.shape[0])*4.0)
                    slots = root_slots(max(int(prefetch_chunks), min(int(reader_workers), 16)))
                    need = slots*int(chunk_size)*per_variant + device_ring_bytes(
                        chunk_variants=int(chunk_size), depth=int(prefetch_chunks), n_samples=int(genotype.shape[0]),
                        n_traits=int(trait_block), covariate_rank=int(0 if covariates is None else np.shape(covariates)[1]),
                        transfer_bytes_per_variant=per_variant)
                    free = gpu_free_bytes([root]).get(root)
                    if free is not None and need > 0.85*free:
                        shared_state['fanout_reason'] = (f'root GPU memory: needs {need/2**30:.1f} GiB, '
                                                         f'{free/2**30:.1f} GiB free')
                        root = None
                shared_state['fanout'] = root

        def shared_subscriber(first):
            """This tile's share of the decode hub, or None when not sharing."""
            if not shared_state['enabled']:
                return None
            from .shared_decode import SharedDecodeHub
            with shared_state['lock']:
                if shared_state['hub'] is None:
                    shared_state['hub'] = SharedDecodeHub(
                        genotype, int(chunk_size), max(int(prefetch_chunks), min(int(reader_workers), 16)),
                        int(reader_workers), subscribers=shared_state['tiles'], variant_range=variant_range,
                        chunk_size_selector=None if empirical_tuner is None else empirical_tuner.control,
                        record_timing=empirical_tuner is not None, fanout_device=shared_state['fanout'])
                return shared_state['hub'].subscriber(first//int(trait_block))
        if full_tiled_output or full_variant_output:
            from .preprocess import _covariate_basis
            from .sumstats_tiled import write_trait_tiled_sumstats,ScanSourceView
            q_matrix = (autotune_basis if autotune_basis is not None else
                        None if covariates is None or covariates.shape[1]==0 else _covariate_basis(covariates))
            rank = 0 if q_matrix is None else q_matrix.shape[1]
            # A panel with missing phenotypes keeps the single-device contract:
            # t from the mean-imputed panel with per-trait df, no variant df
            # sidecar. Tiling or sharding must not change the p-value convention.
            partition_trait_df = (None if not qc['phenotype_missing_cells'] else
                (np.asarray(qc['phenotype_observed_counts'], dtype=np.int64) - rank - 2).tolist())
            offset = 0 if variant_range is None else _resolve_variant_range(variant_range,genotype.shape[1])[0]
            tile_qc = {}
            tile_observed_counts = np.asarray(qc['phenotype_observed_counts'], dtype=np.int64)

            def scan_tile(first,width,tile_device,workers,_variant_span=None):
                # Take this tile's decode share first and give it back if setup
                # fails, so a failing tile cannot hold the shared decoder up.
                shared=shared_subscriber(first) if _variant_span is None else None
                try:
                    yield from _scan_tile(first,width,tile_device,workers,_variant_span,shared)
                finally:
                    if shared is not None:shared.close()  # idempotent after a normal finish

            def _scan_tile(first,width,tile_device,workers,_variant_span,shared):
                # Readers own their handles; metadata and input mappings are
                # read-only. Isolate per-scan profiles/exclusion counters.
                source=ScanSourceView(tile_source)
                if hasattr(source,'decode_workers'):
                    source.decode_workers=workers
                tile_phenotype=np.asarray(phenotype[:,first:first+width],
                                         dtype=np.float32 if resolved_compute_dtype=='float32' else np.float64)
                iterator,_=linear_scan_streaming_chunks(source,tile_phenotype,covariates,
                    chunk_size=chunk_size,device=str(tile_device),compute_dtype=resolved_compute_dtype,
                    reader_workers=workers,prefetch_chunks=prefetch_chunks,compute_p_values=False,
                    compute_log10_p=True,log10_p_dtype='float32',
                    variant_range=variant_range if _variant_span is None else _variant_span,borrow_results=True,return_df=True,_reader_worker_limit=workers,
                    _prevalidated_observed_counts=tile_observed_counts[first:first+width],
                    _prevalidated_covariate_basis=q_matrix,
                    _chunk_size_selector=(empirical_tuner.control if empirical_tuner is not None else
                        None if productive is None else productive.for_partition(
                        str(tile_device),_resolve_variant_range(variant_range if _variant_span is None else _variant_span,int(genotype.shape[1])),
                        (first,first+width))),
                    _chunk_observer=(empirical_tuner if empirical_tuner is not None else
                        calibrator.observer(str(tile_device),
                        variant_range=_resolve_variant_range(variant_range if _variant_span is None else _variant_span,int(genotype.shape[1])),
                        trait_range=(first,first+width),reader_workers=workers,capacity=chunk_size,depth=prefetch_chunks)
                        if calibrator is not None else None if productive is None else
                        productive.stage_observer(str(tile_device),
                            _resolve_variant_range(variant_range if _variant_span is None else _variant_span,int(genotype.shape[1])),
                            (first,first+width))),
                    return_beta=sumstats_fields!='t',
                    _shared_loader=shared)
                try:
                    for chunk in iterator:
                        chunk_offset=offset if _variant_span is None else _variant_span[0]
                        yield (chunk[0]-chunk_offset,chunk[1]-chunk_offset,*chunk[2:])
                    tile_qc[first if _variant_span is None else _variant_span[0]]=dict(getattr(source,'_last_scan_exclusion_counts',{}))
                finally:
                    iterator.close()

            out=mkdir(output_dir)
            tile_metadata={}
            def finalize_tiles():
                if tile_qc:
                    first_qc=tile_qc[min(tile_qc)]
                    if full_variant_output:
                        if any(set(value)!=set(first_qc) for value in tile_qc.values()):
                            raise RuntimeError('Genotype exclusion categories differ across variant shards')
                        genotype._last_scan_exclusion_counts={key:sum(value[key] for value in tile_qc.values()) for key in first_qc}
                    else:
                        if any(value!=first_qc for value in tile_qc.values()):
                            raise RuntimeError('Genotype exclusion counts differ across phenotype tiles')
                        genotype._last_scan_exclusion_counts=first_qc
                if embed_variant_ids:
                    ids_path=out/'sumstats'/'variant_ids.txt'
                    tile_metadata['variant_id_bytes']=_write_variant_ids(ids_path,marker_names)
                    if sumstats_fsync:
                        import os
                        with ids_path.open('rb') as handle:
                            os.fsync(handle.fileno())
            write_started=time.perf_counter()
            writer_settings=dict(n_variants=genotype_shape[1],trait_names=trait_names,
                n_samples=genotype_shape[0],df=genotype_shape[0]-rank-2,
                reader_workers=int(reader_workers),block_bytes=sumstats_block_bytes,
                queue_depth=sumstats_queue_depth or DEFAULT_SUMSTATS_QUEUE_DEPTH,fsync=sumstats_fsync,
                store_beta=sumstats_fields!='t',before_publish=finalize_tiles,
                on_write_progress=output_progress,
                on_writer_open=None if productive is None else productive.register_writer,
                trait_df=partition_trait_df,
                extra_manifest=None if variant_source is None else dict(variant_source=variant_source))
            if full_variant_output:
                from .sumstats_sharded import write_variant_sharded_sumstats
                def scan_shard(first,last,shard_device,workers):
                    return scan_tile(0,phenotype.shape[1],shard_device,workers,
                                     _variant_span=(offset+first,offset+last))
                sumstats_summary=write_variant_sharded_sumstats(out/'sumstats',
                    chunk_size=int(effective_chunk_size),devices=variant_devices,scan_factory=scan_shard,**writer_settings)
            else:
                sumstats_summary=write_trait_tiled_sumstats(out/'sumstats',trait_block=int(trait_block),
                    devices=trait_devices or [str(resolved_device)],scan_factory=scan_tile,**writer_settings)
            sumstats_summary.update(tile_metadata)
            sumstats_summary.update(streaming_timing(
                _phase_prep_done, write_started, time.perf_counter()))
            n_rows=sumstats_summary['cells']
            table=[]
            beta=t_stat=p_value=None
        elif output_dir is not None:
            # Who may borrow the result ring instead of being handed copies.
            # Audited individually, because getting this wrong corrupts output
            # silently rather than raising:
            #   safe   `_drain_linear_chunks`        reads shapes, retains nothing
            # Significant selection gathers owned arrays before advancing the
            # source, so it may also borrow the dense native result ring.
            # Discard-only scans may borrow result buffers. Binary writers
            # retain owned arrays until their background writes complete.
            borrow_results = significance is not None or (sumstats_format == "none"
                              and trait_block is None)
            full_binary_df = (sumstats_format == 'binary' and reduction is None
                              and significance is None
                              and p_value_threshold is None and trait_block is None)
            dense_log10_p = (sumstats_format == 'binary' and reduction is None
                             and significance is None and p_value_threshold is None)

            def _scan(trait_slice=None, device=None, *, source=None, workers=None,
                      observed_counts=None, basis=..., shared_loader=None):
                panel = phenotype if trait_slice is None else phenotype[:, trait_slice]
                if significance is not None:
                    # QC can preserve a mmap or lazy column subset. Convert
                    # only this tile to the actual compute precision.
                    panel = np.asarray(panel, dtype=np.float32 if resolved_compute_dtype == 'float32' else np.float64)
                return linear_scan_streaming_chunks(
                    genotype if source is None else source,
                    panel,
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
                    # `_p_chunk` and ignores it, and the writer recomputes p
                    # for the handful of pairs that clear the threshold. It
                    # presented as a multi-GPU problem (one core pinned, all
                    # cards idle, disk idle, worse as K grows, unaffected by
                    # running one process per card) and it was not one.
                    compute_p_values=False,
                    # Dense binary output stores -log10 P, computed on the device.
                    compute_log10_p=dense_log10_p, log10_p_dtype='float32',
                    variant_range=variant_range,
                    reduction=reduction,
                    borrow_results=borrow_results,
                    return_df=full_binary_df or significance is not None,
                    significance=significance,
                    significance_n_traits=int(phenotype.shape[1]),
                    _reader_worker_limit=workers,
                    _prevalidated_observed_counts=observed_counts,
                    _prevalidated_covariate_basis=basis,
                    _chunk_size_selector=(empirical_tuner.control if empirical_tuner is not None else
                        None if productive is None else productive.for_partition(
                        str(device or resolved_device),_resolve_variant_range(variant_range,int(genotype.shape[1])),
                        (0,int(phenotype.shape[1])) if trait_slice is None else (trait_slice.start,trait_slice.stop))),
                    _chunk_observer=(empirical_tuner if empirical_tuner is not None else
                        calibrator.observer(str(device or resolved_device),
                        variant_range=_resolve_variant_range(variant_range,int(genotype.shape[1])),
                        trait_range=(0,int(phenotype.shape[1])) if trait_slice is None else (trait_slice.start,trait_slice.stop),
                        reader_workers=reader_workers if workers is None else workers,capacity=chunk_size,depth=prefetch_chunks)
                        if calibrator is not None else None if productive is None else
                        productive.stage_observer(str(device or resolved_device),
                            _resolve_variant_range(variant_range,int(genotype.shape[1])),
                            (0,int(phenotype.shape[1])) if trait_slice is None else (trait_slice.start,trait_slice.stop))),
                    return_beta=(sumstats_fields!='t' or reduction is not None or p_value_threshold is not None),
                    _shared_loader=shared_loader,
                    _significance_backend=significance_backend,
                )

            blocked_significance = False
            finalize_reduced = None
            if jagwas_variant_output or significant_variant_output:
                from .linear import linear_scan_multigpu, multigpu_variant_ranges
                ranges = multigpu_variant_ranges(genotype_shape[1], int(effective_chunk_size), len(variant_devices))
                active_variant_devices = variant_devices[:len(ranges)]
                # JAGWAS keeps a joint factor per device; significant pairs are
                # selected cell by cell on each device (stateless).
                if jagwas_variant_output:
                    # Shards share the run's kept-trait decision (collinear traits dropped once).
                    shard_output = dict(reduction_factory=jagwas.spawn,
                                        column_groups=getattr(jagwas, 'column_groups', None))
                else:
                    from .linear import _significant_pairs_iterator
                    # Select pairs on each shard's thread, as trait tiles do;
                    # the chunks carry their own df (return_df).
                    shard_output = dict(significance=significance, significance_n_traits=int(phenotype.shape[1]),
                                        return_df=True, return_beta=sumstats_fields != 't',
                                        _significance_backend=significance_backend,
                                        shard_transform=lambda chunks: _significant_pairs_iterator(
                                            chunks, significance, len(trait_names), None),
                                        # Selection copies the passing pairs before the next
                                        # chunk, so it can read the scan's ring directly.
                                        transform_borrows=True)
                chunk_iterator, q_matrix = linear_scan_multigpu(
                    genotype, phenotype, covariates, devices=variant_devices,
                    chunk_size=int(effective_chunk_size), compute_dtype=resolved_compute_dtype,
                    reader_workers=int(reader_workers), prefetch_chunks=prefetch_chunks,
                    compute_p_values=False, ordered=False, variant_range=variant_range,
                    **shard_output,
                    _chunk_size_selector=(empirical_tuner.control if empirical_tuner is not None else
                        None if productive is None else productive.run),
                    _chunk_observer=empirical_tuner if empirical_tuner is not None else
                        calibrator if calibrator is not None else
                        productive if productive is not None and productive.stage_sample is not None else None,
                    shared_queue_depth=DEFAULT_SUMSTATS_QUEUE_DEPTH if sumstats_queue_depth is None else sumstats_queue_depth,
                    result_queue_registration=None if productive is None or not jagwas_variant_output
                        else productive.register_indexed_result_queue)
            elif trait_block is None:
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

                blocked_significance = significance is not None
                if blocked_significance:
                    from .sumstats_tiled import ScanSourceView
                    from .reduced_output_work import significant_execution_layout
                    layout = significant_execution_layout(int(phenotype.shape[1]), int(trait_block),
                        trait_devices or [str(resolved_device)], int(reader_workers),
                        queue_depth=sumstats_queue_depth)
                    reader_budgets = dict(zip(layout['devices'], layout['readers_per_device']))
                    tile_records = {}
                    observed = np.asarray(qc['phenotype_observed_counts'], dtype=np.int64)

                    def _blocked(offset, width, device=None):
                        # Take the decode share first; give it back on any exit.
                        shared = shared_subscriber(offset)
                        iterator = None
                        try:
                            dev = str(device or layout['devices'][0])
                            source = ScanSourceView(tile_source)
                            workers = reader_budgets[dev]
                            if hasattr(source, 'decode_workers'):
                                source.decode_workers = workers
                            iterator, _ = _scan(slice(offset, offset + width), dev,
                                source=source, workers=workers,
                                observed_counts=observed[offset:offset + width], basis=q_matrix,
                                shared_loader=shared)
                            yield from iterator
                            tile_records[offset] = dict(trait_range=[offset, offset + width],
                                device=dev, reader_workers=workers,
                                exclusions=dict(getattr(source, '_last_scan_exclusion_counts', {})),
                                profile=dict(getattr(source, '_last_scan_profile', {})))
                        finally:
                            try:
                                if iterator is not None:
                                    iterator.close()
                            finally:
                                if shared is not None:
                                    shared.close()

                    def finalize_reduced():
                        records = [tile_records[key] for key in sorted(tile_records)]
                        if len(records) != layout['tiles']:
                            raise RuntimeError('Incomplete significant phenotype tiles')
                        counts = records[0]['exclusions']
                        if any(record['exclusions'] != counts for record in records):
                            raise RuntimeError('Genotype exclusion counts differ across phenotype tiles')
                        genotype._last_scan_exclusion_counts = counts
                        profile = dict(mode='significant_trait_tiles', layout=layout, tiles=records,
                            scope='Per-tile profiles; durations can overlap and are not total wall time')
                        if all('result_payload_bytes' in record['profile'] for record in records):
                            profile['result_payload_bytes'] = sum(record['profile']['result_payload_bytes'] for record in records)
                        genotype._last_scan_profile = profile
                else:
                    def _blocked(offset, width, device=None):
                        iterator, _ = _scan(slice(offset, offset + width), device)
                        return iterator
                if blocked_significance:
                    from .sumstats_indexed import IndexedOutputPartition
                    chunk_iterator = _trait_blocked_significant_chunks(
                        _blocked, significance, int(phenotype.shape[1]),
                        int(trait_block),
                        genotype_shape[0] - (0 if q_matrix is None
                                             else q_matrix.shape[1]) - 2,
                        devices=layout['devices'], queue_depth=layout['queue_depth'],
                        partition_context=None if calibrator is None and productive is None else lambda first,width,dev: IndexedOutputPartition(
                            str(dev or layout['devices'][0]),
                            tuple(_resolve_variant_range(variant_range,int(genotype.shape[1]))),(first,first+width)))
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
                    from .sumstats_indexed import PartitionedIndexedChunk
                    try:
                        for chunk in source:
                            # Dense chunks with -log10 P and df are seven long
                            # too; only selected pairs carry variant indices.
                            if len(chunk) == 7 and significance is not None:
                                payload=(chunk[0]-offset,chunk[1]-offset,chunk[2]-offset,*chunk[3:])
                            else:
                                payload=(chunk[0]-offset,chunk[1]-offset,*chunk[2:])
                            yield chunk.with_payload(payload) if isinstance(chunk,PartitionedIndexedChunk) else payload
                    finally:
                        if hasattr(source,'close'):source.close()

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
                if finalize_reduced is not None:
                    finalize_reduced()
            elif (significance is not None or jagwas is not None or reduction is not None
                  or p_value_threshold is not None):
                from .sumstats_indexed import write_indexed_sumstats,IndexedOutputPartition,COALESCE_BYTES
                kind = ("significant" if significance is not None else "jagwas" if jagwas is not None
                        else "reduced" if reduction is not None else "filtered")
                indexed_observed=output_progress is not None and kind in ('jagwas','significant')
                source_start,source_stop=_resolve_variant_range(variant_range,int(genotype.shape[1]))
                partition_for_range=None
                if indexed_observed and not blocked_significance:
                    partitions=([IndexedOutputPartition(dev,(source_start+lo,source_start+hi),(0,len(trait_names)))
                        for dev,(lo,hi) in zip(active_variant_devices,ranges)]
                        if jagwas_variant_output or significant_variant_output else
                        [IndexedOutputPartition(str(resolved_device),(source_start,source_stop),(0,len(trait_names)))])
                    def partition_for_range(start,end):
                        matches=[p for p in partitions if p.variant_range[0]<=start<end<=p.variant_range[1]]
                        if len(matches)!=1:raise ValueError('Indexed source range must identify one producer')
                        return matches[0]
                indexed_live_progress=None
                if productive is not None and kind in ('jagwas','significant') and indexed_observed:
                    from .indexed_writer_progress import IndexedWriterProgress
                    indexed_live_progress=IndexedWriterProgress(kind)
                    productive.register_indexed_writer(indexed_live_progress)
                n_rows, sumstats_summary = write_indexed_sumstats(
                    out / "sumstats", marker_names, trait_names, genotype_shape[0], chunk_iterator,
                    kind=kind, df=residual_df,
                    # The kept rank is known once the factors are prepared, before the manifest.
                    chi2_df=(lambda: jagwas.degrees_of_freedom) if jagwas is not None else None,
                    extra_manifest=((lambda: dict(jagwas_rank=jagwas.rank_report(trait_names),
                                                  **({} if jagwas_groups is None else dict(groups=jagwas.names))))
                                    if jagwas is not None else None),
                    p_value_threshold=p_value_threshold,
                    variant_metadata=variant_metadata,fsync=sumstats_fsync,
                    store_beta=sumstats_fields != "t", before_publish=finalize_reduced,
                    on_chunk_written=output_progress if indexed_observed else None,
                    partition_for_range=partition_for_range,variant_offset=source_start,
                    live_progress=indexed_live_progress,
                    # Large parts unless the JIT path observes each chunk's write.
                    coalesce_bytes=None if indexed_observed else COALESCE_BYTES,
                    variant_source=variant_source, embed_variant_ids=embed_variant_ids)
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
                    write_variant_ids=embed_variant_ids,
                    store_beta=sumstats_fields != "t",
                    extra_manifest={"genotype_format": genotype_meta.get("genotype_format"),
                                    **({} if variant_source is None else dict(variant_source=variant_source))},
                    borrow_results=borrow_results,
                    store_variant_df=full_binary_df and bool(np.all(trait_df == residual_df)),
                    on_write_progress=output_progress,
                )
            sumstats_summary.update(streaming_timing(
                _phase_prep_done, write_started, time.perf_counter()))
            if blocked_significance:
                sumstats_summary['execution_layout'] = layout
            if jagwas is not None:
                sumstats_summary['jagwas_rank'] = jagwas.rank_report(trait_names)
            if jagwas_variant_output or significant_variant_output:
                sumstats_summary.update(devices=active_variant_devices, genotype_passes=1,
                    execution_layout=dict(partition_axis='variant', variant_ranges=[list(span) for span in ranges],
                        devices=active_variant_devices, reader_workers=int(reader_workers),
                        shared_result_queue_depth=(DEFAULT_SUMSTATS_QUEUE_DEPTH if sumstats_queue_depth is None else sumstats_queue_depth)
                            if len(active_variant_devices) > 1 else 0,
                        reduction_state=('independent factor per active device' if jagwas_variant_output
                                         else 'stateless per-cell selection'), writer_workers=1))
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
        "p_value_threshold": p_value_threshold,
        "reduce": ("jagwas" if jagwas is not None
                   else "significant" if significance is not None
                   else (None if reduction is None else reduction.mode)),
        "significance_threshold": (
            None if significance is None
            else significance.resolved_threshold(len(trait_names))),
        "reduce_top_k": None if reduction is None else reduction.width,
        "jagwas_groups": None if jagwas_groups is None else jagwas.names,
        # host_pages.py: NumPy's MADV_HUGEPAGE advice during the run.
        "numpy_hugepage_advice": numpy_hugepage_advice(),
        "phenotype_outlier_sd": phenotype_outlier_sd,
        "phenotype_outlier_rows": None if outlier_rows is None else int(outlier_rows.sum()),
        # Recorded because it changes how many passes the run made over the
        # genotypes, which is the first thing to check against a wall time.
        "trait_block": None if trait_block is None else int(trait_block),
        # WHICH devices the blocks ran on, not just how wide the blocks were.
        # The 600,000-trait stress run recorded `trait_block: 75000` and
        # nothing else, so "did it use eight cards or one?" could only be
        # argued from the fact that 75,000 is 600,000/8 -- an inference, when
        # the run already knew the answer and simply never wrote it down.
        "variant_devices": (sumstats_summary.get('devices') if full_variant_output or jagwas_variant_output
                            or significant_variant_output else None),
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
    shared = locals().get('shared_state')
    if shared is not None and shared.get('genotype_cache') is not None:
        genotype_meta['genotype_cache'] = shared['genotype_cache'].audit()
    if shared is not None and shared.get('hub') is not None:
        shared['hub'].abandon_untaken()  # scans are done; release any tile never started
        genotype_meta['shared_decode'] = dict(shared['hub'].audit(), enabled=True,
                                              fanout_mode=shared.get('fanout_mode'),
                                              fanout_reason=shared.get('fanout_reason'))
    if empirical is not None:
        genotype_meta['autotune'] = dict(
            method='empirical', chunk_sizes=empirical.get('sizes'), layout=empirical.get('layout'),
            devices_seen=empirical.get('devices_seen'),
            tuner_disabled_reason=empirical.get('tuner_disabled_reason'),
            control_path=empirical.get('control_path'),
            chunk=None if empirical_tuner is None else empirical_tuner.audit())
    if hasattr(genotype, 'backend_used'):
        genotype_meta['decode_backend_used'] = genotype.backend_used
        if hasattr(genotype, 'backend_reason'):
            genotype_meta['decode_backend_reason'] = genotype.backend_reason
    run_metadata.update(genotype_meta)
    if isinstance(genotype, ChunkedGenotype) and 'reader_workers' in genotype_meta:
        # The loader records the worker count it opened the source with; a
        # planner (empirical or detailed) may have changed it for the scan.
        run_metadata['genotype_open_reader_workers'] = genotype_meta['reader_workers']
        run_metadata['reader_workers'] = int(reader_workers)
    result = GWASResult(table=table, run_metadata=run_metadata, qc_summary=qc)
    if output_dir is not None:
        out = mkdir(output_dir)
        streamed = isinstance(genotype, ChunkedGenotype) and run_metadata["results_streamed"]
        if not streamed:
            if sumstats_format == "binary" and (p_value_threshold is not None):
                from .sumstats_indexed import write_indexed_sumstats
                written, sumstats_summary = write_indexed_sumstats(
                    out / "sumstats", marker_names, trait_names, genotype_shape[0],
                    iter([(0, len(marker_names), beta, t_stat, None)]), kind="filtered",
                    df=genotype_shape[0] - covariate_rank_used - 2,
                    p_value_threshold=p_value_threshold,
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
                    # Rows follow the QC-kept markers, not a contiguous input
                    # range, so the IDs are stored.
                    write_variant_ids=True,
                    store_beta=sumstats_fields != "t",
                )
                run_metadata["sumstats_write"] = sumstats_summary
        write_json(run_metadata, out / "run.json")
        write_json(qc, out / "qc.json")
    if calibrator is not None:
        # No record is published if the association writer or required run/QC
        # sidecars fail. Cache I/O failure is reported without losing results.
        run_metadata['initial_calibration']=calibrator.finish(successful=True)
        try:
            write_json(run_metadata,out/'run.json')
        except OSError as error:
            run_metadata['initial_calibration']['report_write_error']=str(error)
    if productive is not None:
        productive.finish(successful=True)
        try:
            write_json(run_metadata,out/'run.json')
        except OSError as error:
            genotype_meta['autotune']['productive']['report_write_error']=str(error)
    return result
