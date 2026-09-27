"""Measure this machine, so the planner can predict it.

This ships WITH the model. `pipeline_model.estimate()` is structural -- its
terms come from the code and the workload shape -- but every term is divided by
a machine rate, and those rates are properties of the host: disk bandwidth,
PCIe, HBM, achievable FP32. On a different box they all change, and the ratios
between them change too, which moves the roofline crossover and therefore which
resource binds.

So a calculator alone is not shippable. What ships is *calculator plus
calibration*: measure the host once, cache the profile, then plan against it.

**Why measured and not looked up.** A spec sheet gives peak numbers that no real
pipeline reaches, and the gap is not a constant -- it depends on transfer size,
on whether memory is pinned, and on the GEMM's shape. Feeding peak FP32 into the
model is what made a hand-estimate predict 60-100 s for a scan that measured
24.26 s: the statistics kernel reaches about 5.6 TFLOPS on 155 design columns
and 29.3 TFLOPS on 540, neither of which is peak.

**And FP32 is not one number -- prefer `measure_gemm_curve`.** A scalar is only
right at the width it was measured at. Measured with CUDA events around the
real kernel: 8.35 TFLOPS at design width 36, 8.71 at 60, 20.57 at 156, 20.28 at
540, a 2.5x spread that the roofline predicts exactly (arithmetic intensity is
about (K + C) / 2, so a narrow design is memory-bound). The model carried one
44.88 TFLOPS calibrated at width 540 and so over-estimated the GPU roughly 2x
at K = 128 and 5x at K = 8. That error survived a long time because the GEMM
rarely binds, so it never showed up end to end -- which is the argument for
checking each term against its own clock rather than checking a total.

**Why not fitted to scan times.** Nothing in this module times a GWAS scan. If
the rates were fitted to scans, the model would reproduce the runs it was fitted
to and predict nothing else -- and a cross-machine check would be circular. Each
rate is an independent microbenchmark, so predicting a scan is a real test.
"""
from __future__ import annotations

import json
import os
import statistics
import tempfile
import time
from pathlib import Path

NEWLINE = bytes([10])
TAB = bytes([9])


def _resident_fraction(path, span_bytes: int) -> float:
    """Share of the first `span_bytes` of `path` still in the page cache.

    NaN when it cannot be determined, which must never be reported as zero --
    "unknown" and "cold" are exactly the distinction a cold-read probe turns
    on, and conflating them is what let a warm measurement pass as cold.
    """
    try:
        import ctypes
        import mmap

        span = min(int(span_bytes), Path(path).stat().st_size)
        if span <= 0:
            return float("nan")
        with open(path, "rb") as handle:
            # ACCESS_COPY, not ACCESS_READ: `ctypes.from_buffer` requires a
            # writable buffer to yield an address, and a read-only mapping
            # refuses. A private copy-on-write mapping still shares the page
            # cache, so mincore reports the real residency and nothing is
            # actually copied.
            mapped = mmap.mmap(handle.fileno(), span, access=mmap.ACCESS_COPY)
            try:
                page = os.sysconf("SC_PAGE_SIZE")
                pages = (span + page - 1) // page
                vector = (ctypes.c_ubyte * pages)()
                address = ctypes.addressof(
                    (ctypes.c_char * span).from_buffer(mapped))
                libc = ctypes.CDLL("libc.so.6", use_errno=True)
                if libc.mincore(ctypes.c_void_p(address),
                                ctypes.c_size_t(span), vector) != 0:
                    return float("nan")
                return sum(value & 1 for value in vector) / pages
            finally:
                mapped.close()
    except Exception:  # noqa: BLE001 - any failure means "unknown"
        return float("nan")


def _median_rate(sizes_bytes: float, seconds: list[float]) -> float:
    """Bytes per second at the MEDIAN time, never the best.

    Best-of-N reports the run where the machine happened to be quiet, which is
    not the rate a scan will see. This project has already retracted one
    conclusion that came from measuring a busy host; medians are the standing
    methodology.
    """
    return sizes_bytes / statistics.median(seconds) if seconds else 0.0


def measure_disk_read(path, block_bytes: int = 64 << 20, blocks: int = 64,
                      repeat: int = 2,
                      workers: tuple = (1, 16)) -> dict:
    """Cold sequential read rate, in bytes/second, at several concurrencies.

    Drops the file's cache between repeats with `posix_fadvise(DONTNEED)`,
    because a warm read measures RAM. If the cache cannot be dropped -- some
    network filesystems ignore it -- the result carries `cache_dropped: False`
    rather than being passed off as a cold number.

    **One stream is the wrong number.** The scan reads with `reader_workers`
    threads, and this storage delivers far more in aggregate than to a single
    reader: measured 1.52 GB/s single-stream against a scan that achieved
    4.02 GB/s (35.93 GB of PGEN in 8.94 s). A profile reporting only the
    single-stream rate under-predicts the read term by ~2.6x, which would make
    the model blame the wrong resource entirely. The planner should use the
    rate at the worker count it intends to run.

    **The span must be large, and this probe used to be far too small.** It read
    1.5 GB (24 blocks of 64 MB) and reported **9.00 GB/s** at 16 readers on a
    host where the scan sustains **5.56 GB/s** over 79 GB. That 1.62x was not
    the scan's inefficiency: a read-only probe over 20 GB gets **6.09 GB/s with
    the very same contiguous-span pattern the small probe used**, and 6.29 GB/s
    with the scan's chunked pattern, and 6.34 GB/s writing into pinned buffers.
    Pattern and pinning make no difference; the SPAN does. A gigabyte and a half
    is short enough to be flattered by drive-level caching and by incomplete
    eviction, and the optimistic number it produced propagated into a
    "the pipeline does not overlap, 1.45-1.84x available" conclusion that was
    entirely an artifact of it.

    **Controlled**, because "contention, not span" was the obvious competing
    explanation: the 20 GB probe ran while another user had a dozen processes
    in D state on the same filesystem. Re-run in exactly those conditions, all
    arms verified 0.0% resident, at 16 workers:

        1.5 GB (the original)   9.41 GB/s
        4 GB                    6.82 GB/s
        12 GB                   6.68 GB/s
        20 GB                   6.44 GB/s

    The small span still reports 9.41 under the same contention where 20 GB
    gives 6.44, so the span is the cause. Single-worker reads are flat at
    ~1.8 GB/s across every span, so whatever the short burst exploits is
    specific to concurrent reads -- a device queue or cache that a 1.5 GB burst
    fits inside and a sustained read does not.

    The default now reads 4 GB per arm, with fewer repeats and fewer worker
    counts so the total work is comparable. 4 GB is a compromise and still runs
    ~6% above the 20 GB asymptote; a caller that cares should pass more blocks.
    Residency is verified rather than assumed, for the same reason:
    `posix_fadvise` is advisory.
    """
    from concurrent.futures import ThreadPoolExecutor

    path = Path(path)
    dropped = True
    total = block_bytes * blocks
    # Never claim to have read more than exists; a short file would otherwise
    # be reported at the nominal span and at a wildly inflated rate.
    available = path.stat().st_size
    if total > available:
        blocks = max(1, available // block_bytes)
        total = block_bytes * blocks
    by_workers = {}
    residency = []

    def drop_cache():
        nonlocal dropped
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        except (AttributeError, OSError):
            dropped = False
            return
        finally:
            os.close(fd)
        # Verified, not assumed. `posix_fadvise(DONTNEED)` is advisory and
        # frequently does not evict; a warm arm reporting itself as cold is
        # how this probe produced 9.00 GB/s where the sustained rate is ~6.
        share = _resident_fraction(path, total)
        if share == share:
            residency.append(share)
            if share > 0.05:
                dropped = False

    def read_span(index, count):
        # A contiguous span per worker, so the streams do not interleave and
        # this measures parallel throughput rather than seek thrash.
        fd = os.open(path, os.O_RDONLY)
        try:
            span = max(blocks // count, 1)
            got = 0
            for b in range(index * span, min((index + 1) * span, blocks)):
                block = os.pread(fd, block_bytes, b * block_bytes)
                if not block:
                    break
                got += len(block)
            return got
        finally:
            os.close(fd)

    for count in workers:
        seconds = []
        for _ in range(repeat):
            drop_cache()
            started = time.perf_counter()
            if count == 1:
                read_span(0, 1)
            else:
                with ThreadPoolExecutor(max_workers=count) as pool:
                    list(pool.map(lambda i: read_span(i, count), range(count)))
            seconds.append(time.perf_counter() - started)
        by_workers[str(count)] = _median_rate(total, seconds)

    peak = max(by_workers.values()) if by_workers else 0.0
    return {"bytes_per_second": peak,
            "bytes_per_second_by_workers": by_workers,
            "single_stream_bytes_per_second": by_workers.get("1", 0.0),
            "cache_dropped": dropped, "bytes": total, "path": str(path),
            # The span is part of the result, not a footnote: the same host
            # reports 9.00 GB/s over 1.5 GB and ~6.1 GB/s over 20 GB, and a
            # consumer cannot tell those apart without knowing which was read.
            "span_bytes": total,
            "resident_fraction_before": (max(residency) if residency
                                         else float("nan"))}


def measure_write(directory: str | Path, block_bytes: int = 64 << 20,
                  blocks: int = 8) -> dict:
    """Durable sequential write rate, including the fsync.

    The fsync is not optional: without it this measures the page cache, and a
    profile that prices the sumstats write at cache speed will under-predict
    every end-to-end number. An earlier calibration harness made write rate a
    REQUIRED argument for exactly this reason.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    payload = b"\0" * block_bytes
    with tempfile.NamedTemporaryFile(dir=directory, delete=True) as handle:
        started = time.perf_counter()
        for _ in range(blocks):
            handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
        elapsed = time.perf_counter() - started
    return {"bytes_per_second": (block_bytes * blocks) / elapsed if elapsed else 0.0,
            "bytes": block_bytes * blocks, "seconds": elapsed, "fsynced": True}


def measure_text_parse(path, max_bytes: int = 64 << 20) -> dict:
    """Single-threaded text parse rate, bytes/second.

    This is the rate that sets `open` -- the variant-metadata parse. It is a
    real term and a large one: parsing the 236 MB `.pvar` costs ~8.5 s cold,
    which at M = 1,000,000 is **76% of the entire end-to-end run** and is
    completely flat in M. A model validated only at full genome would never see
    it, because there it is under 10% and hides inside the noise.

    Measured on the actual metadata file where possible, because the rate
    depends on line length and field count, not just on bytes. It is
    single-threaded on purpose: plink2 parses its `.pvar` serially and so do we,
    and pretending otherwise would price a serial tail at parallel speed.
    """
    path = Path(path)
    fd = os.open(path, os.O_RDONLY)
    try:
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        except (AttributeError, OSError):
            pass
        started = time.perf_counter()
        read = 0
        fields = 0
        while read < max_bytes:
            block = os.pread(fd, 1 << 20, read)
            if not block:
                break
            read += len(block)
            # Split the way a parser does, so this measures parsing and not
            # just reading: a pure read would report storage bandwidth and
            # under-price the term by more than an order of magnitude.
            for line in block.split(NEWLINE):
                fields += line.count(TAB)
        elapsed = time.perf_counter() - started
    finally:
        os.close(fd)
    return {"bytes_per_second": read / elapsed if elapsed else 0.0,
            "bytes": read, "seconds": elapsed, "fields_seen": fields,
            "path": str(path)}


def measure_transfers(device: int = 0, payload_bytes: int = 256 << 20,
                      repeat: int = 7) -> dict:
    """H2D, D2H and device-to-device bandwidth, from PINNED host memory.

    Pinned, because that is what the scan uses: a pageable transfer measures an
    extra staging copy the real pipeline does not pay, and the two differ by
    more than a factor of two on this class of hardware.
    """
    import torch

    dev = torch.device("cuda", device)
    torch.cuda.set_device(dev)
    elements = payload_bytes // 4
    host = torch.empty(elements, dtype=torch.float32, pin_memory=True)
    dev_a = torch.empty(elements, dtype=torch.float32, device=dev)
    dev_b = torch.empty(elements, dtype=torch.float32, device=dev)

    def timed(fn) -> list[float]:
        out = []
        for index in range(repeat + 1):
            torch.cuda.synchronize(dev)
            started = time.perf_counter()
            fn()
            torch.cuda.synchronize(dev)
            if index:           # discard the first, which pays warm-up
                out.append(time.perf_counter() - started)
        return out

    h2d = timed(lambda: dev_a.copy_(host, non_blocking=True))
    d2h = timed(lambda: host.copy_(dev_a, non_blocking=True))
    d2d = timed(lambda: dev_b.copy_(dev_a))
    return {
        "h2d_bytes_per_second": _median_rate(payload_bytes, h2d),
        "d2h_bytes_per_second": _median_rate(payload_bytes, d2h),
        # d2d moves the payload twice (read + write), which is the figure that
        # belongs in a bandwidth term.
        "hbm_bytes_per_second": _median_rate(2 * payload_bytes, d2d),
        "payload_bytes": payload_bytes,
    }


def measure_gemm(n_samples: int, n_traits: int, covariates: int,
                 chunk_variants: int = 4096, device: int = 0,
                 repeat: int = 7) -> dict:
    """Achievable FP32 at the SHAPE the scan runs, not at peak.

    The statistics kernel is `(chunk x N) @ (N x (K + C + 1))` -- tall and
    skinny. Its arithmetic intensity is roughly `(K + C) / 2`, so it sits below
    the device's ridge point at small K and above it at large K, and the
    achieved rate moves by a factor of five across the K range. Measuring at one
    shape and extrapolating to another is the error that produced a 2.5-4x
    over-estimate, so the shape is a parameter here.
    """
    import torch

    dev = torch.device("cuda", device)
    torch.cuda.set_device(dev)
    width = n_traits + covariates + 1
    genotype = torch.randn((chunk_variants, n_samples), device=dev)
    design = torch.randn((n_samples, width), device=dev)
    flops = 2.0 * chunk_variants * n_samples * width
    seconds = []
    for index in range(repeat + 1):
        torch.cuda.synchronize(dev)
        started = time.perf_counter()
        genotype @ design
        torch.cuda.synchronize(dev)
        if index:
            seconds.append(time.perf_counter() - started)
    achieved = flops / statistics.median(seconds) if seconds else 0.0
    return {"gemm_flops_per_second": achieved, "flops": flops,
            "shape": {"chunk_variants": chunk_variants, "n_samples": n_samples,
                      "design_columns": width},
            "arithmetic_intensity_flops_per_byte":
                flops / (chunk_variants * n_samples * 4.0
                         + n_samples * width * 4.0
                         + chunk_variants * width * 4.0)}


def measure_gemm_curve(n_samples: int, covariates: int = 27,
                       trait_counts: tuple = (8, 32, 128, 512),
                       chunk_variants: int = 4096, device: int = 0,
                       repeat: int = 3) -> dict:
    """Achieved FP32 across several design widths, not one.

    **A single number is wrong at every width but the one it was measured at.**
    Measured on the H100 with CUDA events around the real statistics kernel:
    8.35 TFLOPS at design width 36, 8.71 at 60, 20.57 at 156, 20.28 at 540 --
    a 2.5x spread. The roofline predicts exactly that, since arithmetic
    intensity is about (K + C) / 2: a narrow design is memory-bound and a wide
    one approaches the FLOP limit.

    The model used one scalar, 44.88 TFLOPS calibrated at width 540, and so
    over-estimated the GPU roughly 2x at K = 128 and 5x at K = 8. That error
    was invisible in end-to-end checking because the GEMM rarely binds -- which
    is the whole argument for validating term by term.

    Returns `{design_width: flops_per_second}`, which `explain_time` accepts
    directly in place of a scalar.
    """
    curve = {}
    detail = {}
    for traits in trait_counts:
        one = measure_gemm(n_samples, traits, covariates,
                           chunk_variants=chunk_variants, device=device,
                           repeat=repeat)
        width = traits + covariates + 1
        curve[width] = one["gemm_flops_per_second"]
        detail[str(width)] = {"traits": traits, **one}
    return {"gemm_flops_by_width": curve, "detail": detail,
            "covariates": covariates, "n_samples": n_samples}


def measure_pinning(bytes_to_pin: int = 1 << 30, repeat: int = 3) -> dict:
    """Rate at which the host can PIN memory, bytes/second.

    `cudaHostAlloc` is far slower than ordinary allocation -- the pages must be
    locked and registered with the driver -- and the scan pins its whole
    staging ring before steady state begins. That makes it a startup cost
    proportional to `depth * chunk * transfer_bytes_per_variant`, which is
    exactly the shape of the intercept the M sweep measured but the model could
    not explain:

        pgen  1.30 GB pinned  ->  2.10 s intercept
        bed   0.65 GB pinned  ->  0.91 s intercept

    a 2.0:1 ratio of bytes against a 2.3:1 ratio of intercepts. CUDA context
    creation, cuBLAS warm-up and thread-pool startup together measure 0.12 s,
    so they cannot account for it; pinning can.

    Measured rather than assumed because it varies by an order of magnitude
    with the host's memory pressure and hugepage configuration, which is
    precisely the sort of thing that differs between a developer box and a
    cluster node.
    """
    import torch

    # PyTorch CACHES pinned memory in its host allocator, so freeing a buffer
    # and allocating the same size again is served from the cache and times at
    # 0.0 s. Measuring that way reported 29,159 GB/s from times
    # [0.778, 0.0, 0.0] -- an impossible number that only looks plausible if
    # nobody checks the units. Each repeat therefore asks for a DIFFERENT size,
    # which the cache cannot satisfy, and the rate is taken per actual byte.
    rates = []
    held = []
    for index in range(repeat):
        size = bytes_to_pin + index * (bytes_to_pin // 8)
        started = time.perf_counter()
        buffer = torch.empty(size, dtype=torch.uint8, pin_memory=True)
        elapsed = time.perf_counter() - started
        # Keep every buffer alive: releasing returns it to the cache and the
        # next distinct size may then be served from a recycled block.
        held.append(buffer)
        if elapsed > 0:
            rates.append(size / elapsed)
    del held
    return {"bytes_per_second": statistics.median(rates) if rates else 0.0,
            "bytes": bytes_to_pin, "rates": rates,
            "note": "distinct sizes per repeat; torch caches pinned blocks"}


def measure_scan_startup(n_samples: int, n_traits: int = 128,
                         covariates: int = 27, chunk_variants: int = 4096,
                         reader_workers: int = 16, device: int = 0) -> dict:
    """One-time setup a scan pays before steady state, in seconds.

    The M sweep measures a per-variant slope AND an intercept: 2.10 s for PGEN,
    0.91 s for BED. The intercept is real and it is not pipeline geometry --
    `sum(resources) / chunks` gives 0.01 s at 2,181 chunks, three orders of
    magnitude too small. It is one-time setup: CUDA context creation, cuBLAS
    autotuning its first GEMM at this shape, and the reader thread pool.

    Measured here as components rather than fitted to the intercept, because
    fitting would make the model reproduce the sweep it was fitted to and
    predict nothing on a new machine -- and this term differs a lot between
    machines, since cuBLAS warm-up scales with the GPU and context creation
    with the driver.
    """
    import torch
    from concurrent.futures import ThreadPoolExecutor

    dev = torch.device("cuda", device)
    out = {}

    # CUDA context: the first real allocation forces it, so time that.
    started = time.perf_counter()
    torch.cuda.set_device(dev)
    _ = torch.zeros(1, device=dev)
    torch.cuda.synchronize(dev)
    out["cuda_context_seconds"] = time.perf_counter() - started

    # cuBLAS picks an algorithm on first use at a given shape; subsequent calls
    # reuse it. The scan pays this once, at the shape it runs.
    width = n_traits + covariates + 1
    genotype = torch.randn((chunk_variants, n_samples), device=dev)
    design = torch.randn((n_samples, width), device=dev)
    torch.cuda.synchronize(dev)
    started = time.perf_counter()
    genotype @ design
    torch.cuda.synchronize(dev)
    first = time.perf_counter() - started
    started = time.perf_counter()
    genotype @ design
    torch.cuda.synchronize(dev)
    steady = time.perf_counter() - started
    out["cublas_warmup_seconds"] = max(first - steady, 0.0)

    # Reader pool: thread creation plus the first scheduling round.
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=reader_workers) as pool:
        list(pool.map(lambda _: None, range(reader_workers)))
    out["thread_pool_seconds"] = time.perf_counter() - started

    out["total_seconds"] = sum(v for v in out.values())
    out["shape"] = {"chunk_variants": chunk_variants, "n_samples": n_samples,
                    "design_columns": width, "reader_workers": reader_workers}
    return out


def measure_zstd_decode(sample_bytes: bytes | None = None,
                        payload_bytes: int = 8 << 20,
                        level: int = 3, threads: int = 16,
                        repeat: int = 3) -> dict:
    """Aggregate zstd decompression rate, measured on OUTPUT bytes.

    The rate a compressed store has to be priced at. Quoted against the bytes
    PRODUCED, not the bytes read, because that is what the codec's own figures
    mean and what the pipeline then has to move.

    **This term was missing from the model entirely, and it mattered.** With no
    decode cost `explain_time` predicted the zstd hard-call store at 7.63x
    where the measurement gives 1.84x. A model that sees compression's saving
    but not its price can only ever over-promise.

    Measured across `threads`, because the scan decodes in its reader pool and
    a single-thread figure would under-price it by that factor. `zstandard`
    releases the GIL while decompressing, so the threads are real work.

    Pass real genotype bytes as `sample_bytes` where possible: the rate depends
    on how compressible the content is, and random bytes are the pathological
    case -- they barely compress and so decompress almost trivially.
    """
    from concurrent.futures import ThreadPoolExecutor

    try:
        import zstandard as zstd
    except ImportError:
        return {"bytes_per_second": 0.0, "available": False}

    if sample_bytes is None:
        # Structured, not random: mostly one value with occasional others,
        # which is roughly what packed two-bit genotypes look like to a
        # compressor.
        import random
        rng = random.Random(20260914)
        sample_bytes = bytes(
            rng.choice((0, 0, 0, 0, 0, 0, 0, 1, 2, 85, 170))
            for _ in range(payload_bytes))

    blob = zstd.ZstdCompressor(level=level).compress(sample_bytes)
    produced = len(sample_bytes)

    def decode_once(_index: int) -> int:
        return len(zstd.ZstdDecompressor().decompress(blob))

    seconds = []
    with ThreadPoolExecutor(max_workers=threads) as pool:
        for _ in range(repeat):
            started = time.perf_counter()
            list(pool.map(decode_once, range(threads)))
            seconds.append(time.perf_counter() - started)
    median = statistics.median(seconds) if seconds else 0.0
    return {
        "bytes_per_second": (produced * threads) / median if median else 0.0,
        "threads": threads, "level": level,
        "ratio": produced / len(blob) if blob else 0.0,
        "output_bytes": produced, "available": True,
    }


QUIET_LOAD_PER_CORE = 0.25
QUIET_FOREIGN_CPU_PERCENT = 100.0
# Processes blocked in uninterruptible I/O. A handful is ordinary; a dozen
# means the storage is saturated.
QUIET_BLOCKED_PROCESSES = 4


def is_quiet(load_1min, foreign_cpu_percent, cpu_count,
             blocked_processes=None) -> bool:
    """Whether a rate measured now deserves to be called the machine's.

    Every signal here is needed because each alone misses a real case. Load
    average is a decaying mean over a minute, so a burst saturating the box
    right now has barely moved it; foreign CPU is instantaneous and catches
    that. Conversely a churn of short-lived processes can show low
    instantaneous CPU while load stays high, and that box is not quiet either.

    **`blocked_processes` exists because CPU signals are blind to a saturated
    disk, and this module's headline product is a disk rate.** Observed
    directly: a host reading load 17.04 and 429% foreign CPU -- borderline by
    the CPU measures alone -- with **12 processes in D state** from another
    user's training jobs, enough that a test suite on the same filesystem sat
    at 4.8% CPU waiting on I/O. Calibrating storage there and calling the
    result the machine's capability is precisely the mistake this predicate is
    supposed to prevent.

    An unknown signal cannot veto: `/proc/loadavg` and `ps` are both absent on
    some hosts, and refusing to ever call such a machine quiet would mean never
    calibrating on it at all.
    """
    cores = cpu_count or 1
    load_ok = load_1min is None or load_1min < QUIET_LOAD_PER_CORE * cores
    foreign_ok = (foreign_cpu_percent is None
                  or foreign_cpu_percent < QUIET_FOREIGN_CPU_PERCENT)
    io_ok = (blocked_processes is None
             or blocked_processes <= QUIET_BLOCKED_PROCESSES)
    return bool(load_ok and foreign_ok and io_ok)


def host_contention() -> dict:
    """Load and foreign CPU, so a profile records the conditions it was taken in.

    A calibration taken on a busy machine is worse than none: it silently
    reports the contended rate as the machine's capability and the planner then
    mis-sizes everything against it. Observed directly -- the same host measured
    HBM at 2.52 TB/s quiet and 1.05 TB/s while another user held ~35 of 48
    cores, which moved the computed roofline ridge from 17.8 to 43.1 FLOP/byte
    and would have changed which resource the model says binds.

    This cannot be fixed by averaging or by more repeats, because the
    contention is not noise around a true value -- it is a different machine.
    The only honest options are to wait, or to record the conditions and let
    the consumer judge. Both need this measurement.
    """
    result = {"load_1min": None, "foreign_cpu_percent": None, "cpu_count": os.cpu_count()}
    try:
        with open("/proc/loadavg", encoding="ascii") as handle:
            result["load_1min"] = float(handle.read().split()[0])
    except OSError:
        pass
    try:
        import subprocess
        me = subprocess.run(["id", "-un"], capture_output=True, text=True).stdout.strip()
        ps = subprocess.run(["ps", "-eo", "pcpu=,user="], capture_output=True, text=True).stdout
        foreign = sum(float(line.split()[0]) for line in ps.splitlines()
                      if line.split() and line.split()[-1] != me)
        result["foreign_cpu_percent"] = foreign
    except Exception:  # noqa: BLE001 - a missing `ps` must not fail calibration
        pass
    # Processes blocked in uninterruptible I/O, any user. The CPU signals above
    # are blind to a saturated disk, and a disk rate is this module's headline
    # product -- see `is_quiet` for the host that read borderline on CPU while
    # twelve of another user's processes had the filesystem pinned.
    try:
        import subprocess
        states = subprocess.run(["ps", "-eo", "stat="],
                                capture_output=True, text=True).stdout
        result["blocked_processes"] = sum(
            1 for line in states.splitlines() if line.strip().startswith("D"))
    except Exception:  # noqa: BLE001 - a missing `ps` must not fail calibration
        result["blocked_processes"] = None
    result["quiet"] = is_quiet(result["load_1min"],
                               result["foreign_cpu_percent"],
                               result["cpu_count"],
                               result.get("blocked_processes"))
    return result


def calibrate(sample_file: str | Path, scratch_dir: str | Path,
              n_samples: int, device: int = 0,
              metadata_file: str | Path | None = None,
              gemm_shapes: tuple[int, ...] = (1, 8, 32, 128, 512)) -> dict:
    """Measure everything the planner needs, and return a profile.

    `gemm_shapes` is a LADDER of trait counts rather than one value, because the
    kernel crosses the device's roofline ridge somewhere inside that range and
    the crossover is where it is, not where a default guessed. Recording the
    whole ladder lets the planner interpolate to the K it is actually asked for
    instead of assuming one rate holds everywhere.
    """
    import torch

    profile: dict = {
        "torch": torch.__version__,
        "device_name": torch.cuda.get_device_name(device),
        "device_capability": list(torch.cuda.get_device_capability(device)),
        "device_memory_bytes": torch.cuda.get_device_properties(device).total_memory,
        "cpu_count": os.cpu_count(),
        "n_samples": n_samples,
    }
    profile["contention_before"] = host_contention()
    profile["disk"] = measure_disk_read(sample_file)
    if metadata_file:
        profile["text_parse"] = measure_text_parse(metadata_file)
    profile["write"] = measure_write(scratch_dir)
    profile["transfers"] = measure_transfers(device)
    profile["startup"] = measure_scan_startup(n_samples, device=device)
    profile["pinning"] = measure_pinning()
    profile["gemm"] = {
        str(k): measure_gemm(n_samples, k, 27, device=device)
        for k in gemm_shapes
    }
    ridge = (profile["transfers"]["hbm_bytes_per_second"] or 1.0)
    best = max(v["gemm_flops_per_second"] for v in profile["gemm"].values())
    profile["ridge_flops_per_byte"] = best / ridge
    profile["contention_after"] = host_contention()
    # A profile that knows it was taken on a busy box is usable -- one that does
    # not is a trap, because every rate in it reads as a machine capability.
    profile["trustworthy"] = bool(profile["contention_before"]["quiet"]
                                  and profile["contention_after"]["quiet"])
    return profile


def write_profile(profile: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(profile, indent=1))
    return path
