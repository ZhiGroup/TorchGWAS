# Autotune redesign: data movement, online model, planner — 2026-09-24

Status: design for review. It replaces the "empirical segments" tuner
described in [empirical autotune](empirical_autotune_20260923.md) and keeps
that document's measurements as evidence.

Implemented since the first draft (2026-09-24):
- `topk_per_trait` removed (API, CLI, writer, tests); `p_value_threshold`
  stays.
- Significant pairs on variant shards (`variant_devices` with
  `reduce='significant'`). Pairs are selected on each shard's own thread.
- Indexed output published in (variant, phenotype) order (section 1.1).
- Chunk candidates checked against the memory model before the run.

FP32 only: TF32 was tested by the user and rejected.

## 0. Decisions proposed

1. Treat the scan as the blocked product `C = Gᵀ Y` (G: N × M genotypes,
   Y: N × K residualized phenotypes) and choose the partition by bytes moved
   and work repeated, priced by the calculator — not by fixed tile rules.
2. When the phenotype panel fits on every GPU, **split variants, not
   phenotypes**, for significant pairs too (today only full and JAGWAS output
   can). Each variant is decoded, transferred and centred once, on one GPU.
3. When it does not fit (voxel scale), keep phenotype tiles resident and
   **decode the genotype once**: cache it in its compact transfer form (GPU
   memory if it fits beside a tile, else pinned host memory) and replay later
   passes from the cache. Streaming phenotypes past a resident genotype block
   moves 32x more bytes at voxel scale (section 2).
4. Replace segment-level A/B trials with a per-chunk, per-stage model: the
   calculator supplies per-variant costs, the run's own chunks supply
   per-chunk overhead and contention, and the tuner re-plans whenever
   measurement and prediction diverge. No forward/reverse ordering, no fixed
   trial schedule, and every candidate's memory is known before it is tried.
5. Remove `p_value_threshold` and `topk_per_trait` (section 1). That removes
   the only single-writer special case.

## 1. Output path and filters

Multi-GPU output already goes through queues: trait tiles each write their
own store from a per-GPU writer thread, variant shards write shard stores, and
multi-GPU JAGWAS uses a bounded shared queue; manifests join them. Nothing
requires output to be on one GPU.

The exception is the row-filter writer used by two options that are still in
the code on `main` and `jagwas-dev` (no branch removes them):

- `topk_per_trait` (`--topk-per-trait`): keep the k smallest p per phenotype.
- `p_value_threshold` (`--p-value-threshold`): compute full beta/t, then keep
  (variant, phenotype) rows with p <= threshold in an indexed store. This is
  a host-side duplicate of `reduce='significant'`, which selects pairs on the
  GPU against a critical |t| and already supports multi-GPU tiles.

Decision: `topk_per_trait` was stale and is removed. `p_value_threshold` stays
(it is the significant-pair output for full-output runs). Autotune still keeps
it on one GPU; routing it through `reduce='significant'` would lift that.

### 1.1 Output order

Indexed pair stores (significant pairs, p-value filter) and JAGWAS stores are
published in (variant, phenotype) order. Readers iterating parts in manifest
order see sorted rows. The internal per-variant reduction keeps its defined
order (descending |t| within a variant).

- **Parts are producer runs, not chunks.** The writer used to write one part
  file, and fsync it, per chunk: 196 fsyncs per GPU pass on the benchmark,
  about 6 s on lab-2080ti's output disk, all on the one writer thread. Now,
  when no JIT progress observer needs per-chunk durability, chunks are
  coalesced into parts of up to 2^18 rows. Each chunk joins the "lane" whose
  range ends where it starts, so each variant shard or phenotype tile keeps
  its own contiguous run without the writer knowing the producers. This is
  opt-in (`coalesce_rows`), set by `run_linear_gwas` whenever the JIT path
  is not observing per-chunk writes. The writer's default, and every
  observed run, stays one part per chunk, which is what the calculator's
  write model (`reduced_output_work`) describes and tests against.
- **Within each part.** Rows are sorted by (variant, trait) as the part is
  written. Device selection walks blocks of variants × trait strips, so its
  raw order is not sorted. The cost is one `lexsort` per part.
- **Across parts, disjoint variant ranges** (one GPU, variant shards,
  JAGWAS). The manifest lists parts by variant range; no data moves.
- **Across parts, overlapping ranges** (phenotype tiles, and device selection's
  trait strips over one variant block). Only in coalescing mode; per-chunk
  mode keeps parts as written, and the manifest then makes no global-order
  claim. Before
  the manifest is written they are merged in 65,536-variant windows (memory
  bounded by the parts that intersect a window) and rewritten once as large
  `ordered_*.npz` parts. The first version wrote one part per window with
  fsync and cost 6 s for 16k pairs; it now writes a part per 2^18 rows.
  `sumstats_write.ordering` records the time.

A cheaper alternative, if the merge ever matters: when all tiles run at once
with shared decode, they deliver the same chunk within a ring's depth of each
other. The writer could merge per chunk as they arrive, with no second pass.

## 2. Data movement: which operand stays, which one moves

Sizes (float32 phenotypes; genotype in its transfer form):

| Workload | N | K | M | Phenotypes | Genotype, 2-bit | Genotype, int8/uint8 |
|---|---:|---:|---:|---:|---:|---:|
| Frozen H100 panel | 35,365 | 16,385 | 1,048,576 | 2.3 GB | 9.3 GB | 37 GB |
| Voxel (pipeline_model example) | 33,417 | 2,085,000 | 1,048,576 (assumed) | 279 GB | 8.8 GB | 35 GB |

From `pipeline_model.auto_trait_block` (80 GB GPU, chunk 2,048, depth 4;
`benchmarks/voxel_memory_probe_20260924.py`): one GPU holds about 510,000
voxels, and the phenotype block is 68 of the 73 GB the model charges. The host
result ring (17-34 GB per GPU) is not the binding limit at 500 GB host.

Three ways to organize the product:

- **Phenotype-stationary** (today): each GPU holds a phenotype tile; the
  genotype streams past it. Genotype passes P = ceil(K / (tile capacity x
  GPUs)). Every pass today re-reads and re-decodes the genotype.
- **Genotype-stationary** (the phenotype-streaming idea): each GPU holds a
  genotype block, and phenotype tiles stream past it. Every GPU must see all K
  phenotypes for its variants.
- **2-D blocks** (SUMMA-like): a grid of variant × phenotype blocks. The two
  above are its corners.

Bytes over PCIe per GPU, voxel case:

| GPUs | Phenotype-stationary, re-decode | Phenotype-stationary, cached 2-bit genotype | Genotype-stationary |
|---:|---|---|---|
| 8 | 1 pass: 8.8 GB | 8.8 GB (no replay needed) | 279 GB |
| 2 | 3 passes: 26 GB, 3 decodes | 26 GB from host cache, 1 decode (8.8 GB if cached on GPU) | 279 GB |
| 1 | 5 passes: 44 GB, 5 decodes | 44 GB from host cache, 1 decode (8.8 GB if cached on GPU) | 279 GB |

Streaming phenotypes moves 32x more than one genotype pass here, because the
compact genotype (0.25-1 byte per call) is 4-16x smaller per element than a
float32 phenotype and K > M. In general: with the stationary operand filling
the fast memory, both corners move about the same bytes once several rounds
are needed. When the genotype fits in aggregate GPU memory,
phenotype-stationary moves at most one genotype pass more than
genotype-stationary (P·genotype ≤ panel + genotype) and usually far less. It
also keeps residualization and per-trait output contiguous. The planner should still evaluate the 2-D grid, because the
answer depends on N, K, M, encoding and GPU memory. For example, a huge
imputed genotype with a small panel gives P = 1 with phenotypes stationary,
and genotype-stationary is infeasible.

What the idea gets right is that **decode should happen once**. The
mechanism that achieves it cheaply is a compact genotype cache, not phenotype
streaming:

- Pass 1 decodes as today and tees each chunk's transfer buffer (2-bit or
  int8) into the cache: GPU memory when the cache fits beside a tile (voxel
  2-bit: 8.8 GB next to a ~440,000-voxel tile), otherwise pinned host memory,
  otherwise nothing (fall back to re-decode).
- Passes 2..P replay chunks from the cache: no decode, and no PCIe if the
  cache is on the GPU. From host, a pass costs 0.16 s on lab-h100 (54.7 GB/s)
  and 1.3-2.7 s on lab-a100 (6.7 GB/s, or 3.3 behind a shared uplink).
- With the cache on the GPU, the number of passes stops mattering. Total GEMM
  is fixed at 2·N·M·K; only the memory split changes.
- Inside a GPU, phenotypes can additionally be walked in sub-tiles against the
  same genotype chunk. Sub-tile width then bounds only the per-chunk
  product/result buffers (chunk × width), not residency, and becomes a free
  runtime knob like chunk size.

How much this buys depends on the decode/GEMM ratio, which is the
calculator's job to price:

- Voxel scale in FP32 is GEMM-bound. 2·N·M·K = 1.5e17 FLOP is ~2,900
  GPU-seconds at an assumed 50 TFLOP/s, so extra decode passes are a few
  percent.
- Decode-bound jobs on a loaded host are where it matters. On lab-a100, three
  redundant concurrent decodes cost 9 s of a 27 s job (4 tiles, median 26.9 s
  separate vs 17.8 s shared).

### Variant sharding for significant pairs

When the panel fits on every GPU (the frozen panel: 2.3 GB), phenotype tiles
are the wrong split. Each tile still receives, converts and centres every
genotype chunk, so per-chunk GPU work and H2D are duplicated G times, and
shared decode only removes the CPU part. Variant shards give each GPU M/G
variants against the full panel. The indexed significant-pair store already
uses range-relative variant indices (as JAGWAS shards do), so shard stores can
be merged. `variant_devices` used to be refused for `reduce='significant'`. No
statistical reason existed: selection is per cell. It was a guard around an
unwired path, and it is now implemented
(`tests/test_variant_output.py::test_significant_variant_shards_match_single_gpu`:
2- and 3-GPU shards, with a variant range and missing calls, write exactly the
single-GPU pair set).

A pitfall found while measuring: significant pairs are selected on the host
by default (`TORCHGWAS_SIGNIFICANCE_BACKEND=host`). The full chunk × K t-matrix
is copied to the CPU and filtered with NumPy. The first shard version filtered
after the shared result queue, i.e. on one consumer thread for all shards, and
ran slower than one GPU (H100: 35.2 s vs 14.7 s executor). Trait tiles already
filter on each device's worker thread. Shards now do the same
(`linear_scan_multigpu(shard_transform=...)`).

Host selection is itself a scaling problem. Per chunk it moves chunk × K × 4
bytes to the host and scans chunk × K cells on one core: 33 MB at K = 8,192,
but 2 GB at a 500,000-voxel tile. Device selection returns only the passing
pairs. Section 2.2 measures the two.

### JAGWAS

JAGWAS keeps the full phenotype panel and the K x K factor on every GPU and
shards variants; phenotype partitioning is not allowed for it. When the panel
and factor do not fit one GPU, the planner stops with an error.

For full and significant output (not JAGWAS), the planner's order of preference is:
1. Panel fits one GPU → variant shards over the GPUs that pay (one decode,
   no duplication).
2. Panel fits across the GPUs → one round of phenotype tiles, shared decode,
   fan-out only where the uplink is the bottleneck.
3. Otherwise → P rounds, genotype cached after pass 1.

## 2.2 Measured: layouts, selection backend, where the time goes

`benchmarks/layout_profile_20260924.py`: significant pairs, N=20,000,
M=200,000, K=8,192, threshold 1e-5, chunk 1,024, fresh process per run, 2
repeats. Median executor seconds (setup + scan + write):

| Layout | 2080 Ti device | 2080 Ti host | A100 device | A100 host | H100 device | H100 host |
|---|---:|---:|---:|---:|---:|---:|
| 1 GPU | 13.3 | 22.2 | 15.4 | 21.6 | 8.6 | 82.8 |
| 2 tiles | 10.3 | 12.3 | 10.4 | 17.3 | | |
| 4 tiles | 9.5 | 13.0 | 20.6* | 17.2 | | |
| 2 variant shards | 11.2 | 11.1 | 12.5 | 15.3 | 26.7* | 91.9 |
| 4 variant shards | 10.3 | 12.8 | 19.4 | 18.6 | 8.2 | 154.0 |

(* one noisy run. The H100 row is partial, one run per cell; the A100 had other
users' jobs on GPUs 3, 4 and 6 and a load average near 90.)

- **Device selection is never meaningfully slower, and host selection
  collapses under CPU load** (H100: 10x). Autotuned runs now select on the GPU
  unless `TORCHGWAS_SIGNIFICANCE_BACKEND` is set. The global default stays
  `host` because the JIT calculator prices that path. Recommendation: flip
  it.
- **Two shard bugs were found and fixed while measuring.**
  - Selection ran on the one consumer thread after the shared queue.
  - Each shard received owned chunk × K copies of every result: about 40 s of
    system CPU on the 2080 Ti and 11x slower than one GPU on the H100. Shards
    now select on their own thread and read the scan's ring directly.
- **Variant shards vs phenotype tiles is not a fixed rule.** Each shard
  copies the whole residualized panel to its GPU: setup was 2.1-5.4 s per
  shard vs 0.3-1.0 s per tile, and on the 2080 Ti the copy goes through host
  memory. Tiles instead repeat each chunk's transfer and centring on every
  GPU. At M = 200k the panel copy dominates; at genome scale it amortizes.
  The planner has to price both.
- **Multi-GPU gains are small at this size (1.4x at best)** because most of
  the run is not GPU work. A cProfile of the 1-GPU run (24.2 s):
  - **Phenotype QC:** `_phenotype_column_mask`, 7.0 s on one core. About eight
    full passes with full-size temporaries over the 655 MB panel, called
    twice. That is more than the GPU scan, and it grows with the panel (279 GB
    at voxel scale). The residualization right after it already has the
    panel on the GPU and computes the same column moments there.
  - **Per-chunk synchronization in device selection:** 8,012 `.cpu()` waits,
    7.2 s. The scan submits a chunk and waits for its selection before
    submitting the next, so GPU work does not overlap host submission. This
    is the per-chunk floor below.
  - **Residualization:** 2.6 s. **Genotype metadata:** 0.9 s.

### Stage timings per chunk (`benchmarks/stage_model_probe_20260924.py`)

One GPU on lab-2080ti, chunk size switched every 6 chunks to a random choice
of 256-4,096, every chunk measured.

| K | Statistics (GPU) | H2D | Decode CPU | Period at 1,024 |
|---:|---|---|---|---|
| 8,192 | 31.4 µs/variant, intercept ≈ 0, R² 0.997 | 1.7 µs/variant (12 GB/s), R² 0.95 | 16.5 µs/variant, R² 0.86 | 39 ms |
| 2,048 | 11.3 µs/variant, R² 0.984 | 1.7 µs/variant, R² 0.94 | 21.0 µs/variant, R² 0.96 | 10.7 ms |

- **Per-variant stage costs are linear and predictable.** 31.4 µs/variant at
  K = 8,192 is about 10 TFLOP/s against the calculator's FLOP count, and H2D
  matches the step-1 calibration. These are the "easy" prices.
- **How stages combine is not.** At K = 8,192 the period matches the sum of
  the GPU-side stages; at K = 2,048 it matches the max. A max-of-stages
  prediction from two sizes was 28-40% off.
- **Small chunks show a floor of roughly 10-17 ms** that no stage span
  accounts for: the per-chunk synchronization above.
- **Throughput is not monotone in chunk size** (K = 2,048: 1,024 beat 2,048
  and 4,096).

Consequence for section 3: fit the observed chunk period directly as
a + b·c, taking b from the calculator as the starting guess and a from the
run, rather than composing stage spans. And remove the synchronization floor
first, since it currently decides the chunk size.

### Fixes 1 and 2 (implemented, lab-2080ti, 2 repeats, device selection)

1. **Phenotype QC on the scan GPU** (`preprocess._phenotype_qc`): the same
   keep mask, observed counts and missing-cell count, in ~512 MB column
   blocks, float64 reductions. `tests/test_phenotype_qc_device.py` checks it
   against NumPy: constant, all-missing, single-observation and large-mean
   near-constant columns, and infinities. `TORCHGWAS_PHENOTYPE_QC=host`
   restores the CPU path.
2. **Selection one chunk behind, on its own stream**: chunk i+1 is queued
   before chunk i's pairs are collected. It is used only when one extra
   chunk × K result is under 20% of free GPU memory (`TORCHGWAS_SELECTION_LAG`
   overrides). At voxel scale a chunk's GEMM dwarfs the wait anyway.

| Layout | API s before → after | Executor s before → after |
|---|---|---|
| 1 GPU | 20.8 → 15.3 | 13.3 → 12.4 |
| 4 tiles | 15.1 → 12.2 | 9.5 → 9.8 |
| 4 variant shards | 17.5 → 13.1 | 10.3 → 10.2 |

The QC saving is 3-5 s per run, and at voxel scale it grows with the panel.
The re-profiled 1-GPU run (14.9 s, was 24.2 s):
- QC 1.5 s, mostly reading the memory-mapped phenotype.
- The scan (~9 s) now matches its GPU work: statistics 6.1 s plus the
  selection kernels. The remaining `nonzero`/`.cpu()` time is the host
  waiting on a busy GPU, not the GPU waiting on the host.
- Residualization 1.4 s, PGEN metadata 0.8 s.

Fix 2 therefore changed who waits rather than the executor time. On one GPU
the GPU really is the limit. Multi-GPU runs are still bounded by per-GPU
setup and the host work above.

### Fix 3: per-chunk model tuner (implemented, default)

`model_autotune.ModelChunkTuner` (`autotune_options={'tuner': 'segments'}`
restores the old trials). It has the same interface as EmpiricalChunkTuner.

- **Samples.** A sample is one chunk's completion interval on a device,
  taken only between two same-size chunks after the switch has settled. A
  few samples per GPU per size suffice (default 4), so probing costs a few
  percent of the job rather than a quarter.
- **Order.** Sizes are probed in random order. No drift model is assumed.
- **Load.** Each sample carries the reader's scheduler-wait fraction.
  Rates are compared at the median load through a fitted
  log(period) = α_size + λ·w. Hosts without schedstat give no w, and λ
  stays 0.
- **Decision.** Keep the starting size unless another beats it by the
  margin (3%).
- **Drift.** After committing, a load-adjusted EWMA of the committed size's
  rate is tracked against its estimate. A deviation over 15% for 8 samples
  re-probes the neighbouring sizes, reusing the recent committed samples;
  at most 3 times per job.
- **Budget.** Probing takes at most `trial_fraction` (default 30%) of the
  remaining rows; the largest sizes are dropped first. A job too short to
  repay the probe keeps the start size and is re-checked as it runs.
- **Model.** `model` in the audit is period(c) = a + b·c, fitted on the
  measured sizes. It is reported, not used to pick unmeasured sizes, because
  overlap makes stage composition unreliable (above).

Simulation tests (`tests/test_model_autotune.py`): fastest size decided in
under 15% of the job; random 1x/4x load bursts do not decide the size (six
seeds); small gains keep the start; a mid-job change re-probes and
switches; short jobs and two GPUs. All real-GPU autotune and format tests
pass with it.

Measured (`--mode tuner`: 1 GPU, significant pairs, device selection, the
benchmark shape, tuning forced on with `min_job_seconds=0`, 3 repeats,
executor seconds):

| Host | fixed 512 | fixed 1024 | fixed 2048 | model tuner | segment tuner |
|---|---|---|---|---|---|
| lab-2080ti (idle) | 11.5, 11.9, 12.0 | 11.4, 11.6, 13.2 | 12.3, 12.3, 12.6 | 12.8, 13.8, 14.0 | 13.5, 14.0, 14.0 |
| lab-h100 (CPU-loaded) | 3.2, 3.7, 7.1 | 3.5, 5.6, 7.1 | 3.5, 12.2, 57.2 | 3.6, 4.1, 9.6 | 11.0, 33.1, 50.2 |

- **Under load, the model tuner tracks the best fixed size** (median 4.1 s vs
  3.7 s), and the segment tuner collapses (33 s). Its long trials hit the
  load bursts that make fixed 2048 swing from 3.5 to 57 s.
- **On the idle host neither tuner can win at this size.** The sizes differ
  by under 6%, less than the cost of probing plus the ring allocated for
  2048. The model tuner chose 512, the best measured size, and lost about
  15% to fixed 1024/512 overall. With the default `min_job_seconds=20`,
  runs this short skip probing and keep 1024.
- **One re-probe was triggered by per-chunk noise** (±20% at this job size),
  and it re-confirmed 512. λ stayed 0 on lab-2080ti because it exposes no
  scheduler-wait counter.

**Update 2026-09-25: decisions and drift account for noise.** At full scale
on the loaded A100, per-size estimates of the same chunk size moved 3x
between visits (18,262 then 5,889 rows/s at 1024). The fitted period had a
negative intercept, and the tuner switched sizes on noise: autotune took
279.9 s where fixed chunk 1024 took 218-231 s.

- **Switching:** a size is taken only if its gain over the incumbent
  exceeds both the margin and one standard error of the difference of
  medians. Each median's error is 1.253 · 1.4826 · MAD / √n from its own
  samples. The decision audit records the errors.
- **Drift:** the committed size's first 2·`residual_chunks` samples are the
  reference. Each later non-overlapping window of the same size is compared
  with it by a Mann-Whitney rank test, because per-chunk rates are skewed and
  heavy-tailed. A window strikes when |z| > 3 and the median moved more than
  `residual_threshold`, and two strikes in a row are drift.
  - Under ±60% synthetic chunk noise the old EWMA rule re-probed without any
    change. The rank test does not, and a sustained change still re-probes
    (`tests/test_model_autotune.py`).

### Fix 4: variant shards vs phenotype tiles, priced (implemented)

`layout_pricing`, used by autotune for significant pairs when the panel fits
every GPU and the rules chose two or more GPUs (`autotune_options['split']`:
`auto`, `variants`, `traits`). Following section 2, GPU arithmetic is equal
for both splits, so only data movement is priced:

- **Variant shards:** up-front panel copies from the first GPU,
  Σ N·K·4 / peer bandwidth, then max(GEMM/G, (M/G)·N·b / H2D).
- **Phenotype tiles:** N·K·4/(G·H2D) up front, then max(GEMM/G, M·N·b / H2D).

H2D and peer bandwidth are measured at startup with one 64 MB pinned copy to
all chosen GPUs at once, so a shared uplink or other users' traffic shows up
in the price. The GEMM rate is a prior from SM count and architecture. Shards
are chosen only when priced at least 5% cheaper. On the benchmark shape and
H100-like prices both splits are GEMM-bound and tiles are kept. At A100
shared-uplink bandwidth (3.3 GB/s) with K = 4,096, M = 10⁷, the tiles' stream
exceeds their GEMM and shards win (`tests/test_layout_pricing.py`, plus
real 2-GPU runs matching one GPU). The measured shard setup on lab-2080ti
(2-5 s) is larger than the priced copy (~0.3 s): the rest is per-shard
preprocessing, not yet priced.

**Update 2026-09-25: ties go to shards.** Like-for-like runs (device
selection on both layouts, `--mode selection`, K=8,192, 200k variants,
executor s, 2 repeats):

| | 1 GPU | 2 shards | 2 tiles | 4 shards | 4 tiles |
|---|---|---|---|---|---|
| lab-2080ti (busy host) | 13.4-13.8 | 7.1-7.9 | 9.5-9.7 | 5.7-7.3 | 9.0-12.0 |
| H100 (load ~140) | 3.8-3.9 | 1.8-2.9 | 10.2-13.8 | 2.3-2.5 | 5.1-10.4 |

The price had both within 5% (3.29 s vs 3.37 s on lab-2080ti), and the old
tie rule kept tiles. Autotune then took 2 tiles and ran in 10.8 s on
lab-2080ti and 27 s on the H100.

The model has no term for what separates the layouts. Tiles run in lockstep
off one shared decode, and every tile pipeline handles every chunk. Shards
are independent and handle 1/G of the chunks each.

`choose_split` now returns tiles only when they are priced at least 5%
cheaper. That still happens when the up-front panel copies over a slow peer
path decide, which is the wide-panel case in the tests. The earlier
"tiles by default" rested on runs made before the shard setup fix.

### Variant-shard count: setup against work (implemented)

Shard counts used to be min(GPUs, M // (8·chunk)). With 200k variants that
is 12, so H100 JAGWAS autotune took all 7 idle GPUs and ran in 9.0-14.7 s,
where 2 shards took 2.4 s. Each shard adds serial setup (CUDA context,
cuBLAS handle, pinned rings, thread start), while the work splits. So G
minimizes

    s·G + W/G,   W = M × gpu_seconds_per_variant,   s = shard_setup_seconds

- **s** is timed on the second GPU (context plus a small GEMM), which the
  job uses anyway because at least two are kept.
- **Fitted value:** T(G) = a + s·G + W/G fitted to lab-2080ti JAGWAS
  (1, 2 and 4 shards: 43.3, 24.2, 16.2 s) gives s = 1.0 s and W = 40 s.
  The context alone measured 0.8 s there.
- **Resulting counts:** 7 for lab-2080ti JAGWAS (W ≈ 50 s), 2 for H100
  JAGWAS (W ≈ 1.6 s, s ≈ 0.5 s), 3 for lab-2080ti significant pairs
  (W ≈ 7 s).
- **Overrides:** `autotune_options['shard_setup_seconds']` replaces the
  measurement, and 0 restores the variant-count rule.
- **Audit:** `layout['shard_model']` records the inputs.

### Variant-shard setup (fixed)

py-spy over all threads (`benchmarks/pyspy_summary_20260924.py`) put each
shard's setup in the column-blocked panel upload (`native_scan.py:148-149`).
The panel was residualized on the first GPU and downloaded to host. Every
shard then re-uploaded it from pageable memory, after a CPU copy of each
strided column block, with four shards at once. Now the first GPU keeps the
residualized panel (`residualize_and_standardize(keep_on_device=True)`), and
each shard copies its design device to device. This is used only when a second
copy fits on that GPU; otherwise the host path is kept.

| lab-2080ti, device selection, 2 repeats | Setup per shard (s) | Executor (s) |
|---|---:|---:|
| 4 shards, before | 4.3 | 10.2 |
| 4 shards, after | 2.0 | 8.0 |
| 2 shards, after | 0.9 | **7.5** |
| 2 tiles / 4 tiles | 1.1 / 0.3 | 9.3 / 8.7 |
| 1 GPU | 1.2 | 12.4 |

Variant shards are now the fastest layout at this shape on lab-2080ti (2
shards, 19% under 2 tiles). The remaining per-shard setup (context, pinned
rings, design) is not priced yet.

### Host genotype cache for tile rounds (opt-in)

`genotype_cache.CachedFillSource` wraps a native-fill source. The first time
a variant range is decoded, its transfer-form rows are copied into a host
array, and later fills of those rows are a memory copy. It is enabled by
`TORCHGWAS_GENOTYPE_CACHE=1` for tile runs with more tiles than GPUs (rounds),
and only when the cache fits in half of available host memory. Outputs are
identical and later rounds hit the cache (`tests/test_genotype_cache.py`).
Measured (`--mode cache`: 4 rounds, significant pairs, executor s):

| Host | 1 GPU off | 1 GPU on | 2 GPUs off | 2 GPUs on |
|---|---|---|---|---|
| lab-2080ti | 15.4, 15.7, 16.7 | 15.9, 17.7 | 14.7, 16.2 | 15.7, 18.0, 18.0 |
| lab-h100 (1 run) | 23.9 | 23.5 | 28.4 | 55.4 |

No gain. Hardcall PGEN decode is cheap (scan-profile decode waits near
zero), so rounds are GPU-bound. The cache adds its fills (4 GB of fresh host
memory) and a 4 GB memory copy per later round. It stays opt-in, for
expensive decoders (zlib BGEN on the CPU) or a starved CPU. That section 2
argument (decode once) holds only when decode is the bottleneck, and here it
was not.

### JAGWAS stays on the full panel (split removed)

A split-panel JAGWAS (phenotype tiles on several GPUs, t joined per chunk)
was built and then removed. The rule is that JAGWAS runs only when the whole
phenotype panel and its joint-test state fit on every GPU; multi-GPU JAGWAS
shards variants. The planner now enforces that up front: the one-GPU check
counts the K x K float64 factor (3·K²·8 bytes while preparing) as well as the
panel, and stops with a clear error when they do not fit.

### What-if: the JAGWAS join as its own pipeline stage (experiment only)

The question: if the panel were split over G GPUs, could the per-chunk join
of t run as its own stage, overlapped with the next chunk, and cost nothing?
This is not a product path; JAGWAS stays full-panel. The experiment,
`benchmarks/jagwas_join_stage_bench_20260925.py`, is synthetic and
decode-free, so only the GPU side is measured:

- **vshards:** each GPU takes every G-th chunk with the full panel. It runs
  the scoring GEMM and then the block-triangular FP64 projection inline, as
  the product does.
- **split:** GPU b holds traits [s_b, e_b), and every GPU scores every chunk.
  GPU b pulls t_g for g < b from its peers and computes
  Σ_g t_g L⁻¹[b, g]ᵀ (its diagonal block is also sub-blocked) and ||·||².
  - `split_inline`: the join runs on the scoring stream.
  - `split_stream`: the join runs on a second stream, issued one chunk
    behind, with events for the cross-GPU dependencies and for buffer reuse.
- **Tiles:**
  - `equal`: the same width on every GPU.
  - `triangle`: an equal area of L⁻¹ per GPU.
  - `balanced`: equal a·K_b + p·area_b, with a and p measured on GPU 0.

Throughput relative to vshards (1.00), using the best tile choice for each split
row. Settings: 16,384 samples, chunk 2048; H100 384 chunks × 5 repeats, 2080 Ti 48 × 3:

| | vshards, projection on 2nd stream | split_inline | split_stream | split, scoring only vs vshards scoring only |
|---|---|---|---|---|
| H100, 2 GPUs, K=8192 | 1.07 | 0.87 | 0.96 | 0.91 |
| H100, 2 GPUs, K=16384 | 1.03 | 0.85 | 0.88 | 0.94 |
| H100, 4 GPUs, K=8192 | 1.06 | 0.72 | 0.87 | 0.48 |
| H100, 4 GPUs, K=16384 | 1.04 | 0.72 | 0.77 | 0.76 |
| 2080 Ti, 2 GPUs, K=4096 | 1.04 | 0.83 | 0.89 | 0.70 |
| 2080 Ti, 2 GPUs, K=8192 | 1.01 | 0.83 | 0.86 | 0.73 |
| 2080 Ti, 4 GPUs, K=4096 | 1.05 | 0.66 | 0.67 | 0.31 |
| 2080 Ti, 4 GPUs, K=8192 | 1.01 | 0.74 | 0.80 | 0.35 |

What this says:

1. **The stage does overlap.** A second stream recovers 1-20% over the
   inline join, and more with more GPUs, because the join's waits on peers
   stop blocking the next chunk's scoring.
2. **It still does not break even.** Overlap hides latency, not work or
   bandwidth, and the split adds both:
   - Every GPU receives and converts every genotype chunk, so the
     scoring-only split falls to 0.31-0.94 of vshards, worse with more GPUs.
   - GPU b pulls s_b·chunk scores per chunk. On the 2080 Ti (no P2P) this
     goes through host memory.
   - GEMMs are narrower.
   - No tile choice balances both the GEMM (∝ K_b) and the triangle
     (∝ e_b² − s_b²). On H100 `triangle` tiles are the worst choice; on the
     2080 Ti `equal` tiles are.
3. **The loss grows with GPU count** (4 GPUs: 0.67-0.87). The projection is
   the same total work in both layouts, so the split has nothing to trade
   against these costs while the panel fits.
4. **A second stream for the projection inside vshards gains only 1-7%**
   (3 of the 8 rows gain 5-7%). Scoring and projection compete for the same SMs.
   The gain is too small to add a stream and events to the product.

So when the panel fits, variant shards dominate, which matches the rule.

### JAGWAS projection: skip the zero half of L⁻¹ (implemented)

`T = ||L⁻¹ t||²`, and `L⁻¹` is lower triangular, but `JagwasReduction.reduce`
multiplies the full K x K matrix, so half of its FP64 work multiplies zeros.
`jagwas_projection.TriangularJagwasReduction` cuts `L⁻¹` into row blocks
(about K/512 blocks, at most 16). Block b needs only columns up to its own
end, which leaves (B+1)/(2B) of the dense work (0.53 at 16 blocks). The blocks
are views of the one factor, so memory is unchanged. `reduce.py` is
byte-identical (the calculator binds its source).

- **Default (2026-09-25, second pass):** the triangular projection is the
  JAGWAS projection for every run, autotuned or not, and the calculator
  prices it. `TORCHGWAS_JAGWAS_PROJECTION=dense` selects the dense reference
  (`reduce.JagwasReduction`). `sumstats_write.jagwas_projection` records
  which one ran. `reduce.py` stays byte-identical, so the factor calibration
  bound to its sha (`jagwas_preparation`) stays valid; the factor preparation
  is the same code.
- **Calculator changes:**
  - `reduction_tensor_work` traces the class `jagwas_reduction_class()`
    returns. Its cache key carries the projection, and its source binding
    adds `jagwas_projection.py` and `jagwas_blocks.py`. `out=` products are
    counted as GEMMs, and their outputs as writes, not reads.
  - `tensor_service.tensor_stage_service` takes one GEMM triple per product
    and matches each product to its own launches in the census (main kernel,
    plus a separate split-K reduction when present). `gemm` in the result
    sums the products and keeps each one's census under `gemms`.
  - Block b's GEMM is (inner e_b, rows w_b, columns chunk). The captured
    CUTLASS grids confirm the orientation: only it reproduces the grid under
    the identity swizzle (K=512, chunk 128: grid [64, 1, 1], swizzle 4).
  - A one-variant chunk makes each block a GEMV with the vector on the other
    side from dense. It is the same cuBLAS kernel, and `gemm_work` now
    accepts either orientation.
  - New host API prices: `joint_gemm_out_fp64`, `joint_empty_fp64`,
    `joint_slice_view`, `joint_square_inplace_fp64`,
    `joint_sum_fp64_columns` and `joint_unsqueeze_view_fp64`
    (`benchmarks/direct_jagwas_host_primitives.py`).
  - `layout_compute_floor` counts the triangular FP64 FLOPs through the
    torch-free `jagwas_blocks`.
  - The kernel census fixtures were recaptured on A100 with the triangular
    projection (`tests/fixtures/jagwas_projection_geometry.json`,
    `jagwas_chunk_geometry.json`). The statistics launches are unchanged.
    The dense captures are kept as `*_dense.json`.
- **Accuracy:** only the FP64 summation order changes. The FP32 chi² was
  bit-identical to dense in every microbenchmark row and in every end-to-end
  store compared (max |Δ| = 0).

Each block is one GEMM, written straight into its own contiguous rows of a
single K x chunk buffer (`L⁻¹[b, :e_b] · t[:, :e_b]ᵀ`). A single square and
a single column sum follow. The first version launched four kernels per
block (matmul, square, sum, add). On a loaded host with a fast-FP64 GPU,
those launches cost more than the FLOPs they saved.

Projection alone (`benchmarks/jagwas_projection_bench_20260925.py`, dense
and triangular interleaved, median of 9-41 repeats). Each cell is the
speedup over dense at the block count the rule picks (K/512 blocks, at most
16), for chunk 1024 / chunk 4096:

| GPU | K=2048 (4 blocks) | K=4096 (8) | K=8192 (16) | K=16384 (16) | K=32768 (16) |
|---|---|---|---|---|---|
| 2080 Ti (0.5 TF FP64), wall | 1.69 / 1.40 | 1.72 / 1.66 | 1.79 / 1.77 | 1.82 / 1.78 | — |
| H100 (32-62 TF FP64), wall | 1.06 / 1.25 | 1.26 / 1.47 | 1.45 / 1.65 | 1.67 / 1.70 | 1.75 / 1.79 |
| A100 (9-19 TF FP64, busy host), GPU time | 1.06 / 1.03 | 0.90 / 1.43 | 1.42 / 1.66 | 1.72 / 1.91 | 1.74 / 1.80 |

- **Fused vs four-kernel version:** these cells use the fused version.
  - Where launches dominated, the four-kernel version was much worse: A100
    at K=4096 and chunk 1024 measured 0.60, and H100 at K=4096 measured 1.14.
  - FLOP-bound cells did not change.
- **Remaining loss:** one cell loses, A100 at K=4096 and chunk 1024. There
  the dense product takes 2.8 ms, next to an 8.6 ms scoring GEMM, so the
  loss is about 2% of the chunk. The best count differs by up to ±1 power
  of two per GPU and chunk (2080 Ti prefers 32, A100 at small K 4). Picking
  it per device would gain at most about 10% of the projection, and on the
  GPUs where the projection is a large share of the chunk, the rule is
  already within 3% of the best count.
- **Earlier timing artifact:** the first, non-interleaved A100 run showed
  B=1 (the dense kernel) at 0.26x. That was drift on the loaded host, not
  the kernel.

End to end: 200,000 variants x 8,192 traits x 20,000 samples, hardcall
PGEN, chunk 1024, final code, executor seconds per repeat. On lab-2080ti,
another user's job held about 35 of the 48 cores throughout.

| layout | lab-2080ti dense | lab-2080ti triangular | H100 dense | H100 triangular |
|---|---|---|---|---|
| 1 GPU | 70.2, 69.9, 68.1 | 43.3, 43.9, 44.3 | 4.5, 4.0, 3.9 | 4.0, 4.1, 3.6 |
| 2 variant shards | 36.5 | 24.2 | 2.1, 2.8, 2.3 | 2.5, 2.4, 2.6 |
| 4 variant shards | 23.6, 23.1, 23.3 | 16.2, 15.6, 17.6 | — | — |
| autotune (triangular) | | 16.0, 16.4 (4-7 shards) | | 1.8, 2.5, 4.4 (2-3 shards) |

- **lab-2080ti:** the projection is about 54 of the 70 s on one GPU
  (2·M·K² FLOPs at 0.5 TF), and the triangular path removes close to half of
  that (-37% end to end). With 4 shards the rest is setup and decode, so the
  gain is -30%.
- **H100:** FP64 is fast (≈0.5 s of projection per GPU). The difference is
  within the noise of that host (load ~140), which also produced one-off 11
  and 45 s runs.
- **A100 (one GPU, busy host):** dense 10.3 / 14.0 s, triangular
  10.4 / 11.9 s.

### JAGWAS factor: FP64 Gram and a rank cutoff for collinear panels (implemented, 2026-09-26)

The default reduction (`TriangularJagwasReduction`) now forms R with an FP64
Gram and factors it with a rank-revealing (pivoted) Cholesky. The dense
reference keeps the legacy FP32 Gram and strict full-rank Cholesky.

**What R is.** The t-statistics come from the scanned panel: residualised
and standardised in FP32. So their null correlation is exactly that panel's
Gram matrix, and R is computed as the FP64 Gram of those FP32 values, in
sample blocks (`jagwas_blocks.gram_rows`: max(2K, 4096) rows). Cost is
2·N·K² FP64 FLOPs once per device: K=2048, N=35k is ≈17 ms on A100, ≈6 ms
on H100 and ≈0.7 s on lab-2080ti.

**Cutoff.** tol = K·ε·max diag(R), with ε the precision of the scanned panel
and its t-statistics (FP32: 2.4e-4 at K=2048, 1e-3 at K=8192). This is the
LAPACK dpstrf rule, using the data's ε rather than the factorisation's. The
pivot of trait j after the kept traits is 1 − R²_j, so a trait is dropped
when VIF > 1/tol.
- **Rank and df:** T is the quadratic form over the r kept traits, on r df.
  The manifest `df` is r. `manifest['jagwas_rank']` and
  `run_metadata['sumstats_write']['jagwas_rank']` record the method, rank,
  tolerance, smallest pivot, and each dropped trait with its residual
  variance and VIF.
- **Override:** `TORCHGWAS_JAGWAS_PIVOT_TOLERANCE` sets the relative
  tolerance.

**Execution.**
- **Full-rank fast path:** plain Cholesky, L⁻¹, then one check. Every greedy
  pivot is ≥ min_j 1/(R⁻¹)_jj = 1/max_j‖L⁻¹[:, j]‖², so if that exceeds tol
  nothing would be dropped and L⁻¹ is used as is. The check is one column-norm
  reduction and one read of three FP64 scalars (info, max diag, max norm):
  ≈0.25 ms host time on A100.
- **Otherwise:** host LAPACK `dpstrf` picks the kept set, and each device
  factors R[kept, kept] in input order. `reduce` gathers the kept columns
  before the projection.
- **Multi-GPU:** variant shards share one `JagwasRankSelection` (the factory
  is `jagwas.spawn`). The first decision wins, and a device whose own R would
  decide differently adopts it.

**Measured on real panels** (`benchmarks/jagwas_rank_panels_20260925.py`,
`benchmarks/jagwas_rank_error_20260925.py`). These are the 22 panels in
txia2's `batch_seven_plus_torchgwas/phenos` (35,298 samples; the `nceq_*`,
CNN and PCA panels at K=128, the raw `fourier_PE_*` and `graphunet` at
K=86-114). Nothing is dropped on the 17 well-conditioned panels (smallest
1 − R² between 0.0036 and 0.86). The five raw `fourier_PE_*` / `graphunet`
panels have cond(R) 5e7-1e9, one with a negative FP32 eigenvalue.

| panel | K | kept | legacy FP32 Gram, no cut: max rel. T error | legacy FP32 Gram, default cut | FP64 Gram, default cut | FP64 Gram, no cut |
|---|---|---|---|---|---|---|
| fourier_PE_L4_xyz | 90 | 87 | 3.2e-2 | 1.2e-2 | 5.1e-6 | (not PD) |
| fourier_PE_L4_xyz_thick_curv | 86 | 83 | 7.8e-3 | 1.0e-2 | 7.0e-6 | (not PD) |
| fourier_PE_L5_xyz_thick_curv | 92 | 90 | 2.5e-2 | 4.9e-3 | 4.9e-6 | 1.4e-5 |
| fourier_PE_L5_xyz_thick | 97 | 94 | 1.1e-1 | 4.7e-3 | 4.5e-6 | 1.6e-5 |
| graphunet | 114 | 108 | 3.7e-1 | 7.7e-3 | 6.5e-6 | 2.0e-5 |
| CNN (control) | 128 | 128 | 1.6e-5 | 1.6e-5 | 2.2e-7 | 2.2e-7 |

How the errors were measured:
- The error is the maximum over 20,000 null z ~ N(0, R) against the exact R
  of the scanned panel, with z rounded to FP32 as the scan's are.
- The FP32 Gram entries are off by up to 2.4e-7 (≈2ε32). That error is
  almost all Gram rounding: the FP32 residualisation adds nothing
  measurable.
- With the legacy Gram, the maximum error falls as about ε32/τ as the
  tolerance τ rises. It was still 1e-3 at τ = 1e-4, which drops 6-13
  traits. The FP64 Gram removes that source, and what remains is the FP32
  rounding of z.
- The near-duplicate trait in the unit test (true 1 − R² = 1e-8) reads as
  9.8e-9 from the FP64 Gram and 3.6e-7 from the legacy one.

**The cutoff is a heuristic, and t makes it insufficient (measured 2026-09-26).**
K·ε is LAPACK's numerical-rank convention, not a JAGWAS derivation. What a
small pivot λ actually amplifies in T:
- **Rounding of R:** gone with the FP64 Gram.
- **FP32 rounding of z:** ≈ 2·ε_z·√T_v/√λ.
- **The nonlinearity of t in r:** t_j = r_j·√(df/(1−r_j²)). A linear
  dependence among traits (r_c = a·r₁ + b·r₂) holds for r but not for t, so
  a strongly associated variant gives the collinear direction v'·t of about
  z³/2N, and T gains about (z³/2N)²/λ.

`benchmarks/jagwas_nonlinearity_20260926.py`
(`results/jagwas_nonlinearity_20260926.jsonl`) tests this. It injects
genotype effects (strongest single trait at |z| = 0-30, 200 variants per
cell, with effects along a random direction, on a most-collinear trait, or on
a random trait) into the five near-singular real panels and the CNN control.
It compares, on the kept set S the production rank selection picks, against
the exact score form (FP64 √df·r, FP64 R_SS):
- the production statistic (FP32 t through the reduction);
- the FP32 score form √df·r through the same reduction.

Maximum excess T_t − T_ref, worst panel and direction (T_ref median ≈ 90,
119, 197, 491, 994 at z = 0, 5, 10, 20, 30):

| tolerance | z=0 | z=5 | z=10 | z=20 | z=30 |
|---|---|---|---|---|---|
| 0 | 3e4 | 1.3e6 | 4.8e7 | 1.4e9 | 1.4e10 |
| 1e-6 | 1.4 | 32 | 491 | 2.1e4 | 1.9e5 |
| default (K·ε32 ≈ 1.1-1.4e-5) | 0.73 | 5.7 | 69 | 2.3e3 | 2.5e4 |
| 1e-4 | 0.26 | 2.2 | 12 | 285 | 3.3e3 |
| 1e-3 | 0.11 | 0.71 | 3.9 | 53 | 374 |
| 1e-2 (half the traits dropped) | 0.04 | 0.39 | 1.7 | 18 | 88 |
| CNN control (no collinearity) | 0.06 | 0.48 | 2.4 | 31 | 150 |

- **t form:** no tolerance makes strong hits accurate. The control row is
  the ordinary Wald-versus-score gap, which exists on any panel.
- **Score form:** the maximum |T_z − T_ref| over every panel, direction and
  z is 0.017 at the default cutoff, 0.05 at 1e-6 and 0.14 at 1e-7. It is 4-7
  only at tolerance 0, on the singular panels. It is rounding-limited and
  does not grow with effect size.
- **Signal lost by dropping traits** (median T_all − T_S net of the dropped
  df): indistinguishable from 0 at every tolerance ≤ 1e-3. It becomes
  0.4-4 at 1e-2.
- **Score form as a transform:** z = √df·r = t/√(1 + t²/df), element-wise
  from t and the per-variant df. It matches t in the null to O(r²). Unlike t,
  its quadratic form is invariant to adding a redundant trait.

**Implemented (2026-09-26): the default statistic is the score form.**
- **Where:** `TriangularJagwasReduction.reduce` applies z = t·(1 + t²/df)^(−1/2)
  in the scan precision, before the FP64 cast. The FP32 finiteness guards
  now run on t too. It is five element-wise FP32 ops per chunk.
- **Unchanged:** the dense reference (`TORCHGWAS_JAGWAS_PROJECTION=dense`)
  keeps the legacy t form, FP32 Gram and strict factor. The tolerance stays
  K·ε32, which measured ≤ 0.017 error in T up to z=30 and no signal lost.
- **Calculator:** typed host-API keys for the new ops (square, div_, add_,
  rsqrt_, mul_ and the FP32 guards), with matching entries in
  `benchmarks/direct_jagwas_host_primitives.py`. The reduce trace's df
  follows the compute dtype, as the scan's does. The projection and chunk
  kernel census (22 kernels, from 17) and the factor fixtures were
  recaptured on A100.
- **Tests:** added a test that T is unchanged by a redundant trait under a
  z≈20 effect, where the legacy t form changes by more than 0.1%.
  API runs are now checked against independent FP64 OLS references: the
  score form for the default, the t form for dense.

**The cutoff is an accuracy target, and there is one JAGWAS reduction (2026-09-26).**
- **Why not K·ε:** with R exact (the FP64 Gram) and the score form, a small
  pivot amplifies only the rounding of z. K·ε, the numpy/LAPACK
  numerical-rank default, has no relation to that error and scales with K,
  which that error does not.
- **The rule:**
  - First-order error propagation gives the null rms rounding error of T over
    a kept set S as 2·ε_z·√tr(R_S⁻¹), where tr R_S⁻¹ = Σ VIF.
  - ε_z = u·√N is the random-rounding (Higham √n·u) bound on each z, with u
    the scan dtype's unit roundoff.
  - The kept set is the longest prefix of the greedy pivoted-Cholesky order
    whose error is ≤ ΔT\*. That is `TORCHGWAS_JAGWAS_T_ROUNDING`, default
    0.01, about 0.002 in −log₁₀p.
  - tr R_S⁻¹ does not depend on order. So the fast path, which checks
    ‖L⁻¹‖_F on the unpivoted factor with one read of (info, norm), and the
    pivoted selection agree exactly.
  - The pivoted path uses host dpstrf at its FP64 numerical-rank default,
    then cumulative row sums of the leading block of L⁻¹ give every prefix's
    trace.
- **ε_z measured** (`benchmarks/jagwas_score_precision_20260926.py`): the
  null rms rounding of the scan's FP32 z, end to end (cuBLAS GEMM, t
  formula, score transform), grows as √N as the bound says.

  | GPU | N=8,000 | N=35,365 |
  |---|---|---|
  | H100 | 5.8e-7 | 1.2e-6 |
  | A100 | 1.1e-6 | 1.7e-6 |
  | 2080 Ti | 1.1e-6 | 2.4e-6 |

  That is 0.1-0.2 × u·√N. The bound is used as is, with no fitted constant,
  and is 5-10× conservative.
- **On the 22 real panels** (`benchmarks/jagwas_rounding_cutoff_panels_20260926.py`,
  N=35,298, ε_z = 1.1e-5):
  - The 17 well-conditioned panels keep every trait (tr R⁻¹ 129-1.5e4,
    estimated error 0.0003-0.003).
  - The five near-singular panels keep 78/76/81/87/99 of 90/86/92/97/114
    traits (K·ε32 kept 87/83/90/94/108), each at an estimated error just
    under 0.01.
  - The extra traits dropped sit where the injection experiment found no
    measurable signal loss (τ between 1e-4 and 1e-3).
- **FP64 numerical-rank floor (added the same day):** the kept set must also
  satisfy tr(R_S⁻¹) ≤ 1/(K·ε64), R's FP64 numerical rank.
  - Without it, FP64 statistics (rounding limit about 1e25) kept an exactly
    duplicated trait: T stayed accurate, but df counted a direction that
    carries no chi-square.
  - With FP32 statistics the rounding limit (about 2e5) is the binding one,
    so production is unchanged.
  - Every greedy pivot is at least 1/tr, so the fast path stays exactly
    consistent with dpstrf's default rank.
- **Validated through the production reduction:** the injection experiment
  sweeps the target over the five near-singular panels (3,000 variants per
  cell, z = 0-30, `results/jagwas_nonlinearity_D_20260926.jsonl`):

  | ΔT\* | traits kept (of 90/86/97/92/114) | estimated error | measured max \|ΔT\|, z=0 → 30 | signal lost (median, worst) |
  |---|---|---|---|---|
  | ≥ 1 | all (FP64 rank) | ≤ 0.25 | 0.058 → 0.158 | 0 |
  | 0.1 | 89/84/97/92/114 | ≤ 0.091 | 0.017 → 0.063 | 0 |
  | **0.01 (default)** | 78/76/87/81/99 | ≤ 0.0098 | 0.0013 → 0.0062 | ≤ 0.36 |
  | 0.001 | 42-54 | ≤ 0.001 | ≤ 0.0014 | up to 7.4 |

  - The realized error stays within the target at every effect size.
  - The CNN control is at ≤ 0.0016 throughout.
  - A target of 0.1 would keep 10-15 more traits on these panels, at ≤ 0.06
    error and no signal lost.
- **Selector census:** removing the old class changed `reduce.py`, whose hash
  the significant-pair selector census pins. That census
  (`tests/fixtures/device_selection_geometry_cuda0.json`) was recaptured on
  the same A100 GPU. The context and all 12 kernel rows are identical; only
  the source hashes changed.
- **The dense reference is removed:** the legacy t form, FP32 Gram, strict
  factor and `TORCHGWAS_JAGWAS_PROJECTION`. `reduce.JagwasReduction` is the
  single implementation in `jagwas_projection.py`, and it owns the CUDA
  linalg initialisation lock. The calculator lost its dense branches, the
  rank check is (info, ‖L⁻¹‖_F), and the census and factor fixtures were
  recaptured on A100. `benchmarks/jagwas_projection_bench_20260925.py`,
  which compared dense with triangular, is gone. Its results remain in the
  projection section above.

**Computing z from r in the scan (skipping t) was measured, not
implemented.** This is z = √df·gy/√(ss_g·ss_y), instead of the t tail plus
the transform. `benchmarks/jagwas_score_path_bench_20260926.py` gives CUDA-event
stage times per chunk at N=35,365, chunk 2048/4096 and K=512-8192:
- **The GEMM dominates:** the FP32 genotype × panel product is 85-95% of the
  chunk (no TF32), and the FP64 projection is most of the rest.
- **The t tail plus the transform** is 0.2-5 ms of 2-167 ms.
- **The z-from-r path would save** 1.2-2.9% of chunk compute on A100 and
  2.8-4.8% on H100, and less end to end. That is not worth changing every
  statistics backend's return contract and recalibrating the calculator.

**Calculator.**
- The factor ledger traces the blocked Gram: `aten.addmm_` FLOPs are counted,
  and the block casts go to the `cast` phase.
- A new `rank_check` phase is priced: 2K² FP64 multiply-adds, K square roots,
  and a 24-byte read.
- The factor memory floor is the larger of the Gram stage (phenotype, FP64 R,
  one FP64 block) and the solve stage (phenotype + 32K²).
- `FACTOR_PHASES` gains `rank_check`. The fixed-reference bank
  (`benchmarks/direct_jagwas_factor_primitives_20260921.py`) mirrors the new
  phases and binds the ledger's sources. `tests/fixtures/jagwas_factor_phases_cuda{0,2}.json`
  were recaptured on A100.
- The collinear fallback (host dpstrf plus the kept-block factor) is
  data-dependent and listed as unpriced.

### JAGWAS cutoff: eigen truncation (implemented, 2026-09-26, default)

**Decision:** T keeps R's eigen-directions with eigenvalue above 1e-3 of the
largest (numpy `pinv`'s rule), and never more than the 0.01 rounding target
allows (Σ 1/λ over the kept directions plays tr R_S⁻¹). df is the kept count.
- **Settings:** `TORCHGWAS_JAGWAS_RCOND` / `run_linear_gwas(jagwas_rcond=...)`,
  default 1e-3. 0 selects the rounding cutoff over traits from the previous
  section, which `jagwas_min_residual` extends with a VIF threshold.

**Why:** the colleague's collinear imaging panels (fourier_PE_* and graphunet,
86–114 traits, 35k samples) had 662–1,198 loci and hits to P ≈ 1e-1155.
The reference JAGWAS, `pinv(R, rcond=1e-3)` over fastGWA statistics, had
110–240 loci. The evidence, in order:
- **Outliers:** 60–140 samples per panel are extreme in many traits at once.
  They made the low-variance eigen-directions heavy-tailed, with median
  kurtosis 100–1,400 against about 0 for a well-behaved CNN panel.
  - Dropping their rows made those directions Gaussian.
  - Masking or clipping only the extreme values made them heavier-tailed
    still, because it breaks the traits' linear relations for those samples.
- **Stability:** a one-ulp FP32 perturbation of a panel, the difference
  another device or kernel makes, changed which traits greedy pivoting kept,
  because these panels hold near-exact ties (sine/cosine features of equal
  variance).

| cutoff | df under perturbation | T under perturbation |
|---|---|---|
| rounding cutoff over traits | ±1–2 | up to 7% |
| VIF ≤ 100 | ±1–2 | up to 12% |
| eigen 1e-3 or 1e-4 | fixed | within 1e-7 |

- **Result:** with outlier rows excluded and eigen truncation at 1e-3, the
  five panels gave 126–176 loci and 12–17 isolated hits (P < 5e-8 with no
  neighbour at P < 1e-5 within 100 kb). The reference has 56–196 isolated
  hits.

**Projection:** T = ‖Rz‖² with Λ_k^-1/2 U_k' = QR, where R (k × K) is upper
trapezoidal.
- **Cost:** row block [s, e) multiplies columns [s, K), so the projection costs
  what the triangular factor did.
- **Uniqueness:** R'R = U_k Λ_k⁻¹ U_k' whatever basis eigh chose inside a
  repeated eigenvalue.
- **Validated:** over the colleague's 22 groups, the dense and QR projections
  gave byte-identical locus tables.

**Missing phenotypes:** JAGWAS now accepts them. It takes the mean-imputed
panel's t with the scan's common df. z = √df·r then has that panel's Gram as
its null correlation, which is R. The full-output rescaling by √(trait df / df)
is skipped for JAGWAS, because it would leave z a null variance below R's
unit diagonal. `run_linear_gwas(phenotype_outlier_sd=...)` is an opt-in mask:
whole panel rows for JAGWAS, single values for per-trait scans.

**Calculator:** every JAGWAS ledger is per method
(`reduction_tensor_work.jagwas_cutoff_method`: `eigen` or `rounding`), and the
planner prices `TORCHGWAS_JAGWAS_RCOND`'s method.
- **Trace:** cache `jagwas.v3`. The meta trace takes k = K, so the prepare
  phase issues eigh, sqrt, div and QR after the FP64 Gram, and the reduce
  phase installs the upper factor.
- **Arithmetic ledger (conventional counts):**
  - eigh: tridiagonalization (2/3)K³ plus back-transformation K³
    multiply-adds; divide and conquer depends on deflation and is unpriced;
  - the K-value spectrum read: 8K bytes to the host;
  - scaling: K square roots and K² divisions;
  - QR: K³ − K³/3 multiply-adds.
- **Memory floor:** phenotype + correlation + max(eigenvectors + values +
  their scaled copy, scaled copy + R), which is 24K² bytes plus the phenotype.
  The Cholesky path's is 32K².
- **Preparation phases:** `upload, correlation, cast, eigh, spectrum, scale,
  qr`. Banks record their method; `spectrum` is transfer-only, like `upload`.
- **Workspace:** the `xpotrf` workspace census applies only to the rounding
  method. The eigh (`syevd`) and QR (`geqrf`) workspace is reported as
  unresolved until a census exists.
- **Fixtures recaptured on A100:** `jagwas_factor_phases_cuda{0,2}.json` and
  `jagwas_projection_geometry.json`. The chunk census is unchanged (K = 512,
  one block).

### GPUs per idle core: measured demand, not 3 CPUs per GPU (implemented)

The planner used to cap GPUs at `idle cores // 3`. On lab-2080ti, another
user's job left 6 idle cores, and a JAGWAS autotune took 2 GPUs (24.7 s) where
4 fixed shards took 16.2 s. What one GPU needs is decode for the variants it
consumes, plus its shard thread:

    cores per GPU = decode CPU per variant / GPU time per variant + 0.25

- **Decode CPU per variant:** `decode_cpu_seconds_per_variant` decodes 512
  variants through the source's own native fill. The first half warms the
  reader; the second half is timed in process CPU time, so disk waits do not
  count. It takes 19 ms. It needs `native_reader_session` (PGEN, zstd store);
  other sources keep the old rule.
- **GPU time per variant:** `gpu_seconds_per_variant` is 2·N·K over a
  measured FP32 2048³ GEMM rate, plus the JAGWAS projection FLOPs over a
  measured FP64 rate. It takes 0.1 s after CUDA setup. It counts the GEMMs
  only, so it is a floor on GPU time, and the demand errs toward fewer GPUs.
- **The 0.25 host share per GPU:** process CPU in the benchmark records.
  4 JAGWAS shards used 1.48 cores in total and 4 significant-pair shards
  about 2, and roughly 1 core of each is fixed setup.
- **Cap:** GPUs are capped at (idle cores − 1) / cores per GPU, keeping the
  existing floor of two. An explicit `autotune_options['cpus_per_device']`
  restores the fixed rule. `layout['cpu_demand']` records the inputs.

Measured on lab-2080ti, benchmark PGEN, K=8,192, 20,000 samples:
- decode 31 µs CPU per variant;
- GPU: 257 µs per variant for JAGWAS (the scan runs at 216 µs per variant on
  one GPU), 36 µs for significant pairs;
- cores per GPU: 0.37 for JAGWAS, 1.1 for significant pairs.

With 5.25 idle cores, JAGWAS autotune now takes all 7 GPUs and runs in 17.1 s.
Before, it took 2 GPUs and 24.7 s. Fixed 4 shards measured 16.2-20.3 s across runs.

**Supply: fair share, not idle cores (2026-09-25, full scale).** A busy host
is not a full one. Other users' threads time-slice with ours. At full scale
on the A100 (96 cores, load ~82, 6 idle), the idle-core supply allowed 2
significant-pair shards, and autotune ran in 348.1 s. Fixed 4 shards ran in
230.8 s, so they plainly got the decode CPU they needed. The supply is now
what our threads can attain:

    T = 5·G + 1 threads (up to 4 readers and a shard thread per GPU, plus the writer)
    attainable(G) = min(T, T·C / max(C, L + T))      C = affinity cores, L = 1-minute load − 1
    G = the largest count with G·(cores per GPU) + 1 ≤ attainable(G), at least 2

- **A100** (C 96, L 82): 4 shards attain ~19.6 cores and need ~9, so 4.
- **Oversubscribed H100** (load 97 on 48 cores): still the floor of 2.
- **Fallback:** without a load reading, the idle-core supply applies.
- **Readers per GPU:** first sized from this supply. That was superseded
  at full scale: readers now come from measured decode demand with no
  supply cap (see "Full-scale full output across input formats").

Taking a fair share on a shared host is a policy choice: the job competes
with other users' threads rather than using only idle cores.
`autotune_options['cpus_per_device']` restores a fixed cap.


### Full-scale runs (2026-09-25)

Real hardcall PGEN: 22,250 samples × 8,086,101 variants (20.8 GB). Synthetic
phenotypes: 8,192 traits and 10 covariates, the same seed on both hosts.
Executor seconds, one repeat. Earlier sections used 200,000 synthetic variants.

The A100 was shared: GPUs 0, 2, 5 and 7 each sit on their own PCIe root
port and had other users' memory resident at ~0% utilization, and the host
ran at load ~82 on 96 cores. On the H100 another user's jobs started
mid-run, taking the load from 13 to 97-120 on 48 cores; rows marked * ran
under that load.

JAGWAS (chi² identical across all layouts):

| | A100 dense | A100 triangular | H100 dense | H100 triangular |
|---|---|---|---|---|
| 1 GPU | 311.8 | 261.9 | 138.0 | 165.6* |
| 4 variant shards, chunk 1024 | 154.0 | 153.2 | 325.7* | 232.6* |
| autotune (4 shards, chunk 2048) | | 125.3 | | 32.9 (host quiet) |

- **A100, 1 GPU:** the scan is GPU-bound, mostly on the FP32 scoring GEMM
  (2·N·K per variant), and the triangular path saves 16%.
- **A100, 4 shards:** the run is decode-bound. It averaged ~2.3 busy cores,
  at 47 µs of decode CPU per variant (the process CPU matches), so the
  projection does not show. Autotune's chunk 2048 amortizes the per-chunk
  cost: 125 s vs 153 s.
- **H100:** only the quiet-host rows (1 GPU dense and autotune) are
  comparable. The others need a clean rerun.

Significant pairs (p < 1e-5, device selection; about 663,640 pairs):

| | A100 | H100 |
|---|---|---|
| 1 GPU | 364.9 | 120.0* |
| 4 variant shards | 230.8 | 106.3* |
| 4 phenotype tiles | 367.9 | 97.2* |
| autotune, idle-core cap / tile-count shards | 348.1 (2 shards) | 347.1* (2 shards, chunk 512) |
| autotune, fair-share supply | 274.1 (2 shards) | |

- **Pair counts:** they differ by up to 4 between hosts and by 1 between
  layouts. These are borderline p-values under FP32 GEMM rounding, which
  differs by kernel shape and architecture. Each host's runs agree within
  the benchmark's tolerance.
- **H100 rerun (2026-09-25, final code; load ~42-55 on 48 cores after the
  other user reduced their job; GPUs 1-4; one repeat each):**

  | | JAGWAS dense | JAGWAS triangular | significant pairs |
  |---|---|---|---|
  | 1 GPU | 103.5 | 99.8 | 118.9 |
  | 4 variant shards | 29.6 | 32.5 | 184.8 (chunk 1024, 16 readers) |
  | 4 phenotype tiles | | | 335.4 |
  | autotune | | 36.1 (4 shards) | 115.4 (4 shards, chunk 2048, 12 readers) |

  - **JAGWAS:** max |Δchi²| across layouts was 0.003 at chi² ~ 8,192, FP32
    output rounding.
  - **Significant pairs are not GPU- or decode-bound here.** 1 GPU (119 s)
    matches 4 autotuned shards (115 s), about 2x the GPU floor
    (8.09M × 6.5 µs ≈ 52 s), at ~2 busy cores. That points to a serial
    per-chunk host cost of ~15 ms per 1024-chunk (selection transfer, writer,
    GIL). Fixed shards at chunk 1024 with 16 readers do worse (185 s), tiles
    worst (335 s). Profiled and fixed on 2026-09-26: see "Significant pairs:
    the per-chunk host cost". One GPU was GPU-bound; the shards waited on the
    GIL behind the device selector's per-block round trips.
- **Remaining autotune gaps:**
  - JAGWAS on H100: 36 s vs 30-32 s (pre-scan probes and chunk probing on a
    30 s job). Closed on 2026-09-26: 28.6 / 30.0 s against 28.8-29.0 s; see
    "Tuner probing, pinned rings and dense-output multi-GPU".
  - Significant pairs on the loaded A100: autotune 250-293 s vs 209-234 s
    fixed, even with the early probe stop.
- **Paired runs (autotune and its best fixed layout, interleaved in one
  shuffled schedule, same window, A100):**
  - JAGWAS: autotune 134.6 / 139.1 s, 4 shards triangular 151.7 / 152.9 s.
  - Significant pairs: autotune 266.5 / 317.3 s, 4 shards 231.5 / 251.9 s.
  The significant-pair loss is probing on a host whose per-GPU throughput
  moved ~10x within minutes (3k-31k rows/s). A visit to chunk 512 at a slow
  moment ran at 2.5k rows/s and cost ~20 s. A probe visit now ends early once
  its size is slower than the incumbent by more than three standard errors
  and the margin (`ModelChunkTuner._clearly_worse`).
- **What the autotune rows exposed:**
  1. The idle-core CPU cap was too pessimistic on a busy host; the
     fair-share supply fixes it.
  2. Significant pairs took their shard count from the tile rule (2 at
     K = 8,192). They now use the setup-vs-work shard count over every GPU
     the CPU supply allows, and each layout is priced at its own GPU count.

### H100 layout table, quieter host (2026-09-26)

Setup:
- GPUs 1, 3, 4 and 5; load 25-29 on 48 cores; two shuffled repeats.
- Chunk selector for significant pairs.
- Results in `results/empirical_layout_20260923/{jagwas,significant}_full_h100_v3`.

Executor seconds (setup, scan and write):

| | 1 GPU | 4 variant shards | autotune |
|---|---|---|---|
| JAGWAS | 94.2, 94.3 | 30.6, 29.6 | 35.0, 37.7 (4 shards; chunk 4096, 2048) |
| significant pairs | 79.2, 79.5 | 24.4, 25.7 | 26.5, 27.6 (4 shards, chunk 4096) |

- **Outputs:** max |Δchi²| 0.003 and max |Δt| 4e-5; all keys were equal.
- **Autotune vs best fixed (API time):** 1.20× for JAGWAS and 1.11× for
  significant pairs.
- **Significant pairs vs the 2026-09-25 loaded-host table:** 118.9 / 184.8 /
  115.4 s then.

**Where JAGWAS autotune loses time** (`benchmarks/autotune_overhead_20260926.py`,
`results/autotune_overhead_20260926`):
- **Split runs (executor seconds):**
  - fixed 4 shards: 30.0, 28.7;
  - autotune restricted to chunk 1024 (layout path, no probing): 31.9, 29.0;
  - full autotune: 29.1, 36.8.
  - Across all 5 full-autotune runs today: 29.1-37.7 s, median 35.0, against
    28.7-33.0 s for fixed.
- **Before the executor:** API minus executor is 8.2-8.5 s for autotune and
  7.1-8.3 s for fixed. So pre-scan measurement (`measured_gemm_rate`, shard
  setup) costs under 1 s.
  - An nsys comparison of "first GEMM" times misleads here: autotune
    initializes CUDA before the pvar parse, so its trace clock starts earlier.
- **Inside the executor:**
  - Rings are sized for the largest candidate chunk (4096), so pinned
    allocation took 2.9 s against 0.8 s at chunk 1024.
  - Chunk 512 was the slowest size in every decision: 20-69k rows/s per GPU
    against about 100k. Under the tuner's own affine chunk model with a
    nonnegative per-chunk cost, throughput cannot fall with chunk size.
  - Re-probes reached the maximum (3) in 4 of 5 runs. They were triggered by
    host-load drift, which moves every size together.
  - The first decision compares the initial size, measured during start-up,
    against sizes probed later in steady state.

### Full-scale full output across input formats (2026-09-25)

Setup:
- Genotypes: 35,365 samples × 8,931,083 variants (UKB MRI imputed cohort)
  in every input format: BED, BED through the hard-call store, PGEN
  hardcall, PGEN dosage, zstd store and BGEN.
- Phenotypes: K ∈ {512, 2048} synthetic, with 27 covariates.
- Output: full binary (beta + t; 36.7 GB at K=512, 146.5 GB at K=2048),
  cold page cache.
- Arms: `benchmarks/direct_cold_arm.py`. The fixed arm uses the 09-15
  benchmark flags (chunk 4096, 16 workers, prefetch 32, one GPU); the
  autotune arm passes `--autotune`. `benchmarks/full_formats_20260925.sh`
  runs them, deleting outputs after each run.

Wall seconds, fixed / autotune, one repeat; A100 loads 80-95, H100 loads noted:

| format | H100 K=512 | H100 K=2048 | A100 K=512 | A100 K=2048 |
|---|---|---|---|---|
| BED, hard-call store | 30.3 / 33.1 | 86.6 / 82.7 | 461.8 / 477.4 * | 611.6 / 619.7 |
| BED | 28.1 / 29.6 | 84.6 / 83.5 | 189.8 / 190.7 | 557.0 / 591.0 |
| PGEN hardcall | 44.4 / 42.4 | 92.6 / 82.5 (load 95-120) | 195.8 / 252.2 | 609.3 / 622.6 |
| PGEN dosage (load 120-140) | 72.1 / 155.2 | 95.8 / 155.7 | 290.7 / 335.3 | 613.0 / 679.9 |
| zstd (load 135-205) | 96.1 / 303.4 | 163.9 / 271.8 | 291.8 / 317.4 | 631.1 / 681.6 |
| BGEN (load 210-250) | 71.3 / 87.8 | 185.5 / 140.4 | 280.5 / 301.0 | 676.6 / 679.8 |

\* The A100 was under memory pressure from other users' jobs: PSI memory
stalls 30-39%, ~230 compaction stalls/s, 45 GB free, 437 GB of shared
memory held. That made our process spend ~45% of its CPU time in the kernel.
Its disks were ~8% busy.

- **Versus 09-15 (torchGWAS1.1, H100, BED via the hard-call store, load
  70-85):** K=512 took 41.9-47.2 s and K=2048 98.3-114.6 s then. The current
  tree at load ~20 takes 30.3 s and 86.6 s. The runs at matched high load
  (rounds 2-3 of `results/full_benchmark_20260925`) were 67-157 s, so the
  gain at equal load is not established.
- **Autotune vs fixed:** at low H100 load, autotune is within ~10% at K=512
  and 1-11% faster at K=2048, where it picks smaller chunks for the
  writer-bound run. At load 120-250 it lost on PGEN dosage and zstd. The
  cause was the supply cap on readers below.

Autotune changes the full-scale runs forced:

1. **Chunk candidates now go up to 4096** (they were ≤2048). 4096 is what
   the fixed benchmark uses.
   - 8192 was never faster: 1M-variant PGEN fixed runs took 7.0 vs 6.0 s
     warm and 6.4 vs 5.1 s cold.
   - Full-scale autotune runs that committed to 8192 lost 60-500%. Full
     output there is disk-bound: cold input and output share one array, at
     1.4 GB/s each at 4096 and 0.7-0.9 GB/s at 8192.
   - Sequential probes run early, while the page cache still absorbs
     writes, and they favored 8192.
2. **Readers per GPU are sized from measured demand:**
   - the offered load is decode CPU per variant over GPU time per variant,
     divided by the per-thread CPU share;
   - readers are the smallest M/M/c count with P(wait) ≤ 0.2, from 2 to
     16, or 4 when the source cannot be measured (BGEN);
   - there is no supply cap. Time-slicing gives each runnable thread a
     share, so on a busy host more readers get the job more CPU. A supply
     cap held autotune to 2 readers at load 116-182, and it lost 1.7-3.2×;
   - the fitted counts reproduce the measured bests: A100 significant pairs
     → 4 (4 readers 218-231 s, 2 readers 276 s), H100 full output K=512 → 16.
3. **Decode probe:** PLINK sources (BED, with or without a hard-call store)
   are timed on their packed transport (`read_packed_into`), since the GPU
   unpacks them. The probe runs on one GPU too.
4. **Ring depth follows the GPUs in use.** Readers were first sized for the
   GPUs the CPU supply allowed. When fewer are used, the depth, which caps
   decode concurrency, is raised if the rings still fit.
5. **Framed sources get frame-multiple chunks (2026-09-26).**
   - **Symptom:** the A100 rerun's zstd rows were fixed 285 / autotune 461 s
     at K=512, and 606 / 705 s at K=2048.
   - **Cause:** the zstd store has 2,500-variant frames (the hard-call store
     has 2,048), and its pinned-loader fill decodes every frame a chunk
     touches. A frame cut by a chunk boundary is decoded by both chunks,
     through an 88 MB scratch buffer each time. Decode work is
     1 + (F − gcd(C, F))/C, the `pipeline_model` formula: 5.9× at 512, 2.2×
     at 2048 (the A100 K=512 commit) and 1.61× at 4096 (the fixed arm).
     The probes at 512-2048 paid the rest.
   - **The probe:** `decode_cpu_seconds_per_variant` timed 256 rows inside one
     frame, so the whole frame was charged to them. That is 8-10× high, and it
     fed the reader count and the CPU-demand GPU cap.
   - **Fix:**
     - `ZstdGenotype.chunk_alignment_variants` reports the frame, as
       `PlinkBedGenotype` already did for the hard-call store.
     - The default candidates become the nearest frame multiples
       (`frame_aligned_sizes`): {2500, 5000} for the zstd store and
       {2048, 4096} for the hard-call store. This applies when chunks start
       on a frame boundary; explicit chunk sizes are kept.
     - The probe times one whole aligned frame, with its buffer faulted in
       first.
   - **Measured, H100** (`results/frames_20260926`, load 74-86 on 48
     cores, cold, full binary output, one run each), wall seconds:

     | zstd store | fixed 4096 | fixed 5000 (aligned) | autotune (2500 chosen) | before the fix (v2 autotune) |
     |---|---|---|---|---|
     | K=512 | 97.7 | 45.7 | 53.2 | 303.4 |
     | K=2048 | 211.4 | 113.8 | 105.7 | 271.8 |

     An aligned chunk halves the zstd scan, against the 09-15 fixed flags
     (chunk 4096 is 1.61× decode). The hard-call store rows (fixed 4096,
     already aligned) were 36.8 / 50.4 s at K=512 and 93.7 / 100.2 s at
     K=2048: the autotune arm committed early to 2048 with depth 16. That gap
     is split into chunk and depth effects by the grid below.

6. **The ring is double-buffered (2026-09-26).** The planner set the ring
   depth equal to the readers per GPU. Every slot was then being decoded at
   once and the GPU waited on the next fill.
   - **Grid:** H100, hard-call store, full output K=512, 16 readers,
     `results/hcstore_grid_20260926`, scan seconds in two rounds (load
     77-121, then 120-140):

     | chunk × depth (variants in flight) | round 1 | round 2 |
     |---|---|---|
     | 2048 × 16 (32k), the autotune plan | 38.8 | 60.8 |
     | 4096 × 16 (65k) | 28.3 | 45.5 |
     | 2048 × 32 (65k) | 28.2 | (no free GPU) |
     | 4096 × 32 (131k), the fixed flags | 32.4 | 81.8 |

   - **Reading:** chunk size alone does not matter here. The variants in
     flight do.
   - **Rule:** depth is 2 × readers per GPU (a slot per reader being filled,
     as many filled ones waiting) when the rings fit at the largest
     candidate. Otherwise it falls back to readers.
   - **Validation:** H100, `results/depth_20260926`, load 76-93, one run
     each, wall seconds for the fixed flags / autotune. Autotune ran at
     depth 32 (16 readers), or 14 (7 readers) for the hard-call store at
     K=2048.
     - zstd: 128.0 / 60.0 at K=512, 115.3 / 100.8 at K=2048.
     - Hard-call store: 31.3 / 36.7 at K=512 (it was 50.4 at depth 16) and
       84.8 / 83.7 at K=2048.
     - The K=512 hard-call store gap is about 5 s on a 29 s scan: the probes
       commit at about 7%.
     - The fixed zstd K=2048 row varied from 115 to 211 s between sessions at
       similar load. Single runs on this host are noisy.

   A100 rerun with the reader fix (`results/full_formats_20260925_v3`, load
   90-100, before the frame fix), fixed / autotune:
   pgen_hardcall 202.8 / 234.0 (K=512) and 598.5 / 618.3 (K=2048);
   pgen_dosage 337.9 / 383.7 and 1027.5 / 782.6;
   zstd 284.9 / 460.7 and 606.2 / 705.2;
   bgen 290.9 / 271.4 and 683.9 / 672.1.

### Significant pairs: the per-chunk host cost (profiled and fixed, 2026-09-26)

Profiled on the H100 with py-spy (per thread, `--gil`) and nsys (CUDA and OS
runtime) on the full-scale input's first 2M variants: K = 8,192, p < 1e-5,
device selection, chunk 1024.
- Tools: `benchmarks/significant_host_profile_20260926.py` (driver, with a
  `--selector legacy` A/B switch), `nsys_gpu_busy.py`,
  `nsys_idle_attribution.py` (each GPU's idle time against what its launching
  thread was doing) and `speedscope_threads.py`.
- Results: `results/significant_host_profile_20260926`.

**One GPU was GPU-bound.**
- Steady state was 10.7 ms per chunk: the FP32 GEMM 7.7 ms (48.6 TFLOP/s,
  7.5 µs per variant rather than the 6.5 µs assumed above), other statistics
  kernels 1.7 ms, and selection 2.4 ms, mostly overlapped with the next GEMM.
- The ~80 ms idle gaps all fell in the first 1.4 s (input ramp-up).
- The 119 s full-scale run was ~85 s of scan plus 8-14 s of pvar parsing
  (the benchmarks pass no `genotype_cache_dir`), setup and publication.

**Four shards were GIL-bound.**
- Each GPU was 46% busy: 21 ms per chunk against 9.7 ms of GPU work.
- Per GPU, 9-10.5 ms of idle time per chunk. The launching thread spent
  4.5-5.6 ms of it in `pthread_cond_timedwait` (waiting for the GIL), ~1 ms in
  GIL handoff and 2.5 ms running Python.
- `py-spy --gil` showed each shard holding the GIL for only ~0.7 s. Hold time
  is not the cost; the wait to reacquire it after every released call is.
- Cause: the device selector ran 1M-cell blocks (8 per 1024 × 8192 chunk).
  Each block made ~45 Python tensor calls, one `nonzero` and five blocking
  `.cpu()` copies. Per chunk that was ~300 kernel launches, 47 blocking syncs
  and ~360 GIL release/reacquire pairs. Four threads convoyed on them, and on
  a loaded host each wake-up also waits for a core.

**Host selection** (still the default without autotune) costs more on one
GPU. The NumPy predicate took 25 s of a 38 s scan on the main thread, 12.6 ms
per chunk, against a 10.7 ms GPU chunk.

**Fix (`reduce.device_significant_pairs`):**
- **Block size:** a selection block is the whole chunk. It is split only at
  CUDA `nonzero`'s INT_MAX-cell limit
  (`selection_geometry.DEVICE_SELECTION_MAX_CELLS`).
  - Predicate temporaries stay within `PREDICATE_MAX_CELLS` (1M) strips that
    write one bool mask.
  - The mask costs 1 byte per cell, against the 8 bytes per cell that beta
    and t already hold.
- **Predicate:** an invalid variant gets an infinite limit, so the test is
  |t| >= limit & |t| < inf. That is exactly finite & threshold & valid.
- **Transfers:** one `nonzero`, then one packed int32 copy of row, trait,
  beta, t and df. That is 20 bytes per pair instead of 28.
- **Per chunk:** 82 Python tensor calls (was ~360), ~119 launches (was ~300),
  3 blocking syncs (was 49).
- **Output:** byte-identical on 1 GPU. Identical after sorting by
  (variant, trait) with shards, whose write order was never deterministic.

**Measured** (H100, full scale, 8.09M variants, 663,646 pairs, identical
content in every run; scan-and-write seconds, interleaved repeats; metadata
cache on):

| | legacy selector | chunk selector | GPU floor |
|---|---|---|---|
| 1 GPU | 92.8, 80.2 | 77.6, 76.8 | ~75 (9.5 ms × 7,897 chunks) |
| 4 variant shards (free GPUs) | 31.6, 41.5 | 24.7, 25.9 | ~19 |
| 4 variant shards (one GPU shared with another job) | 69.2, 53.9 | 37.9, 43.4 | |

- On the 2M prefix, per-GPU idle with 4 shards fell from 9-10.5 to 1.5-3.4 ms
  per chunk.
- The remaining wait is the statistics path's own ~40 Python calls per chunk.

**Ledgers:**
- `device_significance_work` traces the new source: one packed copy per
  nonempty block.
- The primitive map and the fixed CPU bank gained `where`, `new_empty`,
  `ge(out=)`, `lt`, `&=`, `unbind`, `view(dtype)`, `stack`, the int32 cast and
  the int32 copy.
- The A100 launch census was re-captured
  (`results/device_selection_geometry_chunk_20260926`, now the test fixture).
  All 12 rows matched the independent CPU reference.
- `selection_gpu_work` still refuses blocks over 1M cells: its CUB model was
  checked only there. Pricing whole-chunk blocks needs a census at those
  extents.
- The layout count and launch floors stay valid lower bounds, because they
  count one `nonzero` per selection block.

### Tuner probing, pinned rings and dense-output multi-GPU (2026-09-26)

**Why.** On the H100, JAGWAS autotune lost ~5 s (median 35 vs 30 s) on a
30 s job; the 36 s figure above was a single noisy repeat. Split runs
(`benchmarks/autotune_overhead_20260926.py`) put layout setup and pre-scan
measurement under 1 s, so the loss was inside the scan:
- rings pinned for the largest candidate chunk (4096): 2.9 s of
  `cudaHostAlloc` against 0.8 s at 1024;
- chunk 512 probed although the slowest size in every decision;
- the maximum three re-probes in 4 of 5 runs, triggered by host-load drift;
- a first comparison of the start size (measured during start-up) against
  sizes probed later in steady state.

**Chunk tuner (`model_autotune.ModelChunkTuner`).**
- **Upward first probe.** With period(c) = a + b·c and a per-chunk cost
  a >= 0, per-row time a/c + b cannot fall as the chunk shrinks, so only sizes
  above the start are probed.
- **Start-size revisit.** After the probes the start size is measured again,
  and that visit replaces its warmup samples in the decision. The load fit
  (within-size) keeps every sample.
- **Settling.** After a switch, `depth` completions per device are skipped,
  not one. The first new-size chunks were decoded while the GPU still ran the
  old size and complete back to back: the first size probed read
  178-248k rows/s per GPU against ~92k at its revisit and in steady state.
- **Gated re-probes.** A drift re-probes a neighbour only if the fitted model
  lets it gain more than the margin. For a larger size the bound credits the
  whole observed slowdown to per-chunk cost; a smaller size qualifies only
  when the fit's intercept is negative beyond its error (per-row cost grows
  with chunk size). Declined drifts are recorded (`skipped_reprobes`). A
  slowdown that scales every size alike changes no ranking and re-probes
  nothing.

**Pinned staging on first use.** The input loaders
(`streaming.PinnedDosageLoader`, BED `PinnedPackedBedLoader`) and the native scan's
pinned result ring pin a slot for its first chunk and grow it once, to the
capacity, when a larger chunk first needs it. Scans without a size selector
pin the capacity as before (they use it). A tuned job that never probes above
the start (short, deferred or committed) never pins the 4096-row rings.

**Dense output (reduce=None, binary beta + t).**
- **Measured (H100, full scale, K = 2,048, 132 GB written):** one GPU 120 s,
  four variant shards 45 s. One GPU is held by its single writer (1.1 GB/s);
  four shards bring four writers (2.9 GB/s together). At K = 512 (33 GB), four
  shards 15-19 s with warm input against 25 s on one GPU.
- **Planner rule replaced.** Dense panels under 4,096 traits were pinned to
  one GPU by a lab-2080ti rule (K = 64, writer-bound). The planner now uses
  `output_write_rates`: the scan's own `BinarySumstatsWriter`, one writer and
  then one per candidate GPU, fsync included, on the output filesystem, at
  most 5% of the job's own bytes (skipped below 16 MB per writer). A shard then
  serves a variant in max(GPU time, output bytes / one writer's rate), and the
  shard count minimizes `dense_shard_seconds`: setup·G + max(GPU work / G,
  bytes / min(G × one writer, the measured aggregate)). The same service time
  sizes readers and the CPU cap. Without a measurement the old rule remains.
- **Missing phenotypes.** Trait tiles and variant shards now accept them with
  the single-device contract: t of the mean-imputed panel, per-trait df in the
  manifest, no per-variant df sidecar (`trait_df` in both writers; readers and
  `open_binary_df` accept it). Before, autotune could choose a partition and
  then fail with "Full-output partitioning currently requires complete
  phenotype columns".

**Found on the way (not changed, cause not yet confirmed).** Under `/data`
fsync contention, significant and JAGWAS runs gained a 38-160 s tail after
their last chunk (fixed and autotuned alike).
- These runs did not write one part per chunk. Without a per-chunk durability
  observer (JIT productive or calibration paths), the API coalesces 262,144
  rows per part per producer (`COALESCE_ROWS`), so a full-scale run writes a
  few dozen parts at most.
- Resolved on 2026-09-27 (see "NumPy hugepage faults" below). It was not the
  disk. Writing and fsyncing `variant_ids.npy` (356 MB) takes 0.2 s. The time
  went into converting the 8.09M object IDs to fixed-width unicode while
  publishing: 75-154 s, because every large NumPy allocation faulted into
  failed THP compaction on the host's full NUMA nodes.

**Measured after the changes** (H100; `results/tuner_v3_*`, `tuner_v4_*`,
`dense_layout_*`; executor seconds, interleaved repeats):

Other users loaded the host through the evening: a foreign job on GPU 3 during
`tuner_v3_*`, bursts of fsync latency on `/data` during `tuner_v4_*`, then a
load of ~120 on 48 cores, which stopped the tmpfs rerun. The numbers are
therefore given with the condition each ran under. Every output was checked:
max |Δt| < 1e-5 on sampled rows across all layouts, same NaN pattern.

Dense output, executor seconds (full scale, 8.09M variants):

| | 1 GPU | 4 variant shards | autotune (its choice) |
|---|---|---|---|
| K=512, quiet (`dense_layout_k512`, v2) | 24.7 | 18.7 (warm input) | — |
| K=512, GPU 3 shared (v3) | 36.3, 35.7 | 23.0, 14.8 | 13.1, 14.4 (4 shards) |
| K=512 with 64 traits 2% missing (v3) | 35.3 | — | 25.8 (4 shards) |
| K=512, `/data` fsync contention (v4) | 95.1 | 285.6 | 124.4 (1 GPU) |
| K=2048, quiet (`dense_layout_k2048`, v2) | 120.4 | 45.3, 41.1 | — |
| K=2048, GPU 3 shared (v3) | 97.3, 88.1 | 119.6, 32.8 | 29.1, 70.3 (4 shards) |

- On a normal disk the planner now shards dense output, and it matched or
  beat the fixed four-shard layout in 3 of 4 runs.
- Under disk contention the write probe measured a slow aggregate, and the
  planner kept one GPU. That was right: four shards took 3x longer there.
- The missing-phenotype panel ran sharded end to end.

JAGWAS (K = 8,192), GPUs 4-7 without foreign processes (v4), executor seconds:
- fixed 4 shards: 28.97, 28.87, 28.78;
- autotune: 28.55 and 30.02, both with zero re-probes, plus one run of 194.9.
  That run's last chunk finished at 34 s; the rest was spent after the scan,
  under `/data` contention (see "Found on the way").
- Before the changes, the same comparison was 29.1-37.7 (median 35.0) against
  28.7-33.0 s.

Significant pairs are not reported: the post-scan tail (above) dominated
every run in the contended windows, fixed and autotuned alike.

### Stores, selection, factor workspace and host pages (2026-09-27)

**Index-only stores (`variant_source.py`).** Rows already carry their
variant: dense row i is input variant `variant_offset + i`, and indexed rows
store `variant_index`. A streaming store therefore no longer copies the IDs
(356 MB of `<U11` per indexed store at 8.09M variants) or the variant
metadata. Its manifest records `variant_source` instead: genotype path,
format, reopen options, offset, count, and a digest (sha256 over the count and
4,097 evenly spaced IDs, milliseconds at 8M). `store_variants` resolves IDs
and metadata through the genotype and refuses one whose list no longer matches
the digest. `sumstats_variant_ids=True` still embeds them, and in-memory
genotypes and in-memory results always do. The format is in
`docs/sumstats-format.md`.

**Device selection by default (`significance_backend.py`).** Host selection
moves dense beta and t, 8 B per cell, and filters them on one core: 12.6 ms per
1024 × 8192 chunk on the H100, against a 10.7 ms GPU chunk. Device selection
moves 20 B per passing pair. It moves less while the passing fraction, about
the threshold under the null, stays below 8/20, and every run now uses that
rule. Before, only autotuned runs did. Panels with missing phenotypes still
select on the host.

**Indexed parts close by bytes.** Parts closed at 2^18 rows: 4 MB for
JAGWAS, about 7 MB for significant pairs. They now close at 64 MiB, so the
fsync count is total bytes / 64 MiB plus one per producer. The JIT per-chunk
path is unchanged, because its write model counts chunks.

**Grouped JAGWAS on this line.** Ported from the clump branch: one genotype
pass for several panels (`jagwas_groups`). Each group is residualized as a run
of it alone, and gets its own cutoff, factor, and df column. Group columns
are remapped through phenotype QC's kept columns. On the clump branch, a
dropped trait shifted every later group and the run failed with IndexError.
That branch's streaming preparation did not record the kept columns, and it
is fixed there too (c872213, not deployed). The capacity gate models one group
factored at a time. Autotune still prices a grouped job as one panel of the
total width.

**Eigen factor workspace (`cusolver_memory.jagwas_eigen_factor_workspace`).**
The default factor is `torch.linalg.eigh` followed by `torch.linalg.qr(mode='r')`
on the kept rows. PyTorch 2.5.1 runs these as one cusolverDnXsyevd (vectors,
lower) and one Xgeqrf in place on R, both with NULL params. The untimed census
(`benchmarks/direct_jagwas_eigen_workspace_20260927.py`) queries Xsyevd per K,
and Xgeqrf for every kept row count k in [1, K], because k is only known from
the data. At K = 7, 512 and 2048 the caching allocator's peak exceeds outputs
plus the queried workspace by exactly the 512 B info tensor, for both eigh and
QR.

| K | Xsyevd | Xgeqrf (max over k) | plan setup at N = 22,250 |
|---|---|---|---|
| 2,048 | 0.10 GB | 9.4 MB | 1.01 GB (unchanged: the residualization peak dominates) |
| 8,192 | 1.62 GB | 34.6 MB | 4.34 → 5.96 GB |
| 16,384 | 6.46 GB | 68.2 MB | 12.97 → 19.43 GB |

- Xsyevd asks for about 3·K² FP64, three times the correlation itself. It is
  the largest single term of the eigen factor.
- Xgeqrf does not depend on k on the A100. On the H100 it is not monotone in k
  at K = 16,384, though the maximum is the same, so the census keeps the
  maximum.
- The two workspaces are never live together, so the plan adds the larger.
- The API's capacity gate is a necessary floor and still leaves the workspace
  out. A census covers only its own installation and trait counts.

**Selector launches above 1M cells.** Device selection runs one `nonzero` per
chunk, but `selection_gpu_work` still refused blocks over 1M cells, the old
strip bound. PyTorch 2.5.1's `nonzero` is one CUB count reduce and one flagged
select at any size below INT_MAX, so the guard is now that limit. New A100
censuses at 1024 × 2048, 1024 × 8192 and 4096 × 8192 (both harnesses take
large extents) assign every launch once:

| block | count grid | select tiles | count scratch | select scratch |
|---|---|---|---|---|
| old 1M strip: 256 × 4093 (scratch at 1,048,576) | 256 | 455 | 17,663 B | 4,351 B |
| 1024 × 2048 | 512 | 911 | 17,663 B | 7,935 B |
| 1024 × 8192 | 2,048 | 3,641 | 17,663 B | 29,695 B |
| 4096 × 8192 | 4,320 (saturated) | 14,564 | 17,663 B | 117,247 B |

The count grid is the lesser of ⌈cells / 4096⌉ and CUB's occupancy bound, 108
SMs × 8 × 5 = 4,320. Count scratch is sized for that bound, and select scratch
is (tiles + 32) × 8 B.

**NumPy hugepage faults (`host_pages.py`).** NumPy madvises MADV_HUGEPAGE on
arrays of 4 MB or more. With THP defrag set to `madvise`, a fault in such a
region compacts memory synchronously. On lab-h100, NUMA nodes 1-3 were full
of page cache (0-2 GB free; node 0 had 166 GB). The kernel had counted 811M
compaction stalls, 98% of them failed, and memory PSI sat at about 70%.

| operation | advice on | advice off |
|---|---|---|
| 8.09M object IDs → `<U9` | 75-154 s | 0.42-0.46 s |
| fresh 356 MB array, touched | 19.5-33.4 s | 0.11 s |
| fresh 1 GiB array, touched (node 0 CPUs) | 10.8-120 s | 0.36-0.45 s |
| streaming add, resident 3 × 1 GiB | 16.7-17.9 GB/s | 10.9-15.5 GB/s |

- Resident hugepages gain at most about 1.6x on streaming, and nothing on
  random gathers.
- A run allocates host arrays per chunk, part and publication, so the fault
  cost dominates.
- `run_linear_gwas` now turns the advice off for its duration and restores
  it on return or failure. It is reentrant and counted across threads, and the
  run metadata records the state. `TORCHGWAS_NUMPY_HUGEPAGE=1` keeps NumPy's
  setting.
- PyTorch's CPU allocator does not madvise hugepages by default, and pinned
  buffers are not THP-backed.
- Earlier host timings on this machine carry this cost in any large NumPy
  allocation: benchmark harnesses and CPU reference checks as well as
  torchGWAS.

**Measured after the changes** (`results/rerun_20260927/`, H100 GPUs 4-7,
from the snapshot). Outputs went to tmpfs, except dense K = 2,048, which went
to `/data`. Executor seconds, fresh interleaved processes, two repeats.
Conditions: load 9-16 on 48 cores, another user's job on GPU 3, and memory
PSI back at 50-70% (it had eased to 3% at the start). Every output was
checked: max |Δt| ≤ 9.8e-6 on sampled rows across layouts, the same NaN
pattern, and the same row counts (JAGWAS 8,086,101; significant 663,646).

| job | fixed | auto, one chunk size | autotune |
|---|---|---|---|
| significant pairs (fixed 4 shards) | 26.6, 21.6 | 23.2, 31.8 | 27.3, 25.7 (4096, 0 re-probes) |
| JAGWAS K = 8,192 (fixed 4 shards) | 27.0, 28.8 | 50.1*, 31.7 | 53.3*, 29.1 (4096, 0 re-probes) |
| dense K = 512, tmpfs (1 GPU; 4 shards) | 16.1, 15.3; 17.3, 33.0* | — | 15.9, 53.3* (2 shards) |
| dense K = 512, 2% missing | 15.8; 12.6 | — | 20.1 (2 shards) |
| dense K = 2,048, `/data` (1 GPU; 4 shards) | 402.8, 476.6; 27.0, 125.5 | — | 36.2, 92.5 (2 shards) |

- **No post-scan tail.** Significant and JAGWAS runs end within a second of
  their last chunk. Publication takes milliseconds, and embedding the IDs adds
  0.4-1.1 s. This held under memory pressure too, which is where the 38-160 s
  tails had appeared.
- **The tail A/B could not reproduce the fault.** It ran while PSI had eased to
  about 3%, so turning the advice on cost only 0.4-1.1 s (22-37 s runs, noise
  of the same size). The evidence for the cause remains the microbenchmarks
  above, which were taken under pressure.
- **Starred runs fell in slow host windows.** The JAGWAS pair ran back to back
  at 01:08-01:10. The autotuned run decided at 9.4% of the job, measuring 95k
  rows/s per GPU (the same as in the fast runs) and forecasting a 43 s job;
  throughput then fell after the decision. The dense pair each started after a
  ~70 s gap in which a 33 GB tmpfs output was freed and four CUDA contexts
  closed. With nodes 1-3 full, tmpfs output and pinned rings reclaim page cache
  as they allocate.
- **The disk was 3-4x slower for one writer than on 09-26** (403-477 s against
  88-120 s). Autotune never chose one GPU. It took 2 shards (36 and 93 s), and
  fixed 4 shards ran 27 and 126 s.

Significant pairs on the loaded A100 (`results/rerun_20260927_a100/`; GPUs 1,
2, 6 and 7; load 82-86 on 96 cores; other users' jobs at 100% on GPUs 4-5):

| fixed 4 shards | auto, chunk 1024 | autotune |
|---|---|---|
| 180.4, 146.3 | 193.2, 189.7 | 149.0 (3 re-probes, ends at 2048), 126.2 (4096) |

Under load, larger chunks win, and probing is what finds them. Measured per
GPU in the tuner's decisions: 9.8-33.7k rows/s at 1024, 13.5-33.6k at 2048,
27.6-35.1k at 4096. The A100 runs return
663,642 pairs, except the one at chunk 4096, which returns 663,646 like every
H100 run. Four pairs lie at the p = 1e-5 boundary, where FP32 rounding differs
by GPU and GEMM shape.

### Quiet-window measurements (2026-09-27, midday)

H100 GPUs 4-7, from the snapshot (pre-log10-p code). Load was 3-10, memory PSI
1-6%, IO PSI 4-15%, and another user's job held GPU 3. Three repeats, executor
seconds (`results/quiet_20260927/`).

| job | fixed 1 GPU | fixed 2 shards | fixed 4 shards | autotune |
|---|---|---|---|---|
| JAGWAS K = 8,192 | — | — | 24.6, 24.7, 24.6 | 25.2, 24.9, 26.2 |
| significant pairs | — | — | 20.0, 24.5, 20.3 | 20.7, 34.8, 20.8 |
| dense K = 512, tmpfs | 15.4, 15.3, 15.3 | 10.3, 10.1, 10.8 | 17.0, 13.2, 13.4 | 17.6, 20.4, 19.7 (2 shards) |
| dense K = 512, missing | 16.1, 15.7 | — | 14.1, 13.9 | 19.8, 18.9 (2 shards) |
| dense K = 2,048, `/data` | 91.0, 80.1 | 38.8, 35.9 | 27.9, 28.0 | 37.1, 37.9 (2 shards) |

- **JAGWAS and significant pairs** stay within 0.3-1.6 s of the fixed layout,
  apart from one 34.8 s run.
- **Dense K = 512: the planner chose the right layout.** Two shards were the
  fastest layout at 10.1-10.8 s, yet autotune took 17.6-20.4 s with them. It
  estimated a 17.9 s job, below min_job_seconds (20 s), so it never probed and
  stayed at chunk 1024, where the fixed layouts ran 4096. Its layout record
  also shows readers 26 against 16.
- **Dense K = 2,048: the planner chose the wrong layout.** The write probe
  measured an aggregate of 3.65 GB/s (one writer 2.02 GB/s) and predicted no
  gain beyond two shards. Four shards sustained 4.7 GB/s (132.5 GB in 28 s).
- **The write probe misjudges sustained rates in both directions.** On tmpfs
  it measured 5.74 GB/s for one writer where the run sustained 2.2 GB/s; with
  nodes 1-3 full, sustained tmpfs writes pay page reclaim that a 264 MB probe
  does not. On disk it measured 3.65 GB/s for four writers where the run
  sustained 4.7 GB/s.
- **The CPU-demand cap also misjudged.** It limited the K = 512 run to 2 GPUs
  (10.6 cores per GPU) on a per-variant service time of 0.71 us, where one
  GPU averaged 1.9 us.

Not yet changed: the probe and its short-job rule should be revisited with the
12-byte cells of the log10-p store (below).

### Exact -log10 P from the release (2026-09-27)

The public release stores `neglog10p.f32` beside beta and t (upstream 7c0dd01,
7abf2b1, ed8a5c9, 4a3e410), and this line now does too. The stored format,
openers, indexed fields and in-memory rows follow the release; see
`docs/sumstats-format.md`.

**The device tail.** The release computes the tail with
`upper_tail_log10_from_t_torch`. That function runs 40 continued-fraction
iterations of eager FP64 elementwise work, materializes a chunk-sized
temporary per operation, and evaluates the direct and the reflected fraction
for every cell. Timings on one 8.4M-cell chunk (H100):

| form | ms per chunk | notes |
|---|---|---|
| release, eager | 181 | 1.1 GB transient; about 17x the chunk's GEMM |
| release, `torch.compile` | 31.8 | 5.5 GB transient; 18 s to compile per process |
| one fraction per cell, compiled 4-iteration blocks | 4.1 | 4.7 s to compile per process |
| same, AOTInductor build | 3.3 | loads in 0.7 s, then 0.2 s per further GPU |

- `tails.neg_log10_p_device` evaluates one fraction per cell: the direct
  K(a, 1/2, x), or the reflected K(1/2, a, 1 - x) with the parameters swapped.
- It stays within 1.3e-11 of scipy's `stdtr` at df 22,238, and within 4.4e-7
  of the release's function.
- The stages are built per GPU architecture with `build_device_tail.sh` into
  `.build-libs/device_tail_<key>`. The key covers torch, CUDA, the
  architecture and the stage source.
- Without a build the stages compile per process: about 8 s for two GPUs,
  started at API entry. Calls then share the compile lock, because PyTorch's
  FX-tracing flag is process-wide and a compiled call during another thread's
  compile is refused.
- Loading an Inductor library turns on flush-to-zero in the loading thread,
  because it is linked with `-ffast-math` (crtfastmath). That changed NumPy
  results for subnormal thresholds in the test suite. The loader now restores
  the mode.

**Staging.** Stores stage -log10 P as float32, cast on the device. That is the
stored precision and bit-identical to a host cast, so dense output moves 12
bytes per cell to the host rather than the release's 16. In-memory results
stage FP64.

**Missing phenotypes.** The device scales t and df per trait exactly as the
host does afterwards.

**Models.** The writer, schedule, queue and planner models count the -log10 P
stream: payload, staging, write calls, fsyncs, the indexed part size (now 32
bytes per pair), the pinned ring and the device transient. The
first-principles calculator lists the tail's compute and staging as unpriced.

**Measured cost** (H100, K = 512, 2 variant shards, tmpfs output; same window,
same layout):

| | executor seconds |
|---|---|
| without -log10 P (pre-port snapshot) | 9.6, 10.9 |
| with it, per-process compile | 25.6 |
| with it, AOT build | 18.0, 19.6, 18.1 |

- Stored -log10 P against the host tail on 1.02M sampled cells: max relative
  error 6.0e-8 (float32), with the same NaN pattern.
- The store is 50% larger (49.7 against 33.2 GB), and the run is write-bound:
  the GPUs were about 60% busy.
- On writer-bound dense runs, storing -log10 P costs roughly its share of the
  bytes. GPU-bound runs pay about 3.3 ms per 8.4M cells.

### Grouped JAGWAS in the planner (2026-09-27)

The planner priced a grouped JAGWAS job as one panel of the total width. That
meant 3 x K x K FP64 of factor memory and a K-wide projection per variant, so
it refused grouped jobs whose total K has no single factor that fits a GPU,
though the API's own capacity gate accepted them. It now takes the group sizes:
- `jagwas_factor_bytes`: the largest group's three k x k matrices plus the
  other groups' retained factors;
- `jagwas_projection_flops`: the sum of each group's own projection.

Measured (H100 GPUs 4-7, full scale, K = 8,192 in 22 groups;
`results/grouped_jagwas_autotune_20260927/`):

| | executor seconds |
|---|---|
| autotune, group pricing | 31.0, 31.3 |
| autotune, one-panel pricing | 31.4, 31.2 |
| fixed 4 shards, chunk 1024 | 54.6, 55.7 |

- Both pricings chose 4 shards, since the FP32 scoring dominates at this size.
  The grouped estimate of GPU time per variant is 13% lower.
- Grouped JAGWAS favours larger chunks: 22 projections per chunk make a fixed
  per-chunk cost, which autotune's larger chunks amortize (43% faster than
  chunk 1024).

### min-P, the device tail's host calls, and pricing the scan (2026-09-27, evening)

H100 GPUs 4-7, full scale (22,250 x 8.09M hard calls), fresh interleaved
processes, executor seconds; the host was quiet (load 1.5-5).

**`reduce='min-p'` (user-facing, on request).**
- One indexed row per variant: the smallest-p trait, its beta, t, df and
  -log10 P, all from the device (`min_p.MinPReduction`).
- For a complete panel the df is per variant, so the winner is the largest
  |t| and the tail runs on the winners only.
- With missing phenotypes, df differs by trait, so the winner is ranked by
  the exact tail of every cell.
- It runs on one GPU, on variant shards or on trait tiles. Trait tiles merge
  by float64 -log10 P.
- It lives in its own module because reduce.py's bytes identify the recorded
  selector census (`device_significance_work.source_sha256`). Editing
  reduce.py failed 35 census tests.
- **Check** (`results/min_p_check_*`): on the first 200,000 variants, the
  dense store's argmax matched on 200,000/200,000 for both K = 512 panels, the
  complete one and the one with 64 traits 2% missing. -log10 P agreed within
  2e-6 relative, t was bit-identical, and one GPU matched two shards row for
  row.
- **Convention, flagged and not changed:** for missing phenotypes, full
  output's t is the mean-imputed t times sqrt(trait_df / df). At 250 of 400
  observed samples that is about sqrt(observed / n) below a complete-case fit:
  2.60 against 3.22 (`benchmarks/missing_phenotype_t_check_20260927.py`).
  min-p ranks the p-values full output reports.

**The device tail held the GIL across shards.**
- Each AOT call holds the GIL while it launches, and the block build makes 12
  calls per strip.
- A min-p winner column (4096 x 1) took 0.59 ms on one GPU and 2.89 ms per
  call on each of four (`tail_thread_scaling_20260927.py`). Four-shard min-p
  at K = 8,192 therefore ran 28.2-28.8 s, against 20.4 s for significant pairs.
- The build now also exports the whole tail as one graph. It costs one host
  call (~0.1 ms) but 0.83 ns per cell of GPU time, against 0.38 for the
  blocks: Inductor fuses the 40 iterations with recomputation. The best fusion
  option tried reached 0.56 (`tail_whole_build_options_20260927.py`).
- The block form's host time grows with the strip: 0.51 ms at 4K cells and
  1.44 ms at 4M.
- `prepare_device_tail` times both forms once per device. Each strip takes the
  form with the lower host seconds x (CUDA devices sharing the GIL) + GPU
  seconds per cell x cells. On the H100 that is:
  - one GPU: the whole graph up to about 1M cells, blocks above;
  - four GPUs: the whole graph throughout.
  - A busy A100 measured equal per-cell cost for the two forms and chooses
    the whole graph.
- Result: four-shard min-p at K = 8,192 now runs 20.7-22.0 s against 20.4 s for
  significant pairs, and at K = 512 it runs 7.0-8.5 s, down from 9.6-10.3.

**The planner priced the GEMM alone.**
- `gpu_seconds_per_variant` was 2nK / GEMM rate: 0.47 us per variant at
  K = 512. The scan's own profile shows about 3.4x that. With one GPU and 8
  readers, GPU compute was 1.58-1.65 us for min-p and 1.85 us for dense output,
  and the loader wait 0.23 us (`gpu_pipeline_profile_20260927.py`). The torch
  statistics make about ten elementwise passes over chunk x n besides the
  GEMM.
- The low price inflated CPU demand per GPU (decode 6.7 us of CPU per variant
  over 0.47 us, i.e. 14.6 cores). That held K = 512 min-p to 2 GPUs: 18.4-21.0 s,
  against 9.6-10.3 s on four fixed shards.
- `scan_device_seconds_per_variant` now times the scan's per-chunk device
  work with CUDA events on a synthetic chunk: conversion, statistics, and the
  mode's tail, reduction or selection. It is timed at 256 and 2,048 traits and
  extended linearly in K, in about 50 ms. The price is the larger of that and
  the GEMM.
- H100 prices (us per variant, dense / min-p / significant / JAGWAS):
  - K = 512: 2.26 / 1.85 / 1.86 / 1.77;
  - K = 8,192: 12.4 / 9.65 / 9.71 / 10.7.
- The narrow-panel cap (at most 2 shards below 4,096 traits, from the
  lab-2080ti) now applies only without measurements. With them, the
  setup/work model decides.

**K = 512 on four GPUs is host-bound.**
- Readers x depth factorial for min-p, chunk 4096
  (`host_stage_attribution_20260927.py`):
  - one GPU: r4d4 24.5-25.0, r8d8 21.2, r16d16 22.0, r16d32 22.2-22.7;
  - two shards: r4d4 13.3-13.5, r8d8 12.1-12.8, r16d16 14.8, r16d32 14.7-15.3.
  - Pinning a depth-32 ring of 4096-row slots costs 1.1 s.
- The chunk size dominates. At 1024 the tuner estimated a 27.9 s job, against
  7.1 s fixed at 4096. That fits a shared host cost of about 3.4 ms per chunk:
  1,975 chunks at 4096 take 6.9 s. The device-op issue time accounts for
  0.7-0.9 ms of it on one thread and 2-2.5 ms with four
  (`chunk_host_dispatch_20260927.py`).
- The tuner's warmup and revisit at its 1024 start put about 16% of the job at
  the slow sizes. A job estimated under min_job_seconds never probed away at
  all.
- The tuner now starts at the largest size whose rings fit
  (`autotune_options['start_chunk']` overrides). With a per-chunk cost a >= 0,
  per-row time a/c + b cannot fall as the chunk shrinks. Pinning, the only
  counter-cost, is paid once and bounded by the ring check.

Autotune against fixed four shards (chunk 4096, 4 readers per GPU, depth 4),
executor seconds (`results/autotune_vs_fixed_v2_20260927/`):

| job | fixed 4 | fixed 4, r8d8 | autotune (its layout) | autotune before |
|---|---|---|---|---|
| min-p, K = 512 | 7.0, 8.5 | 7.5, 8.1 | 9.5, 10.5 (4 GPUs, r6d12) | 18.4, 21.0 (2 GPUs) |
| dense, K = 512, `/data` | 11.1, 10.7 | 10.7, 12.0 | 12.9, 11.8 (2 GPUs, r6d12) | 16.6, 23.5 |
| min-p, K = 8,192 | 21.4, 22.2 | 21.4, 21.7 | 22.2, 22.6 (4, r2d4) | 30.3 |
| JAGWAS, K = 8,192 | 24.5, 24.4 | 25.6, 25.1 | 24.3, 24.3 (4, r2d4) | — |
| significant, K = 8,192 | 20.4, 20.4 | 21.0, 21.0 | 19.7, 19.7 (4, r2d4) | — |

**Still open.**
- *Readers and depth at K = 512.* Erlang-C sizes readers for each GPU running
  at its device rate, and depth doubles them. In the host-bound regime the
  shards run slower than that, and more threads or slots cost time: 6 readers
  and depth 12 against 4 and 4.
- *A per-chunk host model.* `T = max(N (s + h/c) / G, N h / c)` with h about
  3.4 ms reproduces one GPU (21.7 s modelled, 21.2 s measured), two (10.8 /
  12.1-12.8) and four (6.7 / 7.0-8.5). The device-op probe measures only about
  a quarter of h, so h itself still needs measuring.
  - That h is not the GIL. py-spy --gil on four min-p shards at K = 512
    (`benchmarks/gil_holders_20260927.py`) found the GIL mostly idle during the
    scan. Most GIL samples were process start (importing api and io); scan-time
    holders, the tail's calls and the indexed writer among them, were under
    10%. The host was loaded (load ~95), so this needs repeating on a quiet
    host, together with a timeline of where one GPU's 2.6 us per variant goes
    beyond its 1.6 us of compute.
- *Pricing every backend* (d9957f2). The device-work probe timed only the
  torch statistics, so native fused kernels and packed BED input fell back to
  the GEMM-only price. `scan_statistics_path` names the path the scan
  dispatches to, and the probe times that path's own kernels. H100, us per
  variant, dense K = 512 / 8,192: torch 2.26 / 12.4; native dosage
  1.28 / 11.2; native packed 1.19 / 11.1; BED 1.20 / 11.1; generic
  1.65 / 12.6.
- *Dense shard count.* The write model priced 2, 3 and 4 shards within 0.4 s
  of each other (10.5, 10.7, 10.9) and chose 2. The unmodelled per-shard host
  cost favours more shards.

### Missing phenotypes: complete-case OLS per trait (2026-09-27, night)

**Status (2026-09-28): opt-in, `missing_phenotype='exact'`.** The user kept
the release convention as `'impute'` and made dropping the subject the
default (see "Missing phenotypes: drop the subject" below). The method below is otherwise
unchanged. The release convention is now a plan too
(`complete_case.ImputedPlan`, rescaling t and returning the pair df), so
both conventions share the per-pair plumbing.

**The release convention was conservative.**
- For a trait with missing values, full output took the t of the
  mean-imputed panel times sqrt(trait_df / df), with pair df
  variant_df x trait_df / df.
- The mean-imputed t is already about the complete-case t. The rescale makes
  it about sqrt(observed / n) of it: 2.60 against 3.22 with 250 of 400
  observed.
- Against FP64 least squares on each trait's observed rows, the error in
  -log10 P (`benchmarks/missing_phenotype_conventions_20260927.py`):
  - the full-scale panel (64 traits, up to 2% missing): -1.3% median,
    -1.7% mean over pairs with -log10 P > 2;
  - a synthetic panel with 30% missing by a covariate: -31% median.
- Without the rescale the imputed t is -0.05% on the first panel, but still
  -5.5% median on the second. Residualizing each trait on its own rows alone
  brings both within 0.13%.

**`'exact'`: each trait on its own observed samples, as plink2's --glm does
per phenotype** (`complete_case.py`).
- Each trait with missing values is residualized on [1, covariates] over its
  own observed rows and zero-filled. The scan's product g^T y is then the
  exact numerator.
- The call's residual sum over the subset is a downdate of full-sample
  quantities over the rows the trait lacks. With Z the scan's covariate
  design [1/sqrt(n), Q]:
  (sum g^2 - sum_M g^2) - u^T (Z^T Z - Z_M^T Z_M)^+ u, where
  u = Z^T g - Z_M^T g_M.
- The df is the pair's own count less rank(Z_S) less 1.
- Traits sharing a missingness pattern share the work. Its cost is
  (missing cells) x (rank + 1) per variant, gathered in blocks of at most
  16,384 cells and 256 patterns.
- It covers every backend: torch, native fused (dosage and packed 2-bit),
  packed BED (unpacked, then the torch statistics) and generic. -log10 P,
  p-values, min-p and significant pairs all use the pair's t and df.
- JAGWAS keeps the imputed panel: its joint test needs one sample set.
- Calls missing inside a subset keep the genotype convention: centred at the
  variant's observed mean, df one lower each.

**Two traps the tests caught.**
- A float32 covariate basis is orthonormal only to ~1e-7. Downdating I
  instead of Z^T Z lifted a subset's null direction above the rank cut. And
  lstsq on the phenotype kept that direction where the call's
  pseudo-inverse dropped it: 63% error in t for a trait whose samples made
  a covariate constant. Both projections now use the plan's one
  pseudo-inverse.
- The tile and shard writers read -log10 P's presence from the chunk's
  length. A chunk sent without its df made them recompute every tail on the
  host: 19 s against 1.9 s for 200k variants.

**Measured.**
- The formula reproduces lstsq to machine precision.
- Full-scale dense output on 300 variants x 64 missing traits
  (`complete_case_full_scale_check_20260927.py`): t within 2.4e-6 absolute,
  -log10 P within 1.8e-6 relative.
- The correction costs 2.2 ms per 4,096-variant chunk for 28k missing cells
  (`complete_case_gather_probe_20260927.py`), 0.9 ms of it the gather. It was
  3.5 ms before its elementwise passes were cut.
- Full scale on four H100s, K = 512, 64 traits up to 2% missing, fixed four
  shards, host load 7-9 (`results/complete_case_cost*_20260927/`), executor
  seconds against the complete panel:
  - min-p: 8.8-8.9 against 6.9;
  - dense: 14.2-16.4 against 10.5-11.4.
- This regime is host-bound (see above). correct() issues ~80 small torch
  ops per chunk, about 0.9 ms of GIL time, and that is likely most of the
  gap.
- Open: fewer host ops (one group per chunk, or a captured graph), and a
  memory-model term for the gather transient (chunk x 16,384 floats, about
  0.5 GB with its temporaries at chunk 4096).

### GIL and process start

- **Correction (2026-09-26):** hold time understates the GIL's cost. With
  device selection, shard threads waited 4.5-5.6 ms per chunk to reacquire it
  (see the section above). The note below measured only hold time.
- `py-spy --gil` on 4 variant shards:
  - Decode readers essentially never hold the GIL (the native decoder releases
    it).
  - The four shard threads together hold it for 1.6 s.
  - The main thread holds it for about 5 s of writer and preprocessing work.
  The GIL is not the limit for these scans. Changing the glue would not help
  here; the per-GPU setup and single-core host work above matter more.
- Process start: `import torch` takes 11.5 s and `import torchgwas.api` 13-17 s
  on lab-2080ti (mostly CPU; libraries on the network filesystem). No API or
  executor timing includes this, but every CLI run pays it, and on short jobs
  it exceeds the scan.

### Missing phenotypes: drop the subject (default, 2026-09-28)

**Decision (user).** `missing_phenotype='drop_subject'` is the default.
- A sample with any missing or outlier-masked phenotype value leaves the
  analysis for every trait. Under this policy `phenotype_outlier_sd` always
  masks whole rows.
- `'impute'` keeps the release convention: the trait mean, t × sqrt(trait_df
  / df), and pair df.
- An exact per-trait method is to follow as an opt-in only. The design for it
  is the phenotype-side Qᵀg correction measured in
  `benchmarks/missing_value_scaling_20260927.py`.

The reasons: real panels rarely miss phenotypes, but they have outliers, and
those drop the whole subject. JAGWAS needs one sample set in any case.

**Exact by construction.** The drop is a sample selection. `select_samples`
is new for PGEN, where it re-decides the transport; BED and BGEN already had
it. Everything downstream then sees a complete n_kept panel: QC,
residualization, df, the JAGWAS factor and the writers.
`tests/test_drop_subject.py` checks, on every backend, that the output equals
a run on files without those samples:
- CPU;
- CUDA Torch;
- native with int8 transport;
- native with packed transport;
- BED on its Torch and native paths.

**The selection is applied on the GPU, not the host.** One thread decoding
4,096 variants at 22,250 samples (`benchmarks/pgen_subset_decode_20260928.py`):

| decode | s / chunk |
|---|---|
| every sample, one pass into the output | 0.022 |
| the reader's subset path (expand, fancy-index, table) | 1.50 |
| the cheapest host gather (one pass, then `np.take`) | 0.128 |

So a selection that keeps at least half the file's samples travels as whole
file-order rows:
- `native_physical_samples`, `native_sample_positions`;
- the unselected columns are marked missing on the device:
  - int8 and dosage rows: `index_fill_` with the transport's sentinel;
  - packed PGEN: OR 0b11 (code 3);
  - native BED: AND and OR to code 01;
- each design row sits at its sample's file position, and every other row is
  zero.

Because the unselected calls are missing, they leave the mean, Σg², the
products and the observed count, so the variant df is exact as well.

A smaller selection, such as a sub-cohort, is still decoded on the host: the
extra transfer and GEMM rows would exceed the kept work. That host path now
uses the one-pass gather, 12× faster than before.

Packed PGEN and the native BED kernel cannot gather at all. They used to
refuse any selection (int8 transport, or the Torch BED path). They now take
whole rows.

`TORCHGWAS_PGEN_WHOLE_ROWS=0` forces the host path, for A/B comparisons.

**Full scale.**
- Setup: 22,250 × 8.09M hard calls, K=512, `reduce='min-p'`, 4 H100s, chunk
  4096, 16 readers, the default int8 transport with Torch statistics.
- 1% of the subjects (222) are given one missing value each.
- Seven interleaved rounds in fresh processes
  (`benchmarks/drop_subject_full_scale_20260928.py`); the host load rose from
  13 to 30 during the run.

| config | executor s, median | best |
|---|---|---|
| complete panel | 9.60 | 7.96 |
| 1% dropped, whole rows (default) | 9.95 | 8.25 |
| 1% dropped, host gather | 34.35 | 32.29 |

- With whole rows, dropping subjects costs nothing measurable at this
  host-bound K.
- Decoding on the host costs 3.6× even with the faster gather.
- The old subset path is about 12× slower again per chunk, which extrapolates
  to roughly 185 s of decode at 16 readers.

**Open.**
- A memory-mapped panel that contains missing values is copied, keeping only
  the kept rows. That is O(n_kept × K) once, and too much for a voxel-scale
  panel.
- Stores without `select_samples` (zstd, the hardcall store) refuse
  `drop_subject` when a value is missing; the error names `'impute'`.
- Grouped JAGWAS drops the union of rows across groups, as outlier masking did
  before.
- `'exact'` on a whole-row selection reads calls at the samples' file
  positions (`CompleteCasePlan.at_positions`); on BED it keeps the gathered
  layout, since the plan gathers calls by sample.

## 3. Tuning during the run

### What the current tuner hides, and what it does not

- Measurement is hidden. Trials run on real, productive chunks; observer
  instrumentation measured +1% (scan 18.1 vs 17.9 s).
- Trials are not free. For about 25% of the job the scan runs sizes that may
  be slower. At the lab-h100 quiet-window rates (512: 104-118K, 2,048: 139K
  rows/s), that costs roughly 2-4% of the job (estimated, not measured).
- Planning is not hidden. `plan_layout` (nvidia-smi, a 0.25 s CPU sample, the
  memory model) takes about 1 s after phenotypes load and before the scan.
  It could overlap phenotype loading and residualization, which take seconds
  at voxel scale.
- Not good enough, for three reasons:
  1. It scores only segment wall throughput, so it needs a quarter of the job
     to separate three sizes under noise. Early-window choices matched the
     whole-job ranking only 0.69-0.91 of the time on the H100 panel.
  2. It ignores the per-chunk stage timings the scan already records.
  3. It commits once and never reacts afterwards. It assumed load drifts
     linearly, which on these servers it does not; load is bursty.

### Online model

Each chunk already yields a `ChunkObservation`: size, reader wall, reader CPU
and scheduler runnable-wait (from `/proc/thread-self/schedstat`), H2D,
conversion, compute and result CUDA spans, and consumer time. The model per
stage s is

    time_s(c) = (a_s + b_s · c) · load_s

- **b_s (per-variant cost), easy to measure.** Priced by the calculator:
  GEMM FLOPs over the device's rate, bytes over H2D bandwidth, decode work
  per variant. It is confirmed from the first chunks' CUDA spans and reader
  CPU time, so no separate calibration run is needed and the measurement is
  hidden in productive work.
- **a_s (per-chunk overhead), hard to measure.** Launches, synchronizations,
  the Python loop, the GIL, selection blocks. Chunk size trades exactly this
  term, so it is fitted from the run by regression over chunks. Warmup chunks
  alternate between two sizes (still productive) to identify it; with b_s
  known, even one size gives a_s.
- **load_s (contention), measured per chunk, not fitted.** CPU:
  (cpu + runnable wait) / cpu for the reader thread. GPU: measured over
  predicted compute (other users' kernels). Storage: reader wall minus CPU
  and wait. Kept as short exponentially weighted averages, because load is
  bursty.

The predicted throughput of chunk size c is c / max over stages of
time_s(c) / parallelism_s. CPU decode parallelism is readers × CPU share;
each GPU stage is serial per GPU; the writer is serial per store.

### Control loop

- Re-plan every few chunks (a few microseconds of arithmetic), and
  immediately when the prediction residual for the bottleneck stage exceeds a
  threshold.
- Switch when the predicted gain clears a margin and repays the switch cost
  (drain of in-flight chunks) within the remaining job. Then check the
  observed rate against the prediction; a miss refits that stage.
- **Memory is known before any switch.** Candidates are filtered by the
  memory model plus the job's current free memory (`mem_get_info` on its own
  context). Rings are allocated for the largest feasible candidate, so
  switching down never allocates. Growing beyond the ring happens only at a
  pass boundary, and only if the model says it fits. Implemented today: the
  start-time filter (chunk sizes whose rings do not fit are dropped, with a
  clear error if none fits).
- Levers, cheapest first:
  - chunk size and sub-tile width: any chunk;
  - reader count: any chunk;
  - fan-out on/off: per chunk inside the hub;
  - decode placement, CPU vs GPU (2-bit PGEN with GPU unpack, BGEN GPU
    decoder): at pass boundaries;
  - cache tier: at pass boundaries.

  When the CPU share collapses, the model moves toward bigger chunks, fewer
  readers or GPU decode. When a GPU slows (foreign kernels), it leaves that
  GPU fewer phenotypes in the next round.

## 4. Planner

The search space is the partition (variant shards / tile rounds / 2-D grid),
GPUs used, tile and sub-tile widths, cache tier, chunk size, readers and
fan-out, about 10^4-10^5 closed-form evaluations of `pipeline_model`. That
is tens of milliseconds in numpy. It runs in a background thread started
before phenotypes load, and it re-runs inside the control loop from updated
prices. If it later grows into real combinatorial work (assigning tiles to
heterogeneous or partly busy GPUs, a bin-packing/scheduling problem), the
candidate grid evaluates as one batched torch expression on an idle GPU; the
formulas are elementwise. Target: planning below 1% of job time and entirely
overlapped.

## 5. Implementation order, each step gated by a measurement

1. Remove `p_value_threshold`/`topk_per_trait`. Gate: existing test suites.
2. Variant shards for significant pairs. Gate: identical pairs to 1 GPU, and
   lab-2080ti/lab-a100 layout benchmarks against 2- and 4-tile shared decode.
3. Online per-stage model and control loop, replacing trial segments; keep
   the old tuner behind an option for A/B. Gate: time to decision, and job
   time against the best fixed chunk, on all three hosts under their normal
   load.
4. Planner over the partition space with the calculator's prices; current
   rules as the fallback without prices. Gate: chooses the measured-best
   layout in the existing layout benchmarks.
5. Compact genotype cache and replayed passes, then in-GPU sub-tiles. Gate:
   voxel-shaped run (K ≫ one GPU) on 1 and 2 GPUs: passes after the first
   show no decode CPU; outputs identical.
6. Shared decode per round when tiles > GPUs, if the cache does not already
   make it moot.

## 6. Questions

- Delete `p_value_threshold` and `topk_per_trait` outright?
- Is TF32/BF16 GEMM ever acceptable for screening? At voxel scale, GEMM
  precision is the largest lever (5-10x) and dwarfs everything above. The
  current design assumes FP32.
