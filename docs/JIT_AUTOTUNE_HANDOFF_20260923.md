# TorchGWAS JIT autotune handoff — 2026-09-23

## Update (2026-09-27, afternoon): in line with main, exact -log10 P stored

Details are in the design doc's "Exact -log10 P from the release" and
"Quiet-window measurements" sections.
- **Upstream commits.** main's commits since the release are on this line:
  cherry-picked where they applied (NumPy genotype input removed,
  `_QUOTE_TRIGGERS`), ported otherwise (log10-p, the top-k tests, the
  predictor accounting).
- **Stores are the release's format.** Dense stores hold `neglog10p.f32` (12
  bytes per cell), and `open_binary_sumstats` returns four values.
- **Device tail.** `tails.neg_log10_p_device`, built ahead of time per GPU
  architecture: run `build_device_tail.sh` once per host type, or it compiles
  per process.
- **Grouped JAGWAS.** The planner now prices factor memory and projection per
  group. Same layout at K = 8,192 in 22 groups (31 s), where fixed chunk 1024
  took 55 s.
- **Dense planner.** Autotune's layout was right at K = 512 but it skipped
  probing and kept chunk 1024 (17.6-20.4 s against 10.1-10.8 s fixed). At
  K = 2,048 the write probe underestimated four writers (2 shards 37 s, 4
  shards 28 s).

## Update (2026-09-27): git repo, stores, selection, workspace, host pages

Details are in the design doc's "Stores, selection, factor workspace and host
pages" section.
- **Repository.** This tree is a git repo on branch `jagwas-dev`, on top of
  origin/main 4a3e410 (one overlay commit, then one commit per change). It
  lacks the release's log10-p output (see the overlay message). It has not
  been pushed.
- **Stores are index-only.** The manifest records `variant_source` (input,
  offset, count, ID digest), and `store_variants` maps rows back to IDs.
  `--sumstats-variant-ids` embeds them.
- **Device selection is the default** while the threshold is below 8/20.
- **Indexed parts** close at 64 MiB.
- **Grouped JAGWAS** is on this line, with the QC remap.
- **Eigen factor workspace is priced.** Xsyevd asks for about 3·K² FP64
  (6.46 GB at K = 16,384), which the plan now adds when the census is in the
  profile.
- **The selector launch model covers whole-chunk blocks.** A100 census up to
  33.5M cells, with a saturated count grid of 4,320.
- **The post-scan tail is explained:** NumPy hugepage faults, not fsync.
  - On lab-h100, NUMA nodes 1-3 were full, and each large NumPy allocation
    waited on THP compaction that failed. 8.09M IDs took 75-154 s to convert
    against 0.46 s, and a fresh 1 GiB array 10.8-120 s against 0.4 s.
  - `run_linear_gwas` now turns NumPy's hugepage advice off
    (`host_pages.py`; `TORCHGWAS_NUMPY_HUGEPAGE=1` keeps it).
  - Check `/proc/pressure/memory` before trusting host timings.

## Update (2026-09-26, late): tuner probing, pinned rings, dense multi-GPU

Details are in the design doc's "Tuner probing, pinned rings and dense-output
multi-GPU" section.
- **Chunk tuner:**
  - The first probe goes upward only.
  - The start size is revisited, so its start-up samples leave the decision.
  - After a switch, `depth` completions per device are skipped.
  - A drift re-probes only when the fitted model bounds the gain above the
    margin (`skipped_reprobes` records the rest).
  - H100 JAGWAS: autotune 28.6 / 30.0 s against 28.8-29.0 s fixed, zero
    re-probes. Before: median 35 s.
- **Pinned staging:** slots are pinned on first use and grow once to the
  capacity. Untuned scans still pin the capacity.
- **Dense output:**
  - The planner prices writers: `output_write_rates` probes the scan's own
    writer, and `dense_shard_seconds` models setup, GPU work and the writes.
  - It now shards K=512/2048 panels, which were one GPU by rule. The rule
    stays only as the fallback without a measurement.
  - At K=2048, one GPU takes 88-120 s (one writer) and four shards 29-45 s.
  - Under disk contention it measured the slow disk and kept one GPU.
  - Tiles and shards accept missing phenotypes (per-trait df, as on one device).
- **Open (resolved 2026-09-27):** the 38-160 s post-scan tail of
  indexed-output runs was NumPy hugepage faults while converting the IDs
  (see the update above), not the disk.

## Update (2026-09-26): significant pairs, per-chunk host cost fixed

Details are in the design doc's "Significant pairs: the per-chunk host cost"
section.
- **One GPU was GPU-bound** at 10.7 ms per 1024 × 8192 chunk (GEMM 7.7 ms).
- **Four variant shards were GIL-bound.** GPUs were 46% busy, and shard
  threads spent 4.5-5.6 ms per chunk waiting to reacquire the GIL. The device
  selector's 1M-cell blocks cost 47 blocking syncs and ~360 Python tensor
  calls per chunk.
- **Fix:** a selection block is now the whole chunk (split only at CUDA
  nonzero's INT_MAX cells): one nonzero and one packed int32 copy. Predicate
  temporaries stay in 1M-cell strips.
- **Full-scale H100 scan, 4 shards:** 24.7-25.9 s, against 31.6-41.5 s before
  (GPU floor ~19 s). 1 GPU: 77 s against 80-93 s. Output was identical.
- **Ledgers:**
  - The source trace, the primitive bank and the A100 launch census follow
    the new selector.
  - The CUB model still covers only blocks up to 1M cells.
- **Benchmark note:** the layout benchmarks pass no `genotype_cache_dir`, so
  every run re-parses the 8.09M-line pvar (8-14 s). Compare `executor_seconds`
  or pass a cache.
- **Not changed:** without autotune the significant backend still defaults to
  host selection. It costs 12.6 ms per chunk on the main thread at K = 8,192,
  against a 10.7 ms GPU chunk.

## Update (2026-09-26, final): eigen truncation is the JAGWAS default

Details are in the design doc's "JAGWAS cutoff: eigen truncation" section.
- **Default cutoff:** T keeps R's eigen-directions above 1e-3 of the largest
  eigenvalue (`TORCHGWAS_JAGWAS_RCOND`, `run_linear_gwas(jagwas_rcond=...)`),
  and still no more than the 0.01 rounding target allows.
  - `jagwas_rcond=0` (or the variable set to 0) gives the rounding cutoff over
    traits. `jagwas_min_residual` / `TORCHGWAS_JAGWAS_MIN_RESIDUAL` adds a
    VIF threshold to it.
- **Why:** on the colleague's collinear imaging panels, a one-ulp FP32
  perturbation changed which traits greedy pivoting kept.
  - Trait dropping moved df by 1-2 and T by up to 7% (rounding cutoff) or 12%
    (VIF <= 100). Eigen truncation kept df fixed and T within 1e-7.
  - The panels' low-variance directions were heavy-tailed because a few
    samples (60-140 per panel) are extreme in many traits at once.
- **Projection:** R from QR of Λ^-1/2 U', upper trapezoidal, so its row
  blocks cost what the triangular factor's did.
  - `projection_gemm_dimensions(..., upper=)` and the reduce trace follow the
    method.
  - Over 22 real groups, the dense and QR projections gave byte-identical
    locus tables.
- **Missing phenotypes:** JAGWAS now accepts them. It takes the mean-imputed
  panel's t with the common df, whose null correlation is that panel's Gram.
- **Outliers:** `run_linear_gwas(phenotype_outlier_sd=...)` is opt-in. It
  masks a sample's whole panel row for JAGWAS and the single value for
  per-trait scans. The clumping pipeline turns it on at 5 SD.
- **Ledgers and pricing are per method** (`reduction_tensor_work.jagwas_cutoff_method`):
  - trace cache `jagwas.v3`;
  - eigen arithmetic ledger (eigh, spectrum read, scaling, QR);
  - memory floor;
  - preparation phases `upload, correlation, cast, eigh, spectrum, scale, qr`,
    with the bank recording its method.
  The planner prices `TORCHGWAS_JAGWAS_RCOND`'s method.
- **Fixtures:** recaptured on A100:
  - `jagwas_factor_phases_cuda{0,2}.json` (eigen bank);
  - `jagwas_projection_geometry.json` (upper-trapezoidal blocks at K=2048).
  The chunk census, all K=512 with one block, is unchanged.
- **Tests:** the full suite shows no new failures against the unmodified tree.
  Both run with the same 11 collection errors, from benchmark modules missing
  from the tree.
- **Not done:**
  - the cuSOLVER workspace census covers `xpotrf` only; the `syevd`/`geqrf`
    workspace is unpriced;
  - grouped JAGWAS (`jagwas_groups`, per-group residualization) exists only in
    the clumping branch.

## Update (2026-09-26, late): FP64 numerical-rank floor; port to the public clumping branch

- **Floor in the kept set.** The kept set must also satisfy
  tr(R_S⁻¹) ≤ 1/(K·ε64). That is R's FP64 numerical rank, the LAPACK
  default.
  - With FP64 statistics the rounding target alone allowed a trace of about
    1e25, so an exact duplicate passed and added a degree of freedom with no
    chi-square behind it.
  - FP32 production is unaffected: its rounding limit is about 2e5.
  - Fixtures were recaptured again on A100; the kernel rows are unchanged.
- **Port.** `~/work/torchGWAS-jagwas-clump` (branch `jagwas-loci-clumping`,
  commit `635a7f5`, local only, not pushed to GitHub) now carries the single
  reduction plus writer and API support. The server tree at
  `/data484_4/zxie3/torchGWAS-jagwas-clump` is synced.
  - The txia2 batch pipeline imports from that server tree, so it now uses
    the new reduction.
  - `run_linear_gwas(jagwas_rank_callback=...)` reports the kept rank before
    the first chunk. The overlapped clumping wrapper used the phenotype
    column count as df; it now takes df from that callback.
  - `ZstdGenotype.chunk_alignment_variants` is included.
  - The branch's own four test failures (numpy genotype input, CLI prep, the
    missing `lowrank/zstore.py`) fail the same way on its unmodified HEAD.

## Update (2026-09-26, night): one JAGWAS reduction, cutoff from a rounding target

- **Single implementation:** `reduce.JagwasReduction` is now the only
  JAGWAS reduction, living in `jagwas_projection.py`: score form, FP64 Gram,
  block-triangular projection. The dense/legacy path and
  `TORCHGWAS_JAGWAS_PROJECTION` are removed.
- **Cutoff:** it keeps the longest greedy (pivoted Cholesky) prefix whose null
  rms rounding error in T, 2·ε_z·√tr(R_S⁻¹) with ε_z = u·√N, is at most
  `TORCHGWAS_JAGWAS_T_ROUNDING` (default 0.01). This replaces the
  K·ε32 pivot tolerance and `TORCHGWAS_JAGWAS_PIVOT_TOLERANCE`.
  - ε_z was measured at 0.1-0.2 × u·√N on H100, A100 and 2080 Ti.
  - On real panels it keeps every trait of the well-conditioned ones, and
    76-99 traits of the near-singular ones.
- **Output:** the manifest's `jagwas_rank` reports the target, ε_z, the
  estimated rounding error, and the dropped traits with residual variance
  and VIF.
- **Fixtures:** the factor, projection census and chunk census fixtures were
  recaptured on A100 again. The `*_dense.json` fixtures are deleted.

## Update (2026-09-26, evening): JAGWAS uses the score statistic

Details are in the design doc's JAGWAS factor section.
- **Why:** the pivot cutoff (K·ε32) turned out to be a heuristic. On
  near-collinear real panels, the t-based quadratic form overstated strong
  hits by up to 69 (z=10) and 2,300 (z=20), because t is nonlinear in r.
  No tolerance fixed it (`benchmarks/jagwas_nonlinearity_20260926.py`).
- **What changed:** the default reduction now forms z = t/√(1 + t²/df)
  (= √df·r) before the projection. The error is then ≤ 0.017 at any effect
  size, and T is invariant to redundant traits. The dense reference keeps
  the legacy t form.
- **Calculator:** new host-API keys and bank entries. The projection and
  chunk census and the factor fixtures were recaptured on A100.
- **Measured, not implemented:** computing z from r inside the scan saves
  1-5% of chunk compute (`benchmarks/jagwas_score_path_bench_20260926.py`),
  because the FP32 GEMM dominates.

## Update (2026-09-26, later): framed sources and ring depth

Items 5 and 6 of the design doc's full-output fixes.
- **Framed sources.** `ZstdGenotype.chunk_alignment_variants` reports the
  store's frame size (2,500 variants in the full-scale store). The default
  autotune candidates become frame multiples via `frame_aligned_sizes`, and
  the decode probe times one whole aligned frame.
  - Chunk 4096 decoded 1.61× the frames and chunk 2048 decoded 2.2×.
  - H100 zstd: autotune 53.2 / 105.7 s (K=512 / K=2048), against the 09-15
    fixed flags' 97.7 / 211.4 s and fixed chunk 5000's 45.7 / 113.8 s.
- **Ring depth.** Depth is 2 × readers per GPU when it fits, instead of
  depth = readers. On the H100 hard-call store at K=512: 2048 × 16 took
  38.8 s, while 2048 × 32 and 4096 × 16 took 28 s.
- **Full suite (lab-2080ti checkout, after the JAGWAS changes):** the same
  63 failures and 11 collection errors as the baseline. One JAGWAS test that
  patched `JagwasReduction.prepare` was updated to patch the class the scan
  runs.

## Update (2026-09-26): JAGWAS factor handles collinear panels

Details and tables are in the design doc section "JAGWAS factor: FP64 Gram
and a rank cutoff".
- **R in FP64.** `TriangularJagwasReduction.prepare` forms R as an FP64 Gram
  of the scanned FP32 panel, in sample blocks. On near-singular real panels,
  the legacy FP32 Gram was the dominant error in T: up to 37%, or 0.4-1.4%
  after the cut. With the FP64 Gram it is ≤7e-6. The dense reference keeps
  the legacy Gram.
- **Rank cutoff.** Pivoted Cholesky with tol = K·ε32·max diag(R). The kept
  traits give df = r, written to the manifest along with
  `jagwas_rank`: dropped traits, residual variance and VIF.
  - The full-rank fast path is the plain factor plus one 3-scalar check.
  - The collinear fallback uses host LAPACK dpstrf.
  - Variant shards share the decision through `jagwas.spawn`.
  - `TORCHGWAS_JAGWAS_PIVOT_TOLERANCE` overrides the relative tolerance.
- **Calculator.**
  - The factor ledger traces the blocked Gram, and a `rank_check` phase was
    added.
  - The memory floor takes the larger of the Gram and solve stages.
  - The factor bench mirrors the new phases.
  - `tests/fixtures/jagwas_factor_phases_cuda{0,2}.json` were recaptured on
    A100. `attach_factor_calibration` binds the ledger's sources:
    `reduce.py`, `jagwas_projection.py` and `jagwas_blocks.py`. Any edit to
    those files needs a recapture (the bench runs in a few seconds on the A100
    GPUs 0 and 2).
- **Tests:** `tests/test_jagwas_projection.py` covers:
  - exact and near collinearity;
  - the Schur complement against an explicit solve;
  - moving the cutoff with the tolerance setting;
  - a χ²_r null, which is rejected as χ²_K;
  - devices sharing one kept set;
  - an API run with a duplicated trait (df 6 of 7, same T as without it),
    single device and variant shards.

## Update (2026-09-25, second pass): triangular projection is the default, calculator follows

- **Triangular JAGWAS projection for every run.** Autotuned or not.
  `TORCHGWAS_JAGWAS_PROJECTION=dense` selects `reduce.JagwasReduction`.
  `reduce.py` is unchanged, so the factor calibration bound to its sha still
  holds.
- **Calculator changes** (details in the design doc):
  - `reduction_tensor_work` traces the running class and binds
    `jagwas_projection.py` and `jagwas_blocks.py`.
  - `tensor_stage_service` prices one GEMM per block, each matched to its
    own launches.
  - `gemm_work` accepts a singleton GEMV in either orientation.
  - Six new JAGWAS host-API prices, re-measured on A100
    (`results/jagwas_host_primitives_triangular_20260925`).
  - `layout_compute_floor` counts the triangular FLOPs.
- **Fixtures:** the kernel census fixtures were recaptured on A100. The dense
  captures are kept as `tests/fixtures/*_dense.json`.
- **Full suite (separate checkout on lab-2080ti):** 4,013 passed. The same
  63 failures and 11 collection errors as the baseline, all pre-existing.
- **Full scale** (22,250 × 8,086,101 hardcall PGEN, K=8,192; A100 load ~82,
  H100 loaded by another user mid-run):
  - A100 JAGWAS: 1 GPU dense 311.8 s → triangular 261.9 s; 4 shards 154 s;
    autotune 125.3 s.
  - A100 significant pairs: 1 GPU 365-370 s, 4 shards 218-231 s,
    4 tiles 368-398 s.
  - H100 JAGWAS autotune: 32.9 s (host quiet); 1 GPU dense 138.0 s.
  - The other H100 rows ran under load 97-120 and need a rerun.
- **Autotune fixes from the full-scale runs:**
  - The CPU supply is our fair share of the host,
    min(T, T·C/max(C, L+T)), not idle cores.
  - Significant pairs size variant shards with the setup-vs-work model over
    every allowed GPU, rather than inheriting the tile count.

## Update (2026-09-25): JAGWAS projection, join-stage experiment, CPU-demand cap

Details and tables are in [autotune design](autotune_design_20260924.md).

- **Triangular JAGWAS projection**
  (`jagwas_projection.TriangularJagwasReduction`):
  - L⁻¹ is lower triangular, so the FP64 projection now does about 53% of the
    dense work: one GEMM per row block, written into one K x chunk buffer.
  - `reduce.py` is unchanged, and the chi² is bit-identical to dense in
    every store compared.
  - Autotuned runs use it. Non-autotuned runs keep dense, because the
    calculator prices dense.
  - `TORCHGWAS_JAGWAS_PROJECTION=dense|triangular` overrides either default.
  - lab-2080ti, K=8,192, 200k variants: 1 GPU 70.2 → 43.3 s; 2 shards
    36.5 → 24.2 s; 4 shards 23.6 → 16.2 s.
- **JAGWAS join as a pipeline stage (experiment only, not in the API).**
  With the panel split over GPUs and the join on its own stream, the split
  reached only 0.67-0.96 of variant shards. The stream recovers 1-20% over
  an inline join. The losses are work and bandwidth, which overlap cannot
  hide:
  - every GPU copies and converts every chunk;
  - peer pulls grow with the triangle;
  - GEMMs are narrower.
  JAGWAS stays full-panel.
- **CPU cap on GPUs from measured demand.** Cores per GPU = decode CPU per
  variant (a 512-variant native-fill probe) / GPU GEMM time per variant
  (measured FP32/FP64 rates) + 0.25. This replaces `idle cores // 3`. On a
  busy lab-2080ti, JAGWAS autotune went from 2 GPUs (24.7 s) to 7 (17.1 s).
  An explicit `cpus_per_device` keeps the old rule.
- **Pricing ties go to variant shards.** Like-for-like runs, with device
  selection on both layouts, found shards faster than tiles on lab-2080ti
  and H100, where the price had them within 5%. `choose_split` keeps tiles
  only when they are priced at least 5% cheaper, as with slow peer copies of
  a wide panel.
- **Shard count balances setup against work.** G minimizes s·G + W/G:
  - W = variants × measured GPU time per variant;
  - s = context plus cuBLAS setup, timed on the second GPU (lab-2080ti
    0.34 s, H100 0.2 s);
  - `autotune_options['shard_setup_seconds']` overrides s, and 0 disables
    the term.
- **Autotune vs best fixed layout (2026-09-25, loaded hosts, executor s).**
  Autotune now matches the best fixed layout on both hosts:

  | | autotune | best fixed | autotune before |
  |---|---|---|---|
  | lab-2080ti, significant | 7.0 / 8.0 | 6.8-7.4 (2 shards) | 10.8 (2 tiles) |
  | lab-2080ti, JAGWAS | 16.0 / 16.4 | 15.6-17.6 (4 shards) | 24.7 (2 GPUs) |
  | H100, significant | 2.3 / 2.3 | 2.2-2.5 (2-4 shards) | 27.3 (2 tiles) |
  | H100, JAGWAS | 1.8 / 2.5 / 4.4 | 2.1-2.8 (2 shards) | 9.0-14.7 (7 shards) |

## Update (2026-09-24): measured redesign, see [autotune design](autotune_design_20260924.md)

- **Top-k output removed.** `p_value_threshold` kept.
- **Significant pairs can shard variants** (`variant_devices`). Selection
  runs on each shard's thread over the scan's ring.
- **Indexed output is (variant, phenotype) ordered.** The writer coalesces
  parts unless the JIT path observes per-chunk writes. That per-chunk mode is
  still the writer default and is what `reduced_output_work` models.
- **Phenotype QC runs on the scan GPU** (7.0 s → 1.5 s on the benchmark;
  it grows with the panel).
- **Device pair selection runs one chunk behind on its own stream.**
  Autotuned runs use device selection; the global default is still host,
  which the JIT pricing assumes.
- **New default tuner: `model_autotune.ModelChunkTuner`.** Per-chunk samples,
  load-adjusted, re-planned on drift; `tuner='segments'` restores the old
  trials. Under load it tracks the best fixed size (H100: 4.1 s vs 3.7 s)
  where segment trials collapsed (33 s).
- **Variant shards vs phenotype tiles are priced from measured transfer**
  (`layout_pricing`, `autotune_options['split']`).
- **Variant shards copy the residualized panel GPU to GPU** instead of
  re-uploading it from pageable host memory per shard. Setup per shard went
  from 4.3 to 2.0 s; 2 shards (7.5 s) are now the fastest lab-2080ti layout
  at the benchmark shape.
- **JAGWAS stays full-panel.** It runs only when the whole panel and its
  factor fit every GPU, and multi-GPU means variant shards. A split-panel
  version was built and removed. The planner now counts the factor and stops
  with a clear error.
- **Host genotype cache for tile rounds** (`TORCHGWAS_GENOTYPE_CACHE=1`,
  opt-in). No gain with hardcall PGEN, because the rounds are GPU-bound.
- **GIL is not the bottleneck** (py-spy `--gil`). `import torch` takes 11.5 s
  on lab-2080ti, which every CLI run pays.
- **Full suite on lab-2080ti:** 4,054 passed. Every one of the 63 failures and
  11 collection errors is outside this work:
  - modules or data missing from this checkout (`benchmarks/direct_*`,
    `scripts/`, `examples/toy`, `lowrank/`);
  - A100/H100-only launch profiles on a 2080 Ti;
  - a JAGWAS factor fixture recorded before the 09-23 `reduce.py` change.

## Update (2026-09-23 evening): a working empirical autotune

The user allowed autotune to depart from first-principles pricing. There is
now a working automatic path that needs no priced profile:
`run_linear_gwas(..., autotune=True)` / CLI `--autotune`, implemented in
[`empirical_autotune.py`](../src/torchgwas/empirical_autotune.py) and wired
into the three scan sites of `api.py` (tiled/sharded dense output, single or
tiled significant pairs, multi-GPU JAGWAS). It chooses GPUs, phenotype tiles
or variant shards and reader workers at startup (live idle GPUs and CPUs plus
the existing memory model), then tunes the chunk size from the job's own
chunks with forward/reverse interleaved trial segments. See
[empirical autotune](empirical_autotune_20260923.md) for rules, measurements
and the user-time-versus-wall-time panel.

Status: the autotune, shared-decode and format tests (62 with test_adaptive_chunks) pass on lab-2080ti, lab-h100 and lab-a100 (tuned output
equals fixed output for dense, 2-GPU significant tiles and 2-GPU JAGWAS
shards). On lab-2080ti it picks the best layout in every tested mode; its
remaining overhead versus the best hand-picked layout is startup (about a
second). The priced detailed/JIT calculator below is unchanged and still the
first-principles path; neither path calls the other.

Since then:

- **Shared decode** ([`shared_decode.py`](../src/torchgwas/shared_decode.py)):
  concurrent phenotype tiles share one host decode pass instead of each
  rereading the genotype. Wall time on lab-2080ti is unchanged (decode is not
  its 4-tile bottleneck) and user CPU drops about 20% at 4 tiles.
- **NVLink fan-out** (`gpu_fanout='uplink'|'nvlink'|'peer'`): each chunk
  crosses PCIe once into the first tile GPU and the others pull it over
  NVLink. On lab-a100, whose GPU pairs share a PCIe uplink, it restores 6.7
  GB/s per GPU where per-GPU copies get 3.3; end to end it stayed within load
  noise on lab-a100 and lab-h100 because the jobs were not transfer-bound, so
  it is opt-in. Shared decode itself cut a 4-tile lab-a100 job from 26.9 s to
  17.8 s (median).
- **Memory**: chunk candidates whose rings do not fit are dropped (the run
  stops with a clear error if none fits), memory-limited tiles are balanced in
  whole rounds over the GPUs, and a fan-out root ring must fit beside the root
  tile.
- **Other input formats**: chunk switching now also works for BED (native
  packed kernel), BGEN with GPU decode (tile slicing), BGEN with CPU decode
  and disk-backed matrices (selector-driven `read_chunk` prefetch). In-memory
  arrays make autotune a no-op. `run.json` → `autotune.control_path` names the
  path used.

Harnesses: `benchmarks/empirical_layout_bench_20260923.py` (layout
comparisons including `--mode shared`, correctness cross-check) and
`benchmarks/usertime_panel_20260923.py` (signal panel). Run them on other
hosts through `~/work/torchGWAS-jagwas-dev-h100` and `-2080ti`, which are
`proj` mappings of the same shared directory with a `- *` `.rsyncignore` so
they can never push.

## Objective and current decision

Tune chunk width from the first written chunks of **one** GWAS job, reusing
immutable independent measurements from prior jobs only when their source,
hardware, execution context, artifact and age bindings remain valid. Later
allow phenotype tile and GPU allocation moves through a bounded
output-inclusive optimization. JAGWAS always keeps the entire phenotype panel
on each variant shard; it must not phenotype-partition. Dense,
significant-pair and JAGWAS output are distinct priced regimes. Cold startup
must admit the starting layout without a whole-job timing-grid search.

**Current status:** the public `initial_chunks.staged_screen` path builds and
prices up to four fixed-partition chunk candidates after first useful written
output, with bounded metadata work, exact unissued coverage, active profile
price audits and stale-frontier rebases. It writes evidence only. It cannot
change a running chunk width: the current envelope is a partial necessary
resource floor, not a finite remaining completion-time interval. Missing or
unbound prices fail only this optional screen; the admitted scientific scan
continues. No production tile or GPU reassignment is implemented.

The frozen H100 two-GPU native PGEN significant-pair job (N=35,365,
M=1,048,576, K=16,385, 27 covariates, empty significant output) was predicted
at 15.862 seconds against a 55.769-second observed median. Loaded CPU
read/decode and selection work exceeded isolated component estimates under
contention. Do not rescale all prices to fit this one job, extrapolate the
first few chunks to the full job, or interpret a partial floor as an elapsed
completion ceiling. See [the frozen diagnosis](frozen_runtime_diagnosis_20260922.md).

## Source, execution and preserved evidence

- Editable WSL source: `/home/x/work/torchGWAS-jagwas-dev` (Windows UNC
  `\\wsl.localhost\Ubuntu\home\x\work\torchGWAS-jagwas-dev`). This project has
  no `.git` directory.
- `.remote`: `lab-a100:/data484_4/zxie3/torchGWAS-jagwas-dev`. Edit the WSL
  copy only, then `proj push`, `proj run`, and `proj pull` selected outputs.
  The adapter is
  `C:\Users\X\.codex\skills\remote-lab-proj\scripts\invoke-proj.ps1`.
- Before each push, verify `.rsyncignore` is 232 bytes. It deliberately
  excludes `results/`. A file path passed to this version of `proj pull`
  becomes a directory destination and fails; pull a result directory. To
  pull JSON inside this excluded directory, temporarily add `+ *.json`
  before `- *` to the local ignore, pull, and restore its exact original
  bytes in a `finally` block. Ordinary push does not delete remote files.
- A100 Python: `/data4012/zxie3/anaconda3/envs/heart/bin/python`;
  `PYTHONPATH=src:/data484_4/zxie3/torchGWAS1.1/.deps:tests:benchmarks`.
  The last broad staged candidate/controller/profile regression passed
  **387 tests in 56.75 seconds** before the transfer probe. The probe itself
  completed on A100 after the UUID normalization change; no new broad
  regression was needed for this standalone diagnostic and documentation.
- Large source scale: the real local-data PGEN header has 8,086,101 variants
  and a 20,838,552,600-byte file on A100. Matched metadata-only staged
  controls measured 4.298–6.496 CPU seconds for eight source steps and
  4.126–7.292 CPU seconds for a two-candidate screen, with peak RSS
  571–579 MiB under shared-server load. Those controls simulated written
  output and used synthetic service rates. A later same-stage large-K
  reuse change reduced a six/twelve-tile screen from median 5.236 to
  1.336 CPU seconds without changing its conditional work. See
  [staged screen](jit_staged_partial_screen_20260923.md) and
  [large-K control](jit_large_k_staged_screen_20260923.md).

## Relevant implementation

| Area | Current code and boundary |
|---|---|
| Public first-chunk controller | [`initial_chunk_autotune.py`](../src/torchgwas/initial_chunk_autotune.py): opt-in `source_staging` plus `staged_screen`; registration does no source schedule or candidate calculation; worker starts after written output and retains `current_evidence`, stale, budget or error status. |
| Native chunk candidates | [`productive_staged_candidates.py`](../src/torchgwas/productive_staged_candidates.py): 1–4 admitted chunk widths including baseline, fixed source suffix, trait ranges and GPU ownership; live dense writer settings or validated reduced-output prices. Requires explicit shared H2D/D2H prices and never invents link capacity. |
| Compact screen | [`productive_staged_screen.py`](../src/torchgwas/productive_staged_screen.py): completed source ledger, exact unissued frontier, bounded candidates and conditional occupancy; includes available dense writer, host/device significant, and full-panel JAGWAS floors. Returns `prediction_complete=False`. |
| Issued/output state | [`productive_output_backlog.py`](../src/torchgwas/productive_output_backlog.py), writer queue and indexed-result observations described in [finite continuation design](jit_finite_continuation_design_20260923.md): useful workload bounds, still missing an atomic complete queue/GPU/final-durability state. |
| Immutable prices | [`detailed_calibration.py`](../src/torchgwas/detailed_calibration.py) and [`productive_staged_work_binding.py`](../src/torchgwas/productive_staged_work_binding.py): execution context, component artifact hashes and price targets; a matching value without an original measured target is reported unbound. `shared_transfer_capacities` is optional and older profiles omit it. |

`initial_chunks.staged_screen` requires `source_staging` and declares
`chunk_sizes`, `occupancy_scenario`, `max_partitions` (at most 16),
`max_unique_records`, `max_chunks_per_partition`, `max_cpu_seconds`,
`max_wall_seconds`, and `max_rebases` (0–2). The staged worker can retry a
stale issued/output revision within the same budget. Its result is not a
switch forecast even when named `current_evidence`. See
[live integration](jit_live_staged_screen_20260923.md) and
[native factory](jit_native_staged_chunk_candidates_20260923.md).

Other hosts: `/data484_4` is shared, so the same checkout runs on lab-h100
and lab-2080ti via `ssh <host> 'cd /data484_4/zxie3/torchGWAS-jagwas-dev && ...'`
with the same interpreter and `PYTHONPATH`. Use each host's local `/data` for
job inputs and outputs. Installed A100-bound profiles (for example the device
selector launch profile) correctly refuse to run elsewhere, so a full suite on
another host shows host-bound failures.

## New shared-transfer diagnostic and what it means

[`shared_transfer_capacity_probe_20260923.py`](../benchmarks/shared_transfer_capacity_probe_20260923.py)
records pinned H2D/D2H single- and concurrent-GPU wall/event spans, physical
UUIDs, load snapshots and CPU affinity. The final UUID-checked A100 run on
physical GPUs 0+7 used CPU 24, 16 MiB copies ×8, two warmups and eight
repeats. Median H2D was 6.62 GB/s on GPU 0, 22.08 GB/s on GPU 7 and
13.10 GB/s aggregate concurrently; D2H was 6.50, 13.85 and 12.99 GB/s.
CPU-0 binding changed GPU 7's medians to 9.46 GB/s H2D and 5.26 GB/s
D2H. Even two CPU-24 observations differed (GPU 7 D2H 19.07 versus
13.85 GB/s). The server had other active GPUs and the selected GPUs held
resident memory. These are diagnostic service observations, **not** safe
shared-capacity ceilings or reusable profile values. The raw artifacts,
topology and caveats are in
[the transfer diagnostic](jit_shared_transfer_diagnostic_20260923.md) and
`results/shared_transfer_diagnostic_20260923/`. The final report's script
SHA-256 matches the pushed source:
`519efdb88074a36f8f732bdc56fb336b3539ccfbb8e7a2fde51de635ee31f032`.

## Next work in order

1. **Transfer calibration: producer and binding done; A100 and native-job
   bracketing open.** See [transfer calibration](jit_transfer_calibration_20260923.md).
   `transfer_calibration.py` measures single, declared-link and whole-set
   pinned H2D/D2H with seeded rounds, alternating calibration/held-out halves,
   foreign-process and NUMA-placement gates, and publishes a
   `transfer_capacity` record with low/median/high scenarios only when all
   checks pass. `attach_transfer_prices` binds one scenario to
   `shared_transfer_capacities`, per-device rates and `shared_links`,
   keeping existing bindings. Qualified records exist for 2080 Ti GPUs 1,2,5
   (1+2 behind one switch: 12.1 GB/s together, the same as one GPU) and
   H100 GPUs 0,3. Records carry the producer's execution context, so run it
   with the job's own input/output paths and CPU affinity. Remaining: an A100
   record (every A100 GPU had a foreign process), and checking bound scenarios
   against transfer service inside held-out native jobs (step 5).
2. **Build a cheap finite continuation after written output.** Starting at
   one issued/output revision, price issued source/GPU work on original
   owners, live producer queues and active results, candidate unissued work,
   selector/reduction service, writer queues, dirty writeback, fsync and
   final publication. Count shared CPU, DRAM, input, output and link
   capacities once across GPUs. Use compact full/tail groups and scenario
   occupancy; do not expand hundreds of thousands of graph nodes in cold
   startup. A conservative conditional completion ceiling is required in
   addition to the existing resource-load floor.
3. **Qualify a chunk-only JIT move.** Compare baseline lower completion with
   candidate upper completion for each declared load/output scenario, after
   charging all staging, screen, measurement, switch, publication and reserve
   costs. Require fresh immutable prices, live memory admission, exact
   unissued coverage and a still-current issue/output revision. If intervals
   overlap or the job outruns the worker, retain the starting chunk width.
4. **Then extend tile/GPU moves.** Prove the new rectangles cover each
   unissued variant–phenotype pair once; keep issued work and queued output
   on original owners. Check memory, transition overhead, shared resource
   throughput and writer ownership. JAGWAS remains full-phenotype
   variant-sharded. Significant-pair large-K tiling must include nonempty
   output scenarios, selector service and NPZ parts.
5. **Falsify against native jobs.** Use matched cold read–scan–read and
   output-inclusive controls on local `/data`, with held-out dense,
   nonempty significant-pair and full-panel JAGWAS jobs. Verify exact
   statistics/output and compare the actual elapsed completion boundary,
   not scan-only or synthetic metadata times. Measure first written-output
   latency separately: the public starting admission still invokes compact
   PGEN memory/header work, while the staged candidate screen is deferred.

No currently measured result justifies enabling a production JIT layout
switch or claiming multi-GPU throughput scaling.
