# Empirical autotune: layout at startup, chunk size from the job itself

`run_linear_gwas(..., autotune=True)` (CLI `--autotune`) is a working automatic
path for chunk size, phenotype tiling / variant sharding and multi-GPU use.
Unlike the detailed/JIT calculator path it needs no priced profile: it measures
the running job. The user explicitly allowed autotune to depart from
first-principles pricing (2026-09-23); the analytical calculator is unchanged
and remains the first-principles tool.

## What it does

**Layout, once at startup** ([`empirical_autotune.plan_layout`](../src/torchgwas/empirical_autotune.py)):

- GPUs: every CUDA device with no foreign compute process, utilization at most
  20% and at least 4 GiB free (`eligible_gpus`, from nvidia-smi so no CUDA
  context is created), unless the caller names devices
  (`autotune_options={'devices': [...]}`, `device='cuda:N'`, or an explicit
  `trait_devices`/`variant_devices`/`trait_block`, which are kept).
- CPUs: idle capacity inside the process affinity, sampled over 0.25 s in a
  thread that overlaps the nvidia-smi call (`usable_cpus`), because these
  servers are shared. At most one GPU per 3 idle CPUs, but never fewer than 2
  GPUs when 2 are free (one noisy sample read a single idle core).
- Partition: JAGWAS always keeps the full panel on each GPU and shards
  variants (wide panels over every GPU, narrow panels over at most 2). Dense
  output shards variants when the panel fits one GPU, otherwise tiles
  phenotypes; narrow dense panels and filtered dense output stay on one GPU.
  Significant pairs tile phenotypes (a tile is a block of phenotype columns,
  e.g. voxels): the tile count starts at 1 below 4,096 traits, otherwise
  `max(2, ceil(K / 32,768))`, and the GPUs used are `min(GPUs, that count)`.
  Tile width is `min(ceil(K / GPUs used), widest tile that fits memory)`.
  When memory forces more tiles than GPUs, each GPU scans its tiles in turn
  (GPU i takes tiles i, i+G, ...) and the count is rounded up to a multiple of
  the GPUs, so no GPU runs a last pass alone (K = 1,000,000 on 4 GPUs with
  60,000 traits fitting: 20 tiles of 50,000 in 5 full rounds rather than 17
  in 5 rounds with 3 GPUs idle in the last). Width is not tuned by
  measurement: every extra tile rereads and retransfers the genotype, and
  GEMM efficiency is flat once a tile is a few thousand traits wide, so the
  widest tile that fits is used.
- Memory: the existing `pipeline_model` memory model (`device_ring_bytes`,
  `host_pinned_bytes`, `auto_trait_block`; within about 5% of measured peaks)
  is evaluated at the largest chunk candidate and the chosen ring depth, so
  every size the tuner may switch to fits. Candidates whose rings cannot hold
  even a 64-trait tile (JAGWAS: the whole panel) on the GPU or in pinned host
  memory are dropped, largest first, and recorded (`chunk_sizes_dropped`); if
  none fits the run stops with the sample count and free memory in the error
  instead of failing at allocation. Explicit `trait_block`/devices bypass this.
- Readers: 4 decode workers per GPU and a prefetch depth of 4 (an A/B showed 8
  per GPU at depth 8 is slower).

**Chunk size, during the scan** (`EmpiricalChunkTuner`):

- Candidates default to 512/1,024/2,048 (4,096 opt-in; multiples of the
  smallest, so switching never creates new tail shapes;
  `AlignedChunkSizeControl`). The ring is allocated once for the largest. The
  scan starts at 1,024.
- Trials start after a 3% warmup unless the measured rate says less than
  `min_job_seconds` (20 s) of scan remains; then they are deferred and
  re-checked every 1% of the job, because early chunks can run 20x faster
  than later ones on a loaded server (a frozen H100 run estimated 15.7 s at
  warmup and took 328 s).
- Each candidate then runs for one segment of real, fully productive chunks,
  in forward-then-reverse order (A B C C B A), so a linear drift in server
  load cancels. All partitions switch together; a segment starts only after
  every active GPU has delivered a chunk issued at the new size, so nothing
  older is in flight. The score is wall rows per second; process user/system
  CPU per row is recorded alongside.
- The best size is kept for the rest of the job unless it beats the starting
  size by less than `min_gain` (2%). The trial budget (warmup + segments +
  rows draining after each switch) may use at most half the job; the largest
  sizes are dropped first, and a job too short even for that keeps 1,024.
  Trials use at most `trial_fraction` (25%) of the rows.

**Shared decode across phenotype tiles** ([`shared_decode.py`](../src/torchgwas/shared_decode.py)):

- Without it every tile opens its own reader and decodes the whole genotype
  again (N tiles = N decode passes of host CPU). With it one
  `PinnedDosageLoader` decodes each chunk once and hands the pinned buffer to
  every tile scan; the buffer returns to the decoder after every tile has
  released it. All tiles see the same chunk sequence, including tuner-chosen
  sizes. On for autotuned runs when the source is decoded on the host (PGEN,
  zstd dosage store, BGEN on the CPU decoder), float32, 2+ tiles and no more
  tiles than GPUs; `autotune_options={'shared_decode': False}` turns it off,
  `TORCHGWAS_SHARED_DECODE=1` turns it on without autotune.
- GPU fan-out (`gpu_fanout`, env `TORCHGWAS_GPU_FANOUT`): by default (`auto`,
  `pcie`) each GPU copies the chunk from pinned host memory. Fan-out instead
  copies each chunk once over PCIe into a small ring (at most 4 chunks) on the
  first tile GPU, and the other GPUs pull it peer-to-peer on dedicated
  streams. `uplink` fans out when NVLink connects the tile GPUs (`nvidia-smi
  topo -m` NV# entries mapped through UUIDs, plus CUDA peer access) and two of
  them sit below one PCIe root port (sysfs PCI tree); `nvlink` whenever NVLink
  connects them; `peer` always. Autotuned runs check that the root ring fits
  beside the root tile's own rings and otherwise fall back to PCIe
  (`fanout_reason`). Opt-in, because it has not paid end to end (below).
  `run.json` → `shared_decode.transfer` records which path ran.

**Input formats.** Chunk switching goes through the source's own chunking:

| Source | Control path | Notes |
|---|---|---|
| PGEN (native hardcall/dosage), zstd dosage store | `native` | full read timing |
| PLINK BED | `packed` | needs the native packed kernel (sm_80/sm_90 build); otherwise fixed chunk |
| BGEN, GPU decoder | `device` | tiles decoded on GPU, cut to the chosen size; cuts never cross a tile |
| BGEN, CPU decoder; disk-backed arrays | `range` | selector-driven `read_chunk` prefetch (`SelectedRangeLoader`) |
| In-memory arrays | none | autotune is a no-op (one in-memory scan) |

Paths other than `native` report delivery without read timing
(`MinimalChunkObservation`), which the empirical tuner accepts; the strict
validation for the older adaptive controller is unchanged.

`run.json` → `autotune` records the layout and why, the CPU measurement, every
segment (size, rows, seconds, rows/s, user and sys CPU per row), the scores,
the choice and the fraction of the job at which it was made; `shared_decode`
records subscribers, chunks and the transfer path.

## Verification

On lab-2080ti (GPUs 5-7) and lab-h100 (GPUs 2-3):

- `tests/test_empirical_autotune.py` (11): simulated jobs check the trial
  schedule, drift cancellation, the min-gain rule, short-job skipping,
  deferred trials and multi-GPU segment starts; layout rules cover JAGWAS,
  full, significant, memory-limited, CPU-limited and filtered output. Real-GPU
  runs compare tuned output with fixed-chunk output: dense (1 GPU),
  significant pairs over two GPUs in phenotype tiles with shared decode
  (default transfer `pcie_per_gpu`), the same with fan-out forced (`peer`)
  and with `nvlink` (fans out on the H100, falls back to PCIe on the 2080 Ti),
  and JAGWAS over two GPUs in variant shards. All match.
- `tests/test_shared_decode.py` (4): every subscriber receives each chunk from
  one decode, fan-out copies each chunk once to the root GPU, a scan that stops
  early does not stall the others, decode errors reach every scan.
- `tests/test_autotune_formats.py` (5): BED (`packed`), BGEN CPU decode
  (`range`), BGEN GPU decode (`device`) and a disk-backed matrix (`range`)
  each switch through all three sizes, commit, record the control path, and
  match the fixed-chunk statistics.
- `tests/test_adaptive_chunks.py` still passes; one expectation there predated
  the `device_selection_shape` block geometry (3-variant chunks now select in
  5 blocks, not 6) and now uses that helper.

## Measured layout comparisons (lab-2080ti)

`benchmarks/empirical_layout_bench_20260923.py`: plink2 `--dummy` genotypes
(N=20,000, M=200,000, native hardcall PGEN on local `/data`, page cache warm),
10 covariates, Gaussian phenotypes. Fresh process per observation, interleaved
order, GPUs 1-7 (GPU 0 had a foreign process). Score: output-inclusive API
wall time. Every configuration's output was compared with the 1-GPU fixed
run: identical keys, max |Δt| 6.2e-6 (significant), max |Δχ²| 2.4e-4 (JAGWAS).

Final rules (v4), 3 repeats each; every output identical to the reference:

| Workload | Autotune choice | Autotune API s | Best fixed API s | Other fixed layouts |
|---|---|---|---|---|
| Significant, K=8,192 | 2 GPUs × 4,096 tiles, chunk 1,024 | 26.0, 29.9, 25.1 | 2 GPUs: 25.4, 27.4, 29.2 | 1 GPU 31-37; 4 GPUs 38-41 |
| JAGWAS, K=256 | 1 GPU | 12.1, 12.9, 12.5 | 2 GPUs: 9.6, 11.1, 12.3 | 1 GPU 11.7-13.1; 4 GPUs 12.4-13.2 |
| Dense, K=64 | 1 GPU | 6.5, 6.8, 7.0 | 1 GPU: 5.7, 5.6, 5.3 | 2 GPUs 5.4-6.8; 4 GPUs 6.1-8.3 |

Significant output is at parity with the best hand-picked layout (median
ratio 0.95) and 20-30% faster than 1 or 4 GPUs. For narrow JAGWAS two variant
shards were 5-10% faster than one GPU, which the narrow-panel rule forgoes.
On the 5-7 s dense job autotune pays about 1.2 s of startup (GPU/CPU checks,
ring sized for 2,048) — negligible on multi-minute jobs, visible here.

Earlier iteration (v2), significant pairs, K=8,192 (threshold 1e-5, 16,371 rows):

| Configuration | API seconds, 3 repeats | Median |
|---|---|---:|
| autotune → 2 GPUs × 4,096-trait tiles, chunk 2,048 (all 3) | 43.4, 29.5, 29.8 | 29.8 |
| fixed 2 GPUs, tiles 4,096, chunk 1,024 (best fixed) | 29.4, 27.3, 28.6 | 28.6 |
| fixed 1 GPU, chunk 1,024 | 37.5, 35.7, 38.8 | 37.5 |
| fixed 4 GPUs, tiles 2,048 | 46.0, 41.2, 39.0 | 41.2 |

Autotune matched the best layout on every repeat and ran within 4% of it
(scan times 17.8-18.5 s versus 17.3-19.2 s); the remainder was startup.

What the earlier iterations taught, each now fixed:

- v1 allowed 4 tiles of 2,048 when 30 cores looked idle (46.5 s): tiles
  reread the genotype, so more of them is slower until per-tile trait work is
  large. Rule now: 1 tile below 4,096 traits, otherwise 2 plus one per 32,768.
- v1 gave 2 readers per GPU from the idle-core count (31.9 s versus 26.8 s at 4
  per GPU): fair scheduling rewards more runnable readers; floor now 4 per GPU,
  and prefetch depth is raised to match (decode concurrency is capped by depth).
- The chunk tuner skipped its trials because it charged every switch a full
  ring of the largest size; budget is now per size, dropping the largest sizes
  before skipping. Very short jobs (warmup-estimated scan < 20 s) skip trials.
- JAGWAS at K=256 and dense output at K=64 were no faster on 2 or 4 GPUs than
  on 1 (decode-bound; shards add setup), so narrow panels stay on one GPU.
- Probing free memory with `torch.cuda.mem_get_info` created a CUDA context on
  every candidate GPU (~4 s over 7 GPUs); nvidia-smi values are used instead.
- `run.json` reported the loader's `reader_workers` (24) instead of the scan's;
  the scan value now wins and the loader's is kept as `genotype_open_reader_workers`.

- An overhead A/B (same 2-GPU layout, chunk 1,024, 3 repeats, identical
  output) isolated the rest: tuner instrumentation +1% (scan 18.1 s versus
  17.9 s), a ring sized for 4,096 +6% (19.0 s), 8 readers per GPU at depth 8
  +8% (19.5 s). Defaults are now 4 readers per GPU at depth 4 and candidates
  512/1,024/2,048 (4,096 opt-in via `chunk_sizes`). A trial set that is cut
  short or skipped falls back to 1,024, never to the smallest size.

## Shared decode (lab-2080ti)

`benchmarks/empirical_layout_bench_20260923.py --mode shared`: significant
pairs, N=20,000, M=200,000 hardcall PGEN, K=8,192, GPUs 1-4, chunk 1,024,
3 interleaved repeats, all outputs identical (max |Δ| 5.7e-6).

| Configuration | API s | Median | User CPU s | Sys CPU s |
|---|---|---:|---|---|
| 1 GPU | 32.8, 42.7, 35.3 | 35.3 | 25-33 | 12-14 |
| 2 tiles, separate decode | 26.0, 26.5, 24.7 | 26.0 | 32-33 | 11-12 |
| 2 tiles, shared decode | 34.3, 26.1, 23.7 | 26.1 | 29-30 (44 in the slow run) | 11-15 |
| 4 tiles, separate decode | 39.8, 36.5, 37.4 | 37.4 | 41-43 | 14-20 |
| 4 tiles, shared decode | 37.9, 36.7, 47.5 | 37.9 | 33 (54 in the slow run) | 13-20 |

Wall time is unchanged: on this host hardcall PGEN decode is not what makes 4
tiles slower than 2, so the tile rule stays as it is. Sharing does remove the
repeated decode work (4 tiles: 33 s user CPU instead of 41-43 s), which
matters on a loaded server and more for expensive decoders (zlib BGEN on the
CPU).

## NVLink fan-out (lab-h100)

`--mode fanout`, same data shape (N=20,000, M=200,000, K=8,192), H100 GPUs 2-5
(NV18 between every pair, each GPU on its own PCIe link), 3 interleaved
repeats, load average 70-100, all outputs identical (max |Δ| 4.4e-5).
Executor seconds (setup, scan and write):

| Tiles | Separate decode | Shared, PCIe per GPU | Shared, NVLink fan-out |
|---|---|---|---|
| 2 | 6.5, 5.4, 6.2 | 5.8, 6.7, 13.5 | 5.8, 6.4, 4.1 |
| 4 | 7.4, 7.5, 7.7 | 6.7, 5.1, 5.5 | 7.8, 6.3, 14.0 |

(13.5 and 14.0 coincide with load spikes.) Fan-out ties at 2 tiles and loses
at 4. A 1,024-variant hardcall chunk is 20 MB, under a millisecond on each
GPU's own PCIe link, so host transfer is not the bottleneck.

## Shared PCIe uplinks and fan-out (lab-a100)

lab-a100 (8 A100-SXM4, NVLink NV12 all-to-all) hangs GPU pairs 0-1, 2-3, 4-5
and 6-7 off one PCIe switch each (sysfs: GPUs 0 and 1 share root port
0000:00:01.1). `benchmarks/fanout_transfer_probe_20260923.py`, 256 MiB pinned
chunks, pipelined like the hub, GB/s delivered to each GPU:

| GPUs | One alone | Each copies from host | One host copy + NVLink |
|---|---:|---:|---:|
| 0+1 (same switch) | 6.7 | 3.3 | 6.7 |
| 0+2 (different switches) | 6.7 | 6.7 | 6.7 |
| 0-3 (two pairs) | 6.7 | 3.3 | 6.2-6.5 |

So per-GPU copies split a shared uplink and fan-out restores the full rate.
End to end (`--mode fanout`, N=20,000, M=200,000, K=8,192, hardcall int8
transport, GPUs 0-3, load average ~83, outputs identical), executor seconds:

| Tiles | Separate decode | Shared, PCIe per GPU | Shared, NVLink fan-out |
|---|---|---|---|
| 2 (GPUs 0+1, shared uplink) | 18.9, 19.6, 18.8 | 19.4, 21.7, 19.1 | 18.5, 19.2, 22.4 |
| 4 (GPUs 0-3) | 26.1, 27.2, 26.9 | 12.1, 19.1, 17.8 | 20.5, 22.3, 25.3 |

Shared decode is a clear win here at 4 tiles (median 26.9 s to 17.8 s; this
host is CPU-loaded, so decoding once matters). Fan-out was not: a profiled
rerun (`benchmarks/fanout_profile_20260924.py`, API seconds) had it ahead
instead, 29.9 and 29.3 s against 33.1 and 32.6 s. Its effect is inside the
load noise because these jobs are not transfer-bound: 20 MB per chunk takes
about 6 ms even at the shared 3.3 GB/s, overlapped with compute. Fan-out
therefore stays opt-in (`gpu_fanout='uplink'` applies it exactly where the
probe shows a gain). A transfer-bound case needs a heavier transport (float32
dosage is 4 bytes per sample), but PGEN dosages are read by pgenlib, which the
shared-decode path does not use; no such case was measured.

## Wall time or user time as the tuning signal (frozen H100 panel)

`benchmarks/usertime_panel_20260923.py` ran the frozen significant-pair job
(N=35,365, M=1,048,576, K=16,385, 2 H100s, tiles of 8,193, empty output) at
chunks 512/1,024/4,096, 5 fresh-process repeats each, cold private inputs,
blocking CUDA events, `NUMPY_MADVISE_HUGEPAGE=0`, while the server's load
average was 100-150. A sampler recorded process user/sys CPU and faults.

| Chunk | Executor seconds, 5 repeats | Median | CV wall/var | CV user/var |
|---|---|---:|---:|---:|
| 512 | 128, 84, 20, 125, 43 | 83.8 | 0.69 | 0.34 |
| 1,024 | 255, 508, 26, 45, 81 | 80.7 | 0.88 | 0.84 |
| 4,096 | 391, 426, 136, 401, 241 | 390.6 | 0.68 | 0.37 |

- The same configuration varied up to 25x between repeats with other users'
  load; comparing sizes across separate runs is unreliable, which is why the
  tuner compares them inside one job with interleaving.
- User time per variant was about half as variable as wall time for two of
  the three sizes, but not load-invariant: at 1,024 it ranged 18.7-171.7
  microseconds per variant across repeats (memory-bandwidth and cache
  contention inflate CPU work, not only waiting).
- Early-window (2-12% of variants) pairwise choices, 75 per metric, agreed
  with the whole-job ranking 0.79 (wall), 0.87 (user), 0.91 (user+sys), and
  with the scan-only ranking 0.69, 0.61 and 0.65. Neither signal dominates.
  The tuner therefore scores wall throughput (the objective) with
  forward/reverse interleaving, and records user and sys CPU per segment.
- The first smoke run without `NUMPY_MADVISE_HUGEPAGE=0` spent 1,139 s in the
  kernel (70 s user) before its first chunk: THP defrag=madvise plus NumPy's
  huge-page advice stalls in direct compaction on fragmented memory. A tuner
  reading user time alone would have missed a 20x slowdown.

## Full-scale frozen H100 workload

Same significant-pair job (N=35,365, M=1,048,576, K=16,385, threshold 1e-5,
158,838 rows), fresh process per observation, 3 interleaved repeats, H100 load
average 100-150. All 12 outputs matched (max |Δt| 8.1e-6).

| Configuration | API seconds (r0, r1, r2) |
|---|---|
| autotune → 2 GPUs × 8,193 tiles | 339, 249, 39 |
| fixed 1 GPU, chunk 1,024 | 340, 274, 198 |
| fixed 2 GPUs × 8,193, chunk 1,024 | 153, 256, 68 |
| fixed 4 GPUs × 4,097, chunk 1,024 | 191, 166, 35 |

Load dominated: the same layout differed 2x inside one repeat (autotune and
fixed 2 GPUs are the same layout). One GPU was slowest in every repeat; 4 tiles
won two repeats. The chunk tuner did not engage in these runs: its warmup
rate estimate (16-17 s) was far below the actual executor time, which led to
the deferred-trial fix above.

A timeline diagnostic (`benchmarks/autotune_timeline_20260923.py`) during a
quieter window:

- autotune chose 2 GPUs × 8,193 tiles; API 35.6 s, scan+write 24.4 s, first
  chunk 4.7 s after scan start, last at 24.7 s. Trials measured 512: 104-118K,
  1,024: 109-130K, 2,048: 139K rows/s and committed to 2,048 (about +17% over
  the frozen plan's 1,024).
- A second run's 0.25 s CPU sample read 1 idle core and chose one GPU (49.8 s);
  the CPU check now never trims below two GPUs.

The CLI path (`python -m torchgwas linear --autotune`, script
`benchmarks/cli_autotune_check_20260923.sh`) on lab-2080ti produced the same
16,371 rows, 2 GPUs × 4,096 tiles, and committed to chunk 2,048.

## Limits and next steps

- Layout is rule-based. The best tile count is host-dependent: 2 tiles on
  lab-2080ti (4 were 40-50% slower, with or without shared decode), while on
  the H100 4 tiles often matched or beat 2. A measured layout probe is the
  next step: run short variant windows with 2 and 4 tiles as real output
  (indexed stores carry range-relative variant indices, so window stores can
  be merged), then continue with the faster layout.
- Narrow panels: JAGWAS uses at most 2 variant shards, dense output 1 GPU.
- Shared decode covers host-decoded sources (PGEN, zstd store, CPU BGEN). BED
  packed reads and GPU-decoded BGEN still read once per tile; they are cheap
  to read, and GPU BGEN decodes on each GPU.
- NVLink fan-out is correct but opt-in. It doubles delivered bandwidth on
  lab-a100 pairs that share an uplink, but no measured job was transfer-bound,
  so end-to-end it stayed within load noise. Turning it on automatically
  should wait for a priced comparison of per-chunk transfer time (bytes x GPUs
  per uplink / uplink bandwidth, from transfer calibration) against compute.
- With more tiles than GPUs (voxel scale), tiles do not share decode: each
  still decodes the genotype itself. Sharing one pass per round of G tiles is
  a natural extension.
- BED chunk switching needs the native packed kernel; without it BED gets the
  layout with a fixed 1,024 chunk.
- Startup adds about a second (nvidia-smi, 0.25 s CPU sample overlapped, ring
  sized for 2,048). Explicit settings still bypass autotune entirely.
