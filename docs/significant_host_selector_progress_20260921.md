# Bounded host significant selection, 2026-09-21

The calculator is still not qualified for automatic reduction tuning. The
completed 60-observation H100 panel selected B1024/T8193/two GPUs at 101.160603 s,
versus the best observed B1024/T4096/two GPUs at 89.471496 s: 13.1% median regret.
No correction was fitted to those association timings.

## Executor change and model agreement

`host_significance.py` now constructs the significance predicate in rectangles
of at most 1,048,576 cells, retaining one mask and one selected result per source
chunk. FP32 statistics compare against FP64 cutoffs rounded upward to FP32;
FP64 statistics retain FP64 cutoffs. This preserves inclusive comparisons for
representable statistics. Flat-index extraction avoids the slow two-dimensional
empty-mask scan. Coordinate buffers are reused, and payload gathering precedes
in-place global variant rebasing. Pair-specific df, arbitrary input strides,
nonempty output, ordering and queue ownership remain intact.

The development host calculator counts cutoff rounding, each predicate block
including tails, flat extraction, coordinate divmod, payload gathers and the
remaining index transformations. It refuses the former primitive bank: the
bank must identify `bounded_flat_v1` and the same predicate limit. The existing
full-retention memory allowance remains conservative; no sparse-occupancy
memory discount was introduced. The source ledger is for native FP32 PGEN,
while the generic executor also preserves FP64 and broadcast/pair-df behavior.

Significant NPZ accounting now includes local ZIP-header rewrites in submitted
page-cache bytes and logical host traffic. Final extent/storage transfer still
counts each header once. Counting-stream tests verify the actual NumPy archive,
its reopened arrays, submitted bytes, CPU/DRAM work and final storage extent.

## Completed controls and tests

H100 job `20260921-164810-615973` completed the DMA and flat-coordinate controls.
There was no consistent large penalty attributable to freshly GPU-written host
buffers. Flat coordinates helped empty/sparse masks, but the coordinate-only
experiment did not establish the full selector's dense-output performance.

The complete-selector controls, job `20260921-170048-622345`, tested 32 cases:
1M/8M cells, main/two-worker execution, CPU-initialized/freshly GPU-written
buffers, and empty/dispersed sparse/clustered sparse/dense output. Seven paired
randomized algorithm repetitions and three calls per repetition were retained.
Every complete result was compared exactly before timing. GPU transfers and
returned-array destruction are outside this CPU-selector boundary. Import
order and actual NumPy hugepage, Torch thread and loaded numerical-pool state
were recorded to match the frozen pipeline's host settings.

Illustrative completed results below are per-worker CPU medians from the
fresh-DMA, two-worker cases; these are component timings, not full-GWAS speedups:

| Cells / output | Original selector | Bounded selector |
| --- | --- | --- |
| 1M / empty | 3.81-3.86 ms | 1.07-1.08 ms |
| 1M / dispersed sparse | 3.86-3.88 ms | 1.23 ms |
| 8M / empty | 40.32-40.33 ms | 7.41-7.48 ms |
| 8M / dispersed sparse | 39.06-39.17 ms | 8.67-8.77 ms |
| 8M / dense | 335.63-335.76 ms | 222.72-222.98 ms |

The full table also retains substantial load variation; notably, some old
8M empty/dense observations were much slower. No unusually large outlier ratio
is used as the performance claim. Dense 1M main-thread performance was similar
in one case; it was not assumed to improve by the sparse-output ratio.

The initial attempt `20260921-165920-619905` stopped during test collection:
the sync filter omitted a helper without the `direct_` filename prefix. No
performance observations were taken. The helper name/filter were corrected,
and the matched run above passed 44 candidate correctness tests before timing.

A100 development job `20260921-170727-625325` passed **138 tests and 15 subtests**
in 88.10 s, with no failures or skips. Tests cover production selector precision,
noncontiguous/broadcast inputs, bounded temporary sizes, retained-array ownership,
calculator geometry and old-bank rejection, actual archive submissions, native
PGEN missing-call/phenotype behavior, threshold-one output, trait blocking and
real native multi-GPU output. Four expected prefetch-depth warnings were kept.
Earlier isolated archive/model tests passed 32 tests in 45.14 s.

A100 job `20260921-171029-626133` completed new fixed generic selector primitives.
All declared bulk-minus-fixed CPU differences were nonnegative; none was clipped.
The raw record is `results/bounded_host_selection_prices_20260921/selection_prices.json`.
All 104 local package hashes matched its remote source capture after pulling.
These are A100-host observations, not H100 prices. Concurrent transfer remains
unqualified, and archive/queue/storage prices must be supplied separately.

## Running full-executor control and remaining work

H100 job `20260921-171303-626976` (PID 2409277 at launch) is the paired cold-input
executor A/B. It preregisters three repetitions of original versus bounded
selection at each of the two finalist tile sizes, for 12 fresh processes.
It reuses the same private local PGEN/phenotypes, cold read/API/read controls,
independent FP64 sampled reference, output durability and 4.5-GiB allocation cap.
Only the tested selector function/helper are installed inside proposed child
processes; main package source and the original frozen benchmark stay intact.
The proposed helper, development linear source, passing test XML, harness and
protocol are hashed. This experiment tests the implementation change; it is
not an all-candidate autotuning qualification.

The root is the main project's `results/significant_bounded_selector_ab_20260921`.
Its final report is pending. Fresh H100 primitive binding, a new frozen ranking
panel, nonempty sustained output, device filtering and additional formats remain
required. The modified development `linear.py` also invalidates old whole-file
geometry audit bindings; recapture/rebind that source before claiming a new
JAGWAS geometry audit. Prior captured observations remain immutable.

## Refreshed source audit and H100 checkout

A100 job `20260921-171742-629147` completed new duration-free geometry captures,
67 passing tests, and the source audit. Its report is
`results/jagwas_bounded_host_source_audit_20260921/replay/audit.json`: 104 source
files, 12 tensor geometry checks, 20 exact significant graph replays, 6 dense
scan replays, 5 trait-shape replays and 3 FP32 service replays. Both linked test
reports have zero failures, errors and skips. All output directories were
pulled; their source hashes match the A100 development checkout and the new
H100 checkout. Old records and fixtures were preserved. The audit now accepts
explicit new geometry paths instead of requiring overwritten fixed filenames.

The new local WSL `/home/x/work/torchGWAS-calculator-h100` project maps to
`lab-h100:/data484_4/zxie3/torchGWAS-calculator-h100`. Its documentation records
job `20260921-172729-636792`, which waits for the A/B job before fresh independent
H100 selector pricing, source/context binding, and a new all-candidate panel.
The first matched A/B pair was 51.565884 s bounded versus 100.090220 s original
at B1024/T4096/two GPUs; the B1024/T8193 bounded observation was 34.427044 s.
These are partial observations, not final speedup or tuning qualification.

## Completed full-executor selector A/B

H100 job 20260921-171303-626976 completed all 12 observations (two fixed
configurations, three randomized paired repeats). With chunk size 1024 and two
GPUs, median scan-and-durable-indexed-write time changed from 138.185170 s to
51.565884 s at tile width 4096 (2.679779x), and from 87.875962 s to 28.132105 s
at tile width 8193 (3.123690x). Every observation checked all 140 independent
FP64 reference cells. This is implementation A/B evidence under the recorded
load/cache context, not all-candidate selection qualification.

The raw observations and final report were pulled to the main checkout at
`results/significant_bounded_selector_ab_20260921`. Fresh independent H100
bounded-selector prices and their source/context binding completed in the
isolated H100 checkout. Full ranking preparation then stopped at publication
because three helper files were missing at relative paths. The completed
geometry collection was retained; exact byte copies of the original helpers
and a checked resume path repair that packaging error without changing prices,
candidate axes, source executor, repetitions or gates. Resumed job
20260921-180931-650226 froze 12 feasible candidates and selected candidate 17
(B1024/T8193/two GPUs, predicted tile executor 15.862022 s) before observations.
The large absolute-time difference from the A/B is not fitted away; the panel
still tests the predeclared selection-regret and timing-boundary criteria.
