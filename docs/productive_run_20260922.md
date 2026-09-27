# Productive execution boundaries for incremental planning

`ProductiveTuningRun` connects the existing aligned chunk control, bounded
planning budget and computational cache to real read reservations and indexed
writer completion. It is an internal execution bridge. Its caller must first
admit the shared allocation envelope, fixed partitions and every permitted
chunk/tail shape. It does not replace admission, calibrate component prices, or
provide the public automatic policy.

## Actual source prefixes

Every partition has a stable id, device, variant interval and phenotype interval.
The selector is called after the pinned reader acquires a free ring slot. At
that point the selected interval is reserved and entered in the prefix, before
its native read is submitted. It therefore includes unsampled reads, prefetch
and submission about to occur. A later chunk-size choice cannot change any of
these intervals. A failed submission fails the scan; its reservation is not
silently removed to make a different history.

The variant-sharded driver binds unique scan partitions before shared
preprocessing. Phenotype tiles with otherwise identical device, variant range
and width require explicit partition-id binding, preventing them from sharing
one source cursor. Every reservation checks its exact contiguous cursor, stop
and allocated capacity. Source reads and writes continue after the tuning
budget closes, using the final admitted size.

Only the early prefix is retained. At the declared prefix limit, tracking
closes tuning and marks the recorded prefix incomplete for further decisions.
It then advances scalar counters without accumulating a full-job range list.
The number of planning decisions, early window and computational cache are
also bounded. Caller admission must include host reserves for this bookkeeping,
cache copies and calculator temporaries; the cache's Python-size estimate is
not an allocator or process-RSS bound.

## Written work is a separate boundary

`write_indexed_sumstats(..., on_chunk_written=...)` reports an immutable
`IndexedChunkWrite` after processing each JAGWAS or significant-pairs writer
block. A material part is reported after its file flush and requested `fsync`.
The event contains source range, surviving rows, physical part-file bytes and
start/completion clocks. It does not claim final manifest or directory
durability. A failed write emits no completion event.

An empty significant-pairs block still represents useful completed work and
can open the planning window. It reports zero rows and bytes, no part name and
no file-fsync claim. First processed output, first material part and first
fsynced part are distinct timestamps. Existing sampled scan events still refer
to downstream iterator/queue acceptance and must not be substituted for these
writer events.

## Bounded planning during the run

No planning callback is invoked before the first writer completion. Each
subsequent call to `planning_step` can evaluate at most one proposed chunk size
using the existing calculator. The callback receives a snapshot of the actual
reserved prefix. It must return baseline and candidate completion forecasts
for equivalent remaining work under the same scientific output request.

The issue frontier is held synchronously during this one calculation, while
already-issued GPU and reader work continues. The callback must not wait for
the executor. This prevents a proposal from silently applying to a different
source prefix, but it can stall future input submission; moving planning here
is not free. The budget charges actual planner-thread CPU and wall time.
Changing size requires a timely result and a forecast gain greater than all
recorded planning wall time plus declared switching, publication and reserve
costs. Structural persistence requires an explicit publication-cost forecast
before cache access. The follow-up `window_forecast_20260922.md` documents this
cost accounting and the remaining-work scenario API. Invalid or failed
optional planning closes tuning and leaves the valid scan running. Late
results, expired windows, exhausted prefixes and fully issued sources cannot
authorize a change. Scientific modes and thresholds are not proposal fields.

The policy supplies forecasts and their evidence. This bridge does not turn
loaded CUDA intervals into independent capacities, implement a Bayesian
posterior, or infer live GPU/queue service state from the source prefix.
Analytical checkpoints remain model states. Mapping observations to an
identified predictive state and choosing which candidate to evaluate remain
necessary for the automatic policy.

## Verification scope

Unit and integration tests cover unsampled/prefetched reservations, exact
partition binding, held-frontier concurrency, mutation isolation, early-window
closure, unprofitable/invalid/late proposals, empty indexed output, requested
file fsync and failed writes. The existing multi-GPU, output and adaptive scan
tests exercise unchanged callers as well.

`direct_productive_run_20260922.py` runs the retained real native-PGEN input on
two A100 GPUs with a full-width JAGWAS factor per GPU and durable indexed output.
It compares alternating-order pairs with fixed chunk sizes and lifecycle
tracking disabled/enabled, after one explicit warmup. Both arms use the same
minimal writer timestamp callback. It reports first completed part and final
output-inclusive time separately. Admission is outside these scan timings;
this is not a measurement of total public-API startup or cold production work.

A separate run changes future chunk sizes using scripted completion forecasts
to exercise the bridge. It checks that the logged reservations equal the
actual delivered source ranges, no variant is omitted or repeated, and all
statistics match the fixed run within the existing FP32 tolerance. Scripted
forecasts are test controls, not calculator accuracy or autotuning benefit.

The first validation passed 96 tests with no failures or skips. The final
regression run passed 105 tests and five subtests, also without failures or
skips (`results/productive_run_v2_20260922`). Expected reader-depth warnings
are retained in the test log.

The final live audit is `results/productive_lifecycle_v2_20260922/report.json`.
Both input and association-output mounts were verified as local XFS on
`/dev/md0`, mounted at `/data`. Reports remain in the shared results directory.
All 121 package-source hashes, the benchmark and three helper hashes, and the
prior execution-report hash match the delivered source. The five immutable
input-file hashes were checked on the server before and after execution.

Seven fixed-size matched pairs gave these warm-fixture measurements:

| Metric | Fixed median | Tracked median | Paired tracked-minus-fixed range |
| --- | ---: | ---: | ---: |
| First completed/fsynced part | 48.03 ms | 48.52 ms | -33.73 to +12.15 ms |
| Output-inclusive completion | 208.71 ms | 179.12 ms | -33.52 to +9.80 ms |

The spread and warming across runs prevent interpreting the lower tracked
median as a speedup or a general overhead bound. The preceding NFS-output
audit (`productive_lifecycle_v1_20260922`) was noisier and is retained as such.
Neither audit evaluates a production-scale planner or public startup cost.

The local-storage transition run retained all 4,097 variants exactly once in
18 parts. Its reservations equaled the delivered source intervals, including
prefetched reads. Two written-block callbacks changed future sizes from 128
to 256 and then 512. The maximum absolute difference from the fixed run was
0.0001220703125. These two callbacks used explicitly scripted forecast values;
their cheap evaluation time must not be reported as analytical planner cost.

Source-generated single-option proposals are now implemented and checked in
jit_proposal_20260922.md. Their initial sufficient-improvement rule uses
conditional model bounds rather than an invented live queue checkpoint.
Public startup selection, automatic candidate choice and cost forecasts,
observation/parameter updates and safe phenotype/GPU reassignment are still
unfinished. JAGWAS continues to forbid phenotype partitioning. Immutable
measurement records and their original observation ages remain governed by
the calibration cache; productive execution snapshots are job-specific.
