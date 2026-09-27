# Bounded forecasts at the productive source frontier

`ProductiveTuningRun.forecast_step` now connects bounded prepared-window
comparisons to actual unissued source work. The callback runs only after useful
writer completion, within the existing held-frontier planning step. Header
construction, model comparisons and forecast arithmetic share that step's CPU
and wall accounting. Existing prefetched reads continue; new reservations wait
until the calculation returns. This is an internal integration path, not the
public automatic tuning policy.

## Exact execution binding

`productive_window_proposal` checks the complete retained reservation history:
partition identities and nonoverlapping association rectangles, contiguous
issued ranges, admitted chunk/tail widths, reader ordinals, source cursors and
the issued revision. Work remaining is derived from those cursors, not from
the number of written rows. Prefetched or submitted reads have already left
the unissued extent even when their output has not been written.

Both compared layouts must include every unfinished partition exactly once,
starting at its held unissued cursor. A finished source partition is excluded
from the unissued total without claiming that its GPU or output work has
finished. Fixed device/phenotype ownership, input dimensions, scientific data,
execution settings and component prices are preserved. The baseline uses the
current chunk size; the candidate uses one common admitted size. The comparison
contract includes a digest of profile fields other than chunk size and the two
duration-free kernel-geometry banks, allowing shape-specific geometry while
detecting changed reader allocations, depth, math settings or prices.

The adapter currently permits at most eight explicit partitions and 128
retained issued ranges by default. Excessive, incomplete or stale inputs fail
optional planning while the valid scan continues. It does not silently drop
remaining tiles to meet these limits. JAGWAS requires variant partitioning and
the complete phenotype panel on every device.

## Unequal remaining extents

Different devices can reserve different numbers of chunks before the callback.
A single remaining-pair count would hide this imbalance. The window forecaster
therefore also accepts exact remaining marker counts for each bound window.
Let h_i be the largest modeled marker span, r_i the actual unissued span, and
W3 the largest horizon's pair count. Its balanced workload scenarios are:

    core     = floor(min_i(W3 * r_i / h_i))
    envelope = ceil(max_i(W3 * r_i / h_i))

The lower time scenario extends the model to the core; the upper time scenario
extends it to the envelope, retaining one fill/drain term and the existing
slope/error/boundary adjustments. Every r_i must cover its modeled horizon,
and the phenotype-weighted sum must equal the exact remaining pair count.
The extrapolation limit applies to the envelope, not just the average ratio.
Dense periodic-writeback checks also use that upper extent.

These are explicit balanced workload scenarios around unequal remaining work.
They are not a proof that fluid graph duration is monotone under work removal,
or hardware-time confidence bounds. Future source work, survivor occupancy,
available resources and imbalance assumptions remain visible caller inputs.
The held cursor is an observed boundary; the GPU/queue/writer continuation
costs are still explicit supplied adjustments, not inferred live checkpoints.

## Applying a result

A stable forecast supplies baseline lower and candidate upper completion
scenarios to the existing productive cost gate. It can change only future
chunk reservations, after repaying cumulative planning, switching, publication
and reserve costs. Unstable slopes or an unseen dense writer regime yield no
change. A delayed, over-budget or invalid calculation also leaves the current
size in place. The run retains a compact audit containing the issued revision,
exact remaining rectangles, scenario extents, forecasts and comparison-contract
hash. It does not retain the three full model reports in every run snapshot.

Memory admission, source/runtime dependency checks, price age/qualification,
actual writer state, future output occupancy and candidate selection still
belong to the caller. The adapter preserves `selection_validated=False`.
It closes a source-position/accounting gap; it does not certify the calculator
for public automatic decisions.

## Verification

Remote job `20260922-071734-956112` passed **168 tests in 40.18 seconds**.
The log is `results/productive_forecast_v2_20260922/tests.log`. Tests cover
nonzero source starts, unequal device progress, prefetch, full source coverage,
unchanged prior reservations, stale cursors, missing/overlapping partitions,
reader ordinals, altered prices/data/settings, unadmitted sizes, source tails,
limits, profitable/unprofitable scenarios, unstable forecasts and the full-panel
JAGWAS rule. A controller test applies a mathematically specified profitable
scenario and verifies that only later read reservations change size. These
synthetic timing controls test decision wiring, not prediction accuracy.
Job `20260922-072211-957663` additionally passed the two newly added
source-generated contract tests in 12.86 s (`profile_binding.log` in the same
directory). They verify that actual comparison construction records changed
decoder allocation and writer price bindings, and that the adapter rejects
those changes before considering a chunk switch. Total: 170 targeted cases.

The same job completed a real two-GPU native-PGEN JAGWAS audit:
`results/productive_forecast_execution_v1_20260922/report.json`. It used 2,049
samples, 4,097 variants, 512 phenotypes and two covariates. The memory-admitted
chunk choices were 128/256, with capacity 256 and initial size 128. Structural
admission built one fixed layout from metadata, with zero payload census passes
and zero runtime candidates evaluated. Source setup, geometry preparation and
admission took 2.947 s and are reported separately from execution.

After the first fsynced result part, the callback inspected this actual held
frontier:

| Device | Complete source partition | Already issued | Unissued variants |
| --- | --- | --- | --- |
| CUDA 0 | [0, 2304) | Four 128-variant chunks | [512, 2304), 1,792 variants |
| CUDA 1 | [2304, 4097) | Four 128-variant chunks | [2816, 4097), 1,281 variants |

Thus 1,024 reserved variants are excluded from the unissued pair total even
though only one part has been written. Remaining work is 3,073 × 512 =
1,573,376 pairs. The three modeled horizons contain 256/512/768 variants per
device, start at those exact unissued cursors, and retain the actual reader
ordinals. Header construction occurs inside the charged callback.

The balanced core/envelope extents are 1,311,744 and 1,835,008 pairs. Their
ratios to the largest modeled horizon are 1.668 and 2.333, showing why using
only the mean remaining-pair ratio would hide the different shard lengths.
The current-size marginal slope changed by 1.19%; the candidate's changed by
21.90%, exceeding the declared 10% stability limit. The controller therefore
kept size 128. No threshold or budget was loosened after observing the refusal.

The complete planning callback cost **1.289 s CPU / 1.314 s wall**. It includes
header construction, all three source/resource comparisons, exact-frontier
binding and forecast arithmetic. It used an explicit audit-only 10 s CPU / 30 s
window budget and a generous 10 s remaining-time precheck; these are not the
production defaults or measured remaining-runtime forecasts. Boundary costs
were an explicit 0–1 s scenario, and component prices were synthetic controls
with captured duration-free kernel geometry. The experiment does not qualify
prediction accuracy, future occupancy, boundary costs or tuning profitability.

Both executions wrote all 4,097 variants exactly once in 33 parts. Reserved
ranges equal delivered ranges, all worker threads closed, and the maximum
statistic difference from the fixed control was zero. The fixed/forecast
execution times were 2.166/1.844 s; first output was 1.832/0.086 s. Their order
and different warm state prevent a speedup comparison. Those times include
output and exclude shared preparation/admission and process imports.

Input and output mounts were verified as XFS `/dev/md0` on `/data`. All five
input hashes stayed unchanged. The pulled report matches all 133 package-source
hashes, eight helper hashes and the benchmark hash. Calculator cost and real
output-inclusive prediction qualification remain open; this audit establishes
the live connection and correct refusal, not automatic tuning readiness.

The subsequent [calculator reuse audit](forecast_calculation_reuse_20260922.md)
profiles this calculation, reduces repeated structural work, and verifies the
updated callback on another actual two-GPU execution. Its matched replay
retains exact model reports while reporting total cost and collection pauses.

The [price-binding integration](price_binding_20260922.md) adds optional original
measurement-age checks before calculation and before a productive proposal may
change future chunks. Declared values are linked to immutable records; unlisted
prices remain explicitly unqualified.
