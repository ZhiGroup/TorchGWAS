# Public tuning after useful output

`run_linear_gwas(..., autotune_profile=..., autotune_config=...)` now accepts an
optional `initial_chunks` object in the detailed configuration. Without it the
existing upfront bounded search is unchanged. This is an experimental bridge to
the analytical calculator, not a qualified automatic performance policy.
The [PGEN index-reuse follow-up](jit_index_reuse_20260922.md) avoids repeated
full-index parses between source opening, admission and productive windows;
it also reports the remaining whole-file metadata costs on a large input.

The starting context, GPU layout and phenotype width are explicit. Metadata-only
admission checks memory for every allowed chunk size and tail before execution;
there is no upfront timing search or PGEN payload census. The first completed
writer event starts a bounded planning window. Each callback can compare one
alternative chunk size using three bounded source horizons. Changes apply only
to future reads, with an exact reserved prefix on each variant/phenotype
partition. Already issued work and all association output are preserved.

Dense output supports fixed phenotype tiles or variant shards. Significant-pair
output uses fixed phenotype tiles and the host selector. JAGWAS uses variant
shards and keeps the complete phenotype panel and factor on every participating
GPU. The new path does not reassign GPUs or change phenotype tile widths during
a scan. The existing upfront calculator still handles layout selection.

## Immutable measurements and freshness

The detailed profile binds declared calculator coefficients to canonical,
immutable component records. Reuse keeps their hashes, observation times and
expiry times. Reading or republishing a profile never renews a measurement.
Source/runtime identity, empirical age and current available memory are checked
separately. New measurements are new records; an old profile is never silently
redirected to a newer coefficient.
The generic cache also [rechecks freshness after reading saved records](cache_read_freshness_20260922.md),
so time spent in cache I/O cannot make an expired observation reusable.
The [NumPy context binding](numpy_context_binding_20260922.md) also distinguishes
same-version binary builds and enabled CPU features, while reusing unchanged
file digests under fresh metadata checks.
The [NUMA policy binding](numa_context_binding_20260922.md) also records the
calling thread's allocation policy, allowed memory nodes and automatic
balancing setting; matching CPU affinity alone is insufficient.

Admission requires valid supplied prices. Before and after each in-job forecast,
the bridge rechecks source/context, original measurement age, exact artifact
identity, PGEN identity and current memory availability. A failure during an
optional planning step retains the current chunk size and allows the scientific
scan to continue. Invalid startup evidence fails before scanning.

Optional [public resident-copy refresh](public_refresh_20260922.md) now schedules
small probe batches in the initial useful chunks. It can reuse a consistent
original record or publish a new stable window and derived profile. Only named
copy coefficients may await this check; other resource probes and fully measured
production profiles remain required work.
Verifying declared fields does not certify unlisted coefficients or runtime
prediction accuracy.

## Configuration

Keep the existing `bounds`, `joint`, `qc_trait_block` and mode-specific component
artifact fields. Add `initial_chunks` with the following explicit fields:

| Field | Meaning |
| --- | --- |
| `context` | One named context in the bound profile |
| `chunk_size` | Starting size from the aligned `bounds.chunks` |
| `partition_axis` | `trait` or `variant`, as allowed by the output mode |
| `trait_block` | One allowed trait width, or `null` for variant sharding |
| `window_markers` | Three increasing capacity-aligned horizons, at most 32 smallest chunks |
| `budget` | `max_steps`, `max_cpu_seconds`, `max_window_seconds` |
| `cost_forecasts` | Remaining time and explicit CPU, wall, switching, publication and reserve costs |
| `forecast_options` | Boundary intervals, relative model error, slope tolerance, extrapolation limit and assumptions |
| `structural_cache_dir` | Optional source-bound structural work cache |
| `resident_copy_refresh` | Optional bounded check/refresh of declared CPU copy coefficients |
| `capacity_scenarios` | Optional named fractions of the bound shared CPU, DRAM, input and output capacities; include one all-ones nominal case |
| `background_planning` | Optional boolean; run a proposal on one planner thread and reject it if the source or written-output frontier advances |

The complete audit example is
[`direct_public_initial_chunks_20260922.py`](../benchmarks/direct_public_initial_chunks_20260922.py).
Its synthetic prices are test controls, not reusable hardware calibration.
`plan_cache_dir` is incompatible with this deferred mode because it would imply
an upfront ranking. Dense output requires an explicit `sumstats_block_bytes`.
Profiles must explicitly bind `return_beta` to the selected transport layout.

At most four declared host-sharing/output-occupancy/capacity combinations are evaluated.
Each capacity fraction must be positive and no greater than one. The fractions
reduce shared graph capacity while keeping component prices, per-device transfer
rates and the profile digest fixed. They are conditional availability assumptions,
not measured live capacity or a correction to loaded runtime. Adding scenarios
increases the synchronous first-chunk calculation cost, which is charged before
any switch. Omitting this field evaluates the nominal capacities only.
Reduced output occupancy is explicitly `empty` or `dense`; it is not inferred by
changing a scientific threshold. Each candidate must repay cumulative measured
planning wall time plus declared switch/publication/reserve costs across all
scenarios. Savings are now [compared within matched scenarios](matched_productive_scenarios_20260922.md),
using the smallest paired saving rather than crossing conditions between
baseline and candidate. Unstable source-window forecasts do not change the size. Planning
currently requires at most eight total partitions and sufficient unissued
source for all three horizons; otherwise it stops while scanning continues.

CPU and wall budgets are checked around each step; they cannot preempt a single
model calculation. By default, forecasting briefly holds new source
reservations while previously issued work continues. With
`background_planning=true`, the writer callback launches one worker and source
reservations continue. An advanced issue or output frontier invalidates the
proposal; the same candidate can be retried at a later completed output, with
both attempts charged. Finalization joins the worker to record its full cost.
See [released-frontier planning](jit_background_planning_20260922.md).
After planning stops, callbacks maintain counters without copying the full
prefix history.

## Output and evidence

`run.json` contains `autotune.productive`: exact bounded reserved prefixes,
complete source cursors, written-output counters, first useful output,
decisions and charged costs, scenario forecasts and original price evidence.
Dense writer progress is successful beta/t writes, not complete-store durability;
indexed part-file fsync is also distinct from final manifest publication.
Output-inclusive API time is reported separately in the audit benchmark.

The first real two-GPU audit used N=2,049, M=4,097 and K=512, native hard-call
PGEN and JAGWAS on server-local `/data`. Control, deferred and same-process reuse
outputs were bit-identical. Both deferred decisions completed model evaluation
and retained size 128 because the candidate marginal costs were unstable.
Planning cost 1.059 s and 0.860 s. The original measured-copy record hash and
observation timestamp were unchanged; its reported age increased. All 136
package source hashes matched the executed report.

Evidence: `results/public_initial_chunks_20260922/execution/report.json`.
The resident NumPy copy coefficient was independently measured; most other
prices were synthetic controls. Different warm states and lack of matched
repetitions preclude any speedup claim. This audit does not qualify calculator
accuracy or profitable tuning on a large production job.

The follow-up audit records preparation time outside the API boundary and first
written output relative to each API call. A separate Python process then loads
the saved profile without collecting another coefficient. It preserves every
bound artifact byte and all original observation/creation/expiry timestamps.
The copy record's age increases from 8.960 s in the previous call to 41.349 s in
the new process. All 4,097 joint outputs remain bit-identical to the control;
the new process also retains chunk size 128 after an unstable forecast.

Evidence: `results/public_initial_chunks_20260922/execution_v2/report.json` and
`results/public_initial_chunks_20260922/separate_process_reuse/report.json`.
The follow-up records 6.735 s of profile/probe preparation, 1.786/1.819 s of
in-process planning and 0.586 s of planning in the separate process. API times
include all output and differ with process/CUDA warm state; they are diagnostic
costs rather than matched performance comparisons.

The regression batch passed 263 tests; its two failures were test expectations
that omitted the dense writer's existing exception wrapper. The corrected
cleanup cases passed all four output routes. Another 80 tests cover public
observation collection, tiled output and immutable price binding. There were no
package source changes between the measured execution and these checks.

Two additional accounting tests execute the complete deferred dense and
significant-pair forecast path with captured tensor launch geometry, real PGEN
header windows and immutable bindings. Their hardware rates are synthetic and
their completed-output event is simulated; the executor tests above separately
exercise actual writes. Together the targeted runs pass 347 distinct checks.
Final logs include `tests_cleanup.log`, `tests_observation_regression.log`,
`tests_real_models_final.log` (significant case) and `tests_dense_model.log`
(corrected dense fixture), alongside the regression batch `tests_final.log`.
