# Calculator, autotuning and multi-GPU handoff

Updated 2026-09-23 after compact source, writer and selected-output work.
This is a continuation guide, not a completion claim.

## Current status

The implementation has an analytical calculator, bounded candidate construction,
public execution bridges, immutable cross-job measurement reuse, and an
experimental controller that changes chunk size after useful output. Dense,
significant-pair and full-panel JAGWAS routes are connected and have correctness
and lifecycle tests. The calculator is **not yet qualified for reliable large-job
runtime prediction or profitable production autotuning**.

The central unresolved result is a frozen H100 panel in which the chosen
candidate was the observed winner, but its predicted executor time was
15.862 s versus an observed median of 55.769 s. Later observations vary widely.
Neither a whole-job correction multiplier nor a favorable component median has
been substituted for missing resource behavior.

The goal tracker is **active**. The user has clarified that the current cold
upfront planner is too slow and wants progress toward just-in-time autotuning.
The statistics-submission diagnostic is complete. No goal has been marked
complete.

The productive controller now has optional bounded GPU-stage observations and
an immutable, source/config/profile-bound history of actual planner CPU/wall
cost. A later job can use the latter only to raise its expected planning cost;
the record cannot replace component prices or approve a switch. Both features
were exercised against real two-GPU JAGWAS output. See
[stage observations](jit_stage_observations_20260922.md) and
[planner-cost history](jit_planning_cost_history_20260922.md). The observed
productive planner still ranges from roughly 0.4 to 2.5 seconds in these
shared-server controls; avoiding a long output-thread callback remains open.
The [planner phase audit](jit_planner_phase_audit_20260922.md) now isolates
source-window, graph-comparison and live-validation costs, and a bounded
per-header source-work cache removes repeated exact PGEN chunk bounds. Matched
header replay showed only a modest 54.8-to-46.9 ms improvement, so the cache
alone does not resolve the callback delay.

The [passive frozen H100 capacity check](frozen_passive_capacity_20260923.md)
completed three further fresh-process scans on the original GPU pair. Every
run passed the same 140 reference checks and cold local-`/data` controls.
Executor times were 22.135, 30.210 and 25.442 s versus the unchanged nominal
15.862 s prediction. The slowest run had the least aggregate observed
runnable wait but more decoder and selector CPU; a scheduler-wait-only or
global timing correction is unsupported. These observations use the older
frozen selector, and the current development selector has a different source
protocol. The diagnostic script and complete reports are preserved in the
separate H100 diagnostic checkout. The model still needs independent current
component prices, loaded first-chunk validation and a complete finite
output-inclusive continuation before an actual JIT switch.

The [H100 decode–selector coupling control](decode_selector_coupling_20260922.md)
completed nine randomized, independent two-second conditions on the frozen
source and real PGEN. Mixed selector CPU service rose by roughly 5–17%, while
throughput fell more sharply and runnable scheduler wait varied even in
selector-only controls. This does not reproduce the worst loaded selector cost
or justify a price correction; it strengthens the case for passive scheduler
and live-capacity observations during a matched full executor scan. GPU 1 was
occupied at the time, so no conflicting original GPU 1+2 rerun was launched.

The [JIT reader CPU/wait audit](jit_reader_cpu_wait_20260922.md) now records
reader-thread CPU, nullable scheduler runnable wait, and probe bracketing wall
time only on reserved productive chunks. Matched cached windows compare the
new CPU field and refresh on CPU-only drift or a missing old protocol field.
The final public two-GPU JAGWAS audit had equal control/deferred/reuse output,
two samples per GPU, and no pending measurements; scheduler wait was ambiguous
on A100 rather than certified zero. This adds live falsification evidence, not
an independent CPU-capacity price or a production tuning qualification.

The [downstream CPU/wait extension](jit_consumer_cpu_wait_20260923.md) samples
the caller thread only while a reserved result is yielded. It distinguishes
CPU work from nullable scheduler wait without changing the existing consumer
wall boundary. Separate probe-wall cost explains why wait can exceed that
narrower consumer interval under load; abandoned yields have null diagnostics.
Its scope includes output-mode-dependent downstream work, not necessarily
durable writing. The immutable per-device calibration key is now v3, and
cache validation rejects old or missing protocol markers. The previous v2
15-run public native-PGEN lifecycle audit matched control output in five
execution modes and reused immutable records in later processes. Three
injected failures published nothing, with all 144 historical source hashes
verified. The current v3 boundary passed 27 targeted CPU/GPU tests and a
separate 15-run public lifecycle rerun with the same output/cache guarantees.
Completed probe-wall brackets reached about 25 ms across three chunks on a
loaded GPU, beyond a 10 ms illustrative reserve. The productive JIT cost gate
now charges the larger of declared reserve and known completed probe-wall
time; a 100-test remote controller/output batch passed on the resulting source.
No loaded interval has been promoted to an independent capacity price.

The [cold JIT diagnosis](jit_cold_planner_diagnosis_20260922.md) found that
fresh context captures repeatedly hash source/library bytes; a stat-bound
digest cache reduced the second productive check and later-job checks. The
controller now uses the already configured planning-cost history directory
for those digests when no other digest directory is supplied. The real
8.1M-variant PGEN header also parsed in 0.079 CPU seconds in a current control,
versus an earlier loaded 11.608-second observation; matched huge-page on/off
controls did not justify changing the parser. Graph and source-window callback
work, plus the H100 absolute prediction gap, remain open.

The [conditional capacity study](jit_capacity_scenarios_20260922.md) separated
shared availability from immutable window/tile service profiles. The first-chunk
controller can now compare optional named shared-resource fractions within its
four-scenario limit, including a required nominal case. This is a charged
sensitivity check, not an inferred live capacity. A read-only frozen H100
experiment gave 15.862 s at nominal CPU capacity and 48.964 s at one eighth
capacity, still below the 55.769 s observed median. It does not resolve the
model gap or qualify a production switch.

The [large-source census and extrapolation audit](jit_census_scale_20260922.md)
profiles the frozen 1.05M-variant preparation and identifies 831,330 tiny
per-record NumPy unique calls. A source-equivalent development change reduced
median CPU for a matched 65,536-variant exact census from 2.027 to 1.370 s;
the full candidate still must not run in cold JIT startup. The productive
controller now skips a graph calculation when its exact remaining/horizon
ratio already exceeds the configured extrapolation cap, while preserving
scheduled component refresh. This avoids wasted callback time but exposes a
remaining scale limitation: short-window extrapolation alone cannot tune an
8M-variant job early under a conservative ratio cap.

The [released-frontier JIT work](jit_background_planning_20260922.md) adds an
opt-in background calculator step after completed output. Source issuance and
later writer events continue while it runs; a changed source/output revision
invalidates the proposal before it can change chunk size. The candidate is
eligible for a later-output retry, with both attempts charged. The default
remains synchronous, and finalization joins any worker to retain a complete
audit. This does not solve the 8M-variant extrapolation gate or the frozen
H100 absolute model error. Read-only source diagnostics showed a cheap
whole-header interval, but the frozen price profile lacks possible 4/5-byte
varint rates; the exact payload census confirmed none occur in that file at
18.795 process CPU seconds. Do not silently assign zero prices to the header
possibilities or run this census before first output.

The [whole-source schedule bound](jit_source_schedule_20260922.md) now counts
every native PGEN chunk-start LD replay using only indexed metadata and cancels
shared primary records and replay starts in a paired chunk-size comparison.
On the real 8.086M-variant A100 source, 128-to-1,024 markers removed 13,572
replays and 57,979,837 indexed read bytes; the full 128-marker schedule took
2.219 s in one ordered diagnostic. A frozen H100 1.05M-variant header schedule
enclosed an exact payload census in 0.580 s versus 18.870 s for that census.
These are **source-work** bounds only. They have not been connected to a finite
GPU/output continuation model and do not lift the early-job extrapolation gate
or resolve the frozen H100 15.862-versus-55.769 s executor prediction gap.
A further frozen H100 paired control priced 1,024-to-4,096 LD decoder work at
at least 0.000179 s saved under its component profile, while omitting read,
control, selector, GPU, output and contention effects. The LD decoder term
alone is far too small to account for the 39.907 s gap; larger chunks were not
admitted as a frozen runtime candidate.
An [individual-chunk expansion study](jit_source_schedule_20260922.md#cost-of-individual-chunk-bounds)
then found the ordinary full 8.1M-variant header window costs 290 CPU seconds
with the default signature cache, or 77 seconds with a 16,384-entry cache.
The new bounded vectorized window matches individual chunk bounds and costs
7.483 CPU seconds on that source at 571,616 KiB peak process RSS. The compact
aggregate remains below one second when its signature calculations are warm.
A narrower-index revision conserved the same work at 516,776 KiB peak RSS;
its CPU timing was confounded by a loaded shared-server run and does not
establish a speed change.
The 369-test related A100 batch passed before that memory change; the final
source passed 63 targeted PGEN and dense/significant graph-equivalence tests.
The vectorized path is a post-output, memory-admitted building block; it is
not yet wired into a whole-job GPU/output graph or public JIT switch.
The [finite continuation design](jit_finite_continuation_design_20260923.md)
specifies the next bounded optimization: use compact whole-source work to
screen layouts, refine only plausible candidates in charged background
steps, and include GPU, transfer, mode-specific output, shared resources,
issued work and drain before a switch. It is a design, not an implementation
or a qualification of the H100 absolute prediction.
The [source-complete read/decode floor](jit_source_resource_floor_20260923.md)
now prices the entire indexed PGEN chunk schedule and composes fixed dense,
significant or full-panel JAGWAS partitions across GPUs without duplicating
shared CPU/DRAM/input capacity. It is a conditional necessary floor for the
source stage only. A 249-test A100 source/resource batch passed before the
final identity recheck; 19 current-source focused tests and an aligned LD
shard-boundary check passed afterward. On the real 8.1M source, the ordered
header, 1,024-marker schedule and floor cost 0.482, 0.903 and 0.185 process
CPU seconds respectively in the final current-source diagnostic. These
shared-server, ordered timings remain post-output/background work. The
diagnostic's primitive rates were synthetic and no pipeline runtime or
profitable switch was inferred.
`native_layout_output_floor` now adds compact dense, significant-pair and
JAGWAS result D2H/output-array payload loads for the same fixed partitions.
Selected-row occupancy is an explicit interval, never inferred from an empty
first chunk. Shared output/D2H capacity is counted once and per-device D2H
limits remain separate. The result is still only a necessary payload floor:
GPU/H2D, selector and writer service, storage commit, queues and final drain
must be added before it can guide a switch.
The later read-only 8.1M-marker compact-admission probe measured 1.114 CPU
seconds for the header and 0.139 for its compact 128-marker memory layout;
three fixed/shifted envelope pairs each cost under 0.003 CPU seconds. This
is an ordered metadata-only component result under uncontrolled shared-server
load, not the complete cold public startup or a speedup measurement.
The [compute/frontier extension](jit_compute_frontier_floor_20260923.md)
counts unpacked native PGEN H2D and matrix products under declared shared and
per-GPU capacity ceilings. A bounded variant/phenotype sweep proves that a
future tile or GPU layout covers the same unissued pairs after the first
written chunk, without reclaiming reserved reads. A compact partial envelope
now composes matching read/decode, H2D/matrix and D2H/output payload floors;
the upper endpoint is only an upper endpoint on this incomplete floor, not a
completion ceiling. The source composer now handles many phenotype tiles with
a sorted overlap check rather than a quadratic pass. The final unequal-cursor
batch passed 265 related remote tests after a 59-test many-tile check. The
controller still cannot relocate live executors or price complete
output-inclusive continuation.
The [compact dense-writer ledger](jit_dense_writer_compact_20260923.md) adds
constant-work exact counts for fixed-chunk beta/t and optional df output,
including staging copies, minimum writes and writeback requests. It can be
attached to the native output floor with explicit writer settings. A matched
8.1M-marker component probe found the same 8,312,511,828 payload bytes and
189,519 staging calls as the expanded ledger while avoiding its per-chunk
event list. The 407-test output/continuation batch passed. This reduces the
cost of a future post-output first pass; public cold planner time and
conditional completion accuracy remain unresolved.
The final 656-test combined source/frontier/compute/output/writer batch passed
as A100 job `20260922-235648-1253302`.

The [multi-GPU transfer extension](jit_multigpu_transfer_floor_20260923.md)
now shares the finite graph's overlapping H2D/D2H link schema and adds one
combined host-memory load for decoder traffic plus both DMA directions. A
71-test A100 related batch passed; fixture placements with shared versus
separate links produced different necessary floors. No measured link bank or
complete GPU/output continuation is yet bound to a JIT switch.

## User requirements to preserve

- Improve the calculator until it is suitable for tuning, then finish chunk and
  phenotype-tile selection and efficient multi-GPU execution. Improvements to
  torchGWAS itself are authorized when analysis exposes wasted work.
- Keep one source/resource-derived analytical calculator, using independent
  component measurements. Do not fit or interpolate a whole-GWAS timing grid.
- Model shared CPU, memory, transfer, storage, GPU and output work. Balancing
  stages is a throughput objective under constraints, not a reason to assume
  linear GPU scaling.
- Include dense output, significant pairs and JAGWAS in tuning costs. Preserve
  scientific thresholds and output semantics while comparing layouts.
- **JAGWAS never partitions phenotypes.** Every active GPU retains the complete
  phenotype panel and factor; parallelism is by variant shards. Reject panels
  that cannot fit the required full-panel configuration.
- Expensive measurements should be spread over the first useful chunks of one
  job. This does not mean the first several separate jobs.
- Persist measurements as immutable records for later jobs. Separate structural
  compatibility, empirical age/drift and live availability. Reading, checking
  or republishing a record must not renew its original observation timestamp.
- Charge planning, validation, probing, switching and publication. A proposed
  switch must repay those costs over remaining work.
- Loaded stage durations are diagnostic observations, not independent service
  capacities. Keep CPU work, elapsed spans, waits and durable output separate.
- Edit in WSL, push, execute remotely, then selectively pull outputs. Do not
  edit remote source or run scientific computation locally. No custom CUDA or
  Triton kernels; CPU C/C++ is permitted. No subagents are authorized.

The sibling `/home/x/work/GWAS_reproducibility` voxel job was inspected for scale
context earlier. Preserve its role as a scale reference; internal unpublished
scientific findings are not public benchmark evidence or a new task scope.

## Authoritative working copies

| Purpose | Local WSL | Remote |
| --- | --- | --- |
| Current development | `/home/x/work/torchGWAS-jagwas-dev` | `lab-a100:/data484_4/zxie3/torchGWAS-jagwas-dev` |
| Frozen comparison package | `/home/x/work/torchGWAS-calculator-h100` | `lab-h100:/data484_4/zxie3/torchGWAS-calculator-h100` |
| Editable diagnostics | `/home/x/work/torchGWAS-calculator-diagnostics-h100` | `lab-h100:/data484_4/zxie3/torchGWAS-calculator-diagnostics-h100` |

The main checkout has no `.git` directory. Its current package contains 144
`.py`/`.cpp` files. The frozen comparison contains 104 such package files. Do not
edit or push the frozen project; diagnostics import it read-only. Historical
source hashes, plans, prices and results must retain their original identities.

Windows can edit through
`\\wsl.localhost\Ubuntu\home\x\work\torchGWAS-jagwas-dev` and the equivalent
diagnostic UNC path. The original `C:\Users\X\Downloads\torchGWAS` working
directory is unavailable; use `C:\Windows` as the shell working directory.

The remote interpreter on both servers is
`/data4012/zxie3/anaconda3/envs/heart/bin/python`.
Main `.remote` supplies `PY` and
`PYTHONPATH=src:/data484_4/zxie3/torchGWAS1.1/.deps:tests`.
Add `:benchmarks` when the selected test/audit imports benchmark helpers.

Read the installed `remote-lab-proj` skill and main `AGENTS.md` before further
work. The skill is at `C:\Users\X\.codex\skills\remote-lab-proj\SKILL.md`.
The main sync allowlist excludes `AGENTS.md`, `.codex`, results and job state.

## Implemented behavior and its limits

### Analytical plans and execution

Bounded searches compose preparation, exact source chunks/tails, decode,
transfer, eager statistics, reduction, queueing, writing and final drain in the
shared execution graph. Candidate memory admission includes explicit reserves
and mode-specific output/workspace costs. Host-sharing and output-occupancy
scenarios are explicit; retained counts are not guessed by changing a threshold.

Relevant entry points include `detailed_autotune.py`, `trait_tiling_plan.py`,
`trait_tiling_model.py`, `significant_host_model.py`, `indexed_schedule.py`,
`execution_graph.py`, and the JAGWAS candidate/planner modules. Public execution
enters through `run_linear_gwas(..., autotune_profile=..., autotune_config=...)`.
See [public JAGWAS tuning](jagwas_public_autotune_20260922.md) and
[bounded JAGWAS construction](jagwas_bounded_space_20260922.md).

Multi-GPU correctness work includes forwarding the aggregate reader budget even
when only one GPU remains active, and serializing first-use CUDA linalg
initialization before concurrent JAGWAS factors. Other backend failures are no
longer mislabeled as non-positive-definite phenotype matrices. Dense output
peer cancellation also preserves the originating writer failure. These are
correctness fixes; they do not establish throughput scaling.

The dense bounded planner now retains compact admission rows for all feasible
priced layouts. A fresh live-memory check can choose an already evaluated
alternative if the planned layout no longer fits, and old caches without these
rows are rebuilt. See [dense live admission](dense_live_admission_20260922.md).

### Tuning during initial useful chunks

`productive_run.py`, `productive_forecast.py`, `run_calibration.py`, and the
measurement/initial-calibration plumbing implement the experimental bridge.
Its configuration is documented in
[public initial chunks](public_initial_chunks_20260922.md). The [index-reuse change](jit_index_reuse_20260922.md) carries the PGEN
header from source open into JIT admission. Current [deferred-base admission](deferred_pgen_admission_20260922.md)
computes exact fixed/shifted reader-memory maxima on the fine grid without
building a full non-LD base array on a cold job. Bounded productive windows
use local predecessor lookup after useful output. An optional, source-bound
v2 cache publishes the full base index only after successful output, then
reuses aligned vectors and dtype-matched lookups in later jobs. It charges
retained bytes and never renews empirical prices. The source reader still
parses the full PGEN header before the first output.

- Starting chunk size, device assignment and phenotype width are explicit.
  Startup uses metadata-only memory admission rather than a timing search or
  full PGEN payload census.
- The first completed writer event admits bounded planning. Changes apply only
  to unissued source ranges; exact reserved prefixes include prefetch work.
- `initial_chunks.background_planning=true` is an experimental opt-in that
  releases the issue frontier during calculator work and discards stale
  proposals. Its worker is joined before final audit, so a short run may pay
  its calculation at the end.
- The current bridge changes chunk size within a fixed layout. It does **not**
  retile phenotypes or reassign GPUs during a scan. Upfront bounded planning
  handles layout selection.
- Dense output supports fixed phenotype tiles or variant shards. Significant
  output uses fixed phenotype tiles with the host selector in this bridge.
  JAGWAS uses variant shards and the complete panel on every GPU.
- Current bounded comparisons use three increasing source horizons, at most
  four declared host/output/capacity scenarios, at most eight total partitions, and explicit CPU,
  wall, step and continuation limits. Unstable extrapolation retains the current
  chunk size.
- The corrected decision compares matching scenarios:
  `min_s(baseline_lower_s - candidate_upper_s) > cumulative_planning_wall + switching + publication + reserve`.
  Every scenario also has to satisfy the remaining-time horizon. See
  [matched scenario accounting](matched_productive_scenarios_20260922.md).

### Immutable reuse and first-chunk refresh

`calibration_cache.py`, `detailed_calibration.py`, `price_binding.py`, and
`binding_digests.py` bind exact records, source/runtime context, artifact bytes,
field values and original ages. Structural records have dependency validity;
empirical records have finite lifetimes. Live memory/contention is not a saved
capacity. New evidence creates new records and profiles rather than modifying
old evidence or silently redirecting an existing binding.

Recent fixes include a post-read freshness check, backwards-clock miss behavior,
NumPy binary/CPU-feature identity, and fresh NUMA policy/allowed-node/balancing
context. Cached file digests still require fresh metadata checks. NUMA collection
is read-only, Linux-specific and fails closed if required policy queries are
unavailable. It describes the calling thread, not every existing worker, VMA
policy or page residency. See [cache freshness](cache_read_freshness_20260922.md),
[NumPy context](numpy_context_binding_20260922.md), and
[NUMA context](numa_context_binding_20260922.md).

The public refresh implementation is deliberately limited to two named
resident-copy CPU fields: `owned_result_copy_scenario.resident_cpu_seconds_per_byte`
and `process_units.numpy_copy_bytes`. It does not refresh GPU, decode, transfer,
storage, allocator or reduction capacities yet.

`resident_copy_refresh.py` and `cpu_service_refresh.py` use two consistency
samples for reuse, or seven samples over four callbacks (2+2+2+1) for a new
window. Drift/spread/temporal checks are declared heuristics, not confidence
bounds. Buffers, checks, I/O and publication are charged, and unfinished or
rejected windows are not published. See [public refresh](public_refresh_20260922.md).

## Verification that matters

The following are distinct evidence scopes; do not combine their test counts
into a claim that the entire current package was rerun as one suite.

| Evidence | Result | Scope |
| --- | --- | --- |
| Background JIT A100 jobs `20260922-211823-1210235`, `20260922-212156-1211189` | 119 controller/run tests and 107 related forecast/binding/history tests passed | Released-frontier staleness, retry and default synchronous compatibility; distinct targeted batches |
| Public background JAGWAS job `20260922-212005-1210462` | Three equal 4,097-row outputs, zero maximum difference | Two background proposals became stale and retained chunk size 128; mostly synthetic prices, no speedup or switch qualification |
| Census/JIT scale job `20260922-205915-1207707` | 319 targeted tests passed | Exact census/range, candidate, bounds, first-chunk gate and issue frontier; legacy `test_pgen.py` was uncollectable due an absent benchmark module |
| Conditional JIT capacity jobs `20260922-204359-1205016`, `20260922-204559-1205351` | 215 focused tests and 3 added real-model targeted tests passed | Shared-capacity sensitivity and controller pairing; synthetic fixture prices, no production switch qualification |
| NUMA/current binding job `20260922-153744-1098179` | 143 tests passed, 25.94 s | Latest targeted package validation; includes cache, detailed/NumPy/NUMA context, price binding, refresh and productive digests |
| Live NUMA binding audit, same job | Modes 2 and 8194 rejected the original saved profile; restoration accepted it without byte changes | Synthetic coefficient tests compatibility, not performance |
| Matched scenario job `20260922-142826-1079293` | 255 tests passed, 86.72 s | Includes 36 new dense/significant/JAGWAS scenario controls and existing public execution/forecast checks |
| Post-read cache freshness job `20260922-142135-1078815` | 267 tests passed, 72.76 s | Expiry during reads, older valid fallback, clock rollback, immutable bytes and integration |
| Public deferred JAGWAS execution | Bit-identical outputs; separate process reused original record/time | Most supplied prices were synthetic; selected chunk size remained 128 |
| Public refresh five-run audit | 16,385 joint rows per run, zero numerical difference; original records preserved through reuse/drift/expiry | Only resident-copy price was independently measured; no speedup established |
| JIT index-reuse split A100 regression | 20 native-reader tests and 88 subtests, plus 180 admission/planner/output tests passed | Includes retained-header reuse and stale-file rejection; separate full-suite qualification remains |
| Public JIT index-reuse execution | Bit-identical control/deferred/reuse JAGWAS output; 0.721/1.449 s proposal steps | Synthetic rates and varying warm state; no production speedup claim |
| Certified JIT-window and vectorized-layout A100 batch | 199 tests passed in 56.02 s | Includes exact per-chunk LD replay extents and stale-file rejection |
| Final public JIT execution, `header_reuse_execution_v2` | Bit-identical control/deferred/reuse JAGWAS output; no switch | 2.832/1.496 s proposal steps under varying load, mostly synthetic rates; no speedup claim |
| Retained-byte public JIT execution, `header_reuse_execution_v3` | Bit-identical control/deferred/reuse JAGWAS output; exact 36-byte fixture base-index charge | 2.100/0.547 s proposal steps; highly variable shared-server timing and mostly synthetic rates |
| JIT index-reuse A100 batch | 200 tests passed in 60.28 s | Prior source, including 32-bit retained-base accounting and live credit |
| Compact JIT admission differential | 48 tests passed; public two-GPU JAGWAS audit bit-identical | Full and compact fixed/shifted memory envelopes agree; audit prices mostly synthetic |
| Large compact envelope comparison | Whole-file and two-shard read/scratch maxima agree exactly on 8,086,101 variants | Ordered component timing only; 2,021,544 bytes of vectors, initial index parse remains expensive |
| Current compact JIT A100 batch | 229 tests passed in 64.31 s | Full and compact memory, header, controller, candidate and output regression |
| Source-bound admission-cache A100 batch | 235 tests passed in 67.99 s | Exact/corrupt/stale cache behavior, config validation, unsupported-type rejection, short-remaining-work gate and admission/controller regression |
| Public admission-cache JAGWAS audit | Control/deferred/reuse output equal, zero maximum difference; first miss stored, second hit | Original component-price observation time preserved; mostly synthetic prices, no runtime speedup claim |
| Large real PGEN structural cache | 8,086,101 variants; ordered earlier build 4.916 CPU s versus load 0.104 CPU s | Historical v1 component evidence; header parse remains separate |
| Deferred-base A100 regression | 238 tests passed in 58.38 s | Final cold/hot cache, source-window, admission/controller and output-mode tests |
| Exact cold-admission comparison | Deferred 0.066/0.890 CPU s versus full-base 3.290/8.069 in one ordered 8.1M-variant run | Same prepared header and identical envelopes; source header was another 11.608 CPU s |
| Aligned v2 cache on large PGEN | After-job publication 5.613 CPU s; hot load 0.091 CPU s; sampled LD predecessors equal | Component costs under shared load; publication occurs after useful output |
| Public deferred-base JAGWAS execution | Control/cold/reuse output equal; cold retained no base, published; reuse hit | Original component-price date preserved, mostly synthetic rates, no whole-job speedup claim |
| Public no-cache JAGWAS execution | Control and both deferred passes matched output; no cold retained base | Default path correctness, mostly synthetic rates, no speedup claim |

Latest pulled structural evidence is in
`results/deferred_pgen_handoff_20260923/`: the large-file and aligned-lookup
reports, public two-GPU JAGWAS audit, parser comparison and six matching
package source hashes. Earlier NUMA artifacts are under
`results/numa_context_binding_20260922/`; their source and test counts are
historical and must not be represented as validation of the current package.

Public lifecycle evidence is under `results/public_initial_chunks_20260922/`
and `results/public_refresh_20260922/execution_v3/`. These verify useful-output
callbacks, immutable reuse and correctness. Their largely synthetic profiles,
unchanged selected sizes and varying warm/load conditions do not qualify
profitable tuning. The refresh audit charged 2.330-5.357 s of planning wall time,
despite only 0.028-0.119 s of copy-probe CPU work; cost accounting must remain.

## Performance investigation: what has and has not worked

The frozen workload is N=35,365, M=1,048,576, K=16,385, C=27, native hard-call
PGEN, chunk 1,024, trait width 8,193, GPUs 1+2, 4.5 GiB/GPU allocator cap and
empty significant-pair output. Genotypes are a real prefix; phenotypes are
synthetic null inputs. Input/metadata are on verified local `/data`. Full scans
use cold-input and read-scan-read controls and 140 sparse reference cells with
exact df. They do not compare every dense statistic.

Frozen plan: `torchGWAS-calculator-h100/results/significant_bounded_ranking_20260921/plan.json`.
Frozen inputs: `/data/zxie3/torchgwas_bench/significant_host_sustained_20260921`.
The 12-candidate, five-repeat panel selected the observed winner; its 15.862 s
prediction still underestimated that candidate's 55.769 s median.

| Investigation | Main finding | Decision |
| --- | --- | --- |
| Loaded stage/profile controls | Excess variable CPU work in decode/selection; summed GPU compute spans much closer to priced work | Diagnose components; do not reprice from loaded intervals |
| Blocking CUDA completion events, three pairs | Result-worker CPU fell sharply, but elapsed-time direction reversed in the final pair | Default unchanged; no reliable speedup claim |
| Host serialization, nine unchanged-price scenarios | Predictions 15.862-17.995 s | Does not explain the observed median; not measured GIL parameters or guaranteed bounds |
| Per-statement selector sampling | Predicate/mask costs and some large recurring minor-fault counts | Fault counts alone do not prove fresh allocation or first touch |
| Private scratch reuse, three pairs | Correct outputs, but executor ratios reuse/current 1.128, 1.108, 1.173 | Prototype rejected; not in production |
| Independent NUMA selector control, 18 observations | 27,648 calls; roughly 6.4-7.0 CPU ms/call with one worker and 7.2-7.5 with two | No consistent binding benefit or reproduction of the large loaded penalty |
| Earlier fused CPU selected-output packing | Correct generic/public execution checks, no whole-pipeline benefit | Benchmark-only prototype |

Scratch reuse still had a later fault burst despite retained buffers. The earlier
fresh-allocation interpretation was corrected. NUMA experiments queried 918
sampled input pages, all on node 1, but do not certify all scratch/output pages.
Raw scheduler counters were sometimes already large at entry; absolute values
must not be described as migration performed by the current run.

Read [frozen diagnosis](frozen_runtime_diagnosis_20260922.md),
[selector stages and rejected reuse](selector_stage_diagnostic_20260922.md),
[NUMA controls](selector_numa_capacity_20260922.md), and
[CPU packing control](selector_pack_control_20260922.md) before repeating these
experiments. The current host selector is `bounded_flat_v2`; the frozen package
uses `bounded_flat_v1`. NumPy remains the default, with a separate explicit
native CPU path. Never change source identities to reuse an old price bank.

## Completed statistics-submission control

Project: `/home/x/work/torchGWAS-calculator-diagnostics-h100` on `lab-h100`.

- Job: `20260922-161156-1107371`.
- Wrapper PID: `3169861`; Python parent PID: `3169889`.
- Log: `.proj/logs/20260922-161156-1107371.log`.
- Script: `qualify_statistics_submission.py`.
- Probe SHA-256: `aafbfdafefd3eaac052f63299e832aee6e2ddf4d4c1207763fe44fc2fdf7683a`.
- Results: `results/statistics_submission_20260922/`.
- The job wrote all 18 observation reports, `complete.json` and `summary.json`.
  Do not rerun it into the same result path; writers use exclusive creation.

The schedule contains three complete randomized repeats of tiny/large source
shapes and GPU 2 alone, GPU 3 alone, and GPUs 2+3 concurrently. Each worker has
eight warmups and 128 measured calls, with a blocking gate every four calls.
CPU submission work, wall spans, waits, allocator retries and event spans are
separate. CUDA event spans include launch gaps and are not kernel-service sums.

The two-case smoke job `20260922-160513-1106587` is terminal. It passed stable
selected-output and exact-df checks on both GPUs. Large-shape submission CPU was
0.671-0.678 ms/call; allocated/reserved peaks were 2,060,235,776/2,367,684,608
bytes per GPU, with zero retry/OOM counters. Only four measured calls per GPU
were used, so no transfer/scaling conclusion follows yet. All 104 frozen package
hashes and the probe hash matched the pulled smoke report.

GPU 1 was occupied by another job; GPU 2 and GPU 3 were idle at launch. The
complete three-repeat control found large-shape isolated submission CPU near
0.60/0.64 ms per call and concurrent CPU near 0.64/0.63 ms per call on GPUs
2/3. There was no consistent large CPU penalty. This generic control therefore
does not explain the 15.862-versus-55.769-second full-executor discrepancy or
qualify the original GPU1+GPU2 plan. Decode, result DMA, host selection, output
and full-pipeline allocation pressure were excluded. All 104 frozen package
hashes and both diagnostic script hashes matched after selective pull. Read
[submission control](statistics_submission_control_20260922.md).

## Immediate continuation

The opt-in released-frontier controller has now been tested against actual
two-GPU JAGWAS output. Its 0.91/0.59-second proposals became stale while source
issue advanced and were correctly rejected. To achieve an early switch on a
large job, add a source-complete continuation that can validate priced
possible decoder units, output regime and shared capacity without a cold
payload census or a thousandfold short-window extrapolation. A future
long-running public run must demonstrate a stable profitable decision and
output equality; this small fixture cannot.

The first-chunk controller now samples completed reduced-output counts and
rejects an empty-only or full-only scenario that the current job has already
refuted. A cheap issue-frontier gate skips price revalidation and model work
when too little unissued source remains. The bounded observation and 148-test
focused A100 verification are documented in
[JIT output feedback](jit_output_feedback_20260922.md). This does not yet
update independent prices or establish a profitable switch.

The JIT controller can now collect bounded native read/decode and CUDA stage
observations from the first productive partition on each GPU, with an explicit
instrumentation cost reserve charged to any chunk decision. The final 150-test
A100 batch and public two-GPU JAGWAS audit are in
[JIT stage observations](jit_stage_observations_20260922.md). These loaded
intervals remain diagnostic; no independent component rate was replaced.

1. Reduce or defer the first-job PGEN header parse before output. Exact cold
   admission now avoids the full base index and took 0.066-0.890 CPU seconds
   versus 3.290-8.069 for the old path in one ordered large-file comparison.
   The initial header still took 11.608 CPU seconds there and varied widely
   across observations. A tested global-assembly parser was not faster; any
   reader redesign must preserve exact file scope, native LD replay,
   source invalidation and memory bounds. The short-remaining-work gate
   already skips source-window construction while retaining live validation.
2. Extend productive independent measurement to the bottleneck services and test a real
   bounded switch after useful output. The current bridge changes chunk size
   only, with fixed phenotype/device assignment. Preserve the exact issued
   frontier, immutable price records and charged callback costs. Dense live-memory
   fallback passed its corrected A100 batch: 120 tests in 68.11 seconds.
3. Investigate coupled decoder/host-selection/result ownership and resource
   load for the H100 discrepancy. Do not replace independent capacities with
   the completed submission control's loaded wall spans.
4. If a source inefficiency is demonstrated, change current development code
   and its work/memory/schedule ledger together, test ownership and numerical
   correctness, then run matched remote controls. Preserve frozen evidence.

## Remaining work before completion

The post-output partial envelope now includes an optional dense-writer service
floor, with a compact exact work ledger and independent finite-graph price
fields. It conserves expanded graph CPU work in a differential test; the
latest 416-test A100 batch passed in job `20260923-001603-1256153`.
See [writer service](jit_dense_writer_service_floor_20260923.md). This is
still only a necessary floor, not a complete continuation or JIT switch rule.
The full-panel JAGWAS path additionally has a compact, source-bound
[host selector floor](jit_jagwas_selector_floor_20260923.md) for explicit
retained-count intervals. It keeps the single indexed consumer's serial CPU
work and shared CPU/DRAM constraints visible without expanding all chunks.
The related A100 batch passed 133 tests in job `20260923-002451-1261842`;
the final shared-DRAM and tail-focused batch passed 20 in
`20260923-002620-1262116`.
The [JAGWAS indexed archive floor](jit_jagwas_archive_floor_20260923.md)
extends that compact continuation with nonempty NPZ part framing, page-cache,
storage and fsync service, plus combined shared CPU/DRAM and one-consumer
loads. Its writer declaration is explicitly durable; qualifying a finite
completion ceiling and a profitable same-job switch remains necessary.
The broader A100 layout/JAGWAS batch passed 141 tests in
`20260923-003326-1262674`, followed by 19 focused tests in
`20260923-003441-1262899` after the sparse-part bound was tightened.
The significant-pair path now has a compact
[indexed archive floor](jit_significant_archive_floor_20260923.md), using the
same independent writer services as the detailed graph. Host selection allows
one part per nonempty source chunk; device selection allows one per nonempty
selection block. Device output also includes the blocking four-byte nonzero
count transfer for every block, even when no pair survives. The explicit
device cell limit binds the source selector's chunk/trait block geometry.
The host significant D2H floor now follows the public scan's unconditional
beta return for significance, including `sumstats_fields='t'`; archive
payload still follows the writer's beta omission.
The matching [host selector floor](jit_significant_host_selector_floor_20260923.md)
uses independent primitive prices and explicit empty/sparse/dense NumPy
regimes to count entire unissued trait tiles without a chunk graph. It is
combined with the indexed archive in shared CPU/DRAM loads, but producer
selection and single-writer archive remain distinct serial chains.
The A100 significant/layout regression passed 154 tests in
`20260923-004814-1263921`; the combined selector/archive follow-up passed
32 in `20260923-004947-1264219`. The focused host/device archive follow-up
passed 38 in `20260923-005617-1273488`. Neither establishes a loaded JIT gain.
The broader selected-output/controller batch passed 229 tests in
`20260923-010029-1274003` with the `benchmarks` primitive bank on `PYTHONPATH`.
After exposing the count bytes in the source ledger, 37 affected tests passed
in `20260923-010221-1277165`.
The compact [device count-transfer floor](jit_device_count_barrier_floor_20260923.md)
also binds the existing independent per-count latency price to the actual
selection-block census. Successive tiles on one GPU serialize; multiple GPUs
can overlap. It is a necessary service floor only; device selector kernels,
host dispatch, the finite queue and loaded calibration remain outstanding.
The corresponding A100 cross-mode regression passed 253 tests in
`20260923-010907-1283679`, including a comparison with the exact device
selector graph. This validates necessary-floor arithmetic, not switch quality.
The final focused validation passed 19 tests in `20260923-011119-1283981`.
The [shared device-selection geometry](jit_device_selection_geometry_20260923.md)
now minimizes blocking nonzero calls under the existing selection-cell cap,
and the calculator follows the same full/tail block shapes. A 512-by-600,000
empty or sparse selector uses 293 blocks instead of 512; matched A100
selector-only controls favored the new shape in every adjacent pair. The
80-test focused regression passed in `20260923-012009-1284623`. This is a
source-level optimization and accounting correction, not a measured full-job
gain or a qualified JIT tile choice.
The follow-up corrected device-selection memory admission to use that exact
shape and bound the geometry helper in saved selector source identities. A
fresh untimed A100 12-case launch census passed source and kernel audits; exact
nonzero allocation captures for the two 512-by-600,000 mask extents yielded a
138,489,856-byte selector-local bound. The [geometry note](jit_device_selection_geometry_20260923.md)
links the captured reports. The updated focused A100 batch passed 76 tests in
`20260923-013929-1307377`, and the final cross-mode regression passed 342
tests in `20260923-014111-1307765`. Old kernel/workspace censuses must be
recaptured for shapes they do not cover; changing a recorded hash is invalid.
The [sparse productive-output scenarios](jit_sparse_output_scenarios_20260923.md)
extend first-chunk comparisons beyond empty/dense endpoints, rebinding one
fine-grid survivor ledger across both chunk sizes. A partial first-chunk part
now rules out an empty/dense-only switch. Fractions and placements are declared
conditional scenarios, not learned future selectivity or fresh price records.
The final A100 productive-planner regression passed 230 tests in
`20260923-015657-1317408`, and the post-edit focused regression passed 12
tests in `20260923-020028-1317646`.
The [productive whole-source floor](jit_productive_source_floor_20260923.md)
now binds compact indexed PGEN read/decode work to the exact post-output
unissued frontier, including alternative chunk, tile and GPU geometry.
It now composes necessary source, H2D/GEMM and output payload floors for
dense, significant and JAGWAS scenarios, but does not switch the public
executor.
The focused A100 accounting regression passed 140 tests in
`20260923-020731-1321659`; it does not qualify service prices or a
profitable layout change.
The [global sparse output ledger](jit_whole_layout_sparse_20260923.md)
then makes conditional retained-pair and JAGWAS row counts invariant to
whole-source tile and GPU repartitioning. It feeds the existing reduced
output floor but does not infer future selectivity or qualify a JIT switch.
The global ledger and output regression passed 262 A100 tests in
`20260923-021817-1322713`; after short-window ledger reuse, 201 focused
tests passed in `20260923-022013-1325839`. The composed post-output floor
passed 134 tests in `20260923-022204-1326005`, followed by 27 focused
tests after the stale-frontier check in `20260923-022325-1326232`.
The [productive output boundary](jit_output_boundary_20260923.md) now
retains bounded per-partition completed output against reserved source
ranges, separates dense matrix and df prefixes, and counts issued pairs
without writer completion. Invalid bound events stop optional chunk-size
planning before model construction. Background binding is synchronized with
writer callbacks. A separate conditional backlog ledger covers the upper
array payload of that issued-but-unwritten work, using the same global
occupancy scenario as unissued work in significant and JAGWAS modes. It is
workload evidence for a future finite continuation, not an elapsed-time
completion bound or a public switch. `productive_partial_floor` can now
attach that synchronized backlog and sum it with unissued array payload
without changing its necessary resource floor. The boundary and real public
output paths passed 147 A100 tests in `20260923-022947-1326771`, 19 in
`20260923-023118-1326991`, and 21 after the backlog/df extension in
`20260923-024139-1329483`. The composed floor/backlog integration then
passed 26 focused tests in `20260923-024453-1329852`.
The [issued-work ledger](jit_issued_work_20260923.md) now inspects only
bounded original chunks without matrix or indexed-part completion and
retains their full conditional PGEN source, H2D and mandatory GEMM workload.
It has no elapsed-time ceiling, and is not yet joined to the finite schedule.
The dense/significant/JAGWAS accounting regression passed 28 A100 tests in
`20260923-025236-1330487`. A read-only server-local 8.1M-variant PGEN probe
counted 63 pending 128-marker chunks in 0.692 wall / 0.683 CPU seconds after
a separate 1.052-second header construction; see the issued-work note and
`results/issued_work_probe_20260923/report.json`. This remains background
calculator cost, not cold startup or executor throughput.
A second ordered read-only probe conserved all work counts with 0.773-second
ledger CPU and recorded exact input/code hashes in
`results/issued_work_probe_20260923/report_v2.json` (job
`20260923-025919-1331906`).
The [checkpoint workload join](jit_checkpoint_ledger_20260923.md) now binds
fixed issued and candidate unissued reports at one revision, counts shared
source/H2D work once, and retains separate old/new GPU assignment. It is a
workload join with explicit missing completion terms, not a time ceiling or
a production tile/GPU switch.
The composed checkpoint and public-output regression passed 30 A100 tests in
`20260923-031627-1337387`; the report also identifies absent mode-specific
selector/writer floors before attempting a completion model.
The [compact GPU shape service](jit_gpu_shape_service_20260923.md) now reuses
the exact tensor calculator for each distinct full/tail candidate shape after
output, with per-device kernel and shared host-dispatch totals in the
checkpoint. It also prices bounded issued chunks on their original GPUs.
Its nine focused A100 tests passed in `20260923-032955-1340677`, and the
preceding broader regression passed 42 tests in `20260923-032722-1340524`.
The final source/shape/checkpoint/mechanistic/JAGWAS regression passed 45
tests in `20260923-033407-1341205`; a real captured-geometry JAGWAS check
matches the compact and expanded tensor-service totals.
The [incremental PGEN source schedule](jit_incremental_source_schedule_20260923.md)
now spreads exact full-source metadata accounting over bounded steps, then
rebases it to the current cursor and a changed chunk or variant-shard span.
`productive_source_floor` accepts the complete staged ledger and reports its
prior planning cost separately. The 8.09M-variant read-only matched control
conserved 63,173 chunks and 15,475 LD restarts; a cold first staged pass was
still around 9 to 12 seconds, so this belongs in charged background work.
The PGEN/source/checkpoint regression passed 235 A100 tests in
`20260923-035322-1346070`. The public controller now has an opt-in,
evidence-only background driver: each completed output event earns one
bounded source step, with a declared host reserve and step/CPU/wall/ledger
limits. It leaves the starting chunk in place because its planning cost is
not yet included in a complete finite continuation. Direct-first A100
controls showed 1M-record staged work took 3.38 s with a large primitive
cache versus 10.87 s with the default cache; the larger cache raised peak
process RSS to about 266 MiB. The production evidence path retains the
default cache until memory and throughput are qualified. The final combined
regression is recorded in the linked incremental-source note.
The [dense writer queue observation](jit_dense_writer_queue_20260923.md)
adds stream-local accepted, staged, queued, active and fully written byte
counts to native completed-prefix events. The output boundary retains one
latest observation per writer directory and binds device-free native events
only when their full source/trait interval has a unique admitted owner.
This exposes output-backlog state without treating asynchronous stream
observations as an atomic execution checkpoint or completion bound.
This is conditional component service, not a calibrated completion ceiling;
the public JIT still does not use it to switch layouts.

- Qualify absolute prediction and selection over realistic workloads and
  resource conditions, including nonempty significant output and dense/JAGWAS
  output-inclusive runs. One successful candidate ranking is insufficient.
- Obtain complete independently measured profiles for the current source;
  lifecycle audits with synthetic coefficients are not production calibration.
- Replace the short-window large-job extrapolation gate with a finite
  continuation calculation. Use the header-only whole-source PGEN schedule for
  source work, then include per-chunk GPU, transfer and mode-specific output
  service, shared-capacity constraints, in-flight work and final drain. Do not
  interpret the one-sided decoder CPU-work bound as a pipeline-time saving.
- Expand bounded in-job measurement/refresh beyond resident-copy CPU service,
  with source-compatible evidence and full cost repayment. Validate an actual
  profitable decision, including preparation and output boundaries.
- Finish/revalidate huge-panel chunk/tile/device selection and multi-GPU shared
  resource behavior. In-job layout changes are not implemented by the current
  initial-chunk bridge. Keep JAGWAS full-panel restrictions.
- Resolve conservative in-job host-memory admission: it compares a full envelope
  with current `MemAvailable`, potentially double-counting this job's already
  resident allocations. Do not credit all RSS; explicit owned-allocation
  accounting is required. Global availability also does not certify a strict
  per-NUMA-node memory budget.
- Keep first-use initialization, unpriced Python/control work, allocator-route
  changes, live load, extrapolation assumptions and boundary costs explicit.
  Do not hide gaps behind a successful validator or an aggregate correction.

## Remote commands and sync precautions

From PowerShell, with working directory `C:\Windows`:

```powershell
$adapter = 'C:\Users\X\.codex\skills\remote-lab-proj\scripts\invoke-proj.ps1'
$diag = '/home/x/work/torchGWAS-calculator-diagnostics-h100'
& $adapter -Project $diag logs '20260922-161156-1107371' -n 30
& $adapter -Project $diag run 'ps -p 3169861,3169889 -o pid=,etimes=,stat=,args='
```

Use `proj run --bg` for long work; preserve returned job IDs/PIDs. Do not
hand-write background SSH launchers. Never kill an unrelated job. The main
`.rsyncignore` is 232 bytes; the diagnostic ignore is 9 bytes (`results/` plus
newline). Verify the expected ignore before every push. A typical push is:

```powershell
wsl.exe -d Ubuntu --cd /home/x/work/torchGWAS-jagwas-dev -- bash -lc 'export RSYNC_RSH="ssh -o ControlMaster=no -o ControlPath=none"; proj push'
```

Results are intentionally excluded by these sync rules. For a selective pull,
save the exact ignore bytes, temporarily blank the local ignore, pull only the
named result directory, and restore the bytes in `finally`. Do not pull the
whole remote root or leave the ignore blank. Do not push during that temporary
state. `proj push` has no deletion flag here; do not use `--mirror`.

Typical A100 targeted-test environment, for reuse only when appropriate:

```bash
export PYTHONPATH=src:/data484_4/zxie3/torchGWAS1.1/.deps:tests:benchmarks
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export OMP_WAIT_POLICY=PASSIVE GOMP_SPINCOUNT=0
export TORCHGWAS_PGEN_LIBRARY=/data484_4/zxie3/torchGWAS-jagwas-dev/.build-libs/libtorchgwas_pgen.so
export TORCHGWAS_PGEN_BACKEND=native TORCHGWAS_PGEN_PACKED=0
export TORCHGWAS_NATIVE_STATS=0 TORCHGWAS_SCAN_PROFILE=0 TORCHGWAS_BLOCKING_EVENTS=0
export TORCHGWAS_SIGNIFICANCE_BACKEND=host TORCHGWAS_HOST_PREDICATE=numpy
```

H100 diagnostics intentionally use a different frozen environment: OMP=4,
OpenBLAS/MKL=1, CPUs 12-19, NumPy huge-page advice off, passive OpenMP waits,
and TF32 off. Let their harness establish it; do not silently substitute the
A100 test settings or promote diagnostic measurements across contexts.

### 2026-09-23 dense writer continuation increment

The native dense writer now emits a compact queue observation with each
completed-prefix callback. Beta, t and optional df accepted/written counters
are captured atomically across the streams of one writer; the lower/upper
pending byte endpoints distinguish staged/queued bytes from an active block
that may be partly written. `ProductiveBoundaryProgress` retains observations
by writer directory. `price_dense_writer_queue_observation` converts one
writer observation into nominal page-cache CPU and storage work using the
existing independent profile coefficients. This is useful partial evidence
for a finite JIT continuation, but the callback, producer, GPU, other writer
directories and durable tail are not one synchronized checkpoint. The
production layout switch remains disabled pending that join and loaded
end-to-end validation. See `docs/jit_dense_writer_queue_20260923.md`.

### 2026-09-23 live writer checkpoint increment

Dense tile and variant-shard writers register weakly with the productive
controller. A background planning step can take two passes over current
writer counters and bound all registered writers at a common anchor using
monotone accepted/written bytes. It rejects an issue/output revision or
writer-identity change during capture and prices the accepted-write upper
workload with the existing independent device profiles. The result is stored
in the planner attempt audit, but does not affect the existing switch rule.
It is still not a full checkpoint: already-issued producer/GPU work,
unregistered future writers, dirty writeback, final fsync and publication
remain outside. See `docs/jit_live_dense_writer_bracket_20260923.md`.

The optional post-output source-staging worker now captures and prices a live
dense-writer bracket after each bounded metadata step. Observation storage is
charged to its retained-ledger cap, and errors are recorded as optional
evidence failures. This gives large jobs first-chunk writer observations
without putting the metadata calculation on the writer callback.
A completed public stage can now be passed directly to
`productive_source_floor`. Its prior-cost report includes the worker's measured
header setup, live-output observation and finalization as well as incremental
metadata work; passing the inner schedule alone remains explicitly narrower.
The same complete ledger is reused for a later chunk-size source rebase. A
public-stage integration regression compares it with a direct schedule without
allowing a whole-header rescan. The 136-test related A100 suite passed in
`20260923-062129-1376227`. This repairs planning-cost accounting for the
future JIT comparison, but does not enable a source-only switch.

The [JAGWAS shared-result queue increment](jit_jagwas_result_queue_20260923.md)
now observes the owned bounded multi-GPU result queue without modifying it,
checks an unchanged issue/output revision, and binds queued result ranges to
issued-but-not-part-written chunks on their original devices. Queued chunks
have completed their producer source/H2D/GPU work. A same-revision calculator
refinement removes those completed producer workloads from the candidate
checkpoint ledger while preserving the issued output backlog. The calculator
rejects combining this refinement with the old full issued-GPU shape estimate.
The staged evidence path also retains exact source reservations across initial
chunks (up to its existing bounded range/window limit), even with more than
eight phenotype partitions; it no longer ends the issue ledger after first
output. The legacy short-window switch remains suppressed in staged mode.
The affected A100 suites passed 182 and 109 tests before the join; the later
join/refinement suite is recorded in the linked note. These are accounting
and output-correctness checks, not a loaded runtime validation.

The [active indexed-writer increment](jit_indexed_writer_state_20260923.md)
now records the exact JAGWAS source range being selected or written by the
synchronous consumer. A two-sided snapshot around the queue anchor can prove
that the same result was already in the writer, then deduct its completed
producer work alongside queued results. The indexed completion event still
controls part completion, including empty chunks. An A100 164-test regression
and a two-test checkpoint follow-up passed; the held writer/queue
integration passed 13 tests. Significant-pair phenotype tiles now also expose their active indexed-writer range; the related
A100 suite passed 270 tests and 5 subtests. This further reduces conservative
JAGWAS issued-work overcount and adds significant output evidence without
asserting an elapsed completion bound.

The [staged source-cache admission update](jit_incremental_source_schedule_20260923.md)
lets a large job request up to 32,768 cached record signatures and 256 source
bounds after first output, with the extra Python-object cache charged to an
explicit host-memory reserve. A prior 8M-variant control found severe thrash
at the old 1,024-signature default; the larger cache shortened staged
metadata accounting in that ordered probe. The public staging path now uses
the requested cache geometry and records hit/miss counts. Its source/controller
regression passed 127 A100 tests. Four ordered A100 local-XFS, 8.086M-variant
metadata probes then matched exact work at both cache sizes. The 32,768-entry
cache reduced signature misses from 124,584 to 16,181; two-run median staged
CPU fell from 7.106 to 3.704 seconds, at about 45 MiB higher cumulative peak
process RSS. Server load and run order confound that comparison; it is not a
loaded whole-job speedup. See the linked note and
`results/jit_source_cache_20260923/`. This makes a post-output full-source
screen cheaper when the job admits the reserve, but no staged layout switch is
enabled.

Next, use the staged source ledger in a bounded whole-pipeline continuation,
price active indexed selection/part service and final metadata/manifest
durability, and observe the significant-pair host selector and upstream queue.
The frozen H100 absolute prediction gap still blocks a production JIT switch.

The [bounded staged partial screen](jit_staged_partial_screen_20260923.md) now joins the completed post-output source ledger, exact unissued frontier, chunk/tile/GPU candidates and mode-specific dense/significant/JAGWAS service floors with one source-stage cost. Two A100 8.1M-variant metadata-only probes completed eight bounded stage steps and two candidates: stage CPU 4.298–6.496 seconds, screen CPU 4.126–7.292 seconds, and peak RSS 571–579 MiB under shared-server load. Its service prices and output events were synthetic; the screen is evidence-only and does not establish a completion-time winner or authorize a switch.

The [device-selector launch floor](jit_device_selector_launch_floor_20260923.md) now adds source-geometry-bound mandatory CUDA launches to the significant-pair screen and joins them with the per-GPU count barrier. Live torch/CUDA/device identity is checked. Remaining selector kernels, finite output continuation and measured absolute calibration are still required before a JIT switch.

The [large-K staged screen and matched cost control](jit_large_k_staged_screen_20260923.md) now cover six versus twelve two-GPU significant-pair phenotype tiles at K=600,000 on the real 8.1M PGEN header. Spread occupancy accepts rare rational fractions up to a one-trillion denominator, and identical span/profile source pricing is calculated once while physical rereads remain charged per tile. A same-stage A/B/B/A control reduced median metadata-screen CPU from 5.236 to 1.336 seconds with identical model fingerprints. The run simulated output events and used synthetic service rates; no switch is qualified.

Final post-change A100 validation passed 312 focused-to-broad calculator/controller tests. The large-K report and script/source SHA-256 values were checked after selective pull; the report is metadata-only and remains evidence rather than a production tuning decision.

### 2026-09-23 first-chunk live staged screen increment

The [first-chunk staged screen integration](jit_live_staged_screen_20260923.md)
adds an explicit evidence-only registration to the public productive
controller. Completing the bounded source stage launches the existing
multi-layout partial screen on a separate worker after real written output.
The worker binds the actual unissued frontier and output boundary, charges its
own CPU/wall work, checks current profile/memory, and labels results stale if
the issue/output token or profile changes. Shutdown suppresses late starts and
joins an active screen. Candidate prices are supplied by the registering
caller and are not automatically bound to the active profile. The related A100
suite passed 117 tests, a broader suite passed 316, and a separate real-screen
controller integration passed. These tests use synthetic typed writer events,
not a loaded scan. No finite completion estimate or layout switch is
enabled; the next step is a cheap rebase at the later live checkpoint and
candidate-price binding, followed by output-inclusive validation.

The [live source-price binding audit](jit_staged_source_price_binding_20260923.md)
now compares every staged candidate's PGEN source coefficients with the active
per-device profile and counts which matched leaves have declared immutable
measurement targets. Mismatched source prices abort the optional screen;
matching but undeclared rates produce `source_prices_unbound`. The saved public
profile in the earlier audit binds resident-copy CPU prices only, so its source
fields cannot yet qualify tuning. The 121-test focused A100 suite passed.
The final broader A100 calculator/controller suite passed 320 tests.
Transfer, GPU, selector, writer and final-output prices, plus the finite
completion model, remain outside this source-only binding check.

The [staged work-price audit](jit_staged_work_price_binding_20260923.md) now
checks available active-profile H2D/D2H/GPU and dense-writer coefficients,
reduced-output per-device service profiles, and exact immutable host
significant/JAGWAS primitive banks where the current reduction record supplies
them. Missing shared transfer ceilings, GPU shape measurement targets and
device-significant count-transfer prices remain explicitly unbound. The live
controller reports `work_prices_unbound` rather than treating profile identity
alone as sufficient. This still needs independently measured missing terms,
finite output drain and loaded whole-job validation before a JIT switch.

The [live staged-screen retry](jit_live_staged_screen_20260923.md) now permits
zero to two extra off-writer rebases after a stale issue/output token. Each
attempt binds a new written-output checkpoint and exact unissued frontier,
regenerates candidates, and reuses the completed source metadata ledger. All
attempts share the registered CPU/wall budget and retain bounded audit
summaries; a finished job, exhausted frontier, or spent budget stops retries.
The focused A100 stage/screen/controller suite passed 120 tests, including a
retry that becomes current and one stopped by its wall budget. This improves
the chance of a usable current screen during fast scans, but no automatic
production candidate factory, finite output-inclusive completion bound,
loaded-job absolute validation, or layout switch is in place.

The staged work audit now checks the actual tensor shape-service price inputs
against immutable profile targets instead of reporting GPU shapes unbound by
construction. An exact compiled-geometry match remains a separate structural
requirement. We also corrected the shared-transfer profile schema: measured
H2D/D2H ceilings belong in optional `context.shared_transfer_capacities`,
not in the legacy CPU/DRAM/input/output `shared_capacities` map. The old
layout planner leaves this new field unused, so its existing candidate
behavior is unchanged. Earlier profiles lack the new measurements and still
report shared H2D/D2H as unbound. The focused A100 schema and binding suite
passed 57 tests; the broader calculator/controller/profile suite passed 379
tests in 71.22 seconds. Independent transfer measurements, their artifact/age
bindings, and loaded multi-GPU output-inclusive validation remain required.

The [native first-chunk candidate factory](jit_native_staged_chunk_candidates_20260923.md)
now gives `initial_chunks.staged_screen` a public opt-in path. Registration
does no cold PGEN schedule work. After written output, it uses the exact
unissued source suffix, active context prices, live dense-writer settings or
the current reduced-output record, and one to four explicitly admitted chunk
widths. The original tile and GPU ownership remain fixed, including full-K
JAGWAS variant shards. A real PGEN/controller integration verified the
background factory after native dense output; a separate test ran constructed
dense candidates through the real staged calculator. The broader A100
calculator/controller/profile suite passed 387 tests in 56.75 seconds. The
path still reports evidence only and never changes a running layout.
Independent shared-link
measurements, memory/output-ownership admission for retile or GPU moves,
finite output-inclusive makespan and loaded selection validation remain.

The [shared pinned-transfer diagnostic](jit_shared_transfer_diagnostic_20260923.md)
now measures single- and simultaneous-GPU H2D/D2H copies with raw CUDA-event
and wall spans, load snapshots, CPU affinity and verified physical GPU UUIDs.
On the shared A100 server, CPU binding changed GPU 7's observed transfer
service substantially, and repeated node-1 D2H observations also differed.
The four immutable raw reports are under
`results/shared_transfer_diagnostic_20260923/`. They are not calibrated
capacity ceilings and have not been installed into a detailed profile.
This leaves the optional live staged screen fail-closed when shared transfer
prices are missing, and the runtime switch remains disabled. The next
transfer calibration must control pinned allocation placement and competing
load, then pass held-out output-inclusive validation together with the
finite issued/queue/output continuation.
