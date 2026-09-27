# JAGWAS executor and calculator progress, 2026-09-21

This work is isolated in `torchGWAS-jagwas-dev` on A100. The main H100 source
and the component prices used by the sustained significant-output experiment
remain frozen. No JAGWAS timing coefficient or tuning-readiness claim is made.

## Executable behavior

`run_linear_gwas(..., reduce='jagwas', variant_devices=[...])` now shares one
phenotype preprocessing pass, creates an independent joint factor per active
device, scans disjoint variant ranges, and feeds one indexed writer through a
bounded shared result queue. Reader workers are a global budget. Small ranges
trim idle devices. Results and metadata use the requested variant extent;
chi-square degrees of freedom use the full retained trait count.

The low-level driver refuses a shared mutable reduction object across workers.
Worker source views have independent counters and profiles. Completed disjoint
counts are aggregated only after successful exhaustion. A worker or writer
failure closes the peers and prevents final manifest publication.

Blocked phenotype QC preserves mmap/lazy input storage. JAGWAS now explicitly
converts the retained panel to the selected compute precision before shared
preprocessing; a float64 source cannot silently change a requested FP32 factor.
Both array and `.npy` inputs were checked in FP32 and FP64. Trait tiling remains
invalid for this joint quadratic form. Missing phenotype columns and impossible
residual rank are refused before factor allocation.

## Source memory and output accounting

For N samples, K retained traits and scalar size d (4 or 8 bytes), the explicit
factorization tensors simultaneously require at least

    d*N*K + (d+24)*K*K bytes.

These are the uploaded phenotype, correlation, FP64 Cholesky factor, identity,
and inverse factor. The API checks this necessary capacity on every active CUDA
device after constant-trait filtering. Available memory includes cached but
unallocated blocks. Passing this check does not establish overall feasibility:
workspace, scan tensors, allocator fragmentation and host memory remain separate.

The shared `eager_scan_memory`/`eager_memory_plan` functions accept
`reduction='jagwas'`. They retain dense statistics work, add the source-traced
FP64 reduction temporaries and persistent factor, and retain narrow outputs
across asynchronous D2H slots. The setup bound includes factor preparation and
the factor surviving design upload. Existing dense-mode defaults are unchanged.
The conservative sum of distinct temporaries is not a measured allocator peak.

`pinned_scan_work(..., reduction='jagwas')` counts the actual narrow result ring:
17 bytes per variant (beta, joint statistic, index, status and residual df), plus
the genotype ring. Status and df are released before the indexed writer, so the
owned result tuple crossing the shared queue contains 12 bytes per variant.
The writer stores 16 bytes per retained variant (int64 index and FP64 chi-square),
plus exact NPY/ZIP headers. Empty chunks produce no part and no part fsync.
Writer array accounting includes the previous iteration's values and index array
surviving assignment, and bounded serialization buffers. Metadata objects,
allocator retention and filesystem memory are not hidden in these payload counts.

Significant and JAGWAS output now share the same NPZ header accounting. Ten exact
replays against the frozen original significant implementation preserved all
reported fields at five row counts, with and without beta.

## Validation and retained evidence

The executor and memory test batch passed **217 tests and 10 subtests**. It includes
BED/native-PGEN, CPU/CUDA, independent FP64 OLS quadratic-form checks, missing and
invalid genotypes on CUDA, chunk/range changes, QC, dtype preservation, factor
ownership, writer/worker failure cleanup, NPZ bytes and memory-layout checks.
The A100 native decoder was compiled on A100; the H100 host's `-march=native`
binary is not portable to this CPU.

Duration-free CUDA captures covered N/B/K = 257/13/7, 2049/128/512 and
4097/512/2048, each in FP32 and FP64 and in prepare/reduce phases. All 12
operation sequences matched the current meta ledger after removing zero-work
meta detach aliases and accounting separately for NumPy wrapping and H2D upload.
Profiler timing fields were discarded. Raw trace files were temporary.

| K | FP32 prepare kernels | FP32 reduce kernels | FP64 prepare kernels | FP64 reduce kernels |
|---|---:|---:|---:|---:|
| 7 | 14 | 17 | 13 | 15 |
| 512 | 26 | 17 | 25 | 15 |
| 2048 | 26 | 17 | 25 | 15 |

On the K=512 and K=2048 captures, observed extra allocated bytes during factor
preparation equaled the explicit live-tensor floor above. Small allocations
included rounding overhead. These observations validate this source term for the
captured cases; they do not replace vendor-workspace or allocator accounting.
FP64 projection used a different compiled GEMM family from the FP32 scan, so the
existing FP32 component rate cannot simply be reused.

Artifacts, relative to this development directory:

- `results/jagwas_variant_checks_v8_20260921/tests.log` and `tests.xml`
- `results/jagwas_geometry_20260921/census.json`
- `results/jagwas_source_audit_20260921/audit.json`
- `docs/jagwas_development_baseline_20260921.json`
- `benchmarks/direct_jagwas_geometry_20260921.py`
- `benchmarks/direct_jagwas_audit_20260921.py`

That audit binds 98 package source files and identifies seven changed/added files
relative to the isolated baseline. Geometry is reused only after verifying the
execution source hashes; subsequent memory-ledger changes are recorded explicitly.

## Shared scheduling and projection service

`indexed_schedule.py` now supplies both significant-pairs and JAGWAS schedules.
JAGWAS requires one explicit shared preprocessing graph and a separate factor
and design preparation graph for each device. Host conversion and finite
filtering run after dequeue under the single consumer, including empty chunks.
The graph retains bounded queue credits, arrival ordering, producer cleanup,
terminal sentinels, shared resources and the final durable-writer drain.
Twenty exact replays against the frozen significant scheduler preserved all
nodes, demands, tokens, FIFO actions and solved results.

The shared tensor service accepts the captured FP64 projection GEMM family and
requires separate scalar FP64 and Tensor Core FP64 capacities. It counts
eight-byte operands, compiled tile padding and the observed split-K work.
The K=2048 capture's in-kernel split-K synchronization remains explicitly
unpriced. Source host APIs are grouped across their ATen operations, so one
`isfinite` API is not charged as several independent Python calls. Typed fixed
host prices must be supplied; the tests use synthetic resource controls, not
calibration or association durations.

`torch_scan_work` can compose this projection after dense statistics and before
D2H for the native FP32 JAGWAS path. Host submission overlaps earlier GPU work;
the component spans are not simply added. Each chunk transfers and copies five
arrays totaling 17 bytes per variant. Its index survives to the consumer, while
status and residual df are released by the finish worker. A matching five-array
32-row finish measurement is mandatory: the dense four-array baseline cannot
be reused. Its copy baseline is 544 bytes, versus 416 for dense K=1 output.
The dense full-output runtime entry point rejects reduced output so it cannot
silently attach the dense writer to a joint scan.

The projection currently begins with an explicitly empty modeled cache because
the two source traces have independent storage namespaces. This boundary
approximation is reported. Scan-to-discard timing excludes preprocessing,
factor preparation and indexed output; their graphs and prices must be supplied
before a complete candidate can be evaluated.

The earlier projection integration batch passed **146 tests**. It checks the actual
finish function, captured projection kernels, separate arithmetic capacities,
host/GPU composition, result ownership, missing-evidence rejection and existing
dense/significant behavior. An earlier component batch passed 109 tests. The
earlier audit binds **99 package source files**, reconciles all 12 CUDA operation
captures and prior tensor-storage ledgers, and exactly replays six dense-scan
accounting cases in addition to the 20 significant schedules.

Additional retained evidence:

- `results/jagwas_scan_work_checks_v3_20260921/tests.log` and `tests.xml`
- `results/jagwas_service_checks_v2_20260921/tests.log` and `tests.xml`
- `results/jagwas_service_audit_20260921/audit.json`
- `benchmarks/direct_jagwas_service_audit_20260921.py`

## Independent host and writer observations

The fixed32 FP64 operation bank matches the actual reducer's operand dtypes and
transposed GEMM RHS. Two bank tests passed. The five-array owned finish probe
uses the source finish function, with the 544-byte baseline and `return_df=False`.
Its one- and four-worker observations and raw controls are retained; they are
context-specific observations, not scan timing coefficients.

The two-array JAGWAS writer now has separate selection, NPZ serialization,
page-cache copy, storage and durable part-commit stages. Selection explicitly
prices conversion, finiteness, nonzero, index offset and gather. Empty chunks
still pay the five host calls but emit no archive. A seekable NPZ rewrites its
local headers: the submitted bytes exceed the final file extent by 125 bytes
per nonempty part. CPU copying uses submitted bytes; storage uses final extent.
The archive schema is checked before accepting primitive prices. Forty-three
writer/reduced-output/schedule checks passed in the initial writer batch.

The CPU-owned writer controls use fixed empty and 1M-element operations, plus
small/large two-array NPZ serialization. They do not establish fresh-DMA cache
behavior, serial CPU fractions, final metadata cost or disk commit service.
Survivor counts remain explicit scenarios, and memory admission must cover all
variants retained.

The GIL meter's positive/negative and interval controls passed in the initial
rerun, but this is now explicitly separate from transfer qualification.
`gil_service_prices` marks total CPU and serial transfer unqualified. The new
comparison validates runtime, affinity, GPU contexts, executable bank/library
identity, sample coverage and raw CPU conservation against unhooked endpoints.
A zero after empty-call subtraction is unresolved measurement resolution, not
a zero-cost API. Its 10% comparison tolerance is a declared diagnostic choice,
not an autotuning accuracy gate.

The first bracket did not support transfer: unhooked single-GPU FP64 GEMM CPU
medians changed from 136.28 to 45.23 microseconds between its endpoints. Three
further brackets alternated two persistent processes after excluded warm-up.
The first two also failed compatibility, with substantial endpoint drift,
particularly for threaded dispatch. The third failed the meter-correction
invariant for `joint_all_bool_rows`; it is rejected without clipping the
correction or deleting a repeat. These observations do not isolate the cause
of variation and must not be used as calibrated scan costs. The raw records,
checkpoints, comparison reports and failure logs are retained.

Evidence:

- `results/jagwas_independent_primitives_20260921/finish/primitives.json`
- `results/jagwas_independent_primitives_v2_20260921/dispatch.json`
- `results/jagwas_independent_primitives_v2_20260921/writer/`
- `results/jagwas_dispatch_bracket_20260921/`
- `results/jagwas_interleaved_dispatch_20260921/`
- `results/jagwas_qualification_checks_v2_20260921/` (62 tests passed)

## Exact factor workspace requests

A size-only query calls the installed cuSOLVER `Xpotrf_bufferSize` for the
single-matrix FP64 lower/default operation selected by PyTorch 2.5.1 CUDA 12.4.
The query allocates no large trait matrix and records no timing. Device requests
for K = 1, 7, 512, 2048, 8192, 16384 were respectively 160, 160, 1152, 16512,
262272 and 1048704 bytes; host workspace was zero for these cases.

Independent small Cholesky allocator checks at K = 7, 512 and 2048 left only
1024, 512 and 512 bytes beyond the rounded output plus queried workspace.
These checks concern that operation, not total setup or reserved GPU memory.
The original query and controls are in `results/jagwas_workspace_20260921/`.

`cusolver_memory.py` accepts only exact queried shapes from the same runtime,
GPU architecture, host, preferred backend and installed library metadata. It
refuses interpolation. `eager_memory_plan` adds the rounded request to factor
preparation when this evidence is supplied and reports the host request
separately. It leaves scan buffers unchanged. Library file metadata is an
installation identity, not a binary digest. Factor info/error tensors,
triangular-solve workspace, driver allocations and allocator reservation remain
explicitly unresolved; the plan is still an incomplete memory candidate.

The completed integration batch passed **266 tests**, with zero failures,
errors or skips. Its audit binds **101 package source files**, verifies all
12 operation/storage captures, and exactly replays 20 significant schedules
and six dense-scan accounting cases. A separate untimed audit composes all
six queried factor workspaces. All 101 local source hashes were checked
against the completed remote record after pulling it.

- `results/jagwas_memory_writer_audit_20260921/tests.xml`
- `results/jagwas_memory_writer_audit_20260921/source_replay/audit.json`
- `results/jagwas_memory_writer_audit_20260921/workspace.json`
- `results/jagwas_memory_writer_audit_20260921/meter_validity_only.json`
- `results/jagwas_memory_writer_audit_20260921/rejected_dispatch_comparison.json`
  (local retained copy of the third interleaved bracket rejection)

## Bounded JAGWAS candidate model

`jagwas_candidate.py` connects the existing native scan/projection work, indexed
writer service and shared execution graph. The candidate contract uses the
executor's full-file, chunk-balanced variant ranges, complete FP32 phenotypes,
owned results, one global reader budget and one global result queue. Every GPU
retains the entire trait panel. Missing tail geometry is rejected rather than
interpolated from a full chunk. Dense and significant defaults retain their
existing borrowed-result contract.

Memory admission counts one shared phenotype preprocessing allocation set,
per-device factor/design storage, the narrow pinned rings, active native reader
workspaces and owned result workers. The global queue carries 12 bytes per
variant (beta/stat/index after status and df release); it is counted once.
The writer and consumer budgets assume all variants survive, regardless of the
selected timing scenario. Vendor host-workspace requests are charged per
factor. Reserves still cover mapped input, metadata, library and allocator
terms outside the explicit array ledger.

Preparation services must be supplied explicitly: a shared preprocessing graph,
a factor/design graph and cleanup for each device, and final publication/commit
service. Their dimensions include N, K, covariate rank and column count, chunk
size and ring depth. This prevents reuse of preparation costs for a different
allocation geometry. The adapter does not fill missing stages with zero or
reuse dense-output writer prices. Independent factor/preparation price
construction is still outstanding; current integration tests supply synthetic
component services, not observed association runtimes.

The development optimizer solves

    minimize over admitted candidates theta: max over declared scenarios s T_graph(theta, s)

subject to the total reader budget and host/per-device memory capacities,
including caller-specified reserves. Candidate count, scenario evaluations and
source-chunk graph expansion have explicit bounds. It can compare chunk sizes,
ring depths, supplied worker contexts and GPU sets. The objective is completion
time under the shared resources; it does not force stage times to be equal.
Scenarios specify empty/full retained output and CPU serialization assumptions.
They are neither statistical confidence intervals nor guaranteed runtime bounds.
The result includes executable JAGWAS variant-device arguments, while
`automatic_selection_ready` and `selection_validated` remain false.

The candidate integration batch passed 119 tests. The latest focused batch
passed 34 tests after extending preparation identity to include covariate
columns, chunk size and ring depth. These check global memory accounting,
full-retention admission under an empty timing scenario, bounded optimization,
exact variant/tail extents, typed projection integration, owned cleanup,
shared storage/PCIe demand, required preparation/finalization and API argument
mapping. Their services are synthetic accounting controls, not timing validation.
The latest audit binds 102 package files and preserves 12 geometry/storage
captures, 20 significant schedules, six dense-work replays and five default
trait-shape replays. All 102 local hashes match the remote audit.

- `results/jagwas_candidate_checks_v2_20260921/tests.xml`
- `results/jagwas_candidate_checks_v4_20260921/tests.xml`
- `results/jagwas_candidate_checks_v4_20260921/source_replay/audit.json`

## Full chunk geometry and factor arithmetic

The duration-free A100 capture in
`results/jagwas_chunk_geometry_20260921/census.json` uses the same N=2049,
K=512, C=2 workload for B=1, 128, 256 and 512. The statistics kernels number
66, 66, 65 and 65; projection has 17 kernels at every width. No elapsed
kernel/association durations or fitted service coefficients are retained.

The one-variant tail exposed two unsupported paths. Statistics uses a captured
FP32 NSP matrix-vector kernel plus an explicit split reduction. Projection uses
a scalar FP64 matrix-vector kernel. The calculator now admits these verified
families, including the exact N2049/B1/W515 NSP case, and rejects unsupported
layouts. FP64 matrix-vector arithmetic uses the independent scalar capacity;
it does not require or substitute Tensor Core throughput. Internal padding and
workspace-layout uncertainty remain explicit logical-floor terms.

`test_jagwas_actual_candidate.py` composes the actual statistics and projection
services without mocking either stage, over all three full chunk sizes and
one/two-device plans. It checks exact tail accounting, narrow D2H payloads,
empty/full-retention output and the bounded optimizer. Prices and preparation
graphs in these checks are synthetic controls, so they do not validate rankings
or elapsed times. The service audit also compares the three previously
supported FP32 full-chunk services exactly against the frozen baseline.

`jagwas_factor_arithmetic_work` adds a source-checked mathematical ledger for
preparation on each active device. Correlation has 2NK^2 useful multiply/add
FLOPs and K^2 normalization divisions. Conventional dense Cholesky has
(K^3-K)/3 multiply/add FLOPs, K(K-1)/2 divisions and K square roots. The dense
triangular solve has K right-hand sides, K^2(K-1) multiply/add FLOPs and K^2
divisions. The ledger retains the NK upload and the 8K^2 persistent factor.
Candidate reports include the number of independent factors and aggregate
upload/arithmetic counts. These are mathematical work counts, not the closed
vendor libraries' issued instructions or runtime service.

An executable scalar Cholesky/triangular-solve reference checks the counts and
numerical factor for both FP32 and FP64 inputs. It reconstructs the symmetric
matrix from the lower triangle, matching Cholesky's input convention; separately
computed FP32 Gram triangles can differ by rounding. The numerical tolerance
remains 1e-12 against the FP64 factor and inverse identities.

Final validation is retained in
`results/jagwas_actual_candidate_v5_20260921/tests.xml` (89 tests, no failures,
errors or skips) and `source_replay/audit.json`. The audit reconciled 12 earlier
operation/storage captures, 20 significant graphs, 6 dense-work cases, 5 default
trait shapes and 3 unchanged FP32 tensor-service cases. All 102 local package
source hashes match the completed remote audit. Main H100 execution remains
frozen while its sustained benchmark runs.

## Constructed preparation and independent FP64 observations

`jagwas_preparation.py` now constructs one residualization graph on the first
GPU and a separate factor/design/pinned-buffer graph for every active GPU.
It reuses the existing source-derived setup phases, including blocked uploads,
downloads, host copies and pageable allocation/release when supplied. Shared
residual storage is released only after all workers finish. Per-device setup
has six JAGWAS ring arrays per slot and explicitly assumes fresh pinned pages.
Common CPU input/count/basis work and final metadata/commit remain mandatory
caller-supplied independent services.

Factor service uses one fixed N=65, K=32, FP32-input reference with six phases:
upload, correlation/normalization, FP64 cast, Cholesky, identity and triangular
solve. Additional arithmetic comes from the source-checked work ledger, while
transfer and logical GPU traffic use independently supplied capacities. Scalar
and Tensor Core library arithmetic are explicit scenarios. Division/sqrt use
scalar-equivalent work; dependency latency and instruction-specific cost remain
unresolved. The model states its blocking-phase approximation: the reference
synchronizes after each phase, while the executor can overlap asynchronous host
submission. It therefore does not claim exact preparation runtime or a bound.

`results/jagwas_factor_primitives_v2_20260921/report.json` retains complete
five-repeat, randomized-phase-order observations on A100 devices 0 and 2.
The CPU intervals include the explicit call-plus-device-synchronize boundary.
No association duration, phenotype/variant timing sweep, or fitted coefficient
was used. The initial capture correctly needed a transfer-only case for the
upload phase, which launches no compute kernel.

Independent FP64 capacities are retained in
`results/jagwas_fp64_lt_capacity_20260921/report.json`. Fixed 4096-square vendor
cuBLASLt operations use explicit FMA versus DMMA numerical-implementation masks,
the first matching vendor heuristic, and at most 64 MiB workspace. Algorithm
flags, untimed kernel names and 32 independent CPU dot products verify the
separation; five paired repeats preserve both rates under recorded load.
Device 0 measured 3.6877e12 scalar and 1.6847e13 Tensor Core useful FLOP/s;
device 2 measured 3.5140e12 and 1.6578e13. These are observations for the selected
algorithms, not peak capacity guarantees or qualification under concurrent GWAS.
No custom CUDA kernel was added.

The earlier DGEMM and GemmEx pedantic-mode attempts are retained in jobs
20260921-163102-609720 and 20260921-163353-610772. Both still launched a captured
FP64 Tensor Core kernel, so neither was accepted as a scalar measurement.
The successful path uses the documented numerical implementation flags:
https://docs.nvidia.com/cuda/archive/12.4.1/cublas/index.html#cublasltnumericalimplflags-t

`attach_factor_calibration` requires matching device/runtime/affinity/thread
context, the executor source hash, complete unique repetitions, conserved
CPU/wall partitions and FLOP/time rates, matching instruction flags and kernel
families, and the CPU reference error check. It reconstructs the medians from
all observations and records observation hashes. Other profile prices remain
untouched. Transfer qualification remains false; successful accounting tests
do not certify throughput rankings.

Final preparation validation is retained in
`results/jagwas_preparation_v4_20260921/tests.xml`: 127 tests passed with no
failures, errors or skips. The source audit retained all prior operation/storage,
significant/dense graph and FP32 service replays. All 103 current package source
hashes match that completed audit. The new checks consume the real two-device
phase/FP64 observations and verify their integration into one/two-device
preparation graphs; remaining synthetic prices keep this an accounting check.

## Remaining work before reduction autotuning

Independent FP64 observations and a source-derived preparation builder now
exist. Concurrent-context transfer of those observations and the factor model's
blocking/roofline approximation remain unvalidated. Host dispatch/finish/writer
observations likewise have unqualified transfer and serial partitions; the
bracket failures above are retained. Common CPU count/basis services and final
metadata/commit still need production pricing. Accurate fixed cold startup is
not the sustained-ranking gate; repeated preparation and output costs must
remain in the relevant executor timing boundary.

Unknown survivor counts are explicit scenarios; memory admission always uses
full retention. Shared host arrays, worker buffers and the global queue are now
counted. Remaining library, allocator, mapped-input and metadata allocations
still require quantified reserves or additional source accounting. Exact
Cholesky workspace requests are integrated for queried runtime contexts.

The public detailed autotuner continues to reject reductions. Sustained ranking
validation, device-side significant selection, nonempty significant output and
additional input formats remain outstanding. The voxel job is a scale reference,
not evidence that a full-rank JAGWAS factor can exist for millions of traits.

## Completed sustained host-selection panel

The frozen H100 job 20260921-131040-516960 completed all 60 observations
(12 candidates, five randomized repetitions). All numerical, checker-overhead,
cold-input, empty-output, minimum-duration and CUDA-cap checks passed. The
predeclared selection target failed: selected candidate 17 (B=1024, T=8193,
two GPUs) had median 101.160603 s, while candidate 16 (B=1024, T=4096, two GPUs)
had median 89.471496 s. Selected median regret was 0.130646159 against the
predeclared 0.05 target. The model predicted 32.359788 s for its selection.
These are scan-and-write executor intervals including repeated tile setup,
not whole-process cold-start times. The result does not justify automatic
selection, and no rate or correction was fitted to these association timings.

The authoritative report is in the main project's
`results/significant_host_sustained_20260921/report.json`; all observation files
were pulled after completion. The prepared independent cache-state and
flat-coordinate controls started as H100 job 20260921-164810-615973 only after
the sustained job's process was confirmed absent and its final report existed.
Main executor source remains unchanged during these controls.

## Subsequent host significant-selector work

The development host significant selector and its work ledger now use bounded
predicates and flat coordinates. See `significant_host_selector_progress_20260921.md`
for exact scope, completed controls, 138 passing tests plus 15 subtests, fresh
independent primitive observations, and the running matched full-executor A/B.
This changes `linear.py`; old whole-file source bindings in geometry audits must
be recaptured/rebound before reusing them as current evidence. The JAGWAS factor
primitive artifacts are unchanged, and automatic reduction tuning remains gated.

## Comparable executor clocks and resource explanation

The development API now retains `scan_and_write_seconds` and adds
`scan_setup_seconds` and `setup_scan_and_write_seconds` to streamed output
metadata. The latter begins immediately after input QC and ends at the same
captured writer-return timestamp. It includes shared preprocessing, all tile
or device preparation, iterator cleanup and requested output publication;
input opening/QC and subsequent run metadata are excluded. Durability still
depends on the requested fsync setting. These overlapping diagnostics are not
extra additive entries in `phase_seconds`.

This fixes a comparison error for JAGWAS: a single active device prepares its
factor before the writer starts, while multiple device workers prepare theirs
after writer entry. The old writer interval therefore charged factor setup
only to the multiple-device case. JAGWAS calculator reports now explicitly name
`setup_scan_and_write_seconds` as the observation metric.

Remote A100 job 20260921-174501-642262 passed 80 tests plus 10 subtests, including
native single/multiple-GPU JAGWAS numerical checks. Controlled-clock API tests
exercise one, two and three active workers and prove factor inclusion, QC
exclusion, and a fixed writer-return endpoint. Artifacts:
`results/jagwas_executor_timing_20260921/{tests.txt,tests.xml}`.

Both JAGWAS and host significant-pair estimates now include `resource_balance`.
For each declared shared resource r, the calculator reports its supplied
capacity C_r, modeled work W_r, capacity time L_r = W_r / C_r, and average
capacity fraction L_r / T_graph. Work includes simulated active waits and only
conditional handoff service that actually ran. These are model quantities,
not measured hardware utilization. Dependency stalls, resource overlap, queue
limits and unpriced terms can keep the graph time above every resource bound.

The bounded optimization remains

    minimize over feasible configurations x: max over supplied scenarios s T_graph(x; s)

where feasibility includes host/device memory, shared reader budgets, allowed
GPU assignments, output semantics and the bounded search axes. Significant
output permits phenotype tiling; full-rank JAGWAS keeps all traits together and
requires K <= N - covariate_rank - 1. Equal stage durations are not an added
constraint: making a faster stage wait cannot improve completion time.

These edits are isolated to the A100 development checkout. H100 A/B and queued
all-candidate validation retain their frozen source and timing definitions.
The API file hash has changed; prior whole-file JAGWAS geometry bindings are
historical until a subsequent current-source audit. No automatic reduction
selection has been enabled by these timing/reporting changes.

Resource accounting passed 78 remote tests in A100 job
20260921-175031-644801. Tests cover CPU/DRAM sharing, dependency gaps, skipped
versus active conditional handoffs, resource-consuming waits, malformed or
partial solutions, and both reduction candidate calculators. Artifacts:
`results/reduction_resource_balance_20260921/{tests.txt,tests.xml}`. The tested
report additions do not alter execution-graph service prices or optimization
scores.

## Explicit JAGWAS scope, 2026-09-21

JAGWAS is restricted to workloads whose complete retained phenotype panel and
joint-test state fit on every active GPU. Its tuning axes include variant
chunk size, variant-device assignment and reader/queue resources. Phenotype
partitioning is excluded, both across successive tiles and across devices.
The API rejects `trait_block` and `trait_devices` before loading inputs or
selecting a device. `variant_devices` remains supported, with the full panel
and an independent full factor on each active GPU.

The calculator already enforces variant-only partitioning and per-device
full-panel/factor memory admission. More GPUs do not combine memory for one
joint factor. Bounded input QC may inspect columns in blocks without dividing
the statistical test. Rank and complete-phenotype requirements still apply.
The early factor capacity check is a necessary lower bound; scan buffers,
library workspace and allocator reserves require additional memory. Large
voxel jobs requiring phenotype tiling remain in significant-pairs mode, which
is a different requested output rather than an automatic fallback.

This API validation change does not open public reduced-output autotuning.
Whole-file API source hashes in earlier geometry artifacts remain historical;
the frozen H100 validation checkout is unchanged.

Remote validation: A100 job 20260921-203519-727948 passed all 103 tests in
95.68 seconds. This covers early partition rejection before input/device setup,
full factors on each GPU, serial and two-GPU agreement with independent FP64
OLS at two chunk sizes, and the existing calculator/memory regressions. The
two warnings report the expected decoder limit imposed by prefetch depth.
Artifacts: `results/jagwas_full_panel_scope_20260922/{tests.log,tests.xml}`.
The broader autotune goal remains active; this change fixes the JAGWAS scope.

## Bounded candidate construction and cold multi-GPU initialization

`docs/jagwas_bounded_space_20260922.md` records the finite JAGWAS candidate
builder and JSON calculator entrypoint. The shared PGEN census now feeds
variant-only candidates with full phenotype/factor state per GPU and separate
statistics/projection geometry. Source-derived preparation follows each host
sharing scenario and is built only after reader/memory and expansion budgets.

Executing emitted configurations exposed a cold CUDA linalg initialization
race. `JagwasReduction.prepare` now initializes the CUDA linalg dispatch with a
tiny locked Cholesky, then permits the full factors to run concurrently. It
also preserves backend/OOM exceptions rather than labeling every runtime error
as a non-positive-definite phenotype matrix. This changes reduce.py's source
binding and introduces a disclosed one-time unpriced startup term. Historical
component evidence is not silently rebound. Public reduced-output autotuning
and cached initial-chunk candidate selection remain unfinished.
