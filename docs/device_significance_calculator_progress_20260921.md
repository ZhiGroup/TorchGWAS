# Device significant-pair calculator progress, 2026-09-21

`device_significance_work.py` now traces the production selector on metadata.
Only two operations are substituted: declared dynamic nonzero output extents
and host-copy requests recorded at their original source positions. Selection
loops, comparisons, gathers and coordinate arithmetic execute from the source
function. The trace uses no real phenotype or GWAS timing. Survivor counts are
explicit inputs; the significance threshold does not determine selectivity.

The independent CUDA census covers three shapes (including a trait axis wider
than a selection block and a 1,048,576-cell limit with a tail), with empty,
sparse, dense and invalid-row scenarios. Every emitted value was compared to a
complete independent CPU reference. A100 job 20260921-175751-646613 produced
`results/device_significance_geometry_20260921`; its harness and all records
were pulled. CUDA allocator probe limit was 512 MiB.

The initial source audit caught an empty-output stride discrepancy: the real
CUDA nonzero output has strides [1,1] for shape [0,2], not [1,0]. Correcting the
metadata substitution yielded 26 passing tests and exact operation/dtype/shape/
stride/transfer reconciliation for all 12 cases in A100 job
20260921-180431-648675 (`results/device_significance_work_v2_20260921`). CPU
NumPy detach aliases are explicitly excluded from this comparison.

The census measures no durations. It confirms one four-byte pinned D2H count
per selection block and five pageable D2H copies per nonempty block, totaling
28 bytes per retained pair. Native status transfer and critical-table setup
sit outside this selector boundary and remain separate costs.

The wrapper-level nonzero ledger follows the [PyTorch 2.5.1 CUDA source](https://github.com/pytorch/pytorch/blob/v2.5.1/aten/src/ATen/native/cuda/Nonzero.cu).
It includes count reduction, blocking count transfer, flagged selection even
when empty, and conditional coordinate scatter. The count barrier precedes
selection/scatter; it is not a barrier after all selection work. Flat indices
alias the coordinate output and the counting iterator is not a materialized
input array. CUB internal passes, prefix arithmetic, workspace and launch
policy remain explicit unresolved terms. This ledger rejects unverified
PyTorch versions. Wrapper arithmetic is not claimed as issued instructions.

A separate paired implementation trial tested two typed packed copies against
the existing five-copy selector. All numerical arrays matched in eight cases
and all seven randomized paired repetitions were retained. Timing was mixed,
including a 4.4x difference on the identical empty-output path at one shape;
that undermines interpreting the apparent gains on other cases. Dense packed
output also increased the observed allocation peak. The production selector
is therefore unchanged. This control is not a calibration bank or a GWAS
performance claim. Artifacts: `results/device_significance_transport_control_20260921`.

Source work is a prerequisite for, not evidence of, accurate device-candidate
ranking. Remaining integration includes independently priced synchronization/
pageable transport, CUB service and workspace, native synchronous-result
scheduling, memory admission and full output-inclusive candidate validation.

The extended nonzero wrapper ledger passed 28 tests and all 12 exact CUDA
operation/transfer replays in A100 job 20260921-181338-651090. Artifacts are
`results/device_significance_work_v3_20260921/{tests.txt,tests.xml,audit.json}`.
The current trace additionally exposes total selector D2H bytes, including
count transport, and GPU logical accesses with the explicit nonzero wrapper
passes. CUB internals and service pricing remain unresolved as stated above.

## Source scheduling and independent phase measurements, 23:59 UTC

The device selector now records source Python API groups as well as ATen
operations. `device_significance_service.device_selection_host_primitives`
resolves each group to a typed, fixed primitive. Unknown groups are refused.
The fixed bank uses 32-element/32x32 tensors; a CPU dispatch audit caught and
corrected an alias probe that initially performed two slices instead.

A100 job `20260921-182846-654139` passed 66 tests, including native multi-GPU
selection correctness. Its audit command initially used the wrong argument;
the corrected audit in job `20260921-183143-654814` matched all 12 stored CUDA
operation/transfer cases. Artifacts: `results/device_significance_service_20260921`.
The earlier host-call trace audit is in `results/device_significance_work_v4_20260921`.

`device_selection_graph` now preserves the source's ordering across CPU calls,
count reduction, four-byte count transfer, flagged selection, coordinate
scatter, and the five owned payload copies. Host prices for blocking calls
must separate CPU work before and after the barrier. Whole blocking-API CPU
measurements cannot silently enter as dispatch plus a second GPU wait.
Pinned-count and pageable-payload transport use separate explicit prices.
Active waiting uses a declared CPU-demand scenario and is included in resource
accounting, separately from dispatch work.

Each selection block exposes yield and resume nodes. Empty output may yield
before flagged selection finishes. Nonempty output waits for all selected
payloads. The shared indexed scheduler connects resume to either synchronous
writing or bounded-queue admission. Pending selection work remains on the
same CUDA stream as the next statistics chunk, while its input H2D may overlap.
Producer cleanup drains the final selection work before another tile begins.
This extends the existing graph and shared writer; it introduces no fitted
GWAS rates or automatic tuning authorization.

Hand-calculated tests exercise a 3-second empty yield with a 13-second GPU tail,
a nonempty yield after five copies, writer overlap, stream ordering across
chunks, dual-device queue admission, and CPU work consumed by blocking waits.
The native and Python graph solvers agree. Job `20260921-183926-656513` passed
70 tests. The broader regression run `20260921-185030-658844` passed 148 tests
covering streamed graphs, multi-GPU scheduling, trait tiling, host significant
output, JAGWAS, resource waits and resource bounds. These suites overlap;
the counts should not be added as unique tests.
Artifacts: `results/device_significance_graph_20260921` and
`results/device_significance_graph_v2_20260921`.

### Primitive observations and their limits

`results/device_significance_primitives_20260921` contains 1,008 observations
for 112 primitive/context combinations: fixed tiny operations and fixed 1M
copy/nonzero operations, on cuda:0, cuda:2, and both concurrently. All nine
repetitions are retained. Caller-thread CPU, API-return wall time and external
GPU drain are separate. Some tiny blocking operations varied from tens of
microseconds to milliseconds; these are raw measurements, not qualified fixed
prices. Process CPU is diagnostic only when threads run concurrently.

A host-only C++ preload probe forwards CUDA copy/synchronization APIs unchanged
and records caller-thread CPU and wall intervals. It launches no kernels and
is enabled only in its dedicated child process. The first build failed because
the pip CUDA runtime headers lacked `crt/host_defines.h`; the successful build
used the installed CUDA 12.9 headers. Header, compiler and library hashes are
stored with the build record; the measured PyTorch runtime is 2.5.1+cu124.

The first probe retained 864 traced invocations. The corrected comparison uses
identical boundary markers in plain and instrumented paths: 1,728 invocations,
of which 864 include internal runtime clocks, over nine randomized repetitions.
Job `20260921-184809-658352` completed this collection. The phase summary
verifies exact CPU/wall interval conservation and preserves negative paired
differences. Every copy/nonzero invocation had the expected D2H byte count and
one stream synchronization. Artifacts:

- `results/device_selection_runtime_probe_20260921`
- `results/device_selection_runtime_probe_v2_20260921`, including `phase_summary.json`
- `results/device_selection_runtime_probe_build_20260921`
- `results/device_selection_runtime_probe_build_v2_20260921`

For the bulk pageable copies on cuda:0, the median recorded CPU intervals were
about 0.993 ms inside `cudaMemcpyAsync` for 4 MiB and 2.232 ms for 8 MiB. The
subsequent stream waits were only about 4.4 and 8.7 microseconds. Treating these
whole API times as fixed Python dispatch would therefore substantially
misattribute transport work. The tiny-copy instrumentation controls changed
paired wall time by approximately 14% and 25%; no automatic price publication
is justified. Markers include their own bridge overhead. The interval after
synchronization also includes immediate result destruction in this probe;
it is not yet an isolated returned-array wrapping price. Production output
ownership and later release require that boundary to be separated before
using these values directly in the source graph.

### Pinned selected-payload proposal

A frozen implementation control replaces five blocking pageable copies with
one 28*S-byte pinned allocation, five nonblocking copies into disjoint typed
views, and one final stream synchronization. Each block remains bounded by
1M pairs (28 MiB payload). It creates no dense host beta/t ring and adds no
custom CUDA code. The proposal is only in the benchmark; `reduce.py` remains
unchanged.

Job `20260921-185302-659400` completed eight cases at B=1024 and K=1024/4096,
with empty, one-per-block, sparse and dense occupancy. All arrays matched the
baseline, retained buffers survived later calls with changed inputs, and the
observed extra GPU allocation peaks matched the baseline. All nine randomized
paired repetitions are retained in `results/device_significance_pinned_control_20260921`.
Dense ratios were about 2.06x and 2.29x, but sparse/tiny cases regressed and the
identical empty branch differed by 2.43x at K=1024. Those controls prevent a
credible general speed claim or default replacement. Pinned memory retention
also needs explicit queue/producer lifetime accounting; 28 MiB per block is
not a total host-memory bound.

### Remaining work

Device selection now has a checked source graph, but still needs independently
priced CUB phases/workspace, source-consistent host allocation/release and
transport services, and integration into native candidate memory admission
and ranking. Uncertain component prices must remain explicit scenarios.
Nonempty output, format coverage and JAGWAS candidate validation remain part
of the full autotune objective. Accurate fixed cold startup is not the gate.

The isolated H100 host-selector ranking job `20260921-180931-650226` remained
active at 23:59 UTC, in its second randomized repetition. Its source, prices,
candidates and predeclared gates remain frozen. No ranking conclusion is
available from the partial panel. Automatic reduction tuning remains gated.

## Native device candidate integration (2026-09-22 00:43 UTC)

This update connects the device selector to the shared calculator. It does not
qualify automatic reduction tuning or publish a new performance claim.

### Memory admission

`nonzero_memory.py` imports duration-free allocator requests for exact installed
PyTorch/CUDA/GPU/library context and exact predicate extent. It requires empty,
one-selected and all-selected controls (two distinct controls for C=1), verifies
the allocation/release order, and refuses interpolation or timestamp fields.
The collector now records the live output allocator block separately from its
requested coordinate extent. At C=1,047,808, a 16,764,928-byte request reused a
16,777,216-byte block. That 12,288-byte excess is allocator reuse, not CUB scratch.

The v3 census covers C=1,32,4093,1047808,1048576 in fourteen occupancy cases.
Private rounded nonzero allocation is 2,048 bytes for the small extents and
23,040 bytes for the two large extents. `device_selection_memory` uses exact
full/tail extents and conservative old/current tensor lifetimes. Its 138,379,264
byte bound covered observed extra selector allocations of 7,336,960 bytes
(empty/sparse), 17,595,904 bytes (invalid) and 34,591,744 bytes (dense) at
B=257,K=4093. This is conservative allocated-memory admission, not allocator
reserved-memory prediction. Explicit reserves remain necessary.

`significant_device_memory` combines that ledger with native design/statistics,
input-only pinned buffers, status/threshold arrays, producer outputs, one global
bounded queue and the single writer's arrays/archive buffers. It assumes every
pair is retained for memory admission, independent of runtime occupancy.
No dense pinned beta/t output ring is charged for this executor.

Jobs/artifacts:

- `20260921-190611-662354`: historical allocator census containing timestamps;
  retained as raw evidence and refused by the strict importer.
- `20260921-191209-663668`: v2 census and 44 tests.
- `20260921-191727-664681`: 120 tests passed; the first audit correctly failed
  on cached output-block excess before that field was added.
- `20260921-192059-665439`: v3 fourteen-case census, 15 focused tests and all four
  selector allocation checks passed. Artifacts are pulled under
  `results/device_nonzero_memory_v3_20260921` and
  `results/device_selection_memory_v3_tests_20260921`.

### Native status and complete conditional candidate graph

`torch_scan_work` accepts the native device-significant result contract. It
charges only B status bytes before selection and refuses dense finish/copy
services or borrowed results in this mode. `native_status_service` separates
copy CPU before/after the blocking barrier, explicit pageable transfer service,
the empty malformed-status check, and both missing/invariant counts. CPU
waiting is an explicit resource scenario. The input lease can be released while
the status call or downstream selection/writer is still active.

`significant_device_runtime` now composes the same setup/decoder/statistics
calculator with status QC, the device selector's count/copy barriers, API trait
index rebasing, a shared bounded queue and one durable indexed writer. Cached
critical values are rounded/uploaded between residualization and native design
creation. Per-device transfers, shared PCIe links, storage, CPU/DRAM and optional
host serialization scenarios use the same resource graph. GPU selection tails
remain ordered with the next statistics chunk even when an empty block yields
before its flagged-select kernel finishes. Small fixed-copy return, count and
selected payload bytes are reported separately.

Selection occupancy and every component service are explicit inputs; the
bridge does not derive p-value selectivity or accept a GWAS runtime fit. Repeated
identical source/service shapes are reused only within one call. Expansion
limits reject oversized requests; scalable full voxel-job optimization still
requires a bounded sustained-throughput objective instead of expanding millions
of chunks. A status-only scan cannot be reported as complete device selection.

Job `20260921-193131-667682` passed 115 tests and exposed conflicting CPU
capacities in one new synthetic composition fixture. Correcting its supplied
capacity left the production guard intact. Job `20260921-193312-668165` then
passed all 116 tests in 18.55 seconds. The candidate bridge plus JAGWAS and tensor
memory regression batch `20260921-193854-669598` passed 93 tests in 59.61 seconds.
A preceding launch named a nonexistent test file and ran no tests; it is not
counted as validation. Final focused validation after adding the pre-expansion
selection limit is recorded separately below.

### API return and ownership measurement

The host-only runtime probe now records an API-return checkpoint before result
destruction. Both plain and instrumented paths use identical markers. The
post-return interval retains destruction and marker-bridge overhead rather than
charging it to returned-array wrapping. New independent primitives cover uint8
status copies and the two NumPy QC operations on fixed 32-row and 1M-row vectors.

Job `20260921-193312-668165` collected 2,592 invocations (1,296 with runtime event
clocks), in 324 batches over nine randomized repetitions. All expected copy
bytes and count/copy barriers passed, and the summary conserves every measured
CPU/wall interval. Raw artifacts and source/library hashes are pulled under
`results/device_selection_runtime_probe_v3_20260921`. Plain API-return CPU
medians were 40.485 microseconds for 32-byte uint8 status copies and 355.7155
microseconds for 1MiB copies; output-release intervals were 3.930 and 5.260
microseconds. These are probe observations with bridge overhead, not qualified
fixed prices or an end-to-end prediction. Nonzero timing and paired controls
remain variable; no measured GWAS duration is used to adjust them.

### Full-scope audit

Done in source: bounded host/device significant output, variant-sharded JAGWAS,
source graphs for shared resources and output ownership, device memory admission,
and the conditional native-to-indexed device candidate bridge.

Still required: resource-derived GPU/CUB service pricing and held-out transfer
checks; qualification of the independent host/transport services; device
candidate ranking with nonempty outputs; additional format and JAGWAS candidate
validation; scalable bounded optimization for very large jobs; and public
autotune wiring only for validated scopes. Fixed cold startup accuracy remains
outside the readiness gate. Device and pinned transport proposals remain
separate from production until supported by matched evidence.

The frozen H100 ranking job `20260921-180931-650226` was progressing through
repeat 3 (the fourth of five repetitions) around 00:41 UTC. Source, candidate
predictions and the predeclared five-percent regret gate are unchanged. No
conclusion is drawn from that partial panel. The broader goal remains active.

Final focused validation: job `20260921-194248-670812` passed all 42 checks in
23.31 seconds after the selection expansion guard. The 93-test bridge/JAGWAS
batch, 116-test native/status regression batch, final tests and v3 probe build
identity have all been pulled. The failed SSH push before this last launch
performed no execution; VPN status and a fresh read-only connection passed,
then the push and single background launch succeeded.

## 2026-09-22 01:19 UTC: source GPU pricing and large-job resource floors

### Exact launch partition, conditional timing prices

`device_selection_gpu_work.py` associates every installed selector launch with
an eager tensor operation or the count, flagged-compaction, and coordinate
phases of CUDA nonzero. It checks the PyTorch 2.5.1/CUDA 12.4 library and device
identity, exact source hash, mask extent, survivor count, and CUB launch policy.
The source model includes the large count-reduction partial buffer, count init,
flagged-selection work even for empty results, padded scan tiles, descriptor
writes and declared predecessor-window reads. Nonempty coordinates require a
separate integer division/remainder capacity. Logical HBM/L2 scenarios remain
attainment scenarios, not guaranteed bounds; CUB retry spinning, shared-memory
scatter and issued-instruction expansion remain unresolved.

The source references are PyTorch 2.5.1 `Nonzero.cu` and NVIDIA CCCL/CUB 2.3.1
`dispatch_reduce.cuh`, `tuning_select_if.cuh`, `agent_select_if.cuh`, and
`single_pass_scan_operators.cuh`. The CUDA 12.4 release notes identify CUB 2.3.1.
The SM80 policy was checked against fresh installed launches. The SM90 policy
is source-derived and has not yet been checked against an H100 selector census.

Job `20260921-195206-676355` captured 12 untimed cases: three shapes, each with
empty, sparse, dense and invalid inputs. All arrays match an independent CPU
reference. Job `20260921-195655-678831` passed 25 GPU-model tests, audited every
launch in those cases, and collected generic integer division/remainder data.
The bulk two-kernel int64 quotient/remainder primitive attained about 39.38
billion element pairs/s. This is not an issued-instruction rate or a qualified
transfer to the compiled coordinate kernel.

Fresh fixed generic resource measurements (`selection_resource_capacity_20260922`)
gave 1.7505 TB/s streaming add, 3.0768 TB/s cache-resident add, 15.637 TFLOP/s
FP32 square GEMM, and 1.108 microseconds per captured fill launch on A100 cuda:0.
These are observations under recorded load. No selector or GWAS timing set them.

### Held-out observations and measurement limitation

`selection_gpu_transfer_v3_20260922` freezes all predictions before observing
the selector, then records five randomized plain/profile pairs for each of the
12 cases, including exact array and installed-launch checks. All 120 calls are
retained. Two earlier harness failures are preserved: a missing plain-path
invocation, and an attempt to overwrite an immutable incremental record. The
v3 harness uses distinct incremental records and completed all cases.

For B=257, K=4093, the HBM scenario predicts GPU service of 74.7/98.9/217.1/142.0
microseconds for empty/sparse/dense/invalid inputs. Profiled sums are
161.8/233.4/331.5/230.4 microseconds. The coordinate phase is overestimated for
dense output while other phases are underestimated; a single multiplier would
conceal that difference.

The paired generic controls (`gpu_instrumentation_control_v2_20260922`) show
that profiling materially changes short GPU execution. Median event spans for
32 captured 1M-element adds were 169 microseconds plain and 827 microseconds
profiled; for 256 captured fills they were 789 versus 2,349 microseconds. The
profiler did not expose all graph kernel events, so missing kernel sums are
recorded as unavailable, not zero. Plain/profile runs share a process, and
persistence of instrumentation effects remains possible. Therefore these
controls neither recover uninstrumented selector kernel time nor qualify a
correction factor. The original prediction is retained unchanged.

### Bidirectional source proposal

`selector_simplify_ab_20260922` tests reuse of one absolute-value tensor and
finite-upper comparisons instead of repeated `abs`/`isfinite` operations. The
proposal deletes the temporary magnitude before yielding. Both implementations
match independent arrays, including NaN, positive/negative infinity, invalid
status, invalid df, and trait-window tails. Nine randomized uninstrumented
repetitions per case retain the final CUDA drain and output destruction.

Runtime evidence is mixed. At B=257, K=4093, the proposal's original/proposal
median ratios are 1.045 empty, 0.937 sparse, 0.978 dense, and 1.014 invalid.
It is retained as a benchmark proposal and **not adopted** into `reduce.py`.
There is no claimed full-GWAS speedup and no source change to the frozen H100
ranking package.

### Constant-size work accounting for very large jobs

`resource_balance.repeated_resource_floor` accepts distinct mandatory source
graph templates and integer multiplicities. A billion full chunks plus one tail
need two templates, not a billion graph copies. For candidate x and supplied
resource scenario s it computes

    L(x,s) = max(max_r sum_j n_j W_jr / C_rs,
                 max_j mandatory_single_copy_dependency_path_j).

This is a necessary lower bound inside the declared service model. Optional
handoff delays and contention-dependent active waits contribute zero to the
floor; their cost from one short solved schedule must not be multiplied as if
it were mandatory work. The full scheduler still accounts for those costs.
Copies are not assumed to drain serially, and FIFO, queue, token and cross-copy
constraints can only raise this bound.

This supports safe elimination of a candidate whose floor already exceeds a
feasible schedule under the same scenario and workload. It does not select the
winner from floors alone. Joint tuning remains minimization of the worst
supplied schedule over bounded chunk, tile, device and worker choices, subject
to memory and shared capacity constraints. Equal isolated stage times are not
a constraint: shared CPU, DRAM, storage and PCIe work must be added before
identifying the limiting resource. Binding repeated full/tail templates to the
large voxel candidate space and constructing scalable upper schedules remain
outstanding.

The broader goal remains active. The complete device path, reduced-output
candidate validation, large-job schedule integration and public reduction
autotune binding are not declared ready by these component checks.

Validation of the new repeated-work floor and the device graph/pricing regressions:
job 20260921-201833-716607 passed all 87 tests in 45.86 seconds. The billion-copy
test retains two graph nodes. No source change was made to the selector itself.

At 01:22 UTC, frozen H100 job 20260921-180931-650226 was still live:
controller PID 2461085, child PID 2518204 running candidate 4, repeat 4. Ten
candidates in the fifth repetition had been reported, and report.json did
not yet exist. No ranking gate has been changed and no partial pass is claimed.
A100 component artifacts and the 87-test log are pulled; its launched jobs
have completed. This goal turn made progress and the full goal remains active.

## Completed frozen H100 ranking

Job 20260921-180931-650226 has completed. The pulled
`torchGWAS-calculator-h100/results/significant_bounded_ranking_20260921/report.json`
records all 60 observations (12 configurations, five randomized repetitions)
and passes every preregistered gate. Selected candidate 17 is also the observed
best: chunk 1024, phenotype tile 8193, cuda:1+cuda:2, median executor time
55.769260546 seconds and selection regret 0.0. Its unchanged predicted time
is 15.862022115 seconds, so this establishes ranking on the declared panel,
not accurate absolute runtime. The plan SHA256 is
`23daffd72c1a6fce051b28fdd1da8f1b654cdec312a6b967b6799aa2aad73723`.

The scope is real native hardcall PGEN with complete synthetic null phenotypes
and empty significant output. Nonempty output, unseen shapes, GPU selection,
JAGWAS and other formats still require separate evidence. All ranking and
bounded host profile artifacts have been pulled; no frozen source was changed.

JAGWAS scope is now explicitly limited to a full phenotype panel and factor on
each GPU, with variant chunking/sharding only. Significant-pairs mode retains
phenotype tiling. A proposed online tuner would use calculator-admitted
configurations, retain early chunk results and update component measurements
before changing future work. That adaptive runtime is not implemented yet.

## Initial-chunk measurement update

The earlier online-tuner proposal has been clarified: distribute measurements
across the first chunks of one real GWAS job, retain their output, save reusable
parameters across jobs, and refresh measurements whose context or age is no
longer valid. `docs/initial_chunk_calibration_20260922.md` records the implemented
chunk mechanism, bounded sampling, immutable parameter cache and automatic
per-GPU validation/refresh. The four-process A100 refresh audit preserved
identical complete joint output and detected a controlled reader delay and its
removal. Both normal runs also differed enough to expand the checks, so this
audit does not demonstrate fewer fresh samples after a cache hit. Measurement
age is anchored to observation rather than publication; structural records
remain reusable while their dependencies match. Public reduced-output
autotuning and bounded candidate-selection integration remain outstanding.
