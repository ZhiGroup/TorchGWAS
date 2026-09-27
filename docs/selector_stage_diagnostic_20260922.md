# Selector costs during the frozen H100 scan

The completed diagnostic identifies two distinct problems: persistent excess
CPU work in predicate/mask scanning, and a recurring minor-page-fault
regime in one worker of one run. The counters alone do not establish its cause.
A single multiplier or unconditional fresh-page
surcharge would conflate them. No loaded duration is installed as a capacity
price. Frozen source, prices and historical observations remain unchanged.

## Method and verification

H100 job `20260922-150755-1093015` completed four fresh-process scans using the
frozen candidate and inputs described in
[the runtime diagnosis](frozen_runtime_diagnosis_20260922.md). Each uses
N=35,365, M=1,048,576, K=16,385, C=27, chunk 1,024, trait width 8,193, two H100s,
native hard-call PGEN, and empty significant-pair output. Input and metadata
remain on verified local `/data`, with cold-input and read–scan–read controls.

The common scan/decode/reference instrumentation is identical in all four
runs. The middle two add statement timers to every sixteenth selector call.
They sample 64 calls from each of two workers per run. The AST-based helper
verifies that removing its timer statements restores the exact original
function structure. Twenty-four generic shape/layout/broadcast cases
check unchanged inputs and identical outputs. Every real scan passed the same
140 reference checks with exact df and emitted zero result rows.

| Run | Executor seconds | Complete API seconds | Process CPU seconds |
| --- | ---: | ---: | ---: |
| Control before | 43.175 | 51.552 | 144.437 |
| Sampled 1 | 42.994 | 52.193 | 116.927 |
| Sampled 2 | 74.943 | 84.442 | 145.389 |
| Control after | 62.600 | 73.999 | 181.012 |

This spread does not identify profiling overhead or establish a performance
gain. The executor interval excludes API opening/QC and final metadata; process
CPU is summed CPU consumption, not additive elapsed time.

## Observed components

The original independent bank predicts approximately 6.756 ms of predicate CPU
and 0.344 ms of empty-mask extraction CPU per 1,024 x 8,193 call. The 8,192-wide
second tile has almost the same predicted work. Sampled loaded observations are:

| Run/worker | Predicate mean / median CPU ms | Flatnonzero mean / median CPU ms | Predicate faults after first sampled call |
| --- | ---: | ---: | ---: |
| Sampled 1 / 0 | 16.024 / 7.788 | 2.090 / 0.714 | 7 total across 63 calls |
| Sampled 1 / 1 | 9.315 / 7.527 | 0.870 / 0.701 | 1 total across 63 calls |
| Sampled 2 / 0 | 17.577 / 7.511 | 1.511 / 0.687 | 12 total across 63 calls |
| Sampled 2 / 1 | 17.790 / 11.736 | 0.907 / 0.669 | 257,988 total across 63 calls |

The last worker incurred 4,095–4,098 minor faults in every sampled call after
the first, whereas the other worker/runs had 0–6. All workers incurred roughly
4,096 faults on their first sampled predicate call. These counts do not prove
fresh allocation or first touch. The recurring faults also cannot explain the
elevated mean CPU costs in workers with almost no later faults. Means substantially
exceed medians in several cases. Replacing total work with a favorable steady
median would discard real cost. Sampling is deterministic, not a statistical
confidence bound or proof that unsampled calls behave identically.

Timers retain CPU and wall time, minor/major faults, and voluntary/involuntary
scheduling counts separately. The untimed residual combines timer bookkeeping,
Python dispatch, return and local cleanup. It is not a clean release-service
measurement and must not become an allocator coefficient.

The before/after generic controls use `numpy.zeros` for statistics and only read
those input values. Warm page-table mappings can still reference shared kernel
zero pages. These rows are descriptive controls, not physically equivalent to
GPU-written arrays and not a replacement for the independent price bank. The
limitation is recorded in the derived summary; original artifacts are preserved.

## Completed scratch-reuse control

`selector_workspace_control.py` tests private per-worker mask, absolute-value
and finite-value scratch reuse. The prototype retains owned selected arrays and
the existing predicate, coordinate and gather semantics. Its intended effect
is to make scratch allocation a bounded startup/growth cost, followed by reuse,
rather than assuming that repeated calls happen to receive resident pages.
It does not address output-array allocation or eliminate shared-resource load.

The fixed follow-up is three counterbalanced current/reuse pairs, fresh process
per scan, started only after the four-run diagnostic completed. Generic checks
include FP32/FP64, strided and read-only inputs, scalar/row/pair df, zero axes,
tail blocks, and preservation of previously returned results after scratch
reuse. All 96 generic comparisons passed, including ownership and preserved
earlier outputs. All six scans passed the same 140 reference checks with exact
df, maximum absolute t error 1.297e-5, and zero emitted rows.

| Run | Executor seconds | Complete API seconds | Process CPU seconds | Selector CPU seconds |
| --- | ---: | ---: | ---: | ---: |
| Current 1 | 43.832 | 52.108 | 131.200 | 19.703 |
| Reuse 1 | 49.447 | 59.930 | 144.576 | 19.647 |
| Reuse 2 | 47.606 | 57.392 | 167.310 | 25.657 |
| Current 2 | 42.955 | 53.702 | 162.569 | 33.660 |
| Current 3 | 54.348 | 63.842 | 168.317 | 40.030 |
| Reuse 3 | 63.771 | 74.006 | 200.337 | 44.870 |

Paired reuse/current executor ratios are 1.128, 1.108 and 1.173; API ratios are
1.150, 1.069 and 1.159. The prototype did not improve throughput in any pair,
and its selector CPU effect was inconsistent. It remains outside production.
Its buffer lookup and shape bookkeeping add work, so this does not establish
that every scratch-reuse implementation must be slower. No regression factor
or loaded service price is installed in the calculator.

Reuse 1 worker 1 still incurred a burst of 3,331 minor faults in a later sampled
call while retaining its scratch arrays. This corrects the earlier attribution
of recurring faults to fresh allocation. Automatic NUMA balancing was enabled
on this host; Linux documents that it periodically unmaps pages and traps
accesses to assess placement. That is a possible mechanism, not a demonstrated
attribution of this burst. Reclaim and page placement also remain unresolved.
See the [kernel NUMA balancing documentation](https://docs.kernel.org/admin-guide/sysctl/kernel.html#numa-balancing).

The finding also exposed a calibration compatibility gap: matching CPU affinity
did not bind memory policy. The separate [NUMA context change](numa_context_binding_20260922.md)
now records that policy without changing it or claiming a performance benefit.
The subsequent [18 independent NUMA controls](selector_numa_capacity_20260922.md)
used physically written resident arrays. They changed fault behavior but did not
reproduce the loaded CPU penalty or establish a consistent binding benefit.

## Artifacts

The separate diagnostic project contains:

- `selector_stage_diagnostic.py` and `summarize_selector_stages.py`.
- `results/selector_stage_diagnostic_20260922/`: four completed scan reports,
  per-statement observations, preserved source, schedule, and derived summary.
- `selector_workspace_control.py`: the follow-up experiment, not a production
  selector replacement.
- `results/selector_workspace_control_20260922/`: the six completed scans,
  generic checks and preserved source from job `20260922-152112-1093605`.

After pulling, all 104 frozen package hashes and both executed diagnostic script
hashes matched. The summarizer also verified the original plan and price-bank
identities. The scratch experiment's executed harness, shared base harness and
104 frozen package hashes also matched the local files. The production selector
and calibration records are unchanged by these diagnostics; the independently
tested NUMA compatibility check is documented separately.
