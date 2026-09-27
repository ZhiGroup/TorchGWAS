# Diagnosing the frozen H100 runtime underestimate

The completed ranking panel selected the observed best candidate, but its
15.862-second prediction was below the 55.769-second observed median. New
diagnostics narrow the problem to transferring isolated CPU service prices into
concurrent execution. They do not justify multiplying all component rates by the
observed whole-job error.

## Preserved implementation and workload

Diagnostics run in the separate project
`/home/x/work/torchGWAS-calculator-diagnostics-h100`, mapped to
`lab-h100:/data484_4/zxie3/torchGWAS-calculator-diagnostics-h100`. They import the
unchanged frozen package and utilities from `torchGWAS-calculator-h100`.
The frozen package, plan, prices, benchmark observations and inputs are not
edited. The inspected 104 package hashes and diagnostic/cold-helper hashes
match the executed artifact.

The workload remains N=35,365, M=1,048,576, K=16,385 with 27 covariates, a real
native hardcall PGEN prefix, complete synthetic null phenotypes, chunk 1,024,
trait width 8,193, and two H100 GPUs. Each GPU has the original 4.5 GiB allocator
cap. CPU affinity, numerical pools, GPU identities and other execution settings
match the frozen profile. Input files and metadata remain on local `/data`.
Each observation runs in a fresh process with verified cold input pages and
read–scan–read controls. Detailed profiling changes only its declared profiling
flag and observational wrappers.

All four completed observations checked the original 140 variant–phenotype
pairs against the saved reference, with exact df and maximum absolute t error
1.29674e-5. Each tile emitted all 1,024 source chunks; significant output remained
empty as declared. These are sparse numerical reference checks, not a comparison
of every dense statistic.

| Observation | Output-inclusive executor | Complete API |
| --- | ---: | ---: |
| Control before | 95.556 s | 107.748 s |
| Profile 1 | 53.437 s | 63.668 s |
| Profile 2 | 28.539 s | 38.832 s |
| Control after | 52.024 s | 65.665 s |

The large variation prevents estimating instrumentation overhead or making a
speedup claim from this sequence. Sequential read controls ranged from 0.796
to 4.949 seconds for the same 6.584 GB. Loaded stage spans are diagnostic
observations and must not be substituted for independent capacities.

## Unchanged model decomposition

Reconstructing the candidate with the original source and prices reproduces
15.862022114765843 seconds exactly. No observed association duration enters this
calculation. Its declared CPU work is 38.802 CPU-seconds plus 60.147 CPU-seconds
of simulated active waits, divided across a supplied capacity of 6.865 cores.
The resulting CPU work/capacity floor is 14.413 seconds. The declared host-memory
traffic is 710.350 GB against 72.129 GB/s, a 9.848-second floor. These conditional
bounds describe the model, not measured utilization.

| Component | Frozen priced work | Profile 1 | Profile 2 |
| --- | ---: | ---: | ---: |
| Decode and input read | 18.579 CPU-s | 32.662 CPU-s | 27.631 CPU-s |
| Host selection, including cutoff lookup | 14.591 CPU-s | 45.322 CPU-s | 26.402 CPU-s |
| GPU kernel service / loaded compute spans, summed over GPUs | 27.654 s | 33.147 s | 30.279 s |

The read/decode wrapper excludes reader construction. The GPU measurements
include stream gaps and overlap other work; none of these rows can be summed
into executor wall time. Native result workers additionally used 31.870 and
34.775 CPU-seconds while handling completion events and result views. Their
spans include waiting and are not independent finish-service measurements.

## Independent controls that reject a simple rate correction

`qualify_selector.py` uses generic arrays, never GWAS input. It compares the
original 256-by-4,096 primitive shape with 1,024-by-8,192 and 1,024-by-8,193
buffers. Controls include resident arrays, a four-buffer pinned ring, recent
GPU-to-host writes, a reused output mask, a fresh output mask, and the complete
selector. Six repetitions per condition retain CPU/wall time and thread page
faults; no coefficient is fitted or published.

Large full-selector median CPU times were 6.82–7.18 ms per call, consistent with
the saved roughly 7 ms estimate. Fresh masks and recent DMA writes did not create
a persistent penalty of the size seen in loaded scans. Some first calls incurred
page faults, but steady median fault counts were near zero. The 8,193-column
tail did not reproduce the earlier large inter-worker timing difference.

`qualify_selector_contention.py` then runs fixed generic selector work in
threads, with three randomized repetitions and 40 calls per worker. One selector
took median CPU times of 6.92–7.21 ms; two concurrent selectors took 7.45–7.95 ms;
two selectors plus two independent 64 MiB copy streams took 8.88–10.63 ms.
This supports a contention effect but does not fully explain the loaded scan's
larger and variable CPU cost. Copy traffic, CPU time, wall spans and faults are
retained separately. The controls are not a new calibrated contention model.

## Completion-event experiment

The result workers' CPU usage makes completion-event policy a concrete executor
candidate. `diagnose_waits.py` runs three counterbalanced pairs with spinning or
blocking completion events, the same full workload, and profiling in both arms.
It keeps source and scientific checks fixed and records the changed flag in the
context. Job `20260922-105519-1007393` completed all six observations; results
are reported below. The default event policy is unchanged.

Artifacts are in the diagnostic project's `results/`:

- `frozen_runtime_diagnostic_v2_20260922/`: four complete real observations.
- `frozen_model_decomposition_20260922.json`: reproduced model and work ledger.
- `selector_buffer_qualification_20260922/`: 162 fixed generic observations.
- `selector_contention_20260922/`: nine concurrent fixed-work experiments.
- `diagnosis_summary_20260922.json`: combined completed evidence.
- `event_wait_control_v1_20260922/`: six completed wait-policy observations.

The initial diagnostic attempt failed before any GWAS execution because a cold
input helper used a relative path. Its log is retained. The second attempt
includes an exact copy of that frozen helper and completed all declared runs.
Calculator accuracy, loaded-resource qualification, nonempty output and broader
workloads remain open; ranking success on the earlier panel does not close them.
# Follow-up: blocking CUDA event waits

The six fresh-process controls in diagnostic job
`20260922-105519-1007393` completed. Each used the frozen input, source,
selected layout, allocator cap, cold-page checks, read–scan–read controls and
the same instrumentation. The only declared execution-setting difference was
`TORCHGWAS_BLOCKING_EVENTS`. Results are in the diagnostic project's
`results/event_wait_control_v1_20260922/`.

| Order | Event wait | Executor seconds | API seconds | Process CPU seconds | Result-worker CPU seconds |
|---|---|---:|---:|---:|---:|
| 1 | spin | 107.022 | 122.425 | 212.235 | 27.865 |
| 2 | block | 33.701 | 44.576 | 97.201 | 0.638 |
| 3 | block | 25.726 | 36.978 | 89.734 | 0.673 |
| 4 | spin | 87.564 | 100.346 | 190.019 | 29.458 |
| 5 | spin | 25.282 | 36.448 | 116.093 | 37.648 |
| 6 | block | 45.127 | 67.953 | 104.902 | 0.806 |

Blocking removed most result-worker CPU consumption and reduced total process
CPU within each adjacent pair. Elapsed times remained variable: the final
pair reversed the direction of the earlier pairs. These observations do not
establish a reliable elapsed-time speedup or an independent capacity price.
The executor spans exclude API startup/finalization; summed worker CPU is not
an additive wall-clock decomposition. All six scans checked the same 140
reference cells, had exact degrees of freedom, maximum absolute t-statistic
error 1.296739267e-5, and zero selected output rows. This does not validate
nonempty significant output or every dense statistic.

## Follow-up: unchanged-price host-serialization sensitivity

Diagnostic job `20260922-145708-1086005` completed nine historical model
calculations. It imported the frozen H100 source and kept candidate 17, the
workload, empty output occupancy and every independent component price fixed.
Only the existing declared host-serialization fraction and section ordering
changed. No new GWAS timing was collected or used to fit a coefficient.

| Declared fraction | Fluid demand | Held section first | Held section last |
| ---: | ---: | ---: | ---: |
| 0 | 15.862 s | 15.905 s | 15.894 s |
| 0.5 | 15.862 s | 15.898 s | 15.927 s |
| 1 | 17.833 s | 17.986 s | 17.995 s |

The zero/fluid calculation exactly reproduces 15.862022114765843 seconds.
None of these nine assumptions explains the earlier 55.769-second executor
median (five observations, 40.006–76.958 seconds). This is a finite sensitivity
check: fractions and orderings were not measured GIL parameters, and the
observed 15.862–17.995 range is not a guaranteed bound over all schedules or
possible resource conditions. Unmodeled work and price transfer remain open.

The fraction applies to modeled eligible host work, not to all process CPU.
Known serialized API/allocator work remains even at fraction zero. Fluid
host-serial work increases from 0.775 to 16.156 CPU-seconds across the tested
fractions. Declared total node CPU work remains 38.802 CPU-seconds. Simulated
active-wait work changes with scheduling: it falls from 60.147 CPU-seconds at
zero/fluid to 30.317 at one/fluid. It must not be held fixed and added to the
other schedules. Held-first/last convert fluid serial demand into an exclusive
token; their report's zero `host_serial.node_work` reflects that representation,
not absence of serialization. The resource-demand ledger does not tabulate
token occupancy.

This diagnostic expands the complete workload into 356,395–518,205 graph nodes.
Candidate preparation took 30.637 seconds, and each graph construction, solve
and ledger traversal took 6.692–25.383 seconds. These costs exclude imports,
artifact loading and publication. They are historical full-graph research
costs, not initial-chunk planner overhead or new calibration parameters. The
productive controller uses bounded source horizons and separate cost gates.

The report and nine immutable per-scenario records are in the diagnostic
project's `results/frozen_host_serial_scenarios_v2_20260922/`, generated by
`host_serial_scenarios.py`. All 104 frozen package hashes, the three frozen
plan/profile/price artifacts and the diagnostic script hash were verified
after pulling. The script also verifies frozen source/artifacts before and
after calculation. Frozen source and empirical records are unchanged.

The preceding attempt, job `20260922-144931-1085183`, stopped on its second
progress write because the shared save helper intentionally creates files
exclusively. Its first scenario and log are retained under the original result
name. The corrected diagnostic saves a separate file for every scenario and
uses a new result directory. No partial attempt is counted as complete.

The subsequent [per-statement selector diagnostic](selector_stage_diagnostic_20260922.md)
completed four further scans. It locates excess CPU cost in predicate/mask
scanning and records a worker-specific recurring page-fault regime, motivating
an isolated private-scratch reuse experiment instead of repricing loaded spans.

The follow-up [statistics submission control](statistics_submission_control_20260922.md)
reviews the main-thread selection dependency and tests isolated versus concurrent
GPU submission on generic tensors. Its smoke checks passed; the full comparison
is pending and supplies no replacement price yet.
