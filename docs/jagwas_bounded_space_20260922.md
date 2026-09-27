# Bounded JAGWAS candidate construction and cold GPU initialization

The calculator can now construct JAGWAS candidates from finite chunk axes and
independently priced device/reader contexts. `jagwas_candidate_space.py` connects
the existing census, memory, preparation and execution-graph components. This
is calculator integration. Public JAGWAS binding is now described separately in
`jagwas_public_autotune_20260922.md`; initial-chunk candidate selection remains
unfinished.

## Planning contract

`prepare_jagwas_candidates` retains the complete phenotype panel on every GPU,
partitions only variants, and shares exact PGEN censuses across compatible
chunk grids. Each shard retains its exact source interval and final chunk.
Unsupported phenotype partitions, insufficient residual rank, borrowed/dense
profiles, dense block coalescing and beta output are refused. JAGWAS projection
geometry is tracked separately from statistics geometry; only the needed
full/tail shapes are copied into each candidate profile.

`bounded_jagwas_plan` uses the existing minimax graph objective over explicit
output-occupancy and host-sharing scenarios. Memory admission assumes every
variant survives. Reader/memory admission and evaluation/chunk expansion
budgets precede preparation graph construction. Missing geometry on a feasible
candidate remains an error; an impossible candidate does not require a capture.

Preparation is constructed once per feasible candidate and host scenario, then
reused across occupancy scenarios. This makes factor/design/pinned allocation
services follow the scenario's host serialization fraction. Shared CPU and
finalization steps are supplied explicitly for each context and host scenario.
The lower-level planner still accepts fixed preparation graphs and identifies
that policy explicitly as `fixed_across_host_scenarios` in its report.

The calculator CLI accepts `model: detailed_jagwas_space` in its JSON request:

```bash
python -m torchgwas.autotune request.json --output plan.json
```

The request has exactly `model`, `workload`, `contexts`, `bounds`, `joint`,
`output`, `prices`, and `preparation_services`. Bounds contain `chunks` and
optional `max_candidates`, `max_candidate_tiles`, and `max_census_chunks`.
There is no phenotype-tile axis. `preparation_services[context][host_scenario]`
contains `library_arithmetic` (scalar/tensor per GPU), `shared_cpu_steps`, and
`finalize`. These are explicit component-model inputs, not complete-GWAS times.
The selected row contains API settings and its required environment; execution
still requires source/context freshness and live resource admission.

## Cold multi-GPU failure found by executing the plans

A100 job `20260921-214114-748769` constructed its plan, then failed in the first
concurrent factor preparation with `lazy wrapper should be called at most once`.
The previous broad `RuntimeError` catch translated that backend failure into a
misleading claim that the phenotype correlation was not positive definite.
The failed request/plan and job log are retained in
`results/jagwas_bounded_execution_v1_20260922` and the corresponding proj log.

PyTorch 2.5.1's [CUDA linear-algebra dispatch source](https://github.com/pytorch/pytorch/blob/v2.5.1/aten/src/ATen/native/cuda/LinearAlgebraStubs.cpp#L38-L50)
loads the library from the first linalg call and checks a shared invocation
counter. The observed concurrent failure is consistent with this first-use
race; [upstream issue 90613](https://github.com/pytorch/pytorch/issues/90613)
reports the same error from threaded CUDA linear algebra.

`JagwasReduction.prepare` now initializes CUDA linalg with a 2x2 FP64 Cholesky
under a process-wide lock before any full phenotype factor. CPU and meta
execution skip this. Successful initialization is reused within the process;
a failed initialization retains its exception and leaves the state retryable.
The actual independent factors run after the initialization lock is released.
Only `torch.linalg.LinAlgError` from the phenotype factorization receives the
specific non-positive-definite explanation. Other runtime and allocation errors
propagate unchanged.

The numerical scan/projection formulas are unchanged. The calculator discloses
the one-time initialization as an unpriced startup term. The source hash of
`reduce.py` changed; older whole-file factor/geometry/calibration bindings remain
historical and are not silently promoted to the current source.

## Evidence

A100 job `20260921-213735-747766` passed 150 tests with zero failures/skips in
159.21 seconds before the cold-initialization fix. These cover the finite
JAGWAS builder, exact shard/tail geometry, shared-census reuse, explicit
per-scenario preparation, missing geometry, memory/evaluation budgets, JSON CLI,
and dense-output planner/API regressions. Artifacts:
`results/jagwas_bounded_space_v1_20260922`.

`benchmarks/direct_jagwas_bounded_execution_20260922.py` executes all six emitted
configurations (B=128/256/512, one/two GPUs) on synthetic native PGEN input with
N=2049, M=1025, K=512 and two covariates. Inputs and metadata are on server-local
/data. Prices are synthetic resource controls and retained geometry is used as
an accounting fixture, so the audit is not a throughput/ranking qualification.
The audit retains all statistics, checks complete unique coverage, compares
across configurations and against independent FP64 OLS/correlation calculations,
and enforces a 2 GiB allocator cap on each GPU.

A100 job `20260921-214637-752429` passed all 110 post-fix tests with zero
failures/skips in 209.39 seconds. This includes two independent fresh Python
processes whose first CUDA linalg calls enter concurrently from two GPUs,
failed-initialization retry/error identity, specific singular-matrix errors,
CPU/meta isolation, and the existing numerical/memory/source-work regressions.
The two warnings describe the explicit reader limit imposed by prefetch depth.
Artifacts: `results/jagwas_cold_init_v1_20260922`.

The same job then completed the six-configuration execution audit in a fresh
process. Artifact: `results/jagwas_bounded_execution_v2_20260922/report.json`.
All candidates retained every one of the 1025 variant statistics. The one- and
two-GPU outputs had identical hashes at each fixed chunk size. Across chunk
sizes, the maximum absolute difference was 0.0001220703125. Thirty-four
independent FP64 checks per candidate had maximum absolute error
0.00013829697138589836 and maximum relative error 2.523471669073157e-7.
The source mount was /dev/md0 XFS at /data. Both GPUs used a 2 GiB allocator cap.
The request, plan, independent reference, indexed output, per-candidate checks
and final source hashes are retained. The first selected configuration was
B=128 on cuda:0+cuda:1; its first factorization therefore exercised the fixed
cold multi-GPU path without a preceding single-GPU scan.

This verifies numerical execution and bounded-model wiring, not sustained
candidate ranking. No association durations were used to fit the calculator.
The original failed artifact remains available; the successful audit uses a
new directory and source snapshot. The H100 ranking checkout was not changed.
