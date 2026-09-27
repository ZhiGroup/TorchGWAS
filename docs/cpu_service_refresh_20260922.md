# Immutable evidence with productive CPU refresh

This extends [price binding](price_binding_20260922.md) with a bounded fresh
check before reusing one independent CPU coefficient. An observation remains
immutable; its eligibility for reuse can change.

The [t-only transport correction](t_only_transport_20260922.md) supplies a
concrete output-layout invalidation case: the old four-array finish measurement
cannot be silently reused for a three-array result path.

The [public initial-chunk bridge](public_initial_chunks_20260922.md) reuses bound
immutable prices and checks their original ages during productive scans.
[Public resident-copy refresh](public_refresh_20260922.md) now schedules these
small probe batches after useful output and publishes derived immutable profiles.

| Evidence | Reuse rule |
| --- | --- |
| Structural work, geometry and device properties | Exact relevant dependencies; no empirical age renewal or expiry |
| Measured CPU/GPU, transfer and storage service | Original observation time, producer lifetime and a possibly shorter caller lifetime; compatible measurement context |
| Available memory and current contention | Fresh live observation for each admission; cached values are audit evidence only |

The new `CpuServiceRefresh` implements the fresh-check policy for fixed-work CPU
samples. It does not infer an independent capacity from loaded scan intervals,
implement probes for the other resources, or enable a complete public autotuner.

## Incremental measurement

Construction performs no cache I/O or probing. The integration invokes
`advance()` inside `ProductiveTuningRun.planning_step`, after actual written
output. Each call makes at most two samples; each sample in this audit is four
copies of an 8 MiB resident NumPy array. Raw thread CPU time is additive across
the four copies. The coefficient is the median of seven repeat means per byte.
Allocation and first touch of the two 8 MiB arrays, sample validation, record
lookup, publication and profile binding are all charged to the productive step.
The buffers are released on completion, failure or scan exit. The audit's host
reserve covers these additional buffers; this is not a general probe-memory
admission policy.

A compatible, unexpired cached record receives two fresh samples. The declared
drift heuristic rejects a median ratio outside [0.5, 2] when its CPU difference
exceeds 2 microseconds per fixed work batch. A consistent check retains the exact
old coefficient, record hash and observation time. Fresh validation samples are
reported separately and never renew its lifetime. A newer record published by
another job cannot substitute for the particular record being checked.

On drift or expiry during validation, those two samples become the start of a
seven-sample replacement, completed on later productive callbacks. Missing,
expired or incompatible records also require seven samples. A stable refresh
therefore uses batches of 2, 2, 2 and 1 after four written chunks. Publication
uses the earliest actual sample time, not the time the last callback completes.

Incomplete and failed windows are not published. A complete window with a
maximum/minimum sample ratio above four and an absolute difference above the
declared tolerance is refused as unstable. A window that has aged out before
publication is also refused. If publication itself outlasts the lifetime, the
complete immutable artifact can remain in the cache, but it cannot become a
usable decision result. Subsequent binding and final decision checks enforce
the same original lifetime.

The subsequent [temporal-window check](cpu_window_stability_20260922.md) also
compares early and recent sample medians using the declared drift tolerance.
A window may therefore be refused despite passing its overall spread limit;
matching a cached median alone cannot authorize reuse of inconsistent samples.

These thresholds detect substantial changes; they are not statistical
confidence bounds. In particular, passing this check does not establish the
five-percent synthetic model-error scenario used in the integration audit.

## Cost and scientific correctness

The existing productive budget charges every callback and includes their
cumulative wall cost in the later switch decision. No useful scan chunk is
discarded or rerun for calibration. CPU/window overruns stop optional tuning
while the admitted source scan continues. Samples are cooperatively bounded;
an individual native operation cannot be forcibly preempted by this controller.

The two-GPU native-PGEN/JAGWAS audit uses the complete 512-phenotype panel on each
variant shard, 2,049 samples, 4,097 variants and two covariates. Chunk sizes are
128 and 256. The active target remains
`owned_result_copy_scenario.resident_cpu_seconds_per_byte`; every model call
asserts that source-derived additional copy bytes consume this coefficient.
All other component prices remain synthetic controls. Forecasting begins after
a later written part and uses three bounded windows at the actual issue
frontier. Their lengths adapt to the remaining source; insufficient source
stops optional forecasting.

The budget is explicitly an audit setting: up to five steps, 10 seconds of
planning-thread CPU and a 30-second early window, with a 10-second remaining
horizon precheck and 10-millisecond reserve. It is not a production policy or
a measured remaining-runtime forecast. The calculation still uses declared
0–1 second in-flight boundary scenarios. This work establishes evidence reuse
and measurement scheduling, not prediction accuracy or beneficial autotuning.

## Validation

Remote job `20260922-083405-968430` passed 108 tests in 22.57 seconds. Job
`20260922-083616-968940` passed another 120 in 17.38 seconds. The 228 distinct
cases cover immutable caching, binding, dependency changes, changed rates in
both directions, expiry during check and publication, failed and noisy probes,
concurrent publication, cumulative cost, continued source coverage, loaded
observation separation and analytical forecast integration.

Test logs are `results/cpu_service_refresh_v1_20260922/tests.log` and
`results/cpu_service_refresh_broad_v1_20260922/tests.log`.

## Verified productive runs

| Job | Evidence action | Samples per callback | Age at decision | Measurement callbacks | Forecast callback | Cumulative cost plus reserve |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| First | Publish seven-sample record | 2, 2, 2, 1 | 1.70 s | 0.376 s | 1.436 s | 1.822 s |
| Reuse | Consistent check; retain original | 2 | 103.63 s | 0.229 s | 0.956 s | 1.195 s |
| Shorter lifetime | Expired; publish new record | 2, 2, 2, 1 | 0.61 s | 0.201 s | 0.480 s | 0.690 s |
| Follow-up | Consistent check; retain original | 2 | 221.65 s | 0.099 s | 0.764 s | 0.873 s |

The first, reuse and follow-up jobs used the same record
`3b0a29694eafa15f65e7225f06bed3d76d6e8e0a8c6bc29ff3624cb34a4b6dcd`,
observed at Unix time 1790084091.540527. Its coefficient was
4.967086613178045e-10 CPU seconds per byte. The two fresh checks produced median
ratios of 0.700 and 0.611, both within the declared factor-two band. They kept
the old coefficient; they did not relabel their measurements as a fresh copy of
the original record. These checks allow substantial variation and do not
certify precise runtime predictions.

The third job requested a 10-second lifetime. The original record was expired
under that restriction, so the job published
`96691b5f99f7d1a7b82eff8ba18369d3f8339591d253b46dc0fd9d0806e67b4a`,
observed at 1790084214.9924803, with coefficient 1.9630619883538155e-10.
The follow-up's 300-second request could not extend this new record's producer
lifetime of 10 seconds. It checked the still-eligible original record instead.
The original artifact remained byte-identical. The live runs exercise checked
reuse and expiry refresh; out-of-band drift refresh is covered by controlled
tests, not claimed as an observed live event.

Every forecast used the declared coefficient in 12 positive source-derived
copy-service calculations, 48 checks total. All eight fixed/forecast executions
wrote all 4,097 variants exactly once in 33 fsynced parts and had zero statistic
difference from their control. All decisions kept chunk size 128 because the
candidate marginal-cost stability check failed. JAGWAS retained all 512 traits
on each GPU throughout.

| Job | Shared preparation/admission | Fixed / forecast execution | Fixed / forecast first written part |
| --- | ---: | ---: | ---: |
| First | 2.396 s | 1.934 / 2.287 s | 1.488 / 0.182 s |
| Reuse | 1.730 s | 1.642 / 1.680 s | 1.266 / 0.131 s |
| Shorter lifetime | 1.207 s | 0.922 / 0.850 s | 0.776 / 0.053 s |
| Follow-up | 3.296 s | 1.661 / 1.997 s | 1.421 / 0.064 s |

Execution includes writing and excludes shared preparation/admission and
imports. Fixed controls precede forecast executions and have different warm
states; these values support no speedup comparison. The tuning cost is large
relative to these deliberately small jobs. Useful production behavior still
requires qualified component prices and a cost-aware decision to skip or defer
work when it cannot repay that cost.

Reports are `results/productive_drift_{first,reuse,expiry,followup}_v1_20260922/report.json`.
The first job ran in `20260922-083405-968430`, reuse and expiry in
`20260922-083616-968940`, and follow-up in `20260922-083810-969960`.
Verification job `20260922-083940-970069` matched all 135 package-source hashes,
the benchmark and helper hashes, all five inputs, exact record contents and
file hashes. It checked callback sample prefixes, original ages, coefficient
consumption, cumulative cost and the eight completed outputs. Its report is
`results/productive_drift_verification_v1_20260922.json`, produced by
`benchmarks/verify_productive_drift_20260922.py`. The selectively pulled copy is
`results/productive_drift_verification_v1_20260922/report.json`.
