# Original measurement ages in calculator decisions

Detailed calibration previously checked source, execution context and artifact
hashes, but those checks did not establish the age of most embedded prices.
The new optional v2 profile links declared calculator values to exact immutable
component records. A new profile binding time cannot renew those measurements.

## Binding contract

`bind_detailed_profile(..., price_bindings=[...])` creates
`torchgwas.detailed_calibration.v2`. Each declaration names an artifact, empirical
kind, measurement name, source/execution/protocol dependencies, an optional
shorter maximum age, and exact paths from record values to context values.
For example, one resident-copy coefficient may supply
`contexts[0].profiles['cuda:0'].owned_result_copy_scenario.resident_cpu_seconds_per_byte` and the
corresponding CPU service coefficient for another GPU's worker.

The checker requires an artifact already present in the profile's hash map,
the original content-addressed record, exact target values, matching source and
execution dependencies, and an explicit measurement protocol. It accepts CPU,
GPU, transfer or storage component records. Loaded stage observations and live
memory/contention observations cannot substitute for these records. Targets
cannot overlap. Declarations are bounded to 128 records and 512 target paths.

Freshness uses the original observation timestamp and the minimum of the
producer's lifetime and any caller restriction. Loading, rebinding, a later
publication, or a newer record in the cache cannot extend that lifetime or
replace the artifact named by an existing profile. A final clock check includes
time spent loading the remaining records. Refresh creates a new record and
requires a newly bound profile; the previous record remains unchanged.

The audit reports only `declared_targets_verified`. Unlisted prices, the
scientific validity of a measurement protocol, current contention and prediction
accuracy remain unqualified. Existing v1 profiles are readable and explicitly
report `undeclared` price freshness. This is not a claim that all existing
profiles are ready for automatic selection.

## Execution checks

The detailed autotune bridge checks declared prices at construction, after
input QC, and again after planning and live memory admission. Cached plan hits
go through the same checks. Its output audit retains the actual record hashes,
observation times, expiry and ages used by the final validation.

Geometry and setup-profile completion preserve the existing declarations.
Adding untimed geometry or new setup references cannot silently drop the
original price evidence. Any new unlisted coefficients remain unqualified.

`ProductiveTuningRun.forecast_step(..., price_profile=profile)` adds checks inside
the charged productive step. Before model construction it checks the declared
records. After constructing the proposal it verifies the profile hash, source
identity and each device's fixed profile values, then rechecks original ages.
Chunk size and untimed geometry may vary; productive evidence targets must be
fixed per-device profile parameters. A changed, expired or substituted record
stops optional tuning and leaves valid source reservations on the current size.
The existing scientific scan can continue. Caller-side admission still owns
fresh source/runtime context and available-memory validation.

## Productive measurement and later-job reuse

The real integration audit uses native PGEN, two A100s, N=2,049, M=4,097,
K=512 and two covariates. JAGWAS retains the complete phenotype panel on each
GPU. Inputs, result parts and the component cache are on XFS `/dev/md0` at
`/data`; reports are in the shared project results directory.

The first written part opens a measurement/reuse step. On a cache miss, a
bounded independent probe copies an 8 MiB resident uint8 array four times in
each of seven repeats. Its coefficient is the median thread-CPU seconds per
copied byte. It is a generic copy CPU-service measurement during the job, not
DRAM capacity, a loaded association-stage rate, or a fit to GWAS runtimes.
The fresh-copy fraction is explicitly zero in this scenario; the experiment
does not qualify that assumption or extrapolation from 8 MiB copies to small
result arrays. The active `owned_result_copy_scenario` coefficient is checked
in every modeled window: additional copy CPU time must equal the source-derived
additional byte count times the measured price.
Allocation, warmup, measurement, lookup, publication and binding costs all
remain inside the charged step.

The next written part triggers three bounded source-window comparisons at the
actual held frontier. Their price profile is checked by the productive bridge.
Cumulative tuning cost includes both steps. All other component prices remain
synthetic controls; the experiment validates evidence reuse and decision wiring,
not complete prediction accuracy or tuning profitability. It retains the earlier
explicit 10 s CPU / 30 s audit window and 0–1 s boundary-cost scenario.

The first job uses a 300 s maximum age. The next job reuses its record. A third
job requests a 10 s maximum age to exercise a real expired-cache refresh. The
output-inclusive timing and immutable-record checks are recorded with the
corresponding reports rather than inferred from microbenchmark latency.

Remote job `20260922-080500-965808` passed 241 targeted tests in 37.51 s.
The log is `results/price_binding_tests_v3_20260922/tests.log`. Coverage includes
record and value substitution, original age, shortened lifetimes, expiration
during QC/planning, cached analytical plans, productive no-switch behavior,
continued source coverage and preservation through profile extension.

The preliminary v1–v3 execution reports verified persistence and age, but put
the coefficient in `process_units.numpy_copy_bytes`, which this JAGWAS scenario
bypasses. They therefore do not establish consumption of the measurement by
the window model. The v4 audit corrects the target and checks the actual model
service. A source-level sensitivity test also verifies that tripling the active
resident coefficient triples additional copy service and increases modeled
chunk-finish service by the corresponding amount.

## Verified v4 results

Job `20260922-081224-966141` passed all nine productive price-binding tests in
15.56 s, including the new source-consumption sensitivity case. Together with
the preceding suite this covers 242 distinct targeted cases. It then completed
three real audits using the corrected active coefficient:

| Job | Evidence action | Original age at binding / decision | Measurement step | Forecast step |
| --- | --- | ---: | ---: | ---: |
| First | Measure and publish | 0.20 / 1.59 s | 0.322 s | 1.358 s |
| Reuse | Same immutable record | 26.58 / 27.81 s | 0.116 s | 1.187 s |
| Shorter lifetime | Expired; publish new record | 1.01 / 1.80 s | 1.129 s | 0.774 s |

The first two jobs used record
`f6e1806990beb3b7d719f0f99131722d6afdc8e2ec8e8062209bb090dda4907f`,
observed at Unix time 1790082791.048596. The third created
`49c2a8d2e2ab4e92c4763b70838846aa13dfc2a465a17b0f5558aa5dbc366b49`,
observed at 1790082844.8612676. The old artifact's content hash remained
unchanged. The copy coefficients were 6.360e-10 and 7.774e-10 CPU s/byte,
respectively. Freshness and exact reuse do not establish stability under load.

The active coefficient was verified in 12 model calls per job, 36 total, with
positive source-derived copy work in every call. Both declared device targets
were checked at the final decision. All six fixed/forecast executions wrote
4,097 variants exactly once in 33 fsynced parts, with zero statistic difference
from their corresponding control. All forecast decisions retained size 128
because marginal-cost stability failed. Cumulative tuning cost, including the
declared 10 ms reserve, was 1.689, 1.313 and 1.913 s.

| Job | Shared preparation/admission | Fixed / forecast execution | Fixed / forecast first written part |
| --- | ---: | ---: | ---: |
| First | 1.981 s | 2.603 / 2.385 s | 1.967 / 0.178 s |
| Reuse | 3.313 s | 2.344 / 2.025 s | 1.886 / 0.140 s |
| Refresh | 3.036 s | 1.974 / 2.345 s | 1.517 / 0.175 s |

Execution includes result writing and excludes shared preparation/admission and
process imports. Controls precede forecast executions, so their different warm
states prevent any throughput comparison. These small jobs do not justify the
measured tuning overhead; the generous audit budget tests the wiring only.

Reports are `results/productive_price_{first,reuse,refresh}_v4_20260922/report.json`.
`results/productive_price_verification_v4_20260922.json` records the immutable
record, expiry, numerical consumption and output checks. Verification matched
all 134 current package-source hashes, the benchmark/helper hashes and all five
input hashes. The audit test log is `results/price_binding_tests_v4_20260922/tests.log`.

The remaining work is to collect and qualify the other independent coefficients,
connect bounded refresh policy to the public automatic path, and validate the
calculator's output-inclusive predictions and tuning benefit. This change does
not make synthetic or unlisted prices measured.

The subsequent [productive CPU refresh](cpu_service_refresh_20260922.md) adds
fresh drift checks and distributes replacement measurement across initial
useful chunks, while retaining the original age of reused evidence.
