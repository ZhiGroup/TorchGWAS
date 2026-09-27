# Reducing productive forecast calculation cost

The productive forecast calculation now reuses bounded structural decoder
counts and constructs source windows without copying an encoded description
that it will immediately replace. Its numerical reports and decision remain
unchanged. These changes reduce repeated calculator work; they do not qualify
component prices or enable the public automatic tuning policy.

## What is reused

`PgenHeaderWork` keeps an LRU of at most 1,024 record signatures by default.
Each entry depends on sample count, record encoding form, byte length and
whether the decoder expands the record to int8. LD base-only replay is a
distinct entry from ordinary expansion. Setting `max_cached_signatures=0`
disables reuse. Each aggregate result owns its interval lists, and file identity
is still checked on every requested range. The cache contains source-count
intervals, no measured rates, elapsed times or live capacity. It lasts only as
long as the header object and does not renew any persisted measurement's age.

`prepared_source_window` copies retained data/profile fields and builds only
the requested encoded range. Each result remains independently mutable,
including its geometry bank. A memory-only admitted template must supply the
admission's `input_file_identity`; an embedded source identity is checked too.
Path, dimensions, fixed partition and window-work limits are checked. A later
replacement file cannot become the admitted input just because a new header
object was constructed for it.

## Attribution and matched replay

The initial profile replayed the actual source frontier from
`results/productive_forecast_execution_v1_20260922/report.json`: two devices,
128/256-variant chunk choices and 256/512/768-variant horizons. It used the same
synthetic component prices and asserted exact model reports, proposal and
forecast-audit equality. The profile identified repeated record-bound
calculations and container copying as substantial costs. The native event
solver was a small part of that initial profile and was left unchanged.
Attribution is in `results/forecast_cost_profile_v1_20260922/`.

The matched benchmark alternates order, constructs a new header and planning
cache for each replay, and charges all three horizons. Control disables the
signature cache and uses the earlier copy-then-replace window construction.
The reuse arm uses both changes. Source-meta initialization precedes both arms,
as it did in the original admission. This is calculator replay, not GWAS
throughput or prediction validation.

The first four pairs were noisy: control/reuse median wall time was
312/326 ms. The follow-up used 12 alternating pairs and observed garbage
collection without disabling it or changing its thresholds. All observed
pauses remain included in the reported times.

| Measure, median of 12 replays per arm | Control | Reuse |
| --- | ---: | ---: |
| Complete calculation wall time | 452.0 ms | 353.4 ms |
| Planning-thread CPU time | 433.7 ms | 337.7 ms |
| Source-window construction, all horizons | 178.3 ms | 79.4 ms |
| Actual decoder-signature calculations per replay | 2,623 | 115 |

The matched median wall reduction was 21.8%. Wall ranges were 273–782 ms
for control and 183–700 ms for reuse. Generation-2 collection took 489 ms
across two control collections and 676 ms across three reuse collections;
several slow replays contained large collection pauses. These samples do not
establish a tight latency bound or a uniform improvement on every invocation.
The cache retained 115 signatures and served 2,508 hits in each reuse replay.
Every replay retained exact numerical reports, proposal and compact audit.

Artifacts are `results/forecast_reuse_v1_20260922/` and
`results/forecast_reuse_v2_20260922/`. The latter includes per-replay collection
events and separate instrumented profiles. Script, helper, input-report and
package hashes identify each run; instrumented profile times are not the
uninstrumented benchmark results.

## Actual output-inclusive execution

Remote job `20260922-073921-960175` also ran the updated callback on two A100s
with native PGEN input: N=2,049, M=4,097, K=512 and two covariates, JAGWAS with
the full panel on each GPU. Inputs and output were on XFS `/dev/md0` at `/data`.
The fixed and forecast executions each wrote all 4,097 variants exactly once
in 33 fsynced parts, with identical statistics and all worker threads closed.

The callback cost 348 ms CPU / 356 ms wall. Its observed frontier had five
chunks issued on CUDA 0 and three on CUDA 1; source starts were 640 and 2,688.
It modeled the actual remaining ranges, not the earlier replay's four/four
split. The candidate's marginal slope changed 17.6%, above the fixed 10%
stability limit, so the controller retained chunk size 128.

Shared preparation/admission took 1.308 s. Fixed/forecast execution took
0.842/0.549 s, with first durable output at 0.699/0.048 s. Execution includes
writing and excludes shared preparation and process imports. The ordered
warm-state difference prevents a GWAS speedup comparison, including against
the older 1.314 s callback run. The live run retains synthetic prices, explicit
0–1 s boundary scenarios and the earlier generous audit-only planning budget.
Report: `results/productive_forecast_reuse_execution_v1_20260922/report.json`.

Remote job `20260922-073654-959895` passed 155 targeted tests in 26.62 s.
These cover cached versus uncached counts across every supported record form,
LD restart and sample/tail boundaries, eviction, returned-result mutation,
replacement inputs, window admission and work limits, writer composition and
productive forecast binding. Log: `results/forecast_reuse_tests_v1_20260922/tests.log`.

The remaining calculator work includes expensive component/geometry copies
and source tracing. The present callback still exceeds the default 50 ms
planning CPU budget and is too expensive for this tiny scan. This optimization
does not change that budget or relax the existing no-switch checks. Independent
price qualification, live boundary costs and the automatic policy remain open.
