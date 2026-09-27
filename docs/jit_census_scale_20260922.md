# Large-source preparation and JIT extrapolation gate

The read-only frozen H100 candidate with 35,365 samples, 1,048,576 variants
and 16,385 phenotypes took 77.376 wall seconds under `cProfile` to construct.
The earlier uninstrumented preparation took 52.292 seconds in a separate
capacity-sensitivity run; these are not matched timing comparisons. The new
profile attributes 77.307 instrumented seconds to the exact PGEN payload
census, including 1,048,806 record visits and 831,330 per-record `np.unique`
calls. Candidate enumeration itself was negligible. The diagnostic script and
frozen source/artifacts were checked before and after execution; the frozen
package was not edited. Evidence is under the diagnostic project's
`results/frozen_candidate_prepare_profile_20260922/`.

The development census now computes end positions once, uses `np.diff` with
`prepend`, and counts validated 1–5-byte varint lengths with `np.bincount`.
It removes a duplicate predicate pass, an `np.r_` temporary and a sort/unique
for every difflist record. On the same 65,536-variant H100 prefix, three old
and three new observations produced exactly the same full structural report
SHA-256 (`728efd281664773437d8ddafa0e7da9050277aebce9805741f1e82de399ab8f9`),
including 64 chunk reports. Median process CPU was 2.027 seconds old and
1.370 seconds new, a 32% reduction in this bounded control. The schedule and
original observation times are in
`results/census_numpy_calls_20260922/report.json` in the diagnostic project.
This is a census improvement, not an end-to-end GWAS speedup or a replacement
for deferred source-window modeling.

The current JIT path already avoids the full payload census before output.
It compares three short header-only horizons after output and imposes a
caller-declared maximum extrapolation ratio. A large source can exceed that
limit even when the largest modeled horizon fits every unfinished partition.
For illustration, 8,086,101 remaining variants divided by a 4,096-marker
horizon is about 1,974; a cap of 8 cannot support a remaining-job forecast.
The public controller now detects this exact issue-frontier condition before
switch-proposal price validation or graph construction, records the ratio, and continues the
scientific scan at its admitted chunk size. Scheduled independent component
refresh still runs first so it may produce reusable evidence. This gate avoids
a knowingly futile synchronous callback; it does not make long-range JIT
tuning possible.

The next architectural requirement is a scale-aware continuation model:
aggregate bounded source evidence over the first chunks or validated header
metadata, check source/output regime changes, and combine it with independent
component prices and finite fill/drain work. Simply raising the extrapolation
cap while keeping a fixed short-window error allowance would not establish
predictive accuracy. Huge-phenotype tile changes and live GPU reassignment
also still require an executor/output handoff design.

The updated A100 development source passed 319 targeted tests in job
`20260922-205915-1207707`, covering PGEN census/range views, candidate
construction, header bounds, the public first-chunk controller and exact
issue-frontier behavior. The older `tests/test_pgen.py` could not collect
because its imported `benchmarks.benchmark_pgen_converter` file is absent from
this checkout; it was excluded from the passing batch. This limitation is
separate from the census change.
