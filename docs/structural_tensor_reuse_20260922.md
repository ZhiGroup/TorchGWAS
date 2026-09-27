# Reusing tensor-work calculations across jobs

`StructuralTensorWorkCache` saves the duration-free operation ledgers used by
the shared scan calculator. These include eager statistics and JAGWAS factor
preparation/projection traces. They contain shapes, dtypes, operation counts,
logical storage accesses and source provenance. They contain no measured
association durations, service rates or live resource availability.

Records use `CalibrationParameterCache`'s immutable `source_work` artifacts.
Dependencies include the relevant package source hashes, Python version,
PyTorch version/build, CUDA version, default dtype, exact typed shape/mode
arguments, and tracing-function implementations/defaults. Changes cause a
miss. Dynamic wrappers, closures or nonportable requests bypass persistent
reuse. The implementation fingerprint uses explicit bytecode, constants and
signature data: raw Python marshal output proved unstable after executing a
nested trace and was replaced following a failing regression test.

The calculator still prices loaded work with its current supplied independent
service parameters and compiled geometry. The empirical records behind those
prices retain their separate validity and expiry requirements. Structural
reuse cannot renew their measurement times or certify a hardware capacity.
Survivor counts, source-window observations and live memory/contention are
not cached as tensor work.

## Lifecycle and overhead

The optional cache is explicit and scoped:

```python
cache = StructuralTensorWorkCache(cache_directory)
try:
    with cache.activate():
        comparison = compare_prepared_windows(baseline, candidate, **common)
    cache.publish(successful=True)
finally:
    cache.close()
```

The caller must charge construction and lookup to its calculation budget.
New ledgers are staged, without writing files inside the calculation. Publish
only after successful work, and include its cost in the output-inclusive run
accounting. A loaded record is never republished, touched to renew its age, or
overwritten. Source changes prevent publication. Cache I/O failure permits
ordinary computation and is reported rather than invalidating GWAS results.

Retained state defaults to at most 16 entries and an 8 MiB Python memory
estimate. Entries hold immutable JSON, which gives each consumer an independent
mutable ledger on decoding. This avoids repeated deep copies and recursive
sizing of thousands of loaded objects. In-process priced-component caching
remains separate and still binds all supplied service and geometry inputs.
These retention limits do not constitute a preemptive filesystem-I/O deadline.

`ProductiveTuningRun(..., structural_cache_dir=...)` owns this lifecycle for
incremental planning. It constructs and activates the cache only inside an
admitted planning step after a typed writer-completion event. Actual
initialization/loading CPU and wall time are charged by the existing planning
budget. An initialization overrun cannot authorize a chunk-size change.
New records remain bounded in memory after the tuning window closes, until
successful run completion publishes them. Failed runs discard staged work.
Publication CPU/wall time and errors appear in the final audit. The caller's
cost forecast must include this deferred work and reserve the cache's memory.

This adds reuse to the internal productive controller, not a complete public
automatic policy. `initial_calibration` continues to collect and validate
bounded early-chunk observations without selecting a new layout itself.
Remaining-work forecasts, live continuation/output state and prediction
qualification still govern whether a proposed change is worthwhile.

## Verification scope

Tests cover cross-job reuse without age refresh, changed source/library/dtype/
shape/mode/implementation dependencies, mutation isolation, changed price
recalculation, bounded retention, corruption, failed publication, source drift,
and controller initialization/finalization. Existing source-count, geometry,
memory, window-composition and planning-budget regressions are also exercised.

`direct_structural_tensor_reuse_20260922.py` starts separate control, populate
and reuse processes. Each compares chunk sizes 128 and 512 over the same
prepared windows in five modes: dense, empty/sparse/dense significant output,
and JAGWAS. Four structural records cover two statistics and two projection
shapes. The audit requires all predictions to remain exactly equal and the
reuse process to leave the saved artifacts byte-for-byte unchanged.

`direct_structural_tensor_profile_20260922.py` alternates four uncached/cached
pairs in one process and separately profiles a cached comparison. Cache
construction is included in these paired times. This helps distinguish the
cost of the implementation from variation between separate processes. The
profiler's instrumented elapsed times are not performance measurements.

Both audits use real PGEN headers and captured launch geometry but synthetic
service prices. They measure calculator overhead, not GWAS runtime, selection
accuracy or achieved speedup. Fixture construction, header indexing and process
imports are outside the comparison timer; initialization and publication costs
are explicitly reported by the separate-process audit.

## Final remote evidence

Job `20260922-052226-923392` completed with **221 tests passing in 65.93 s**.
The test log is `results/structural_tensor_cache_v4_20260922/tests.log`.
The separate-process audit is
`results/structural_tensor_reuse_v4_20260922/report.json`; paired timings and
diagnostic profiles are in `results/structural_tensor_profile_v3_20260922`.
The report's 130 package files, five helpers/geometry files and benchmark hash
matched the tested revision. Input and cache artifacts reside on `/data`, XFS
`/dev/md0`; project source remains on the existing NFS checkout.

All five mode comparisons were exactly equal across control, population and
reuse processes. The reuse process loaded all four structural records, missed
none, published none, and left all artifact contents unchanged. Its retained
trace-cache estimate was 185,574 bytes. Changed independent GPU prices were
separately tested to force repricing while reusing the same source work.

Timing did **not** establish a reliable speedup or compliance with the default
50 ms total planning CPU budget. In the four alternating matched pairs, control
CPU time ranged from 178.71 to 403.98 ms and reuse from 180.68 to 255.77 ms;
reuse was faster in two pairs and slower in two. First control wall time was
837.39 ms, demonstrating that early costs cannot be omitted from accounting.
The separate-process mode comparisons also varied substantially. These are
observations under the measured server conditions, not hardware bounds.

Cache initialization cost 3.06 ms CPU / 29.17 ms wall for population and
2.58 ms CPU / 20.68 ms wall for reuse. Publishing four new records cost
33.13 ms CPU / 57.32 ms wall. The reuse process's no-new-record publication
check cost 1.86 ms CPU / 1.98 ms wall. These costs must be added where they sit
outside a comparison's own timer; caching must not hide them after the GWAS
timer or reset the productive tuning window.

The implementation removes repeated Python-object traversal from structural
reuse, but graph composition, priced-component copying, pricing itself and
live binding remain. Whether to spend more than the conservative default
planning budget must follow a worthwhile remaining-work forecast. This patch
does not silently increase that budget or enable an unqualified automatic
switch. The next policy work must account for remaining work, already-issued
ranges, output continuation, planning/publication cost and uncertainty together.

The follow-up `window_forecast_20260922.md` adds bounded remaining-work
scenarios and requires an explicit publication forecast in the productive
controller when persistence is enabled. Its decision cost includes previous
planning steps; the cache no longer relies on the caller folding publication
into an unspecified switching estimate.
