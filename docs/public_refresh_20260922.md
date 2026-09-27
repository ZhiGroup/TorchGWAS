# Refreshing a saved copy coefficient during a public GWAS run

The optional public initial-chunk tuner can now collect or check an independent
resident-copy CPU coefficient after useful output. An expired named coefficient
does not block memory-admitted scientific work. It remains ineligible for timing
decisions until a bounded check or measurement window completes.

This extends [public initial-chunk tuning](public_initial_chunks_20260922.md) and
the [CPU refresh controller](cpu_service_refresh_20260922.md). It does not collect
GPU, transfer, decode, allocator, storage or reduction-service rates. Those
supplied coefficients and their required evidence remain unchanged.

## Configuration and admission

Add this object to `autotune_config["initial_chunks"]`:

```python
autotune_config["initial_chunks"]["resident_copy_refresh"] = {
    "cache_dir": "/data/project/calibration/parameters",
    "profile_dir": "/data/project/calibration/profiles",
    "binding_indexes": [0],
    "max_age_seconds": 3600.0,
    "expected_cpu_seconds": 0.04,
    "expected_wall_seconds": 0.1,
}
```

The indexes refer to declared `price_bindings` in the supplied immutable
profile. Every target in a selected binding must be an existing
`owned_result_copy_scenario.resident_cpu_seconds_per_byte` or
`process_units.numpy_copy_bytes` field of an explicit device profile. Other
fields, loaded-stage observations, live capacities and mixed unsupported targets
are rejected. The two expected costs are caller forecasts, not measured defaults.

Startup still checks source, hardware/runtime context, artifact bytes, exact
target values, all unselected measurement ages and the mode-specific reduction
artifact. Only the selected bindings may await an age/drift check. Their status
is `pending_refresh`, not verified. Normal calculator validation still rejects
their expired prices. Memory admission adds 32 MiB for the two 8 MiB probe arrays,
the equality-check temporary and probe bookkeeping while preserving the user's
existing reserve. There is no probe allocation, sample or cache lookup at startup.

## Useful work, measurements and publication

Each completed writer callback can perform at most one charged planning step.
The copy probe uses four `numpy.copyto` operations per sample with pre-touched
uint8 arrays. A callback takes at most two samples. Thread CPU time prices the
fixed copy work; the controller separately charges buffer allocation, validation,
cache I/O, context checking and publication wall/CPU time. The coefficient is not
inferred from a loaded scan interval or identified as DRAM bandwidth.

A compatible saved record receives two check samples. A consistent result keeps
the original coefficient, record hash, observation time, lifetime and profile
binding timestamp. It does not create a younger measurement. A missing, expired
or drifting record requires seven samples, collected as 2, 2, 2 and 1 across four
callbacks. Drift and spread use the existing declared heuristics: a twofold rate
ratio beyond the absolute CPU tolerance, and at most a fourfold sample spread.
These are heuristics rather than statistical confidence bounds.
The later [temporal-window check](cpu_window_stability_20260922.md) additionally
requires consistency between early and recent samples before reuse/publication.

The live/source checks run before and after sampling, including a guard before
any measurement is published. A context change detected after the last sample
cannot publish a new record or profile. A complete accepted measurement creates
a new immutable record and a new content-addressed profile. Selected coefficients
in future model windows are explicitly rebound to that profile. Old evidence
and caller objects are never overwritten. Unneeded predecessor probe-artifact
dependencies are removed from the derived profile; the original files remain.

Measurements and runtime comparisons use the same cumulative tuning budget.
Four refresh callbacks plus one candidate comparison require at least five
admitted steps. A probe does not change chunk size; a later comparison must repay
all charged steps plus declared switching, publication and reserve costs. A
failed/noisy/late probe or an exhausted budget retains the current size and
allows scientific output to finish. Partial windows are not published, and probe
buffers are released on completion, failure, budget stop or API exit.

The audit is in `run.json` under
`autotune.productive.resident_copy_refresh`. It contains original record ages,
the callback/sample sequence, fresh checks, publication identity and buffer
cleanup. `autotune.active_profile_sha256` distinguishes a derived profile from
the original admission profile. A caller can load the published profile path in
a later job with the same refresh configuration.

## Validation and limitations

Controlled tests cover deferred expiry without weakened integrity checks,
two-sample reuse, drift, expiry, partial/failed measurements, context changes
before publication, charged callbacks and unchanged source coverage. Real
writer-failure tests also exposed and fixed a peer-cancellation race in both
dense tile and variant-shard writers: the first originating failure is now
recorded before cancellation, so a faster cancelled peer cannot hide it.

The first real audit used N=2,049, M=16,385, K=512, two A100 GPUs and server-local
native hard-call PGEN, with JAGWAS output. A first job collected seven samples over
four useful callbacks. The next job's stable check differed by 2.141x and produced
a new record; both GWAS outputs were bit-identical to the explicit-size control.
That audit's original test expectation incorrectly demanded reuse and stopped
after recording the valid drift response. Its evidence remains in
`results/public_refresh_20260922/execution/report.json`.

A subsequent audit rejected a window whose samples ranged from 7.93 to 38.51 ms
of CPU time, exceeding the fourfold spread limit. It published no new profile,
kept chunk size 128 and still completed identical GWAS output. That evidence is
in `results/public_refresh_20260922/execution_v2/report.json`. The audit driver
now records both reuse and justified refresh/rejection instead of demanding one
outcome from a shared machine. It never widens thresholds to make a run pass.

The final fixed five-run audit completed collection, consistent reuse, a later
drift replacement and expiry replacement. Every run emitted all 16,385 joint
rows with zero numerical difference from the explicit-size control; all prior
record/profile bytes remained unchanged. The measured record used by `check`
kept exactly the hash and observation time from `refresh`. `check_again` detected
a fresh/cached ratio of 0.281 and published a replacement; `expired` waited until
that replacement's original lifetime ended and collected a new window.

| Run | Samples | Probe CPU time | Total planning wall time | Result |
| --- | ---: | ---: | ---: | --- |
| refresh | 7 | 0.119 s | 5.357 s | New record after four callbacks |
| check | 2 | 0.038 s | 2.330 s | Same original record and timestamp |
| check_again | 7 | 0.028 s | 2.982 s | New record after detected drift |
| expired | 7 | 0.098 s | 4.402 s | New record after expiry |

All decisions retained size 128. The scope includes the complete metadata,
validation, refresh and model cost; it is not a speedup experiment. The extra
time is mostly outside the copy samples and motivates reducing repeated
validation and calculation work before a performance claim. The control took
3.563 s; output-inclusive API times of the four tuned calls were 9.607, 4.174,
5.618 and 6.354 s under differing warm/load conditions.

Final evidence: `results/public_refresh_20260922/execution_v3/report.json`.
All 137 package source hashes and the benchmark hash matched the executed
version. The final targeted regression run passed 92 tests, including the new
pre-publication context guard and both deterministic peer-error tests; earlier
broader validation results are retained beside it.

The independent resident-copy coefficient is the only newly measured price in
these execution audits; the other prices are explicit synthetic controls.
Realistic component qualification, calculator accuracy, and profitable tuning
on large production jobs are not established by these lifecycle tests.

The subsequent [productive digest reuse audit](productive_digest_reuse_20260922.md)
addresses repeated source/library hashing while preserving live checks and
original measurement ages.
