# Reusing file digests during productive tuning

Public initial-chunk tuning now uses the existing bounded `BindingDigestCache`
for repeated source/library identity checks. It still reads all source and
library bytes during startup admission. Cache initialization and the first
productive lookup occur after useful output, inside the charged planning step.

This addresses the repeated validation work observed in the
[public refresh audit](public_refresh_20260922.md). It changes neither the
resident-copy measurement protocol nor the calculator's physical resource model.

## Configuration and validity

Digest reuse defaults to enabled when `initial_chunks.resident_copy_refresh`
provides a `cache_dir`, or when `initial_chunks.structural_cache_dir` is set.
An explicit `initial_chunks.binding_digest_cache_dir` overrides that location.
With none of those directories, the tuner keeps its original full-byte checks.
Set `initial_chunks.reuse_binding_digests = False` to require full-byte checks
throughout a job. These options affect only productive callbacks; startup
admission continues its existing full-byte validation.

Each lookup rechecks requested/resolved path, device, inode, byte count, mtime and
ctime before and after using a digest. Host and boot identity bind saved records.
Files changed within the previous second bypass reuse, and their original bytes
are rechecked before using the result or publishing a cache entry. Every
productive validation additionally rechecks all tracked identities. The existing
cache bounds remain eight groups and 512 files per job. This trusts filesystem
change metadata; it does not authenticate against adversarial metadata changes.

Only file content hashes are cached. Runtime thread counts, affinity, flags,
device/driver identity, mounts, currently available CPU/GPU memory, original
measurement age, declared coefficient values and component artifact contents
are still checked. A changed/expired profile stops optional tuning while
scientific work retains its admitted chunk size. An unavailable optional digest
cache falls back to reading file bytes at initialization or lookup; exceeding
its bounds or detecting changed files stops the optional tuning step.

New digest records are published only after successful scientific completion and
another file-identity check. Failure, changed files or a failed publication never
overwrites a record. Publication errors do not invalidate completed scientific
output. All digest state is closed on API success or failure. A digest cache hit
does not renew any empirical measurement's observation time or lifetime.

## Cost accounting and validation

`autotune.productive.validation` reports the number and cumulative wall/CPU cost
of productive checks, including unsuccessful checks. These costs are already
inside `productive.planning`, not additional time to add to that total.
`autotune.productive.binding_digests` reports cache hits, byte hashes, publication,
cleanup, and measured final publication wall/CPU time. Final cache publication
occurs after the productive steps, so it is reported separately; the caller's
explicit `publication_seconds` forecast remains part of every switch payback
decision. Startup binding/admission cost remains separately reported.

The targeted batch passed 110 tests. Coverage includes cross-job digest reuse,
unchanged measurement bytes and age, disabled reuse, fresh source/library/input
checks, live memory/settings changes, price expiry and artifact changes,
post-check file changes, failure cleanup and publication errors. Existing tests
also cover same-size/restored-mtime edits, symlink retargeting, coarse/recent file
timestamps, corruption, bounded retention, and fresh physical GPU/runtime data.

The real audit uses `benchmarks/direct_public_refresh_20260922.py
--digest-comparison`: a control and initial refresh followed by three
counterbalanced pairs of full-byte and cached-digest jobs. Every job processes
the same server-local PGEN with N=2,049, M=16,385, K=512 on two A100s and compares
all JAGWAS results. Only the resident-copy CPU coefficient is measured; other
prices remain explicit synthetic controls. This audit measures tuning overhead
and output preservation, not calculator accuracy or profitable tuning.

The completed audit emitted all 16,385 rows in each of eight runs, with zero
numerical difference from the explicit-size control. All decisions retained
chunk size 128. All previously saved measurement/profile bytes were preserved.

| Pair | Full-byte checks: mean wall time per check | Cached checks: mean wall time per check | Full-byte / cached sample counts |
| --- | ---: | ---: | ---: |
| 1 (full-byte first) | 273.2 ms | 68.0 ms | 2 / 7 |
| 2 (cached first) | 258.0 ms | 138.0 ms | 2 / 7 |
| 3 (full-byte first) | 226.5 ms | 72.0 ms | 2 / 7 |

The median of these per-run check means was 258.0 ms with full-byte validation
and 72.0 ms with cached digests. All cached comparison jobs had three disk hits,
39 in-memory hits, and zero source/library bytes hashed in productive callbacks.
The first refresh populated those three groups by hashing 141 files totaling
63,201,356 bytes. Its final digest publication cost 0.556 s; later no-new-record
publications cost 0.028–0.034 s, including final file checks.

The full-byte jobs each reused their original copy coefficient after two
samples. Each cached comparison job detected drift and collected seven samples;
therefore it performed 14 checks versus five for each full-byte job. This
variation, shared-machine load and different warm states prevent interpreting
the API durations as a matched full-GWAS speedup. Per-run productive planning
times were 2.338/2.655/2.576 s (full-byte) and 2.077/3.145/1.519 s (cached).
Expensive work beyond source/library hashing still remains.

Evidence is in `results/productive_binding_digests_20260922/execution/report.json`
and `results/productive_binding_digests_20260922/tests.log`. All 137 package
source hashes and the benchmark hash matched the executed version. The report
retains all outcomes and the original coefficient/profile identities; the
comparison does not loosen the existing drift or noisy-window thresholds.
