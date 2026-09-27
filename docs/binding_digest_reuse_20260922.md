# Reusing source and library digests across jobs

Initial calibration now optionally reuses SHA-256 digests of stable package and
loaded numerical/native-library files. Each job still reads its current device
identity, driver, CPU topology, thread pools, execution flags and storage mounts.
Digest reuse changes the cost of constructing the dependency binding; it does
not make empirical measurements permanent or replace live resource admission.

The option is enabled with `initial_calibration` and can be disabled explicitly:

```python
initial_calibration={
    "cache_dir": "/data/calibration/torchgwas",
    "reuse_binding_digests": False,  # Hash source/library bytes on every binding.
}
```

Standalone `source_identity()` and `execution_context()` continue to hash full
file contents unless the caller supplies a `BindingDigestCache` object.

## Validity and publication

Records live below `cache_dir/binding-digests-v1`, use the structural
`source_work` kind, and carry their actual SHA-256 values. Dependencies include
the host and boot identities plus every file's resolved path, device, inode,
size, mtime and ctime. The same group can be requested in a different order.
Each lookup checks file identities before and after use. End-of-job checks also
retain the original requested paths, so changing a symlink target invalidates
the binding even when its previous target remains unchanged.

Files with recent (less than one second old) or future change timestamps cannot
use a saved digest. Their bytes are read, and the first digest is retained for
a final content check. This covers same-size edits within one filesystem clock
tick. Ordinary content changes, file replacements, source-set changes and host
reboots prevent reuse. Like the existing input metadata binding, this trusts
filesystem change metadata; it is not authentication against deliberately
concealed changes or a guarantee about remote filesystem cache coherence.

A cache instance retains at most eight groups and 512 files/requested paths.
New entries remain in memory until the scientific job succeeds and its bindings
still match. Failed jobs publish nothing. Cache reads/writes tolerate missing,
corrupt or unavailable records; failed lookups read file bytes normally.
Publication never overwrites a previous record, and a hit does not write a new
record or renew any timestamp. Filesystem operations are not forcibly bounded
in wall time, and these retention limits are not an RSS guarantee.

The public audit includes `initial_calibration.binding_digests`: hit/miss and
byte-read counters, publication results, errors, and publication wall/CPU time.
Its cache is closed before the report, so retained-entry counts are zero there.
Binding and finalization are still charged to the full API call. This cache
does not change empirical lifetimes, the early useful-chunk sampling budget,
output thresholds, or the JAGWAS full-phenotype requirement.

## Why this change

The component profile in `results/binding_profile_v1_20260922/report.json`
measured three bindings with the public API libraries already imported. At that
revision, each binding hashed 131 package files (1,831,300 bytes) and four
library files (61,279,619 bytes). Library hashing took 0.147–0.266 seconds;
source hashing took 0.029–0.420 seconds. The warmed cProfile attribution placed
0.189 of 0.209 seconds in file hashing. Those instrumented observations identify
work worth removing; they are not a GWAS throughput measurement.

Source and reports reside on project NFS; scientific inputs, output and the
calibration cache use server-local `/data`. Source-binding timings must not be
interpreted as a local-storage throughput ceiling.

## Verification and binding cost

Remote job `20260922-065406-950360` passed **165 tests in 25.55 seconds**.
Log: `results/binding_digests_v2_20260922/tests.log`. Tests cover ordinary edits,
restored mtime, file replacement, host/boot changes, recent/coarse timestamps,
symlink changes during and after hashing, immutable cross-job reuse, corruption,
I/O failures, bounded retention, failed-job publication refusal, strict opt-out,
and fresh device/thread/storage context. The coarse-timestamp test reproduces a
gap caught in the first remote test run; completion now checks the original
content for initially recent files instead of relying only on their stat data.

`results/binding_digest_reuse_v1_20260922/report.json` compares the full binding
against cached file digests with public API libraries imported in both arms.
The population process hashed 136 files (132 source files plus four libraries),
63,121,347 bytes in total, and published three groups. The later process loaded
all three groups, hashed no source/library bytes, and published no replacement.
All complete source hashes and execution contexts matched exactly. The saved
records remained byte-identical through reuse and six additional warmed pairs.

Each warmed pair alternated execution order. The timing includes construction,
source/context binding, final source/file checks and deferred cache publication;
imports and scientific work are excluded. CPU is calling-thread CPU time.

| Six warmed repetitions | Full byte hashing | Digest reuse |
| --- | ---: | ---: |
| Median wall seconds | 0.376 | 0.232 |
| Wall range, seconds | 0.249–0.675 | 0.141–0.265 |
| Median CPU seconds | 0.331 | 0.171 |

All six paired wall differences favored reuse, by 0.015–0.534 seconds. This is
a measured reduction in binding work on this host, not a GWAS speedup or a
guaranteed overhead bound. The three first bindings in separate processes were
1.901 s for control, 0.584 s for population and 0.800 s for reuse; these single,
ordered observations include CUDA initialization and do not isolate a cold
startup effect. Fresh metadata/mount/topology checks still incur measurable work.

The report verifies `/data` input/cache mounts as XFS `/dev/md0` and the package
mount as NFS4. All 132 delivered package hashes and the benchmark hash were
checked against the report after pulling it.

The same job completed the public API audit in
`results/indexed_output_partition_execution_v2_20260922/report.json`: twelve
fresh processes plus one injected writer-callback failure. The fixture uses
2,049 samples, the nonzero source interval [129, 1154), 512 phenotypes and
128-variant chunks. Scientific input/output and cache are on local `/data`.

| Output configuration | Result rows | Later-run output-count record |
| --- | ---: | --- |
| Significant pairs, four 128-trait tiles on two GPUs | 5,264 | Original record, age unchanged |
| Empty significant pairs, same tiling | 0 | Original record, age unchanged |
| Significant pairs, complete panel on one GPU | 5,264 | Original record, age unchanged |
| JAGWAS, two variant shards with all traits per GPU | 1,025 | Original record, age unchanged |

Pair/variant identities and significant-pair df are identical to their controls;
statistics match within the existing FP32 tolerance (rtol 6e-5, atol 3e-4).
Producer/source coordinates, complete partition coverage and empty output were
checked against actual parts. Every later run loaded all three digest groups,
read no source/library bytes for binding, published no new digest record, and
preserved prior cache artifact hashes. The failure case left every cache record
unchanged and closed all torchGWAS worker threads. All 132 package hashes and
the public-audit script hash match the pulled report.

The public report includes first completed output and output-inclusive API wall
time, excluding imports. Single-run API times ranged from 3.10 to 10.30 seconds
and are too variable to infer whole-job speedup. Binding ranged from 0.068 to
0.824 seconds. Existing read-only NumPy and prefetch-limited decoder warnings
appeared; input hashes were unchanged. This completes verification of optional
digest reuse, not runtime-prediction qualification or public automatic tuning.
