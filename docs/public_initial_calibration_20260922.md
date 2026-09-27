# Productive measurement reuse through the public API

`run_linear_gwas` accepts an optional `initial_calibration` dictionary. It
collects observations from useful early chunks and writes immutable evidence
after the association output and required run/QC sidecars succeed. A later run
checks matching evidence against a short new productive sample. There is no
separate benchmark run or candidate search in this option.

```python
from torchgwas import run_linear_gwas

result = run_linear_gwas(
    "/data/cohort/input.pgen",
    "/data/cohort/phenotype.npy",
    "/data/cohort/covariates.npy",
    pgen_mode="hardcall", device="cuda:0", chunk_size=1024,
    reader_workers=4, prefetch_chunks=4,
    reduce="significant", significance_threshold=5e-8,
    output_dir="/data/results/job_2",
    initial_calibration={
        "cache_dir": "/data/calibration/torchgwas",
        "max_age_seconds": 300,
    },
)
audit = result.run_metadata["initial_calibration"]
```

This entry point currently requires an explicit native hardcall `.pgen` file,
aligned arrays or array-file inputs, CUDA FP32, an explicit positive chunk
size, and binary output. It supports full output, phenotype-tiled full and
significant output, and variant-sharded full/JAGWAS output. JAGWAS still requires
the entire phenotype panel on each active GPU. It cannot be combined with the
existing upfront `pipeline_profile` or `autotune_profile` selectors.

The supplied configuration stays fixed. Collection does not certify memory
admission, choose a chunk/tile/device layout, or estimate independent hardware
capacities from overlapping pipeline spans. Automatic budgeted adaptation
remains separate work.

## Reuse and freshness

The shared parameter cache distinguishes three kinds of evidence:

| Kind | Reuse rule |
| --- | --- |
| Structural parameters and work/geometry | Matching declared dependencies; no elapsed-time expiry |
| Empirical capacities and stage observations | Matching dependencies and an age limit measured from the original observation |
| Available memory and contention | Fresh live checks; historical observations are audit evidence only |

This public option creates stage-observation records. They contain read/decode
wall intervals and optional CUDA event spans, including device reduction work
where applicable. CUDA spans can contain scheduling gaps. Consumer suspension
does not establish host selection service, writer service, or durable-output
completion. The output mode and thresholds are dependencies and never change
as part of measurement collection.

Each record binds package source hashes, numerical/native libraries, GPU and
host execution context, input identities, the output request, and exact variant
and phenotype ranges. Only the first partition on each GPU is observed. Later
phenotype tiles do not reopen its sample budget or reset the shared wall window.
Different configurations or ranges need separate matching evidence.

File identity uses resolved path, device/inode, size, mtime and ctime. This avoids
an extra complete pass through a huge phenotype panel and is not a content-hash
guarantee against changes hidden by filesystem metadata. Input/source changes
detected before publication prevent saving the baseline. In-memory inputs use
a job-specific identity and deliberately cannot reuse measurements across jobs.

Records are content-addressed and never overwritten. A consistent cache hit
writes a separate validation record and retains the original baseline's age.
Detected drift expands the new sample within its original budget and, if it
completes, writes a new baseline. An incomplete/failed run publishes no baseline.
Optional cache I/O failures are reported without discarding completed results.

Stable source/library digests can also be reused with fresh filesystem identity
checks. Current devices, settings and mounts are still read each job. Set
`reuse_binding_digests=False` in this dictionary to force full byte hashing;
see [digest reuse](binding_digest_reuse_20260922.md) for validity, publication
and measured cost. This structural optimization does not renew empirical age.

## Budgets and timing boundaries

Defaults are two validation samples per GPU on a cache hit, at most eight samples
per GPU on a miss or drift, one warmup chunk, sampling every fourth subsequent
chunk, a shared ten-second reservation window, and a 50 ms decision CPU budget
per GPU. Already reserved reads finish normally; the wall window bounds new
measurement reservations rather than forcibly interrupting CUDA or file I/O.
CUDA event measurement can be disabled with `cuda_events=False`.

The controller CPU budget excludes context binding, cache I/O and CUDA event
costs. The public report records `binding_seconds`, `cache_lookup_seconds`, and
`finalization_seconds` separately. Identity binding precedes the scan and can be
material for a short job. No claim of zero startup cost is made. The existing
`runtime_seconds` boundary precedes final cache publication; measure the complete
API call when comparing end-to-end overhead.

The public report also includes `output_progress`, using actual writer
completion notifications. Dense output reports the first written beta/t prefix;
df flushing and final fsync can follow it. Reduced output separately reports the
first indexed-part fsync. See [dense productive output](dense_productive_output_20260922.md)
for the boundaries and tiled/sharded coverage. These per-run timestamps do not
become independent writer capacities or renew cached observation ages.

The numerical validation and drift heuristic are documented in
[initial chunk calibration](initial_chunk_calibration_20260922.md). The separate
[productive-run bridge](productive_run_20260922.md) and
[analytical proposal bounds](jit_proposal_20260922.md) do not yet turn this option
into an automatic public JIT policy.

## Validation on 2026-09-22

Job `20260922-014419-847079` passed all 140 targeted tests in 26.47 seconds,
covering dependency mismatch, expiry, immutable reuse, shared-window limits,
cache failures, early rejection, and existing dense/significant/JAGWAS planning
contracts. The test log is
`results/public_initial_calibration_v1_20260922/tests.log`.
Job `20260922-015435-851236` additionally passed all 80 adaptive-chunk,
multi-GPU scheduling and productive-run execution tests in 15.57 seconds;
see `execution_tests.log` in the same directory. Total: 220 passed, no failures
or skips in these two targeted suites.

Job `20260922-014649-847780` ran twelve fresh processes: an uninstrumented
control, initial collection, and later cache lookup in each of four modes. The
native PGEN fixture had N=2,049, M=4,097, K=512 and two covariates. Genotype and
metadata, phenotype/covariate inputs, association outputs and the measurement
cache were on local `/data` XFS `/dev/md0`. Reports remained in the shared project
results directory. Chunks were 128 variants, prefetch depth two, and the two
tiled modes used 128-trait tiles on CUDA devices 0 and 1. JAGWAS sharded variants
across those devices and retained all 512 traits on each.

| Output | Control API seconds | First collection | Later cache lookup |
| --- | ---: | ---: | ---: |
| Full, one GPU | 2.185 | 3.467 | 2.501 |
| Full, phenotype tiles, two GPUs | 4.786 | 6.091 | 3.991 |
| JAGWAS, variant shards, two GPUs | 4.655 | 5.488 | 4.617 |
| Significant pairs, phenotype tiles, two GPUs | 4.762 | 5.507 | 3.637 |

All corresponding numerical arrays and selected pair identities were identical,
with maximum absolute difference zero. Every later process found the expected
baseline record and all original record bytes remained unchanged. All seven
per-GPU validation windows detected interval drift and expanded to the audit's
three-sample maximum, publishing new records. The consistent-hit path is covered
by deterministic tests; this hardware audit did not observe that path.

The changing intervals include H2D ratios from 0.07 to 6.16 and significant-mode
result-transfer ratios above 40. These loaded event spans cannot identify a
corresponding hardware-capacity change. They must remain separate from reusable
component capacities and transient-state estimates.

Context binding took 0.63–1.64 seconds. First fsynced indexed parts appeared at
4.213/5.131/4.425 seconds for JAGWAS and 3.163/4.003/1.898 seconds for significant
pairs, in control/first/later order. Dense first-part durability was not
instrumented. These are one run per condition, with uncontrolled warming/load,
and establish lifecycle correctness rather than speedup or an overhead bound.

The full report and the exact audit script used for these twelve processes are
in `results/public_initial_calibration_execution_v1_20260922`. Its 124 package
source hashes and archived script hash were verified against the local source
and pulled artifact. The current benchmark additionally exercises injected
writer and metadata failures.

Job `20260922-015251-850660` verified a writer failure after twelve completed
JAGWAS parts and a required metadata-write failure after all 33 parts. Neither
published a calibration record, and every prior cache record retained its exact
bytes. The two reports in
`results/public_initial_calibration_failures_v1_20260922` bind the same 124
package source files and the updated benchmark script.

Subsequent startup-admission changes have a separate source-bound audit in
[adaptive startup](adaptive_start_20260922.md). The 124-source-file results
above describe their recorded revision; they are not current-source calibration
artifacts after those changes.

The later [context-binding audit](context_binding_20260922.md) measures and
reduces GPU identity-query overhead while retaining content hashes, live
settings checks, immutable records and original observation ages. Its timings
have a separate source-bound report and do not revise the historical results
above.
