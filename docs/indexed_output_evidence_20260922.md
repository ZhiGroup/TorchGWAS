# Partition-bound evidence from reduced output

Completed indexed output now carries the producer's device, complete source
variant interval and post-QC phenotype interval. Its source-file chunk range
is separate from the output file's rebased range. These coordinates are needed
to interpret survivor counts when different GPUs or phenotype tiles feed one
writer. Inferring the phenotype tile from surviving rows fails for empty
chunks and can mix distinct workloads.

## Execution path

`IndexedOutputPartition` is an immutable producer identity.
`PartitionedIndexedChunk` carries it through the significant-pairs queue while
preserving tuple-style access and sharing the existing result arrays. The tile
producer attaches the identity before publishing to the shared queue. Variant
range rebasing changes payload indices and keeps that identity. Closing the
rebasing iterator also closes its underlying source, including on writer error.

For unpartitioned significant output and variant-sharded JAGWAS, source bounds
uniquely identify the configured producer; the public API supplies that
resolver to the writer. JAGWAS metadata requires the complete phenotype panel.
Conflicting or out-of-range producer identities are rejected before a part is
written. Existing callers without partition metadata still receive ordinary
completion events, but those events cannot supply partition-bound evidence.

`IndexedChunkWrite` reports its original output-relative start/end plus
`source_variant_range` and `partition`. Empty chunks retain the same identity
even though they create no file. Material events still occur after file flush
and requested fsync. This is completed output, not queue acceptance, final
manifest/directory durability or an independent writer-capacity measurement.

## Bounded collection and cross-job reuse

The existing public `initial_calibration` option collects survivor-count bins
only after typed writer completion. It retains at most `max_chunks_per_device`
bins per GPU and 128 total, within the existing shared early-run wall window.
Collection then reduces to scalar progress bookkeeping. Later phenotype tiles
may be unobserved; the report never presents a first-tile sample as coverage
of the full phenotype panel.

Each bin contains the source chunk range, producer partition/device, phenotype
range, retained rows, physical part bytes and file-fsync state. Bounds, counts,
output mode and overlap are checked against the prepared request. An invalid
or overlapping sample is excluded from publication while completed scientific
output remains intact. No extra genotype/phenotype pass or association run is
performed to obtain these counts.

Records use the existing immutable `CalibrationParameterCache`:

- Kind: `stage_observations`; name: `initial_output_survivors.v1`.
- Dependencies: source, software/hardware execution context, file input
  identities and the complete scientific/execution request, including threshold.
- Lifetime: the configured empirical maximum age, measured from the oldest
  included output completion, rather than end-of-job publication.

Lookup is lazy, after the first eligible completed output. Its cost is reported
separately. Fresh bins are compared in canonical coordinate order, so different
GPU arrival order does not create a false change. Identical fresh/cached bins
reuse the old record without rewriting it or renewing its age. Different or
previously unobserved bins create a new measurement record and leave the old
record intact. This difference is not itself a hardware drift diagnosis.
Changed dependencies or expired observations require fresh evidence.

Publication occurs only after successful completion and unchanged input/source
bindings. An optional cache failure is reported without invalidating completed
associations. `run.json` exposes the bounded bins and reuse/publication result
under `initial_calibration.output_occupancy`.

These are empirical output observations, not immutable structural constants.
Their presence does not establish future selectivity, a writer throughput
ceiling, accuracy of the remaining-work model or permission to change a layout.
The calculator still needs an explicit occupancy scenario for unobserved work,
and public automatic configuration changes remain unfinished.

## Verification

Remote job `20260922-062420-945685` finished with **186 tests and five subtests
passing in 26.97 s**. The log is
`results/indexed_output_partition_v2_20260922/tests.log`. Coverage includes
empty/concurrent tiles, immutability, source-coordinate rebasing, invalid
producer identities, overlap refusal, bounded sampling, expiry, changed
threshold/source/counts, unchanged-record reuse, publication failures and
existing native multi-GPU significance execution. Four expected reader-depth
warnings remain in the log. A prior persistence lifecycle test used a timing
gate sensitive to cold trace initialization; it now tests lifecycle with a
nonbinding test budget. Controlled-clock overrun tests and production defaults
remain unchanged.

The public execution report is
`results/indexed_output_partition_execution_v1_20260922/report.json`.
`direct_indexed_output_partition_20260922.py` ran 12 fresh processes: control,
first measurement and reuse for each of four output configurations. All scan
the nonzero source interval [129, 1154), with 2,049 samples, 1,025 variants,
512 phenotypes and chunk size 128. This exercises a source tail and prevents
output-relative positions from being mistaken for source indices. Unit tests
also cover phenotype-tile tails.

| Configuration | Execution | Result rows per run | Retained evidence bins |
| --- | --- | ---: | ---: |
| Significant pairs, threshold 0.01 | Four 128-trait tiles on two GPUs | 5,264 | 6 |
| Significant pairs, threshold 1e-300 | Four 128-trait tiles on two GPUs | 0 | 6 |
| Significant pairs, threshold 0.01 | Complete panel on one GPU | 5,264 | 3 |
| JAGWAS | Two variant shards, complete panel per GPU | 1,025 | 6 |

Every measured/reuse output has identical pair/variant identities and
significant-pair df to its control, with statistics within the existing FP32
tolerance. Actual part
contents were checked against event coordinates. Each producer's chunk ranges
are contiguous and cover its entire configured partition exactly once. All
36 empty chunks per empty-output run retain their producer identity without
creating parts. The four reuse processes reference the original output-count
record hashes, publish no replacement count record, and leave all existing
record bytes unchanged.

A thirteenth process injected a writer-callback failure during a rebased,
two-GPU significant run. All torchGWAS worker threads closed and every cache
artifact remained unchanged. This checks the new metadata/rebasing path's
failure cleanup, not just a successful output path.

Input, association output and cache mounts were verified as XFS `/dev/md0` on
`/data`. The benchmark verified all five input hashes and the source identity
before/after execution; the delivered 131 package-source hashes and benchmark
hash match the report. Reports are in the shared project results directory.

The report includes first completed output and output-inclusive API times,
excluding process imports. These single runs on the shared server vary too
much to infer speedup or a general overhead bound. Initial context binding
took 0.267–1.452 s in the observed runs; output-count cache lookup took
1.58–6.21 ms. Binding is still material pre-scan work and needs component-level
profiling before choosing an optimization. These observations do not establish
that the automatic planning budget is met or that forecast accuracy is qualified.
