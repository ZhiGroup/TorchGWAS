# Initial-chunk calibration and reusable evidence, 2026-09-22

The requested behavior is to spread expensive tuning measurements over the
first chunks of one GWAS job, keep their association results, and reuse valid
parameters in later jobs. Tuning must have a bounded early window rather than
continue changing configurations throughout the scan. JAGWAS retains the full
phenotype panel and joint factor on every active GPU; only variants may be
partitioned. Significant-pairs mode can still partition phenotypes.

## Startup and decision contract

There must be no fixed upfront tuning phase. Correctness checks and conservative
memory admission precede the first useful association; a complete candidate
search does not. Start from valid cached evidence or a conservative feasible
configuration. Measurement and planning then proceed incrementally over a small
number of promising alternatives during the early productive chunks. Planning
CPU, elapsed delay, instrumentation, switching and exploration all consume
resources and must be budgeted. Existing blocking public candidate selection
is a separate opt-in path, not the implementation of this JIT contract.

The intended decision model is finite-horizon Bayesian adaptive control with
switching costs, initially approximated with a small parameter/state filter and
a constrained, cost-aware decision rule. Keep the existing analytical resource
calculator as the predictive model. Reusable component parameters and transient
conditions have different roles in y_t = h(a_t, x_t, theta, z_t) + noise:
a_t is the configuration; x_t includes exact source work, shapes, output volume
and pipeline state; theta holds reusable parameters; z_t describes transient
conditions. This is a proposed model, not an implemented or validated posterior.
The current median-ratio drift check is only a heuristic.

Updates require an identifiable observation model. Loaded CUDA spans contain
scheduling gaps, and consumer acceptance is not durable output completion.
Neither can be relabeled as an independent hardware capacity. A posterior
update or refreshed measurement must be a new immutable record with its input
record identities, protocol, dependencies and original observation times.
Posterior publication must not rejuvenate expired contributing measurements.

Candidate decisions minimize expected remaining completion time, including
planning, measurement, switching and extra time spent exploring. Stop when
expected remaining savings cannot repay those costs or the early tuning budget
expires. Short jobs may warrant no exploration. Compare equivalent completed
GWAS work and the requested output; raw chunk latency and survivor counts are
not throughput objectives. Output mode and significance threshold remain fixed
scientific requirements. Only changes admitted at safe execution boundaries
are actions; JAGWAS always retains its full phenotype panel and factor per GPU.

Validation must measure time to first useful output and incremental end-to-end
overhead against matched runs with tuning disabled, separately from productive
sampling duration, planner CPU and throughput improvement. Earlier callback CPU
figures do not establish total JIT overhead or a production throughput benefit.

## Implemented execution support

`adaptive_chunks.ChunkSizeControl` chooses from explicit chunk sizes. The
native pinned reader advances one contiguous source cursor and consults the
control only after acquiring a free slot. Already-issued reads keep their
ranges. Input, device and result rings retain the fixed allocated capacity;
DMA leases and borrowed-result lifetimes are unchanged. The caller must admit
that ring capacity and the workspace/transient allocation envelope for every
permitted chunk shape and tail, including allocator reserves. Admission of
each separately allocated candidate is not automatically admission of their
shared allocation envelope.

The shared allocation model and an aligned-size transition control are now
implemented in adaptive_candidate.py and AlignedChunkSizeControl. They retain
the largest ring, account for mixed old/new storage classes and shifted PGEN
LD restart buffers, and enumerate required full/tail geometry. See
adaptive_candidate_20260922.md for exact scope and numerical verification.

Explicit schedules can now be regrouped from a fine source census without
rereading the input. The calculator charges actual LD restarts, transfers and
writer work at each actual chunk boundary. See scheduled_census_20260922.md.
Analytical checkpoints now preserve completed and in-flight model work while
scoring unissued future ranges; see candidate_continuation_20260922.md. These
checkpoints are model state, not measurements of live executor queues.
Bounded computational reuse and a per-step cost gate are described in
planning_session_20260922.md. Connecting them to the public productive-run
lifecycle and a validated decision policy remains open.

The internal productive-run bridge now records every reserved source interval,
including unsampled prefetch, and opens its planning window only after an
indexed writer completes useful work. It can hold the issue frontier during
one bounded analytical proposal and apply a profitable timely size change to
future reads. See productive_run_20260922.md for the two-GPU execution audit,
writer-fsync boundaries and the remaining automatic-policy/public integration.
The source-calculator callback in jit_proposal_20260922.md now produces one
future-size proposal from that exact prefix, with conditional remaining-work
bounds and no upfront candidate-grid search.

The internal `linear_scan_streaming_chunks` hooks `_chunk_size_selector` and
`_chunk_observer` require native CPU fill, CUDA FP32 dosage statistics and an
explicit ring capacity. Unsupported paths fail before phenotype preprocessing.
The multi-GPU driver shares preprocessing and gives each worker the same
thread-safe control, while observations identify their device and exact source
range. Joint factors remain distinct full-width factors on their GPUs.

`InitialChunkMeasurements` reserves a bounded number of samples per GPU at read
issue. Warmup, sampling stride and an early wall-clock window spread those
samples over real production chunks. Both completed and in-flight samples
count against the limit. Once the limit is reached, no additional measurement
CUDA events are recorded. The window deadline stops new reservations; it does
not cancel previously issued work.

Measurements contain reader wall time, submission/delivery/consumer timestamps,
output bytes, and optional H2D/conversion/statistics-and-reduction/result-stream
CUDA event spans. They introduce no profiler and no extra CUDA stream drain.
These intervals overlap and include scheduling gaps. They are not independent
hardware capacity measurements and must not be added as isolated stage costs.
Device significant output reports one observation after all selection blocks
for a source chunk; its result-stream transfer span is unavailable rather than
zero. Consumer resumption can mean asynchronous queue acceptance, so it is not
a durable-write timestamp. A generator closed at a yield does not report that
chunk as fully consumed.

## Cache validity

`calibration_cache.CalibrationParameterCache` publishes content-addressed JSON
records atomically without overwriting earlier records. A refresh appends a new
record and preserves the evidence used by earlier decisions.

- Structural records (`device_properties`, `source_work`, `kernel_geometry`,
  `workspace`) have no age expiry, but require identical explicit dependencies.
- Empirical capacities and `stage_observations` require an explicit maximum age.
  A later caller may shorten this age but cannot extend the producer's expiry.
- `available_memory` and `contention` may be retained for audit, but lookup always
  requires a new live observation.

Bindings must include the relevant source/library hashes, physical device,
shape, math and thread settings, storage path/mount, and measurement protocol.
The cache does not infer missing dependencies or validate the scientific
meaning of a value. Corrupt, modified, future-dated, expired or differently bound
records are misses. A valid cache hit is evidence, not candidate admission or
proof that shared-resource performance has remained stable. Early observations
are still needed to detect changes within an empirical record's lifetime.

`InitialChunkMeasurements.publish` saves completed observations with explicit
provenance and expiry. Pending coverage remains labeled pending. The two-job
audit stores the shorter cache-hit validation window separately so it does not
replace the fuller original measurement window.

## Automatic checks and bounded refresh

`initial_calibration.InitialCalibrationController` now manages a separate
empirical baseline for each GPU. With a valid cached baseline it initially
reserves two matching chunks per device. If their component intervals are
consistent, measurement stops for that device. Drift or incomparable geometry
expands the same window to at most eight samples per device. A missing or
expired baseline starts with the eight-sample budget. These counts, a ten-second
early reservation window, the sampling stride, and the comparison tolerances
are configurable. Extending a check never resets the count or window deadline.
Already-issued samples finish normally; all useful GWAS output is retained.

The comparison matches exact device/source-range/capacity tuples, checks output
bytes/block counts, and compares median reader and CUDA event spans. Its default
twofold ratio and two-microsecond absolute tolerance are declared drift
heuristics, not confidence intervals or validated performance thresholds. Event
spans can vary with host scheduling and resource contention. A comparison with
no result-transfer metric leaves it unavailable. Consumer suspension is not
used as a hardware or writer service rate.

The controller also records cumulative callback thread CPU time. The default
50 ms budget stops new reservations once exhausted. This is a cooperative limit:
a callback already running and in-flight completions can exceed it. It excludes
CUDA event overhead, cache lookup/publication, and candidate planning, so it is
not a bound or measurement of total tuning overhead. No profiler or extra CUDA
stream drain is introduced.

After the complete scan and consumer succeed, `finish(successful=True)` appends
new per-device baselines only for completed collection/refresh windows. A
consistent short check is published under a separate validation name and does
not renew the baseline's age. Incomplete or failed windows do not publish a new
baseline. Optional cache-write I/O failures are reported without discarding the
GWAS results. Historical records remain immutable; there is no invalidation
tombstone or automatic cache pruning.

Empirical age starts at observation time, not publication time. Finishing a
long job cannot make its early samples fresh again. Cache lookup prefers the
newest observation among valid records, even if an older measurement was
published later by a longer job. This also preserves compatibility with older
records that used publication time as observation time. The controller itself
requires explicit observation time for reusable component windows.

## Validation

Remote A100 validation, all artifacts pulled:

- Job `20260921-204634-730786`: 76 initial adaptive/native/JAGWAS checks passed.
- Job `20260921-205710-734242`: 138 execution, sampling, cache and significant
  output checks passed, including variable chunks across native PGEN LD replay,
  ring reuse, early close and error cleanup, full JAGWAS factors on two GPUs,
  multiple selected-output blocks, bounded reservations, expiry, context changes
  and immutable publication.
- Job `20260921-210057-735010`: 65 calculator/source-work regressions passed.
- Job `20260921-210027-734915`: two fresh processes ran a synthetic native-PGEN
  fixture with N=2049, M=32768, K=129, three covariates, chunk 128 and cuda:0+cuda:1.
  Input and metadata are on local `/data` (`/dev/md0`, XFS). Each process retained
  all 32768 joint statistics from 256 chunks. The first saved 16 observations;
  the second manually used the cache-hit policy and took four fresh
  observations. Output SHA256 was identical across processes:
  `9c41fded8ca6be507ce35fb8c68e5adc2050809cfcc7c22d29c098fd2a138c18`.
  Forty-one joint statistics matched an independent FP64 OLS/correlation
  calculation, maximum absolute chi-square error `4.139426076221753e-05`.

The fixture establishes correctness and cache lifecycle, not throughput,
ranking, overhead, or public autotune readiness. A concurrent calculator test
job used the same host; observed component spans remain loaded-context samples.
The read-only NumPy mmap warning came from existing preprocessing, which reads
the fixture without mutating it; the numerical output and fixture identity
checks passed.

Artifacts are under `results/adaptive_chunks_v1_20260922`,
`results/initial_chunk_measurements_v2_20260922`,
`results/adaptive_chunk_model_regressions_20260922` and
`results/initial_chunk_cache_audit_v1_20260922`. The audit harness is
`benchmarks/direct_initial_chunk_measurements_20260922.py`.

### Automatic refresh audit

Job `20260921-212304-742161` passed 75 cache/controller/native-execution tests
(31.09 s), then completed four fresh processes on the same local-PGEN fixture.
Artifacts: `results/initial_refresh_controller_v2_20260922` and
`results/initial_refresh_audit_v2_20260922`. Harness:
`benchmarks/direct_initial_refresh_20260922.py`.

| Process | Prior empirical cache | Final state on both GPUs | Samples, total | Callback CPU time |
|---|---|---|---:|---:|
| First | Miss | Collected | 16 | 12.57 ms |
| Second | Hit | Refreshed | 16 | 14.05 ms |
| Controlled reader delay | Hit | Refreshed | 16 | 7.43 ms |
| Delay removed | Hit | Refreshed | 16 | 9.14 ms |

The second normal run also differed sufficiently to trigger expansion, so this
audit does not demonstrate a reduction in the number of fresh samples. Unit
controls cover the consistent-cache two-sample path and selective per-GPU
refresh. The delay process injected 30 ms inside native fill as a test-only
change in wall service, producing reader-span ratios of 22.46 and 30.08 against
the cached validation ranges. Removing it produced ratios of 0.065 and 0.146.
This verifies refresh decisions, not storage or CPU capacity. The audit cache
is isolated and must not be used as production capacity calibration.

All four processes retained all 32768 statistics from 256 source chunks,
produced the same output SHA256 reported above, and matched the independent
FP64 reference with the same maximum absolute error. Every process completed
its bounded windows without exhausting the callback CPU budget. These timings
cover callback CPU only; total measurement overhead was not benchmarked.

Job `20260921-212836-743882` then passed all 31 cache/controller tests (1.10 s),
including two new out-of-order publication cases: explicit observation times
and a mixed legacy cache. Artifacts:
`results/initial_refresh_freshness_v3_20260922`. This final cache-only change
prefers observation time when selecting among still-valid records. The v2
four-process artifacts retain their original source hashes; no numerical or
scan code changed afterward.

## Remaining work toward the active goal

The public API now exposes bounded productive observation collection and cache
validation through `initial_calibration`; see
[public initial calibration](public_initial_calibration_20260922.md). It uses an
explicit configuration and does not choose or switch it. The public API also
supports opt-in significant-host and JAGWAS candidate
selection using an immutable component-price record and live memory admission;
see `significant_public_autotune_20260922.md` and
`jagwas_public_autotune_20260922.md` for their separate validation. The
initial-chunk controller still needs to connect its evidence to public candidate
selection and safe transitions. Loaded component spans are not independent
hardware capacities and must not be used as fitted whole-GWAS timing tables.
Choosing or changing phenotype tiles and GPU assignments requires appropriate
safe scheduling boundaries and full coverage accounting; the chunk mechanism
alone does not implement those decisions. The current controller does not
choose or switch an execution configuration.

The H100 ranking checkout remains frozen. A100 source hashes changed, so earlier
whole-source profile bindings must not be presented as current-source evidence.
The broad calculator, autotune, tiling and multi-GPU goal remains active.
