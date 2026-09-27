# First-chunk downstream CPU and wait observations

The productive chunk observer now records caller-thread CPU time and nullable
Linux scheduler runnable wait while a result is yielded to its downstream
iterator. This is the same bounded, useful-work window that already records
reader CPU/wait and optional CUDA events. The observer samples only reserved
chunks. The existing `consumer_seconds` wall boundary is unchanged: the new
probes sit outside it and therefore do not move its start or end timestamps.
`consumer_probe_wall_seconds` separately records their bracketing wall cost.
On a loaded host, sampled runnable wait can exceed the narrower consumer
wall interval because preemption during a probe is included in the wait
counter. The two values must not be subtracted to infer active service.

The downstream iterator differs by output mode. It may select significant
pairs, submit queue work, or perform synchronous output; it is not a pure
selector or writer-service measurement. Runnable wait is null when the kernel
counter is unavailable or ambiguous, or when a generator resumes on a different
thread. It is never inferred to be zero from an unchanged ambiguous counter.
Abandoned yields receive null CPU, wait and probe-wall fields. Thread CPU
from normally resumed yields is comparable across matching
source ranges, but is not itself an independent CPU-capacity price. The
observation also does not establish indexed-part fsync or dense-file durability.

`InitialCalibrationController` stores this protocol as
`initial_chunk_components.v3:<device>` and requires the protocol marker on a
cached record. Earlier v1/v2 records remain immutable and cannot be reused as
if they contained the current probe boundaries. The first two matching chunks
validate a candidate
cached window; a changed CPU interval triggers the existing bounded refresh.
Observation timestamps are the original first read in the collected window,
even if publication occurs at job end. Structural dependencies and empirical
age checks still govern cache reuse separately.

The GPU integration test passed all six dense/JAGWAS/significant-pairs by
CUDA-events combinations on lab-a100 (job `20260922-222850-1237617`, 6 passed
in 39.06 s). A separate CPU/unit batch passed 22 tests in 30.57 s (job
`20260922-222631-1237342`); a follow-up 17-test calibration batch passed
in 0.61 s (job `20260922-223526-1238966`). These validate probe lifecycle
and cache fields at the previous v2 boundary. The current v3 boundary passed
27 targeted CPU/GPU cases in 26.63 s (job `20260922-224629-1242333`). None
of these tests establish a throughput gain or authorize a production switch.

The historical v2 native PGEN public lifecycle audit (job `20260922-223210-1238331`,
`results/consumer_probe_v2_20260923/report.json`) completed 15 fresh-process
control/first/second runs over dense single-GPU, dense trait tiles, dense
variant shards, JAGWAS and significant-pairs trait tiles. All first/second
outputs matched their control arrays exactly, including selected-pair
identities. Every second run found its prior per-device record; prior record
bytes stayed unchanged. Three injected writer/metadata failures published no
new record and left earlier records unchanged. All 144 recorded package-file
hashes matched its recorded development source. Consumer CPU was present in the
bounded observations, while scheduler wait was sometimes null as intended.
All later windows refreshed on interval drift under varying shared-server
load. The audit establishes output/cache lifecycle correctness, not timing
improvement or an independent component-capacity calibration.

The v3 public rerun (job `20260922-224728-1242504`,
`results/consumer_probe_v3_20260923/report.json`) again completed 15 normal
runs and three injected failures. Across all five modes, control versus first
and later output had maximum array difference zero; all later per-device
windows found their cache records, original bytes were unchanged, failures
published no record, and all 144 source hashes matched the recorded revision.
Every first-window observation had a probe-wall field; nullable wait occurred
as designed. A three-chunk probe bracket reached about 25 ms on one loaded
significant-pairs GPU, exceeding the illustrative 10 ms JIT measurement
reserve. This is elapsed time around the probes, including possible scheduler
delay and overlap, not intrinsic probe CPU service.

The productive JIT cost gate now uses the larger of its declared measurement
reserve and the sum of completed reader/consumer probe-wall brackets when
checking a proposal or refresh step. It records that known cost in the final
stage audit. This avoids undercharging an already observed loaded probe, while
the declared reserve still covers unobserved future samples and CUDA event
bookkeeping. Summing overlapping GPU brackets is deliberately conservative;
neither quantity certifies a worst-case future instrumentation cost. The
current JIT controller and delivery regression batch passed 100 remote tests
in 63.58 s (job `20260922-225849-1245084`), including real GPU output paths
and both declared-reserve/observed-probe cost cases.

The next calculator step is a finite continuation over the exact unissued
source, combining the indexed PGEN schedule with per-chunk GPU, transfer and
output service, shared-capacity constraints, in-flight work and final drain.
Loaded CPU/wait observations should challenge that model under concurrency,
not replace its independently measured service prices. The frozen H100
absolute prediction gap remains unresolved.
