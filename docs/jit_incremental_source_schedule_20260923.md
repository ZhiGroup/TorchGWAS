# Incremental PGEN source accounting for JIT planning

`IncrementalPgenSchedule` inspects one bounded, chunk-aligned PGEN header
segment per `advance()` call. It reads index metadata only. After all
segments, `finish()` conserves the same complete source-unit intervals,
indexed read bytes, LD restart entries, variable-integer constraints and
chunk count as one direct `PgenHeaderWork.schedule_bounds()` call. The staged
primary record work can then `rebase(cursor, chunk_size, stop=...)` to an
actual later unissued cursor and variant-shard endpoint. Only a bounded
partial edge and the new chunk-start LD replays are inspected; the primary
work already accumulated across earlier chunks is reused. This is important
because source issuance continues while background planning runs, and a
chunk-size or GPU partition proposal must model the current remaining range.

`productive_source_floor(..., staged_source=...)` now accepts a completed
staged ledger. It checks the PGEN identity and exact frontier coverage, then
reuses the ledger for each unique candidate variant span. Trait tiles may
reread the same range, and variant shards may have different endpoints. The
report separates prior staged CPU/wall cost from the current rebase call.
Passing the completed public `ProductiveSourceStage` now reports its full
measured worker-step cost, including header setup, optional live-output
observation and finalization. Passing the lower-level
`IncrementalPgenSchedule` alone still reports only its inner calculation;
a JIT controller must charge any omitted outer work once per job before a
decision.
The public initial-chunk controller does not yet use this ledger for a
switch. The completed-stage/full-worker-cost source rebase and related source,
stage, incremental-schedule and public-controller regressions passed 136 A100
tests in `20260923-062129-1376227`. No source floor or short-window estimate
is a complete pipeline-time ceiling.

The public controller now has an **opt-in evidence-only** `source_staging`
path. Each useful output event grants one background `advance()` step; the
writer callback only queues the work. PGEN scheduling imports and header
construction occur inside the first background step. The path bounds records
per step, step count, cumulative thread CPU, elapsed window and retained
Python ledger. It reserves declared extra host memory at admission. If the
job ends while a step is in flight, finalization joins that single step and
discards queued
ones. The audit records every step's cursor, CPU, wall time and retained
ledger size. The source identity is checked before and during calculation.
The reserve is an admission envelope, not a hard operating-system RSS limit.
Because this work is not yet charged in a complete payback comparison,
enabling `source_staging` disables short-window layout switches. The usual
chunk tuning behavior is unchanged when staging is absent. A later controller
can use the completed ledger only after charging its observed cost.

The read-only A100 probe used the existing 8,086,101-variant,
22,250-sample hardcall PGEN at
`/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen` (20,838,552,600
bytes). At B=128, the complete schedule has 63,173 chunks and 15,475 LD
restarts. Staged/direct/staged matched every work field, excluding only the
descriptive scope string. In the 1,048,576-record-step control, the first
eight steps took 9.04 s total, with a 1.59 s maximum step; the later warm
staged pass took 0.84 s, and a direct pass between them took 1.37 s. Header
construction was separately 0.92 s. The first staged pass reached 159 MiB
peak RSS, versus 285 MiB after the direct pass, but peak RSS is cumulative
within that process. These are metadata-calculation costs under different
cache/load states, not genotype scanning throughput or proof that staging is
faster. They show that a cold whole-source calculation cannot be put in the
first-output callback.

A later control rebased an already staged 8M ledger from cursor 8,320 to a
4,040,000-variant shard at B=256 and matched a direct schedule exactly.
Its first staged pass took 11.80 s and the warm pass 0.82 s; rebase took
1.33 s versus 0.94 s direct after eliminating a duplicated replay-signature
pass. The direct calculation remains faster when the header and primitive
cache are warm. Rebase preserves primary work distributed across earlier job
progress and avoids a fresh whole-source primary walk at a changed frontier;
whether that reduces total planning wall time remains unproven.

Smaller steps are not automatically cheaper. A separate 131,072-record
control needed 62 steps; its first and second staged passes took 38.73 and
44.79 s, with maximum individual steps of 1.27 and 1.46 s. The intervening
direct pass took 2.75 s. This reflects the current bounded signature/bounds
cache and per-segment work under that run's load; it rules out simply
reducing the step size to make the total tuning cost negligible. A controller
must choose a step budget using measured planner costs and stop when the
expected benefit cannot repay them.

Direct-first controls separate some cold-load effects. With the default
1,024-signature/64-bound cache, a direct pass took 4.00 s, the later first
1M-step staged pass took 10.87 s, and the warm direct pass took 0.24 s
(`20260923-035735-1346651`). Raising only the probe cache to 32,768
signatures/256 bounds gave 2.39 s direct-first, 3.38 s across eight staged
steps (maximum 0.73 s), and 0.15 s direct after; peak process RSS after the
first stage was 266 MiB (`20260923-035836-1346780`). At 131,072 records
per step with that large cache, 62 staged steps still took 10.60 s (maximum
0.61 s), with 266 MiB peak after the first stage and 289 MiB by the warm
stage (`20260923-035856-1346847`). The cache held 16,181 unique signatures
in these controls. These process-order comparisons do not isolate the
staging CPU from other runtime work, and the production evidence-only path
currently retains the default smaller header cache. The optional public `source_staging` settings now accept
`max_cached_signatures` (at most 32,768) and `max_cached_bounds` (at most 256).
The staging header uses those values only after the first useful output.
Admission charges the declared extra host reserve; above the old 1,024/64
cache, each additional entry requires a conservative 16 KiB reservation on
top of `max_retained_bytes`. For example, a 32,768/256 cache with a 1 MiB
retained-ledger cap needs a 500 MiB `extra_host_reserve_bytes` setting. This
is an explicit planning envelope, not a hard RSS or Python allocator bound.
Each staged step and the final audit record cache hits, misses and capacity.
The source and public-hook regression passed 127 A100 tests in
`20260923-055510-1373733`; the exact 500 MiB admission edge passed a focused
follow-up in `20260923-055742-1373916`.

Four new A100 runs used the same local-XFS 8,086,101-variant PGEN and B=128,
1,048,576 records per step. All four matched 63,173 complete chunks and 15,475
LD replays exactly. Run order was large/default/default/large cache, each in a
fresh process. The 32,768/256 cache had 16,181 signature misses and 110,383
hits in both runs; the default 1,024/64 cache had 124,584 misses and 1,980
hits. First staged-pass thread CPU seconds were 4.151, 8.611, 5.600 and 3.258
in that order, giving two-run medians of 3.704 s (large) and 7.106 s
(default). Peak process RSS at the end of the first stage was 200/206 MiB
for large and 157/162 MiB for default. That is cumulative process RSS, not
isolated Python cache allocation. Header construction ranged from 0.218 to
1.119 CPU seconds; it is excluded from the staged-pass CPU values. These
ordered metadata-only measurements are confounded by shared load and cache
state. They support the signature-thrashing diagnosis but do not measure a
whole GWAS speedup or qualify a JIT switch. Reports are in
`results/jit_source_cache_20260923/{large_first,default_second,default_third,large_fourth}.json`;
remote jobs were `20260923-055807-1374016`, `20260923-055833-1374082`,
`20260923-055905-1374168`, and `20260923-055933-1374278`.

`benchmarks/incremental_pgen_schedule_probe.py` records input identity,
source-code hashes, per-step wall/CPU time, peak process RSS and exact-work
checks. A100 job `20260923-034150-1343400` produced the original 1M-step
report; job `20260923-035403-1346245` produced the rebased-shard control.
Job `20260923-035442-1346352` compared the 128k-step setting.
The complete PGEN/source/checkpoint regression passed 235 A100 tests in
`20260923-035322-1346070`. Small mixed-form tests cover LD boundaries,
source tails, later cursors, changed chunk grids and shard endpoints; they
also prove that the productive source floor avoids a direct whole-header
rescan when a complete staged ledger is supplied.

The first public-hook run `20260923-040538-1349658` passed 303 tests and
exposed 31 older `initial_chunk_autotune` tests whose synthetic writer events
did not carry the now-required producer/output identity, or whose snapshot
omitted the written-event count. Those fixtures were corrected without
weakening the output-boundary checks. The 38 affected cases passed in
`20260923-041024-1350317`. The combined A100 regression passed 335 tests
in `20260923-041358-1350758`, including real PGEN accounting through the
public writer callback. The final import-deferral and worker-failure
regression passed 336 A100 tests in `20260923-041715-1351211`.

Next, the controller needs to charge this nonblocking source stage and use
its rebased result in a **complete** finite continuation. The current
calculator still lacks synchronized producer/queue state, mode-specific final
service, a calibrated completion ceiling and a qualified layout transition.
The frozen H100 absolute prediction gap remains open.
