# Productive stage observations for JIT tuning

`initial_chunks.stage_observations` now attaches the existing bounded native
chunk observer to the first admitted scan partition on each GPU. It begins at
actual source-read issue, after startup admission, and stops reserving samples
at its declared count or wall-window limit. Its optional CUDA events time H2D,
conversion, statistics/reduction and result transfer; the read interval covers
input and decode together. Later phenotype tiles on the same GPU do not get a
second measurement window. The source scan and its numerical output are the
work being measured; no separate genotype probe or upfront tuning grid runs.

The later [reader CPU/wait extension](jit_reader_cpu_wait_20260922.md) adds
bounded thread-CPU and nullable scheduler diagnostics to these same productive
chunks and compares thread CPU across matched cached windows.

The option declares `max_chunks_per_device` (at most 32), `warmup_chunks`,
`stride`, `max_window_seconds`, `cuda_events`, and a positive
`measurement_reserve_seconds`. The last value is added to every productive
planning cost gate. It is a caller-declared allowance for instrumentation,
not an independently measured upper bound; CUDA event/bookkeeping overhead
still needs qualification before using these observations to authorize a
production performance claim. Loaded read and CUDA intervals can reveal model
drift or scheduling variation, but they cannot replace independent CPU, GPU,
transfer or storage capacity prices. This step does not publish a new price
record or renew the observation age of an existing one.

The final A100 regression passed 150 focused tests across productive control,
JAGWAS and significant-pairs paths. Four profitable-switch integration cases
also exercised the observer in dense trait tiles, dense variant shards,
significant pairs and JAGWAS, with result equality and per-GPU CUDA samples.

The public two-GPU JAGWAS audit is
`results/jit_stage_feedback_20260922/public_jagwas_v4/report.json`. Control,
deferred and reuse outputs were equal (maximum absolute difference zero).
Each deferred pass wrote 33 indexed chunks, retained 32 bounded output
observations, and collected two native stage observations per GPU with no
pending reservation. The cost gate included the declared 0.01-second
measurement reserve; actual total tuning costs were 2.96 and 2.03 seconds.
Neither pass switched chunk size. API times of 2.55, 7.64 and 4.03 seconds
were collected under variable shared-server conditions and mostly synthetic
component prices, so they establish no throughput gain or loss from the
observer.

Next, independently priced bottleneck services need source-compatible
incremental refresh and immutable cache bindings. Productive stage intervals
can then test whether those prices transfer under concurrency. A switch must
repay measured planning cost and a qualified instrumentation allowance;
phenotype tile and device moves require separate admission and a bounded
layout transition. JAGWAS remains full-panel on every GPU.
