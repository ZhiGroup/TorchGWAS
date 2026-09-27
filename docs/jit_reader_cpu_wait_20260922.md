# Reader CPU and scheduler evidence from productive JIT chunks

The bounded first-chunk observer now samples three additional quantities on the
actual native reader thread: `read_cpu_seconds` from `thread_time()`, the
per-thread Linux scheduler runnable-wait counter when interpretable, and the
wall time spent in the probe outside the existing read interval. Unreserved
chunks keep the original `read_into` call. The immutable stage snapshot declares
`thread_cpu_schedstat_v1` so a later consumer can distinguish this protocol.

The scheduler counter is deliberately nullable. On the A100 host the global
`sched_schedstats` switch read 0 and small sampled counters did not advance;
that cannot prove zero runnable wait. On H100, a per-thread counter was nonzero
even while the switch read 0, so the switch alone cannot determine whether a
sample is usable. The observer retains a changing counter and reports `null`
when an unchanged counter is ambiguous. These are loaded reader-thread
diagnostics, not independent CPU capacity or decode price measurements. The
probe wall span can include preemption and is not an exact instrumentation CPU
cost.

The existing cache-validation controller now compares matched reader-thread
CPU spans alongside read wall and CUDA spans. A CPU-only drift can trigger its
bounded refresh. A prior window without this field and a fresh window with it
are incomparable and refresh instead of silently treating the old schema as
current. The comparison still requires identical device, source range, chunk
capacity, payload shape, source/config dependencies and observation age. Reads
do not renew an immutable baseline. This remains a drift heuristic; it does not
change a chunk size or overwrite an independently measured resource rate.

The source-matched public two-GPU JAGWAS audit is
`results/jit_reader_cpu_wait_20260922/public_jagwas_v3/report.json` (pulled
locally). Control, deferred and reuse outputs were equal with maximum absolute
difference zero. Both deferred passes delivered two sampled chunks per GPU and
left no pending reservation. Reader-thread CPU across those samples was
0.709–1.857 ms; scheduler wait was `null` for all eight samples under the
counter-availability rule. The probe's bracketing wall time summed to 5.09 ms
in the deferred pass and 2.70 ms in reuse, below the declared 10 ms measurement
allowance for this particular window. That comparison excludes CUDA event and
observer overhead, so it does not certify the allowance for other workloads.
The shared-server API times (1.682, 2.393 and 1.417 seconds) do not establish a
throughput gain or loss. The audit used the existing mostly synthetic component
prices and made no profitable chunk-size switch.

The final A100 regression passed 179 focused tests spanning adaptive chunks,
cached component validation, the public initial-chunk controller, and run
calibration. The JAGWAS audit ran on the final reader-probe source before the
cache-comparison extension; the latter does not enter that public execution
path.

This closes one evidence gap for just-in-time validation, not calculator
qualification. The original H100 two-GPU full executor remains unmatched while
one GPU is occupied by another process. A passive matched full scan and
independent resource measurements are still needed before turning loaded
scheduler or reader spans into a bounded CPU-availability scenario or using
them to authorize a production tuning decision.
