# Passive capacity check on the frozen H100 executor

The nominal frozen calculator selected B=1,024, trait width 8,193 and two
H100s, but predicted 15.862 s against a historical observed median of
55.769 s. A new read-only diagnostic ran three fresh processes with the same
frozen source, plan, input, phenotype, GPU pair, 4.5 GiB per-GPU allocator
cap, numerical settings and 140 reference checks. The external 100 ms monitor
read each child thread's Linux `schedstat` and the child cgroup's CPU
statistics. No frozen package, price record or production executor changed.
The input plus metadata lived on local `/data` ext4; each process passed its
cold-page checks and retained read–scan–read controls. All three produced zero
significant pairs, exact df and maximum sampled absolute t error 1.297e-5.

| Run | Executor s | API s | Process CPU s | Monitored runnable wait, aggregate s | Decode CPU s | Host selector CPU s | Result completion CPU s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 22.135 | 30.699 | 106.986 | 28.531 | 21.897 | 19.605 | 46.717 |
| 2 | 30.210 | 41.730 | 126.711 | 17.112 | 25.169 | 29.862 | 48.594 |
| 3 | 25.442 | 34.561 | 109.214 | 24.332 | 24.514 | 22.854 | 38.848 |

The slowest new run had the *least* aggregate runnable wait, while its selector
CPU increased by 10.257 s and decoder CPU by 3.272 s over the fastest run.
This weakens a scheduler-wait-only explanation. It does not prove that the
selector alone caused the extra eight executor seconds: decoder, GPU,
completion, overlap and shared-resource load all varied. The previously frozen
independent prices put decode/read at 18.579 CPU s and selection at 14.591
CPU s. Every new loaded observation exceeded those isolated prices, so an
early-job JIT check should challenge their transfer under concurrent work.
Loaded intervals are not replacement component prices. The frozen
package uses the earlier `bounded_flat_v1` selector; the current development
checkout uses `bounded_flat_v2`, so these CPU totals cannot be transferred
directly to a new production profile. The 15.862 s nominal
prediction still misses the new 25.442 s median by 9.580 s; the historical
55.769 s median was not reproduced under current shared load.

The external monitor observed 77 child TIDs across each complete process,
including import and API preparation. Its summed schedstat runtime therefore
exceeds the API process-CPU boundary. The reported runnable-wait counter is
aggregated across overlapping threads, may omit short-lived threads and cannot
be added to executor wall time. Linux exposed cgroup CPU usage but no
`nr_throttled`/`throttled_usec` fields for this session, so the report neither
confirms nor rules out quota throttling. The monitor's single `affinity` field
was sampled before the child applied the frozen CPU affinity; the scan's
execution context separately checked that affinity against the frozen plan.
Cgroup pressure and runnable wait show some contention but cannot establish a
stable available-core fraction or a conditional runtime ceiling.

The executed script is
`/home/x/work/torchGWAS-calculator-diagnostics-h100/passive_capacity_monitor.py`.
Remote job `20260923-061200-1375300` completed; pulled reports are under that
project's `results/frozen_passive_capacity_20260923/`. The completion record,
all three scan/monitor reports and the local script match their recorded source,
plan and diagnostic hashes. This evidence directs the next calculator work to
mode-specific loaded service validation in the first chunks and a finite
output-inclusive continuation, rather than a global timing multiplier. No
production JIT switch is qualified by these three runs.
