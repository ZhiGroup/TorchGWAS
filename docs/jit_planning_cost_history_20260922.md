# Immutable planning-cost history for productive JIT tuning

`initial_chunks.planning_cost_history` is an optional, bounded record of CPU and
wall time actually spent in a productive analytical proposal. The controller
reads it only after a completed output event and the cheap source-frontier gate.
It can only **raise** the next job's declared `expected_cpu_seconds` and
`expected_wall_seconds` before the cooperative planning budget admits a step.
It never supplies an association rate, corrects a resource price, or approves
a layout switch. The current job also charges history construction and lookup
wall time as a reserve and declares a separate publication allowance. Actual
publication time is audited after a successful run.

The immutable empirical record binds source-file identity, complete calculator
profile hash, workload and output geometry, reduction, partitions, admitted
chunks, forecast horizons/options, host and occupancy scenarios, and enabled
cache/instrumentation modes. A configured maximum age limits reuse. A lookup
does not renew the source observation timestamp. A failed or incompatible
record becomes a miss and the current declared forecast remains in force. At
most 32 productive steps can be recorded; the cache is advisory, so a failure
to publish it cannot invalidate a completed scientific scan.

The A100 focused suite passed 153 tests, including staleness, changed-source
binding, immutable lookup, and bounded-step checks. The public two-GPU JAGWAS
audit on the final source, with both structural reuse and stage observations,
is `results/jit_planning_cost_20260922/public_jagwas_final_v2/report.json`.
Control, deferred and reuse outputs were identical. The first productive
proposal took 2.878 s wall; the later job reused that original observation and
raised its expected step wall from 0.020 to 2.878 s. Its actual proposal took
1.121 s, including 0.431 s of current validation. Both tuned jobs collected
four stage observations and made no switch. The control, deferred and reuse API
times were 4.59, 6.57 and 3.55 s under variable shared-server load. These
times are not a measured speedup claim.

An earlier public two-GPU JAGWAS
audit with stage observations is
`results/jit_planning_cost_20260922/public_jagwas_v1/report.json`. Control,
deferred, and reuse outputs were equal (maximum absolute difference zero).
The first productive proposal took 0.652 s wall and 0.312 s of the two live
validation checks. Its next-job history hit raised the expected wall from
0.020 to 0.652 s, preserving the first observation time; that second proposal
actually took 2.535 s wall, including 1.441 s of validation. Both jobs
collected four GPU-stage observations and made no chunk switch. Shared-server
load and mostly synthetic component rates prevent a throughput conclusion.

The separate structural-cache control at
`results/jit_planning_cost_20260922/public_jagwas_structural_v1/report.json`
also preserved exact output. Its cold and repeat proposal steps took 0.636
and 0.395 s, respectively; repeat validation took 0.167 s, and history still
raised the declared next-step wall from 0.020 to 0.636 s. These are two
unmatched runs under variable load, not a causal cache-speed estimate.

This addresses the hidden cost of deciding whether to run the calculator; it
does not remove a long admitted callback from a productive writer thread. The
next controller improvement should either produce a sufficiently cheap upper
bound on repayable gain before that callback or run a bounded proposal without
holding output delivery, then revalidate the source frontier and evidence
before any switch. A cross-job cost history cannot certify current load.
