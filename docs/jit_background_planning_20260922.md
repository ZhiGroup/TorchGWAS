# Released-frontier productive planning

The public first-chunk controller now accepts `initial_chunks.background_planning=true`.
It remains opt-in. After a completed writer event and the existing horizon,
memory and refresh gates, it starts one bounded proposal in a planner thread.
The output callback returns after thread launch; future source reservations and
writer events can proceed during calculator work. The synchronous path remains
the default. A worker is joined before final publication so its cost, cache
state and decision are included in the run audit; a very short scan may still
pay that tail time.

`ProductiveTuningRun.planning_step(release_issue_frontier=True)` captures the
exact reserved-prefix revision and written-output count after lazy cache setup.
It releases the issue lock during analytical work and compares those values,
the current size and run state before applying a proposal. Any change makes
the calculation stale and **cannot** switch the running scan. Its wall and CPU
cost remain charged. The public controller retains that candidate size for a
later completed output event, ahead of untried sizes. A stale result with no
later output cannot be retried. Budget skips and actual planning errors stop
optional tuning. Only one planner worker may run at a time; the final join and
worker errors are audited. Current chunk, phenotype ownership and GPU
assignment remain fixed except for an admitted future chunk-size change.

This addresses a writer-callback stall, not source-scale extrapolation. The
same three short horizons, declared scenario costs and maximum extrapolation
still gate switching. The earlier 8.1M-variant example would fail the
conservative extrapolation cap early, before this worker starts. A successful
background calculation is not a qualified prediction or proof of payback under
contended production load. Scheduling the calculator in a separate thread can
also compete for host resources; planning time remains in the payback test.

Read-only H100 source diagnostics clarify the scale problem. Current
development `PgenHeaderWork` aggregated 1,048,576 frozen variants in 0.658
wall seconds after a 0.010-second header parse and 0.017-second signature
enumeration. It found 20,388 `(record form, byte length)` signatures and
4,228,442,839 payload bytes. Its possible decoder-CPU interval could not be
fully priced: the frozen immutable profile has no `uleb4` or `uleb5` primitive
rates. The partial identified interval was 6.111–9.301 seconds and must not be
called a whole-decoder interval. The current development exact payload census
of the same source took 21.058 wall / 18.795 process CPU seconds and found
797,888,026 one-byte, 29,553,671 two-byte, 6,714 three-byte, and zero four-
or five-byte varints. Those zeros cannot be inferred from the header alone
because the supported decoder permits overlong integers. The exact census is
too expensive for cold upfront planning and would reread the payload. These
are structural diagnostics, not a GWAS timing fit or new service capacities.
The frozen source and plan/profile artifacts were checked unchanged. The
diagnostic project's pulled records are
`results/full_header_envelope_v2_20260922/report.json` and
`results/exact_source_envelope_20260922/report.json`.

The next substantial JIT step is a source-complete continuation with priced
possible decoder units or source-bound exact counts acquired during useful
chunks, plus output-regime and shared-capacity checks. The current frozen
H100 executor prediction still has a large absolute error, so a profitable
production switch remains unqualified.

On lab-a100, the final focused controller/run suite passed 119 tests (job
`20260922-211823-1210235`). Another 107 forecast, binding, price and history
tests passed (job `20260922-212156-1211189`). The public two-GPU JAGWAS audit
(job `20260922-212005-1210462`) wrote 4,097 joint results in each of its
control, deferred and same-process reuse runs; maximum absolute output
difference was zero. Both optional background proposals finished without
errors but were stale: issue revisions advanced 11→33 and 7→33 before their
commit checks. The run kept chunk size 128, and recorded 0.91 and 0.59 seconds
of planning wall time. The fixture has mostly synthetic service prices and
varying warm/load states; its 4.11/1.89/1.04-second API times are not a matched
speed comparison. The executed benchmark SHA and all 144 package source hashes
matched the pulled report at
`results/background_jit_execution_20260922/audit/report.json`.
