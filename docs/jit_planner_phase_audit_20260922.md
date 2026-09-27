# Productive calculator callback phase audit

The public two-GPU JAGWAS controller now records bounded phase timings inside
its productive analytical callback and up to eight successful live-validation
phase rows. They are loaded elapsed spans, **not** independent CPU, GPU,
transfer or storage capacities. The timing calls are inside the charged
planning step. The output and source semantics do not depend on them.

In `results/jit_planner_phases_20260922/public_jagwas_v1/report.json`, the
cold proposal took 2.297 s: 0.619 s in two live evidence checks, 0.860 s in
source windows and 0.757 s in schedule comparisons. Its repeat proposal took
0.780 s, including 0.326 s validation, 0.087 s source windows and 0.333 s
comparisons. Control and tuned outputs were equal, with no switch. Server load
varied, so these are diagnostic partitions of each callback, not a causal
warm-cache speedup.

The source-window construction repeated exact PGEN `(start, stop)` bounds
across three horizons and two chunk sizes. A per-header, at-most-64-entry
default cache now reuses immutable structural bounds, checks file identity on
each lookup, and gives callers independent copies. The public JAGWAS audit at
`results/jit_planner_bounds_cache_20260922/public_jagwas_v1/report.json`
recorded 18 hits and 18 misses in each proposal, equal control/tuned output,
and no switch. Its cold and repeat source-window phases were 0.438 and 0.045
seconds, but cannot be compared causally to a different shared-server run.

The matched header-only replay at
`results/jit_header_bounds_cache_20260922/paired_v1/report.json` alternated
cache off/on for eight trials on the same public input and exactly the first
proposal's two-shard window sequence. Median window time was 54.8 ms without
and 46.9 ms with the bounded cache; every returned work ledger matched. The
improvement is real in this isolated replay but small relative to full planner
cost. It excludes context checks, graph comparisons and GWAS execution.

The later public audit at
`results/jit_validation_phases_20260922/public_jagwas_v1/report.json` identifies
the other recurring costs. In its repeat proposal, two validation checks took
0.445 s; within them, execution-context capture took about 0.168 s,
profile/price validation 0.150 s, and digest rechecks 0.092 s. The schedule
comparison took 0.356 s. The cold pass under different load took 1.558 s
overall, with 0.534 s validation, 0.273 s source windows and 0.676 s
comparison. Results again matched exactly and no switch occurred.

The next useful change needs to reduce the context and graph costs while
retaining a current check before expensive modeling and a final check before
switching. Removing the first full check failed the stale-evidence guard's
fast-failure tests and was reverted. A second avenue is speculative
background planning with exact frontier revalidation, but simply moving the
current callback to another thread while holding the issue lock would still
stall source submission. A late proposal must never switch using an obsolete
issue frontier. Neither avenue is yet implemented or qualified.
