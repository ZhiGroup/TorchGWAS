# Dense live-resource admission

The dense calculator now retains a compact admission row for every feasible,
priced candidate in its bounded search. The rows contain each candidate's
memory request, devices, API settings and worst declared scenario cost; they
omit the scenario execution graphs. A cached structural ranking remains
unchanged. Before execution, fresh host and per-device available memory filters
these rows. Among the survivors, the planner applies the original
`max_slowdown_fraction` and fewest-device preference to the new live best.

Previously only the planned dense winner was checked. If another process used
memory after planning, that choice could fail even though an already priced
layout with a smaller envelope fit. Significant-pair and JAGWAS plans already
retained their ranked feasible candidates, so their admission path is unchanged.

Old dense cache entries without the admission list are recalculated once; they
are not interpreted as proof that no alternative fits. This is structural
reuse, not a new empirical observation. Calibration prices keep their original
timestamps and undergo the existing context/freshness checks. The audit records
the planned and live-selected candidate indices and rejected live capacities.

The fix does not make global `MemAvailable` a per-NUMA-node bound. During an
initial-chunk run the later host check can still double-count already resident
allocations, because the current execution bridge lacks exact ownership of
those host pages. That guard remains conservative. The change also does not
establish runtime prediction or multi-GPU throughput accuracy.
