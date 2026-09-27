# Written-output boundary for first-chunk JIT

The productive controller now records bounded, per-partition output progress
alongside exact reserved source ranges. Indexed significant/JAGWAS events
identify their producer and absolute source chunk; a bound snapshot lists
reserved chunks that have not completed the indexed writer. Empty indexed
chunks still count as completed work. Dense writer events identify a
per-partition beta/t prefix in store coordinates, which the observer maps
back to source variant coordinates; the df sidecar prefix remains separate.
The bound report also counts issued variant–phenotype pairs not yet written
per partition, without treating them as completed GPU or queue service.
`productive_output_backlog` now carries those issued-but-unwritten partitions
into a conditional payload ledger. Dense beta/t bytes use the fixed pending
pair count, with an optional df sidecar count. Significant-pair and JAGWAS
pending chunks use the same global
association-coordinate occupancy scenario as unissued candidate layouts, so
splitting a chunk does not change its modeled survivors. The result is an
upper array-payload workload at the held checkpoint, not an upper elapsed
time: part writes or one dense array can already be ahead of the event.
When supplied to `productive_partial_floor`, this ledger is added only to
the unissued array-payload work upper; it never raises the necessary resource
floor or licenses a chunk-size decision.

This makes the post-output frontier more precise than a global written-event
counter, without retaining a job-sized history. Malformed or unbound optional
events mark the boundary invalid and appear in the audit; they do not change
the scientific scan. An invalid bound event now prevents the current
short-window planner from changing chunk size; valid forecasts still use
explicit caller boundary adjustments. A future finite continuation must
require a valid boundary, account for issued chunks not yet written, and
separately bound work in GPU streams, result queues, dense writeback/fsync,
and final manifest/directory publication. A written prefix is not proof of
durable completion.
