# Active indexed-writer state for productive reduced-output checkpoints

Synchronous indexed writing for JAGWAS and significant pairs now has a
compact optional state:
waiting for a result, selecting, emitting a part, writing final metadata,
publishing the manifest, or published. It retains only one source range and
producer identity while a result is active, plus monotonically increasing
completion counters. It holds no result arrays. The state is cleared after a
part has finished optional fsync and before the existing chunk-completion
callback; an empty chunk clears it without creating a part. A failed chunk
write marks the state failed. The scientific output path is unchanged when
there is no productive observer.

The productive controller registers that writer weakly. For multi-GPU
JAGWAS, it snapshots the writer before and after an atomic shared-queue
snapshot and checks the source-issue/output revision around the sequence. If the same unique active
source range spans both writer snapshots, that result was already consumed by
the writer at the queue anchor. The queue/frontier binding verifies its exact
issued range, full phenotype panel and device, and excludes it from the queued
set. The issued-work refinement can then deduct source read/decode, H2D and GPU
work for both queued and stably active results. Chunks outside those observed
states retain full conservative producer work. The checkpoint ledger preserves
output-array backlog for every unfinished part, including the active result.

This is a workload refinement, not a completion-time prediction. Selection,
NumPy part encoding, possible fsync, result arrays in transit, final variant
metadata, manifest/directory publication, shared-capacity contention and
unissued work still need finite service bounds. An active writer result is
never treated as part-complete before its existing callback. The
writer's metadata/publication phases are observable but not yet priced. Significant-pair phenotype tiles retain their producer identity in the
writer state, but the host selector and its upstream queue remain unobserved.
The significant writer snapshot is read-only evidence and is not yet joined
to a finite completion-time calculator.

The A100 regression `20260923-053829-1372074` passed 164 tests including a
held `np.savez` write and empty-part behavior. The queue-aware checkpoint
follow-up `20260923-054036-1372287` passed 2 tests. The held
writer/queue anchor integration passed 13 tests in
`20260923-054135-1372387`. These tests do not validate loaded throughput
or the frozen H100 absolute runtime model.

A significant-pair writer can obtain its producer identity from a typed
`PartitionedIndexedChunk`; it need not invent a full-panel range resolver.
The source-staging worker captures its active tile and publication phase at
a stable issue/output revision. No selector or writer service is inferred
from the phase alone. The related A100 regression
`20260923-054509-1372690` passed 270 tests and 5 subtests, including tiled significant output.
