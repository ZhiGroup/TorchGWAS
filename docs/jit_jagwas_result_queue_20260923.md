# JAGWAS shared-result queue in the JIT checkpoint

The productive JAGWAS multi-GPU path now registers its one owned, bounded,
unordered result queue before the producer threads start. A staging worker can
copy the queue references under its mutex at one anchor time and count result
ranges and resident NumPy payload bytes without reading or modifying results.
A source-issue/output revision token is checked on both sides of that capture.
After useful indexed output, the held exact source prefix is bound to completed
part events, and queued ranges are matched to issued-but-not-part-written
chunks and their admitted device owners. Empty completed parts count as
completion. A result currently delivered to the indexed writer is outside the
queue; it remains among the issued chunks outside the observed queue.

`refine_issued_jagwas_with_queue` splits the existing full-chunk issued-work
upper ledger at the same revision. Queued chunks have completed source read,
decode, H2D and GPU production of their host result. Chunks outside the queue
retain the full conservative producer workload because they may be upstream or
inside the active writer. The split conserves each workload count, and queued
host payload is recorded separately. The candidate checkpoint ledger accepts
this bound queue state as an optional JAGWAS-only refinement: its combined
source/H2D/GPU upper workload then excludes completed queued producer stages,
while its issued output-array backlog remains. It rejects a simultaneous
full issued-GPU shape estimate, which would double count those stages. It is an accounting refinement, not an
elapsed-time upper bound: indexed selection, part encoding/fsync, final
manifest publication, future source work and shared-capacity contention still
need a finite continuation schedule.

The earlier source-staging path ended the productive issue ledger after one
output. It now suppresses the legacy short-window layout switch while retaining
issued ranges across the first chunks, until the normal source range or wall
window budget ends retention. It also bypasses the old eight-partition
short-window gate; large phenotype-tile layouts can accumulate the same
read-only staging evidence. This makes a later first-chunk decision joinable to
an exact issue frontier. It does not enable a layout switch yet.

The 2026-09-23 A100 regression `20260923-051447-1370110` passed 182 tests for
the initial queue observation; `20260923-051814-1370374` passed 109 tests after
the retained-frontier change. The extended join/refinement regression
`20260923-052406-1370793` passed 165 tests; the queue-aware checkpoint
regression `20260923-052733-1371189` passed 154 tests. The subsequent staged-worker test on an actual PGEN header passed in
`20260923-053131-1371677` (12 targeted tests). No loaded multi-GPU speedup or absolute runtime calibration is
claimed from these state tests.

The [active indexed-writer increment](jit_indexed_writer_state_20260923.md)
now brackets the result consumed between queue pop and the post-fsync
`on_chunk_written` callback. A result outside the queue remains conservative
unless the same active writer range is observed on both sides of the queue
anchor. The completion event still reports fsynced part bytes (including empty
chunks with no part); `variant_ids.npy`, optional variant metadata, final
manifest and directory publication happen only after the iterator drains.
These phases need separate service and durability accounting before any finite
completion guarantee.
