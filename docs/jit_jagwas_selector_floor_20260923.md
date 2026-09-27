# Compact full-panel JAGWAS selector floor

`native_layout_jagwas_selection_floor` prices the existing host JAGWAS
selection primitives for fixed unissued variant shards. It requires a current
PGEN identity, full phenotype panel on every shard, explicit retained-variant
intervals and the same independently measured per-primitive CPU prices used by
the finite graph. It needs constant work per shard: full chunk count, a
possible short tail, and the retained-count interval. It does not expand the
8-million-variant source into a chunk graph.

Each chunk always converts FP32 to FP64, checks finiteness, runs nonzero, and
dispatches index-add and gather. Retained counts only affect the latter two
unit terms and nonzero output traffic. If a scenario permits both empty and
nonempty chunks, the lower/upper selector-work endpoints use the cheaper/more
expensive corresponding nonzero primitive independently for each chunk. That
interval is deliberately conservative because a shard-level retained total
does not determine which chunks are nonempty. Exact empty and full scenarios
select one primitive. The selector floor is the maximum of shared CPU work,
logical DRAM work and the sum of CPU service through the public single
consumer. The high endpoint bounds this partial selector floor only.
The partial envelope also adds selector logical DRAM bytes to decode and both
DMA directions under the same declared host-memory capacity.

The differential tests compare the compact CPU/DRAM totals with the existing
`jagwas_host_selection_service` on empty, full and mixed per-chunk counts. A
first test found and corrected an accidental double charge of retained-only
gather units. The output-partition digest and source geometry are checked
before the floor enters `native_layout_partial_envelope`.
The A100 layout/JAGWAS regression batch passed 133 tests in job
`20260923-002451-1261842`; after the shared-DRAM addition and short-tail
fixture, the focused batch passed 20 tests in job `20260923-002620-1262116`.

This floor excludes GPU projection, indexed archive serialization, the
bounded result queue, final commit and current in-flight work. It does not
qualify a chunk-size or GPU switch. The public JIT still needs a source-matched
finite completion bound and a measured profitable same-job transition.
