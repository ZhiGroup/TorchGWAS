# Native PGEN compute floor and unissued-work coverage

`native_layout_compute_floor` adds compact H2D and matrix-product work to the
existing read/decode and output-payload floors. It is restricted to the
explicitly unpacked, hardcall int8 PGEN path (`TORCHGWAS_PGEN_PACKED=0`). For
each fixed partition with `M` unissued markers, `N` samples, `K_t` phenotypes
and `R` covariate-Q columns, the native scan transfers `M N` genotype bytes
and performs a FP32 product of `2 M N (K_t + R + 1)` FLOPs. The extra column is
the intercept constructed by `native_scan.py`. JAGWAS additionally projects
the complete phenotype panel in FP64, `2 M K²` FLOPs per variant shard. The
module sums repeated genotype uploads for phenotype tiles, counts one shared
H2D ceiling across GPUs, and retains per-device H2D and matrix ceilings. These
are necessary resource-load floors *if the supplied capacities are valid
service ceilings*. They omit conversion, all other statistics/reduction
kernels, phenotype/setup transfer, selectors and final drain. An independent
profile and live admission are still required; no default peak rate is guessed.

For the real 8,086,101-marker, 22,250-sample source, one unpacked pass alone
uploads 179,915,747,250 genotype bytes. At 128 traits and 27 Q columns, its
scan GEMM performs 56,133,713,142,000 FP32 FLOPs; full-panel JAGWAS adds
264,965,357,568 FP64 projection FLOPs. Splitting the 128 traits into two tiles
would upload the genotype twice, while two disjoint variant shards would not.
These are operation counts, not runtime predictions or a tile recommendation.
The source-layout composer now admits up to 10,000 declared partitions by
default and checks trait/variant overlap after sorting intervals. This permits
many phenotype tiles for voxel-scale panels without an O(tiles²) overlap pass;
the explicit maximum still bounds calculator work.

`unissued_frontier` reads one active `ProductiveTuningRun` snapshot after actual
written output. Its cursor follows every reserved read, including prefetch and
in-flight reads. The only movable work is each partition's suffix from that
cursor to its end. It first checks that the starting partitions cover exactly
the declared full job variant range and phenotype panel.
`bind_layout_to_frontier` then sweeps variant boundaries and checks exact
phenotype coverage on each strip, rejecting missing or duplicated
variant–phenotype pairs and any attempted reassignment of the issued prefix.
Different phenotype tiles may have different cursors; their unissued
rectangles remain distinct, so a proposal that rewinds one tile is rejected.
The proof is bounded by a declared active-rectangle-visit budget, source
identity, issued revision and writer-event count. It does not move a live
reader, result queue, phenotype factor or writer, or certify a stale snapshot
at application time. JAGWAS retains the full phenotype panel per shard.

`native_layout_partial_envelope` binds the exact unissued coverage proof to
matching source, compute and output reports. Its interval is the maximum of
the three necessary stage floors under their declared capacities. The high
endpoint is an upper endpoint on this *partial floor*, not an upper bound on
completion time. This makes a compact lower-bound screen possible without
expanding every source chunk into a graph. It still cannot rank feasible
candidates or authorize a switch: the calculator needs complete GPU, selector
and writer service, queued work and final drain, plus a conditional baseline
completion ceiling and held-out executor qualification.

The A100 batch `20260922-233808-1250414` passed 251 source, output, compute
and productive-frontier tests on the first implementation. The subsequent
sample-binding and full-job-coverage changes passed 254 related tests (job
`20260922-234126-1250973`). The first partial-envelope binding passed 261
related tests (job `20260922-234348-1251230`). The final dense,
empty-first-part significant and full-panel JAGWAS integrations passed 263
related tests (job `20260922-234514-1251508`). After the many-tile composer
change, 59 focused source/layout/continuation tests passed (job
`20260922-234622-1251691`). A final unequal-tile-cursor check passed in the
265-test related batch `20260922-234808-1251908`. These are structural and
output-bound checks with fixture capacities, not held-out runtime
qualification.

The [shared transfer-link extension](jit_multigpu_transfer_floor_20260923.md)
uses the existing finite graph's `shared_links` schema for overlapping H2D
and D2H link constraints. The partial envelope additionally combines source
DRAM work with mandatory genotype DMA reads and result DMA writes against
one shared host-memory capacity. These are still necessary conditional
resource floors, not completion-time ceilings. The final related batch passed
71 tests on A100 (job `20260923-000737-1255007`).
