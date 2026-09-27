# Source-complete read/decode floors for JIT continuation

`PgenHeaderWork.schedule_bounds` already enumerates all native PGEN chunk-start
LD replays from indexed metadata, without decoding payload. The new
`native_schedule_source_floor` prices its whole unissued schedule with the
existing independent decoder units and buffered-read primitives. It retains
the exact source identity, marker range, chunk count and read bytes. Decoder
CPU is a conditional interval; buffered-read CPU includes the minimum one
fixed call per chunk, with short-read retries unpriced. Logical DRAM work uses
the same equation as the per-chunk scan model,
including the extra LD base traffic. A missing possible decoder-unit price
fails closed.

For a declared shared-capacity scenario, the source floor is the maximum of
aggregate CPU/DRAM/input work divided by their capacities and a reader-worker
chain bound. The worker term adds minimum aggregate read and decode service,
then divides by at most `min(depth,decode_workers,chunks)` concurrently active
readers. The paired schedule conservation check and live file identity are
revalidated before pricing. These are necessary source-stage bounds in the
declared analytical model, not upper/lower elapsed-time bounds for a GWAS run.

`native_layout_source_floor` composes one such schedule per fixed executor
partition. It adds physical source rereads for disjoint phenotype tiles,
counts shared CPU/DRAM/input capacity once across GPUs, and sums reader-worker
chains for tiles serialized on the same GPU. Variant shards require distinct
devices and complete phenotype panels; JAGWAS rejects phenotype partitions.
The caller must still bind the rectangles to its exact issued/unissued
frontier and admit memory. Neither function prices GPU, transfers, selection,
output, in-flight state or final drain, so neither can authorize a switch.

The focused A100 source-floor tests passed 7 cases (job
`20260922-230634-1245992`). A broader source/resource/continuation batch
passed 249 tests in 28.54 s (job `20260922-230748-1246213`) before the
final identity recheck. The current source passed 19 focused cases (job
`20260922-231701-1247228`) and the subsequent aligned LD shard-boundary
case in a separate remote pytest command. These include
exact PGEN-census enclosure, source mutation rejection, shared-capacity
accounting, physical trait-tile rereads, and full-panel JAGWAS validation.

The read-only local `/data` 8,086,101-variant diagnostic is
`results/large_schedule_source_floor_v2_20260923/report.json` (job
`20260922-231843-1247527`, source hash recorded at execution). It conserved the prior
indexed byte/replay counts for 1,024 and 4,096 markers. Header construction
took 0.482 process CPU seconds; the ordered first 1,024-marker schedule took
0.903 CPU seconds and its source floor 0.185 CPU seconds. The later
4,096-marker schedule reused cached primary work and took 0.0215 CPU seconds;
its source floor took 0.0351 CPU seconds. Peak process RSS was 273,056 KiB.
These ordered timings are sensitive to shared-server load and cache state;
they are not matched speedup measurements. Primitive rates in this
diagnostic are synthetic arithmetic fixtures, so only calculator construction
cost and source conservation are evidence. The current cold path would still
be too slow in a first-output callback; run source-complete work in a charged
background step and reuse structurally bound results in later jobs.

The later [compute/frontier extension](jit_compute_frontier_floor_20260923.md)
adds mandatory H2D/matrix work and exact unissued association coverage. A
finite continuation still needs complete GPU and dense/significant/JAGWAS
output service, issued work and final drain. It must also confront the frozen
H100 absolute prediction error before a production performance switch is
justified.

The first output-regime extension is `src/torchgwas/layout_output_floor.py`.
It uses the fixed source partitions to count native result D2H payload and
mandatory output-array payload in O(partitions), with a shared output/D2H
capacity and a capacity per device. Dense results carry per-tile status and
df transfer; significant device selection carries status, selected pairs and
the four-byte blocking nonzero count for every selection block,
while host selection carries the dense result; JAGWAS carries 17 bytes per
marker before the indexed writer. Selected-row intervals must be declared and
default to the full valid range. ZIP headers, selector service, writer
CPU, H2D, GPU, queues and durable commit are still unpriced. Its payload
resource floor is not a whole-run prediction or switch criterion.
The combined source/layout/output regression passed 226 cases on A100 (job
`20260922-232424-1248276`); the final small readability change passed all 10
output cases again (job `20260922-232846-1249018`).

For the real 8,086,101-marker source with a 128-phenotype panel, the exact
native payload arithmetic is illustrative: dense beta+t output arrays require
8,280,167,424 bytes, with 8,320,597,929 bytes of result D2H payload including
per-marker status and df. JAGWAS transports 137,463,717 result bytes and has
at most 129,377,616 indexed array payload bytes, depending on retained valid
markers. Device-selected significant pairs range from 8,086,101 status bytes
and zero output-array bytes to 28,988,672,085 D2H bytes and 28,980,585,984
output-array bytes if every pair survives. These are byte counts, not timing
predictions; the wide significant interval is why first-chunk occupancy must
remain a scenario rather than a silent global extrapolation.
