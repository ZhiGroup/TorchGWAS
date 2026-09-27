# Device significant-count transfer floor

`native_layout_significant_device_count_floor` carries the existing
`device_selection_graph` price of one blocking four-byte CUDA nonzero count
transfer into the compact, source-bound layout calculator. For `S` selection
blocks on a GPU, the necessary serial count-transfer time is
`S * (latency_seconds + 4 / bytes_per_second)`. Successive phenotype tiles
on that GPU are summed; different GPUs can overlap, so the layout floor is
the maximum per-GPU sum. The block census uses the current selector's source
chunk, trait-column and variant-strip geometry without expanding all blocks.
The [shared geometry helper](jit_device_selection_geometry_20260923.md)
now chooses a minimum-block bounded shape for each full or tail chunk.

The D2H output floor already charges the four bytes per block to shared and
per-GPU bandwidth. The new floor adds only the per-transfer latency chain;
its bytes are not added to host-memory or D2H load again. It requires an
independent count-transfer price for every active GPU and binds the exact
output partition digest and selection cell limit before entering the partial
envelope. An empty block still incurs a nonzero count and host barrier.

This is one necessary service floor, not the selector's complete execution
time. Count kernels, flagged selection, host API calls, active CPU waits,
selected-payload copies, queueing and final drain still need a finite
continuation and loaded validation before JIT can compare layouts safely.

The A100 selected-output, selector-graph, writer, JAGWAS and productive-JIT
regression passed 253 tests in job `20260923-010907-1283679`. The compact
count floor was also checked against the exact source selector graph.
The final price-snapshot and envelope validation passed 19 focused tests in
`20260923-011119-1283981`.
