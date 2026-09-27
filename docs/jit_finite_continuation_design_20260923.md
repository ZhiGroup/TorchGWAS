# Finite, bounded continuation for JIT chunk and layout decisions

The current productive planner compares three short source windows and rejects
an extrapolation that exceeds its declared ratio. That is appropriately
conservative, but an 8M-variant job cannot be decided from a few early chunks.
The whole-header PGEN schedule now supplies a cheap, source-complete work
ledger after the first useful output; it does not supply GPU or output time.
The following is the next calculator implementation target, not a switch
qualified by current measurements.

Let a discrete layout `x` select a chunk size, phenotype tile (where allowed),
variant/trait partition, active GPU set and fixed queue/decode settings from
memory-admitted candidates. For JAGWAS, every partition has trait range
`[0,K)` and a complete phenotype factor; only variant ranges are sharded.
Each candidate must cover every unissued association exactly once. Issued
ranges, in-flight buffers, existing writer staging and first output are a
common checkpoint, not work that can be reassigned retrospectively.
The [output-boundary audit](jit_output_boundary_20260923.md) now records
completed indexed source chunks or dense beta/t and df prefixes per
partition against that issue frontier. It does not yet observe GPU streams,
producer queues or final durability. The [dense writer queue observation]
(jit_dense_writer_queue_20260923.md) now records stream-local accepted,
staged, queued, active and fully written bytes at completed-prefix callbacks;
the streams of one writer are captured atomically, but the scan, GPU and other writers are not.
A [live two-pass bracket](jit_live_dense_writer_bracket_20260923.md) now bounds
accepted-minus-written bytes across registered writers at a common instant
and prices the accepted-write upper workload after useful output. The conditional
output-backlog ledger
also counts upper array payload for issued chunks without writer completion,
using the same global occupancy scenario as future unissued work. Partial
writes mean this payload ceiling is not an elapsed-time ceiling.
The [issued-work ledger](jit_issued_work_20260923.md) now carries whole
original chunks without observed matrix/part completion into bounded
read/decode/H2D/GEMM counts. It intentionally overcounts any work already in
flight and does not price a completion ceiling.
The [checkpoint join](jit_checkpoint_ledger_20260923.md) binds issued and
candidate unissued ledgers to one revision and survivor scenario. It preserves
original GPU ownership for issued chunks even when future work is reassigned.
Its summed counts cover only a subset of required operations, so the join
still cannot authorize a runtime comparison. The later
[JAGWAS queue join](jit_jagwas_result_queue_20260923.md) and
[active indexed-writer observation](jit_indexed_writer_state_20260923.md)
can now identify queued and stably active JAGWAS results as producer-complete
at one issue/output revision. The checkpoint ledger deducts their completed
source/H2D/GPU workload while retaining output backlog. Significant-pair
phenotype tiles now report the active indexed-writer range, but their upstream
host selector and queue are not yet joined. None of these observations prices
the remaining writer/final-publication service or supplies a finite elapsed
completion ceiling.

For each explicit source/output/capacity scenario `s`, form a finite ledger of
per-chunk read, decode, H2D, statistics/reduction, D2H or device selection,
host selection, queue and writer work. Indexed PGEN metadata supplies exact
record bytes and conditional decoder-unit intervals for every chunk-start LD
replay. The existing independently measured primitives price GPU, transfer,
CPU, DRAM and storage work. Dense output bytes and block/writeback thresholds
are deterministic; significant-pair survivors require bounded occupancy
scenarios that include the observed first-chunk regime; JAGWAS output has at
most one row per marker because invalid statistics are omitted. Do not infer
future significant-pair counts from a null or an
empty early window alone.

For a candidate, a resource-load floor is
`L_s(x) = max_r W_r(x,s)/C_r(s)`, augmented by each GPU's own serial/dependency
chain. Shared CPU, DRAM, input, output and PCIe-link capacities appear once
across all GPUs; duplicating a shared capacity per shard would manufacture
scaling. A finite graph or a conservative serial schedule gives a conditional
model ceiling `U_s(x)` including queued work and final drain. Neither endpoint
is a hardware-time guarantee when component prices or live availability are
wrong. The optimizer minimizes the scenario completion envelope; balancing
pipeline stages emerges from the shared-resource makespan, rather than from
forcing equal measured stage wall spans.

The first pass uses the compact whole-source schedule and aggregate resource
loads to reject dominated candidates without expanding thousands of Python
chunk dictionaries. If a candidate can plausibly win, expand bounded groups
of exact chunks in background after output, with a memory and CPU budget.
Group boundaries must include source tails, LD restart changes, dense
writeback thresholds and output-occupancy changes. Stop refinement when the
conditional decision interval separates or the budget expires. If it does
not separate, retain the current layout. The measured 8.1M-header vectorized
expansion still costs seconds and hundreds of MiB, so it cannot sit in a cold
first-output callback by default.

The [source-stage floor](jit_source_resource_floor_20260923.md) now implements
the indexed read/decode and shared multi-GPU source-work portion of this first
pass. `native_layout_output_floor` adds a constant-size payload ledger for
native dense, host/device significant-pair and JAGWAS output on the same fixed
partitions. It counts shared output and D2H capacities once, retains per-device
D2H constraints, and requires explicit intervals for selected rows (defaulting
to the full valid range). Device-selected output also counts the blocking
four-byte CUDA nonzero count transfer per selection block, including empty
blocks. It excludes ZIP/framing, selector service, writer CPU, storage commit
and intermediate shared links. Its endpoints bound only necessary payload
resource loads, not elapsed completion. Neither floor
is used to reject or apply candidates yet: complete shape-specific GPU and
H2D service, output service and a conditional baseline completion ceiling are
still needed for a sound whole-pipeline comparison.

The [post-output source adapter](jit_productive_source_floor_20260923.md)
now binds that compact read/decode floor to the exact unissued frontier and
composes a necessary source/compute/output resource floor.
The [incremental source schedule](jit_incremental_source_schedule_20260923.md)
now distributes exact header work over bounded post-output steps and rebases
it to later source cursors, chunk sizes and variant-shard endpoints. The
productive source adapter can use that staged work without a full direct
rescan, while reporting its earlier CPU/wall cost separately. The public JIT
controller can now drive bounded steps after useful output when explicitly
enabled, with an evidence-only audit. It does not yet charge these steps in a
finite continuation or use them for a switch; enabling the evidence path
suppresses uncharged short-window switches.
The [whole-layout sparse ledger](jit_whole_layout_sparse_20260923.md)
assigns conditional retained output to the same global association
coordinates across alternative tile and GPU layouts. These are inputs to
the finite continuation, not a public JIT decision.

The [compute and frontier extension](jit_compute_frontier_floor_20260923.md)
adds the explicitly unpacked hardcall PGEN H2D payload, FP32 scan GEMM and
FP64 JAGWAS projection as necessary per-device/shared capacity loads. It also
checks that a proposed layout covers exactly the unissued association pairs
from a productive snapshot. `native_layout_partial_envelope` binds the three
necessary stage floors to this exact unissued coverage. Their maximum remains
only a partial lower envelope; no production switch uses it.

The [compact GPU shape service](jit_gpu_shape_service_20260923.md) now
aggregates the existing exact tensor component prices over full and tail
chunks after output, with one price per distinct device/shape. The productive
checkpoint checks candidate and original issued ownership and retains both
conditional shape services separately from necessary floors and work counts.
Producer stage state is partly observed for JAGWAS; a finite queue/output continuation remains missing.

The [compact dense-writer ledger](jit_dense_writer_compact_20260923.md) now
counts exact fixed-chunk array payload, staging copies, minimum writes and
writeback requests without materializing all chunk events. The optional
df-sidecar payload enters the output floor only when the caller supplies
explicit matching writer settings. These work counts still need independently
measured service and a finite continuation that includes queued and in-flight
output. The compact dense-writer service below supplies the existing finite
graph's conditional CPU, page-cache, writeback and fsync prices.

The [shared-link extension](jit_multigpu_transfer_floor_20260923.md) carries
the finite graph's overlapping PCIe-link declarations into the compact H2D
and D2H floors and charges source DRAM plus both DMA directions to one host
memory ceiling. This makes the first-pass lower floor sensitive to GPU
placement and shared host traffic without assigning a link rate from device
identity. It remains insufficient to rank or switch complete layouts.

The [compact dense-writer service](jit_dense_writer_service_floor_20260923.md)
adds source-bound CPU-equivalent writer work and per-GPU serial write/close
chains to the partial envelope, using the existing finite-graph price fields.
The writer work and profile must be current, and a legacy copy-price fallback
explicitly reports its missing per-call price. This strengthens one necessary
floor without providing a finite completion ceiling; live dense-writer queue observations remain separate.

The [JAGWAS selector floor](jit_jagwas_selector_floor_20260923.md) now counts
the fixed per-chunk host primitives for full-panel variant shards, plus a
bounded retained-count scenario, without expanding the entire source. One
consumer serializes that host work across GPUs. The indexed archive floor
below supplies selected-writer service; live queue and active-writer state are now observable but not yet priced as a finite continuation.

The [indexed archive floor](jit_jagwas_archive_floor_20260923.md) adds
nonempty NPZ part counts, framing, independent archive service, storage and
fsync for explicit JAGWAS retained intervals. It combines selector and writer
work under their single consumer, and sums their distinct CPU and host-memory
traffic with source and DMA work at the shared capacities. This is still a
necessary partial load, not a qualified upper completion bound.

The [significant archive floor](jit_significant_archive_floor_20260923.md)
counts the NPZ service for both host- and device-selected pairs. Host selection
can write one part per nonempty source chunk; device selection can write one
per nonempty selection block. The compact bound uses the source selector's
chunk/block geometry and retained-pair intervals, including optional beta,
framing, owned df copy, durable storage, fsync and shared CPU/host-memory
loads. The device selector's GPU kernels, blocking synchronization service and
host dispatch remain unpriced. A finite queue continuation remains necessary.

The [host significant selector floor](jit_significant_host_selector_floor_20260923.md)
now uses the active NumPy nonzero protocol and existing primitive prices to
bound empty, sparse and dense work for complete trait tiles without chunk
expansion. It joins the host indexed archive under shared CPU/DRAM loads,
while retaining per-GPU producer selection and the writer's separate single
consumer. A finite queue/in-flight continuation and qualified completion
ceiling remain necessary before the public JIT may compare layouts.

The [device count-transfer floor](jit_device_count_barrier_floor_20260923.md)
now prices the blocking four-byte nonzero transfer latency for every device
selection block. It sums successive tiles on each GPU and takes the maximum
across GPUs; the bytes were already charged to shared/per-GPU D2H capacity.
The count kernel, rest of the selector and finite queue remain unresolved.
The [device selection geometry change](jit_device_selection_geometry_20260923.md)
also reduces blocking calls and possible NPZ parts for some huge trait tiles;
source and calculator use one block-shape helper. Matched A100 selector-only
controls support this shape, but do not establish whole-job throughput or a
JIT tile decision.
The memory admission path now uses the same shape, and the current A100
geometry and allocation censuses cover the worked 512-by-600,000 selector
case. This repairs one necessary safety bound; it does not provide selector
service, output-inclusive continuation or an admissible JIT switch.

For a permitted chunk-only move, apply it only to the held unissued frontier.
Require every scenario's baseline lower completion minus candidate upper
completion to exceed charged planning CPU/wall, measurement, publication,
switching and reserve cost. Require fresh profile bindings, live memory
admission and an output regime that first-chunk observations have not refuted.
An optional background proposal is discarded if source issuance or writer
progress changes its checkpoint before application. Tile or GPU reassignment
needs a separate bounded transition and output-ownership proof; the current
public JIT bridge does not perform either move.

The next test should first falsify this continuation against held-out native
PGEN controls, then compare exact output and profitable switches in dense,
nonempty significant and full-panel JAGWAS jobs. The frozen H100
15.862-second model prediction versus 55.769-second observed median remains
an absolute-calibration failure; resource-load floors or first-chunk loaded
spans alone must not be used to paper over it.

The [bounded staged partial screen](jit_staged_partial_screen_20260923.md) now implements the first post-output multi-candidate comparison against one completed source stage and exact unissued frontier. It includes the available mode-specific writer/selector floors and reports stage cost once. The real 8.1M-header metadata-only probe shows this work can be spread over first output events, but it still lacks a finite output-inclusive completion ceiling and cannot authorize a layout switch.

The [compact device-selector launch floor](jit_device_selector_launch_floor_20260923.md) adds mandatory CUDA launch service to the device-significant partial envelope, including empty nonzero blocks and conditional coordinate scatter. It joins per-GPU count-transfer latency but does not price all selector kernels or a finite completion ceiling.

The [large-K matched metadata control](jit_large_k_staged_screen_20260923.md) demonstrates bounded six/twelve-tile significant-pair screening at 600,000 phenotypes. Reusing one identical source-price result per span/profile reduced same-stage median screen CPU from 5.236 to 1.336 seconds without changing conditional work. This solves duplicate calculator work, not output-inclusive completion or live switching.
