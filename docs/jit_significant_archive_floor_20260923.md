# Significant-pair indexed archive floor

`native_layout_significant_archive_floor` prices the existing single-writer
NPZ service for host- and device-selected significant pairs. Each nonempty
source chunk creates at most one host-selected part; each nonempty device
selection block creates at most one device-selected part. The report requires
explicit retained-pair intervals, a durable fsync setting, current PGEN
identity and the same independent archive/page-cache/storage/NumPy-copy
prices as the detailed finite graph. It is constructed after useful output;
it does not infer future survivors from an early empty chunk.

For a tile of `M` unissued variants, `Ktile` phenotypes and chunk width `B`,
each host part retains at most `min(B,M)*Ktile` pairs. The device path uses
the source selector's bounded trait columns and variant strips, determined
by explicit `device_selection_max_cells`. Its maximum part size and possible
block count are computed without enumerating chunks. For retained interval
`[Hlo,Hhi]`, nonempty parts lie between `ceil(Hlo/max_part_cells)` and
`min(Hhi,possible_parts)`. These endpoints bound array payload, NPY/ZIP framing and
rewritten local headers. The model adds the writer's owned df copy, page-cache
copy, durable file bytes and per-part fsync. It retains the one-consumer
serial work across GPUs and combines archive CPU and logical host-memory
traffic with source and DMA work under shared capacities. Construction takes
constant work per phenotype tile, independent of the number of chunks.

The device output floor also counts the blocking four-byte CUDA nonzero
count transfer for every selection block, including empty blocks. GPU selector
kernels, count synchronization latency and host dispatch remain separate;
this report only bounds the archive service. The high
endpoint bounds a partial modeled writer floor, not completion time or a
JIT-switch benefit. Queue stalls, live in-flight output, filesystem throttling
and final metadata publication remain outside it.

The public significant scan requests beta even when `sumstats_fields='t'` and
the NPZ writer omits beta. `native_layout_output_floor` now keeps the full
beta-plus-t D2H payload in that case while counting only t and indexing fields
in the file payload. Previously its host significant D2H floor followed the
writer's `store_beta` flag and understated this transfer.

On A100, the focused host/device output, archive and partial-envelope batch
passed 38 tests in job `20260923-005617-1273488`. These are ledger and
equation checks, not a measured whole-job speedup.
