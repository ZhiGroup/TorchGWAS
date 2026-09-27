# Compact JAGWAS indexed archive for a post-output JIT pass

The JAGWAS writer stores two uncompressed NPZ arrays in every nonempty
variant chunk: absolute int64 variant indices and FP64 chi-square values.
`native_layout_jagwas_archive_floor` prices this exact schema with the
existing finite graph's independent archive call, payload-byte, page-cache,
storage-byte and fsync services. The output layout must explicitly declare
durable `jagwas_writer_fsync=True`. No retained count is inferred from the
significance threshold or from an early empty chunk.

For each variant shard with `M` unissued variants, chunk width `B`, and
declared retained interval `[Hlo,Hhi]`, the number of nonempty NPZ parts lies
between `ceil(Hlo/B)` and `min(Hhi,ceil(M/B))`. Each file has `16H` array
payload bytes plus two NPY headers and ZIP framing. The smallest and largest
possible nonempty per-part framing are evaluated once with the existing
indexed-part ledger. This bounds file bytes, submitted header rewrites,
archive CPU, logical DRAM traffic and fsync calls in constant work per shard.
It also retains the one-consumer serial service bound: archive CPU/DRAM,
storage and fsync for every nonempty part.

The partial envelope now sums distinct decoder, selector and archive CPU
work under one shared CPU capacity, and distinct decoder, DMA, selector and
archive traffic under one shared host-memory capacity. JAGWAS host selection
and archive work also share one indexed consumer. Storage bytes include NPZ
framing in the archive floor; the older output payload floor still counts
only the 16-byte-per-retained-variant arrays.

The model remains a conditional partial floor. The upper endpoint bounds
only modeled archive work, not completion time. It omits queue stalls,
filesystem throttling, metadata publication, current in-flight output, and
the unpriced first-use behavior needed for a qualified finite continuation.
It cannot authorize a JIT switch by itself.

On A100, the broader layout/JAGWAS regression batch passed 141 tests in job
`20260923-003326-1262674` (two prefetch-depth warnings). After the sparse
occupancy cap was added, the focused archive/envelope batch passed 19 tests
in job `20260923-003441-1262899`. These validate the ledger and existing
service equations; they are not whole-job performance measurements.
