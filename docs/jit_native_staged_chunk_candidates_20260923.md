# Native first-chunk candidate construction

The public `initial_chunks.staged_screen` option now registers a background
screen while the starting layout is prepared. Registration builds no PGEN
schedule and runs no candidate calculator. After useful written output and
completion of the bounded source-metadata ledger, the worker constructs up to
four chunk-width candidates from the current exact unissued frontier. Their
partition IDs, source suffixes, phenotype ranges and GPU ownership stay fixed;
each candidate changes only the future chunk width. The source, transfer,
GPU and output prices come from the active detailed context. Dense output
uses the live registered native writer's block, queue, borrow, fsync,
writeback and df-sidecar settings. Host significant and full-panel JAGWAS
use the already validated reduced-output primitive record. JAGWAS retains
the complete phenotype panel on each variant shard.

The configuration declares an explicit list of one to four admitted chunk
sizes including the starting size, one output occupancy scenario, limits for
partitions, unique PGEN records, chunks per partition, total CPU and wall
screen time, and zero to two frontier retries. It requires `source_staging`.
Missing shared H2D/D2H context prices or a missing live dense writer fail
only the optional post-output screen and are retained as screen evidence;
the scientific scan continues at its admitted starting layout. The
candidate factory does not infer bus capacity from per-GPU rates.

This establishes a real public route from first chunks to a measured,
frontier-bound calculator screen, while still keeping its result as evidence.
It does not change a running chunk size. Phenotype retile, GPU reassignment,
memory and writer-ownership admission, queued/in-flight completion, final
output durability, independent shared-transfer calibration and loaded-job
selection validation remain separate requirements for a JIT switch.
