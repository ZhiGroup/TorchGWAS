# Compact source admission for the first useful chunks

The public `initial_chunks` path now admits the fixed starting GPU/phenotype
layout and every permitted future chunk size from compact, read-only PGEN index
vectors. The ordinary full per-chunk layout remains available for audits and
other callers. This changes memory admission only: decoder work and runtime
pricing still require bounded real source windows after output starts.

For a fine-grid start, the index vectors hold its exact contiguous payload,
LD-base read prefix and restart scratch. A cumulative payload array gives the
sum for any capacity-wide interval. Fixed-capacity reader memory takes maxima
over regular starts; adaptive reader memory takes maxima over **every** aligned
fine-grid start in each admitted variant shard. Both retain separate read and
scratch maxima, matching the existing grow-only reader allocation model. The
compact route never substitutes a sampled window for an exact memory bound.
Dense output can keep fixed trait tiles or variant shards, significant pairs
keep fixed trait tiles, and JAGWAS keeps the complete phenotype panel on each
variant-sharded GPU.

The [server-local comparison](../results/compact_jit_admission_20260922/report.json)
used the 8,086,101-variant, 22,250-sample hard-call PGEN, fine size 128 and
capacity 1,024. Its 63,173 logical fine chunks occupied 2,021,544 bytes of
compact vectors. Fixed and shifted read/scratch maxima and total payload were
identical to the full layout for the whole file and both variant shards.
Within that one ordered process, compact/full fine-layout construction took
6.144/8.671 CPU seconds; whole-file envelope calculation took 0.018/0.633
seconds. Initial index parsing took 13.201 CPU seconds. These component timings
depend on load, cache and allocation state and do not establish a whole-job
speedup or an acceptable cold first-output latency. The benchmark script and
three package source hashes match the pulled report.

The small-fixture differential batch passed 48 tests across fixed and shifted
reads and dense, significant and JAGWAS startup. A public two-GPU JAGWAS audit
under `results/public_initial_chunks_20260922/compact_execution_v4/` then
wrote bit-identical control/deferred/reuse outputs; it kept the initial chunk
size with mostly synthetic prices. Its deferred/reuse startup validation took
1.05/0.21 seconds and productive planning took 2.23/0.55 seconds. Shared-server
warm/load states differed, so those times do not isolate the compact change.
The consolidated A100 regression batch passed 229 tests in 64.31 seconds.

Cold startup still includes full PGEN index parsing and full-array validation;
the compact representation removes Python row materialization and repeated
per-chunk envelope loops. Next work should reduce or defer those remaining
full-file checks, then validate real output-inclusive JIT decisions with
independently measured component prices and an observed remaining-work gate.
