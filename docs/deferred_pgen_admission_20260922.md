# Deferred PGEN base lookup in public JIT admission

The public `initial_chunks` path now admits exact PGEN reader memory without
building a full non-LD base-position array on a cold job. Source opening still
parses the PGEN header first; this change removes a second whole-file pass in
the planner. It does not complete runtime-calculator calibration or dynamic
phenotype/device reassignment.

For each fine-grid chunk start, admission reads the already prepared record
offsets and lengths and finds the preceding non-LD base. It resolves ordinary
starts with at most 16 vector steps, then walks only unresolved starts locally
under a one-million-step cap; a pathological file falls back to the exact full
index. The fixed and shifted reader-memory maxima remain exact. The first
productive source window checks file identity and scope again, then uses local
LD predecessor lookup for its bounded starts. A cache hit can instead reuse
the retained base index. The short-remaining-work gate runs before constructing
any source window.

An optional structural cache is still keyed to source identity, fine-grid
size, dimensions and implementation hash. A cold miss retains only 2,021,544
bytes of fine-grid vectors on the tested large file. After successful output,
it derives the full base index and atomically publishes a bounded binary
artifact. The current v2 format aligns each read-only vector and the uint32
base array to eight bytes. The hot lookup passes a uint32 scalar to NumPy.
Both matter: the former unaligned artifact and a Python-integer search made
NumPy copy or promote a 24 MB base array for each predecessor query. The
entire retained artifact is charged to host memory on a hit; a cold miss
charges its pending vectors. Original component-price observation times are
never renewed by this structural cache.

On the real local 22,250-sample, 8,086,101-variant PGEN, the final ordered
component run measured 8.069 and 3.290 CPU seconds for the old full-base
admission, versus 0.890 and 0.066 for exact deferred-base admission. The
fine-grid vectors and fixed/shifted envelopes agreed exactly at capacities
128 and 1,024. Initial header parsing separately took 11.608 CPU seconds.
Those timings vary on the shared server and do not establish whole-job
speedup. After-job cache publication took 5.613 CPU and 5.868 wall seconds;
aligned reuse loading took 0.091 CPU seconds and reproduced the vectors
and sampled LD predecessors. In a separate ordered lookup check, 1,000
aligned dtype-matched searches took 0.0019 CPU seconds while one search
with a Python integer took 0.188 CPU seconds. These are diagnostic
component measurements, not independent throughput capacities.

The final public two-GPU JAGWAS audit produced identical control, cold
deferred and reuse output, with zero maximum numerical difference. The cold
job retained zero base-index bytes at admission and published the artifact
after completion; the reuse job loaded a 36-byte base index on this small
fixture. Both proposals cited the same original component-price observation
time. API times were 6.383, 5.515 and 4.546 seconds under varying warm and
server-load conditions. Most supplied prices were synthetic and no
profitable switch occurred, so these times do not qualify a production
autotuner. A separate public run without a structural cache also kept
zero base-index bytes at admission and matched its control output in both
deferred passes. Its API times varied with warm state and load; it is
correctness evidence for the default cold path, not a speedup comparison.
The relevant A100 regression batch passed 238 tests, including dense,
significant and JAGWAS admission/controller coverage.

Selectively pulled evidence is in
`results/deferred_pgen_handoff_20260923/`: `large_report.json`,
`predecessor_lookup.json`, `public_jagwas.json`,
`public_jagwas_no_cache.json`, `parser_comparison.json` and `source_sha256.json`. Six current package
source hashes matched that remote report after the pull. The rejected
global-assembly parser experiment produced identical header arrays but did
not improve the original parser in an ordered comparison, so the PGEN reader
was left unchanged.

The remaining cold-start floor is the source header parse, which materializes
the full record-length and offset arrays before the native reader opens.
Reducing it requires a reader/index design change while preserving exact
scope, LD replay, source invalidation and memory bounds. Separately, the
calculator still needs independently measured full-pipeline services and
output-inclusive qualification before it should authorize automatic chunk,
tile or multi-GPU decisions.
