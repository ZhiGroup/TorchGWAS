# Source-bound PGEN admission cache

Historical v1 evidence. Current behavior and v2 validation are in
[deferred PGEN admission](deferred_pgen_admission_20260922.md).

This change reduces repeat-job startup work for the public `initial_chunks`
path. It does not qualify the analytical runtime calculator or finish JIT
autotuning.

When `initial_chunks.structural_cache_dir` is set, memory admission can load
the exact compact PGEN reader envelope and non-LD base positions from a prior
successful job. A miss still uses the header-only calculation. The first job
holds the fine-grid vectors in memory and publishes them **after successful
output**, so disk publication does not delay the first useful chunk. The cache
never contains duration prices, measured rates, available memory, or a proposed
chunk size. Existing component measurements retain their original observation
timestamps.

The artifact key binds the resolved PGEN path, byte size, mtime, ctime, device
and inode; sample and variant counts; fine-grid size; and the implementation
hash of `pgen_memory_layout.py`. A one-second timestamp-stability gate prevents
reuse immediately after an input edit. The bounded binary artifact has a
manifest and SHA-256 payload checksum, is written through an atomic replacement,
and stores at most eight entries of at most 256 MiB each. Corrupt, oversized,
stale and source-changed entries become cache misses. Cached arrays are
read-only views of the verified immutable bytes. The complete retained blob,
including manifest and fine-grid vectors, is charged to host memory; a miss
charges its retained base array and pending vectors. Live host/GPU capacity and
the full input/profile bindings are still checked before execution and during
productive decisions.

A direct A100 check used the real local
`/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen`:
N=22,250, M=8,086,101, 20,838,552,600 file bytes, fine grid 128.
The retained base index was 24,322,956 bytes; fine vectors were 2,021,544
bytes; the artifact was 26,344,746 bytes. Cached vectors and base positions
were byte-equivalent to the fresh build, and fixed and shifted reader
envelopes matched at capacities 128 and 1,024. In the final ordered component
run, the fresh compact build took 4.916 CPU seconds, publication 0.498 CPU
seconds after completion, and reuse load 0.104 CPU seconds. The separately
parsed header took 3.112 CPU seconds. Earlier ordered observations varied
considerably under shared-server load, so these are component costs, not a
whole-job speedup or a bound on cold latency. The source/header parse remains a
material cold-start cost.

The public two-GPU JAGWAS audit produced identical control, first deferred and
reuse outputs (maximum numerical difference zero). The first deferred job
reported a cache miss and successful publication; the next reported a hit.
Both productive proposals referenced the same original component-price
observation time. The selected chunk size remained unchanged under largely
synthetic prices; API times were 3.359, 1.862 and 1.101 seconds for control,
deferred and reuse under varying warm/load conditions. They are not
performance-comparison evidence.

Remote evidence was selectively pulled to
`results/admission_cache_handoff_20260922/large_pgen.json`,
`public_jagwas.json`, and `source_sha256.json`. The final relevant A100 regression batch passed 235 tests, including
unsupported-record-type rejection and a short-remaining-work JIT gate. The
public execution audit ran against the final controller/cache implementation.
An earlier attempt in this final source stopped before execution because its
resident-copy samples failed the declared drift check; a fresh result directory
held the successful audit. No unstable price was accepted.

Remaining cold-start work is distinct: the original PGEN header parser
materializes record lengths and offsets for the source reader before the
planner can admit memory. On the large input, its measured CPU time varied
from about 3 to 12 seconds across ordered observations; line attribution
showed most cost in writes to those large arrays. The cache does not eliminate
that first-job parse. A future optimization needs exact file-scope and native
reader correctness, per-source invalidation, and memory charges. The JIT bridge
still changes only chunk size after useful output; phenotype retile/device
reassignment and independently priced full-pipeline calibration remain open.
