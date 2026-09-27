# Transfer, GPU and output price binding for staged JIT evidence

The first-chunk live screen now audits the work prices in every supplied
candidate after its source-price check. Per-GPU H2D, D2H, FP32 and JAGWAS FP64
ceilings must match the active device profile. Shared H2D and D2H ceilings
must come from positive `h2d` and `d2h` values in the optional
`context.shared_transfer_capacities` field; if that field is absent,
the audit reports it as unbound rather than deriving a bus rate from GPU
identity. The output capacity may be lower than its context value as a named
availability scenario, but cannot exceed it. Dense-writer copy, zeroing,
writeback, executor and fsync coefficients are also compared with the active
profile. A different available value aborts only the optional screen.
If an active context declares shared transfer links, a candidate must retain
them. A multi-GPU candidate without declared link topology is marked unbound.

For reduced output, the audit compares the active per-device selector/archive
profile fields used by JAGWAS and significant pairs. Host significant and
JAGWAS selection/archive primitive banks may additionally match the exact
value of the already loaded, age-checked reduction record. That record's
content and artifact SHA-256 identities are retained in the screen audit.
Device-selected significant archives can use that same exact record, while
their count-transfer prices still need a separate binding. A supplied GPU
shape profile must match the active fixed profile, and every supplied compiled
geometry row must occur exactly once in the active bank. The CPU fraction,
GPU resource rates and host dispatch primitive prices consumed by the tensor
shape service must have immutable measurement targets. Compiled geometry is
checked by exact identity rather than treated as a timing measurement. An
absent shape service is reported explicitly.

The audit records matched and declared leaves, bounded examples of missing
targets, and rates absent from the active context. A current live screen with
matching but undeclared work prices reports `work_prices_unbound`; a stale
issue/output frontier or profile retains the stale status. This is still a
conditional partial resource screen, not an output-inclusive completion
estimate. Source, work and reduction-record identity checks cannot establish
available capacity under load or authorize a chunk, tile or GPU switch.
