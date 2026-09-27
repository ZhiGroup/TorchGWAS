# Source-price binding at the live staged screen

The live first-chunk screen now audits its candidate PGEN source prices against
the one active detailed-profile context before it prices the whole unissued
header. For each candidate partition, the decoder units, CPU fraction,
prefetch/decode-worker settings, available CPU, DRAM and input rates must
exactly match that device's profile. Optional buffered-read CPU prices must
also match when the candidate supplies them; a candidate cannot omit an active
buffered-read price. The screen's shared source
capacities may be lower than the context's values as an explicit availability
scenario, but may not exceed them. A mismatch stops the optional screen and is
recorded without failing the scientific scan.

Matching values are only the first check. The audit also expands those source
rates into individual price leaves and checks whether each is covered by a
declared immutable measurement target in the already validated active profile.
It reports matched, declared and unbound leaf counts and at most 24 example
unbound paths. `source_prices_unbound` is a distinct live result status when
the issue/output checkpoint and profile remain current but source targets are
missing. A profile change or source/output advance retains its stronger stale
status. Measurement records and their original ages are checked by the
controller's existing detailed-profile validation before and after the
screen; this audit does not renew an old record.

The saved `results/public_initial_chunks_20260922/execution_v2/profile.json`
is a historical illustration: its one declared price binding targets the two
devices' resident NumPy copy coefficients, not the source coefficients this
staged screen consumes. Its source-value agreement alone therefore cannot
qualify a JIT decision. Current server measurements and bindings may differ;
the live audit uses the active profile rather than this saved file.

This increment covers only the source portion. H2D/GPU, selector, indexed
archive, dense writer, storage, final publication and a finite completion
model still need their own current price bindings and validation. The focused
A100 source-stage/controller/screen suite passed 121 tests, and the final
broader calculator/controller suite passed 320 tests. The tests use a
real PGEN header and synthetic typed output events; they establish the
binding logic and lifecycle, not loaded throughput or a profitable switch.
