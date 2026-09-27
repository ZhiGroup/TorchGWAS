# Sparse reduced output in first-chunk JIT comparisons

The existing productive comparison priced only empty or fully retained
significant-pair and JAGWAS output. Sparse output can incur archive and fsync
work per nonempty part while transferring far fewer result rows than the dense
case, so its service need not lie between those endpoint timings.

An `initial_chunks` run may now declare a sparse `occupancy_scenarios` value,
for example `{"retained_fraction": [3, 1000], "placement": "spread"}`.
The rational fraction is conditional retained association pairs for significant
output, or retained variant rows for JAGWAS. `spread` distributes retained
cells across the smallest admitted chunk grid; `clustered` groups them within
repeating denominator-sized association periods. Both now use deterministic
global association coordinates, so a longer planning window, a later
unissued cursor, or an alternate phenotype tile sees the same scenario.
Empty and dense string scenarios remain available.

Partial survivor bins are emitted on the smallest admitted chunk grid, at
most 256 bins across the bounded productive comparison. Both chunk candidates
then sum the same fine bins into their source chunks. This preserves total
retained associations while allowing the indexed archive model to change part
count, framing, CPU, storage and fsync work with chunk size. JAGWAS bins still
carry the complete phenotype panel, with at most one retained row per variant.

A completed partial output chunk now prevents an empty/dense-only productive
forecast from authorizing a switch. First-chunk output observations can refute
an inadequate scenario set, but do not establish future selectivity or renew
an immutable component price. Sparse fractions and placements are explicit
conditional assumptions. The three short-window extrapolation limit and the
unresolved loaded absolute-calibration error still prevent a qualified
large-job JIT decision; this change supplies one missing output case for that
continuation.

The final A100 reduced-output and productive-planner regression passed 230
tests in job `20260923-015657-1317408`. The subsequent focused run passed
12 tests in `20260923-020028-1317646`, including a later issue cursor and
the partial-output guard. A broader API run reached eight unrelated failures
because the `examples/toy` fixture is absent from this checkout and remote
mirror; it did not exercise this scenario code.
