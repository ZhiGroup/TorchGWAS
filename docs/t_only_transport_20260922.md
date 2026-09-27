# T-only result transport and calibration identity

`sumstats_fields='t'` now omits beta from the native unreduced host-result path.
This applies to dense output, variant shards, phenotype tiles and host-selected
significant pairs. The statistic calculation still computes beta internally;
the change removes its pinned destination, asynchronous device-to-host copy,
owned NumPy copy and host significant-pair gather. Statistical thresholds,
missingness handling, variant coordinates and result-file schemas are unchanged.

The low-level streaming interface defaults to `return_beta=True`. JAGWAS and
other joint reductions retain their existing result layouts. JAGWAS still needs
the complete phenotype panel on each GPU and partitions variants only. The
experimental device-significance path still transfers selected beta values;
it does not receive the host-transfer savings described here. Packed BED can
omit the returned field, but its transfer implementation is unchanged.

For B variants and K traits, native dense result transfer changes from
`B * (8*K + 5)` to `B * (4*K + 5)` bytes: t remains FP32, status is one byte,
and residual degrees of freedom remain FP32 per variant. Each ring slot loses
one `4*B*K` pinned allocation. Owned results have three arrays instead of four.
Host significant selection gathers one statistic array instead of two, with
24 rather than 28 bytes per selected pair in its owned payload.

The analytical ledgers count these changes in transfer, control calls, pinned
allocation, owned copy, host selection, phenotype tiling and adaptive memory
admission. Beta's GPU intermediate remains counted. Its previous local can
survive evaluation of the next statistics call, so one previous beta storage
also remains in the memory budget even though asynchronous D2H retention is
removed. Memory budgets still require reserves for the stated unresolved terms.

## Measurement reuse

The [immutable cache and refresh policy](cpu_service_refresh_20260922.md) keeps
the original observation date and requires compatible source, runtime and
measurement context. This optimization is an example of why a cached number
cannot be reused solely because its GPU or parameter name matches.

A t-only finish needs its own independently measured three-array service.
`result_finish_prices(..., return_beta=False)` checks the exact native source,
runtime, layout and paired meter controls. The tiny owned baseline is 288 bytes
for 32 rows and one trait; borrowed output has no owned copy. A four-array
416-byte baseline is rejected. Candidate preparation forwards the requested
output choice without relabeling old finish evidence. Missing or incompatible
evidence prevents that timing calculation; memory-only preparation can still
proceed without claiming a timing prediction.

`benchmarks/direct_calculator_result_lifecycle_primitives.py --omit-beta`
collects the new fixed-size finish primitive. This does not certify all other
calculator parameters or establish prediction error for a complete GWAS job.

## Verification

Remote A100 tests cover statistics and degrees of freedom, missing and invariant
variants, owned and borrowed result-ring reuse, two-GPU source coverage,
significant selection including empty output, and old-layout evidence rejection.
The first batch passed 166 tests; a broader regression batch passed 407; the
final t-only batch passed 19. These batch counts overlap. The broader tests also
cover immutable cache lookup, original age, expiry, drift, interrupted refresh,
public calibration, phenotype tiling, variant shards and JAGWAS accounting.

The public-output audit is `benchmarks/direct_t_only_20260922.py`. Its baseline
forces beta transfer while keeping the same t-only writer, and compares the
persisted t statistics, df and selected coordinates with the optimized path.
It records actual payload counts and output-inclusive API wall time. Its first
chunk timestamp means scientific results available to the consumer, not durable
file publication. A single pair per output layout with differing warm states
does not establish a speedup.

On the fixed synthetic PGEN fixture (2,049 samples, 4,097 variants, 512 traits),
all eight public API runs completed. Persisted t statistics and degrees of
freedom were exactly equal in each baseline/optimized pair. Selected variant
and trait coordinates were also exactly equal for significant output.

| Output route | Baseline result D2H bytes | T-only result D2H bytes |
| --- | ---: | ---: |
| One GPU, full panel | 16,801,797 | 8,411,141 |
| Two GPU variant shards | 16,801,797 | 8,411,141 |
| Two GPUs, three phenotype tiles | 16,842,767 | 8,452,111 |
| Host significant pairs, same phenotype tiles | 16,842,767 | 8,452,111 |

Every route saved 8,390,656 result-transfer bytes. Tiled scans also transfer
status and df on each phenotype pass, explaining their larger totals. These
counts exclude genotype transfer and preparation. The finish probe passed its
source/layout/meter checks separately for one and four workers; this is not a
whole-model accuracy qualification.

The retained report is
[`results/t_only_20260922/execution_v2/report.json`](../results/t_only_20260922/execution_v2/report.json).
Raw fixed-size finish observations are in `results/t_only_20260922/finish/`.
Remote jobs were `20260922-090610-978822` (first tests),
`20260922-090737-979081` (finish probe), `20260922-091007-979677`
(407 regression tests), and `20260922-091242-981137` (19 final tests and eight
completed public runs). The first audit attempt stopped on a benchmark-wrapper
argument error before scanning; the corrected audit used fresh output paths.

The eight output-inclusive API times ranged from 0.465 to 2.462 seconds. Their
order and first-use effects are visible in the report and are too large for a
performance conclusion. No runtime model was fitted to these runs.

This change does not complete the public just-in-time tuning policy or qualify
all component prices for selecting new chunks, phenotype tiles or GPU layouts.
