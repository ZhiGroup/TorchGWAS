# Output-inclusive bounded window comparisons

`window_model.prepared_window_runtime` connects the header-window scan model
to the existing dense and indexed output schedules. It prices a bounded set
of prepared windows through payload drain and fsync, including shared CPU,
DRAM, input/output storage and declared transfer-link capacities. It does not
expand the rest of a large GWAS job or use an association timing grid.

Dense output includes beta/t, per-tile df streams, coalescing, writer copies,
page-cache/writeback work and fsync. Significant-pair output uses the existing
host selector and an explicit survivor count for every chunk; zero survivors
still incur dense transfers and selection work. JAGWAS retains its full
phenotype panel and projection work, transfers reduced results and writes
indexed parts. Indexed multi-GPU layouts use the existing bounded global queue
and single writer schedule. Queue acceptance is not counted as durable output.

The production-supported partition rules are enforced: dense output permits
trait or variant partitions, significant output uses trait partitions, and
JAGWAS uses variant partitions with the complete phenotype panel on every
device. Variant layouts require one window per active-device shard. Source
chunks, records, windows, dense output blocks and writeback actions have
explicit preflight limits. A final graph-node limit checks the constructed
graph; it is not a preemptive construction-time bound.

`compare_prepared_windows` checks that two layouts cover exactly the same
variant–phenotype rectangles. It supports different chunk sizes, trait tiles
and device assignments while keeping output mode, threshold and writer
settings common. Coverage uses bounded interval sweeps, without constructing
a variant-by-phenotype matrix. Different tile counts may legitimately require
different numbers of df streams; different chunk counts may require different
numbers of NPZ headers even when retained rows are identical.

Reduced comparisons use one common survivor-count ledger, bound to the input
identity, mode, threshold and total phenotype count. Counts can be combined
across whole evidence bins. A bin can be split exactly only if every cell or
no cell survived; otherwise finer evidence is required. This prevents a
configuration comparison from silently assuming uniform selectivity. A live
caller must additionally bind that ledger to its phenotype/covariate inputs
and scientific settings. This arithmetic API does not establish those
content identities or forecast occupancy in unobserved variants or tiles.

## Verification and measured calculation cost

Remote job `20260922-045456-913782` passed **153 tests in 23.05 seconds**.
The suite includes window composition, the prior header/scan model, existing
dense and reduced calculators, and immutable calibration-cache validation.
It checks writer/selector service equality against the exact-census models,
shared-resource work conservation, coverage and partition restrictions,
survivor rebinning, changed evidence, and rejection before costly expansion.
The log is `results/prepared_window_model_v3_20260922/tests.log`.

The arithmetic audit in
`results/prepared_window_arithmetic_v2_20260922/report.json` evaluates 15
two-GPU scenarios at chunk sizes 128, 256 and 512: dense output, significant
output with empty/sparse/dense survivors, and JAGWAS. It then evaluates five
equivalent-work comparisons twice, starting with a fresh calculation cache
for each mode. The source fixture has 2,049 samples and 4,097 variants;
each device processes a 512-variant window. Dense/significant scenarios use
two 512-trait partitions; JAGWAS uses the complete 512-trait panel on each
variant shard. GPU launch geometry is captured, but all service prices are
synthetic controls. No association run is timed by this audit.

Actual `numpy.savez` file sizes agree with every nonempty modeled indexed
payload, and empty results produce no part. Dense beta/t/df payload is checked
analytically. Genotype and output artifacts are on `/data`, XFS `/dev/md0`.
All 129 package files, five helpers/geometry files and the benchmark hash were
verified against the delivered source. The audit is evidence for accounting
and calculation overhead, not prediction accuracy or measured GWAS speedup.

| Output | Fresh comparison CPU | Repeated comparison CPU |
|---|---:|---:|
| Dense | 58.80 ms | 17.93 ms |
| Significant, empty | 58.69 ms | 23.40 ms |
| Significant, sparse | 75.14 ms | 24.01 ms |
| Significant, dense | 58.98 ms | 20.78 ms |
| JAGWAS | 69.54 ms | 23.71 ms |

Comparison timing includes validation, survivor rebinning, graph construction
and solution. It excludes caller header construction and parameter binding.
Each retained comparison cache holds two shape components, approximately
218 KB for dense/significant or 304 KB for JAGWAS. These are in-process
calculation caches, not empirical calibration records. Their reuse neither
validates nor renews any measured service price.

## Integration boundary

These are isolated prepared windows with empty initial queues/writer staging
and a final payload drain. They omit phenotype/factor setup, already-issued
work, existing writer state, API metadata and directory publication. Decoder
endpoint simulations are scenarios, not proven elapsed-time bounds. The
finite-window difference must not be multiplied by remaining work and treated
as a certified live-switch saving.

Fresh comparisons still exceed the default 50 ms total planning CPU budget.
The existing productive-run budget would reject an over-budget result. A
public automatic policy still needs dependency-valid structural reuse or
explicitly charged incremental construction, an observed continuation state,
an accountable remaining-work forecast, live memory admission and measured
prediction qualification. There is no automatic configuration change in this
patch.

Saved calibration follows the existing separate rules: structural records
require matching dependencies; empirical records retain their original
observation time and expiry; memory availability and contention are read
again live. Early-chunk validation is separate evidence and cannot make an
old observation newer. A cached loaded GPU interval remains a loaded interval,
not an independently measured hardware capacity.

The follow-up `structural_tensor_reuse_20260922.md` implements immutable
cross-job reuse of the duration-free tensor ledgers and connects it to the
productive controller. Its final tests and timing limitations are documented
there; reuse does not by itself establish that the planning budget is met.
