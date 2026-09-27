# Significant-output gathers and independent timing qualification

The significant-pair selector now gathers contiguous t/beta payloads using
flat indices before converting them into row/column coordinates. A df array
broadcast across traits uses only the row index. General pair/trait df and
noncontiguous payloads retain two-dimensional gathers without copying the
original dense arrays. Selected outputs still own their storage and survive
input-ring reuse.

NumPy remains the default predicate; the optional CPU C++ predicate is
unchanged. Selector identities advance to `bounded_flat_v2` and
`native_row_flat_v2`, rejecting prior gather banks. The source ledger uses
separate `matrix_gather_flat` and `df_gather_row` prices in execution order.
Both count 16 logical bytes per selected FP32 element: int64 index, FP32 read
and FP32 write, compared with 24 bytes for the old two-index gather. Savings
are 16 logical bytes per retained t-only pair, or 24 with beta. This is not a
physical DRAM measurement. Output bytes and conservative dense memory
admission remain unchanged.

## Controlled comparison

A100 job `20260922-123515-1035146` compared the frozen v1 selection sequence
with the proposed sequence on generic buffers containing finite values,
NaN/infinity, sparse survivors and dense survivors. Selected arrays matched
exactly by dtype, shape and byte hash. Each selector case used six
counterbalanced observations on a four-buffer pinned ring. Separate gather
primitives used seven repetitions. No GWAS duration was fitted.

| Shape | Retention | Fields | v1 median CPU ms | v2 median CPU ms |
| --- | --- | --- | ---: | ---: |
| 256 x 4,096 | Sparse | t | 7.925 | 5.536 |
| 256 x 4,096 | Sparse | beta+t | 10.391 | 9.755 |
| 256 x 4,096 | Dense | t | 96.504 | 76.015 |
| 256 x 4,096 | Dense | beta+t | 96.012 | 76.548 |
| 1,024 x 8,193 | Sparse | t | 35.899 | 28.954 |
| 1,024 x 8,193 | Sparse | beta+t | 52.480 | 39.240 |
| 1,024 x 8,193 | Dense | t | 318.227 | 269.343 |
| 1,024 x 8,193 | Dense | beta+t | 966.696 | 817.692 |

These are selector CPU observations, not end-to-end GWAS speedups. Large cases
retain about 3% and 97% of pairs. Input restoration, validation and output
destruction are outside both timing boundaries.

## Numerical and public execution validation

The regression suite passed 339 cases. Two additional cases initially failed
because their ndarray input guard also intercepted the assertion helper's
flattening of selected output. Restricting that guard to the selector boundary
fixed the test, and both cases passed: 341 checked cases across the two runs.
Production code did not change for the test correction. Coverage includes
FP32/FP64 cutoffs, nonfinite values, empty dimensions, block tails, independent
beta/t layouts, general df broadcasting, ownership, source traffic,
legacy-bank rejection, prepared windows and public tuning.

Job `20260922-124313-1036959` completed 16 public GPU scans using the local
`/data/zxie3/torchgwas_adaptive_candidate_fixture_v1_20260922` fixture:
N=2,049, M=4,097, K=512, chunk 128, three reader workers, prefetch three and a
durable indexed writer. Both predicates ran with v1/v2 gathers.

| Case | Layout | Fields | Pairs in each of four runs |
| --- | --- | --- | ---: |
| Sparse, threshold .02 | One GPU | beta+t | 41,742 |
| Sparse, threshold .02 | Two GPUs, trait width 193 | t | 41,742 |
| Dense, threshold 1 | Two GPUs, trait width 193 | beta+t | 2,097,664 |
| Empty, threshold 1e-30 | Two GPUs, trait width 193 | t | 0 |

Persisted coordinates, t, df and optional beta matched exactly across the four
implementations within each case. Timings varied and each condition ran once;
these scans establish output parity, not aggregate throughput gains. Input,
package and harness hashes are retained in the report.

## Price transfer remains unqualified

Fresh primitive banks retain all observations and their original start/finish
timestamps. They are new artifacts, not edits to old measurements or automatic
profile installations. The collector verifies NumPy core identity and CPU
features before and after measurement.

An initial probe mismatch was corrected: `np.nonzero` supplied interleaved
row-index views, whereas production's flatnonzero/divmod indices are
contiguous. The revised collector constructs contiguous indices and records
their strides. Job `20260922-125353-1039223` recollected both banks and a
held-out panel. Transfer checks require the correct index layout and matching
source, threads, affinity, NumPy settings and native binary. No coefficient
is fitted to the held-out observations.

The corrected probe's held-out 1,024 x 8,193 t-only selector results are:

| Predicate | Retention | Predicted CPU ms | Observed median CPU ms |
| --- | --- | ---: | ---: |
| NumPy | Empty | 56.878 | 66.739 |
| NumPy | Sparse | 68.597 | 61.674 |
| NumPy | Dense | 322.210 | 517.008 |
| Native | Empty | 16.654 | 21.684 |
| Native | Sparse | 23.752 | 34.689 |
| Native | Dense | 183.200 | 451.564 |

Large dense service remains substantially underestimated. Earlier unqualified
banks had markedly different rates even for unchanged primitives. Favorable
rows and matching identities do not establish stable capacity. These banks
must not be relabeled as qualified full-pipeline models.

A further locality control, job `20260922-125017-1038650`, collected 192
randomized observations of one fixed 1M-cell predicate. Reused buffers had two
untimed warm calls; rotating views spanned 64 MiB of input and 16 MiB of masks
per allocator, without forced cache eviction. NumPy pageable medians were
4.171/4.016 ms for reused/rotating buffers; pinned medians were 3.164/4.584 ms.
Native medians were 0.809/1.536 ms pageable and 0.842/1.732 ms pinned.
Locality matters but does not explain the whole rate discrepancy. All calls
began/ended on CPU 14, with zero major faults and few minor faults. Input
mappings reported NUMA node 0. This does not establish absence of external
interference, cache eviction or frequency changes.

Next qualification should check stability/drift of independent probes during
useful initial chunks and price allocation/page state in large selected arrays.
It must not fit a whole-job slowdown or assume a single cache explanation.
Frozen H100 ranking evidence is unchanged; full calculator readiness is open.

The subsequent [fused packing control](selector_pack_control_20260922.md)
reduced logical gather/rebase passes but did not establish a whole-job benefit
in 24 completed two-GPU scans. It remains experimental; production packing and
calibration prices are unchanged.

The subsequent [allocation-state audit](selector_allocation_state_20260922.md)
records substantial first-touch and returned-array release costs, including
faults in allocations with no live mmap increase. The calculator now exposes
the corresponding source extents and unpriced transfer explicitly. No page
surcharge is added to allocation-inclusive prices.

## Artifacts

- `results/selector_gather_control_20260922/`
- `results/selector_gather_checks_20260922/`
- `results/selector_gather_checks_v2_20260922/`
- `results/selector_gather_gwas_20260922/`
- `results/selector_gather_numpy_prices_v2_20260922/`
- `results/selector_gather_native_prices_v2_20260922/`
- `results/selector_gather_heldout_v2_20260922/`
- `results/selector_gather_price_transfer_v2_20260922/`
- `results/predicate_locality_control_20260922/`

Earlier banks/panel/transfer artifacts without `_v2` remain historical,
unqualified observations. Source and audit scripts are in `src/torchgwas/`
and `benchmarks/`.
