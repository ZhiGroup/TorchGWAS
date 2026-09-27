# Significant-output pricing: source branches and immutable observations

The host calculator now prices empty, sparse and dense boolean index extraction
separately. NumPy 2.2.6 first counts retained entries, skips extraction when
there are none, uses its sparse search at occupancy <=10%, and otherwise uses a
branchless loop. The previous bank applied the largest nonempty per-cell rate
to both branches. In these independent observations that was the all-true
rate: 5.716 ms per 1M cells, versus 0.688 ms for stride-64 survivors.
Source: [NumPy 2.2.6, PyArray_Nonzero](https://github.com/numpy/numpy/blob/v2.2.6/numpy/_core/src/multiarray/item_selection.c#L2688-L2769).

`numpy_nonzero_work.py` declares the branch and version-bound protocol.
Full-run planning, productive-window planning and public reduction-price
loading all check this protocol. Missing or old pooled prices are rejected.
Ordinary GWAS execution is unchanged; this timing protocol is qualified only
for the inspected NumPy 2.2.6 implementation. Unsupported NumPy versions need
source/protocol qualification before this calculator can price them.

The correction reuses the original independent empty, all-true and stride-64
observations. It does not fit the held-out selector timings. The original
artifact and both observation timestamps remain unchanged. The derived bank
records the original artifact hash and measurement-source hashes separately
from the current calculator-source hashes. It is a diagnostic derivation,
not a newly certified calibration profile or a renewal of empirical freshness.

For the held-out 1,024 x 8,193 t-only selector, CPU time was:

| Output occupancy | Original prediction | Corrected prediction | Observed median |
|---|---:|---:|---:|
| Empty | 20.12 ms | 20.12 ms | 23.56 ms |
| Sparse, 262,620 pairs | 78.52 ms | 38.29 ms | 28.76 ms |
| Dense, 8,141,207 pairs | 524.45 ms | 524.45 ms | 498.78 ms |

The large sparse error falls from 173% high to 33% high. On the 1M-cell sparse
control the corrected estimate is 25% low; the small dense estimate remains
52% high. Branch fidelity improves the model but does not establish a tight
prediction interval or whole-pipeline accuracy. Within-branch pattern,
allocator history, first-touch cost and concurrent resource use still matter.
Memory admission continues to assume dense retained output.

The first transfer audit rejected an actual PyTorch thread-count mismatch:
environment declarations matched, but the held-out process had one thread and
the saved probe had four. The accepted comparison uses an explicit four-thread
rerun with matching NumPy huge-page settings, affinity, selector sources and
native binary. Its original timestamps were not replaced by the rerun date.

Remote job `20260922-120853-1028888` passed 204 focused regression tests, then
completed the derivation and checked that both source observation artifacts
were byte-identical afterward. Tests cover the 10% boundary, invalid counts,
distinct branch services, version/protocol rejection, public and productive
planning, existing binding guards and optional native selection.

Artifacts:

- `results/nonzero_branch_checks_20260922/pytest.txt`
- `results/nonzero_branch_repricing_20260922/report.json`
- `results/nonzero_branch_repricing_20260922/derived_prices.json`
- `results/native_host_predicate_prices_20260922/selection_prices.json`
- `results/native_host_predicate_matched_v2_20260922/report.json`
- `benchmarks/reprice_nonzero_branches_20260922.py`

## Completed large H100 diagnostic

Diagnostic job `20260922-113847-1023606` completed eight fresh-process scans at
the frozen N=35,365, M=1,048,576, K=16,385 layout. The two factors were the
NumPy/native selector and spinning/blocking CUDA completion waits. The native
arm injected the separately hashed production selector into the frozen scan;
the frozen source, inputs, plan and original results were not overwritten.
Each run retained cold-input checks, read-scan-read controls and 140 numerical
reference cells with exact df and maximum absolute t error 1.296739267e-5.
All runs had empty significant output.

| Selector | Wait | Executor seconds, two observations | Median selector CPU seconds |
|---|---|---|---:|
| NumPy | Spin | 77.31, 142.93 | 85.92 |
| NumPy | Block | 102.75, 76.32 | 60.47 |
| Native | Spin | 80.58, 67.66 | 33.45 |
| Native | Block | 93.37, 111.49 | 44.38 |

The native selector reduced the measured selector CPU work, but overall scan
times varied substantially and do not establish an end-to-end speedup. Loaded
worker CPU and CUDA event spans remain diagnostics, not independent capacity
prices. The selector remains opt-in and the event default is unchanged.
The diagnostic project's `results/native_selector_factorial_v1_20260922/`
contains every run plus `summary.json`. Current full-model runtime accuracy
and robust profitable autotuning remain unqualified.

Subsequent [gather work](selector_gather_20260922.md) changes the selector
protocol again, verifies persisted output and records new independent prices.
The branch derivation above remains historical evidence with its original
observation dates; it is not automatically valid for the new implementation.
