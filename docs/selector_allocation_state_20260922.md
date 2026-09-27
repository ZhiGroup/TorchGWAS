# Selector allocation state and first-touch evidence

The calculator now records source allocation extents for the host significant
selector and explicitly exposes unpriced allocation transfer in both full and
prepared-window reports. No page surcharge or fitted slowdown is installed.
Prediction and selection qualification remain false.

## Source accounting

For a chunk with `B*K` cells and `R` retained pairs, FP32 row-df contiguous
selection returns two int64 coordinate arrays, FP32 t and df, and optional
FP32 beta: `(24 + 4*return_beta)*R` data bytes. The mask adds `B*K` bytes.
The current outer trait rebasing makes one owned int64 copy, another `8*R`
bytes, and adds its offset in place. Variant-index rebasing is also in place.
The initial note counted `16*R` from the old compound cast/add expression;
subsequent allocator tracing showed that NumPy elided its large temporary.
That earlier count was conservative source accounting, not observed
allocation. The [rebase correction and boundary census](trait_rebase_20260922.md)
document both the small-output saving and the absence of a large-array saving
on this installation.

The new ledger identifies each extent, allocating primitive and logical
lifetime. It counts nonzero's underlying storage even if NumPy returns a view.
Zero-length output has zero data bytes but still has fixed dispatch/metadata
cost. Cutoff and bounded predicate scratch remain in their primitive prices;
this is not a complete allocator-event trace or a peak-memory replacement.
Extents alone establish neither fresh pages nor the thread that eventually
releases an array. Conservative memory admission is unchanged.

## Independent observations

Job `20260922-131928-1049986` completed on lab-a100 with affinity 12-19,
PyTorch four threads, NumPy huge-page advice disabled, and passive OpenMP
waiting. It first collected new NumPy/native primitive banks with per-call
thread faults, user/system CPU counters, net mallinfo2 state and separately
timed destruction. It then ran seven paired first/warm copies of one fixed
64-MiB generic NumPy allocation, followed by generic held-out selectors on a
four-buffer pinned input ring. No genotype or GWAS duration entered pricing.

The copy control observed 16,384 fresh-page faults per repetition (16,386 in
the first), and zero faults on every warm rewrite. The median signed
first-minus-warm CPU increment was 2.906 microseconds/page, or 47.611 ms per
64 MiB. Raw signed observations were neither clipped nor selected. First
necessarily precedes warm, so the difference can also contain cache/order
effects; this single-size observation does not certify a universal page price.

Held-out t-only observations, medians of seven repetitions:

| Shape | Retention | Predicate | Priced CPU ms | Selector CPU ms | Minor faults |
| --- | --- | --- | ---: | ---: | ---: |
| 256 x 4096 | stride 31 | NumPy | 6.008 | 8.926 | 0 |
| 256 x 4096 | all pairs | NumPy | 47.732 | 46.880 | 4,064 |
| 1024 x 8193 | stride 31 | NumPy | 47.717 | 79.837 | 447 |
| 1024 x 8193 | all pairs | NumPy | 381.541 | 480.781 | 51,214 |
| 256 x 4096 | stride 31 | Native | 3.072 | 4.682 | 0 |
| 256 x 4096 | all pairs | Native | 40.169 | 43.692 | 4,064 |
| 1024 x 8193 | stride 31 | Native | 23.708 | 57.594 | 3,603 |
| 1024 x 8193 | all pairs | Native | 320.518 | 417.582 | 51,178 |

This timing boundary starts with supplied critical values and ends when
`select_host_pairs` returns. It includes destruction of internal temporaries,
but excludes outer trait-index rebasing and destruction of returned arrays.
The predictions use the same boundary. NumPy/native outputs were checked by
shape, dtype and byte hash outside the timed calls.

Every large dense call had a net four live mappings after selection, totaling
192.04 MiB, which disappeared when its result was released. Returned data
alone covers 49,158 pages; the separate all-fresh CPU scale from the fixed
control is 142.850 ms. This is **not added to the prediction**: existing
primitive prices already contain allocation/page work, the gather price uses
the maximum declared pattern, and reference/candidate residency differs.
Adding it blindly would double count and would over-correct this particular
residual. Large sparse predictions still miss substantial service as well.

Large dense returned-array destruction separately cost 52.82 ms (NumPy) and
46.32 ms (native), compared with about 4.5 ms for small dense output. Small
dense selections had no net live mmap increase yet still faulted on about
4,064 pages per call. Thus an arena route does not establish warm residency.
Net mallinfo2 counters do not identify every transient allocation, and these
observations do not isolate the allocator mechanism responsible for arena
faults. User/system counters are diagnostic and may have coarse resolution.

## Reuse and qualification

All three new observation artifacts preserve their original start/finish
timestamps, package identity, harness identity and runtime context. They are
new records; previous banks are untouched. The read-only summary verifies the
three input hashes before and after analysis and publishes zero price records.
Structural identity is necessary but not enough to reuse an empirical service
as a capacity estimate. Allocator history and residency require an explicit
state assumption or a fresh check; they must not become immutable facts merely
because one job observed them.

Before numerical integration, resident kernel work, first touch and destruction
need compatible independent measurement boundaries. The reference contribution
must be removed before adding a candidate page term. Release placement also
needs the actual producer/consumer lifetime. Current evidence supports that
investigation, not a whole-job correction factor or automatic qualification.

## Validation and artifacts

Job `20260922-132353-1052464` completed 225 tests in 33.29 seconds, covering
significant model graphs, prepared windows, source branching and production
selector correctness. New tests compare ledger extents with actual returned
arrays for empty, sparse and dense output, with and without beta, and check
that both planning paths disclose unpriced residency transfer.

The summary records both the observation source and the subsequent analysis
source. Only `significant_host_work.py`, `significant_host_model.py` and
`window_model.py` differ: allocation metadata/disclosure was added after data
collection; selector execution and pricing equations did not change. All 140
current package hashes and all three observation artifact hashes were verified
against the pulled summary.

- `results/selector_page_numpy_prices_20260922/selection_prices.json`
- `results/selector_page_native_prices_20260922/selection_prices.json`
- `results/selector_first_touch_control_20260922/report.json`
- `results/selector_first_touch_summary_20260922/report.json`
- `results/selector_allocation_checks_20260922/pytest.txt`

Frozen H100 evidence is unchanged. These CPU controls are not an end-to-end
GWAS speedup claim, and full calculator readiness remains open.

The follow-up [resident and default-address prefault controls](selector_resident_controls_20260922.md)
separate allocator callback service from unchanged NumPy numerical kernels.
They confirm substantial page-touch cost but fail stable-price/held-out
qualification; no resulting coefficients are installed automatically.
