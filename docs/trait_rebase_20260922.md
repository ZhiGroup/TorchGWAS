# Explicit single-copy trait rebasing and calculator correction

The significant-pair producer now copies block-local trait indices into one
owned int64 array and adds the tile offset in place. Both host and device
selection use this common producer. The calculator now prices `index_cast`
followed by `inplace_index_add`, and its allocation ledger counts one rebased
array. Inputs may be borrowed, strided or read-only; they are never modified.

## Allocation evidence corrected the initial expectation

The earlier expression was `trait_index.astype(np.int64) + offset`. Counting
those two Python operations as two material arrays was too conservative for
large outputs on the tested installation. NumPy 2.2.6
[tries temporary reuse in array addition](https://github.com/numpy/numpy/blob/v2.2.6/numpy/_core/src/multiarray/number.c#L215-L225).
Its [elision implementation](https://github.com/numpy/numpy/blob/v2.2.6/numpy/_core/src/multiarray/temp_elide.c#L261-L328)
requires an eligible unique, writable, owned numeric temporary and compatible
callers/casting. Non-debug builds use a 256-KiB eligibility threshold.

The actual source iterator was compared against a frozen pre-change copy,
using the forwarding allocator control. Four counterbalanced observations
per row count produced identical allocation totals within each condition:

| Retained indices in the nonempty tile | Old requested data bytes | New requested data bytes | Difference |
| ---: | ---: | ---: | ---: |
| 0 | 36 | 34 | 2 |
| 7 | 146 | 89 | 57 |
| 1,024 | 16,418 | 8,225 | 8,193 |
| 32,767 | 524,306 | 262,169 | 262,137 |
| 32,768 | 262,186 | 262,177 | 9 |
| 32,769 | 262,194 | 262,185 | 9 |
| 1,048,576 | 8,388,650 | 8,388,641 | 9 |
| 8,389,632 | 67,117,098 | 67,117,089 | 9 |

Each control also emits an empty first tile. Totals include NumPy data buffers
for that tile and small scalar temporaries; they do not include Python object
headers, allocator rounding or RSS. No durations are retained in this census.
The observed crossover agrees with NumPy's 256-KiB threshold for int64 arrays.

Thus the change removes a material temporary for small outputs and makes
single-copy behavior explicit for all sizes. It does **not** save another
64-MiB array for the tested 8,389,632-index case: that temporary was already
elided. The earlier expectation of that saving was withdrawn after tracing.

## Calculator consequences

An independent `rows + offset` probe retains a reference to `rows`, whereas
the old compound production expression could consume its temporary in place.
Identical dtype and shape therefore did not guarantee identical allocation
behavior. The explicit in-place production operation now matches the separate
in-place primitive, without depending on interpreter/build/size elision.

For retained count R, the outer rebase has one 8R-byte data allocation.
The two logical read/write passes still contribute 32R source bytes; an
allocation saving is not presented as elimination of a numeric pass. The
conservative host/device producer memory ledgers remove one worst-case int64
temporary per current/previous iteration. The tighter memory estimate is an
accounting improvement, not a measured large-array RSS reduction.

Device-selection price banks must provide `inplace_index_add`; an old
allocating-add price is rejected rather than relabeled. Existing measurements
remain immutable. Allocating-add prices can still describe operations that
actually allocate, including other unchanged reduction paths. Remaining
allocator-state and full-pipeline qualification gaps are not resolved by this
change.

## Correctness and timing boundaries

Job `20260922-135731-1070602` passed 330 tests, including independent source
pricing, host/device selection, prepared windows, multi-GPU queues and new
borrowed-index ownership cases. Four existing warnings describe tests that
request three readers with prefetch depth two. Job
`20260922-140345-1074020` passed another 72 memory/adaptive-tuning tests and
completed the allocation census. Both jobs are terminal.

The first job also completed 24 paired public GPU scans: old/current rebasing
for NumPy-host, native-host and device selection in each of four cases.
The local `/data` fixture has N=2,049, M=4,097 and K=512, with chunk 128,
trait width 193, three readers, prefetch three, queue depth two and durable
indexed output. Coordinates and every persisted statistic field matched
exactly within all twelve before/after comparisons:

| Case | Devices | Fields | Retained pairs per run |
| --- | --- | --- | ---: |
| Tiled sparse, threshold .02 | cuda:1 | beta+t | 41,742 |
| Tiled sparse, threshold .02 | cuda:1, cuda:2 | t | 41,742 |
| Tiled dense, threshold 1 | cuda:1, cuda:2 | beta+t | 2,097,664 |
| Tiled empty, threshold 1e-30 | cuda:1, cuda:2 | t | 0 |

Separate seven-repeat CPU controls did not establish a speedup. Default-
allocator median rebasing CPU times were old/new 4.789/5.875 ms at 1,048,576
indices and 71.469/73.336 ms at 8,389,632 indices. Shared-host variation
remains substantial; public scans had one observation per condition. No
end-to-end throughput gain is claimed.

## Artifacts

- `results/trait_rebase_checks_20260922/pytest.txt`
- `results/trait_rebase_memory_checks_20260922/pytest.txt`
- `results/trait_rebase_control_20260922/report.json`
- `results/trait_rebase_allocation_census_20260922/report.json`
- `results/trait_rebase_gwas_20260922/report.json`
- `benchmarks/legacy_trait_rebase_20260922.py`

Persisted outputs are under `/data/zxie3/torchgwas_trait_rebase_20260922` on
lab-a100. The frozen pre-change function hash is
`e46cf5c54f97beb2c8420c6512cd2fa909eae475cbee1ca1d7e2bbe3912e296b`.
All 140 current package source hashes and all three public-run harness hashes
were verified against the pulled report. Old measurement artifacts and the
frozen H100 project were not changed. Full calculator readiness remains open.
