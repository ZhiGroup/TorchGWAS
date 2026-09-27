# Experimental fused selected-result packing

The CPU packing prototype remains experimental. It preserved every tested output
but did not provide a consistent whole-GWAS improvement. Production source and
calibration prices are unchanged.

## Proposed work reduction

`benchmarks/selector_pack_control_20260922.cpp` combines selected-value gathers,
flat-index division/remainder, row-df gathering and global variant-index rebasing
in one allocation-free CPU loop. The Python harness retains the existing bounded
predicate, inclusive FP32 cutoff handling and `numpy.flatnonzero`. It allocates
the same owned coordinate and value arrays and calls the loop with the GIL
released. This is CPU C++; it introduces no CUDA or Triton kernel.

The experimental contract is contiguous FP32 t and optional beta, row-broadcast
FP32 df (including scalar and reversed row stride), and an overflow-safe
nonnegative variant offset. The production fallback for arbitrary array layouts
is not replaced. For R retained pairs and b=0/1 for optional beta, the replaced
steps' logical source accesses decrease from (72+16b)R to (40+8b)R bytes. These
are source-level access counts, not measured DRAM traffic. Predicate work,
flatnonzero, output allocations and array destruction remain.

## Generic component controls

A100 job `20260922-143622-1081055` completed 188 exact checks: 20 boundary/layout
cases and 168 timed comparisons over two shapes, three retention patterns,
optional beta, seven randomized paired repeats and two implementations.

| Shape | Retention | Fields | Median paired CPU ratio, fused/current |
| --- | --- | --- | ---: |
| 256 x 4,096 | Empty | t | 0.963 |
| 256 x 4,096 | Empty | beta+t | 0.991 |
| 256 x 4,096 | Sparse | t | 1.014 |
| 256 x 4,096 | Sparse | beta+t | 0.772 |
| 256 x 4,096 | Dense | t | 0.796 |
| 256 x 4,096 | Dense | beta+t | 0.430 |
| 1,024 x 8,193 | Empty | t | 1.186 |
| 1,024 x 8,193 | Empty | beta+t | 0.906 |
| 1,024 x 8,193 | Sparse | t | 1.031 |
| 1,024 x 8,193 | Sparse | beta+t | 0.933 |
| 1,024 x 8,193 | Dense | t | 0.862 |
| 1,024 x 8,193 | Dense | beta+t | 0.948 |

The timed call includes allocation and selection but excludes destruction of
returned arrays. Destruction is recorded separately after correctness checking.
For the large dense arrays, its median CPU cost was 74.1/55.3 ms for current/fused
t-only and 63.8/63.4 ms for beta+t. Allocation history differed between cases;
these values are observations, not transferable allocator prices. Removing
logical passes does not remove the large first-touch or release costs.

## Completed public API comparison

A100 job `20260922-144111-1083456` completed all 24 fresh-process observations:
four cases, three counterbalanced paired repeats and two implementations.
Inputs and output were on verified `/data` XFS storage (`/dev/md0`). Input
hashing and sequential read controls warmed the input; CUDA was warmed before
the API boundary. The read controls include hashing and do not establish an I/O
ceiling. Source, input, binary, execution-context and harness identities are
retained in the report.

The fixed native hard-call PGEN workload used N=2,049, M=16,385, K=512 and C=2;
chunk width 4,096; trait width 256; cuda:1 and cuda:2; four total readers;
prefetch two; queue depth two; and one durable indexed writer. Both arms used
the same native host predicate. Only the post-predicate packing implementation
changed. This is a small control workload, not voxel-scale qualification.

| Case | Output pairs/run | Median paired API ratio | Median paired scan-and-write ratio | Median paired process CPU ratio |
| --- | ---: | ---: | ---: | ---: |
| Empty, t, threshold 1e-30 | 0 | 1.408 | 1.255 | 1.317 |
| Sparse, beta+t, threshold .02 | 173,779 | 1.251 | 1.102 | 1.186 |
| Dense, t, threshold 1 | 8,389,120 | 1.007 | 0.999 | 1.066 |
| Dense, beta+t, threshold 1 | 8,389,120 | 1.255 | 1.237 | 1.080 |

All ratios are fused/current; values below one favor the prototype. Each is the
median of three within-repeat ratios, not a ratio of separate medians. API time
includes input opening/QC, preparation, scanning, durable output and final
metadata. The writer's scan-and-write interval excludes its separately reported
setup. Process CPU sums CPU consumption and is not elapsed time.

All persisted coordinates, t, df and optional beta were bit-identical across
the six runs of each case. The empty case checks zero-length outputs. Whole-job
timings varied substantially; for example, current dense beta+t API runs took
1.58, 1.83 and 3.46 seconds. Even after removing setup, sparse beta+t was slower
in all three prototype comparisons. These observations support declining the
change, not claiming a universal regression or fitting a slowdown coefficient.

The harness's optional `phase_breakdown` field is null because the API calls
that metadata `phase_seconds`. The explicitly captured writer boundaries and
API/process CPU clocks are present. The original report/harness are preserved;
the null field is not used in the comparisons above.

## Evidence and decision

- `results/selector_pack_control_20260922/report.json`: generic observations,
  original timestamps, 188 checks, source/compiler/binary identity.
- `results/selector_pack_gwas_control_20260922/report.json`: 24 API observations,
  exact persisted-field hashes, preregistered order and timing boundaries.
- Corresponding scripts and CPU source are in `benchmarks/`.

All 140 package source hashes and the recorded harness/library hashes were
verified against the pulled local artifacts. No production change is promoted,
and neither component nor GWAS observations are installed as calibration prices.
The source work reduction alone is insufficient evidence for a faster pipeline.
