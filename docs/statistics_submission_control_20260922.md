# Statistics submission control

The source review did not identify a missing main-thread selection dependency
in the scan graph. The independent GPU submission comparison is complete. It
does not support assigning the large full-scan runtime gap to the isolated
statistics API submission step. No coefficient or execution policy has changed.

## Ordering reviewed

In `native_scan.py`, the main iterator fetches its next input, resolves and
yields the oldest pending result when the ring is full, then submits the next
GPU operations. Host significant selection runs in the consumer of that yield.
The source cannot submit into the reused slot until that selection returns.

`execution_graph.py` already makes `host_start:i` depend on `consume:i-depth`
after the initial ring fill. It also keeps the previous host submission as a
dependency and waits for the prefetched input before resolving the old result.
`indexed_schedule.py` places host significant selection before queue admission
and `consume:i`. Input-buffer release remains an independent background path.
This review covers those ordering relationships, not every modeled dependency
or service price.

The frozen decomposition prices 1.240 CPU-seconds of eager API work over 2,048
source chunks. The measured iterator `next()` CPU span includes other work and
instrumentation, so its larger value cannot simply replace the API coefficient.
The new control times the unchanged eager statistics function and input
conversion on generic resident tensors, keeping explicit waits separate.

## Fixed experiment and completed smoke check

`qualify_statistics_submission.py` imports the unchanged frozen H100 package.
The full schedule crosses tiny and large shapes with isolated GPU 2, isolated
GPU 3, and concurrent GPUs 2+3, with three randomized complete repeats. Each
fresh process warms eight calls and measures 128 calls per participating GPU.
A blocking gate every four calls bounds queued computation. CPU inside the
submission wrapper, host wall spans, gate waits and CUDA event spans are
recorded separately. The event spans include launch gaps and are not sums of
kernel service.

The shapes are `(N,B,K,C)=(128,32,32,27)` and
`(35365,1024,8193,27)`. Inputs are deterministic generic device tensors with
declared missing calls. Selected outputs must remain unchanged across repeated
evaluation, with exact expected df. This checks repeatability of the primitive;
it is not an independent OLS validation or a GWAS result comparison.

Smoke job `20260922-160513-1106587` completed the two concurrent-GPU cases,
using four measured calls per GPU after warmup:

| Shape | GPU | Mean submission CPU ms | Mean submission wall ms |
| --- | ---: | ---: | ---: |
| Tiny | 2 | 0.675 | 0.718 |
| Tiny | 3 | 0.671 | 0.721 |
| Large | 2 | 0.678 | 0.725 |
| Large | 3 | 0.671 | 0.710 |

Large-shape peak allocated memory was 2,060,235,776 bytes and reserved memory
2,367,684,608 bytes per GPU, below the 4.5 GiB allocator cap. Both GPU allocator
retry and OOM counters stayed zero. The smoke checks passed, but four measured
calls do not establish concurrency scaling or a capacity price.

Full job `20260922-161156-1107371` completed all 18 planned fresh-process
observations: three repeats of six shape/device conditions. The schedule,
128 samples per worker, stable selected results, exact df, unchanged runtime
identity and zero allocator retry/OOM deltas passed the summary checks.
Startup library loading from shared storage is outside the primitive interval.

| Shape | GPU | Isolated median repeat-mean CPU ms/call | Concurrent median repeat-mean CPU ms/call |
| --- | ---: | ---: | ---: |
| Tiny | 2 | 0.47 | 0.53 |
| Tiny | 3 | 0.44 | 0.53 |
| Large | 2 | 0.60 | 0.64 |
| Large | 3 | 0.64 | 0.63 |

For the large shape, paired concurrent/isolated CPU ratios were
1.28/0.95/1.03 on GPU 2 and 1.19/0.99/0.96 on GPU 3. The concurrent wall
spans were sometimes longer, but include scheduling and waiting. The original
model priced 1.240 CPU-seconds of eager API work over 2,048 calls, about
0.605 ms/call. This generic control's large-shape CPU service is of similar
order and its concurrency effect is inconsistent; it cannot account for a
roughly 40-second full-executor underestimate. That comparison does not prove
the original isolated price transfers inside a loaded decoder/selector/writer
pipeline. CUDA event spans include launch gaps and remain diagnostic only.

## Limits and evidence

GPU 1 was occupied by another job, so the independent control uses idle GPUs
2 and 3. It does not qualify the original GPU1+GPU2 plan. It also excludes
input decoding, result DMA, host selection, writer work and their buffer
lifetimes. The large smoke case uses less memory than the full GWAS pipeline;
the nominal allocator cap alone does not reproduce full-pipeline pressure.

The separate diagnostic project contains `qualify_statistics_submission.py`,
`summarize_statistics_submission.py`,
`results/statistics_submission_smoke_20260922/`, and the complete
`results/statistics_submission_20260922/` including `summary.json`. The executed probe hash is
`aafbfdafefd3eaac052f63299e832aee6e2ddf4d4c1207763fe44fc2fdf7683a`.
All 104 frozen package hashes and both probe/summary script hashes matched
local files after pulling. The summarizer checks the complete schedule,
source/runtime identity, call counts, selected-output stability and allocator
evidence. Maximum allocated/reserved device memory was 2.060/2.368 GB.

This follows the [NUMA selector control](selector_numa_capacity_20260922.md).
The calculator's large-job accuracy and profitable production tuning remain
unqualified.
