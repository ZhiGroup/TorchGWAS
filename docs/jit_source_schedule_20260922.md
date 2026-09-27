# Whole-source PGEN schedule work for JIT planning

The cold productive planner used three short source windows. On an 8,086,101-variant job, the actual remaining extent exceeds its conservative eightfold extrapolation cap, so a post-output callback cannot consider a chunk switch. Raising the cap would convert a small-window slope into an unsupported whole-job estimate. This work adds a bounded, header-only description of the **actual remaining chunk schedule** for the native hard-call PGEN reader. It is a source-work component for a later finite continuation model, not a switch decision.

`PgenHeaderWork.schedule_bounds(start, stop, chunk_markers)` aggregates the indexed primary records once and includes the LD base decoded at each independent chunk start. It reads no genotype payload. Explicit record, signature and chunk budgets reject unbounded requests. It retains the source identity, exact indexed read/decode byte counts, and conditional intervals for decoder primitive counts. The intervals assume valid supported record payloads; they are not an exact payload census. `decoder_work` rejects this schedule kind as an exact census.

`paired_schedule_source_difference(baseline, candidate, prices)` cancels the shared primary records **and shared LD restart positions** before comparing schedules. It validates source identity, chunk geometry, replay entries and conservation. Its decoder CPU interval requires a price for every possible primitive. If an unpriced primitive has a nonpositive count difference throughout, the function can still give a one-sided upper bound on candidate-minus-baseline decoder CPU. This is a source CPU-work bound, not a completion-time bound: different chunk sizes also change GPU submission, transfer, result handling and output timing.

For a later finite graph that needs individual chunk bounds,
`PgenHeaderWork.vectorized_window` groups all record signatures once, then
distributes their counts to chunks and adds each chunk's LD base replay. It
returns the same typed per-chunk work as `window`; the ordinary small-window
path remains in use for current productive horizon comparisons. Budgets cap
records, chunks, distinct signatures and signature–chunk pairs. This path is
for a charged background calculation after useful output, not cold startup.
It still creates many chunk objects and is not itself a whole-pipeline graph
or a tuning decision.

## Remote validation

On lab-a100, 165 focused tests passed after the paired-source implementation (`20260922-213936-1214706`). A later chunk-geometry validation patch passed its 10 targeted tests (`20260922-214043-1214866`). Synthetic PGEN cases compare all schedule intervals against exact per-chunk payload censuses, including LD restarts and paired differences; tests also cover budgets, stale input and tampered replay entries.

The vectorized path passed 61 targeted tests (`20260922-220516-1228072`),
including 60 combinations of sample count, source range and chunk size that
matched every field of independent `window` chunk bounds. A broader A100
batch passed 369 tests (`20260922-220934-1229196`), including exact graph
equality for dense and significant-output windows. After the narrower-index
change, 63 targeted source/output-graph tests passed
(`20260922-221433-1229970`). The broader batch predates that final memory
change; it must not be counted as a full rerun of the final source hash.

A read-only lab-h100 diagnostic of the frozen 35,365-sample, 1,048,576-variant source found 1,024 chunks and 230 LD base replays at chunk size 1,024. The header schedule took 0.580 s wall/CPU, versus 18.870 s wall and 18.866 s CPU for the exact payload census. All source-unit intervals enclosed the exact counts; indexed read bytes (4,230,062,817) and decoder input bytes (4,229,874,414) matched. The frozen price profile omits possible 4/5-byte ULEB rates, so the conditional decoder CPU interval is unavailable even though the payload census found no such integers. Artifact: `torchGWAS-calculator-diagnostics-h100/results/full_schedule_bounds_20260922/report.json`. That probe preceded the later paired-difference changes.

On lab-a100, a fresh run on the local `/data` PGEN with 8,086,101 variants and a 20,838,552,600-byte source built these full schedules after a 0.241 s header parse. The timing is an ordered diagnostic on a shared server, not an estimated saving:

| Chunk markers | Chunks | LD replays | Schedule wall | Indexed read bytes |
| ---: | ---: | ---: | ---: | ---: |
| 128 | 63,173 | 15,475 | 2.219 s | 20,884,425,671 |
| 1,024 | 7,897 | 1,903 | 0.181 s | 20,826,445,834 |
| 4,096 | 1,975 | 423 | 0.013 s | 20,820,180,640 |

For 128→1,024 the paired calculation took 0.764 s, removed 13,572 LD replays and 57,979,837 indexed read bytes. For 1,024→4,096 it took 0.175 s, removed 1,480 replays and 6,265,194 read bytes. Every source-unit difference interval was nonpositive, so with no primitive prices the one-sided decoder CPU-work increase bound is zero; this does **not** estimate the saved seconds. The report and exact code hashes are in `results/large_schedule_bounds_v3_20260922/report.json` (job `20260922-214205-1214955`). Earlier ordered v1/v2 timings should not be combined into a speedup comparison.

A separate read-only lab-h100 comparison used the frozen two-GPU decoder primitive prices on the frozen 1.05M-variant PGEN (`20260922-214807-1219426`). For 256→1,024 markers it removed 646 LD replays and 4,401,634 indexed read bytes; the paired decoder CPU-work increase has a conditional upper bound of **−0.000697 s**. For 1,024→4,096 markers it removed 174 replays and 1,239,709 bytes; the corresponding upper bound is **−0.000179 s**. Both bounds retain missing ULEB4/5 prices; those unpriced units can only increase the saving for these aligned coarser schedules. They exclude read service, per-chunk control, selector, GPU, output and live contention, and 4,096 markers have not been admitted as a feasible frozen executor layout. The header calculations used no payload read, preserved the frozen plan/profile/source hashes, and are recorded in `torchGWAS-calculator-diagnostics-h100/results/paired_schedule_source_20260922/report.json`. The small priced LD decoder savings do not explain the frozen 39.907 s prediction gap.

### Cost of individual chunk bounds

On the same 8.1M-variant A100 PGEN, ordinary `window` expansion of all 7,897
1,024-marker chunks took 290.368 s wall / 289.754 s CPU with the default
1,024-entry signature cache (`20260922-215504-1227035`). The source had 12,950
distinct primary signatures; this run recorded 3,094,850 signature-cache
misses. Its next 4,096-marker arm stopped at the explicit 2,048-signature
per-chunk budget, so this job did not publish a final JSON report. The logged
1,024-marker counts matched the compact schedule; do not interpret the
incomplete job as a successful two-arm comparison.

A separate complete 1,024-marker expansion with a 16,384-entry cache took
77.074 s wall / 76.917 s CPU, with 13,904 misses and 210,328 KiB peak process
RSS (`results/large_window_cache_20260923/report.json`, job
`20260922-220019-1227636`). The vectorized path took **7.496 s wall / 7.483 s
CPU**, with 571,616 KiB peak RSS. It exactly conserved all 7,897 chunks,
1,903 LD replays, indexed read/decode bytes and every source-unit interval
against the compact aggregate (`results/large_vectorized_window_20260923/report.json`,
job `20260922-220719-1228627`). These ordered diagnostics show a CPU-time
improvement at a larger peak-memory cost; they are not controlled executor
speedups. Compact `schedule_bounds` remains the cheaper choice when a planner
only needs aggregate source work. Individual bounds should be materialized
only when a finite graph needs them and the host scratch budget admits them.

A follow-up implementation narrows the per-variant chunk/signature index to
`uint32` when its admitted range fits and releases scratch arrays after their
last use. Its real-panel run conserved the same source work with 516,776 KiB
peak RSS (`results/large_vectorized_window_v2_20260923/report.json`, job
`20260922-221526-1230238`). It took 31.790 CPU seconds while the compact
control in that same job rose to 4.345 CPU seconds, versus 0.836 seconds in
the previous job. The shared server and remote filesystem were disturbed; no
algorithmic timing direction is inferred from these unmatched runs. The peak
RSS reduction is directly observed, but a larger phenotype job still needs
explicit host scratch admission before materializing the full chunk list.

## Continuation needed for a switch

The intended bounded optimization has a finite admitted set
`X = {(chunk markers, phenotype tile, device assignment)}` and declared
resource/output scenarios `S`. For each `x` and `s`, construct a finite
continuation graph from the exact unissued variant frontier and the current
in-flight/writer state. Its resource demand includes the header source interval
`U_p(x)`, exact indexed read bytes `R_p(x)`, per-chunk GPU and transfer work,
mode-specific result selection and durable output. A necessary throughput floor
is `max_r W_r(x,s)/C_r(s)` over shared CPU, DRAM, input/output and device
resources, together with dependency critical paths. A feasible finite schedule
or a separately qualified conditional continuation model must supply an upper
completion estimate; the resource floor alone cannot rank candidates. The
decision compares the current lower and candidate upper in each **same**
scenario and requires the smallest gain to exceed measured planning, probing,
switch and publication costs. Memory and full-panel JAGWAS rules constrain
`X` before comparison. During one running job, the present execution bridge
can only vary chunk markers within its fixed admitted tile/device assignment.

The next planner must combine source work with independently priced per-chunk CPU handoffs, host–device transfer, GPU statistics, reduction and dense/significant/JAGWAS output. A full finite schedule should account for shared CPU, DRAM, input/output storage and device capacity, unequal partition tails, current in-flight work and final drain. It must compare bounded baseline and candidate completion under the **same** occupancy and capacity scenario, then repay measured planning and switch costs. The header schedule can replace the small-window source extrapolation without running an exact census before useful output. It should be evaluated on a background worker after the first written chunk, and immutable cross-job records must bind source bytes, decoder implementation, price protocol and original observation age separately from live resource availability.

The frozen H100 calculator's 15.862 s prediction versus 55.769 s observed median remains unresolved. Source schedule correctness does not qualify that whole-pipeline model or justify increasing the extrapolation cap. Full-panel JAGWAS retains variant sharding only; phenotype tiling is for modes that permit phenotype partitioning.
