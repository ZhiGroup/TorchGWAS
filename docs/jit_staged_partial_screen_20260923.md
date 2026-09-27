# Bounded post-output screen for staged JIT layouts

The experimental `productive_staged_partial_screen` compares up to four explicit chunk/phenotype-tile/GPU layouts against one exact unissued association frontier. It reuses one completed `ProductiveSourceStage` ledger and reports its worker CPU/wall cost once. Each layout has independently supplied source, compute, output, and mode-specific service prices. Source capacity is shared across GPU partitions, and every candidate must cover the held frontier exactly. JAGWAS candidates must keep the entire phenotype panel in each partition and may shard variants.

The screen adds the native dense writer service, host significant selection and archive, device significant count transfer, mandatory selector launches and archive, or JAGWAS selection and archive as appropriate. It includes conditional occupancy for reduced output, so empty early significant chunks cannot silently establish a sparse whole-job regime. The returned envelope contains only necessary resource/serial-service floors. Remaining device significant selection kernels, queued/in-flight work, final durability, live capacity, and a finite completion ceiling are still missing. `prediction_complete` and `selection_validated` remain false; the screen cannot apply a chunk, tile, or GPU change.

The candidate cap is four by default. Source records, partitions, chunks per partition, CPU, and wall also have explicit bounds. The CPU/wall budget is cooperative after a candidate: it cannot interrupt one source rebase. Malformed layout coverage is checked before the expensive schedule walk. The output includes each source-floor calculation cost, each partial-floor cost, the one prior-stage cost, and the complete/stop status.

## Large-source metadata probe

The two A100 probe reports are [final.json](../results/staged_partial_screen_probe_20260923/final.json) and [retry1.json](../results/staged_partial_screen_probe_20260923/retry1.json), produced by [staged_partial_screen_probe.py](../benchmarks/staged_partial_screen_probe.py) from the server-local `/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen` (8,086,101 variants; 20,838,552,600 file bytes; `/data` on `/dev/md0` xfs). The probe simulated eight completed dense output events to trigger staged planning. It performed no genotype scan, GPU computation, or output write. All service prices were synthetic; these are planner-cost measurements, not throughput or a candidate winner. The final report's source-code hashes match the local checkout.

Both runs completed eight bounded metadata steps and two exact-coverage candidates using `staged_primary_rebase`, retaining one prior-stage charge. The first run used 6.496 CPU / 6.567 wall seconds for the stage and 4.126 CPU / 5.311 wall seconds for the screen, with 584,692 KiB peak RSS. The final-code run used 4.298 CPU / 4.318 wall seconds for the stage and 7.292 CPU / 9.332 wall seconds for the screen, with 593,228 KiB peak RSS (579 MiB). The baseline 128-marker/single-GPU source-floor calculation took 0.904 then 1.104 CPU seconds; the 1,024-marker/two-GPU phenotype-tile source-floor calculation took 0.274 then 0.431 CPU seconds. The screen also priced the dense writer service in both layouts. The spread under shared-server load makes a fixed planning-time estimate unsafe. These times do not establish whether eight real first chunks provide enough wall time; the stage and screen must remain off the output callback.

The final focused A100 suite passed 26 tests across staged screening, source floors and the partial envelope. Tests cover dense, host and device significant output, full-panel JAGWAS, two-GPU partitions, global retained-pair invariance, empty device-selection count barriers, staged source reuse, bounded stopping and invalid frontier rejection.
## Next implementation gate

Run the screen in the already opt-in background planning path after real written output, with current immutable component-price bindings and the actual output boundary. Rebase or discard stale source/output revisions. Establish a finite, output-inclusive continuation for both the current and proposed layouts, including active/queued writer work, device significant selection service and final drain. Only then can a bounded optimization compare completion envelopes plus planner/switch cost. A live switch requires fresh memory admission and measured agreement on held-out native PGEN jobs in dense, significant-pair and full-panel JAGWAS modes.
The [compact device-selector launch floor](jit_device_selector_launch_floor_20260923.md) now charges mandatory CUDA count/flagged-selection launches and bounded scatter launches, joined serially with the blocking count transfer per GPU. It still omits other selector kernels and cannot establish completion time.

The [large-phenotype metadata probe](jit_large_k_staged_screen_20260923.md) covers 600,000 phenotypes in six versus twelve two-GPU significant-pair tiles on the 8.1M-variant PGEN. It exposed and repaired the sparse-fraction cap and repeated source-price calculation. A same-stage A/B/B/A control reduced median screen CPU from 5.236 to 1.336 seconds with identical conditional floors; this does not qualify a runtime winner.

The final broader A100 layout/occupancy/controller suite passed 312 tests after the large-K sparse-fraction and source-price-reuse changes.

The [first-chunk live screen hook](jit_live_staged_screen_20260923.md) can now
launch this same evidence-only calculation from the public productive
controller after the staged source worker completes. It uses a bound written
output checkpoint and labels stale issue/output or profile state. Candidate
price binding, finite completion and a live switch remain separate work.
