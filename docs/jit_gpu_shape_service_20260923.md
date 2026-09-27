# Compact GPU shape service after first output

`native_layout_gpu_shape_service` prices the unissued native PGEN scan with
the existing source-faithful `_shape_component` calculator. Each partition
contributes its complete chunks and at most one tail. Distinct
`(device, B, K)` shapes are priced once; a repeated tile or shard multiplies
the same shape service by its exact chunk count. The result records per-GPU
serial kernel service and the shared host-dispatch CPU work separately. It
keeps JAGWAS on the complete phenotype panel and requires exact compiled
statistics/projection geometry. A bounded shape count is checked before any
tensor trace is priced.

This is opt-in through `productive_partial_floor(shape_profiles=...)` after
the useful-output frontier. `productive_checkpoint_ledger` checks the shape
report against the same source identity, reduction, sample count, covariate
rank, chunk size and candidate partition geometry. It does not blend this
conditional service estimate into the necessary resource floors or the
issued-work upper counts. The report retains profile hashes so a caller can
compare its price inputs with the active context; hashing alone
does not establish freshness. The shape report is not used by the public JIT
switch yet.

`productive_issued_gpu_shape_service` applies the same pricing to bounded
full chunks issued under the original layout but lacking output completion.
It keeps the original device and full JAGWAS panel, even when the candidate
moves future work to another GPU. The checkpoint can retain both optional
shape reports. The issued service is conservative replay work: some of it
may already have executed when the checkpoint was taken. Significant-pair
output accepts the scan's host-selection (`None`) or device-selection
(`device_significant`) mode while preserving the output scenario label. The
productive composer checks that the supplied scan mode matches its actual
host or device selector backend before accepting a shape report.

The aggregation allocates O(partitions plus distinct shapes) objects instead
of O(chunks) GPU graph objects. Pricing a new tensor shape can still take
meaningful CPU time, especially for a huge phenotype tile; it belongs in a
charged background planning step or a fresh cross-job structural cache hit,
never a cold first-output callback. The existing first-chunk measurements
may inform loaded-service scenarios, but they cannot overwrite the saved
immutable hardware and geometry evidence.

This fills conditional GPU statistics/projection shape terms for future and
issued chunks when both optional reports are supplied. Transfer,
result preparation, selector and writer service, queue occupancy, final
durability, and loaded shared capacities must be composed into a finite
continuation before ranking chunk/tile/GPU moves. The H100 absolute
calibration gap remains unresolved.

A100 targeted run `20260923-032955-1340677` passed nine shape, source and
checkpoint tests after the issued-service join. A broader regression before
that extension passed 42 tests in `20260923-032722-1340524`. The shape tests
verify tile ownership, full/tail counts, JAGWAS full-panel binding, exact
geometry refusal and source-bound joining; their component prices are
synthetic accounting controls, not calibration.
An additional A100 test (`20260923-033232-1341028`) matched the compact
JAGWAS per-device kernel sums against the existing expanded scan calculator
using real captured launch geometry and synthetic resource prices.
The final focused and calculator regression passed 45 A100 tests in
`20260923-033407-1341205`, including host/device significant scan-mode
binding and the source-faithful JAGWAS comparison.
