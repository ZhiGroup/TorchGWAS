# One bounded workload ledger at a productive JIT checkpoint

`productive_checkpoint_ledger` now joins three reports at the same source
identity, issue revision and written-output count: necessary work on a
candidate's exact unissued association rectangles, conservative full-chunk
work on original producers without matrix/part completion, and conditional
array payload for issued output without writer completion. The issued source
report and output backlog are bound to the held productive snapshot; the
candidate source report separately proves exact unissued coverage.

The join preserves **ownership**. A candidate may assign future chunks to a
new GPU, but already-issued chunks stay charged to their original GPU. It
adds indexed read bytes, unpacked H2D bytes, mandatory FP32/FP64 matrix-product
FLOPs and array payload only across the disjoint issued/unissued source
ranges. Shared totals are summed once, while GPU work remains separated by
device. The conditional survivor scenario must be the same for pending
indexed output and candidate future output. An issue or output revision
change invalidates the join.

The report deliberately exposes the two scopes separately. Unissued work is
necessary; entire issued chunks are an upper *nominal workload* because some
stages may already have completed. The summed subset workload cannot be used
as a resource-load lower floor or an elapsed completion upper. Issued D2H,
complete statistics kernels, queue state, writer staging/fsync, final
manifest/directory publication and loaded shared capacity still need
independent service and a finite schedule. The public first-chunk controller
does not call this join or change tile/GPU ownership yet.
The report carries the actual stage-floor audit and lists any missing
mode-specific selector or writer floor separately from the larger set of
missing completion-time services. Even a complete set of floors would not
turn them into an elapsed upper bound.

The A100 checkpoint, source and public-output regression passed 30 tests in
`20260923-031627-1337387`, including a candidate that moves only future
dense chunks to another GPU and a full-panel, two-shard JAGWAS occupancy case.
