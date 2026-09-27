# Compact dense-writer work for post-output JIT planning

`compact_binary_output_work` counts the native dense writer's fixed-chunk
beta/t and optional variant-df streams without building chunk or write-event
lists. It derives staging-copy slices from each chunk's payload, block size
and the `gcd` period at which a chunk ends on a block boundary. It counts
borrowed chunks separately, handles the final short chunk, and derives
`sync_file_range` submit/wait/drop counts from total stream bytes. The result
includes exact array payload, staging bytes and calls, minimum write calls,
closing queued bytes, fsync calls and allocated staging bytes. It is constant
work per stream and does not infer service times.

`native_layout_output_floor` can now take complete explicit
`dense_writer_options`. For each fixed unissued dense partition it adds the
compact writer ledger and counts the optional df array in the output payload
floor. This is useful for a whole-source first pass after useful output. The
existing default floor remains a smaller beta/t payload floor when writer
settings are absent. Significant-pair and JAGWAS output still use their own
selected-row intervals; the dense writer option is rejected for those modes.

On A100, 404 source-expansion equivalence and writeback tests passed in job
`20260922-235306-1252652`. The output-floor integration and three-mode
partial-envelope batch passed 407 tests in job `20260922-235418-1252857`.
The final combined source/frontier/compute/output/writer batch passed 656 tests
in job `20260922-235648-1253302`, including auto-block and disabled-writeback
cases.
An ordered metadata-only 8,086,101-marker, 128-trait, 128-marker-chunk probe
used `borrow_chunks=False`, 16 MiB blocks and a df sidecar. Compact and
expanded ledgers both counted 8,312,511,828 payload bytes, 189,519 staging
copy calls and 525 minimum writes. Compact construction took 0.000092 process
CPU seconds and 10,752 KiB peak process RSS; the expanded construction took
0.142853 CPU seconds and 44,544 KiB peak RSS. Jobs
`20260922-235524-1253070` and `20260922-235525-1253098` ran separately under
uncontrolled server load, so these are component observations, not a measured
end-to-end startup improvement.

The detailed finite graph still expands writer events when it needs their
order, and the public cold planner does not yet use this compact ledger to
select a layout. The [compact writer service](jit_dense_writer_service_floor_20260923.md)
now prices its exact counts with the existing independent finite-graph CPU,
page-cache, writeback and fsync fields after useful output. Dirty-page
throttling, short-write retries, queue stalls, manifest/sidecar metadata and
live writer state remain outside that conditional lower floor. Neither the
ledger nor its priced floor is a completion ceiling or permission to switch
tiles or GPUs.
