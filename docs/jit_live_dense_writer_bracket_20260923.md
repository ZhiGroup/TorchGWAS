# Live dense-writer bracket for JIT accounting

The native dense tile and variant-shard helpers register each active writer
with the productive controller before scanning. The registry holds weak
references, so a completed tile does not stay resident merely to support
planning. A planner captures each active writer twice, with one anchor time
between the passes. It checks the source-issue and completed-output token
before, between and after the passes, and rejects a changed token or writer
identity. With background planning enabled, the two passes run on the
planner worker; each brief controller-lock acquisition can delay a writer
callback. Synchronous planning remains on the callback thread.

For a stream, let `A0,A1` be accepted bytes before and after the anchor, and
`W0,W1` be completely written bytes. These counters never decrease. At the
common anchor, logical pending bytes therefore lie in
`[max(0,A0-W1), A1-W0]`. The intervals can be summed across registered
writers because every first sample precedes the same anchor and every second
sample follows it. This is a cross-writer interval without locking all output
threads. The upper endpoint also bounds bytes still needing `os.write` among
already accepted output. The lower logical endpoint is not a physical-write
lower bound: an active block may already be partly written while its
completed-block counter has not advanced.

The calculator prices that upper endpoint using the existing independent
page-cache CPU and storage coefficients for each writer's device. The public
short-window planner records the capture and priced workload in its attempt
audit. It does not use them to authorize a layout switch. The bracket covers
only writers registered and alive during both passes; it does not include
future output, work inside issued source/GPU chunks, dirty bytes previously
written, writeback calls, fsync, metadata publication or shared-capacity
contention. A stable issue/output token does not freeze GPU work inside an
issued chunk. These remain obligations for a finite completion continuation.

With `source_staging` enabled, each bounded metadata step also takes this
compact live bracket on the staging worker and records its independent writer
price when the capture and profile are valid. A capture failure or retained
observation-budget limit is reported without failing the scientific scan.
The source-stage CPU and wall ledger includes the capture cost. This route
still collects evidence only; it does not invoke the old short-window switch.

The A100 registration/controller/writer regression passed 198 tests in
`20260923-045410-1368746`. The two-pass bracket passed 200 tests in
`20260923-045909-1368968`. After adding the staged-worker observation and
bounded-retention check, the affected source-stage/controller/pricing suite
passed 109 tests in `20260923-050718-1369535`. These establish state and
accounting behavior, not a loaded JIT speedup or absolute model calibration.

The retained-object cap includes captured observations, but it does not prove
a peak-memory bound for temporary NumPy arrays during a metadata step. That
working-set admission remains necessary before a production layout switch.

Reduced-output jobs do not have a dense writer. Their indexed part completion
is already in `ProductiveBoundaryProgress`, but JAGWAS shared-result queue
occupancy and significant-pair selection/part-writing in flight still need
mode-specific observations. JAGWAS continues to require the full phenotype
panel in every variant partition.
