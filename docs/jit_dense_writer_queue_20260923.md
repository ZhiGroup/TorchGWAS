# Dense writer backlog observations for JIT continuation

The native dense writer now attaches a compact queue observation to each
completed common beta/t row-prefix event. Each active array stream reports
bytes accepted into staging or the pass-through queue, bytes completely
written, and the current staged, queued and active bytes. The exact
stream-local invariant is `accepted - written = staged + queued + active`.
The active term counts its entire block until the final successful write,
including when `os.write` has made partial progress. Thus staged plus queued
bytes are still awaiting writes at that observation, while accepted minus
fully written bytes is a conservative upper amount. The same observation
exposes these endpoints as `pending_write_bytes_interval`, together with
the stream's block size, queue depth and queued block count.

A snapshot holds the short state locks for all streams of one writer at
the same time, including during a completed-prefix callback. It is atomic
across that writer's beta, t and df counters, but not across the scan
producer, GPU work, other writer directories, filesystem writeback or the
final fsync. A queue snapshot is an
observation of already accepted output, not a prediction of all issued or
unissued output and not a synchronized execution checkpoint. A later finite
continuation must join it with the held issue/output revision, price any
remaining writer service, and include the final durable drain.

`ProductiveBoundaryProgress` retains only the latest observation per writer
directory within its bounded early-event window. An event's partition id
identifies the completed row prefix; a shared writer may also hold data from
another partition. Native dense events have no producer-device label, so the
boundary now infers ownership only when one admitted variant/trait partition
contains the entire completed interval. A prefix crossing variant shards or
matching more than one owner remains invalid for JIT planning. Explicitly
mismatched device labels remain invalid. This lets the existing single-GPU
dense path bind real writer progress without inventing an owner for ambiguous
multi-GPU output.

The queue counters are installed only when write-progress observation is
enabled. They do not alter output content, borrowing, staging, block size,
writeback, fsync or manifest policy. The calculator does not yet turn these
observations into a completion ceiling or a layout switch. Validation must
include loaded dense output and final drain, not just a zero-backlog unit
test.

The A100 writer/boundary/controller/checkpoint regression passed 172 tests in
`20260923-043125-1363724`; the final interval and multiwriter follow-up
passed 39 focused tests in `20260923-043356-1364029`. These checks include
short writes, staging and pass-through paths, t-only and df-sidecar stores,
and actual native prefix callbacks. They establish byte-accounting and
ownership behavior, not loaded writer throughput or a JIT speedup.

The companion `price_dense_writer_queue_observation` accepts one atomic
writer-local snapshot and the existing independent dense-writer service
profile. It prices the lower and upper accepted-byte endpoints with page-cache
CPU and storage-transfer coefficients; it also reports each stream's nominal
serial page-cache work at the profile's writer CPU fraction. This is a priced
workload interval under fixed coefficients, not a timing confidence interval.
The lower endpoint omits an active block because an `os.write` may have partly
completed it. The report deliberately does not add writeback-call, final
fsync or already-dirty bytes, and it cannot be summed with a producer boundary
from a different instant to assert a completion bound. Joining it to a held
source/GPU/output revision remains necessary for a layout decision.

The atomic snapshot regression passed 172 A100 tests in
`20260923-043944-1367592`. The new priced accepted-work checks and affected
writer/boundary/service tests passed 44 tests in
`20260923-044259-1367851`. These are correctness checks, not a calibrated
loaded writer-throughput measurement.

A [live two-pass bracket](jit_live_dense_writer_bracket_20260923.md) now
combines the stream counters across active writers at one common anchor,
checking source-issue and output-event revisions around the capture. It
prices an upper accepted-write workload with the existing device profiles;
this still omits unaccepted output and the in-flight source/GPU and durable
output tails.
