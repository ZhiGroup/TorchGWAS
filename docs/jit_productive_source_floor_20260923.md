# Whole-source floor at the productive JIT frontier

`productive_source_floor` joins the existing indexed PGEN whole-schedule
reader model to the exact unissued association frontier after written output.
It accepts explicit chunk, tile and device layouts, checks the rectangle
coverage and JAGWAS full-panel rule before any whole-source header walk, and
caps unique records, chunks and partitions. Repeated trait tiles share one
header schedule calculation but each incurs its own physical source work.
The result includes the source floor, coverage proof and charged CPU/wall
calculation time.

`productive_partial_floor` now composes that source result with the
existing H2D/GEMM and native-output payload calculators at the same
frontier. For significant-pair and JAGWAS output it supplies the declared
global sparse counts to the output calculator; dense output has no occupancy
assumption. Optional independently priced stage-service floors can be
included through the same existing envelope. The result is a necessary
whole-job resource floor, with source, compute, output and survivor ledgers
kept separately for audit.

This removes one obstacle to evaluating large jobs without extrapolating three
short windows: the read/decode ledger can cover the whole unissued source
using header metadata and bounded work. It is intentionally not connected to
the public switching bridge yet. The composed floor includes mandatory
GPU matrix work and transfers, but misses some statistics/reduction kernels,
selection and writer service unless independently supplied, plus live
in-flight queues and final drain. Its upper endpoint is an upper bound on the
*necessary partial floor*, not a completion ceiling. It cannot authorize a
chunk, tile or GPU move.

The A100 focused regression passed 140 tests in job
`20260923-020731-1321659`. It includes exact post-output coverage, physical
trait-tile rereads with one metadata pass, full-panel JAGWAS variant shards,
and rejection of wrong coverage or excessive source records before the
whole-source walk. The fixtures use synthetic source-service prices, so this
is accounting and safety evidence, not throughput calibration.
The source/compute/output composition passed 134 focused A100 tests in
`20260923-022204-1326005`; the final stale-frontier guard and sparse
accounting tests passed 27 tests in `20260923-022325-1326232`.

The next implementation step is a finite completion continuation
at the same frontier. It must price remaining stage service and shared capacities once,
include in-flight buffers and writer progress, and provide a qualified
baseline completion lower and candidate completion upper under each
declared source/output/capacity scenario. The measured H100 absolute
prediction gap remains a separate calibration blocker.
