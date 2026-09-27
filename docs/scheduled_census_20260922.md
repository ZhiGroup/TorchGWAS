# Exact source schedules for adaptive chunk sizes

The calculator previously inferred every boundary from the fixed chunk capacity.
After a 128-to-256-to-512 change, this can omit native PGEN LD restarts, count too
few active reader initializations and charge the wrong number of output parts.
It now accepts an explicit contiguous schedule while retaining the allocated
capacity in the execution profile.

## Accounting contract

`pgen_work_census.scheduled_census(encoded, capacity, chunk_ranges)` regroups a
regular fine census. Every start and interior end must align to that census;
the final file tail may be shorter. Ranges must be contiguous, nonempty and no
larger than capacity. The input and returned records do not alias. Source and
output chunk counts have a finite budget. No source file read or elapsed-time
measurement enters regrouping.

Source record counts and instruction units add across the selected range.
Native LD replay is retained only at actual read starts, preserving the exact
base record and contiguous prefix bytes. Per-chunk validation conserves aggregate
payload, record forms, decoder units and replay records. Malformed or incomplete
layouts fail instead of being reweighted over an assumed uniform source.

`torch_scan_work` uses each actual shape, requiring its exact compiled geometry.
Its blocks carry file-global ranges and row counts. Reader count is capped by
the actual number of chunks, and initialization/control accounting is carried
through the complete schedule, rather than restarting once per chunk.

Dense writer staging and write events follow the actual input chunks, including
the first payload's automatic block-size choice and the separate df stream.
Host significant-pair and JAGWAS writer services likewise use actual row counts,
survivor scenarios and indexed-part counts. Expansion budgets count actual
chunks before constructing their execution graphs. JAGWAS still uses complete
phenotype panels and one factor per active GPU.

Fixed-run memory entrypoints deliberately reject explicit schedules. Admission
must first use `adaptive_candidate_memory` on the regular fixed-capacity
candidate and its entire admitted size set. Separate peaks from smaller static
candidates are not a bound on mixed old/new allocation lifetimes. This runtime
extension does not relax memory, shape or scientific output requirements.

## Verification

The regression suite compares regrouped chunks against independent direct
PGEN parses for plain, one-bit and difflist LD bases, shifted source extents,
irregular transitions and sample/file tails. It checks source conservation,
malformed schedules, immutable source views, worker initialization, automatic
dense writer block sizes, full/significant/JAGWAS payloads, graph budgets and
the adaptive-memory guard. Existing regular-grid and calibration-cache tests
are included in remote validation.

A100 job `20260921-235233-807617` passed 361 tests, with no failures, errors or
skips, in 98.52 seconds. Artifacts were pulled to
`results/scheduled_census_v2_20260922`. The earlier v1 run had two erroneous test
expectations about reusing a cold versus cached tile graph; these were corrected
without changing runtime accounting.

The replay harness `benchmarks/direct_scheduled_census_20260922.py` consumes the
saved actual ranges from the prior two-GPU native PGEN execution audit. It
independently reparses every range and compares modeled indexed-part counts and
bytes with that run's physical output files. Its component prices and setup
are synthetic controls. It provides source/output accounting evidence, not a
measured runtime prediction or a throughput qualification.

The same job replayed the N=2049, M=4097, K=512, C=2 execution with fixed 512-row
buffers and actual 128/256/512/1-row work. All 17 physical indexed parts matched
the modeled 74,392 bytes and all 4,097 variants. The actual source ranges require
11,866,725 read bytes, versus 7,169,612 on the incorrectly assumed regular grid:
4,697,113 additional bytes. There were four LD restarts on cuda:0 and six on
cuda:1, whereas both regular grids had none. The two shards retained two and one
reader initializations, respectively. This includes contiguous read prefixes;
the decoder's base-only replay units are accounted separately.

The report is in `results/scheduled_execution_replay_v1_20260922/report.json`.
After pulling, all 118 package source hashes and the replay harness hash matched
the local checkout. Its original execution report and input identities remain
recorded separately; the old GPU execution is not relabeled as a current-source
execution. All five source input/metadata hashes matched on the server.

## Remaining live decision work

Explicit ranges describe work, not live executor state. Scoring a subrange from
time zero would still charge new setup/readers and omit outstanding queues,
leases and issued chunks. Public JIT decisions need a continuation model that
carries this state, leaves already-issued ranges untouched and evaluates only
the unpaid work. Full candidate enumeration must not block the first useful
chunk. Planning and information collection must be incremental and cost-aware,
as specified in initial_chunk_calibration_20260922.md. No automatic switching
policy or Bayesian posterior is claimed by this accounting change.

Model-state continuation and one-candidate future-range construction are now
implemented; see candidate_continuation_20260922.md. Their checkpoints preserve
modeled live queues and nominal running service. They are not measurements of
actual executor progress, and the public incremental decision loop remains open.
