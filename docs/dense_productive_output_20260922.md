# Productive output progress for dense scans

Dense output now has an optional `BinarySumstatsWriter.on_write_progress`
callback. `DenseWriteProgress` reports the new row prefix whose stored beta/t
arrays have completed their `write(2)` calls. It joins the per-array byte
prefixes, so queue submission, one array finishing early, and partially written
rows do not trigger a false completion. It handles short writes, staging-block
coalescing, and direct queuing of owned arrays without changing the write block,
flush or fsync policy.

The event is immutable and reports store coordinates, trait range, newly written
statistic bytes, completion time and directory. Tiled/sharded wrappers attach
device identity and translate trait or variant coordinates to the root store.
Progress can split or combine scan chunks. The tracker retains at most three
byte counters and one row prefix; it does not retain a full-job chunk list.

## Timing boundary

The first dense event means beta/t statistics have been written to the operating
system. It does not mean the store is durable or published. The per-variant df
file can still be coalescing many small rows, so its written prefix is reported
separately. Neither df flushing nor an extra fsync is forced merely to obtain a
tuning timestamp. Complete association output still requires the successful
writer close and manifest publication.

`ProductiveTuningRun.output_written` now accepts these events in addition to
indexed JAGWAS/significant completion events. The same bounded planning window
can start after material dense output. Dense ranges are counted as events and
statistic cells; they are not labelled as scan chunks or fsynced parts. Source
reservation coverage, admitted chunk sizes, cost gates and fixed partition
ownership remain unchanged.

Callbacks execute on writer threads and must not wait for executor progress.
Exceptions fail the writer and prevent a completion manifest. Reader iterators
are explicitly closed when the public dense writer exits, including an
asynchronous writer failure. Stream callbacks are released after join to avoid
retaining a closed writer's staging buffers through reference cycles.

The optional callback duration is reported as `progress_callback_seconds` in
the dense writer summary, separate from storage-facing write/fsync service.
It sums time on multiple writer threads and can include lock contention; it is
not the callback's added end-to-end wall time or CPU-only cost. With callbacks
disabled there are no progress hooks or locks on the write path.

## Public calibration audit

The existing `initial_calibration` option now connects to real writer progress
for ordinary dense output, phenotype tiles, variant shards, significant pairs
and JAGWAS. Its `output_progress` report includes first processed/written work,
first material output, first indexed-part fsync where applicable, and bounded
aggregate output counts. Times are relative to API entry. Final report rewriting
and function return occur after the reported snapshot; use an outer API timer
for complete end-to-end comparisons.

These are per-run observations, not independent storage capacities. They do not
overwrite reusable parameter records or renew their observation dates. Dense
first-statistic-write and reduced first-part-fsync measurements must retain their
different labels when compared. Empty significant chunks count as processed
work but not material output.

This closes the dense-output notification gap for productive planning; it does
not yet connect an automatic public chunk-selection policy. Source-model
proposal construction, cost-qualified decisions and changes to tile/GPU
assignment remain separate unfinished parts of the goal. JAGWAS still retains
the full phenotype panel on each active GPU.

## Verification

Job `20260922-031923-878584` passed 144 tests without failures or skips in
44.21 seconds (`results/dense_productive_output_v1_20260922/tests.log`). The new
tests exercise real asynchronous writes, staging, partial rows and short writes,
delayed arrays, callback failures, exact output values, ownership coordinates
and productive-run budget activation. Existing dense/indexed writer, public
calibration and phenotype-tile regressions also passed.

Job `20260922-032112-879080` completed the fresh-process public audit on the
N=2,049, M=4,097, K=512, C=2 native-PGEN fixture. Fifteen successful processes
covered ordinary dense output, phenotype-tiled dense output, variant-sharded
dense output, JAGWAS and significant phenotype tiles: one control, initial
measurement and later cache reuse for each. Every reported numeric/index/df
comparison had zero difference from its matching control. All nine later-run
GPU windows found their matching cached baseline; the original record bytes
and observation times were preserved.

Three additional processes injected JAGWAS writer failure, JAGWAS metadata
failure, and dense writer failure. Each left the prior cache unchanged and
published no new record. The audit also checked that no scan/writer thread
remained active after these failures. The dense failure happened after three
completed write-prefix events, while the reader could still be active.

`results/dense_public_progress_v1_20260922/report.json` contains the aggregate
evidence and per-process reports. Input, association output and the cache were
on verified local XFS `/dev/md0` under `/data`; reports were saved in the shared
project. All 127 package file hashes and the benchmark hash matched the
delivered source.

These runs establish correctness and timing boundaries, not a speedup or an
overhead bound. Each condition had one fresh-process run, and both control and
measured arms retained the same writer telemetry. Calibration context binding
took 1.14–3.34 seconds before useful output, whereas cache lookup took
milliseconds. Binding remains a material startup cost on these small jobs and
is the next issue to investigate before automatic public JIT integration.
