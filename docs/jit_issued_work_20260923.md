# Issued source and GPU workload at a JIT output checkpoint

The post-output issue frontier identifies chunks already reserved by each
original producer. The written-output boundary identifies which of those
chunks have completed an indexed part, or the common beta/t prefix for dense
output. `productive_issued_work` now intersects those two records and inspects
only issued chunks whose matrix/part completion has not been observed. A dense
prefix can end inside a source chunk; the ledger conservatively keeps that
entire chunk. An empty indexed chunk with a completed writer event is removed.

For each retained original chunk, the existing PGEN header inspector counts
read bytes, conditional decoder-unit intervals, LD base replay, unpacked
int8 H2D bytes and the mandatory FP32 scan product. JAGWAS also counts the
full-panel FP64 projection. The ledger preserves the original device and
phenotype range, so significant-pair trait tiles may count independent source
reads and JAGWAS variant shards retain the full phenotype panel. It caps the
inspection at 16 partitions, 64 pending chunks, one million physical source
records and 256 distinct source signatures by default. These budgets are
checked before any source chunk is inspected.

This is a deliberately conservative **full-chunk workload** for stages whose
exact progress is unknown. Some or all read, decode and GPU work may have
already completed. Decoder-unit intervals are conditional on valid payloads;
they are not measured service or an elapsed-time upper bound. The report does
not resolve GPU streams, result queues, partially written NPZ parts, fsync or
manifest publication. It is a bounded input to the future finite schedule,
not a production chunk/tile/GPU switch.

On lab-a100, 28 focused tests passed in `20260923-025236-1330487`,
covering dense partial prefixes, significant-pair trait tiles and JAGWAS
variant shards as well as existing output-boundary and whole-source checks.
A read-only index probe used the server-local XFS file
`/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen`
(22,250 samples, 8,086,101 variants, 20,838,552,600 bytes). It fabricated
64 early 128-marker reservations and one written dense chunk, leaving 63
pending chunks (8,064 markers). Header construction took 1.052 wall / 1.042
CPU seconds; the separate issued-work ledger took 0.692 wall / 0.683 CPU
seconds. Process peak RSS rose from 148,312 to 152,152 KiB during the
ledger. The report counted 17,745,472 indexed read bytes, 179,424,000
nominal H2D bytes and 47,367,936,000 mandatory FP32 product FLOPs.
Artifact: `results/issued_work_probe_20260923/report.json` (job
`20260923-025447-1331516`). These timings include neither PGEN payload
reads nor an actual GWAS run. The roughly 0.7-second ledger belongs in
charged background work, not a cold startup or blocking writer callback.
A second fresh-process probe on the same file and geometry took 0.317 wall /
0.312 CPU seconds for the header and 0.773 wall / 0.773 CPU seconds for the
ledger; every workload count matched the first run. Its
`results/issued_work_probe_20260923/report_v2.json` records the file identity,
interpreter and SHA-256 hashes of the probe, boundary, header and issued-work
code (job `20260923-025919-1331906`). These are two ordered observations on a
shared server, not a speedup comparison. The different header times reinforce
the need to keep measurements and planner costs separate from immutable
source/architecture parameters.
