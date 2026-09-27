# Decoder work before payload measurement

`PgenHeaderWork` bounds native int8 PGEN decoder operations from record types
and lengths. It reads the index once and accepts budgeted variant ranges through
`bounds(start, stop)`. A range represents one independently issued chunk. Its
LD restart includes the contiguous bytes read back to the most recent non-LD
record, but only that base record is decoded; it is not expanded to int8.

The bound has its own schema, `torchgwas.pgen_header_work_bounds.v1`. It is not
an exact census and cannot be passed to `decoder_work` or `torch_scan_work` as
one. Memory-only startup layouts also remain separate. Bounds are useful when
matching exact work evidence is absent, without performing a full payload
census before the job produces results.

## Derivation and scope

For a difflist with E entries, let G = ceil(E/64), D = E-G, and
A = G*id_bytes + max(G-1,0) + ceil(E/4). After removing the one-bit fixed prefix,
the record length L obeys

```
minimum_header_length(E) + A + D <= L <= 5 + A + 5*D.
```

Both sides increase monotonically with E. Two binary searches give conservative
entry-count endpoints without iterating over the samples. They imply intervals
for group counts, integer count V = 1+D, and integer bytes T = L-A. For example,
the number of one-byte integers is at least max(0, 2*V_lower-T_upper). Each
possible integer length is bounded using both the one-byte and five-byte
extremes. The operation bounds also include unknown high bits in the one-bit
tail. Packed copies, fills, inversion, base updates and expansion are fixed by
the sample count and record forms and use the exact calculator's shared ledger.

These bounds apply conditionally to valid basic hardcall records with 1–5-byte
integer encodings and no unused trailing payload, matching the exact decoder
work model's scope. Nonminimal encodings within five bytes are covered. Header
inspection cannot establish payload validity. A native decoder accepting a
malformed padded record is not evidence that the work model supports it.

`price_header_work` uses the existing independent primitive prices. It also
constrains total integer cost by V and T, so the five marginal histogram upper
bounds cannot invent five full sets of integers. A price missing for an operation
that is only *possibly* present still makes complete pricing unavailable.

`native_decoder_service_bounds` uses the same logical traffic approximation as
`torch_scan_work`, including base refresh, packed staging, int8 write allocation
and the separation of buffered-read CPU/copy work from decoding. It returns
read service and CPU/DRAM work intervals. It does not model setup, GPUs, writing
or complete pipeline runtime. Independent primitive prices remain empirical
approximations; these intervals are not hardware time guarantees. Solving two
endpoint graphs does not establish makespan bounds for a contended pipeline.

## Reuse and incremental planning

Header-derived work is structural evidence. Store it with the existing
`CalibrationParameterCache` under `source_work`, binding the input filesystem
identity, range, schema, header/parser and decoder-ledger source identities.
Use the normal input-stability checks at publication and reuse. Filesystem
identity detects ordinary edits; it is not a full-content hash. A changed
dependency requires a new record. The original record remains unchanged.

CPU, GPU, transfer and storage prices are separate empirical records with their
own measurement protocol, relevant hardware/software dependencies and maximum
age. Reusing a price never renews its observation time. Structural work from an
unchanged file need not be discarded merely because an empirical price expired.
Free memory and current contention must still be queried live. Loaded productive
CUDA event intervals remain observations rather than independent capacities.

Each range request defaults to at most 65,536 records and 256 distinct
(record-form, record-length) signatures, plus at most one replayed base. Exceeding
a budget refuses that optional calculation. These are work limits, not wall-time
guarantees; the productive planner must also charge actual CPU and wall time.
Index construction itself is O(file variants) metadata work and must be charged.
No automatic range enumeration occurs. A later productive exact census or
cached exact record may replace a bound in a new plan without altering the
earlier evidence record.

This component is not yet wired into automatic public JIT decisions. It removes
the need to demand exact compressed work for every unobserved range, but bounded
whole-job proposal construction, profitable decisions and public controller
integration still need work. JAGWAS continues to require the full phenotype panel
per active GPU; neither this component nor cached evidence changes output mode,
significance threshold, phenotype ownership or GPU assignment.

## Verification

The first remote pass (`20260922-024136-866807`) passed 196 tests, including
header-I/O guards, comparisons with exact payload censuses across supported
record forms and LD starts, actual native decoding, integer/group boundaries,
and immutable structural-cache reuse.

The final package passed 297 tests without failures or skips in 52.44 seconds
(`20260922-024732-868449`, `results/pgen_header_bounds_v2_20260922/tests.log`).
This includes the shared scan memory equations, tighter integer-cost
constraints, CPU-limited and DRAM-limited scenarios, separate buffered-read
accounting, and existing adaptive admission/JAGWAS/census regressions.

The audit in `results/pgen_header_bounds_v3_20260922/report.json` checks all
59 ranges formed by chunk sizes 128, 256 and 512 on the N=2,049, M=4,097 native
PGEN fixture. Every exact payload count was enclosed, and read/decode byte
extents and base updates agreed. The input mount was verified as local XFS
`/dev/md0` under `/data`. Index setup took 0.826 ms; individual range calculations
took 4.18 ms median and 6.07 ms maximum with CPU affinity 12–19. These are
fixture-specific header arithmetic costs, not a GWAS speedup or a large-file
startup bound. The independent payload censuses used to check the results are
excluded from those step timings.

All 127 package files, the native decoder source and the audit script were
hash-checked against the delivered source. The report also records the native
library hash. No scientific runtime was fitted from this audit, and no
automatic tuning decision or numerical association change was introduced.
