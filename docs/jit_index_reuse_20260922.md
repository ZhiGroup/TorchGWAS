# Reusing the PGEN index for initial-chunk tuning

The experimental `initial_chunks` route now carries one file-identity-bound
PGEN header from source backend selection through the native readers, memory
admission and the first productive calculator window. A changed file is
rejected before a retained header is reused. Direct callers without a retained
header still parse and validate their own index. After admission, the public
bridge releases the fine-grid layout and retains only the compact source
template and header required for bounded post-output windows. It rechecks live
host and GPU availability before scanning.

This removes duplicate full-index parsing; it does not remove the initial index
parse or the full-file memory-layout scan. Neither path reads genotype payloads
for planning. Admission now also passes its validated, read-only non-LD-base
array to the first calculator window. That window checks the same header object
and current file identity, then prices only the bounded range it needs. A
caller without the admission certificate still performs full-index validation.
The retained base locators use 32-bit PGEN variant indices; their exact byte
count is added to public host-memory admission. Later live checks credit only
that already allocated array, not the process's entire RSS.

On the server-local `/data/zxie3/torchgwas_pgen_benchmark/hardcall_full.pgen`
(8,086,101 variants; 22,250 samples), one fresh-process component breakdown
recorded 10.130 wall/10.109 CPU seconds for `read_header`, 3.523/3.523 seconds
for a 128-marker memory layout using the retained header, and 4.216/4.214
seconds for the still-full-file `PgenHeaderWork` validation using that header.
Separate standalone operations on the same retained header took 1.235 CPU
seconds for distinct type validation, 0.346 for nonempty lengths and 2.006 for
the full non-LD-base locator array. These are ordered observations under one
server load, not cold/warm matched controls or an end-to-end speedup estimate.
Another ordered process took 0.017 seconds to construct a native reader from
the prepared header and 14.495 seconds to construct a fresh reader; its index
and work-bound timings varied substantially from the component breakdown.

After the admission certificate and vectorized chunk-start lookup were added,
a second fresh process measured 5.058 wall/5.055 CPU seconds for the initial
index read, 1.578/1.553 seconds for the 63,173-row 128-marker memory layout,
0.000408/0.000408 seconds to construct the certified work window, and
6.632/6.619 seconds for the same constructor revalidating the retained header
without the certificate. The two constructors returned identical bounded
windows. This is an ordered component comparison under one server load, not a
matched cold-input or full-GWAS speedup comparison.

The reproducible [component harness](../benchmarks/direct_jit_index_reuse_20260922.py)
wrote `results/jit_index_reuse_20260922/report_compact.json` after converting
the retained locators to the production 32-bit representation. It records
24,322,956 retained bytes, 9.587 CPU seconds for the index read, 7.377 for
memory layout, 0.000416 for the certified window constructor and 8.495 for
revalidation using the same open header. The script and three source hashes
match the local code. An earlier run of the harness recorded 10.299 CPU
seconds for memory layout; another ordered check recorded 1.553. This
variation is precisely why the remaining cold pass cannot yet be called cheap
or its benefit inferred from one timing.

The first public two-GPU JAGWAS audit with N=2,049, M=4,097 and K=512 wrote
bit-identical control, deferred and reuse results after this change. The audit
was rerun under `results/public_initial_chunks_20260922/header_reuse_execution/`.
Its deferred/reuse API times were 2.470/2.279 seconds and charged proposal
steps 0.721/1.449 seconds. Both kept chunk size 128 after unstable forecasts.
Preparation took 14.596 seconds outside the API boundary, and supplied rates
were mostly synthetic. These are correctness and overhead observations, not a
production autotuning benefit. The split remote regression run passed 20
native-reader tests (88 subtests) and 180 planner/output tests. The later
certificate and vectorized layout passed 199 targeted tests, including exact
LD restart/payload extents over seven chunk widths.
After the retained-byte accounting change, 60 focused controller execution,
cost and live-capacity cases passed, including the exact host-credit boundary.
An initial broad run exposed controller-only fixtures lacking an admission
memory record; that optional test path now defaults to zero retained bytes.
The final consolidated remote batch passed 200 admission, header-work,
controller and output tests in 60.28 seconds.

The public audit was then rerun after both changes, with the same 4,097-variant
fixture and new exclusive output paths under
`results/public_initial_chunks_20260922/header_reuse_execution_v2/`.
All three runs again wrote bit-identical JAGWAS results. Control, deferred and
reuse API times were 3.274, 6.488 and 4.218 seconds; the two charged proposal
steps were 2.832 and 1.496 seconds and neither applied a switch. Separate
source-open preparation took 3.975 seconds. This run had different load and
warm state from the earlier audit, so the timings do not isolate the index
change or support a speedup claim.

The final public rerun with 32-bit retained-base accounting is under
`results/public_initial_chunks_20260922/header_reuse_execution_v3/`.
Control, deferred and reuse results were again bit-identical. The fixture
retained 36 bytes of base indices and admitted 570,280,088 host bytes including
that exact amount. API times were 8.872, 5.382 and 1.175 seconds; proposal
steps took 2.100 and 0.547 seconds and applied no switch. This large timing
variation across audits on the shared server makes them unsuitable for an
index-reuse speedup estimate.

Next, profile the complete admission path on a realistic large-job context.
The remaining fine-grid layout and mixed-size memory envelope are still
whole-file work, and synthetic public-audit prices cannot qualify a profitable
decision. A compact, mathematically equivalent source-memory envelope may
reduce cold startup further, but must preserve exact LD restart and source
memory upper bounds. Planning steps still need a cheap remaining-work gate so
a short job does not spend more time evaluating a change than it could save.
