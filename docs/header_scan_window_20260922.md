# Bounded source-window scan pricing

The incremental proposal path previously needed an exact full-file payload
census and an execution graph containing every remaining chunk. This change
connects header-only decoder work to the shared scan calculator, allowing a
finite source window to be priced before such a full-job expansion. It is a
component of the forthcoming proposal path, not a completed tuning policy.

`PgenHeaderWork.window(start, stop, chunk_markers)` retains original file
coordinates, input identity and full file variant count. Its default limits
are eight chunks and 65,536 records. It rejects excessive ranges before any
per-chunk work and reads no genotype payload. The existing header index still
costs O(file variants) to construct; that cost has not disappeared or been
moved into a claim of free startup.

`torch_scan_header_work(data, profile, endpoint='lower'|'upper')` feeds these
typed intervals into the same implementation as the exact-census scan
calculator. It retains shape-specific tensor work, JAGWAS projection, H2D/D2H,
result ownership, finish/control work and disclosed uncertainties. Only the
decoder CPU service selects an interval endpoint; read extent and logical
decoder traffic remain the header-derived values. `issued_chunks` advances
the existing first-reader/ring-reuse accounting for a continuation window.
This is the calculator's initialization model, not observed live worker state.

The exact-census entry points explicitly reject these window objects. A
decoder-service interval does **not** imply that simulations of its endpoints
bound pipeline elapsed time under contention. A window also does not establish
that later variants have the same distribution of encoded work. The returned
audit states both limitations. No association duration is accepted or fitted.

The model keeps reader initialization proportional to the original file's
variant count, including for a short window. It preserves base-only LD replay,
contiguous prefix reads, and extra reader scratch. Input changes invalidate an
already-open header index. As with other cached source work, reuse by a later
job requires matching current input and implementation dependencies.

The exact-census production calculator still follows its original path through
the shared implementation. Public GWAS execution and automatic configuration
selection have not been changed by this patch.

## Verification boundary

Tests compare every non-decoder block field against exact payload censuses,
check that exact decoder CPU/traffic falls within the intervals, reject source
binding/budget errors, and verify JAGWAS projection with captured real geometry
and a one-variant tail. Reader continuation and exact-census rejection are also
covered. The arithmetic audit uses the local `/data` PGEN fixture with 2,049
samples and 4,097 variants, 512 phenotypes and 2 covariates, at chunk sizes 128,
256 and 512. It checks beginning, interior and final-tail windows.

The audit's primitive prices are deliberately synthetic, and exact payload
censuses are offline comparison oracles outside the timed window path. Its
timings measure calculator overhead, not GWAS runtime, calibrated hardware
capacity, prediction accuracy, or autotuning benefit. The benchmark records
package, helper, captured-geometry, input and script hashes.

Final remote job `20260922-041744-901542` passed **153 tests in 64.62 seconds**
and checked all nine windows/two endpoints against the exact-count oracle.
Every non-decoder block field matched; exact decoder CPU/traffic was enclosed.
The test log is `results/header_scan_window_v3_20260922/tests.log`; the arithmetic
report is `results/header_scan_window_arithmetic_v2_20260922/report.json`.
The report's 128 package files, benchmark, helpers and captured-geometry hashes
were verified against the delivered source.

On this small fixture, constructing the header index took 0.87 ms. Two-chunk
header windows took 8.58–11.46 ms; single-variant tail windows took 0.70–1.42 ms.
First pricing of each new shape took 58.79–90.76 ms wall and 58.67–86.24 ms CPU.
Subsequent pricing with a retained shape ledger took 3.61–9.52 ms wall. Four
shape ledgers occupied an estimated 605,692 bytes in the bounded session cache.
These are single-process observed costs, not overhead bounds for larger jobs.

The first-shape costs exceed the current default 50 ms total planning CPU
budget. Therefore this change must not be wired to an automatic decision under
that default and declared successful: it would discard an over-budget result.
The integration must reuse valid structural work or explicitly allocate and
charge worthwhile construction across early production chunks. It must also
include header/index work in the charged step, not just the price arithmetic.

## Remaining integration

Window blocks now compose with the existing mode-specific output and shared
resource models; see `prepared_window_model_20260922.md` for verification and
the isolated-window boundary. The live proposal still needs
an explicit forecast for unobserved source/output work, known already-issued
work, and the cost of calculating and switching. JAGWAS must retain its complete
phenotype panel. Memory admission, parameter freshness and the productive
early-chunk budget remain prerequisites to any public automatic switch.
