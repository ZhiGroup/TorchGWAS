# Separating selector numerical work from allocation and paging

The experiment now measures resident numerical work without substituting a
different NumPy operation. This improves the measurement boundary but does
not yet qualify production timing prices. No result from this audit has been
published to the reusable calibration cache or installed in an autotuner.

## Controls and ownership

`benchmarks/selector_allocator_control.c` uses NumPy's documented
[data-memory handler interface](https://numpy.org/doc/2.2/reference/c-api/data_memory.html).
It is an explicitly built, experiment-only extension with four modes:

1. Default allocator, with no wrapper.
2. A forwarding wrapper around that same default allocator.
3. A fixed 512-MiB resident arena, initialized before measurements.
4. The default allocator, with output pages touched in the allocation callback
   before NumPy's numerical operation starts using them.

The wrapper records allocation/free callback CPU time separately. Mode 4 also
records page-touch CPU time and thread faults within that callback. Numerical
residual CPU is the outer operation's CPU minus measured callback CPU, rather
than a first-minus-warm correction added to allocation-inclusive prices.
Callback subtraction does not remove every wrapper bookkeeping instruction.
Prefaulting changes cache state and mode 3 changes addresses; these remain
experimental boundaries, not interchangeable production capacity estimates.

Arena reset is forbidden while any array allocated from it remains alive,
including base arrays retained by views after the default handler is restored.
The extension rejects nested scopes, unsupported reallocation and capacity
overflow, records foreign-thread/invalid-free errors, and is limited to one
owning thread per process. The arena is never used by the public GWAS executor.
The mode 4 callback touches only newly allocated uninitialized or zero-filled
storage; it does not alter input arrays. All accepted calls have balanced
ownership counters and zero allocator errors.

## Measurement design

Job `20260922-133430-1058313` measured modes 1-3 above; job
`20260922-134355-1061612` repeated the controls with mode 4 added. Both are
terminal. Fixed generic primitive controls use 256 x 4096 cells, nine named
operation families and declared empty/dense/stride-64 or divisor patterns.
Each pattern has seven repeats, each containing four bulk calls or eight
fixed/empty calls. Mode order is randomized within each repeat.

Unchanged primitives are measured once and shared by the NumPy/native
predicate scenarios; only the predicates have distinct controls. This avoids
treating separate noisy measurements of an identical gather as different
backend capabilities. Values, dtypes and output shape/hash parity were
checked outside timing.

Held-out selectors use a four-buffer pinned input ring, shapes 256 x 4096 and
1024 x 8193, and empty, stride-31 and all-pair retention. Each shape/density
has seven randomized repetitions across all modes and both predicates.
Inputs, output destruction, validation and arena initialization are outside
the selector numerical timing boundary. Returned-array destruction is
recorded separately. No genotype input, whole-GWAS timing or held-out duration
is used to derive a coefficient.

## What the evidence establishes

The arena control greatly reduces paging. In the first experiment, large
dense residual CPU medians were 209.894 ms with the NumPy predicate and
164.904 ms with the native predicate. Almost all corresponding arena calls
had zero minor faults. However, one reference NumPy predicate observation
incurred 1,281 faults despite using the arena. Residency cannot simply be
assumed from the allocator label.

The default-address prefault control in the second experiment provides a
separate comparison. For 1024 x 8193 all-pair output:

| Predicate | Numerical residual CPU ms | Allocation-callback page-touch CPU ms | Median faults inside page touch |
| --- | ---: | ---: | ---: |
| NumPy | 190.925 | 280.420 | 51,187 |
| Native | 180.811 | 262.623 | 51,174 |

These are medians of separate observations/fields; they must not be added as
though they came from one particular call. Page-touch cost includes walking
the pages and changing cache state, not just kernel fault-handler service.
This experiment confirms that treating allocation/page work as part of a
single transferable per-element kernel price is inadequate. It does not
establish a universal per-page price or a production speedup.

Fixed mode 4 primitive prices predicted the large dense numerical residuals
at 230.247/208.912 ms (NumPy/native), about 21%/16% above observation. Other
cases failed more strongly: large sparse observed/predicted ratios were
1.94/2.13, whereas small empty/sparse ratios ranged from 0.28 to 0.46.
All results, including unfavorable cases, remain in the summary.

The current CPU refresh window-stability rule rejected three mode 4 reference
windows: NumPy predicate bulk, native predicate bulk and the fixed empty
flatnonzero call. These diagnostics therefore do not provide an admissible
stable price bank. The summaries report warnings and publish zero records;
the stability tolerances were not widened to pass these observations.

Read-only host observations after the first run found load averages around
90 on the 96-logical-CPU A100 server and high utilization across most of the
12-19 affinity set. Automatic NUMA balancing was enabled and its system-wide
hint-fault counters were increasing. These later observations are consistent
with interference but do not prove the cause of any specific probe sample or
establish load during an earlier measurement. No system settings or other
users' processes were changed. The forwarding/default median ratio itself
varied substantially, so wrapper overhead is also not universally qualified.

## Validation, preservation and next use

The final extension passed 17 tests covering value parity, zero-sized arrays,
view/base lifetimes, reset refusal, calloc, overflow, unsupported realloc,
scope restoration and separate prefault accounting. Without the explicitly
built extension, 21 ordinary model tests passed and its test module skipped
cleanly. Production package code did not change in this experiment.

All observation files keep original timestamps and hashes. The read-only
summaries verify that their inputs remain byte-identical. The mode 4 probe's
exact C source, Python harness, measurement helper and compiled extension are
preserved in `results/selector_allocator_probe_sources_20260922/`, after each
was checked against the original observation hashes. The pulled archive,
input report and all 140 current package source hashes were verified locally.

Principal artifacts:

- `results/selector_resident_allocator_control_20260922/report.json`
- `results/selector_resident_allocator_stability_20260922/report.json`
- `results/selector_prefault_allocator_control_20260922/report.json`
- `results/selector_prefault_allocator_summary_20260922/report.json`
- `results/selector_cpu_state_20260922/report.json`
- `results/selector_allocator_control_checks_v3_20260922/pytest.txt`
- `results/selector_allocator_optional_checks_20260922/pytest.txt`
- `results/selector_allocator_probe_sources_20260922/manifest.json`

The new controls make resident-operation, page-touch and release boundaries
testable. Before using their coefficients for tuning, the fixed probes must
pass stability and interference controls and transfer to held-out work within
an explicit error criterion. Production allocation lifetimes and loaded
multi-GPU throughput still require validation. Frozen H100 evidence and
previous immutable calibration records remain unchanged.

The later [trait-rebase audit](trait_rebase_20260922.md) used these allocator
controls to identify NumPy's temporary elision, correct the calculator's
allocation/primitive accounting and make the producer's single-copy behavior
explicit. It does not supply a qualified replacement price bank.
