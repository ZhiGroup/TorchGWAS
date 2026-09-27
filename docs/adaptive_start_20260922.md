# Starting productive work without a runtime search

`adaptive_start.prepare_adaptive_start` builds and admits one supplied starting
layout using the existing memory calculator. The first work chunk can be smaller
than the retained ring capacity. Memory admission covers every reachable allowed
chunk/tail shape and mixed old/new result lifetime. The returned partitions bind
directly to `ProductiveTuningRun`.

The caller provides one validated device/reader context, a fixed partition axis
and phenotype tile width, initial and allowed chunk sizes, current resource
budgets, workspace properties, and explicit reserves. The helper
does not enumerate or rank layouts, benchmark associations, or simulate runtime
graphs. Reader budgets, device/host limits and changed inputs are checked before
returning the admitted layout. Missing or ambiguous timing geometry is reported
and prevents ranking, but does not block memory-admitted productive work.
JAGWAS retains the complete phenotype panel on every active device;
its starting layout can partition variants only. Significant-pairs admission
retains the worst-case dense survivor memory allowance.

## Header-only reader memory

The existing exact `pgen_work_census.census` is unsuitable as a mandatory JIT
startup operation: it reads compressed genotype payloads to count difflist
entries and variable integers. Calling it once is still a full extra payload
pass, potentially expensive on a large job.

The new `pgen_memory_layout` reads only the header and record locators. For each
possible smallest-grid chunk start it records the contiguous payload extent,
the distance back to the most recent non-LD record, and the additional packed
base/relative-offset scratch. Coarser chunks sum their constituent payload
lengths and retain only the first start's LD obligation. The existing adaptive
reader envelope then checks every possible aligned start at the largest retained
capacity. It still accounts for grow-only native input and scratch allocations.

Tests compare these extents with the exact payload census at every start on a
fixture containing ordinary and inverted LD records. A guarded file reader
raises if any admission read reaches genotype payloads; mmap is also forbidden
in that test. Thus the no-payload-read claim is checked at the I/O boundary,
rather than inferred from a function name or timing.

Metadata layouts have their own explicit type and omit decoder operation counts.
The decoder and scan runtime calculators, and future-candidate constructor,
reject them as timing evidence. Startup returns `source_layout` for memory and
`source_census=None`. Exact work counts must come from matching cached evidence
or subsequent productive collection. An absent count is never replaced with
zero or a made-up price. A separate
[header work interval component](pgen_header_work_bounds_20260922.md) can now
bound unobserved decoder work without reading payloads. Its typed intervals
are not substituted for exact census records; automatic public planning with
them remains an integration task.

Header parsing and memory arithmetic still precede productive scanning and have
real cost. The helper reports admission duration, one header pass and zero
payload-census passes. The header/layout allocations remain subject to the host
reserve; the source allocation ledger is not a complete RSS or allocator bound.

## Integration boundary

This is the missing startup admission component for deferred public tuning.
It does not yet select an initial layout automatically, validate an empirical
price artifact, attach itself to `run_linear_gwas`, or change GPU/trait ownership
during a job. The execution bridge must bind source/input/library identities,
recheck live resources, attach the productive chunk controller and route writer
completion events to budgeted planning. Exact model work for unobserved source
regions and profitable bounded proposal construction remain necessary before
automatic public decisions are ready.

The existing public `initial_calibration` option continues to collect and save
early productive intervals, as described in
[public initial calibration](public_initial_calibration_20260922.md). The new
startup component does not reinterpret those loaded intervals as independent
hardware capacities.

## Verified execution

The initial admission implementation passed 99 tests but used the exact payload
census; it was replaced before being offered as a JIT startup path. The header-
only version passed 130 admission, candidate-space and continuation tests
(`20260922-021703-859355`). The final early runtime guards and nonblocking
missing-timing-geometry behavior passed 155 focused admission, decoder, census
and schedule tests (`20260922-022212-860782`,
`results/adaptive_start_v3_20260922/tests.log`). No failures or skips occurred.

That final job also ran the two-GPU JAGWAS productive lifecycle audit on the
N=2,049, M=4,097, K=512, C=2 native-PGEN fixture. Input and association output
mounts were verified as local XFS `/dev/md0` under `/data`. The new helper built
one layout with one header pass, zero payload-census passes and zero evaluated
runtime candidates. Memory admission took 0.419 seconds, including first-use
Python work; this is meaningful overhead for such a small scan, not a zero-cost
startup claim.

The returned layout retained a 512-variant ring and began with 128-variant work.
Two scripted productive transitions, 128→256→512, completed all 4,097 variants
exactly once in 18 indexed parts. Maximum absolute FP32 statistic difference
from the fixed-128 control was 0.0001221, within the declared audit tolerance.
The initial associations were retained and every transition applied only to
future unissued reads. JAGWAS retained all 512 phenotypes on both GPUs.

Three alternating matched pairs gave median first-fsynced-part times of
110.33 ms fixed and 112.30 ms tracked; output-inclusive medians were 625.28 and
562.49 ms. Variability was substantial, so these establish lifecycle and output
correctness rather than a throughput gain or an overhead bound. The transitions
used scripted forecasts, not automatic analytical decisions. Admission and
fixture context construction were outside those matched scan timings.

`results/adaptive_start_execution_v2_20260922` contains the startup, memory and
execution reports. All 126 package source hashes, the benchmark and three helper
hashes were checked against the delivered files. The broad calculator/JIT/
chunking/tiling/multi-GPU goal remains active; the integration boundary above
has not been declared complete.
