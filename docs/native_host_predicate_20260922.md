# Optional fused CPU significance predicate

`TORCHGWAS_HOST_PREDICATE=native` selects an ordinary compiled C++ CPU loop
for contiguous native-endian FP32 statistics with scalar or per-variant
thresholds. Build it explicitly before a job:

```sh
PYTHONPATH=src python -m torchgwas.native_host_predicate --build
TORCHGWAS_HOST_PREDICATE=native PYTHONPATH=src python your_scan.py
```

The default remains `numpy`. A requested native build that is missing fails
with the build command; scans never invoke a compiler. Unsupported statistic
dtypes, noncontiguous or unaligned arrays, and per-cell/per-trait thresholds
retain the NumPy implementation. The native implementation releases the GIL,
has no shared scratch space, uses no CUDA or ISA-specific intrinsics, and is
compiled without fast-math. A boolean bitwise AND permits compiler
vectorization while preserving finite-value and inclusive-threshold rules.
The original scalar/per-row FP64 cutoff is rounded upward to FP32 before
comparison. Selected coordinates and payloads still own their arrays.

## Calculator and measurement identity

The enabled implementation is `native_row_flat_v1`, distinct from the NumPy
`bounded_flat_v1`. Full-run planning and initial-chunk window planning both
reject a price bank for the other implementation. The native predicate needs
its own independently observed `predicate_native` CPU primitive. It uses one
native call per source chunk and source logical traffic `5*B*K + 4*B` bytes
for FP32 input, a boolean mask, and one threshold per row. The NumPy bounded
predicate retains its block-call count and `25*B*K` logical bytes. These are
source access counts, not physical DRAM transaction measurements.

The supported native scan produces eligible FP32 row-threshold arrays. General
fallback layouts are not priced as if they ran the fused loop. Existing memory
admission remains conservatively valid for the larger NumPy temporary set;
this change does not claim a tighter allocator bound.

Build directories are keyed by source, binding, architecture and flags. A
completed directory is atomically published and reused without replacing its
bytes or metadata. The execution context also records the SHA256 of the
actual loaded binary. Changed source or binary invalidates that binding.
Immutable CPU measurement records continue to retain their original dates;
enabling this implementation cannot relabel a NumPy price as a native price.
The generic primitive collector now calls the production predicate, records
the native binary identity when enabled, and saves original observation start
and finish timestamps. Its raw observations do not certify a whole model.

## Completed verification

Remote A100 job `20260922-111532-1017300` passed 205 tests. Job
`20260922-111846-1020025` passed another 49 prepared-window tests. These cover
threshold rounding, NaN/infinity handling, empty and tail shapes, FP64 and
strided fallbacks, t-only ownership, concurrent calls, build reuse, changed
bindings, and rejection of legacy prices in both planning paths. The compiler
diagnostic confirms vectorization of the inner loop. Retained evidence is in
`results/native_host_predicate_checks_20260922/`.

Eight public GPU runs then compared the NumPy and native selectors using the
existing local `/data` synthetic PGEN fixture: N=2,049, M=4,097, K=512,
chunk=128. Two-GPU runs used cuda:1/cuda:2 with phenotype tile width 193,
three reader workers, prefetch depth three and a durable indexed writer.

| Output case | Fields | Selected pairs per backend | Native calls |
|---|---|---:|---:|
| One GPU, threshold .02 | beta+t | 41,742 | 33 |
| Two GPUs, threshold .02 | t | 41,742 | 99 |
| Two GPUs, threshold 1 | beta+t | 2,097,664 | 99 |
| Two GPUs, threshold 1e-30 | t | 0 | 99 |

Every persisted coordinate, t statistic, df and optional beta matched exactly
within its pair. All 330 intended native calls used the fused path, including
tail chunks and the final phenotype tile. The report records input digests,
all 139 package source digests, actual binary/context identities and output
paths. The pulled source digests match the editable WSL source. See
`results/native_host_predicate_gwas_20260922/report.json` and
`benchmarks/native_host_predicate_gwas_20260922.py`.

These small public runs have one observation per configuration and substantial
startup/order variation; they do not show an end-to-end speedup. The frozen
large H100 calculator and its results were not modified.

## Performance observations and remaining qualification

The initial scalar prototype and production control are separate experiments.
The initial prototype's whole-selector helper also uses a simpler df gather;
its whole-selector timings must not be attributed solely to the predicate.
Its large validation calls proved unusually expensive: a separate 8M-element
integer `np.testing.assert_array_equal` control took 63.256 seconds, with stack
dumps inside NumPy's NaN/infinity checks. Validation is outside the reported
timing interval. The subsequent production control uses exact dtype, shape
and byte comparisons for its selected arrays.

The production control in `20260922-112445-1021156` overlapped the older
prototype validation on different CPU affinities and cannot establish an
isolated capacity. Both completed. Job `20260922-111846-1020025` also completed
its later production control and independent primitive collector. Subsequent
controls corrected a thread-count mismatch before comparing saved prices to
held-out observations. All results remain separate artifacts; no original
observation timestamp was renewed.

The completed comparison exposed a source-level pricing error: empty, sparse
and dense NumPy index extraction need distinct costs. The calculator now uses
the inspected branch protocol, with a diagnostic repricing from the original
measurements. See [branch pricing and completed H100 controls](nonzero_branch_pricing_20260922.md).
Those eight H100 runs support reduced selector CPU work but do not establish
an overall scan speedup. Full-model runtime accuracy, multi-GPU contention and
profitable productive switching remain separate qualification requirements.
