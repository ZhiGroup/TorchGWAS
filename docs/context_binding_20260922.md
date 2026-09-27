# Calibration context binding, 2026-09-22

The production-chunk calibration cache still records immutable observations,
matches source/input/execution dependencies, and ages empirical evidence from
its original observation time. This change reduces the cost of checking the
GPU part of that binding; it does not cache live device identities or change
measurement freshness.

## Change

`gpu_identity.py` queries the installed NVIDIA management library through a
small ctypes binding. It resolves each handle by the physical UUID supplied by
PyTorch, checks the returned UUID, and reads the PCI address and driver version.
CUDA and NVML indices are never assumed to match. The init/shutdown pair is
balanced, including partial-query failures. Missing libraries/symbols and query
failures fall back to the existing `nvidia-smi` query. Both paths validate and
normalize the PCI address into the same context field.

The ABI follows NVIDIA's [PCI structure](https://docs.nvidia.com/deploy/nvml-api/api/structnvmlPciInfo__t.html),
[device queries](https://docs.nvidia.com/deploy/nvml-api/api/group__nvmlDeviceQueries.html),
and [initialization/cleanup](https://docs.nvidia.com/deploy/nvml-api/api/group__nvmlInitializationAndCleanup.html)
interfaces. No Python NVML package or custom GPU kernel is required.

Canonical CUDA-name validation is now shared through this small module;
standalone context binding no longer imports the full analytical planner just
to validate `cuda:N`. The public API already imports most numerical machinery,
so the standalone import saving must not be counted as a public API speedup.

At this revision, package source and loaded BLAS/OpenMP and PGEN-library contents were hashed.
Thread pools, affinity, settings, mounts and hardware identity are still read
for each binding. No source/library hash cache or metadata-only replacement
was added in this NVML change. Subsequent optional
[digest reuse](binding_digest_reuse_20260922.md) retains fresh file identity
checks and current context reads. Source changes invalidate older bindings.
Live free memory/contention remain separate from reusable measured parameters.
The later [NumPy identity correction](numpy_context_binding_20260922.md) also
binds the NumPy core binary and runtime CPU feature state; these were not
covered by the original version-string check.

## Attribution and comparison

Remote A100 job `20260922-033936-885141` first profiled the original binding.
Its report is `results/context_binding_baseline_v1_20260922/report.json`.
Across three fresh processes and one repeat per process, the `nvidia-smi`
subprocess alone took 0.349–0.721 seconds. Standalone first-call imports also
contributed, particularly the otherwise unnecessary planner/linear import.

Job `20260922-034518-886896` then compared the two identity paths with the
public API imported before timing in **both** arms. Three fresh processes per
path each made an initial binding and one same-process repeat. All 12 complete
contexts were exactly equal, and direct NVML results were independently checked
against the subprocess results in each process. No NVML arm invoked the
`nvidia-smi` fallback within its timed bindings.

| Execution-context wall seconds | Subprocess path | Direct NVML path |
|---|---:|---:|
| First call, range | 1.600–1.720 | 0.570–0.958 |
| First call, median | 1.611 | 0.871 |
| Same-process repeat, range | 0.673–1.140 | 0.187–0.247 |
| Same-process repeat, median | 0.881 | 0.218 |

The table excludes the separately timed package source hash and Python imports;
it includes CUDA initialization on the first call, content hashing of loaded
libraries, mount/topology queries, and instrumentation overhead. Source hashing
varied from 0.031 to 0.671 seconds in this comparison and remains material.
These are six instrumented fresh processes on a shared server, not a controlled
GWAS speedup or a guaranteed startup-overhead bound. The project/source lives
on NFS; association input/output is on local XFS `/dev/md0`. This experiment
does not measure storage throughput.

The comparison artifact is `results/context_binding_nvml_v1_20260922/report.json`,
with package/script hashes and full contexts. Its benchmark is
`benchmarks/direct_context_binding_20260922.py`.

## Validation and remaining work

The same remote job passed **125 tests** in 28.77 seconds. Tests cover physical
UUID ordering, fresh driver changes, wrong handles, missing symbols, each NVML
query failure, balanced cleanup, fallback identity normalization, cache age,
immutable record reuse, source/context invalidation, and calculator selection.
Log: `results/context_binding_nvml_tests_v1_20260922/tests.log`.

Job `20260922-034909-888181` verified the public API with the final 128 package
files: five layouts, each with control/first-cache/later-cache fresh processes,
plus three injected failures. Dense beta/t arrays, JAGWAS chi-square values,
and significant-pair indices/beta/t/df matched the corresponding controls with
maximum difference zero. All nine later per-GPU windows found their original
baseline; those records retained their exact bytes and observation times.
The single-GPU dense window was consistent and published a separate validation
record; the other eight detected loaded-interval drift and published refreshed
evidence. This does not establish a change in independent hardware capacity.

| Public layout | API wall, first/later cache run (s) | First completed output, first/later (s) |
|---|---:|---:|
| Dense, one GPU | 4.006 / 3.636 | 3.065 / 2.992 |
| Dense, phenotype tiles, two GPUs | 3.210 / 4.939 | 1.918 / 3.245 |
| Dense, variant shards, two GPUs | 3.999 / 4.793 | 3.478 / 4.270 |
| JAGWAS, variant shards, two GPUs | 4.631 / 4.133 | 3.756 / 3.336 |
| Significant pairs, phenotype tiles, two GPUs | 4.076 / 4.974 | 1.944 / 1.752 |

Fixture: 2,049 samples, 4,097 variants, 512 phenotypes, 2 covariates; chunk 128,
prefetch 2, reader budget 4, tile width 128 where applicable, significant-pair
threshold 0.01. Initial-calibration binding took 0.428–0.795 seconds across
these ten calibrated API calls. The API wall includes successful output and
calibration publication, excludes imports, and is externally measured. Dense
first-output events follow common beta/t writes and can precede df flushing
and fsync; reduction events follow part-file fsync. Each condition ran once,
so the table is a lifecycle audit, not an end-to-end overhead/speedup estimate.

The injected JAGWAS writer failure after 12 parts, metadata failure after 33
parts, and dense writer failure after 3 progress events all left prior cache
records unchanged, published no new calibration records, and left no active
torchGWAS worker threads. Report and source/script hashes:
`results/nvml_public_calibration_v1_20260922/report.json`.
Both fixture and association-output mounts were verified as local XFS
`/dev/md0`; reports live in the NFS project directory.

Finally, job `20260922-035814-892143` repeated context binding with
`CUDA_VISIBLE_DEVICES=1,0`. Both logical device identities correctly followed
the opposite physical GPU, with matching NVML/subprocess results. Artifact:
`results/context_binding_nvml_v1_20260922/remapped.json`.

This optimization does not complete automatic tuning. The public API collects
bounded useful-chunk evidence but still needs the incremental analytical
decision policy and profitable, memory-admitted chunk/layout/device switches.
JAGWAS remains restricted to full-phenotype execution per active GPU; significant
and dense outputs retain their distinct tiling and output costs.
