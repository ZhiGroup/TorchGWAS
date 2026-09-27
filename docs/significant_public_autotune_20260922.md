# Public significant-pairs autotuning and immutable calibration

`run_linear_gwas` now dispatches opt-in detailed autotuning for
`reduce='significant'` to the bounded host-selector calculator. It selects
variant chunk size, phenotype-tile width and a supplied one/multiple-GPU
context, then uses the indexed significant-pair writer. The selected tile
cannot accidentally route the result to the dense full-output writer.

```python
run_linear_gwas(
    'input.pgen', 'phenotype.npy', 'covariates.npy',
    pgen_mode='hardcall', reduce='significant',
    significance_threshold=0.001,
    autotune_profile='profile.json', autotune_config='bounds.json',
    output_dir='output',
)
```

This bridge currently supports native int8 hardcall PGEN, complete aligned
FP32 phenotypes in file order, the host significance selector, and durable
binary indexed output. Explicit chunk/tile/device overrides conflict with the
autotuner. Both `beta+t` and `t` output fields are supported. Dense output
coalescing, variant partitioning in significant mode, and device-side
significance selection are refused by this bridge. Explicit non-autotuned
execution paths retain their existing behavior.

The configuration contains `bounds`, `joint`, `qc_trait_block`,
`significant_host_prices`, and optional `plan_cache_dir`. Bounds contain
`chunks`, `trait_blocks`, and optional finite search budgets. `joint` requires
named occupancy and host-sharing scenarios, CPU/host/device budgets, explicit
reserves, and CUDA workspace profiles. A QC width no larger than the smallest
candidate tile keeps input checking bounded before execution selection.

## One immutable record per bound selector price bank

`significant_host_prices` is the path of a record published by
`CalibrationParameterCache.store` with:

- Kind `cpu_capacity`, name `significant_host_components`.
- Dependencies exactly `source_sha256` and `execution_context`, equal to the
  independently bound detailed profile and the current execution context.
- Value containing `host_selector`, `predicate_max_cells`, `prices`, `archive`
  and `queue_cpu_seconds`. These are independent selector, archive and queue
  component prices, not complete-association timings.
- Explicit measurement provenance, `observed_unix_seconds`, and a positive
  `max_age_seconds`. For a bank combining several measurements, use the oldest
  observation time and a lifetime appropriate to every constituent.

The record path and file SHA256 must be included in the detailed profile's
`component_artifacts`. Publication is atomic and does not overwrite an older
record. The reader verifies the record's own content identity, exact dependency
match, observation/publication times, and producer-declared age limit. Legacy
records lacking observation time cannot be used for this binding. The profile
binding timestamp does not renew a measurement.

The price record is checked during controller construction, after input QC,
and again after planning immediately before execution. An expired or changed
record stops the run before a writer is created. Publishing a newer record in
the cache does not silently substitute it for the artifact selected by the
profile; new evidence needs an explicit new binding. A normal later job can
reuse the same still-valid record without measuring it again.

This age check currently applies to the bound significant-host service bank.
The rest of the detailed hardware profile retains its existing source/context
checks; identity alone does not certify current contention or rate accuracy.
The initial-chunk drift controller is still separate from public candidate
selection and is not automatically activated by these arguments.

## Cached ranking with fresh admission

The deterministic plan-cache key includes reduction mode, significance
threshold, profile and configuration hashes, workload, output settings, and
stable input identity. A threshold change cannot reuse another threshold's
plan. Cache hits still check calibration age, current source/context, input
identity and live memory.

All significant candidates admitted under the configured budgets remain ranked
in the cached plan. Each run reads current host and relevant GPU availability,
then picks the first ranked candidate that still fits. The cached ranking is
not modified when a lower-memory candidate is selected. The audit records the
original preferred candidate, rejected candidates with their live-memory
reason, and the actual selection. A later job with more available memory can
select the original preferred candidate without rerunning the calculator.

Memory admission still uses all-pairs retention, including output queues and
explicit reserves. Storage/CPU/GPU contention requires fresh evidence outside
these memory checks; this is not a reservation against concurrent processes.

## Validation status

A100 job `20260921-223756-773154` passed 135 tests in 42.15 seconds before the
live-memory fallback extension. These cover public mode dispatch and real CPU
indexed output, both output field sets, threshold-specific plan reuse, record
immutability and expiry, source/context binding, and dense API/cache regressions.
Artifacts: `results/significant_public_autotune_v1_20260922`.

Job `20260921-224419-775673` passed all 137 final integration tests in 43.53
seconds (zero failures/skips), including live-memory fallback without changing
the cached ranking. It then ran
`benchmarks/direct_significant_public_autotune_20260922.py` in fresh processes.
Both public API runs completed on native local PGEN with N=2049, M=1025,
K=4097 and C=2. Each selected B=256, T=512 and cuda:0+cuda:1 under an explicit
181786666-byte per-GPU allocator cap. Eight tiled configurations were admitted;
the two full-panel proposals were rejected by their memory bounds.

| Public run | PID | Plan cache | Calibration age at admission | Retained pairs |
|---|---:|---|---:|---:|
| First | 1517394 | Miss | 31.835 s | 4220 |
| Second | 1523755 | Hit | 55.076 s | 4220 |

Both used record
`d2fcfac0069aa8abdd7e4da926fcc17cb00ef19288a4679e0116830d1e2fca3c`,
with exactly the same observation timestamp (1790048880.032279). The record file
was unchanged; reading it did not renew its age. Pair identities and integer
df exactly matched independent FP64 OLS, with beta expressed per residual
phenotype standard deviation. Maximum absolute t error was 2.274298e-6.
The sorted pair SHA256 was
`2052f5f9c9a67dcd1cbb7ffe662cd09ac1d2713e83f8877c2b80393f299bc079`
in both processes.

The separate expired-evidence control initially used a relative artifact key
in its manually assembled profile, so it stopped at artifact binding before
reaching the age check. Job `20260921-225025-777606` reran only that control and
the final verifier after canonicalizing the path. It confirmed rejection for
expired observation age before the output directory existed. No package source
or completed numerical run was changed for this harness correction. The final
report records the before/after harness hashes and which phases used each.
The exact original harness is preserved, with a verified matching hash, at
`benchmarks/snapshots/significant_public_before_expiry_path_fix_20260922.py`.

Pulled artifacts are `results/significant_public_autotune_v2_20260922` and
`results/significant_public_execution_v1_20260922`. The final verifier confirmed
that all five fixture files were unchanged. All 117 package-source hashes and
the corrected benchmark hash matched the local checkout after retrieval.

These component prices are explicit synthetic controls, with freshly captured
duration-free CUDA geometry. The audit proves public API wiring, numerical
output, cross-process plan/record reuse and expiry behavior. It does not
establish measured throughput or ranking quality, and the API audit continues
to report selection/runtime validation as false.

JAGWAS remains restricted to a complete phenotype panel per GPU; its separate
public binding is described in `jagwas_public_autotune_20260922.md`.
Automatic initial-chunk measurement,
refresh-driven candidate selection, and further format/device-selector coverage
also remain toward the full goal.
