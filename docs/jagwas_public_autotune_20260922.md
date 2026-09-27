# JAGWAS public autotuning and immutable component reuse

The public run_linear_gwas API now connects reduce="jagwas" with autotune_profile
and autotune_config to the bounded JAGWAS calculator. The search varies chunk
size and calibrated device/reader contexts. Every active GPU retains the complete
phenotype panel and its own factor. Neither phenotype tiling nor phenotype
partitioning is accepted for JAGWAS.

This is an opt-in execution bridge for explicitly supplied independent component
profiles. Binding a profile does not certify its measurement protocol, candidate
ranking, or absolute runtime predictions. The audit continues to report both
selection_validated and runtime_prediction_validated as false.

## Configuration

The configuration has bounds, joint, qc_trait_block, and jagwas_services,
with optional plan_cache_dir. Bounds contain chunks and optional finite
max_candidates, max_candidate_tiles, and max_census_chunks. There is no
trait_blocks or partition_axes field.

The joint map requires host and output-occupancy scenarios, the aggregate reader
budget, host/device capacities, explicit reserves, and a workspace profile for
each calibrated device. Optional search budgets are max_candidates,
max_scenario_evaluations, and max_source_chunks. The shared calculator
constructs preparation separately for each candidate and host scenario.

The qc_trait_block setting bounds input validation work only. It does not
partition the joint statistic. After QC the retained trait count must fit the
residual phenotype rank. The planner performs full-panel memory admission and
repeats live host/device checks before output starts; the API additionally checks
the necessary factor tensor capacity on the actual selected devices.

Output is native FP32 scan with indexed JAGWAS statistics and durable writes.
The effective output contract is store_beta=False even when the API's generic
sumstats_fields default is beta+t: JAGWAS emits variant index and chi-square
statistic. Dense block coalescing and explicit execution overrides are refused.
The significance_threshold argument is refused for autotuned JAGWAS, which emits
every valid joint statistic; that option belongs to significant-pair selection.

## Cross-job calibration

The jagwas_services path names an exact immutable CalibrationParameterCache record:

~~~python
cache.store(
    "cpu_capacity", "jagwas_components",
    {
        "writer_prices": {
            "prices": indexed_writer_primitive_bank,
            "archive": archive_services,
            "queue_cpu_seconds": queue_services,
        },
        "preparation_services": preparation_by_context_and_host_scenario,
    },
    dependencies={
        "source_sha256": profile_source_hashes,
        "execution_context": bound_execution_context,
    },
    provenance=measurement_provenance,
    observed_unix_seconds=measurement_time,
    max_age_seconds=producer_lifetime,
)
~~~

The preparation-service map is the same explicit map used by
bounded_jagwas_plan: per-device scalar/tensor factor arithmetic, common CPU
steps, and finalization steps. GPU factor/setup rates remain separate component
inputs in the detailed profile. This record is not a replacement for those
independent measurements or their validation.

The profile must bind the component artifact's canonical absolute path and raw
file SHA-256. The record reader verifies content identity, component name,
dependencies, provenance, observation/publication order, and producer lifetime.
Checks occur at construction, after input QC, and after planning. Reuse does
not renew the original measurement timestamp. A later record does not silently
replace the exact record named by this job's profile.

Structural records in the general cache have dependency validity without an age
limit. Empirical records require an explicit lifetime. Available memory and
contention are live observations, not reusable capacities. New evidence is
published as a new record; older records are retained.

Plan caching preserves the original analytical ranking. Each job chooses the
first ranked candidate that still fits live memory and records any rejected
choices. A later job can recover the original winner when memory becomes
available. The chosen row and original planned index are both audited.

## Execution fix found by the public audit

The single-device fallback in linear_scan_multigpu previously forwarded an
aggregate reader limit only when the caller originally requested more than one
GPU. With one requested GPU, a source opened before tuning could keep its
original decoder preference, such as 24 workers.

The driver now forwards an explicit reader budget for every active device
count. Direct driver calls retain source preferences while enforcing the cap.
The public autotune adapter also updates its path-backed source's decode
preference after selection, so a higher selected reader count is not capped by
the source's initial default. Two real CUDA tests
exercise a source preferring 24 workers with a one-worker budget and four
prefetch slots: one requested GPU over several chunks, and two requested GPUs
collapsing to one active shard. Both must construct a one-worker native loader
and match independent FP64 joint statistics. A third public API test starts
from a source preferring one worker, selects two, and checks that the native
loader actually uses two workers while preserving the numerical result.

## Verification and remaining work

A100 job 20260921-230527-783928 passed all 122 tests in 132.16 seconds before
the reader-budget fix. Coverage included immutable binding/expiry, source and
mode mismatch, cached live-memory fallback, bounded QC with full-panel output,
significant/dense API regression, JAGWAS candidate construction and preparation.
The artifacts are results/jagwas_public_autotune_v1_20260922.

The first public execution audit, job 20260921-230819-784820, completed four
fresh Python processes on N=2049, M=1025, K=512, C=2 server-local native PGEN
input. It retained every statistic; all four outputs were byte-identical and
matched sampled independent FP64 OLS/correlation calculations within 0.000139
absolute error. A second job reused the same immutable record and disk plan,
and an expired control was rejected before output creation. Those artifacts
remain historical in results/jagwas_public_execution_v1_20260922, because the
audit exposed the reader-budget issue subsequently fixed above.

The final source passed all 154 tests in 19.93 seconds in A100 job
20260921-231706-791999, with zero failures or skips. These include both real
single-active-device reader-cap tests and the public reader-preference update
test. Two warnings report intentional three-reader/two-prefetch limits in the
existing numerical controls. Test artifacts:
results/jagwas_public_autotune_v3_20260922.

The same final-source job completed all seven audit phases. Results are in
results/jagwas_public_execution_v3_20260922/report.json. Four fresh processes
(PIDs 1902317, 1907812, 1910929, 1913629) produced byte-identical arrays for all
1025 variants. The unconstrained context search selected B=128 on cuda:0+cuda:1;
the single-only and dual-only controls also selected B=128. Each GPU had a
physical 2 GiB allocator cap and the corresponding modeled budget.

The repeated unconstrained run reused plan key
7ad73bb71e825128440fd8fc3512769d6eb82629b9c8e4cb270e9d5a1333ca8a
and calibration record
c7f37de3a74ca1b4e64e6880a5fcde60cb2637c9405acc6007e113dfcf2050ca.
Its observation timestamp remained 1790050666.5665438, while its age increased
from 51.9734201431 to 66.5013391972 seconds. The expired control was rejected
before its output directory existed.

The maximum absolute difference from the sampled independent FP64 reference
was 0.00013829697138589836 in every run. The verifier confirmed that all five
input/metadata files were unchanged. After pulling the artifacts, all 117
package-source SHA-256 values and the benchmark SHA-256 matched the local
checkout. The earlier v1/v2 artifacts retain their own original bindings.

The audit uses synthetic service prices and retained duration-free kernel
geometry. It validates execution, cache reuse and numerical results, not
throughput or ranking quality. First-chunk measurements still need to drive
public candidate selection and safe transitions; the existing internal
measurement controller does not yet make those decisions. The full calculator,
autotune and multi-GPU goal remains active.
