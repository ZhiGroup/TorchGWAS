# Bounded significant-pairs candidate construction

`significant_candidate_space.py` connects finite chunk and phenotype-tile axes
to the existing source/resource significant-host calculator. Each context
supplies its own device set, reader allocation, independent component prices
and shared resource capacities. The builder does not multiply shared CPU,
storage or PCIe capacity when adding a GPU.

The objective is `min_x max_s T_graph(x, s)` over the admitted finite candidates
`x` and explicitly declared survivor/host-sharing scenarios `s`. Candidates must
fit host and per-device memory budgets, reserve requirements, reader limits and
finite graph-expansion budgets. The graph includes setup, repeated genotype
passes for phenotype tiles, native reader lifetimes, host significance
selection, the shared result queue, and one durable indexed-part writer.
Balancing isolated stage durations is not imposed as an equality: shared
resource demands and graph dependencies determine the achievable completion
time. Reports retain resource demand and bottleneck information.

Significant-pairs mode partitions phenotypes and preserves every complete
variant pass and final chunk. It uses the native borrowed dense ring followed
by the host selector; this is not a device-reduction price substitution.
Survivor counts are supplied as scenarios, never inferred from the p-value
threshold. Memory admission always assumes all pairs survive. An infeasible
full phenotype panel is rejected before requiring its compiled geometry.
JAGWAS remains a separate full-panel, variant-partitioned mode.

The JSON calculator entry point accepts `model: detailed_significant_host_space`
with `workload`, `contexts`, `bounds`, `joint`, `output`, `prices` and optional
`significance_threshold`. Bounds contain `chunks`, `trait_blocks` and optional
finite `max_candidates`, `max_candidate_tiles`, `max_census_chunks`. The emitted
API settings retain `reduce: significant`, the threshold, output fields and
required host-selector environment. Dense block coalescing and variant
partitioning are not supported by this significant-host bridge.

## Correctness fix in scenario enumeration

The earlier planner concatenated occupancy and host labels with `:`. Distinct
pairs such as (`dense:x`, `y`) and (`dense`, `x:y`) collided. A later combination
could replace the worst scenario and change the selected candidate. Internal
keys now retain the pair. JSON reports contain a list of separate `occupancy`,
`host` and `estimate` fields. Every scenario score must be finite and positive,
so a hidden NaN cannot disappear inside an aggregate maximum.

Identical source scan-work calculations are reused within one plan across
repeated tiles and scenarios. The scoped cache is discarded when the plan
returns or fails; prices are not carried across independent planning calls.
This keeps exact scores while avoiding repeated source tracing. No speedup
claim is made from this structural reduction alone.

## Validation

A100 job `20260921-220014-756216` passed 85 tests in 215.95 seconds with no skips.
The suite covers bounded construction, the scenario collision changing the
winner, invalid scores/resources, full-panel memory rejection before geometry,
within-plan reuse with unchanged scores, JSON output, and related dense/JAGWAS
candidate and significant archive regressions. Pulled artifacts:
`results/significant_bounded_space_v1_20260922`.

The execution audit is `benchmarks/direct_significant_bounded_execution_20260922.py`.
Its synthetic input has N=2049, M=1025, K=4097 and two covariates. An explicit
per-device allocation budget rejects the full panel while admitting 512- and
1024-phenotype tiles. Chunk sizes are 128/256 on one/two GPUs. Fresh captures
retain duration-free CUDA launch geometry. All service prices are synthetic
controls: this audit checks numerical wiring, not measured ranking or capacity.

The first attempt, job `20260921-220825-760367`, stopped before execution because
the harness supplied the workload-sized setup reference instead of the required
independent tiny reference [32,1,2]. The calculator rejected it. The corrected
harness uses the exact reference-shape contract; no production source changed
for that correction. Its rerun, job `20260921-221028-761058`, exposed an unhandled
compiled tail: N=2049, B=1, K=1, C=2 uses the FP32 NSP GEMV specialization with
32 lanes, block [8,32,1], grid [1,1,8], and a separate reduction kernel. The
previous NSP model recognized only the K=512 layout.

`tensor_service.py` now recognizes this additional captured shape and exact
specialization/launch pair. Its matrix arithmetic and intermediate storage
remain explicitly labeled logical floors; internal padding, physical operand
rereads and synchronization are unresolved. No duration or fitted rate was
added. The captured fixture preserves the source artifact hash and job ID in
`tests/fixtures/significant_singleton_n2049_b1_k1.json`. Negative tests reject
unobserved dimensions, layouts, splits, signatures and reduction launches.
The package source hash changes with this calculator correction; earlier
profile bindings remain historical. Job `20260921-221353-764148` passed all 42
geometry regressions in 91.25 seconds (zero skips), then ran the third numerical
audit attempt. The first configuration retained exactly the independent
reference's 4220 pair identities, but the beta check exposed another harness
error: its reference beta was in raw phenotype units, while the API reports
beta per population standard deviation of the covariate-residualized phenotype.
The reference now divides the independent raw OLS coefficient by that residual
standard deviation. Its t statistic, df and selected pair set are unchanged.
No executor or production statistic was altered for this harness correction.
Artifacts from that attempt and the geometry tests are retained under
`results/significant_bounded_execution_v3_20260922` and
`results/significant_singleton_geometry_v1_20260922`.

Job `20260921-222425-766640` completed all eight configurations with the corrected
independent effect-size reference. The candidate axes were B=128/256,
T=512/1024, and cuda:0 versus cuda:0+cuda:1. Every candidate retained exactly
4220 pairs at p<=0.001, with no duplicated identities, identical integer df,
and the planted pair in both final axes (variant 1024, trait 4096). All pair
identities matched an independent full FP64 OLS calculation. The output pair
SHA256 was `2052f5f9c9a67dcd1cbb7ffe662cd09ac1d2713e83f8877c2b80393f299bc079`.

Across all configurations, maximum absolute t error was 2.577071e-6 and maximum
absolute standardized-beta error was 1.098652e-7. Maximum t difference between
configurations was 3.814698e-6. This is numerical agreement within FP32 rounding;
bit identity of floating-point values across every chunk/tile size is not
claimed. The closest independent statistic was 1.910215e-6 from the significance
boundary, and all eight executions still produced the exact reference pair set.

The allocation cap was 147215908 bytes per GPU (about 140.4 MiB), including a
64 MiB modeled reserve. The largest admitted tiled bound was 126471196 bytes;
the smallest full-panel bound was 167960620 bytes. Both full-panel proposals
were rejected, and every admitted tiled proposal executed under the matching
physical CUDA allocator cap. This is deliberately imposed memory pressure on
a moderate fixture, not a production voxel-scale throughput benchmark.

The input and companion metadata were generated together under
`/data/zxie3/torchgwas_significant_plan_fixture_v4_20260922`; `findmnt` confirmed
`/dev/md0`, XFS, mounted at `/data`. Pulled artifacts are under
`results/significant_bounded_execution_v4_20260922`, including the complete
request, plan, geometry, reference, per-candidate outputs and final report.
All 117 recorded package-source hashes and the benchmark hash matched the local
checkout after retrieval. The finite planner chose candidate 8, but synthetic
prices do not support a measured speed or ranking claim for that choice.

## Integration boundary

The guarded public bridge now connects this significant-host calculator to
opt-in `run_linear_gwas` autotuning, with immutable selector-price binding,
mode/threshold-specific plan caching, live memory admission and indexed writer
routing. See `significant_public_autotune_20260922.md` for its current contract
and separate fresh-process validation. JAGWAS public binding is described in
`jagwas_public_autotune_20260922.md`. The initial-chunk measurement controller is implemented internally
but does not yet choose a public execution candidate automatically. Loaded stage
spans are drift evidence, not independent capacity measurements or a whole-GWAS
timing fit.

Cross-job evidence follows `calibration_cache.py`: structural records require
matching dependencies, empirical records require valid observation age, and
memory/contention observations must be refreshed live. New evidence is
published without modifying older records. Out-of-order job completion cannot
make older measurements take precedence over newer observations. See
`initial_chunk_calibration_20260922.md` for the cache/controller audit.
