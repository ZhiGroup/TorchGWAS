# JAGWAS plus locus clumping handoff

This branch starts from the current public TorchGWAS `main` branch and adds one
lab workflow wrapper. TorchGWAS performs the current joint-trait reduction
(`reduce="jagwas"`); the wrapper converts its indexed chi-square output into the
harmonized input consumed by the validated local FUMA-style clumper.

## The joint statistic and its degrees of freedom

`reduce="jagwas"` (`torchgwas.jagwas_projection.JagwasReduction`) computes
T = z'R⁻¹z over the traits it keeps, chi-square on r degrees of freedom:

- **z:** the score form, z = t / sqrt(1 + t²/df) = sqrt(df)·r. It is linear
  in the phenotype, so a trait that is a linear combination of others adds
  nothing to T. A quadratic form of t itself does not have that property.
  On near-collinear panels it overstated strong hits by hundreds of chi-square
  units.
- **R:** the FP64 Gram matrix of the scanned (residualised, FP32) phenotype.
- **Kept traits:** the longest greedy pivoted-Cholesky prefix meeting two
  conditions:
  - its estimated null rounding error in T, 2·u·sqrt(N)·sqrt(tr R_S⁻¹), is at
    most `TORCHGWAS_JAGWAS_T_ROUNDING` (default 0.01);
  - it stays within R's FP64 numerical rank.

  Well-conditioned panels keep every trait. Collinear panels drop the
  redundant ones with a warning.
- **df = r:** the manifest's `df` is the kept count r, and `jagwas_rank`
  records the dropped traits with their residual variance and VIF.

The overlapped wrapper learns r through `run_linear_gwas(jagwas_rank_callback=...)`,
which fires after the factor is prepared and before the first chunk. So
p-values use r, never the column count. `pipeline.json` reports `jagwas_df` and
`jagwas_rank`.

## Established discovery inputs

The colleague workflow reads the original BGEN on the network drive:

```text
BGEN:        /path/to/all_filtered.bgen
sample IDs:  /path/to/inputs/T1_discovery_sample_order.npy
phenotype:   /path/to/inputs/T1_pheno_v2_discovery.npy
covariates:  /path/to/inputs/T1_discovery_covarC.npy
MAF:         /path/to/inputs/maf_discovery.npy
LD tables:   /path/to/ld_tables
clumper:     /path/to/local_clumping/fuma_clump.py
```

`--sample-ids` is essential: it selects the discovery cohort from the full
BGEN and fixes the genotype rows to the same order as the phenotype and
covariate arrays. Do not substitute an array in another sample order.

## Minimal full run

Build the optional direct BGEN decoders once in a fresh checkout:

```bash
bash build_bgen_cpu.sh
bash build_direct_bgen.sh
```

Then, from the project root on the GPU server:

```bash
export PYTHONPATH=src:/path/to/deps
export PYTHONHASHSEED=0
python \
  scripts/run_jagwas_clumping.py \
  --genotype /path/to/all_filtered.bgen \
  --genotype-format bgen \
  --sample-ids /path/to/inputs/T1_discovery_sample_order.npy \
  --phenotype /path/to/inputs/T1_pheno_v2_discovery.npy \
  --covariates /path/to/inputs/T1_discovery_covarC.npy \
  --maf /path/to/inputs/maf_discovery.npy \
  --genotype-cache-dir /path/to/bgen-metadata \
  --clump-cache /path/to/clump-data/bgen_discovery_clump_cache \
  --output-dir results/discovery_jagwas_clump \
  --bgen-decode-backend gpu \
  --reader-workers 8 \
  --prefetch-chunks 16 \
  --device cuda:0
```

The two cache arguments avoid reparsing the 8.9-million-row BGEN index and
rebuilding static clumping columns on every run. They are reusable only with
this exact variant order and filtering configuration. If `--clump-cache` is
omitted, the wrapper builds one inside the output directory before the first
scan.

The wrapper deliberately uses the public TorchGWAS API rather than copying the
JAGWAS calculation. Its default path follows the optimized legacy workflow:
it starts a separate clumping process before the scan, preloads local LD
tables, streams each narrow JAGWAS result chunk through `/dev/shm`, and runs
chromosome-level clumping as soon as that chromosome has passed through the
scan. It writes:

- `torchgwas/`: the ordinary TorchGWAS run record and QC metadata;
- `loci/`: local-clumping results;
- `pipeline.json`: row counts, overlapped-stage timings and post-scan wait.

The default does not materialize whole-genome JAGWAS statistics. Add
`--sequential` only when an indexed JAGWAS result and a standalone
`jagwas_clump_input.npz` are specifically needed; that diagnostic path writes,
then rereads, the complete scan and does not overlap clumping.

`PYTHONHASHSEED=0` makes string ordering in the existing clumper reproducible.
The wrapper excludes the extended MHC interval during harmonization, matching
the established analysis. Use `--include-mhc` only for an explicitly different
analysis.

The complete network-BGEN validation retained 1,141,204 rows at P <= 0.05 and
produced 2,539 independent significant SNPs, 748 lead SNPs and 477 genomic risk
loci. All five locus tables were byte-for-byte identical to the prior
sequential result. On the shared A100 server during that validation, the scan
took 242.1 seconds and the overlapped workflow waited another 20.5 seconds
after the scan, for 262.7 seconds total. The server was shared and these times
are an integration measurement, not a hardware benchmark.

## Several phenotype groups in one genotype pass

Each run reads the whole BGEN once. On the shared server that pass is the run.
The 102 GB network BGEN takes about 10 minutes to scan, and harmonization plus
clumping add about 1.5 minutes. Rereading it for the next group does not come
from memory: the page cache on that busy host evicted a once-read region within
a few minutes, so every pass went back to the file server.

Groups on the same samples and covariates can share that pass. Replace
`--phenotype` with one `--phenotype-group NAME=PATH` per group:

```bash
--phenotype-group CNN=$BATCH/phenos/CNN.npy \
--phenotype-group graphunet=$BATCH/phenos/graphunet.npy \
--output-dir $BATCH/runs
```

- **Scan:** TorchGWAS scans the concatenated phenotypes once. Each group gets
  its own joint test (`run_linear_gwas(jagwas_groups=...)`,
  `jagwas_projection.JagwasGroups`) with its own correlation, kept traits and df.
- **Clumping:** in the default overlapped mode every group has its own
  clumping worker, and all of them finish together after the scan.
- **Outputs:** each group's `loci/` and `pipeline.json` go to
  `OUTPUT_DIR/NAME`, which is the layout a per-group run with
  `--output-dir OUTPUT_DIR/NAME` produces. The scan record and the combined
  `pipeline.json` are at the top level.
- **Trait names:** a `NAME.traits.txt` sidecar next to the `.npy`, one name
  per line, names the traits in the rank reports.

A grouped run matches separate runs:

- **Residualization:** each group is residualized by its own call, exactly as
  a run of that group alone would do it. The per-column arithmetic is the same
  either way, but kernels chosen for a different width round differently. On a
  near-collinear group, that rounding alone changed which traits the rank
  cutoff kept: three of fourier_PE_L4_xyz's twelve dropped traits, and
  graphunet kept 99 traits instead of 98. The joint statistic moved by a median
  1.6% and by up to 32%.
- **Kept traits and df:** with per-group residualization they are identical to
  the separate runs.
- **Statistic:** it differs only by the scan's rounding of t, a median
  relative difference below 1e-5 on the colleague's panels.
- **Loci:** overlapped and `--sequential` grouped runs produced byte-identical
  locus tables.

`--sequential` writes one indexed scan with a chi-square column per group
(the manifest's `groups` and `df` list them). `--reuse-scan` with any subset of
the groups reclumps them, for example at another `--lead-p`, with
`--clump-workers` concurrent clumping processes.

## Short integration check

Add this option to the command above:

```text
--variant-range 0:25000
```

and use a separate output directory such as `results/smoke_jagwas_clump`.
This validates BGEN loading, discovery-sample selection, JAGWAS output,
harmonization, and clumping, but its locus count is not a genome-wide result.

Use `--reuse-scan` to reuse an existing compatible `torchgwas/sumstats`
directory and rerun only harmonization and clumping.

## Optional hard-call zstd input

A current-format hard-call zstd store is available at this prefix:

```text
/path/to/clump-data/discovery
```

Its 22,250 stored sample rows exactly match
`T1_discovery_sample_order.npy`. To use it, omit `--sample-ids`,
`--bgen-decode-backend` and `--genotype-cache-dir`, and replace the genotype,
MAF and clumping-cache arguments with:

```bash
--genotype /path/to/clump-data/discovery \
--genotype-format zstd \
--maf /path/to/clump-data/discovery.maf.npy \
--clump-cache /path/to/clump-data/zstd_discovery_clump_cache
```

The bare value passed to `--genotype` is a store prefix. Its `.zst`, index,
manifest, sample, and variant-metadata sidecars must remain together. This
store contains hard calls; use the network BGEN command above when dosage
genotypes are required. Do not use the BGEN clumping cache with this filtered
store: their variant axes differ. For another zstd store, omit `--clump-cache`
on the first run to build an aligned cache inside that run's output directory,
then reuse that directory explicitly in later runs.
