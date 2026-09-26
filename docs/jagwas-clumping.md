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

## Standard QC in the wrapper

`run_jagwas_clumping.py` applies two steps to every new scan by default:

- **Phenotype outlier rows (`--phenotype-outlier-sd 5`):** for each phenotype
  panel (or each `--phenotype-group`), residualize every trait on the
  covariates and standardize it. Any sample beyond 5 SD in any trait of that
  panel has its whole row for that panel set to missing.
  - TorchGWAS's phenotype-missingness path mean-imputes those rows. The
    genotype and every other group keep the sample, so there is no union of
    exclusions across groups.
  - `OUTPUT_DIR[/NAME]/excluded_samples.txt` lists the samples, and
    `pipeline.json` reports `excluded_samples`.
- **Trait dropping at VIF > 100 (`--jagwas-min-residual 0.01`):** a trait is
  kept only while at least 1% of its variance is its own, given the traits kept
  before it in the greedy pivoted order. The rounding cutoff still applies on
  top.

Pass `0` to either option to switch it off. `--jagwas-rcond` replaces the trait
threshold with eigen truncation. `--reuse-scan` applies neither step, because it
reuses a scan made with its own settings.

Why these are the defaults: the colleague's collinear imaging panels
(fourier_PE_* and graphunet) looked inflated. They had 662–1,198 loci and hits
down to P ≈ 1e-1155, against 110–240 loci in the reference JAGWAS.

- **Cause: outlier samples.** About 60–140 samples per panel carry many traits
  beyond 5 SD, likely failed image processing. They made the panels'
  low-variance directions heavy-tailed, with median kurtosis 100–1,400. Removing
  their rows made those directions Gaussian.
- **The whole row goes, not the one value.** The traits are near-linear
  combinations of one another, so masking or clipping only the extreme value
  breaks those relations for the sample. Doing that made the low-variance
  directions heavier-tailed still.
- **The raised trait cutoff removes most isolated hits that remain.** In a
  test run that removed the outlier samples from the scan, the rounding cutoff
  still left 94–257 isolated hits per panel (P < 5e-8 with no neighbour at
  P < 1e-5 within 100 kb). VIF ≤ 100 left 18–56, fewer than the reference's
  56–196.
- **Result with the standard defaults:** in one pass over all 22 groups, the
  five collinear panels gave 121–175 loci, top −log10 P of 79–135, and 1–3% of
  loci resting on a single SNP. The reference has 110–240 loci and top
  −log10 P of 59–131.
- **Full-rank panels barely move:** CNN went from 801 to 798 loci, and the
  mesh and nceq panels changed by up to about 40 loci each.
- **Library default unchanged:** TorchGWAS's own default (`run_linear_gwas`
  without `jagwas_min_residual`) is still the rounding cutoff alone.

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

Each run reads the whole BGEN once, and on the shared server that pass is most
of the run.

- **Scan time:** six single-group scans of the 102 GB network BGEN on
  the H100 host took 588, 521, 430, 309, 466 and 147 s. Harmonization plus
  clumping added about 75 s.
- **Cache:** the spread is page cache luck. The busy host evicted a
  once-read region within a few minutes, so a reread usually went back to the
  file server; the 147 s scan followed another pass immediately.

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

**Measured on the H100 host:** all 22 of the colleague's groups (2,655 traits),
with the standard QC, took 4 min 53 s of wall time for the whole command.

| phase | time |
|---|---|
| imports | about 16 s |
| reading the 22 panels in parallel, plus the outlier screen | 11 s |
| the one genotype scan | 253 s |
| clumping after the scan | 13 s |

- **Scan time varies:** across runs, the scan took 191–465 s, depending on how
  much of the BGEN the busy host still had cached. The GPU was shared with
  another user's training job at 99% utilization.
- **Before:** the same groups as separate runs had taken the colleague
  27,076 s of run time. A single group takes about 11 minutes from a cold file.
- **Read each panel once:** reading every panel twice, one file after another,
  had cost 266 s of startup on that host. It is now read once, 8 at a time, with
  byte-identical results.

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
  locus tables for all 22 groups. Against five whole-genome single-group runs,
  locus counts were identical. P values differed by at most 0.003 in log10.
  The only other differences were near-ties: adjacent SNPs in tight LD
  swapping as lead or independent SNP, and one SNP at P = 4.996e-8.

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
