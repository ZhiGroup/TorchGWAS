# JAGWAS plus locus clumping handoff

This branch starts from the current public TorchGWAS `main` branch and adds one
lab workflow wrapper. TorchGWAS performs the current joint-trait reduction
(`reduce="jagwas"`); the wrapper converts its indexed chi-square output into the
harmonized input consumed by the validated local FUMA-style clumper.

## Established discovery inputs

The colleague workflow reads the original BGEN on the network drive:

```text
BGEN:        /data4012/zxie3/all_filtered.bgen
sample IDs:  /data484_4/zxie3/torchGWAS1.1/T1_discovery_sample_order.npy
phenotype:   /data484_4/zxie3/torchGWAS1.1/T1_pheno_v2_discovery.npy
covariates:  /data484_4/zxie3/torchGWAS1.1/T1_discovery_covarC.npy
MAF:         /data484_4/zxie3/torchGWAS1.1/maf_discovery.npy
LD tables:   /data/zxie3/ld_tables
clumper:     /data484_4/zxie3/local_clumping/fuma_clump.py
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

Then, from the project root on `lab-a100`:

```bash
export PYTHONPATH=src:/data484_4/zxie3/torchGWAS1.1/.deps
export PYTHONHASHSEED=0
/data4012/zxie3/anaconda3/envs/heart/bin/python \
  scripts/run_jagwas_clumping.py \
  --genotype /data4012/zxie3/all_filtered.bgen \
  --genotype-format bgen \
  --sample-ids /data484_4/zxie3/torchGWAS1.1/T1_discovery_sample_order.npy \
  --phenotype /data484_4/zxie3/torchGWAS1.1/T1_pheno_v2_discovery.npy \
  --covariates /data484_4/zxie3/torchGWAS1.1/T1_discovery_covarC.npy \
  --maf /data484_4/zxie3/torchGWAS1.1/maf_discovery.npy \
  --genotype-cache-dir /data/zxie3/torchgwas-bgen-metadata \
  --clump-cache /data484_4/zxie3/torchGWAS-jagwas-clump-data/bgen_discovery_clump_cache \
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
/data484_4/zxie3/torchGWAS-jagwas-clump-data/discovery
```

Its 22,250 stored sample rows exactly match
`T1_discovery_sample_order.npy`. To use it, omit `--sample-ids`,
`--bgen-decode-backend` and `--genotype-cache-dir`, and replace the genotype,
MAF and clumping-cache arguments with:

```bash
--genotype /data484_4/zxie3/torchGWAS-jagwas-clump-data/discovery \
--genotype-format zstd \
--maf /data484_4/zxie3/torchGWAS-jagwas-clump-data/discovery.maf.npy \
--clump-cache /data484_4/zxie3/torchGWAS-jagwas-clump-data/zstd_discovery_clump_cache
```

The bare value passed to `--genotype` is a store prefix. Its `.zst`, index,
manifest, sample, and variant-metadata sidecars must remain together. This
store contains hard calls; use the network BGEN command above when dosage
genotypes are required. Do not use the BGEN clumping cache with this filtered
store: their variant axes differ. For another zstd store, omit `--clump-cache`
on the first run to build an aligned cache inside that run's output directory,
then reuse that directory explicitly in later runs.
