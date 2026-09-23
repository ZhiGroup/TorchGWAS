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
  --output-dir results/discovery_jagwas_clump \
  --bgen-decode-backend gpu \
  --reader-workers 8 \
  --prefetch-chunks 16 \
  --device cuda:0
```

The wrapper deliberately uses the public TorchGWAS API rather than copying the
JAGWAS calculation. It writes:

- `torchgwas/`: the ordinary TorchGWAS run record, QC, and indexed JAGWAS
  statistics;
- `jagwas_clump_input.npz`: harmonized, MAF-filtered clumping input;
- `loci/`: local-clumping results;
- `pipeline.json`: row counts and stage timings.

`PYTHONHASHSEED=0` makes string ordering in the existing clumper reproducible.
The wrapper excludes the extended MHC interval during harmonization, matching
the established analysis. Use `--include-mhc` only for an explicitly different
analysis.

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
