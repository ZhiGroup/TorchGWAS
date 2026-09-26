#!/usr/bin/env bash
# The colleague's phenotype groups in one genotype pass, with the pipeline's
# standard QC (5-SD phenotype rows set to missing per group, traits dropped at
# VIF > 100). See docs/jagwas-clumping.md.
#
#   DEVICE=cuda:0 bash examples/run_batch_grouped.sh
#   OUT=/some/dir DEVICE=cuda:0 bash examples/run_batch_grouped.sh
#
# Every .npy in $BATCH/phenos becomes a group named after the file. Each
# group's loci/, pipeline.json and excluded_samples.txt go to $OUT/<group>/.
set -euo pipefail
PY=python
export PYTHONPATH=/path/to/TorchGWAS/src:/path/to/deps
export PYTHONHASHSEED=0
ROOT=/path/to/shared
PREP=$ROOT/torchgwas/torchgwas_jagwas_nceq_35k_cnnonly
BATCH=${BATCH:-$ROOT/torchgwas/batch_seven_plus_torchgwas}
OUT=${OUT:-$BATCH/runs_grouped}
DEVICE=${DEVICE:-cuda:7}

# Not GROUPS: that name is a bash builtin and silently ignores assignment.
GROUP_ARGS=()
for pheno in "$BATCH"/phenos/*.npy; do
  GROUP_ARGS+=(--phenotype-group "$(basename "$pheno" .npy)=$pheno")
done

cd /path/to/TorchGWAS
"$PY" -u scripts/run_jagwas_clumping.py \
  --genotype "$ROOT/UKB_bgen/step4_hetqc_fusionN_rsid.bgen" --genotype-format bgen \
  --sample-ids "$PREP/sample_order_35k_fusionN_eid_eid.npy" \
  --covariates "$PREP/covar_design_35k_fusionN.npy" \
  --maf /path/to/inputs/maf_discovery.npy \
  --genotype-cache-dir "$PREP/bgen_meta_cache_rsid" \
  --clump-cache "$PREP/full_rsid_bgen/clump_cache" \
  --ld-dir /path/to/local_clumping/ld_tables \
  --bgen-decode-backend gpu --reader-workers 8 --prefetch-chunks 16 \
  --device "$DEVICE" \
  "${GROUP_ARGS[@]}" \
  --output-dir "$OUT"
