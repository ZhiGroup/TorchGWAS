#!/usr/bin/env bash
# JAGWAS + locus clumping for a batch of phenotype groups in one genotype pass,
# with the pipeline's standard QC (docs/jagwas-clumping.md).
#
# usage: run_jagwas_batch.sh -o OUTPUT_DIR [-d DEVICE] [-n] PHENOTYPES... [-- EXTRA_OPTIONS...]
#
#   PHENOTYPES   any mix of
#                  a directory  every *.npy in it is a group named after its file
#                  a .npy file  one group named after the file
#                  NAME=PATH    one group with that name
#                Each (samples x traits) array's rows follow SAMPLE_IDS; an
#                optional NAME.traits.txt beside it names the traits.
#   -o DIR       output directory; each group's loci/, pipeline.json and
#                excluded_samples.txt go to DIR/<group>/
#   -d DEVICE    CUDA device (default cuda:0)
#   -n           print the command instead of running it
#   EXTRA_OPTIONS  passed to run_jagwas_clumping.py, e.g. --lead-p 1e-9,
#                  --phenotype-outlier-sd 0 (keep every row),
#                  --jagwas-rcond 0 (rounding cutoff only), --sequential
#
# The other inputs default to the 35k fusionN discovery set; override any of
# them from the environment: CHECKOUT PREP DEPS CLUMPING_DIR GENOTYPE
# SAMPLE_IDS COVARIATES MAF GENO_CACHE CLUMP_CACHE LD_DIR PY.
# run_jagwas_batch.env beside this script, if present, is sourced first and
# holds the site's paths (untracked; the environment still wins). Anything
# left unset falls back to a /path/to placeholder.
#
# example:
#   bash /path/to/TorchGWAS/scripts/run_jagwas_batch.sh \
#     -o /path/to/shared/torchgwas/new_batch/runs -d cuda:0 \
#     /path/to/shared/torchgwas/new_batch/phenos
set -euo pipefail

SITE_ENV="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/run_jagwas_batch.env"
# shellcheck source=/dev/null
[[ -f $SITE_ENV ]] && . "$SITE_ENV"
CHECKOUT=${CHECKOUT:-/path/to/TorchGWAS}
PREP=${PREP:-/path/to/shared/torchgwas/torchgwas_jagwas_nceq_35k_cnnonly}
DEPS=${DEPS:-/path/to/deps}
CLUMPING_DIR=${CLUMPING_DIR:-/path/to/local_clumping}
PY=${PY:-python}
GENOTYPE=${GENOTYPE:-/path/to/shared/UKB_bgen/step4_hetqc_fusionN_rsid.bgen}
SAMPLE_IDS=${SAMPLE_IDS:-$PREP/sample_order_35k_fusionN_eid_eid.npy}
COVARIATES=${COVARIATES:-$PREP/covar_design_35k_fusionN.npy}
MAF=${MAF:-/path/to/inputs/maf_discovery.npy}
GENO_CACHE=${GENO_CACHE:-$PREP/bgen_meta_cache_rsid}
CLUMP_CACHE=${CLUMP_CACHE:-$PREP/full_rsid_bgen/clump_cache}
LD_DIR=${LD_DIR:-/path/to/local_clumping/ld_tables}

usage() { sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

OUT="" DEVICE=cuda:0 DRY_RUN=0
while getopts "o:d:nh" option; do
  case $option in
    o) OUT=$OPTARG ;;
    d) DEVICE=$OPTARG ;;
    n) DRY_RUN=1 ;;
    h) usage 0 ;;
    *) usage 2 ;;
  esac
done
shift $((OPTIND - 1))
[[ -n $OUT ]] || { echo "error: -o OUTPUT_DIR is required" >&2; usage 2; }

# Not GROUPS: that name is a bash builtin and silently ignores assignment.
GROUP_ARGS=()
add_group() {
  [[ -f $2 ]] || { echo "error: phenotype file not found: $2" >&2; exit 2; }
  GROUP_ARGS+=(--phenotype-group "$1=$2")
}
while (($#)) && [[ $1 != -- ]]; do
  if [[ -d $1 ]]; then
    shopt -s nullglob
    files=("${1%/}"/*.npy)
    shopt -u nullglob
    ((${#files[@]})) || { echo "error: no .npy files in $1" >&2; exit 2; }
    for file in "${files[@]}"; do add_group "$(basename "$file" .npy)" "$file"; done
  elif [[ $1 == *=* ]]; then
    add_group "${1%%=*}" "${1#*=}"
  else
    add_group "$(basename "$1" .npy)" "$1"
  fi
  shift
done
[[ ${1:-} == -- ]] && shift
((${#GROUP_ARGS[@]})) || { echo "error: no phenotype groups given" >&2; usage 2; }

COMMAND=("$PY" -u "$CHECKOUT/scripts/run_jagwas_clumping.py"
  --genotype "$GENOTYPE" --genotype-format bgen
  --sample-ids "$SAMPLE_IDS" --covariates "$COVARIATES" --maf "$MAF"
  --genotype-cache-dir "$GENO_CACHE" --clump-cache "$CLUMP_CACHE" --ld-dir "$LD_DIR"
  --clumping-dir "$CLUMPING_DIR"
  --bgen-decode-backend gpu --reader-workers 8 --prefetch-chunks 16 --device "$DEVICE"
  "${GROUP_ARGS[@]}" --output-dir "$OUT" "$@")

echo "$((${#GROUP_ARGS[@]} / 2)) group(s) -> $OUT" >&2
if ((DRY_RUN)); then
  printf '%q ' "${COMMAND[@]}"; echo
  exit 0
fi
export PYTHONPATH=$CHECKOUT/src:$DEPS
export PYTHONHASHSEED=0
cd "$CHECKOUT"
exec "${COMMAND[@]}"
