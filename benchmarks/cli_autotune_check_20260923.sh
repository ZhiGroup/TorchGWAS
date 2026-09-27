#!/usr/bin/env bash
# End-to-end check of `torchgwas.cli linear --autotune` on a benchmark dataset.
# Usage: cli_autotune_check_20260923.sh <data_dir> <output_dir>
set -euo pipefail
data=$1
out=$2
rm -rf "$out"
export TORCHGWAS_PGEN_BACKEND=native TORCHGWAS_PGEN_PACKED=0 TORCHGWAS_BLOCKING_EVENTS=1 NUMPY_MADVISE_HUGEPAGE=0
"$PY" -m torchgwas linear \
    --genotype "$data/input.pgen" --phenotype "$data/phenotype.npy" --covariates "$data/covariates.npy" \
    --pgen-mode hardcall --compute-dtype float32 --reduce significant --significance-threshold 1e-5 \
    --autotune --autotune-options '{"min_job_seconds": 0}' --output-dir "$out"
"$PY" - "$out" <<'PY'
import json, sys
run = json.load(open(sys.argv[1] + '/run.json'))
tune = run['autotune']
chunk = tune['chunk'] or {}
print('rows', run['n_result_rows'], 'layout', tune['layout'].get('why'),
      'readers', run['reader_workers'], 'prefetch', run['prefetch_chunks'])
print('chunk', chunk.get('state'), chunk.get('reason'), chunk.get('choice'),
      [(s['size'], round(s['rows_per_second'])) for s in chunk.get('segments', [])])
PY
