#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

while [ ! -f real_case/t1_128_smoke/genotype.npy ]; do
  sleep 5
done

PYTHONPATH=src bash real_case/t1_128_smoke/run_smoke_test.sh > real_case/t1_128_smoke/smoke_test.log 2>&1
